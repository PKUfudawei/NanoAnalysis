#!/usr/bin/env python3
import awkward as ak
import numpy as np
import os, yaml, random, correctionlib, uproot

from coffea import processor, lumi_tools
from coffea.analysis_tools import PackedSelection, Weights
from coffea.lookup_tools import extractor
from coffea.lookup_tools.dense_lookup import dense_lookup
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory
from coffea.nanoevents.methods.base import NanoEventsArray
from .GenMatch import GenMatch


class Processor(processor.ProcessorABC):
    def __init__(self, mode: str, param_dir: str, outdir: str) -> None:
        super().__init__()
        self.event = None
        self.weight = None  # should be coffea.processor.Weights while self.weight.weight() should be an array
        self.tag = {}
        self.object = {}
        self.variables = {}
        self.cutflow = {}

        self.param_dir = os.path.abspath(param_dir)
        with open(os.path.join(self.param_dir, 'channel.yaml'), 'r', encoding='utf-8') as f:
            self.channel_info = yaml.safe_load(f)
        with open(os.path.join(self.param_dir, 'golden_JSON.yaml'), 'r', encoding='utf-8') as f:
            self.golden_JSON = yaml.safe_load(f)

        if mode.split('_')[0] not in ['data', 'mc']:
            raise ValueError("Processor.__init__(): mode must start with 'data' or 'mc'")
        if mode.split('_')[1] not in self.golden_JSON.keys():
            raise ValueError(f"Processor.__init__(): mode must contain {set(self.golden_JSON.keys())}")
        if mode.split('_')[2] not in set(ak.flatten(self.channel_info.values())):
            raise ValueError(f"Processor.__init__(): mode must end with {set(ak.flatten(self.channel_info.values()))}")
        self.mode = mode  # = '$type_$year_$channel'
        self.sample_type = self.mode.split('_')[0]
        self.year = self.mode.split('_')[1]
        self.channel = self.mode.split('_')[2]
        self.outdir = os.path.abspath(outdir)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def lumi_mask(self, events: NanoEventsArray) -> NanoEventsArray:  # only applied on data
        if self.sample_type == 'mc':  # usually skipped cuz type is restricted to data before executing this function
            return events
        elif self.sample_type == 'data':
            if not os.path.exists(self.golden_JSON[self.year]):
                self.golden_JSON[self.year] = os.path.join(self.param_dir, self.golden_JSON[self.year].split('/')[-1])
            lumi_mask = lumi_tools.LumiMask(self.golden_JSON[self.year])
            data_mask = lumi_mask(events.run, events.luminosityBlock)
            return events[data_mask]

    def add_weight(self) -> ak.Array:
        if self.weight is None:
            self.weight = Weights(len(self.event), storeIndividual=True)
        if self.sample_type == 'data' or len(self.event) == 0:
            return self.weight.weight()
        with open(os.path.join(self.param_dir, 'pile-up.yaml'), 'r', encoding='utf-8') as f:
            self.pile_up = yaml.safe_load(f)

        self.weight.add('genWeight', weight=np.sign(getattr(self.event, 'genWeight', ak.ones_like(self.event.event))))
        correction = correctionlib.CorrectionSet.from_file(self.pile_up['json'][self.year])[self.pile_up['key'][self.year]]
        weight = correction.evaluate(np.array(self.event.Pileup.nPU), "nominal")
        weightUp = correction.evaluate(np.array(self.event.Pileup.nPU), "up")
        weightDown = correction.evaluate(np.array(self.event.Pileup.nPU), "down")
        self.variables['event_pu_nominal'] = weight
        self.variables['event_pu_up'] = weightUp
        self.variables['event_pu_down'] = weightDown
        self.weight.add('pile-up', weight=weight, weightUp=weightUp, weightDown=weightDown)

        return self.weight.weight()

    def pass_cut(self, name: str, cut: ak.Array) -> ak.Array:
        if self.sample_type == 'data':
            self.cutflow[name] = ak.sum(cut)
        elif self.sample_type == 'mc':
            self.cutflow[name] = ak.sum(np.sign(self.event.genWeight[cut]))

        # update events and all objects after passing cut
        self.cutflow['final'] = self.cutflow[name]
        self.event = self.event[cut]
        for obj in self.object:
            self.object[obj] = self.object[obj][cut]
        for var in self.variables:
            self.variables[var] = self.variables[var][cut]
        return cut

    def triggered(self, method: str = 'any') -> ak.Array:
        if method not in ['any', 'all']:
            raise ValueError("Processor.triggered(): level must be in ['any', 'all']")
        with open(os.path.join(self.param_dir, 'triggers.yaml'), 'r', encoding='utf-8') as f:
            self.triggers = yaml.safe_load(f)

        trigger_cuts = PackedSelection()
        for t in self.triggers[self.year]:
            if t in self.event.HLT.fields:
                trigger_cuts.add(name=t, selection=self.event.HLT[t])
            else:
                raise ValueError(f"{t} is not in event.HLT")

        if method == 'any':
            return trigger_cuts.any(*trigger_cuts.names)
        elif method == 'all':
            return trigger_cuts.all(*trigger_cuts.names)

    def filtered(self, method: str = 'all') -> ak.Array:
        if method not in ['any', 'all']:
            raise ValueError("Processor.filtered(): level must be in ['any', 'all']")
        with open(os.path.join(self.param_dir, 'filters.yaml'), 'r', encoding='utf-8') as f:
            self.filters = yaml.safe_load(f)

        filter_cuts = PackedSelection()
        for f in self.filters[self.year]:
            if f in self.event.Flag.fields:
                filter_cuts.add(name=f, selection=self.event.Flag[f])
            else:
                print(f"{f} is not in event.Flag")

        if method == 'any':
            return filter_cuts.any(*filter_cuts.names)
        elif method == 'all':
            return filter_cuts.all(*filter_cuts.names)

    def b_tag(self, level: str = 'medium', reconstruct: bool = False) -> ak.Array:
        if level not in ['loose', 'medium', 'tight']:
            raise ValueError("Processor.b_tag(): level must be in ['loose', 'medium', 'tight']")

        raw_AK4jet = self.event.Jet
        # Working Points -- loose: 0.0490, medium: 0.2783, tight: 0.7100
        # refer to https://gitlab.cern.ch/groups/cms-btv/-/wikis/SFCampaigns/UL2018
        WP = {'loose': 0.0490, 'medium': 0.2783, 'tight': 0.7100}
        self.tag['b-jet'] = (
            (raw_AK4jet.btagDeepFlavB > WP[level]) &
            (self.object['AK8jet'].delta_r(raw_AK4jet) > 0.8 + 0.4)
        )
        if reconstruct:
            self.object['b-jet'] = self.event.Jet[self.tag['b-jet']]
        return self.tag['b-jet']

    def extra_AK4jet_tag(self, reconstruct: bool = False) -> ak.Array:
        raw_AK4jet = self.event.Jet
        # Working Points -- loose: 0.0490, medium: 0.2783, tight: 0.7100
        # refer to https://gitlab.cern.ch/groups/cms-btv/-/wikis/SFCampaigns/UL2018
        self.tag['extra_AK4jet'] = (
            self.object['AK8jet'].delta_r(raw_AK4jet) > 0.8 + 0.4
        )
        if reconstruct:
            self.object['extra_AK4jet'] = self.event.Jet[self.tag['extra_AK4jet']]
        return self.tag['extra_AK4jet']

    def muon_tag(self, reconstruct: bool = False) -> ak.Array:
        raw_muon = self.event.Muon  # (event, muon)
        self.tag['muon'] = (  # (event, boolean)
            # high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)
            (raw_muon.highPtId == 2) &
            (raw_muon.tkRelIso < 0.1) &  # Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt
            (abs(raw_muon.eta) < 2.4) &
            (raw_muon.pt > 20)  # I don't use `muon_corrected_pt` coming from ROOT.RoccoR
        )
        if reconstruct:
            self.object['muon'] = self.event.Muon[self.tag['muon']]
        return self.tag['muon']

    def electron_tag(self, reconstruct: bool = False) -> ak.Array:
        raw_electron = self.event.Electron  # (event, electron)
        self.tag['electron'] = (  # (event, boolean)
            (raw_electron.cutBased_HEEP == True) &  # cut-based HEEP ID
            (abs(raw_electron.eta) < 2.5) &
            (raw_electron.pt > 20)
        )
        if reconstruct:
            self.object['electron'] = self.event.Electron[self.tag['electron']]
        return self.tag['electron']

    def photon_correction(self, photon):
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaSFJSON
        with open(os.path.join(self.param_dir, 'uncertainty/photon_json.yaml'), 'r', encoding='utf-8') as f:
            photon_json = yaml.safe_load(f)
        correction = correctionlib.CorrectionSet.from_file(photon_json['SF'][self.year])
        if self.year.startswith('2016'):
            json_year = self.year + 'VFP'
        else:
            json_year = self.year

        photon['SF_nominal'] = correction['UL-Photon-ID-SF'].evaluate(json_year, 'sf', 'Medium', photon.eta,  photon.pt)
        photon['SF_up'] = correction['UL-Photon-ID-SF'].evaluate(json_year, 'sfup', 'Medium', photon.eta,  photon.pt)
        photon['SF_down'] = correction['UL-Photon-ID-SF'].evaluate(json_year, 'sfdown', 'Medium', photon.eta,  photon.pt)

        return photon

    def photon_tag(self, reconstruct: bool = False) -> ak.Array:
        raw_photon = self.event.Photon  # (event, photon)
        self.tag['photon'] = (  # (event, boolean)
            (raw_photon.pt > 200) &
            ((abs(raw_photon.eta) < 1.4442) | ((abs(raw_photon.eta) > 1.566) & (abs(raw_photon.eta) < 2.4))) &
            raw_photon.mvaID_WP90 &  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariatePhotonIdentificationRun2
            (~raw_photon.pixelSeed) &
            (raw_photon.electronVeto == True)
        )
        if reconstruct:
            self.object['photon'] = self.event.Photon[self.tag['photon']]
        return self.tag['photon']

    def AK8jet_correction(self, AK8jet):
        if len(AK8jet) == 0:
            return AK8jet

        extract = extractor()
        uncertainty_dir = os.path.join(self.param_dir, 'uncertainty', self.year)
        for f in os.listdir(uncertainty_dir):
            if 'AK8' not in f or 'UncertaintySources' in f:
                continue
            extract.add_weight_sets([f'* * {os.path.join(uncertainty_dir, f)}'])
        extract.finalize()
        evaluator = extract.make_evaluator()

        jec_names = dir(evaluator)
        jec_inputs = {name: evaluator[name] for name in jec_names}
        jec_stack = JECStack(jec_inputs)
        name_map = jec_stack.blank_name_map
        name_map['JetPt'] = 'pt'
        name_map['JetMass'] = 'mass'
        name_map['JetEta'] = 'eta'
        name_map['JetPhi'] = 'phi'
        name_map['JetA'] = 'area'
        name_map['ptGenJet'] = 'pt_gen'
        name_map['ptRaw'] = 'pt_raw'
        name_map['massRaw'] = 'mass_raw'
        name_map['Rho'] = 'rho'
    
        AK8jet['is_real'] = (~np.isnan(ak.fill_none(AK8jet.matched_gen.pt, np.nan)))*1
        AK8jet["pt_raw"] = (1 - AK8jet.rawFactor) * AK8jet.pt
        AK8jet["mass_raw"] = (1 - AK8jet.rawFactor) * AK8jet.mass
        AK8jet["pt_gen"] = ak.values_astype(ak.fill_none(AK8jet.matched_gen.pt, 0), np.float32)
        AK8jet["rho"] = ak.broadcast_arrays(self.event.fixedGridRhoFastjetAll, AK8jet.pt)[0]
        corrected_AK8jet = CorrectedJetsFactory(name_map, jec_stack).build(AK8jet).compute()
        AK8jet["pt"] = corrected_AK8jet.pt
        AK8jet["pt_nominal"] = corrected_AK8jet.pt
        AK8jet["mass"] = corrected_AK8jet.mass
        AK8jet["mass_nominal"] = corrected_AK8jet.mass

        self.AK8jet_corrections = set()
        for i in corrected_AK8jet.fields:
            if i.startswith("JES") or i.startswith("JER"):
                self.AK8jet_corrections.add(i)
                AK8jet[f"pt_{i}_up"] = corrected_AK8jet[i].up.pt
                AK8jet[f"pt_{i}_down"] = corrected_AK8jet[i].down.pt
                AK8jet[f"mass_{i}_up"] = corrected_AK8jet[i].up.mass
                AK8jet[f"mass_{i}_down"] = corrected_AK8jet[i].down.mass

        return AK8jet

    def AK8jet_tag(self, reconstruct: bool = False) -> ak.Array:
        raw_AK8jet = self.event.FatJet  # (event, fatjet), >=1 AK8 jet per event
        if self.sample_type == 'mc':
            raw_AK8jet = self.AK8jet_correction(raw_AK8jet)

        self.tag['AK8jet'] = (  # (event, boolean)
            (raw_AK8jet.pt > 250) &
            (abs(raw_AK8jet.eta) < 2.6) &
            (raw_AK8jet.msoftdrop > 30) &  # Corrected soft drop mass with PUPPI
            (raw_AK8jet.jetId & 2 > 0) &
            (raw_AK8jet.jetId & 4 > 0)
            # Jet ID flags bit1 is loose (always false in 2017 since it does not exist),
            # bit2 is tight, bit3 is tightLepVeto
        )
        if reconstruct:
            self.object['AK8jet'] = raw_AK8jet[self.tag['AK8jet']]

        return self.tag['AK8jet']

    def HEM_tag(self) -> ak.Array:  # jet shape: (event, jet), jet = AK8jet or AK4jet
        with open(os.path.join(self.param_dir, '2018_HEM_correction.yaml'), 'r', encoding='utf-8') as f:
            self.HEM_parameters = yaml.safe_load(f)
        if self.year != '2018':
            return ak.full_like(array=self.event.event, fill_value=True, dtype=bool)

        if self.sample_type == 'data':
            event_in_HEM = (
                (self.event.run >= self.HEM_parameters['2018']['RunC']['start']) &
                (self.event.run <= self.HEM_parameters['2018']['RunD']['end'])
            )
        elif self.sample_type == 'mc':  # RunC & RunD is 63.2% of 2018
            event_in_HEM = ak.Array([random.random() < 0.632 for _ in range(len(self.event))])

        jet_in_HEM = (
            (self.object['AK8jet'].eta > self.HEM_parameters['eta']['min'] - 0.4) &  # deltaR=R_jet/2
            (self.object['AK8jet'].eta < self.HEM_parameters['eta']['max'] + 0.4) &
            (self.object['AK8jet'].phi > self.HEM_parameters['phi']['min'] - 0.4) &
            (self.object['AK8jet'].phi < self.HEM_parameters['phi']['max'] + 0.4)
        )
        photon_in_HEM = (
            (self.object['photon'].eta > self.HEM_parameters['eta']['min']) &
            (self.object['photon'].eta < self.HEM_parameters['eta']['max']) &
            (self.object['photon'].phi > self.HEM_parameters['phi']['min']) &
            (self.object['photon'].phi < self.HEM_parameters['phi']['max'])
        )

        return ~(event_in_HEM & (ak.any(jet_in_HEM, axis=1) | ak.any(photon_in_HEM, axis=1)))

    def calculate_PU_SF(self):
        with open(os.path.join(self.param_dir, 'uncertainty/pile-up_files.yaml'), 'r', encoding='utf-8') as f:
            pile_up_files = yaml.safe_load(f)

        PUfunc = {}
        with uproot.open(os.path.join(self.param_dir, 'uncertainty', pile_up_files['mc'][self.year])) as f:
            PUfunc['mc'] = dense_lookup(f['pileup'].values(), f['pileup'].axis(0).edges())
        for k in ('nominal', 'up', 'down'):
            with uproot.open(os.path.join(self.param_dir, 'uncertainty', pile_up_files['data'][self.year][k])) as f:
                PUfunc[f'data_{k}'] = dense_lookup(f['pileup'].values(), f['pileup'].axis(0).edges())

        nTrueInt = self.event.Pileup.nTrueInt
        nMC = PUfunc['mc'](nTrueInt+1)
        nData = PUfunc['data_nominal'](nTrueInt)
        nData_up = PUfunc['data_up'](nTrueInt)
        nData_down = PUfunc['data_down'](nTrueInt)
        return np.divide(nData, nMC), np.divide(nData_up, nMC), np.divide(nData_down, nMC)

    def store_variables(self, vars: dict):
        self.variables['nAK4jet'] = ak.num(self.event.Jet, axis=1)
        self.variables['nExtraAK4jet'] = ak.sum(self.extra_AK4jet_tag(reconstruct=True), axis=1)
        self.variables['event_No.'] = getattr(self.event, 'event', ak.ones_like(self.event.event))
        self.variables['photon-jet_deltaR'] = self.object['AK8jet'].delta_r(self.object['photon'])
        self.variables['MET+photon_mT'] = np.sqrt(
            2 * self.object['photon'].pt * self.event.MET.pt * (1 - np.cos(self.object['photon'].phi - self.event.MET.phi))
        )
        self.variables['MET+AK8jet_mT'] = np.sqrt(
            2 * self.object['AK8jet'].pt * self.event.MET.pt * (1 - np.cos(self.object['AK8jet'].phi - self.event.MET.phi))
        )
        self.variables['nMuon'] = ak.sum(self.muon_tag(reconstruct=False), axis=1)
        self.variables['nElectron'] = ak.sum(self.electron_tag(reconstruct=False), axis=1)

        if self.sample_type == 'mc':
            self.variables['PUWeight_nominal'], self.variables['PUWeight_up'], self.variables['PUWeight_down'] = self.calculate_PU_SF()
            self.variables['LHEScaleWeight'] = self.event.LHEScaleWeight
            self.variables['LHEPdfWeight'] = self.event.LHEPdfWeight
            for i in self.AK8jet_corrections:
                for j in ('up', 'down'):
                    self.object['AK8jet']['pt'] = self.object['AK8jet'][f'pt_{i}_{j}']
                    self.object['AK8jet']['mass'] = self.object['AK8jet'][f'mass_{i}_{j}']
                    self.object['photon+jet'] = self.object['photon'] + self.object['AK8jet']
                    self.variables[f'photon+jet_mass_{i}_{j}'] = self.object['photon+jet'].mass
            self.object['AK8jet']['pt'] = self.object['AK8jet']['pt_nominal']
            self.object['AK8jet']['mass'] = self.object['AK8jet']['mass_nominal']

        self.object['photon+jet'] = self.object['photon'] + self.object['AK8jet']

        for obj, vars in vars.items():
            for var in vars:
                if obj != 'event':
                    array = getattr(self.object[obj], var)
                elif obj == 'event' and '_' in var:
                    array = self.event[var.split('_')[0]]['_'.join(var.split('_')[1:])]
                elif obj == 'event' and var == 'genWeight':
                    array = np.sign(getattr(self.event, 'genWeight', ak.ones_like(self.event.event)))
                elif obj == 'event' and var == 'weight':
                    array = self.add_weight()  # sign(genWeight) * PU_weight
                elif obj == 'event':
                    array = self.event[var]
                if f'{obj}_{var}' not in self.variables:
                    self.variables[f'{obj}_{var}'] = array

        for field in self.event.FatJet.fields:
            if 'ParticleNet' in field or 'deep' in field or 'inclParT' in field or 'particleNet' in field or 'btag' in field:
                self.variables[f'AK8jet_{field}'] = self.object['AK8jet'][field]

    def to_parquet(self, array: ak.Array) -> None:
        start_event = self.event[0]
        name = f'{start_event.run}_{start_event.luminosityBlock}_{start_event.event}'
        ak.to_parquet(array=array, destination=os.path.join(self.outdir, f'{name}.parq'))

    def preselect_HGamma(self):
        # at least pass one trigger
        self.pass_cut(name='triggered', cut=self.triggered(method='any'))

        # pass all needed flags
        self.pass_cut(name='filtered', cut=self.filtered(method='all'))

        # Muon veto
        # self.pass_cut(name='muon-veto', cut=(ak.sum(self.muon_tag(reconstruct=False), axis=1)==0))

        # Electron veto
        # self.pass_cut(name='electron-veto', cut=(ak.sum(self.electron_tag(reconstruct=False), axis=1)==0))

        # Photon == 1
        self.pass_cut(name='photon', cut=(ak.sum(self.photon_tag(reconstruct=True), axis=1) == 1))
        if self.sample_type == 'mc':
            self.object['photon'] = self.photon_correction(self.object['photon'])

        # AK8 jet >=1
        self.pass_cut(name='AK8jet', cut=(ak.sum(self.AK8jet_tag(reconstruct=True), axis=1) > 0))

        # HEM filter
        self.pass_cut(name='2018_HEM_correction', cut=self.HEM_tag())

        # Photon-Jet cleaning, a very special part so no function definition here
        pj_pair = ak.cartesian({'photon': self.object['photon'], 'AK8jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_index_pair = ak.argcartesian({'photon': self.object['photon'], 'AK8jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_dr = pj_pair.photon.delta_r(pj_pair.AK8jet)
        pj_clean = (pj_dr > 1.1)
        photon_index, jet_index = pj_index_pair.photon[pj_clean], pj_index_pair.AK8jet[pj_clean]

        self.object['photon'] = self.object['photon'][photon_index]
        self.object['AK8jet'] = self.object['AK8jet'][jet_index]
        self.pass_cut(name='photon-jet_cleaning', cut=(ak.sum(pj_clean, axis=-1) > 0))

        Higgs_score = self.object['AK8jet']['inclParTMDV1_HbbvsQCD']
        if Higgs_score.ndim == 1:
            return ak.fill_none(Higgs_score, value=False)

        Higgs_candidate = ak.argmax(Higgs_score, axis=1, keepdims=True, mask_identity=False)
        Higgs_candidate = ak.mask(Higgs_candidate, mask=ak.firsts(Higgs_candidate) >= 0, valid_when=True)
        # self.variables['AK8jet_rankToHiggsMass'] = ak.argsort(ak.argsort(abs(self.object['AK8jet'].msoftdrop - 125), axis=1), axis=1)
        self.object['AK8jet'] = self.object['AK8jet'][Higgs_candidate]
        # self.variables['AK8jet_rankToHiggsMass'] = ak.firsts(self.variables['AK8jet_rankToHiggsMass'][Higgs_candidate], axis=1)

        self.object['photon'] = ak.firsts(self.object['photon'], axis=1)  # exactly 1 photon per event
        self.object['AK8jet'] = ak.firsts(self.object['AK8jet'], axis=1)  # exactly 1 Higgs candidate per event

        # b veto
        self.pass_cut(name='b-veto', cut=(ak.sum(self.b_tag(reconstruct=True, level='medium'), axis=1) == 0))

        return self.cutflow['final']

    def process(self, events: NanoEventsArray) -> dict:
        # initialize
        if self.sample_type == 'data':
            self.event = self.lumi_mask(events)
            self.cutflow['n_events'] = len(self.event)
        elif self.sample_type == 'mc':
            self.event = events
            self.cutflow['n_events'] = ak.sum(np.sign(self.event.genWeight))

        # pass pre-selection
        N_preselect = self.preselect_HGamma()
        if N_preselect == 0:
            return self.cutflow

        # store variables
        self.store_variables(vars={
            'AK8jet': {'pt', 'eta', 'phi', 'mass', 'msoftdrop'} | (set(self.object['AK8jet'].fields) - set(self.event.FatJet.fields)),
            'photon': {'pt', 'eta', 'phi', 'mass', 'cutBased', 'sieie'} | (set(self.object['photon'].fields) - set(self.event.Photon.fields)),
            'event': {'MET_pt', 'MET_phi', 'genWeight', 'weight'},
            'photon+jet': {'pt', 'eta', 'phi', 'mass'},
        })

        # gen-macthing
        if self.sample_type == 'mc':
            gen = GenMatch(events=self.event)
            if self.channel in self.channel_info['signal']:
                self.variables.update(gen.HGamma())
            elif self.channel in self.channel_info['fake_photon']:
                self.pass_cut(name='no prompt photon', cut=gen.all_fake_photon())
            elif self.channel in self.channel_info['true_photon']:
                self.pass_cut(name='any prompt photon', cut=gen.any_true_photon())

        # store output
        if self.cutflow['final'] > 0:
            self.to_parquet(array=ak.Array(self.variables))
        return self.cutflow

    def postprocess(self):
        pass
