#!/usr/bin/env python3
import awkward as ak
import numpy as np
import os, yaml, random, correctionlib

from coffea import processor, lumi_tools
from coffea.analysis_tools import PackedSelection, Weights
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
        self.cuts = PackedSelection()
        self._accumulator = processor.defaultdict_accumulator()

        self.param_dir = os.path.abspath(param_dir)
        with open(os.path.join(self.param_dir, 'triggers.yaml'), 'r', encoding='utf-8') as f:
            self.triggers = yaml.safe_load(f)
        with open(os.path.join(self.param_dir, 'golden_JSON.yaml'), 'r', encoding='utf-8') as f:
            self.golden_JSON = yaml.safe_load(f)
        with open(os.path.join(self.param_dir, 'filters.yaml'), 'r', encoding='utf-8') as f:
            self.filters = yaml.safe_load(f)
        with open(os.path.join(self.param_dir, '2018_HEM_correction.yaml'), 'r', encoding='utf-8') as f:
            self.HEM_parameters = yaml.safe_load(f)
        with open(os.path.join(self.param_dir, 'channel.yaml'), 'r', encoding='utf-8') as f:
            self.channel_info = yaml.safe_load(f)
        with open(os.path.join(self.param_dir, 'pile-up.yaml'), 'r', encoding='utf-8') as f:
            self.pile_up = yaml.safe_load(f)

        if mode.split('_')[0] not in ['data', 'mc']:
            raise ValueError("Processor.__init__(): mode must start with 'data' or 'mc'")
        if mode.split('_')[1] not in self.triggers.keys():
            raise ValueError(f"Processor.__init__(): mode must contain {set(self.triggers.keys())}")
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
            lumi_mask = lumi_tools.LumiMask(self.golden_JSON[self.year])
            data_mask = lumi_mask(events.run, events.luminosityBlock)
            return events[data_mask]

    def add_weight(self) -> ak.Array:
        if self.weight is None:
            self.weight = Weights(len(self.event), storeIndividual=True)
        if self.sample_type == 'data' or len(self.event) == 0:
            return self.weight.weight()

        self.weight.add('genWeight', weight=np.sign(getattr(self.event, 'genWeight', ak.ones_like(self.event.event))))

        correction = correctionlib.CorrectionSet.from_file(self.pile_up['json'][self.year])[self.pile_up['key'][self.year]]
        weight = correction.evaluate(np.array(self.event.Pileup.nPU), "nominal")
        weightUp = correction.evaluate(np.array(self.event.Pileup.nPU), "up")
        weightDown = correction.evaluate(np.array(self.event.Pileup.nPU), "down")
        self.weight.add('pile-up', weight=weight, weightUp=weightUp, weightDown=weightDown)

        return self.weight.weight()

    def pass_cut(self, name: str, cut: ak.Array, final: bool = False) -> ak.Array:
        self.cuts.add(name=name, selection=cut)

        # calculate cutflow
        event_flag = self.cuts.all(*self.cuts.names)
        if self.sample_type == 'data':
            self.cutflow[name] = ak.sum(event_flag)
        elif self.sample_type == 'mc':
            self.cutflow[name] = ak.sum(np.sign(self.event.genWeight[event_flag]))

        # update events and all objects after passing cut
        if final:  # if it is final cut, let's make projection to drop unwanted events and objects
            self.cutflow['final'] = self.cutflow[name]
            self.event = self.event[event_flag]
            self.cuts = PackedSelection()
            for obj in self.object:
                self.object[obj] = self.object[obj][event_flag]
            for var in self.variables:
                self.variables[var] = self.variables[var][event_flag]
        else:  # if it is intermediate cut, keep event size but fill unwanted events and objects with None
            self.event = ak.mask(array=self.event, mask=event_flag)
            for obj in self.object:  # keep size: (event, object) in intermediate process
                self.object[obj] = ak.mask(self.object[obj], mask=event_flag, valid_when=True)
            for var in self.variables:
                self.variables[var] = ak.mask(self.variables[var], mask=event_flag, valid_when=True)
        return event_flag

    def triggered(self, method: str = 'any') -> ak.Array:
        if method not in ['any', 'all']:
            raise ValueError("Processor.triggered(): level must be in ['any', 'all']")

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

    def b_tag(self, level: str = 'medium', reco: bool = False) -> ak.Array:
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
        if reco:
            self.object['b-jet'] = self.event.Jet[self.tag['b-jet']]
        return self.tag['b-jet']

    def extra_AK4jet_tag(self, reco: bool = False) -> ak.Array:
        raw_AK4jet = self.event.Jet
        # Working Points -- loose: 0.0490, medium: 0.2783, tight: 0.7100
        # refer to https://gitlab.cern.ch/groups/cms-btv/-/wikis/SFCampaigns/UL2018
        self.tag['extra_AK4jet'] = (
            self.object['AK8jet'].delta_r(raw_AK4jet) > 0.8 + 0.4
        )
        if reco:
            self.object['extra_AK4jet'] = self.event.Jet[self.tag['extra_AK4jet']]
        return self.tag['extra_AK4jet']

    def muon_tag(self, reco: bool = False) -> ak.Array:
        raw_muon = self.event.Muon  # (event, muon)
        self.tag['muon'] = (  # (event, boolean)
            # high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)
            (raw_muon.highPtId == 2) &
            (raw_muon.tkRelIso < 0.1) &  # Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt
            (abs(raw_muon.eta) < 2.4) &
            (raw_muon.pt > 20)  # I don't use `muon_corrected_pt` coming from ROOT.RoccoR
        )
        if reco:
            self.object['muon'] = self.event.Muon[self.tag['muon']]
        return self.tag['muon']

    def electron_tag(self, reco: bool = False) -> ak.Array:
        raw_electron = self.event.Electron  # (event, electron)
        self.tag['electron'] = (  # (event, boolean)
            (raw_electron.cutBased_HEEP == True) &  # cut-based HEEP ID
            (abs(raw_electron.eta) < 2.5) &
            (raw_electron.pt > 20)
        )
        if reco:
            self.object['electron'] = self.event.Electron[self.tag['electron']]
        return self.tag['electron']

    def photon_tag(self, reco: bool = False) -> ak.Array:
        raw_photon = self.event.Photon  # (event, photon), >=1 photon per event
        self.tag['photon'] = (  # (event, boolean)
            (raw_photon.pt > 200) &
            ((abs(raw_photon.eta) < 1.4442) | ((abs(raw_photon.eta) > 1.566) & (abs(raw_photon.eta) < 2.4))) &
            raw_photon.mvaID_WP90 &  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MultivariatePhotonIdentificationRun2
            (~raw_photon.pixelSeed) &
            (raw_photon.electronVeto == True)
        )
        if reco:
            self.object['photon'] = self.event.Photon[self.tag['photon']]
        return self.tag['photon']

    def AK8jet_tag(self, reco: bool = False) -> ak.Array:
        raw_AK8jet = self.event.FatJet  # (event, fatjet), >=1 AK8 jet per event
        self.tag['AK8jet'] = (  # (event, boolean)
            (raw_AK8jet.pt > 250) &
            (abs(raw_AK8jet.eta) < 2.6) &
            (raw_AK8jet.msoftdrop > 30) &  # Corrected soft drop mass with PUPPI
            (raw_AK8jet.jetId & 2 > 0) &
            (raw_AK8jet.jetId & 4 > 0)
            # Jet ID flags bit1 is loose (always false in 2017 since it does not exist),
            # bit2 is tight, bit3 is tightLepVeto
        )
        if reco:
            self.object['AK8jet'] = self.event.FatJet[self.tag['AK8jet']]
        return self.tag['AK8jet']

    def HEM_tag(self) -> ak.Array:  # jet shape: (event, jet), jet = AK8jet or AK4jet
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

    def store_variables(self, vars: dict):
        for obj in vars.keys():
            for var in vars[obj]:
                if obj != 'event':
                    array = getattr(self.object[obj], var)
                elif obj == 'event' and '_' in var:
                    array = self.event[var.split('_')[0]]['_'.join(var.split('_')[1:])]
                elif obj == 'event' and var == 'genWeight':
                    array = np.sign(getattr(self.event, 'genWeight', ak.ones_like(self.event.event)))
                elif obj == 'event' and var == 'weight':
                    array = self.add_weight()  # sign(genWeight) * PU_weight

                self.variables[f'{obj}_{var}'] = array

        for field in self.event.FatJet.fields:
            if 'ParticleNet' in field or 'deep' in field or 'inclParT' in field or 'particleNet' in field or 'btag' in field:
                self.variables[f'AK8jet_{field}'] = self.object['AK8jet'][field]

    def to_parquet(self, array: ak.Array) -> None:
        tokens = self.event.behavior["__events_factory__"]._partition_key.split('/')
        name = '_'.join([(t if 'Events' not in t else 'Events') for t in tokens])
        ak.to_parquet(array=array, where=os.path.join(self.outdir, f'{name}.parq'))

    def preselect_HGamma(self):
        # at least pass one trigger
        self.pass_cut(name='triggered', cut=self.triggered(method='any'))

        # pass all needed flags
        self.pass_cut(name='filtered', cut=self.filtered(method='all'))

        # Muon veto
        # self.pass_cut(name='muon-veto', cut=(ak.sum(self.muon_tag(reco=False), axis=1)==0))

        # Electron veto
        # self.pass_cut(name='electron-veto', cut=(ak.sum(self.electron_tag(reco=False), axis=1)==0))

        # Photon == 1
        self.pass_cut(name='photon', cut=(ak.sum(self.photon_tag(reco=True), axis=1) == 1))

        # AK8 jet >=1
        self.pass_cut(name='AK8jet', cut=(ak.sum(self.AK8jet_tag(reco=True), axis=1) > 0))

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

        Higgs_score = ak.sum([self.object['AK8jet'][f] for f in self.event.FatJet.fields if f.startswith('inclParTMDV1_probH')], axis=0, mask_identity=True)
        if Higgs_score.ndim == 1:
            return ak.fill_none(Higgs_score, value=False)

        Higgs_candidate = ak.argmax(Higgs_score, axis=1, keepdims=True, mask_identity=False)
        Higgs_candidate = ak.mask(Higgs_candidate, mask=ak.firsts(Higgs_candidate) >= 0, valid_when=True)
        self.variables['AK8jet_rankToHiggsMass'] = ak.argsort(ak.argsort(abs(self.object['AK8jet'].msoftdrop - 125), axis=1), axis=1)
        self.object['AK8jet'] = self.object['AK8jet'][Higgs_candidate]
        self.variables['AK8jet_rankToHiggsMass'] = ak.firsts(self.variables['AK8jet_rankToHiggsMass'][Higgs_candidate], axis=1)

        self.object['photon'] = ak.firsts(self.object['photon'], axis=1)  # exactly 1 photon per event
        self.object['AK8jet'] = ak.firsts(self.object['AK8jet'], axis=1)  # exactly 1 Higgs candidate per event
        self.object['photon+jet'] = self.object['photon'] + self.object['AK8jet']

        # b veto
        # final=True means to drop events not passing all selections
        final_cut = self.pass_cut(name='b-veto', cut=(ak.sum(self.b_tag(reco=True, level='medium'), axis=1) == 0), final=True)
        
        # AK4 jets
        self.variables['nAK4jet'] = ak.num(self.event.Jet, axis=1)
        self.variables['nExtraAK4jet'] = ak.sum(self.extra_AK4jet_tag(reco=True), axis=1)

        self.store_variables(vars={
            'AK8jet': {'pt', 'eta', 'phi', 'mass', 'msoftdrop'},
            'photon': {'pt', 'eta', 'phi', 'mass', 'cutBased', 'sieie'},
            'event': {'MET_pt', 'MET_phi', 'genWeight', 'weight'},
            'photon+jet': {'pt', 'eta', 'phi', 'mass'},
        })

        # Additional vars by special computing
        self.variables['event_No.'] = getattr(self.event, 'event', ak.ones_like(self.event.event))
        self.variables['photon-jet_deltaR'] = self.object['AK8jet'].delta_r(self.object['photon'])
        self.variables['MET+photon_mT'] = np.sqrt(
            2 * self.object['photon'].pt * self.event.MET.pt * (1 - np.cos(self.object['photon'].phi - self.event.MET.phi))
        )
        self.variables['MET+AK8jet_mT'] = np.sqrt(
            2 * self.object['AK8jet'].pt * self.event.MET.pt * (1 - np.cos(self.object['AK8jet'].phi - self.event.MET.phi))
        )
        self.variables['nMuon'] = ak.sum(self.muon_tag(reco=False), axis=1)
        self.variables['nElectron'] = ak.sum(self.electron_tag(reco=False), axis=1)

        return final_cut

    def process(self, events: NanoEventsArray) -> dict:
        # initialize
        if self.sample_type == 'data':
            self.event = self.lumi_mask(events)
            self.cutflow['n_events'] = len(self.event)
        elif self.sample_type == 'mc':
            self.event = events
            self.cutflow['n_events'] = ak.sum(np.sign(self.event.genWeight))

        # process
        final_cut = self.preselect_HGamma()

        # gen-macthing
        if self.sample_type == 'mc':
            gen = GenMatch(events=self.event)
            if self.channel in self.channel_info['signal']:
                self.variables.update(gen.HGamma())
            elif self.channel in self.channel_info['fake_photon']:
                final_cut = self.pass_cut(name='no prompt photon', cut=gen.all_fake_photon(), final=True)
            elif self.channel in self.channel_info['true_photon']:
                final_cut = self.pass_cut(name='any prompt photon', cut=gen.any_true_photon(), final=True)

        # store output
        if any(final_cut):
            self.to_parquet(array=ak.Array(self.variables))
        return self.cutflow

    @property  # transform method into attribute and make it unchangable to hide _accumulator
    def accumulator(self):
        return self._accumulator

    def postprocess(self, accumulator):
        return super().postprocess(accumulator)
