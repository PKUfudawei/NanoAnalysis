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


class TriggerProcessor(processor.ProcessorABC):
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
        self.variables['pass_cut'] = self.variables['pass_cut'] * cut
        if self.sample_type == 'data':
            self.cutflow[name] = ak.sum(self.variables['pass_cut'])
        elif self.sample_type == 'mc':
            self.cutflow[name] = ak.sum(np.sign(self.event.genWeight[self.variables['pass_cut']>0]))

        # update events and all objects after passing cut
        self.cutflow['final'] = self.cutflow[name]
        return self.variables['pass_cut']

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

        photon_in_HEM = (
            (self.object['photon'].eta > self.HEM_parameters['eta']['min']) &
            (self.object['photon'].eta < self.HEM_parameters['eta']['max']) &
            (self.object['photon'].phi > self.HEM_parameters['phi']['min']) &
            (self.object['photon'].phi < self.HEM_parameters['phi']['max'])
        )

        return ~(event_in_HEM & ak.any(photon_in_HEM, axis=1))

    def store_variables(self, vars: dict):
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

    def to_parquet(self, array: ak.Array) -> None:
        start_event = self.event[0]
        name = f'{start_event.run}_{start_event.luminosityBlock}_{start_event.event}'
        ak.to_parquet(array=array, destination=os.path.join(self.outdir, f'{name}.parq'))

    def preselect_HGamma(self):
        # at least pass one trigger
        self.pass_cut(name='triggered', cut=self.triggered(method='any'))

        # pass all needed flags
        self.pass_cut(name='filtered', cut=self.filtered(method='all'))

        # Photon == 1
        self.pass_cut(name='photon', cut=(ak.sum(self.photon_tag(reconstruct=True), axis=1) == 1))
        if self.sample_type == 'mc':
            self.object['photon'] = self.photon_correction(self.object['photon'])

        # HEM filter
        self.pass_cut(name='2018_HEM_correction', cut=self.HEM_tag())

        self.object['photon'] = ak.firsts(self.object['photon'], axis=1)  # exactly 1 photon per event

        # Additional vars by special computing
        self.variables['event_No.'] = getattr(self.event, 'event', ak.ones_like(self.event.event))
        self.variables['MET+photon_mT'] = np.sqrt(
            2 * self.object['photon'].pt * self.event.MET.pt * (1 - np.cos(self.object['photon'].phi - self.event.MET.phi))
        )
        
        self.store_variables(vars={
            'photon': {'pt', 'eta', 'phi', 'mass', 'cutBased', 'sieie'} | (set(self.object['photon'].fields) - set(self.event.Photon.fields)),
            'event': {'MET_pt', 'MET_phi', 'genWeight', 'weight'},
        })

        return self.variables['pass_cut']

    def process(self, events: NanoEventsArray) -> dict:
        # initialize
        if self.sample_type == 'data':
            self.event = self.lumi_mask(events)
            self.cutflow['n_events'] = len(self.event)
        elif self.sample_type == 'mc':
            self.event = events
            self.cutflow['n_events'] = ak.sum(np.sign(self.event.genWeight))
        self.variables['pass_cut'] = ak.ones_like(self.event.event)

        # process
        final_cut = self.preselect_HGamma()

        # gen-macthing
        if self.sample_type == 'mc':
            gen = GenMatch(events=self.event)
            if self.channel in self.channel_info['signal']:
                self.variables.update(gen.HGamma())
            elif self.channel in self.channel_info['fake_photon']:
                final_cut = self.pass_cut(name='no prompt photon', cut=gen.all_fake_photon())
            elif self.channel in self.channel_info['true_photon']:
                final_cut = self.pass_cut(name='any prompt photon', cut=gen.any_true_photon())

        # store output
        if any(final_cut):
            self.to_parquet(array=ak.Array(self.variables))
        return self.cutflow

    def postprocess(self):
        pass
