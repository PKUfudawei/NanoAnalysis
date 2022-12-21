import awkward as ak
import numpy as np
import os
import json
import random

from coffea import processor, lumi_tools
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods.base import NanoEventsArray
from .GenMatch import GenMatch


class Processor(processor.ProcessorABC):
    def __init__(self, mode: str, JSONdir: str, outdir: str) -> None:
        super().__init__()
        self.event = None
        self.tag = {}
        self.object = {}
        self.variables = {}
        self.cutflow = {}
        self.cuts = PackedSelection()
        self._accumulator = processor.defaultdict_accumulator()
        
        self.JSONdir = os.path.abspath(JSONdir)
        self.outdir = os.path.abspath(outdir)
        
        if mode.split('_')[0] not in ['data', 'mc']:
            raise ValueError("Processor.__init__(): mode must start with 'data' or 'mc'")
        self.mode = mode  # = '$type_$year_$channel'
        self.sample_type = self.mode.split('_')[0]
        self.year = self.mode.split('_')[1].replace('APV', '')
        self.channel = self.mode.split('_')[2]
        
        with open(os.path.join(self.JSONdir, 'triggers.json'), 'r', encoding ='utf-8') as f:
            self.triggers = json.load(f)
        with open(os.path.join(self.JSONdir, 'golden_JSON.json'), 'r', encoding ='utf-8') as f:
            self.golden_JSON = json.load(f)
        with open(os.path.join(self.JSONdir, 'filters.json'), 'r', encoding ='utf-8') as f:
            self.filters = json.load(f)
        with open(os.path.join(self.JSONdir, '2018_HEM_correction.json'), 'r', encoding ='utf-8') as f:
            self.HEM_parameters = json.load(f)
    
    def lumi_mask(self, events: NanoEventsArray) -> NanoEventsArray:  # only applied on data
        if self.sample_type=='mc':  # usually skipped cuz type is restricted to data before executing this function
            return events
        elif self.sample_type=='data':
            lumi_mask = lumi_tools.LumiMask(self.golden_JSON[self.year])
            data_mask = lumi_mask(events.run, events.luminosityBlock)
            return events[data_mask]
        
    def pass_cut(self, cutName: str, cut: ak.Array, final: bool = False) -> ak.Array:
        if self.cuts.names and len(cut) != len(self.cuts.all()):
            self.cuts = PackedSelection()
        self.cuts.add(cutName, cut)
        
        # calculate cutflow
        event_flag = self.cuts.all(*self.cuts.names)
        if self.sample_type == 'data':
            self.cutflow[cutName] = ak.sum(event_flag)
        elif self.sample_type == 'mc':
            self.cutflow[cutName] = ak.sum(np.sign(self.event.genWeight[event_flag]))
        
        # update events and all objects after passing cut
        if final:  # if it is final cut, let's make projection to drop unwanted events and objects
            self.event = self.event[event_flag]
            for obj in self.object:
                self.object[obj] = self.object[obj][event_flag]
        else:  # if it is intermediate cut, keep event size but fill unwanted events and objects with None
            self.event = ak.mask(array=self.event, mask = event_flag)
            for obj in self.object:  # keep size: (event, object) in intermediate process
                self.object[obj] = ak.mask(self.object[obj], mask = event_flag)
        return event_flag
    
    def triggered(self, level: str = 'any') -> ak.Array:
        if level not in ['any', 'all']:
            raise ValueError("Processor.triggered(): level must be in ['any', 'all']")
        elif level == 'any':  # pass any trigger
            return ak.any([self.event.HLT[t] for t in self.triggers[self.year] if t in self.event.HLT.fields], axis=0)
        elif level == 'all':  # pass all triggers
            return ak.all([self.event.HLT[t] for t in self.triggers[self.year] if t in self.event.HLT.fields], axis=0)
    
    def filtered(self, level: str = 'all') -> ak.Array:
        if level not in ['any', 'all']:
            raise ValueError("Processor.filtered(): level must be in ['any', 'all']")
        elif level == 'any':  # pass any filter
            return ak.any([self.event.Flag[t] for t in self.filters[self.year] if t in self.event.Flag.fields], axis=0)
        elif level == 'all':  # pass all filters
            return ak.all([self.event.Flag[t] for t in self.filters[self.year] if t in self.event.Flag.fields], axis=0)
    
    def b_tag(self, level: str = 'tight', reco: bool = False) -> ak.Array:
        if level not in ['loose', 'medium', 'tight']:
            raise ValueError("Processor.b_tag(): level must be in ['loose', 'medium', 'tight']")
        
        raw_AK4jet = self.event.Jet
        # Working Points -- loose: 0.0490, medium: 0.2783, tight: 0.7100
        # refer to https://gitlab.cern.ch/groups/cms-btv/-/wikis/SFCampaigns/UL2018
        WP = {'loose': 0.0490, 'medium': 0.2783, 'tight': 0.7100}
        self.tag['b-jet'] = (raw_AK4jet.btagDeepFlavB > WP[level])
        if reco:
            self.object['b-jet'] = self.event.Jet[self.tag['b-jet']]
        return self.tag['b-jet']
    
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
            (raw_photon.mvaID_WP90 > 0.2) &
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
            return ak.Array([True for _ in range(len(self.event))])
         
        if self.sample_type == 'data':
            event_in_HEM = (
                (self.event.run >= self.HEM_parameters['2018']['RunC']['start']) &
                (self.event.run <= self.HEM_parameters['2018']['RunD']['end'])
            )
        elif self.sample_type == 'mc':  # RunC & RunD is 63.2% of 2018
            event_in_HEM = ak.Array([random.random()<0.632 for _ in range(len(self.event))])
            
        jet_in_HEM = (
            (self.object['AK8jet'].eta > self.HEM_parameters['eta']['min'] - 0.4) &
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
                if obj!='event':
                    array = getattr(self.object[obj], var)
                elif obj=='event' and '_' in var:
                    array = self.event[var.split('_')[0]]['_'.join(var.split('_')[1:])]
                elif obj=='event' and var=='genWeight':
                    array = getattr(self.event, var, ak.ones_like(self.event.run))
                    
                self.variables[f'{obj}_{var}'] = array
    
    def to_parquet(self, array: ak.Array) -> None:
        output_dir = os.path.abspath(self.outdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        tokens = self.event.behavior["__events_factory__"]._partition_key.split('/')
        name = '_'.join([(t if 'Events' not in t else 'Events') for t in tokens])
        
        ak.to_parquet(array=array, where=os.path.join(output_dir, f'{name}.parq'))
    
    def preselect_HGamma(self):  # __ in prefix means private method
        # at least pass one trigger
        self.pass_cut(cutName='triggered', cut=self.triggered(level='any'))
        
        # pass all needed flags
        self.pass_cut(cutName='filtered', cut=self.filtered(level='all'))
        
        # b veto
        self.pass_cut(cutName='b-veto', cut=(ak.sum(self.b_tag(reco=True, level='tight'), axis=1)==0))

        # Muon veto
        # self.pass_cut(cutName='muon-veto', cut=(ak.sum(self.muon_tag(reco=False), axis=1)==0))
        
        # Electron veto
        # self.pass_cut(cutName='electron-veto', cut=(ak.sum(self.electron_tag(reco=False), axis=1)==0))
        
        # Photon == 1
        self.pass_cut(cutName='photon', cut=(ak.sum(self.photon_tag(reco=True), axis=1)==1))
        
        # AK8 jet >=1
        self.pass_cut(cutName='AK8jet', cut=(ak.sum(self.AK8jet_tag(reco=True), axis=1)>0))
        
        # HEM filter
        self.pass_cut(cutName='2018_HEM_correction', cut=self.HEM_tag())
        
        # Photon-Jet cleaning, a very special part so no function definition here
        pj_pair = ak.cartesian({'photon': self.object['photon'], 'AK8jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_index_pair = ak.argcartesian({'photon': self.object['photon'], 'AK8jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_dr = pj_pair.photon.delta_r(pj_pair.AK8jet)
        pj_clean = (pj_dr > 1.1)
        photon_index, jet_index = pj_index_pair.photon[pj_clean], pj_index_pair.AK8jet[pj_clean]
        self.object['photon'] = self.object['photon'][photon_index]
        self.object['AK8jet'] = self.object['AK8jet'][jet_index]
        self.object['AK8jet'] = self.object['AK8jet'][ak.argmin(abs(self.object['AK8jet'].msoftdrop - 125), axis=1, keepdims=True)]
        
        # final=True means to drop events not passing all selections
        final_cut = self.pass_cut(cutName='photon-jet_cleaning', cut=(ak.sum(pj_clean, axis=-1)>0), final=True)
        self.object['photon'] = self.object['photon'][:, 0]
        self.object['AK8jet'] = self.object['AK8jet'][:, 0]
        self.object['photon-jet'] = self.object['photon'] + self.object['AK8jet']
        """
        # Return vars of objects after pre-selection
        self.object['heaviest_jet'] = self.object['AK8jet'][ak.argmax(self.object['AK8jet'].msoftdrop, axis=1, keepdims=True)][:, 0]
        self.object['leading_photon'] = self.object['photon'][ak.argmax(self.object['photon'].pt, axis=1, keepdims=True)][:, 0]
        self.object['photon-jet'] = self.object['heaviest_jet'] + self.object['leading_photon']
        delta_phi = abs(self.object['heaviest_jet'].phi - self.object['leading_photon'].phi)
        delta_phi = ak.min([delta_phi, 2 * np.pi - delta_phi], axis=0)
        final_cut = self.pass_cut(cutName='photon-jet_delta_phi', cut=(delta_phi>2.9), final=True)
        """
        
        self.store_variables(vars={
            'AK8jet': {'pt', 'eta', 'phi', 'mass', 'msoftdrop'},
            'photon': {'pt', 'eta', 'phi', 'mass'},
            'b-jet': {'pt', 'eta', 'phi', 'mass'},
            'event': {'MET_pt', 'genWeight', 'event'},
            'photon-jet': {'pt', 'eta', 'phi', 'mass'},
        })
        
        # Additional vars by special computing
        # self.variables['photon-jet_deltaR'] = self.object['photon'].delta_r(self.object['AK8jet'])
        self.variables['photon-jet_deltaR'] = self.object['AK8jet'].delta_r(self.object['photon'])

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
        if any(final_cut) and self.sample_type == 'mc':
            gen = GenMatch(events=self.event)
            if self.channel == 'ZpToHGamma':
                self.variables.update(gen.ZpToHGamma())
            elif self.channel in ['QCD', 'TTJets', 'WJetsToQQ', 'ZJetsToQQ']:
                final_cut = self.pass_cut(cutName='no prompt photon', cut=gen.all_fake_photon(), final=True)
            elif self.channel in ['GJets', 'TTGJets', 'WGamma', 'ZGamma']:
                final_cut = self.pass_cut(cutName='any prompt photon', cut=gen.any_true_photon(), final=True)
                
        # store output
        if any(final_cut):
            self.to_parquet(array=ak.Array(self.variables))
        return {self.mode: {k: float(v) for (k, v) in self.cutflow.items()}}
    
    @property  # transform method into attribute and make it unchangable to hide _accumulator
    def accumulator(self):
        return self._accumulator
    
    def postprocess(self, accumulator):
        return super().postprocess(accumulator)
