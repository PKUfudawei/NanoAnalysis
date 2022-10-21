import numpy as np
import awkward as ak
import os

from coffea import processor
from coffea.nanoevents.methods import candidate
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods.base import NanoEventsArray
from coffea.nanoevents.methods.nanoaod import FatJetArray, GenParticleArray
from .GenMatch import GenMatch

class Processor(processor.ProcessorABC):
    def __init__(
        self, machine: str, outdir: str, channel: str,
        cutValue: dict={
            'deltaR': {'min': 1.1},
        }, triggers: list=['Photon175', 'Photon165_R9Id90_HE10_IsoM'],
    ) -> None:
        super().__init__()
        self.triggers = triggers
        self.event = None
        self.tag = {}
        self.object = {}
        self.variables = {}
        self.cutflow = {}
        self.outdir = os.path.abspath(outdir)
        if machine not in ['local', 'condor']:
            raise ValueError("Processor.__init__(): machine must be in ['local', 'condor']")
        self.machine = machine
        self.channel = channel
        self.cutValue = cutValue
        self.cuts = PackedSelection()
        self._accumulator = processor.defaultdict_accumulator()
    
    
    def pass_cut(self, cutName: str, cut: ak.Array, mask: bool=True) -> None:
        self.cuts.add(cutName, cut)
        self.cutflow[cutName] = self.cuts.all(*self.cuts.names)
        if mask: ## update all objects after passing cut
            self.event = ak.mask(array=self.event, mask=self.cutflow[cutName])
            for obj in self.object: ## keep size: (event, object) in intermediate process
                self.object[obj] = ak.mask(self.object[obj], mask=self.cutflow[cutName])
        else: ## make projection after all cuts to reduce event size
            self.event = self.event[self.cutflow[cutName]]
            for obj in self.object: 
                self.object[obj] = ak.flatten(self.object[obj][self.cutflow[cutName]], axis=1)
    
    
    def triggered(self, level: str='any') -> ak.Array:
        if level not in ['any', 'all']:
            raise ValueError("Processor.passTriggers(): level must be in ['any', 'all']")
        elif level == 'any':
            return ak.sum([self.event.HLT[t] for t in self.triggers if t in self.event.HLT.fields], axis=0) > 0 
        elif level == 'all':
            raise ValueError("Processor.passTriggers(level='all') not finished yet")## not finished yet
            ## pass all triggers
            
        
    def b_tag(self, level: str='tight', reco: bool=False) -> ak.Array:
        if level not in ['loose', 'medium', 'tight']:
            raise ValueError("Processor.b_veto(): level must be in ['loose', 'medium', 'tight']")
        
        raw_AK4jet = self.event.Jet
        # Working Points -- loose: 0.0490, medium: 0.2783, tight: 0.7100
        # refer to https://gitlab.cern.ch/groups/cms-btv/-/wikis/SFCampaigns/UL2018
        WP = {'loose': 0.0490, 'medium': 0.2783, 'tight': 0.7100} 
        self.tag['b'] = (raw_AK4jet.btagDeepFlavB > WP[level]) 
        if reco:
            self.object['b'] = self.event.Jet[self.tag['b']]
        return self.tag['b']
    
    
    def muon_tag(self, reco: bool=False) -> ak.Array:
        raw_muon = self.event.Muon # (event, muon)
        self.tag['muon'] = ( # (event, boolean)
            # high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)
            (raw_muon.highPtId == 2) & 
            (raw_muon.tkRelIso < 0.1) & # Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt
            (abs(raw_muon.eta) < 2.4) & 
            (raw_muon.pt > 20) # I don't use `muon_corrected_pt` coming from ROOT.RoccoR
        )
        if reco:
            self.object['muon'] = self.event.Muon[self.tag['muon']]
        return self.tag['muon']
    
    
    def electron_tag(self, reco: bool=False) -> ak.Array:
        raw_electron = self.event.Electron # (event, electron)
        self.tag['electron'] = ( # (event, boolean)
            (raw_electron.cutBased_HEEP == True) & # cut-based HEEP ID
            (abs(raw_electron.eta) < 2.5) &
            (raw_electron.pt > 20)
        )
        if reco:
            self.object['electron'] = self.event.Electron[self.tag['electron']]
        return self.tag['electron']
    
    
    def photon_tag(self, reco: bool=False) -> ak.Array:
        raw_photon = self.event.Photon # (event, photon), >=1 photon per event
        self.tag['photon'] = ( # (event, boolean)
            (raw_photon.mvaID_WP90 > 0.2) &
            (raw_photon.pt > 200) &
            (abs(raw_photon.eta) < 2.4)
        )
        if reco:
            self.object['photon'] = self.event.Photon[self.tag['photon']]
        return self.tag['photon']
    
    
    def AK8jet_tag(self, reco: bool=False) -> ak.Array:
        raw_AK8jet = self.event.FatJet # (event, fatjet), >=1 AK8 jet per event
        self.tag['AK8jet'] = ( # (event, boolean)
            (raw_AK8jet.msoftdrop > 40) & # Corrected soft drop mass with PUPPI
            (raw_AK8jet.pt > 250) & 
            (abs(raw_AK8jet.eta) < 2.4) & 
            (raw_AK8jet.jetId&2 > 0)
            # Jet ID flags bit1 is loose (always false in 2017 since it does not exist), 
            # bit2 is tight, bit3 is tightLepVeto
        )
        if reco:
            self.object['AK8jet'] = self.event.FatJet[self.tag['AK8jet']]
        return self.tag['AK8jet']
    
    
    def store_variables(self, vars: dict):
        for obj in vars.keys():
            if obj!='event':
                self.variables.update({obj+'_'+var: getattr(self.object[obj], var) for var in vars[obj]})
            else:
                self.variables.update({
                    'event_'+var: self.event[var.split('_')[0]]['_'.join(var.split('_')[1:])] 
                    for var in vars['event']
                })
    
    
    def to_parquet(self, array: ak.Array) -> None:
        if self.machine not in ['local', 'condor']:
            raise ValueError("Processor.__init__(): machine must be in ['local', 'condor']")

        output_dir = os.path.abspath(self.outdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if self.machine=='local':
            tokens = self.event.behavior["__events_factory__"]._partition_key.split('/')
            name = '_'.join([(t if 'Events' not in t else 'Events') for t in tokens ])
        elif self.machine=='condor':
            name = 'output'
        
        ak.to_parquet(array=array, where=os.path.join(output_dir, f'{name}.parq'))
    
      
    def preselect_HGamma(self) -> ak.Array: ## __ in prefix means private method
        ## at least pass one trigger
        self.pass_cut(cutName='triggered', cut=self.triggered(level='any'))
        
        ## b veto
        self.pass_cut(cutName='b-veto', cut=(ak.sum(self.b_tag(level='tight'), axis=1)==0))

        ## Muon veto
        self.pass_cut(cutName='muon-veto', cut=(ak.sum(self.muon_tag(), axis=1)==0))
        
        ## Electron veto
        self.pass_cut(cutName='electron-veto', cut=(ak.sum(self.electron_tag(), axis=1)==0)) 
        
        ## Photon >=1
        self.pass_cut(cutName='photon', cut=(ak.sum(self.photon_tag(reco=True), axis=1)>0)) 
        
        ## AK8 jet >=1
        self.pass_cut(cutName='AK8jet', cut=(ak.sum(self.AK8jet_tag(reco=True), axis=1)>0))
        
        ## Photon-Jet cleaning, a very special part so no function definition here
        pj_pair = ak.cartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_index_pair = ak.argcartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_dr = pj_pair.photon.delta_r(pj_pair.jet)
        pj_clean = (pj_dr > self.cutValue['deltaR']['min'])
        photon_index, jet_index = pj_index_pair.photon[pj_clean], pj_index_pair.jet[pj_clean]
        ## exactly 1 pair of photon and jet passed jet-cleaning requirement
        self.object['photon'] = self.object['photon'][photon_index]
        self.object['AK8jet'] = self.object['AK8jet'][jet_index]
        self.object['photon-jet'] = self.object['photon'] + self.object['AK8jet']
        self.pass_cut(cutName='photon-jet_cleaning', cut=(ak.sum(pj_clean, axis=-1)==1), mask=False) 
        
        ## Return vars of objects after pre-selection
        self.store_variables(vars={
            'AK8jet': {'pt', 'eta', 'phi', 'mass', 'msoftdrop'},
            'photon': {'pt', 'eta', 'phi', 'mass'},
            'event': {'MET_pt'},
            'photon-jet': {'pt', 'eta', 'phi', 'mass'},
        })

        ## Additional vars by special computing
        self.variables['photon-jet_deltaR'] = self.object['photon'].delta_r(self.object['AK8jet'])

        return self.cuts.all(*self.cuts.names)
    

    def process(self, events: NanoEventsArray) -> dict:
        ## initialize
        self.event = events
        
        ## process
        event_cut = self.preselect_HGamma()
        cutflow = {k: int(ak.sum(v)) for (k,v) in self.cutflow.items()}
        if all(event_cut==False):
            self.to_parquet(array=ak.Array({}))
            return cutflow
        
        ## gen-macthing
        if self.channel == 'ZpToHGamma':
            gen_match = GenMatch()
            self.variables.update(gen_match.ZpToHGamma(self.event))
        
        ## store output
        self.to_parquet(array=ak.Array(self.variables))
        return cutflow
    
    
    @property ## transform method into attribute and make it unchangable to hide _accumulator
    def accumulator(self):
        return self._accumulator
    
    
    def postprocess(self, accumulator):
        return super().postprocess(accumulator)