import numpy as np
import awkward as ak
import uproot
import os

from coffea import processor
from coffea.nanoevents.methods import candidate
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods.base import NanoEventsArray
from coffea.nanoevents.methods.nanoaod import FatJetArray, GenParticleArray
from .GenMatch import GenMatch

class Processor(processor.ProcessorABC):
    def __init__(self, triggers: list=[]) -> None:
        super().__init__()
        self.triggers = triggers
        self.object = {}
        self.variables = {}
        self._accumulator = processor.defaultdict_accumulator() ## useless
    
    def __to_parquet(self, arrays: dict) -> None:
        output_dir = os.path.abspath(os.path.join('..', 'output'))
        tokens = self.object['event'].behavior["__events_factory__"]._partition_key.split('/')
        output_file = '_'.join([(t if 'Events' not in t else 'Events') for t in tokens ])
        ak.to_parquet(arrays, where = output_dir+'/'+output_file)


    def __preselect_HGamma( ## _ in prefix means private method
            self, events: NanoEventsArray, deltaR_cut: float=0.8,
            variables: dict={
                'AK8jet': {'pt', 'eta', 'phi', 'mass', 'msoftdrop'},
                'photon': {'pt', 'eta', 'phi', 'mass'},
                'event': {'MET_pt'},
            }, triggers: list=['Photon175', 'Photon165_R9Id90_HE10_IsoM']
        ) -> ak.Array:
        self.object['event'] = events
        if triggers:
            self.triggers = triggers
        event_cut = ak.sum([events.HLT[trigger] for trigger in self.triggers if trigger in events.HLT.fields], axis=0)>0
        
        ## Muon
        muons = self.object['event'].Muon # (event, muon)
        muon_cut = ( # (event, boolean)
            # high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)
            (muons.highPtId == 2) & 
            (muons.tkRelIso < 0.1) & # Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt
            (abs(muons.eta) < 2.4) & 
            (muons.pt > 20) # I don't use `muon_corrected_pt` coming from ROOT.RoccoR
        )
        self.object['muon'] = muons[muon_cut]
        event_cut = event_cut * (ak.sum(muon_cut, axis=-1)==0) ## 0 muon

        ## Electron
        electrons = self.object['event'].Electron # (event, electron)
        electron_cut = ( # (event, boolean)
            (electrons.cutBased_HEEP == True) & # cut-based HEEP ID
            (abs(electrons.eta) < 2.5) &
            (electrons.pt > 20)
        )
        self.object['electron'] = electrons[electron_cut]
        event_cut = event_cut * (ak.sum(muon_cut, axis=-1)==0) ## 0 electron
        
        ## Photon
        photons = self.object['event'].Photon # (event, photon)
        photon_cut = ( # (event, boolean)
            (photons.mvaID_WP90 > 0.2) &
            (photons.pt > 200) &
            (abs(photons.eta) < 2.4)
        )
        self.object['photon'] = photons[photon_cut]
        event_cut = event_cut * (ak.sum(photon_cut, axis=-1)>0) ## >=1 photon

        ## AK8 jet
        AK8jets = self.object['event'].FatJet # (event, fatjet)
        AK8jet_cut = ( # (event, boolean)
            (AK8jets.msoftdrop > 30) & # Corrected soft drop mass with PUPPI
            (AK8jets.pt > 250) & 
            (abs(AK8jets.eta) < 2.4) & 
            (AK8jets.jetId&2 > 0)
            # Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto
        )
        self.object['AK8jet'] = AK8jets[AK8jet_cut]
        event_cut = event_cut * (ak.sum(AK8jet_cut, axis=-1)>0) ## >=1 AK8 jet
        
        ## Photon-Jet cleaning
        photon_jet_pair = ak.cartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        photon_jet_index_pair = ak.argcartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        
        photon_jet_dr = photon_jet_pair.photon.delta_r(photon_jet_pair.jet)
        photon_jet_clean = (photon_jet_dr > deltaR_cut)
        photon_index, jet_index = photon_jet_index_pair.photon[photon_jet_clean], photon_jet_index_pair.jet[photon_jet_clean]
        
        self.object['photon'] = self.object['photon'][photon_index] ## may exists repetition
        self.object['AK8jet'] = self.object['AK8jet'][jet_index] ## may exists repetition
        event_cut = event_cut * (ak.sum(photon_jet_clean, axis=-1)==1) ## exactly 1 pair
        
        ## Return vars of objects after pre-selection
        for obj in ['event', 'photon', 'AK8jet']:
            self.object[obj] = self.object[obj][event_cut]
            if obj=='event':
                self.variables.update({
                    obj+'_'+var: self.object[obj][var.split('_')[0]]['_'.join(var.split('_')[1:])] for var in variables[obj]
                })
            else:
                self.variables.update({obj+'_'+var: self.object[obj][var] for var in variables[obj]})

        self.variables['deltaR_jet_phton'] = ak.flatten(self.object['photon'].delta_r(self.object['AK8jet']), axis=-1)

        return event_cut

    def process(self, events: NanoEventsArray):
        event_cut = self.__preselect_HGamma(events=events)
        events = events[event_cut]
        gen_match = GenMatch()
        self.variables.update(gen_match.HGamma(events))
        self.__to_parquet(arrays=self.variables)
        
        return {}
    
    @property ## transform method into attribute and make it unchangable to hide _accumulator
    def accumulator(self):
        return self._accumulator
    
    def postprocess(self, accumulator):
        return super().postprocess(accumulator)