import numpy as np
import awkward as ak
import uproot

from coffea import processor
from coffea.nanoevents.methods import candidate
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods.base import NanoEventsArray
from coffea.nanoevents.methods.nanoaod import FatJetArray, GenParticleArray
from .GenMatch import GenMatch

class Processor(processor.ProcessorABC):
    def __init__(self) -> None:
        super().__init__()
        self.object = {}
        self.variables = {}
        self._accumulator = processor.defaultdict_accumulator() ## useless
    
    def __preselect_HGamma( ## _ in prefix means private method
            self, events: NanoEventsArray, deltaR_cut: float=0.8,
            variables: dict={
                'AK8jet': {'pt', 'eta', 'phi', 'mass', 'msoftdrop'},
                'photon': {'pt', 'eta', 'phi', 'mass', 'msoftdrop'},
                'event': {'MET_pt'},
            }, trigger: str='Photon175'
        ):
        self.object['event'] = events
        event_cut = events.HLT[trigger]
        
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
        event_cut *= (ak.sum(muon_cut, axis=-1)==0)

        ## Electron
        electrons = self.object['event'].Electron # (event, electron)
        electron_cut = ( # (event, boolean)
            (electrons.cutBased_HEEP == True) & # cut-based HEEP ID
            (abs(electrons.eta) < 2.5) &
            (electrons.pt > 20)
        )
        self.object['electron'] = electrons[electron_cut]
        event_cut *= (ak.sum(muon_cut, axis=-1)==0)
        
        ## Photon
        photons = self.object['event'].Photon # (event, photon)
        photon_cut = ( # (event, boolean)
            (photons.mvaID_WP90 > 0.2) &
            (photons.pt > 200) &
            (abs(photons.eta) < 2.4)
        )
        self.object['photon'] = photons[photon_cut]
        event_cut *= (ak.sum(photon_cut, axis=-1)>0)

        ## AK8 jet
        AK8jets = self.object['event'].FatJet # (event, fatjet)
        AK8jet_cut = ( # (event, boolean)
            (AK8jets.msoftdrop > 30) & # Corrected soft drop mass with PUPPI
            (AK8jets.pt > 250) & 
            (abs(AK8jets.eta) < 2.6) & 
            (AK8jets.jetId&2 > 0)
            # Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto
        )
        self.object['AK8jet'] = AK8jets[AK8jet_cut]
        event_cut *= (ak.sum(AK8jet_cut, axis=-1)>0)
        
        ## Photon-Jet cleaning
        photon_jet_pair = ak.cartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        photon_jet_pair_index = ak.argcartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        
        photon_jet_dr = photon_jet_pair.photon.delta_r(photon_jet_pair.jet)
        photon_jet_clean = (
            ak.num(photon_jet_dr > deltaR_cut, axis=1)==1
        )
        event_cut *= (ak.sum(photon_jet_clean, axis=-1)==1)
        
        photon_index, jet_index = photon_jet_pair_index[photon_jet_clean].photon, photon_jet_pair_index[photon_jet_clean].jet
        self.object['photon'] = self.object['photon'][photon_index]
        self.object['AK8jet'] = self.object['AK8jet'][jet_index]
        
        ## Return vars of objects after pre-selection
        for obj in ['event', 'photon', 'AK8jet']:
            self.object[obj] = self.object[obj][event_cut]
        
        self.variables.update({
            var: (
                self.object[obj][var] if obj!='event' else self.object[obj][var.split('_')[0]][var.split('_')[1:]]
            ) for (obj, var) in variables.items()
        })
        
        return event_cut

    def process(self, events: NanoEventsArray):
        event_cut = self.__preselect_HGamma(events)
        events = events[event_cut]
        gen_match = GenMatch()
        self.variables.update(gen_match.HGamma(events))
        
        return self.variables
    
    @property ## transform method into attribute and make it unchangable to hide _accumulator
    def accumulator(self):
        return self._accumulator
    
    def postprocess(self, accumulator):
        return super().postprocess(accumulator)