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
    def __init__(
        self, triggers: list=['Photon175', 'Photon165_R9Id90_HE10_IsoM'],
        cut: dict={
            'deltaR': {'min': 0},
        }, outdir: str=os.path.join('..', 'output')
    ) -> None:
        super().__init__()
        self.triggers = triggers
        self.object = {'event': None}
        self.variables = {}
        self.cutflow = {}
        self.outdir = outdir
        self.cut = cut
        self.oldCut = None
        self.newCut = 'raw'
        self._accumulator = processor.defaultdict_accumulator()
    
    def _passCut(self, cutName: str, cut: ak.Array):
        self.oldCut, self.newCut = self.newCut, cutName
        self.cutflow[self.newCut] = ak.fill_none(array=self.cutflow[self.oldCut] * cut, value=False)
        self.object['event'] = ak.mask(array=self.object['event'], mask=self.cutflow[self.newCut])
    
    def _to_parquet(self, arrays: dict) -> None:
        output_dir = os.path.abspath(self.outdir)
        tokens = self.object['event'].behavior["__events_factory__"]._partition_key.split('/')
        output_file = '_'.join([(t if 'Events' not in t else 'Events') for t in tokens ])
        ak.to_parquet(arrays, where = output_dir+'/'+output_file+'.parq')
        

    def __preselect_HGamma( ## _ in prefix means private method
            self, events: NanoEventsArray,
            variables: dict={
                'AK8jet': {'pt', 'eta', 'phi', 'mass', 'msoftdrop'},
                'photon': {'pt', 'eta', 'phi', 'mass'},
                'event': {'MET_pt'},
                'photon-jet': {'pt', 'eta', 'phi', 'mass'},
            }
        ) -> ak.Array:
        ## initialize
        self.object['event'] = events
        self.cutflow['raw'] = ak.Array([True for _ in range(len(self.object['event']))])
        
        ## triggers
        self._passCut(
            cut=ak.sum(
                [self.object['event'].HLT[trigger] for trigger in self.triggers if trigger in self.object['event'].HLT.fields], axis=0
            ) > 0, cutName='trigger'
        ) ## pass any trigger
        
        ## b veto
        raw_AK4jet = self.object['event'].Jet
        b_tagging = (raw_AK4jet.btagDeepFlavB > 0.2783) # Working Points -- loose: 0.0490, medium: 0.2783, tight: 0.7100
        # refer to https://gitlab.cern.ch/groups/cms-btv/-/wikis/SFCampaigns/UL2018
        self._passCut(cutName='b-veto', cut=(ak.sum(b_tagging, axis=-1)==0)) ## b-veto

        ## Muon
        raw_muon = self.object['event'].Muon # (event, muon)
        muon_cut = ( # (event, boolean)
            # high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)
            (raw_muon.highPtId == 2) & 
            (raw_muon.tkRelIso < 0.1) & # Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt
            (abs(raw_muon.eta) < 2.4) & 
            (raw_muon.pt > 20) # I don't use `muon_corrected_pt` coming from ROOT.RoccoR
        )
        self._passCut(cutName='muon', cut=(ak.sum(muon_cut, axis=-1)==0)) ## 0 muon

        ## Electron
        raw_electron = self.object['event'].Electron # (event, electron)
        electron_cut = ( # (event, boolean)
            (raw_electron.cutBased_HEEP == True) & # cut-based HEEP ID
            (abs(raw_electron.eta) < 2.5) &
            (raw_electron.pt > 20)
        )
        self._passCut(cutName='electron', cut=(ak.sum(electron_cut, axis=-1)==0)) ## 0 electron
        
        ## Photon
        raw_photon = self.object['event'].Photon # (event, photon), >=1 photon per event
        photon_cut = ( # (event, boolean)
            (raw_photon.mvaID_WP90 > 0.2) &
            (raw_photon.pt > 200) &
            (abs(raw_photon.eta) < 2.4)
        )
        self.object['photon'] = raw_photon[photon_cut]
        self._passCut(cutName='photon', cut=(ak.sum(photon_cut, axis=-1)>0)) ## >=1 photon

        ## AK8 jet
        raw_AK8jet = self.object['event'].FatJet # (event, fatjet), >=1 AK8 jet per event
        AK8jet_cut = ( # (event, boolean)
            (raw_AK8jet.msoftdrop > 30) & # Corrected soft drop mass with PUPPI
            (raw_AK8jet.pt > 250) & 
            (abs(raw_AK8jet.eta) < 2.4) & 
            (raw_AK8jet.jetId&2 > 0)
            # Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto
        )
        self.object['AK8jet'] = raw_AK8jet[AK8jet_cut]
        self._passCut(cutName='AK8jet', cut=(ak.sum(AK8jet_cut, axis=-1)>0)) ## >=1 AK8 jet
        
        ## Photon-Jet cleaning
        pj_pair = ak.cartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_index_pair = ak.argcartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_dr = pj_pair.photon.delta_r(pj_pair.jet)
        pj_clean = (pj_dr > self.cut['deltaR']['min'])
        photon_index, jet_index = pj_index_pair.photon[pj_clean], pj_index_pair.jet[pj_clean]
        self._passCut(cutName='photon-jet_cleaning', cut=(ak.sum(pj_clean, axis=-1)==1)) ## exactly 1 pair photon-jet
        
        ## final event-cut
        final_cut = self.cutflow['photon-jet_cleaning']
        self.object['event'] = self.object['event'][final_cut]
        self.object['photon'] = self.object['photon'][photon_index][final_cut] ## shape=(event, photon), 1 photon per event
        self.object['AK8jet'] = self.object['AK8jet'][jet_index][final_cut] ## shape=(event, AK8jet), 1 AK8jet per event
        self.object['photon-jet'] = self.object['photon'][:, 0] + self.object['AK8jet'][:, 0] ## shape=(event, )
        
        ## Return vars of objects after pre-selection
        for obj in variables.keys():
            if obj=='event':
                self.variables.update({
                    obj+'_'+var: self.object[obj][var.split('_')[0]]['_'.join(var.split('_')[1:])] for var in variables[obj]
                })
            else:
                self.variables.update({ obj+'_'+var: getattr(self.object[obj], var) for var in variables[obj] })

        ## Additional vars by specific computing
        self.variables['photon-jet_deltaR'] = ak.flatten(self.object['photon'].delta_r(self.object['AK8jet']), axis=-1)
        
        return final_cut

    def process(self, events: NanoEventsArray):
        event_cut = self.__preselect_HGamma(events=events)
        if all(event_cut==False):
            return {k: ak.sum(v) for (k,v) in self.cutflow.items()}
        events = events[event_cut]
        gen_match = GenMatch()
        self.variables.update(gen_match.HGamma(events))
        self._to_parquet(arrays=self.variables)
        
        return {k: ak.sum(v) for (k,v) in self.cutflow.items()}
    
    @property ## transform method into attribute and make it unchangable to hide _accumulator
    def accumulator(self):
        return self._accumulator
    
    def postprocess(self, accumulator):
        return super().postprocess(accumulator)