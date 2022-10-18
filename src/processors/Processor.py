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
        cut: dict={
            'deltaR': {'min': 1.1},
        }, triggers: list=['Photon175', 'Photon165_R9Id90_HE10_IsoM'],
    ) -> None:
        super().__init__()
        self.triggers = triggers
        self.object = {'event': None}
        self.variables = {}
        self.cutflow = {'raw': None}
        self.outdir = os.path.abspath(outdir)
        if machine not in ['local', 'condor']:
            raise ValueError("Processor.__init__(): machine must be in ['local', 'condor']")
        self.machine = machine
        self.channel = channel
        self.cut = cut
        self.prevCut = None
        self.nextCut = 'raw'
        self._accumulator = processor.defaultdict_accumulator()
    
    def passCut(self, cutName: str, cut: ak.Array) -> None:
        self.prevCut, self.nextCut = self.nextCut, cutName
        self.cutflow[self.nextCut] = ak.fill_none(array=self.cutflow[self.prevCut] * cut, value=False)
        self.object['event'] = ak.mask(array=self.object['event'], mask=self.cutflow[self.nextCut])
        
    def passTriggers(self, level: str='any') -> None:
        if level not in ['any', 'all']:
            raise ValueError("Processor.passTriggers(): level must be in ['any', 'all']")
        elif level == 'any':
            self.passCut(cutName='trigger',
                cut=ak.sum([self.object['event'].HLT[t] for t in self.triggers if t in self.object['event'].HLT.fields], axis=0) > 0 
            ) ## pass any trigger
        elif level == 'all':
            return  ValueError("Processor.passTriggers(level='all') not finished yet")## not finished yet
            ## pass all triggers
            
        
    def b_veto(self, level: str='tight') -> None:
        if level not in ['loose', 'medium', 'tight']:
            raise ValueError("Processor.b_veto(): level must be in ['loose', 'medium', 'tight']")
        
        raw_AK4jet = self.object['event'].Jet
        WP = {'loose': 0.0490, 'medium': 0.2783, 'tight': 0.7100} # Working Points -- loose: 0.0490, medium: 0.2783, tight: 0.7100
        # refer to https://gitlab.cern.ch/groups/cms-btv/-/wikis/SFCampaigns/UL2018
        b_tagging = (raw_AK4jet.btagDeepFlavB > WP[level]) 
        self.passCut(cutName='b-veto', cut=(ak.sum(b_tagging, axis=-1)==0)) ## b-veto
    
    def to_parquet(self, array: ak.Array) -> None:
        if self.machine not in ['local', 'condor']:
            raise ValueError("Processor.__init__(): machine must be in ['local', 'condor']")

        output_dir = os.path.abspath(self.outdir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if self.machine=='local':
            tokens = self.object['event'].behavior["__events_factory__"]._partition_key.split('/')
            name = '_'.join([(t if 'Events' not in t else 'Events') for t in tokens ])
        elif self.machine=='condor':
            name = 'output'
        
        ak.to_parquet(array=array, where=os.path.join(output_dir, f'{name}.parq'))
        
    def __preselect_HGamma(self, ## __ in prefix means private method
            obj_vars: dict={
                'AK8jet': {'pt', 'eta', 'phi', 'mass', 'msoftdrop'},
                'photon': {'pt', 'eta', 'phi', 'mass'},
                'event': {'MET_pt'},
                'photon-jet': {'pt', 'eta', 'phi', 'mass'},
            }
        ) -> None:
        ## triggers
        self.passTriggers(level='any') ## at least pass one trigger
        
        ## b veto
        self.b_veto(level='tight')

        ## Muon
        raw_muon = self.object['event'].Muon # (event, muon)
        muon_cut = ( # (event, boolean)
            # high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)
            (raw_muon.highPtId == 2) & 
            (raw_muon.tkRelIso < 0.1) & # Tracker-based relative isolation dR=0.3 for highPt, trkIso/tunePpt
            (abs(raw_muon.eta) < 2.4) & 
            (raw_muon.pt > 20) # I don't use `muon_corrected_pt` coming from ROOT.RoccoR
        )
        ## 0 muon requirement satisfied
        self.passCut(cutName='muon', cut=(ak.sum(muon_cut, axis=-1)==0))
        

        ## Electron
        raw_electron = self.object['event'].Electron # (event, electron)
        electron_cut = ( # (event, boolean)
            (raw_electron.cutBased_HEEP == True) & # cut-based HEEP ID
            (abs(raw_electron.eta) < 2.5) &
            (raw_electron.pt > 20)
        )
        ## 0 electron requirement satisfied
        self.passCut(cutName='electron', cut=(ak.sum(electron_cut, axis=-1)==0)) 
        
        
        ## Photon
        raw_photon = self.object['event'].Photon # (event, photon), >=1 photon per event
        photon_cut = ( # (event, boolean)
            (raw_photon.mvaID_WP90 > 0.2) &
            (raw_photon.pt > 200) &
            (abs(raw_photon.eta) < 2.4)
        )
        ## the photons from events matching all above cuts
        self.object['photon'] = raw_photon[photon_cut]
        ## >=1 photon requirement satisfied
        self.passCut(cutName='photon', cut=(ak.num(self.object['photon'], axis=1)>0)) 
        
        ## AK8 jet
        raw_AK8jet = self.object['event'].FatJet # (event, fatjet), >=1 AK8 jet per event
        AK8jet_cut = ( # (event, boolean)
            (raw_AK8jet.msoftdrop > 30) & # Corrected soft drop mass with PUPPI
            (raw_AK8jet.pt > 250) & 
            (abs(raw_AK8jet.eta) < 2.4) & 
            (raw_AK8jet.jetId&2 > 0)
            # Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto
        )
        ## the AK8jets from events matching all above cuts
        self.object['AK8jet'] = raw_AK8jet[AK8jet_cut]
        ## >=1 AK8 jet requirement satisfied
        self.passCut(cutName='AK8jet', cut=(ak.num(self.object['AK8jet'], axis=1)>0)) 
        
        ## Photon-Jet cleaning
        pj_pair = ak.cartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_index_pair = ak.argcartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_dr = pj_pair.photon.delta_r(pj_pair.jet)
        pj_clean = (pj_dr > self.cut['deltaR']['min'])
        photon_index, jet_index = pj_index_pair.photon[pj_clean], pj_index_pair.jet[pj_clean]
        ## exactly 1 pair of photon-jet passed jet-cleaning requirement
        self.passCut(cutName='photon-jet_cleaning', cut=(ak.sum(pj_clean, axis=-1)==1)) 
        
        ## final event-cut
        final_cut = self.cutflow[self.nextCut]
        ## drop unwanted objects by projection
        self.object['event'] = self.object['event'][final_cut]
        self.object['photon'] = self.object['photon'][photon_index][final_cut] ## shape=(event, photon), 1 photon per event
        self.object['AK8jet'] = self.object['AK8jet'][jet_index][final_cut] ## shape=(event, AK8jet), 1 AK8jet per event
        self.object['photon-jet'] = self.object['photon'][:, 0] + self.object['AK8jet'][:, 0] ## shape=(event, )
        
        ## Return vars of objects after pre-selection
        for obj in obj_vars.keys():
            if obj=='event':
                self.variables.update({
                    obj+'_'+var: self.object[obj][var.split('_')[0]]['_'.join(var.split('_')[1:])] for var in obj_vars[obj]
                })
            else:
                self.variables.update({ obj+'_'+var: getattr(self.object[obj], var) for var in obj_vars[obj] })

        ## Additional vars by specific computing
        self.variables['photon-jet_deltaR'] = ak.flatten(self.object['photon'].delta_r(self.object['AK8jet']), axis=-1)

        return final_cut

    def process(self, events: NanoEventsArray) -> dict:
        ## preprocessing
        self.object['event'] = events
        self.cutflow['raw'] = ak.Array([True for _ in range(len(self.object['event']))])
        
        ## processing
        event_cut = self.__preselect_HGamma()
        cutflow = {k: int(ak.sum(v)) for (k,v) in self.cutflow.items()}
        if all(event_cut==False):
            self.to_parquet(arrays={})
            return cutflow
        
        ## gen-macthing
        if self.channel == 'ZpToHGamma':
            gen_match = GenMatch()
            self.variables.update(gen_match.ZpToHGamma(self.object['event']))
        
        ## storing output
        self.to_parquet(array=ak.Array(self.variables))
        return cutflow
    
    @property ## transform method into attribute and make it unchangable to hide _accumulator
    def accumulator(self):
        return self._accumulator
    
    def postprocess(self, accumulator):
        return super().postprocess(accumulator)