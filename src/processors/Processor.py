import awkward as ak
import numpy as np
import os

from coffea import processor, lumi_tools
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods.base import NanoEventsArray
from .GenMatch import GenMatch


class Processor(processor.ProcessorABC):
    def __init__(
        self, outdir: str, mode: str,
        triggers: dict = {
            '2016': ['Photon175', 'Photon165_R9Id90_HE10_IsoM'],
            '2017': ['HLT_Photon200'],
            '2018': ['HLT_Photon200'],
        }
    ) -> None:
        super().__init__()
        self.triggers = triggers
        self.event = None
        self.tag = {}
        self.object = {}
        self.variables = {}
        self.cutflow = {}
        self.outdir = os.path.abspath(outdir)
        if mode.split('_')[0] not in ['data', 'mc']:
            raise ValueError("Processor.__init__(): mode must start with 'data' or 'mc'")
        self.mode = mode  # = '$type_$year_$channel'
        self.sample_type = self.mode.split('_')[0]
        self.year = self.mode.split('_')[1]
        self.channel = self.mode.split('_')[2]
        self.cuts = PackedSelection()
        self._accumulator = processor.defaultdict_accumulator()
    
    def lumi_mask(self, events: NanoEventsArray) -> NanoEventsArray:  # only applied on data
        golden_JSON = {
            '2018': '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt',
            '2017': '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Legacy_2017/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt',
            '2016': '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt',
        }
        if self.sample_type=='mc':  # usually skipped cuz type is restricted to data before executing this function
            return events
        elif self.sample_type=='data':
            lumi_mask = lumi_tools.LumiMask(golden_JSON[self.year])
            data_mask = lumi_mask(events.run, events.luminosityBlock)
            return events[data_mask]
        
    def pass_cut(self, cutName: str, cut: ak.Array, final: bool = False) -> ak.Array:
        if self.cuts.names and len(cut) != len(self.cuts.all()):
            self.cuts = PackedSelection()
        self.cuts.add(cutName, cut)
        
        # calculate cutflow
        event_flag = self.cuts.all(*self.cuts.names)
        if self.mode.startswith('data'):
            self.cutflow[cutName] = ak.sum(event_flag)
        elif self.mode.startswith('mc'):
            self.cutflow[cutName] = ak.sum(self.event.genWeight[event_flag])
        
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
            raise ValueError("Processor.passTriggers(): level must be in ['any', 'all']")
        elif level == 'any':  # pass any trigger
            return ak.any([self.event.HLT[t] for t in self.triggers[self.year] if t in self.event.HLT.fields], axis=0)
        elif level == 'all':  # pass all triggers
            return ak.all([self.event.HLT[t] for t in self.triggers[self.year] if t in self.event.HLT.fields], axis=0)
            
    def b_tag(self, level: str = 'tight', reco: bool = False) -> ak.Array:
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
            (raw_photon.mvaID_WP90 > 0.2) &
            (raw_photon.pt > 200) &
            ((abs(raw_photon.eta) < 1.4442) | ((abs(raw_photon.eta) > 1.566) & (abs(raw_photon.eta) < 2.5)))
        )
        if reco:
            self.object['photon'] = self.event.Photon[self.tag['photon']]
        return self.tag['photon']
    
    def AK8jet_tag(self, reco: bool = False) -> ak.Array:
        raw_AK8jet = self.event.FatJet  # (event, fatjet), >=1 AK8 jet per event
        self.tag['AK8jet'] = (  # (event, boolean)
            (raw_AK8jet.msoftdrop > 40) &  # Corrected soft drop mass with PUPPI
            (raw_AK8jet.pt > 250) &
            (abs(raw_AK8jet.eta) < 2.4) &
            (raw_AK8jet.jetId & 2 > 0)
            # Jet ID flags bit1 is loose (always false in 2017 since it does not exist),
            # bit2 is tight, bit3 is tightLepVeto
        )
        if reco:
            self.object['AK8jet'] = self.event.FatJet[self.tag['AK8jet']]
        return self.tag['AK8jet']
    
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
        
        ak.to_parquet(array=ak.Array(array), where=os.path.join(output_dir, f'{name}.parq'))
    
    def preselect_HGamma(self):  # __ in prefix means private method
        # at least pass one trigger
        self.pass_cut(cutName='triggered', cut=self.triggered(level='any'))
        
        # b veto
        self.pass_cut(cutName='b-veto', cut=(ak.sum(self.b_tag(reco=False, level='tight'), axis=1)==0))

        # Muon veto
        self.pass_cut(cutName='muon-veto', cut=(ak.sum(self.muon_tag(reco=False), axis=1)==0))
        
        # Electron veto
        self.pass_cut(cutName='electron-veto', cut=(ak.sum(self.electron_tag(reco=False), axis=1)==0))
        
        # Photon >=1
        self.pass_cut(cutName='photon', cut=(ak.sum(self.photon_tag(reco=True), axis=1)>0))
        
        # AK8 jet >=1
        self.pass_cut(cutName='AK8jet', cut=(ak.sum(self.AK8jet_tag(reco=True), axis=1)>0))
        """
        # Photon-Jet cleaning, a very special part so no function definition here
        pj_pair = ak.cartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_index_pair = ak.argcartesian({'photon': self.object['photon'], 'jet': self.object['AK8jet']}, axis=1, nested=False)
        pj_dr = pj_pair.photon.delta_r(pj_pair.jet)
        pj_clean = (pj_dr > self.cutValue['deltaR']['min'])
        photon_index, jet_index = pj_index_pair.photon[pj_clean], pj_index_pair.jet[pj_clean]
        # exactly 1 pair of photon and jet passed jet-cleaning requirement
        self.object['photon'] = self.object['photon'][photon_index]
        self.object['AK8jet'] = self.object['AK8jet'][jet_index]
        self.object['photon-jet'] = self.object['photon'] + self.object['AK8jet']
        # final=True means to drop events not passing all selections
        self.pass_cut(cutName='photon-jet_cleaning', cut=(ak.sum(pj_clean, axis=-1)==1), final=True)
        """
        # Return vars of objects after pre-selection
        self.object['heaviest_jet'] = self.object['AK8jet'][ak.argmax(self.object['AK8jet'].msoftdrop, axis=1, keepdims=True)][:, 0]
        self.object['leading_photon'] = self.object['photon'][ak.argmax(self.object['photon'].pt, axis=1, keepdims=True)][:, 0]
        self.object['photon-jet'] = self.object['heaviest_jet'] + self.object['leading_photon']
        delta_phi = abs(self.object['heaviest_jet'].phi - self.object['leading_photon'].phi)
        delta_phi = ak.min([delta_phi, 2 * np.pi - delta_phi], axis=0)
        final_cut = self.pass_cut(cutName='photon-jet_delta_phi', cut=(delta_phi>2.9), final=True)
        self.store_variables(vars={
            'heaviest_jet': {'pt', 'eta', 'phi', 'mass', 'msoftdrop'},
            'leading_photon': {'pt', 'eta', 'phi', 'mass'},
            'event': {'MET_pt', 'genWeight'},
            'photon-jet': {'pt', 'eta', 'phi', 'mass'},
        })
        
        # Additional vars by special computing
        # self.variables['photon-jet_deltaR'] = self.object['photon'].delta_r(self.object['AK8jet'])
        self.variables['photon-jet_deltaR'] = self.object['heaviest_jet'].delta_r(self.object['leading_photon'])

        return final_cut
    
    def process(self, events: NanoEventsArray) -> dict:
        # initialize
        if self.mode.startswith('data'):
            self.event = self.lumi_mask(events)
            self.cutflow['n_events'] = len(self.event)
        elif self.mode.startswith('mc'):
            self.event = events
            self.cutflow['n_events'] = ak.sum(self.event.genWeight)

        # process
        final_cut = self.preselect_HGamma()
        
        # gen-macthing
        if any(final_cut) and 'mc' in self.mode:
            if 'ZpToHGamma' in self.mode:
                self.variables.update(GenMatch().ZpToHGamma(self.event))
            elif 'QCD' in self.mode:
                final_cut = self.pass_cut(cutName='no prompt photon', cut=GenMatch().QCD(self.event), final=True)
            elif 'GJets' in self.mode:
                final_cut = self.pass_cut(cutName='any prompt photon', cut=GenMatch().GJets(self.event), final=True)
                
        # store output
        if any(final_cut):
            self.to_parquet(array=self.variables)
        return {self.mode: {k: float(v) for (k, v) in self.cutflow.items()}}
    
    @property  # transform method into attribute and make it unchangable to hide _accumulator
    def accumulator(self):
        return self._accumulator
    
    def postprocess(self, accumulator):
        return super().postprocess(accumulator)
