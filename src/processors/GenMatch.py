import numpy as np
import awkward as ak

from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods.base import NanoEventsArray
from coffea.nanoevents.methods.nanoaod import FatJetArray, GenParticleArray

class GenMatch():
    def __init__(self):
        self.PDGID = {
            "d": 1,
            "u": 2,
            "s": 3,
            "c": 4,
            "b": 5,
            "t": 6,
            "b'": 7,
            "t'": 8,
            "e-": 11,   "e+": -11,
            "ve": 12,   "ve~": -12,
            "m-": 13,   "m+": -13,
            "vm": 14,   "vm~": -14,
            "ta-": 15,  "ta+": -15,
            "vt": 16,   "vt~": -16,
            "ta'-": 17, "ta'+": -17,
            "vt'": 18,  "vt'~": -18,
            "g": 21, ## or 9?
            "a": 22,
            "Z": 23,
            "W+": 24,   "W-":-24,
            "H": 25,
            "Z'": 9906663
        }
        self.GEN_FLAGS = ["fromHardProcess", "isLastCopy"]
        self.particle = {}
        self.childs_of = {} ## childs shorter than children
        self.genVars = {}
        
    def _update(self, name: str, genPart: GenParticleArray, variables: set={'pt', 'eta', 'phi', 'mass', 'pdgId'},
                flatten: bool=True, maxNum: int=1, axis: int=-1, clip: bool=False): ## default genPart shape: (event, particle)
        self.particle[name] = genPart  ## shape: (event, particle) 
        self.childs_of[name] = ak.flatten(genPart.children, axis=2) if flatten else genPart.children 
        ## shape: (event, child_particle) if flatten else (event, particle, child_particle) 
        
        return {
            f'gen_{name}_{var}': ak.pad_none(array=genPart[var], target=maxNum, axis=axis, clip=clip) for var in variables
        }
        
    def HGamma(self, events: NanoEventsArray):  
        """Get gen-info. and the wanted variables of H, WW, W_children and gamma"""
        
        self.genVars["Z'"] = self._update( ## shape: (event, particle)
            genPart = events.GenPart[
                (events.GenPart.pdgId == self.PDGID["Z'"]) &
                events.GenPart.hasFlags(self.GEN_FLAGS)
            ], name="Z'", flatten=True, maxNum=1
        )
        
        self.genVars["H"] = self._update( ## shape: (event, H)
            genPart = self.childs_of["Z'"][
                (self.childs_of["Z'"].pdgId == self.PDGID["H"]) &
                self.childs_of["Z'"].hasFlags(self.GEN_FLAGS)
            ], name="H", flatten=True, maxNum=1
        )
        
        self.genVars["a"] = self._update( ## shape: (event, gamma)
            genPart = self.childs_of["Z'"][
                (self.childs_of["Z'"].pdgId == self.PDGID["a"]) &
                self.childs_of["Z'"].hasFlags(self.GEN_FLAGS)
            ], name="a", flatten=True, maxNum=1
        )
        
        self.genVars["WW"] = self._update( ## shape: (event, WW)
            genPart = self.childs_of["H"][ ## shape: (event, H_childs)
                (ak.num(self.childs_of['H'].pdgId, axis=-1) == 2) &
                ak.all(abs(self.childs_of['H'].pdgId) == self.PDGID["W+"], axis=-1) &
                self.childs_of["H"].hasFlags(self.GEN_FLAGS)
            ], name="WW", flatten=True, maxNum=2
        )
        
        self.genVars["WW_childs"] = self._update( ## shape: (event, WW_childs)
            genPart = self.childs_of["WW"][ ## shape: (event, WW)
                (ak.num(self.childs_of["WW"].pdgId, axis=-1) == 4) &
                self.childs_of["WW"].hasFlags(self.GEN_FLAGS)
            ], name="WW_childs", flatten=True, maxNum=4
        )
        
        HWW_decay_mode = (
            ak.Array([0 for _ in range(len(events))]) + 
            1 * ak.num(abs(self.childs_of['WW_childs'].pdgId)==11, axis=-1) + # num. of electrons in WW_childs
            2 * ak.num(abs(self.childs_of['WW_childs'].pdgId)==13, axis=-1) + # num. of muons in WW_childs
            4 * ak.num(abs(self.childs_of['WW_childs'].pdgId)==15, axis=-1) + # num. of tauons in WW_childs
            8 * ak.num(abs(self.childs_of['WW_childs'].pdgId)<=6, axis=-1)    # num. of quarks in WW_childs
        )
        
        return {
            **self.genVars["Z'"], **self.genVars["H"], **self.genVars["a"], **self.genVars["WW"], **self.genVars["WW_childs"], 
            'HWW_decay_mode': HWW_decay_mode, 
        }
        
        
        
        
        
        
        