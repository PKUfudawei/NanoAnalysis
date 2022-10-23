import awkward as ak

from coffea.nanoevents.methods.base import NanoEventsArray
from coffea.nanoevents.methods.nanoaod import GenParticleArray


class GenMatch():
    def __init__(self) -> None:
        self.PDGID = {
            "d": 1,
            "u": 2,
            "s": 3,
            "c": 4,
            "b": 5,
            "t": 6,
            "b'": 7,
            "t'": 8,
            "e-": 11, "e+": -11,
            "ve": 12, "ve~": -12,
            "m-": 13, "m+": -13,
            "vm": 14, "vm~": -14,
            "ta-": 15, "ta+": -15,
            "vt": 16, "vt~": -16,
            "ta'-": 17, "ta'+": -17,
            "vt'": 18, "vt'~": -18,
            "g": 21,  # or 9?
            "a": 22,
            "Z": 23,
            "W+": 24, "W-": -24,
            "H": 25,
            "Zp": 9906663
        }
        self.object = {}
        self.childs_of = {}  # childs shorter than children
        self.genVars = {}
        "statusFlags().isLastCopyBeforeFSR()                  * 16384 +"
        "statusFlags().isLastCopy()                           * 8192  +"
        "statusFlags().isFirstCopy()                          * 4096  +"
        "statusFlags().fromHardProcessBeforeFSR()             * 2048  +"
        "statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +"
        "statusFlags().isHardProcessTauDecayProduct()         * 512   +"
        "statusFlags().fromHardProcess()                      * 256   +"
        "statusFlags().isHardProcess()                        * 128   +"
        "statusFlags().isDirectHadronDecayProduct()           * 64    +"
        "statusFlags().isDirectPromptTauDecayProduct()        * 32    +"
        "statusFlags().isDirectTauDecayProduct()              * 16    +"
        "statusFlags().isPromptTauDecayProduct()              * 8     +"
        "statusFlags().isTauDecayProduct()                    * 4     +"
        "statusFlags().isDecayedLeptonHadron()                * 2     +"
        "statusFlags().isPrompt()                             * 1      "
        self.GEN_FLAGS = ["fromHardProcess", "isLastCopy"]
        
    def updateParticle(
        self, name: str, genPart: GenParticleArray, flatten: bool = True,
        variables: set = {'pt', 'eta', 'phi', 'mass', 'pdgId'},
        maxNum: int = 1, axis: int = -1, clip: bool = False
    ) -> dict:  # default genPart shape: (event, particle)
        
        self.object[name] = genPart  # shape: (event, particle)
        self.childs_of[name] = ak.flatten(self.object[name].children, axis=2) if flatten else self.object[name].children
        # shape: (event, child_particle) if flatten else (event, particle, child_particle)
        
        return {
            f'gen_{name}_{var}': ak.pad_none(
                array=self.object[name][var], target=maxNum, axis=axis, clip=clip
            ) for var in variables
        }
        
    def ZpToHGamma(self, events: NanoEventsArray) -> dict:
        """Get gen-info. and the wanted variables of H, WW, W_children and gamma"""
        
        self.genVars["Zp"] = self.updateParticle(  # shape: (event, particle)
            genPart = events.GenPart[
                (events.GenPart.pdgId == self.PDGID["Zp"]) &
                events.GenPart.hasFlags(self.GEN_FLAGS)
            ], name="Zp", flatten=True, maxNum=1
        )
        
        self.genVars["H"] = self.updateParticle(  # shape: (event, H)
            genPart = self.childs_of["Zp"][
                (self.childs_of["Zp"].pdgId == self.PDGID["H"]) &
                self.childs_of["Zp"].hasFlags(self.GEN_FLAGS)
            ], name="H", flatten=True, maxNum=1
        )
        
        self.genVars["a"] = self.updateParticle(  # shape: (event, gamma)
            genPart = self.childs_of["Zp"][
                (self.childs_of["Zp"].pdgId == self.PDGID["a"]) &
                self.childs_of["Zp"].hasFlags(self.GEN_FLAGS)
            ], name="a", flatten=True, maxNum=1
        )
        
        self.genVars["WW"] = self.updateParticle(  # shape: (event, WW)
            genPart = self.childs_of["H"][  # shape: (event, H_childs)
                (ak.num(self.childs_of['H'].pdgId, axis=-1) == 2) &
                ak.all(abs(self.childs_of['H'].pdgId) == self.PDGID["W+"], axis=-1) &
                self.childs_of["H"].hasFlags(self.GEN_FLAGS)
            ], name="WW", flatten=True, maxNum=2
        )
        
        self.genVars["WW_childs"] = self.updateParticle(  # shape: (event, WW_childs)
            genPart = self.childs_of["WW"][  # shape: (event, WW)
                (ak.num(self.childs_of["WW"].pdgId, axis=-1) == 4) &
                self.childs_of["WW"].hasFlags(self.GEN_FLAGS)
            ], name="WW_childs", flatten=True, maxNum=4
        )
        
        HWW_decay_mode = (
            ak.Array([0 for _ in range(len(events))]) +
            1 * ak.num(abs(self.childs_of['WW_childs'].pdgId)==11, axis=-1) +  # num. of electrons in WW_childs
            2 * ak.num(abs(self.childs_of['WW_childs'].pdgId)==13, axis=-1) +  # num. of muons in WW_childs
            4 * ak.num(abs(self.childs_of['WW_childs'].pdgId)==15, axis=-1) +  # num. of tauons in WW_childs
            8 * ak.num(abs(self.childs_of['WW_childs'].pdgId)<=6, axis=-1)     # num. of quarks in WW_childs
        )
        H_a_pair = ak.cartesian({'H': self.object['H'], 'a': self.object['a']}, axis=1, nested=False)
        self.genVars['event'] = {
            'gen_H_a': ak.flatten(ak.num(H_a_pair, axis=1)==1, axis=-1),
            'gen_deltaR_H_a': ak.flatten(ak.pad_none(H_a_pair.H.delta_r(H_a_pair.a), target=1, axis=1), axis=-1),
            'gen_HWW_decay_mode': HWW_decay_mode,
            'gen_HWW_a': ak.flatten((HWW_decay_mode>0) * ak.num(H_a_pair, axis=1) == 1, axis=-1),
            'gen_MET_pt': events.GenMET.pt,
        }
 
        return {
            **self.genVars["Zp"], **self.genVars["H"], **self.genVars["a"], **self.genVars["WW"],
            **self.genVars["WW_childs"], **self.genVars['event'],
        }
