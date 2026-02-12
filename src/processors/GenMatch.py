#!/usr/bin/env python3
import awkward as ak
import numpy as np

from coffea.nanoevents.methods.base import NanoEventsArray
from coffea.nanoevents.methods.nanoaod import GenParticleArray, GenParticle
from coffea.analysis_tools import PackedSelection


class GenMatch():
    def __init__(self, events: NanoEventsArray) -> None:
        self.PDGID = {
            "d": 1, "d~": -1,
            "u": 2, "u~": -2,
            "s": 3, "s~": -3,
            "c": 4, "c~": -4,
            "b": 5, "b~": -5,
            "t": 6, "t~": -6,
            "b_prime": 7,
            "t_prime": 8,
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
        self.event = events
        self.particle = {}
        self.childs_of = {}  # childs shorter than children
        self.gen_cuts = PackedSelection()
        self.variables = {}
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
        # https://github.com/cms-sw/cmssw/blob/master/DataFormats/HepMCCandidate/interface/GenStatusFlags.h

    def find_lastcopy(self, particle: GenParticle):  # not in use at present
        def iterate(particle):
            if particle.hasFlags('isLastCopy'):
                return particle
            else:
                childs = particle.children[particle.children.pdgId == particle.pdgId]
                for p in childs:
                    return iterate(p)
        return iterate(particle=particle)

    def update_particle(
        self, name: str, candidates: GenParticleArray, cut: ak.Array, flatten: bool = True,
        num: int = 1, variables: set = {'pt', 'eta', 'phi', 'mass', 'pdgId'},
    ) -> None:  # default particles shape: (event, particle)
        self.particle[name] = ak.pad_none(candidates[cut & (ak.sum(cut, axis=1) == num)], target=num, axis=1)  # shape: (event, particle)
        self.childs_of[name] = ak.flatten(self.particle[name].children, axis=2) if flatten else self.particle[name].children
        # shape: (event, child_particle) if flatten else (event, particle, child_particle)
        if len(variables) > 0:
            self.variables.update({f'gen_{name}_{var}': self.particle[name][var] for var in variables})

    def gluon(self) -> np.array:
        return

    def quark(self) -> np.array:
        return

    def top(self) -> np.array:  # single top process
        gen_cuts = PackedSelection()

        self.update_particle(  # select one t or t~
            candidates=self.event.GenPart, cut=(
                (abs(self.event.GenPart.pdgId) == self.PDGID["t"]) &
                self.event.GenPart.hasFlags(["fromHardProcess", "isLastCopy"])
            ), name='t', num=1, flatten=True, variables=set()
        )
        gen_cuts.add('N_top == 1', ak.num(self.particle['t'], axis=1) == 1)
        gen_cuts.add('top has 2 childs', ak.num(self.childs_of['t'], axis=1) == 2)

        self.update_particle(  # select one W with the same sign as top
            candidates=self.childs_of['t'], cut=(
                (self.childs_of['t'].pdgId == self.PDGID['W+'] * np.sign(self.genVars['t']['pdgId'])) &
                self.childs_of['t'].hasFlags(["fromHardProcess"])
            ), name='W', num=1, flatten=True, variables=set()
        )
        gen_cuts.add('N_W == 1', ak.num(self.particle['W'], axis=1) == 1)

        self.update_particle(  # select one b with the opposite sign as top
            candidates=self.childs_of['t'], cut=(
                (self.childs_of['t'].pdgId == self.PDGID['b~'] * np.sign(self.genVars['t']['pdgId'])) &
                self.childs_of['t'].hasFlags(["fromHardProcess"])
            ), name='b', num=1, flatten=True, variables=set()
        )
        gen_cuts.add('N_b == 1', ak.num(self.particle['b'], axis=1) == 1)

        return gen_cuts.all(*gen_cuts.names)

    def HGamma(self) -> dict:  # ZprimeToHGamma
        self.update_particle(  # shape: (event, particle)
            candidates=self.event.GenPart, cut=(
                (self.event.GenPart.pdgId == self.PDGID['Zp']) &
                self.event.GenPart.hasFlags(['fromHardProcess', 'isLastCopy'])
            ), name='Zp', flatten=True, num=1, variables={'pt', 'eta', 'phi', 'mass', 'pdgId'}
        )
        self.gen_cuts.add('N_Zp == 1', ak.num(self.particle['Zp'], axis=1) == 1)
        self.gen_cuts.add('Zp to H Gamma', (
            (ak.num(self.childs_of['Zp'], axis=1) == 2) &
            (ak.sum(self.childs_of['Zp'].pdgId == self.PDGID['a'], axis=1) == 1) &
            (ak.sum(self.childs_of['Zp'].pdgId == self.PDGID['H'], axis=1) == 1)
        ))

        self.update_particle(  # shape: (event, gamma)
            candidates=self.childs_of['Zp'], cut=(
                (self.childs_of['Zp'].pdgId == self.PDGID['a']) &
                self.childs_of['Zp'].hasFlags(['fromHardProcess', 'isLastCopy'])
            ), name='a', flatten=True, num=1, variables=set()
        )

        self.update_particle(  # shape: (event, H)
            candidates=self.event.GenPart, cut=(
                (self.event.GenPart.pdgId == self.PDGID['H']) &
                self.event.GenPart.hasFlags(['fromHardProcess', 'isLastCopy'])
            ), name='H', flatten=True, num=1, variables=set()
        )
        self.gen_cuts.add('H to 2 childs', ak.num(self.childs_of['H'], axis=1) == 2)
        self.gen_cuts.add('H to WW', ak.sum(abs(self.childs_of['H'].pdgId) == self.PDGID['W+'], axis=1) == 2)
        self.gen_cuts.add('H to bb', ak.sum(abs(self.childs_of['H'].pdgId) == self.PDGID['b'], axis=1) == 2)
        self.gen_cuts.add('H to cc', ak.sum(abs(self.childs_of['H'].pdgId) == self.PDGID['c'], axis=1) == 2)
        self.gen_cuts.add('H to ZZ', ak.sum(abs(self.childs_of['H'].pdgId) == self.PDGID['Z'], axis=1) == 2)
        self.gen_cuts.add('H to tautau', ak.sum(abs(self.childs_of['H'].pdgId) == self.PDGID['ta-'], axis=1) == 2)
        self.gen_cuts.add('H to gammagamma', ak.sum(abs(self.childs_of['H'].pdgId) == self.PDGID['a'], axis=1) == 2)

        self.update_particle(  # shape: (event, W)
            candidates=self.event.GenPart, cut=(  # shape: (event, H_childs)
                (abs(self.event.GenPart.pdgId) == self.PDGID['W+']) &
                self.event.GenPart.hasFlags(['fromHardProcess', 'isLastCopy'])
            ), name='W', flatten=True, num=2, variables=set()
        )

        self.update_particle(  # shape: (event, W_childs)
            candidates=self.childs_of['W'], cut=(  # shape: (event, WW)
                self.childs_of['W'].hasFlags(['fromHardProcess'])
            ), name='WW_childs', flatten=True, num=4, variables=set()
        )
        self.gen_cuts.add('WW to 4 childs', ak.num(self.childs_of['W'], axis=1) == 4)

        HWW_decay_mode = (
            ak.Array([0 for _ in range(len(self.event))]) +
            1 * ak.sum(abs(self.particle['WW_childs'].pdgId) == self.PDGID['e-'], axis=-1) +  # num. of electrons in WW_childs
            2 * ak.sum(abs(self.particle['WW_childs'].pdgId) == self.PDGID['m-'], axis=-1) +  # num. of muons in WW_childs
            4 * ak.sum(abs(self.particle['WW_childs'].pdgId) == self.PDGID['ta-'], axis=-1) +  # num. of tauons in WW_childs
            8 * ak.sum(abs(self.particle['WW_childs'].pdgId) <= self.PDGID['t'], axis=-1)     # num. of quarks in WW_childs
        )

        base = ('N_Zp == 1', 'Zp to H Gamma')
        base_H2children = (*base, 'H to 2 childs')
        return {
            **self.variables,
            'gen_ZpToHGamma': self.gen_cuts.all(*base),
            'gen_ZpToH(WW)Gamma': self.gen_cuts.all(*base_H2children, 'H to WW', 'WW to 4 childs'),
            'gen_ZpToH(bb)Gamma': self.gen_cuts.all(*base_H2children, 'H to bb'),
            'gen_ZpToH(cc)Gamma': self.gen_cuts.all(*base_H2children, 'H to cc'),
            'gen_ZpToH(ZZ)Gamma': self.gen_cuts.all(*base_H2children, 'H to ZZ'),
            'gen_ZpToH(tautau)Gamma': self.gen_cuts.all(*base_H2children, 'H to tautau'),
            'gen_ZpToH(gammagamma)Gamma': self.gen_cuts.all(*base_H2children, 'H to gammagamma'),
            'gen_HWW_decay_mode': HWW_decay_mode,
        }

    def ZGamma(self) -> dict:  # GluGluToZGamma
        self.update_particle(  # shape: (event, Z)
            candidates=self.event.GenPart, cut=(
                (self.event.GenPart.pdgId == self.PDGID['Z']) &
                self.event.GenPart.hasFlags(['fromHardProcess', 'isLastCopy'])
            ), name='Z', flatten=True, num=1, variables=set()
        )
        self.gen_cuts.add('Z to 2 childs', ak.num(self.childs_of['Z'], axis=1) == 2)
        self.gen_cuts.add('Z to bb', ak.sum(abs(self.childs_of['Z'].pdgId)==self.PDGID['b'], axis=1) == 2)
        self.gen_cuts.add('Z to cc', ak.sum(abs(self.childs_of['Z'].pdgId)==self.PDGID['c'], axis=1) == 2)
        self.gen_cuts.add('Z to qq', ak.sum(abs(self.childs_of['Z'].pdgId)<=self.PDGID['s'], axis=1) == 2)

        self.update_particle(  # shape: (event, gamma)
            candidates=self.event.GenPart, cut=(
                (self.event.GenPart.pdgId == self.PDGID['a']) &
                self.event.GenPart.hasFlags(['fromHardProcess', 'isLastCopy'])
            ), name='a', flatten=True, num=1, variables=set()
        )
        self.gen_cuts.add('One photon', ak.num(self.particle['a'], axis=1) == 1)

        return {
            **self.variables,
            'gen_GluGluToZGamma': self.gen_cuts.all('Z to 2 childs', 'One photon'),
            'gen_GluGluToZ(bb)Gamma': self.gen_cuts.all('Z to 2 childs', 'One photon', 'Z to bb'),
            'gen_GluGluToZ(cc)Gamma': self.gen_cuts.all('Z to 2 childs', 'One photon', 'Z to cc'),
            'gen_GluGluToZ(qq)Gamma': self.gen_cuts.all('Z to 2 childs', 'One photon', 'Z to qq'),
        }

    def all_fake_photon(self) -> ak.Array:
        return ~ak.any(self.event.GenPart[abs(self.event.GenPart.pdgId) == 22].hasFlags(['isPrompt']), axis=1)

    def any_true_photon(self) -> ak.Array:
        return ak.any(self.event.GenPart[abs(self.event.GenPart.pdgId) == 22].hasFlags(['isPrompt']), axis=1)
