import numpy as np
import awkward as ak
import uproot

from coffea import processor
from coffea.nanoevents.methods import candidate
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.analysis_tools import PackedSelection
from coffea.nanoevents.methods.base import NanoEventsArray
from coffea.nanoevents.methods.nanoaod import FatJetArray, GenParticleArray

class Processor(processor.ProcessorABC):
    def __init__(self) -> None:
        super().__init__()
    
    def _HGamma(self, events: NanoEventsArray, trigger: str):
        # Trigger
        events = events[events[trigger]]
        pass
        
    
    
    def process(self, events: NanoEventsArray):

        
        pass