import awkward as ak
import uproot
import numpy as np
import pandas as pd
import os
import time
import json

from coffea import processor
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

import sys
sys.path.append("..")
from src.processors.Processor import Processor

base = '/stash/user/fudawei/samples/mc/2018'
basedir = {d: os.path.join(base, d) for d in os.listdir(base)}
samples = {s: [] for s in basedir}
for s in basedir:
    for (current_path, dirs, files) in os.walk(basedir[s]):
        for f in files:
            if f.endswith('.root'):
                samples[s].append(os.path.join(current_path, f))
print(samples)

cutflow = {}
t0 = time.time()
                  
def parallelProcess(samples, name):
    global cutflow, t0
    cutflow[name] = processor.run_uproot_job(
        fileset={name: samples[name]},
        treename="Events",
        processor_instance=Processor(outdir=f'../output/{name}'),
        executor=processor.futures_executor,
        executor_args={"schema": NanoAODSchema, "workers": 24}, # running on $workers cpu cores
    )
    cutflow['time_'+name] = (time.time() - t0)/60
    with open('../output/cutflow.txt', 'w') as file:
        json.dump(cutflow, file)

parallelProcess(samples=samples, name='ZpToHGamma')
parallelProcess(samples=samples, name='GJets')
parallelProcess(samples=samples, name='WJetsToQQ')
parallelProcess(samples=samples, name='ZJetsToQQ')
parallelProcess(samples=samples, name='QCD')

