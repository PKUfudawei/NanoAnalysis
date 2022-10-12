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

samples = {
    'HGamma': {},
    'GJets': {}
}

base = '/data/pubfs/fudawei/samples/mc/2018'
basedir = {d: os.path.join(base, d) for d in os.listdir(base)}
samples = {s: [] for s in basedir}
for s in basedir:
    for (current_path, dirs, files) in os.walk(basedir[s]):
        for f in files:
            if f.endswith('.root'):
                samples[s].append(os.path.join(current_path, f))
                                        

cutflow={}

t0 = time.time()
cutflow['ZpToHGamma'] = processor.run_uproot_job(
    fileset={'ZpToHGamma': samples['HGamma']},
    treename="Events",
    processor_instance=Processor(outdir=os.path.join('..', 'output', 'HGamma')),
    executor=processor.futures_executor,
    executor_args={"schema": NanoAODSchema, "workers": 36}, # running on $workers cpu cores
)
t1 = time.time()

cutflow['GJets'] = processor.run_uproot_job(
    fileset={'GJets': samples['GJets']},
    treename="Events",
    processor_instance=Processor(outdir=os.path.join('..', 'output', 'GJets')),
    executor=processor.futures_executor,
    executor_args={"schema": NanoAODSchema, "workers": 36}, # running on $workers cpu cores
)

t2 = time.time()

cutflow['HGamma time'] = '%.2f mins'%((t1-t0)/60)
cutflow['GJets time'] = '%.2f mins'%((t2-t1)/60)

with open('../output/cutflow.txt', 'w') as file:
    json.dump(cutflow, file)
