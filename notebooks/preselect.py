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
for (current_path, dirs, files) in os.walk(top='/data/pubfs/fudawei/samples/HGamma', topdown=True):
    ## topdown = True/False mean breadth/depth first
    ## current_path is the absolute path
    description = current_path.split('/')[-1]
    if description.endswith('GeV'):
        samples['HGamma'][description] = [os.path.join(current_path, f) for f in files if f.endswith('.root')]
        
for (current_path, dirs, files) in os.walk(top='/data/pubfs/fudawei/samples/GJets', topdown=True):
    ## topdown = True/False mean breadth/depth first
    ## current_path is the absolute path
    description = current_path.split('/')[-1]
    if description.startswith('HT'):
        samples['GJets'][description] = [os.path.join(current_path, f) for f in files if f.endswith('.root')] 

cutflow={}

t0 = time.time()
cutflow['HGamma'] = processor.run_uproot_job(
    fileset=samples['HGamma'],
    treename="Events",
    processor_instance=Processor(outdir=os.path.join('..', 'output', 'HGamma')),
    executor=processor.futures_executor,
    executor_args={"schema": NanoAODSchema, "workers": 36}, # running on $workers cpu cores
)
t1 = time.time()

cutflow['GJets'] = processor.run_uproot_job(
    fileset=samples['GJets'],
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
