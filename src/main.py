#!/usr/bin/env python3
import sys
import time
import argparse

from coffea import processor
from coffea.nanoevents import NanoAODSchema
from processors.Processor import Processor


def parse_commanline():
    parser = argparse.ArgumentParser(description='Script to check if each condor job is done')
    parser.add_argument('-f', '--file', help='To specify file path',)
    parser.add_argument('-m', '--machine', help='To specify running on which machine', choices=('local', 'condor'))
    parser.add_argument('-o', '--outdir', help='To specify output directory', default='./output')
    args = parser.parse_args()
    return args

def main():
    if len(sys.argv)<3:
        raise ValueError('main() needs three arguments as file, machine and outdir by -f, -m and -o respectively')
    args = parse_commanline()
    
    t0 = time.time()
    cutflow = processor.run_uproot_job(
        fileset={'input': [args.file]},
        treename="Events",
        processor_instance=Processor(outdir=args.outdir, machine=args.machine),
        executor=processor.futures_executor,
        executor_args={"schema": NanoAODSchema, "workers": 1}, # running on $workers cpu cores
    )
    cutflow['time'] = (time.time() - t0)/60
    print(cutflow)
        
if __name__ == "__main__":
    main()
