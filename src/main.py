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
    parser.add_argument('-c', '--channel', help='To specify gen-matching mode', default='ZpToHGamma')
    args = parser.parse_args()
    return args

def main():
    if len(sys.argv)<4:
        raise ValueError('main() needs three arguments as file, machine, outdir, channel by -f, -m, -o, -c respectively')
    args = parse_commanline()
    
    t0 = time.time()
    run = processor.Runner(
        executor = processor.FuturesExecutor(compression=None, workers=1),
        schema = NanoAODSchema,
        savemetrics = True,
        xrootdtimeout = 60*30,
        #chunksize = 100_000,
        #maxchunks = None,
    )
    cutflow, metrics = run(
        fileset = {'input': [args.file]},
        treename = 'Events',
        processor_instance = Processor(outdir=args.outdir, machine=args.machine, channel=args.channel),
    )
    print(f'===> Time for processing: {(time.time() - t0)/60} mins')
    print('===> Metrics:\n', metrics)
    print('===> Cutflow:\n', cutflow)
        
if __name__ == "__main__":
    main()
