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
    parser.add_argument('-e', '--environment', help='To specify running on which environment', choices=('local', 'condor'))
    parser.add_argument('-o', '--outdir', help='To specify output directory', default='./')
    parser.add_argument('-m', '--mode', help='To specify $type_$year(_$channel) mode', default='mc_2018_ZpToHGamma')
    parser.add_argument('-n', '--ncpu', help='To specify the number of CPUs', default=1)
    args = parser.parse_args()
    return args


def main():
    if len(sys.argv)<4:
        raise ValueError('main() needs three arguments as file, environment, outdir, mode by -f, -e, -o, -m respectively')
    args = parse_commanline()
    
    t0 = time.time()
    run = processor.Runner(
        executor = processor.FuturesExecutor(compression=None, workers=args.ncpu),
        schema = NanoAODSchema,
        savemetrics = True,
        xrootdtimeout = 60 * 30,
        # chunksize = 100_000,
        # maxchunks = None,
    )
    cutflow, metrics = run(
        fileset = {'input': [args.file]},
        treename = 'Events',
        processor_instance = Processor(outdir=args.outdir, environment=args.environment, mode=args.mode),
    )
    print(f'===> Time for processing: {(time.time() - t0)/60} mins')
    print('===> Metrics:\n', metrics)
    print('===> Cutflow:\n', cutflow)


if __name__ == "__main__":
    main()
