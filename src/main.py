#!/usr/bin/env python3
import os
import sys
import time
import argparse
import awkward as ak
import json

from coffea import processor
from coffea.nanoevents import NanoAODSchema
from processors.Processor import Processor


def parse_commanline():
    parser = argparse.ArgumentParser(description='Script to check if each condor job is done')
    parser.add_argument('-f', '--file', help='To specify file path',)
    parser.add_argument('-o', '--outdir', help='To specify output directory', default='./')
    parser.add_argument('-m', '--mode', help='To specify $type_$year(_$channel) mode', default='mc_2018_ZpToHGamma')
    parser.add_argument('-n', '--ncpu', help='To specify the number of CPUs', default=1)
    args = parser.parse_args()
    return args


def main():
    if len(sys.argv)<2:
        raise ValueError('main() needs three arguments as file, mode by -f, -m respectively')
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
    stats, metrics = run(
        fileset = {'input': [args.file]},
        treename = 'Events',
        processor_instance = Processor(outdir=args.outdir, mode=args.mode),
    )
    result = []
    for i in os.listdir(args.outdir):
        if i.endswith('.parq')>0:
            result.append(ak.from_parquet(i))
    if len(result)>1:
        result = ak.concatenate(result, axis=0)
    elif len(result)==1:
        result = result[0]
    else:
        result = ak.Array({})
    ak.to_parquet(result, where='./output.parq')
    with open('./stats.json', 'w', encoding ='utf-8') as f:
        json.dump(stats, f)
    
    print(f'===> Removed {os.system("rm -rf *Events*.parq")} intermediate parquets')
    print(f'===> Time for full-processing: {(time.time() - t0)/60} mins')
    print('===> Metrics:\n', metrics)
    print('===> Statistics:\n', stats)


if __name__ == "__main__":
    main()
