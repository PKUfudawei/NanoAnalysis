#!/usr/bin/env python3
import os
import sys
import time
import argparse
import awkward as ak
import pickle
import uproot


from coffea import processor
from coffea.nanoevents import NanoAODSchema
from processors.Processor import Processor


def parse_commandline():
    parser = argparse.ArgumentParser(description='Main script to run preprocessors and processors')
    parser.add_argument('-f', '--file', help='To specify file path', default=None)
    parser.add_argument('-p', '--param_dir', help='To specify input parameters directory', default='./parameters/')
    parser.add_argument('-o', '--outdir', help='To specify output directory', default='./')
    parser.add_argument('-m', '--mode', help='To specify $type_$year_$channel mode', default='mc_2018_ZpToHGamma')
    parser.add_argument('-n', '--ncpu', help='To specify the number of CPUs', default=1)
    args = parser.parse_args()
    return args


def main() -> None:
    # prepare
    uproot.open.defaults["xrootd_handler"] = uproot.MultithreadedXRootDSource
    if len(sys.argv) < 2:
        raise ValueError('main() needs two necessary arguments as mode, file by -m, -f respectively')
    args = parse_commandline()

    if args.file is None:
        for f in os.listdir('.'):
            if f.endswith('.root'):
                file = f
                break
    else:
        os.system(f"rsync {args.file} .")
        file = os.path.join('.', args.file.split('/')[-1])

    # run
    t0 = time.time()
    print(f'===> Running on file: {file}')

    ## processing
    run = processor.Runner(
        executor=processor.FuturesExecutor(compression=None, workers=args.ncpu),
        schema=NanoAODSchema,
        savemetrics=True,
        xrootdtimeout=60 * 30,
        # chunksize = 100_000,
        # maxchunks = None,
    )
    stats, metrics = run(
        fileset={'input': [f'{file}']}, treename='Events',
        processor_instance=Processor(mode=args.mode, outdir=args.outdir, param_dir=args.param_dir),
    )

    ## post-processing
    result = []
    for i in os.listdir(args.outdir):
        if i.endswith('.parq'):
            result.append(ak.from_parquet(os.path.join(args.outdir, i)))

    if len(result) > 1:
        result = ak.concatenate(result, axis=0)
    elif len(result) == 1:
        result = result[0]
    else:
        result = ak.Array({})

    ak.to_parquet(result, where=os.path.join(args.outdir, 'output.parq'))
    with open(os.path.join(args.outdir, 'stats.pkl'), 'wb') as f:
        pickle.dump(stats, f)

    print()
    print(f'===> Removed {os.system("rm -rf *Events*.parq")} intermediate parquets')
    print(f'===> Time for full-processing: {(time.time() - t0)/60} mins')
    print('===> Metrics:\n', metrics)
    print('===> Statistics:\n', stats)


if __name__ == "__main__":
    main()
