#!/usr/bin/env python3
import os, sys, time, argparse, uproot, yaml
import awkward as ak


from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from processors.Processor import Processor
from processors.TriggerProcessor import TriggerProcessor
NanoAODSchema.warn_missing_crossrefs = False


def parse_commandline():
    parser = argparse.ArgumentParser(description='Main script to run preprocessors and processors')
    parser.add_argument('-f', '--file', help='To specify file path', default=None)
    parser.add_argument('-p', '--param_dir', help='To specify input parameters directory', default='./parameters/')
    parser.add_argument('-o', '--outdir', help='To specify output directory', default='./')
    parser.add_argument('-n', '--name', help='To specify output name', default='out')
    parser.add_argument('-m', '--mode', help='To specify $type_$year_$channel mode', default='mc_2018_ZpToHGamma')
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
        if ':' in args.file:
            os.system(f"xrdcp {args.file} .")  
            file = args.file.split('/')[-1]
        else:
            file = args.file

    # run
    t0 = time.time()
    print(f'===> Running on file: {file}')

    ## processing
    events = NanoEventsFactory.from_root(
        {file: 'Events'}, schemaclass=NanoAODSchema, mode='eager',
    ).events()
    p = Processor(outdir=args.outdir, mode=args.mode, param_dir=args.param_dir)
    cutflow = p.process(events)
    stats = {k: float(v) for (k, v) in cutflow.items()}

    ## post-processing
    result = []
    for i in os.listdir(args.outdir):
        if i.endswith('.parq'):
            result.append(ak.from_parquet(os.path.join(args.outdir, i)))

    if len(result) >= 1:
        result = ak.concatenate(result, axis=0)
        ak.to_parquet(result, destination=os.path.join(args.outdir, f'{args.name}.parquet'))
    else:
        with open(os.path.join(args.outdir, f'{args.name}.parquet'), 'w') as f:
            pass

    with open(os.path.join(args.outdir, f'{args.name}.yaml'), 'w', encoding='utf-8') as f:
        yaml.dump(stats, f)

    print()
    print(f'===> Removed {os.system("rm -rf *_*_*.parq")} intermediate parquets')
    print('===> Time for full-processing: %.2f mins'%((time.time() - t0)/60))
    print('===> Statistics:\n', stats)


if __name__ == "__main__":
    main()
