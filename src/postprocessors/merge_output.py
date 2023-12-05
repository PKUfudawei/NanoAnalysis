#!/usr/bin/env python3
import os
import pickle
import argparse
import awkward as ak
import glob
import yaml


def parse_commanline():
    parser = argparse.ArgumentParser(description='Do cross-section reweighting on files')
    parser.add_argument('-d', '--directory', help='To specify file directory', default='../../condor/output/')
    parser.add_argument('-t', '--sample_type', help='To specify jobs in mc/ or data/', default='*')
    parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
    parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
    parser.add_argument('-j', '--job', help='Specify which job to be submitted', default='*')
    args = parser.parse_args()
    return args


def main():
    args = parse_commanline()

    dirs = os.path.join(args.directory, args.sample_type, args.year, args.channel, args.job)
    dirs = set(glob.glob(dirs))

    for current_path in dirs:
        mode = '_'.join(current_path.split('/')[-4:-1])
        print(f'=> Start postprocessing files in {current_path}')
        os.system(f"rm -rf {os.path.join(current_path.replace('output', 'log'), '*')}")

        if os.path.exists(os.path.join(current_path, f'{mode}.parq')):
            print(f'==> Jump {current_path}')
            continue

        cutflow = {}
        arrays = []
        for f in os.listdir(current_path):
            if os.path.getsize(os.path.join(current_path, f)) == 0 or (not f.endswith('.pkl')) or ('-' not in f):
                continue

            with open(os.path.join(current_path, f), 'rb') as file:
                stats = pickle.load(file)
            for (key, value) in stats.items():
                cutflow[key] = cutflow.get(key, 0) + value
                continue
            if 'final' in stats and stats['final'] != 0:
                arrays.append(ak.from_parquet(os.path.join(current_path, f.replace('pkl', 'parq'))))

        print(f'==> Merging parquet and pickle files in {current_path}')
        cutflow = {key: int(value) for (key, value) in cutflow.items()}
        if len(cutflow.keys()) > 0:
            with open(os.path.join(current_path, f'{mode}.yaml'), 'w', encoding='utf-8') as file:
                yaml.dump(cutflow, file)
        output_file = os.path.join(current_path, f'{mode}.parquet')
        if len(arrays) > 0:
            ak.to_parquet(ak.concatenate(arrays, axis=0), output_file)

        print(f'===> Finish merging parquet and pickle files in {current_path}')
        os.system(f"rm -rf {current_path}/*.pkl {current_path}/*.parq")
        os.system(f"mv {current_path}/{mode}.parquet {current_path}/{mode}.parq")


if __name__ == "__main__":
    main()
