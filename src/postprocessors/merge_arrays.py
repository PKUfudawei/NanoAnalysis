#!/usr/bin/env python3
import os, argparse, glob, yaml
import awkward as ak


def parse_commanline():
    parser = argparse.ArgumentParser(description='Do cross-section reweighting on files')
    parser.add_argument('-i', '--indir', help='Specify input directory of arrays', default='../../condor/output/')
    parser.add_argument('-t', '--sample_type', help='Specify jobs in mc/ or data/', default='*', choices=('data', 'mc', '*'))
    parser.add_argument('-y', '--year', help='Specify jobs in which year', choices=('2018', '2017', '2016pre', '2016post'), default='*')
    parser.add_argument('-c', '--channel', help='Specify jobs in which channel', default='*')
    parser.add_argument('-j', '--job', help='Specify which job to be submitted', default='*')
    parser.add_argument('-o', '--outdir', help='Specify output directory of arrays', default=None)
    args = parser.parse_args()
    return args


def merge(cutflow, array_list, yaml_file, parquet_file, merge_input):
    if not os.path.exists(yaml_file) or os.path.getsize(yaml_file) == 0:
        return cutflow, array_list
    merge_input.append(yaml_file)

    with open(yaml_file, 'r', encoding='utf-8') as f:
        new_cutflow = yaml.safe_load(f)
    for (key, value) in new_cutflow.items():
        cutflow[key] = cutflow.get(key, 0) + int(value)
    if os.path.exists(parquet_file) and os.path.getsize(parquet_file) != 0 and new_cutflow.get('final', 0) > 0:
        array_list.append(ak.from_parquet(parquet_file))
        merge_input.append(parquet_file)

    return cutflow, array_list, merge_input

def merge_files(indir, outdir):
    print(indir)
    mode = '_'.join(indir.split('/')[-4:-1])
    print('\tStart postprocessing files:')
    os.system(f"rm -rf {os.path.join(indir.replace('output', 'log'), '*')}")

    cutflow = {}
    arrays = []
    merge_input = []

    for yml_file in glob.glob(os.path.join(indir, '*.yml')):
        parq_file = yml_file.replace('yml', 'parq')
        cutflow, arrays, merge_input = merge(cutflow, arrays, yml_file, parq_file, merge_input)
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    else:
        merged_yaml = os.path.join(outdir, f'{mode}.yaml')
        merged_parquet = os.path.join(outdir, f'{mode}.parquet')
        cutflow, arrays, merge_input = merge(cutflow, arrays, merged_yaml, merged_parquet, merge_input)

    if len(cutflow.keys()) > 0:
        with open(merged_yaml+'.new', 'w', encoding='utf-8') as file:
            yaml.dump(cutflow, file)
    if cutflow.get('final', 0) == 0:
        print('\tNo events passed final cut!')
    elif len(arrays) > 0:
        print('\tMerging arrays and cutflows...')
        ak.to_parquet(ak.concatenate(arrays, axis=0), merged_parquet+'.new')
        print(f'\tMegred files are stored in {outdir}')
    

    print('\tFinished.')
    for f in merge_input:
        os.remove(f)
    os.rename(merged_yaml+'.new', merged_yaml)
    os.rename(merged_parquet+'.new', merged_parquet)


def main():
    args = parse_commanline()
    if args.outdir is None:
        args.outdir = args.indir
    input_dirs = os.path.join(args.indir, args.sample_type, args.year, args.channel, args.job)

    for input_dir in set(glob.glob(input_dirs)):
        output_dir = os.path.join(args.outdir, *input_dir.split('/')[-4])
        merge_files(indir=input_dir, outdir=output_dir)


if __name__ == "__main__":
    main()
