#!/usr/bin/env python3
import os
import yaml
import argparse
import awkward as ak
import glob


def parse_commanline():
    parser = argparse.ArgumentParser(description='Do cross-section reweighting on files')
    parser.add_argument('-d', '--directory', help='To specify file directory', default='../../condor/output/')
    parser.add_argument('-t', '--sample_type', help='To specify jobs in mc/ or data/', default='*', choices=('mc', 'data'))
    parser.add_argument('-y', '--year', help='To specify jobs in which year', choices=('2018', '2017', '2016pre', '2016post'), default='*')
    parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
    parser.add_argument('-j', '--job', help='Specify which job to be submitted', default='*')
    parser.add_argument('-p', '--param_dir', help='Specify input parameters directory', default='../parameters')
    args = parser.parse_args()
    return args


def cross_section_reweighting(file, lumi, x_section, n_events):
    """
    lumi: fb^-1
    x_section: pb
    """
    if os.path.getsize(file) == 0 or len(ak.from_parquet(file)) == 0:
        print(f'\tEmpty file: {file}')
        return
    array = ak.from_parquet(file)
    array['event_final_weight'] = array['event_weight'] * x_section * lumi * 1000 / n_events
    ak.to_parquet(ak.Array(array), file)
    return array.event_final_weight


def main():
    args = parse_commanline()
    with open(f'{args.param_dir}/luminosity.yaml', 'r', encoding='utf-8') as f:
        LUMI = yaml.safe_load(f)
    with open(f'{args.param_dir}/cross-section.yaml', 'r', encoding='utf-8') as f:
        X_SECTION = yaml.safe_load(f)

    dirs = os.path.join(args.directory, args.sample_type, args.year, args.channel, args.job)
    dirs = set(glob.glob(dirs))

    for folder in dirs:
        print(folder)
        print('\tStart postprocessing files')
        n_raw_events = 0
        sample_type, year, channel, dataset = folder.split('/')[-4:]

        stats_file = os.path.join(folder, '_'.join(folder.split('/')[-4:-1]))+'.yaml'
        output_file = stats_file.replace('yaml', 'parq')
        if not os.path.exists(stats_file) or not os.path.exists(output_file) or os.path.getsize(stats_file) == 0 or os.path.getsize(output_file) == 0:
            print(f'\tEmpty directory!')
            continue

        with open(stats_file, 'r', encoding='utf-8') as f:
            stats = yaml.safe_load(f)
        n_raw_events += stats['n_events']
        print(f'\tn_raw_events={n_raw_events}')

        if sample_type == 'mc':
            cross_section_reweighting(output_file, LUMI[year], eval(str(X_SECTION[channel][dataset])), n_raw_events)
            print(f'\tFinished, lumi={LUMI[year]}, x-section={X_SECTION[channel][dataset]}pb')
        elif sample_type == 'data':
            cross_section_reweighting(output_file, 1, 1, 1000)
            print(f'\tFinished, weight = 1 for data')


if __name__ == "__main__":
    main()
