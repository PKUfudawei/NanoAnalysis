#!/usr/bin/env python3
import os
import yaml
import argparse
import awkward as ak
import glob


def parse_commanline():
    parser = argparse.ArgumentParser(description='Do cross-section reweighting on files')
    parser.add_argument('-d', '--directory', help='To specify file directory', default='../../condor/output/')
    parser.add_argument('-t', '--sample_type', help='To specify jobs in mc/ or data/', default='mc')
    parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
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
        print('!!! File zero size:', file, '!!!')
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

    for current_path in dirs:
        print(f'=> Start postprocessing files in {current_path}')
        n_raw_events = 0
        dataset = current_path.split('/')[-1]
        channel = current_path.split('/')[-2]
        year = current_path.split('/')[-3]

        stats_file = os.path.join(current_path, '_'.join(['mc', year, channel])) + '.yaml'
        output_file = stats_file.replace('yaml', 'parq')
        if not os.path.exists(stats_file) or not os.path.exists(output_file) or os.path.getsize(stats_file) == 0 or os.path.getsize(output_file) == 0:
            print('!!! File zero size in: ', current_path, '!!!')
            continue

        with open(stats_file, 'r', encoding='utf-8') as f:
            stats = yaml.safe_load(f)
        n_raw_events += stats['n_events']
        print(f'==> Finish parsing yaml files with n_raw_events={n_raw_events}')

        cross_section_reweighting(output_file, LUMI[year], eval(str(X_SECTION[channel][dataset])), n_raw_events)
        print(f'===> Finish postprocessing parquet files with lumi={LUMI[year]}, x-section={X_SECTION[channel][dataset]}')


if __name__ == "__main__":
    main()
