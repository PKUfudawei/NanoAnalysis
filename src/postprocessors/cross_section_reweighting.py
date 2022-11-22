#!/usr/bin/env python3
import os
import json
import argparse
import awkward as ak
import numpy as np
import pandas as pd
import glob


def parse_commanline():
    parser = argparse.ArgumentParser(description='Do cross-section reweighting on files')
    parser.add_argument('-d', '--directory', help='To specify file directory', default='../../condor/output/')
    parser.add_argument('-t', '--type', help='To specify jobs in mc/ or data/', default='mc')
    parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
    parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
    parser.add_argument('-j', '--job', help='Specify which job to be submitted', default='*')
    args = parser.parse_args()
    return args


def cross_section_reweighting(file, lumi, x_section, n_events):
    if os.path.getsize(file) == 0:
        print('!!! File zero size:', file, '!!!')
        return
    if len(pd.read_parquet(file))==0:
        return
    array = ak.from_parquet(file)
    array['event_weight'] = np.sign(array.event_genWeight) * x_section * lumi / n_events
    ak.to_parquet(ak.Array(array), file)
    return array.event_weight


def main():
    args = parse_commanline()
    with open('../json/luminosity.json', 'r', encoding ='utf-8') as f:
        LUMI = json.load(f)
    with open('../json/cross_section.json', 'r', encoding ='utf-8') as f:
        X_SECTION = json.load(f)
    
    dirs = os.path.join(args.directory, args.type, args.year, args.channel, args.job)
    dirs = set(glob.glob(dirs))
    
    for current_path in dirs:
        print(f'===> Start postprocessing files in {current_path}')
        n_raw_events = 0
        dataset = current_path.split('/')[-1]
        year = current_path.split('/')[-3]
        files = os.listdir(current_path)
        for f in files:
            if f.endswith('.json'):
                if os.path.getsize(os.path.join(current_path, f)) == 0:
                    print('!!! File zero size:', os.path.join(current_path, f), '!!!')
                    continue
                with open(os.path.join(current_path, f), 'r', encoding ='utf-8') as f:
                    stats = json.load(f)
                n_raw_events += list(stats.values())[0]['n_events']

        print(f'Finish parsing json files with n_raw_events={n_raw_events}')
        for f in files:
            if f.endswith('.parq'):
                cross_section_reweighting(os.path.join(current_path, f), LUMI[year], X_SECTION[dataset], n_raw_events)
        print(f'Finish postprocessing parquet files with lumi={LUMI[year]}, x-section={X_SECTION[dataset]}')


if __name__ == "__main__":
    main()
