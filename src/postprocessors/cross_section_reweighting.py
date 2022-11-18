#!/usr/bin/env python3
import os
import json
import argparse
import awkward as ak


def parse_commanline():
    parser = argparse.ArgumentParser(description='Do cross-section reweighting on files')
    parser.add_argument('-d', '--dir', help='To specify file directory', default='../../condor/output/mc')
    args = parser.parse_args()
    return args


def cross_section_reweighting(file, lumi, x_section, n_events):
    array = ak.from_parquet(file)
    if len(array)==0:
        return
    array['event_weight'] = array.event_genWeight * x_section * lumi / n_events
    ak.to_parquet(ak.Array(array), file)
    return array.event_weight


def main():
    args = parse_commanline()
    with open('../json/luminosity.json', 'r', encoding ='utf-8') as f:
        LUMI = json.load(f)
    with open('../json/cross_section.json', 'r', encoding ='utf-8') as f:
        X_SECTION = json.load(f)
    
    for (current_path, dirs, files) in os.walk(args.dir):
        if (len(current_path.split('/'))<5) or (current_path.split('/')[-5] != 'output'):
            continue
        
        print(f'===> Start postprocessing files in {current_path}')
        n_raw_events = 0
        dataset = current_path.split('/')[-1]
        year = current_path.split('/')[-3]
        for f in files:
            if f.endswith('.json'):
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
