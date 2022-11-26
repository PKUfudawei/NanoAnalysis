#!/usr/bin/env python3
import argparse
import os
import glob


def parse_commanline():
    parser = argparse.ArgumentParser(description='Script to check if each condor job is done')
    parser.add_argument('-d', '--directory', help='To specify base directory', default=os.path.abspath('filelists'))
    parser.add_argument('-t', '--type', help='To specify jobs in mc/ or data/', choices=('data', 'mc', '*'), default='*')
    parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
    parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
    parser.add_argument('-j', '--job', help='Specify which job to be submitted', default='*')
    parser.add_argument('-m', '--mode', help='Specify mode')
    args = parser.parse_args()
    return args


def main() -> None:
    args = parse_commanline()
    filelists = glob.glob(os.path.join(args.directory, args.type, args.year, args.channel, args.job + '.txt'))
    print(filelists)
    for filelist in filelists:
        with open(filelist, 'r') as f:
            files = f.read().splitlines()
        files = [f.replace('root://cms-xrd-global.cern.ch/', '/eos/cms') for f in files]
        count = 0
        for f in files:
            os.system(f'../src/main.py -m {args.mode} -f {f}')
            os.system(f"mv output.parq {filelist.replace('filelist', 'output').replace('.txt', f'/{count}.parq')}")
            os.system(f"mv stats.json {filelist.replace('filelist', 'output').replace('.txt', f'/{count}.json')}")
            count += 1
        
   
if __name__ == "__main__":
    main()
