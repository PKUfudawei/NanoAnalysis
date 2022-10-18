#!/usr/bin/env python3
import argparse
import os
import glob

def parse_commanline():
    parser = argparse.ArgumentParser(description='Script to check if each condor job is done')
    parser.add_argument('-t', '--type', help='To specify jobs in mc/ or data/', choices=('data', 'mc', '*'), default='*')
    parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
    parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
    parser.add_argument('-j', '--job', help='Specify which job to be submitted', default='*')
    args = parser.parse_args()
    return args

def main() -> None:
    args = parse_commanline()
    ## 4506743.36 from QCD_HT1000to1500_TuneCP5
    ## 4506748.45 from QCD_HT300to500_TuneCP5
    os.system(f"""
        ## myschedd bump
        for jdl in submit/{args.type}/{args.year}/{args.channel}/{args.job}.submit
            do condor_submit $jdl
        done
    """)
    
if __name__ == "__main__":
    main()