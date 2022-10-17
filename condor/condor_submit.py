#!/usr/bin/env python3
import argparse
import os

def parse_commanline():
    parser = argparse.ArgumentParser(description='Script to check if each condor job is done')
    parser.add_argument('-t', '--type', help='To specify jobs in mc/ or data/', choices=('data', 'mc', '*'), default='*')
    parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
    parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
    args = parser.parse_args()
    return args

def main() -> None:
    args = parse_commanline()
    
    os.system(f"""
        ## myschedd bump
        for jdl in submit/{args.type}/{args.year}/{args.channel}/*.submit
            do condor_submit $jdl
        done
    """)
    
if __name__ == "__main__":
    main()