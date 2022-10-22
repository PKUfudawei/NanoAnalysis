#!/usr/bin/env python3
import argparse
import os
import subprocess

def parse_commanline():
    parser = argparse.ArgumentParser(description='Script to check if each condor job is done')
    parser.add_argument('-t', '--type', help='To specify jobs in mc/ or data/', choices=('data', 'mc', '*'), default='*')
    parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
    parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
    parser.add_argument('-j', '--job', help='Specify which job to be submitted', default='*')
    parser.add_argument('-b', '--bump', help='Determine whether to change myschedd', choices=('True', 'False', 'true', 'false'), default='true')
    args = parser.parse_args()
    return args

def main() -> None:
    args = parse_commanline()
    if args.bump.capitalize()=='True':
        print("====> Jump to the best schedd")
        schedd = subprocess.check_output(f"myschedd bump", shell=True, encoding='utf-8').split("'")[1]
    else:
        schedd = ''
    ## myschedd set $schedd
    os.system(f"""
        myschedd show
        echo "###########################################################################"
        echo "====> Status fo current schedd"
        condor_status -schedd {schedd}
        echo "###########################################################################"
        echo "====> Wanna submit condor jobs on this schedd?"
        
        while true          # interactive interface
        do
            read -n 1 -p "Press ENTER to submit condor jobs or CTRL+c to exit> " input
            if [[ -z $input ]];then
                break   
            fi
            echo ""
        done
        echo ""

        for jdl in submit/{args.type}/{args.year}/{args.channel}/{args.job}.submit
            do condor_submit $jdl
        done
        """)
    
if __name__ == "__main__":
    main()
