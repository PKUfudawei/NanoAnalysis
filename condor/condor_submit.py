#!/usr/bin/env python3
import argparse, os, glob


def parse_commandline():
    parser = argparse.ArgumentParser(description='Script to check if each condor job is done')
    parser.add_argument('-t', '--type', help='To specify jobs in mc/ or data/', choices=('data', 'mc', '*'), default='*')
    parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
    parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
    parser.add_argument('-j', '--job', help='Specify which job to be submitted', default='*')
    parser.add_argument('-b', '--bump', help='Determine whether to change myschedd', choices=('True', 'False', 'true', 'false'), default='True')
    parser.add_argument('-l', '--lxplus', help='If in lxplus environment', choices=('True', 'False', 'true', 'false'), default='True')
    args = parser.parse_args()
    return args


def main():
    args = parse_commandline()
    if args.bump.capitalize() == 'True':
        print('#' * 150)
        os.system('condor_status -schedd')
        print('#' * 150)
        number = input("====> Which schedd you wanna use? ")
        if number.isdigit() and args.lxplus.capitalize()=='True':
            print(f"====> myschedd set bigbird{number}.cern.ch")
            os.system(f"myschedd set bigbird{number}.cern.ch")
        elif args.lxplus.capitalize()=='True':
            os.system("myschedd show")

    os.system("which condor_submit")
    input("====> Are you ready to submit condor jobs?")

    for jdl in glob.glob(f'submit/{args.type}/{args.year}/{args.channel}/{args.job}.submit'):
        log = jdl.replace('.submit', '/*').replace('submit', 'log')
        os.system(f"rm -rf {log}")
        os.system(f"condor_submit {jdl}")


if __name__ == "__main__":
    main()
