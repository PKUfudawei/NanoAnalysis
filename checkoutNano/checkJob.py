import argparse
import os
import glob

parser = argparse.ArgumentParser(description='Script to check if each condor job is done')
parser.add_argument('-d', '--directory', help='To specify base directory', default=os.path.abspath('./jobs'))
parser.add_argument('-t', '--type', help='To specify jobs in mc/ or data/', choices=('data', 'mc'), default='*')
parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
parser.add_argument('-v', '--version', help='To specify jobs in which nanoAOD version', default='*')
args = parser.parse_args()

print("========> Processing...")
jobs = '/'.join([args.directory, args.type, args.year, args.channel, args.version, '*', '*'])
condor_undone = set(glob.glob(jobs+'.jid'))
condor_resubmit = set(i.replace('jid', 'sub') for i in condor_undone)

if condor_undone:
    for i in condor_resubmit:
        os.system(f'condor_submit {i}')
else:
    print("All condor jobs are done!")
