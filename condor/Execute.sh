#!/bin/bash
export HOME=`pwd`
echo "===> Home Directory:"; echo $HOME
source /cvmfs/cms.cern.ch/cmsset_default.sh

## configure python environment
echo "===> Initializing CMSSW_12_5_5"; cmsrel CMSSW_12_5_5; cd CMSSW_12_5_5/src; cmsenv; cd -
echo "===> Python3 information: "; which python3; which pip3
echo "===> Python installing/upgrading modules"
pip3 install coffea

## print environment info
printf "===> Start time: "; /bin/date
printf "===> Job is running on node: "; /bin/hostname
printf "===> Job running as user: "; /usr/bin/id; voms-proxy-info
printf "===> Job is running in directory: "; /bin/pwd

## execute main.py
python3 src/main.py -m $1 -p $2

## transfer files via condor
#xrdcp -f *.parquet root://eosuser.cern.ch/$3