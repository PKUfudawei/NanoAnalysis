#!/bin/bash
## configure python environment
echo "===> Initializing CMSSW_11_3_4"
cmsrel CMSSW_11_3_4
cd CMSSW_11_3_4/src
cmsenv
echo "===> Python3 version: `python3 -V`"
printf "===> Back to "; cd -
echo "===> Python installing/upgrading modules"
pip3 install -U coffea uproot awkward numpy numba

## print environment info
printf "===> Start time: "; /bin/date 
printf "===> Job is running on node: "; /bin/hostname 
printf "===> Job running as user: "; /usr/bin/id 
voms-proxy-info
printf "===> Job is running in directory: "; /bin/pwd

## execute main.py
echo "`python3 src/main.py -f $1 -m $2`"

## transfer files via condor
#xrdcp -f *.parquet root://eosuser.cern.ch/$3