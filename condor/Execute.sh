#!/bin/bash
## configure python environment
echo "===> Initializing CMSSW_11_3_4"
cmsrel CMSSW_12_2_4
cd CMSSW_12_2_4/src
cmsenv
echo "===> Python3 version: `python3 -V`"
printf "===> Back to working directory"; cd -
echo "===> Python installing/upgrading modules"
pip3 install coffea qiskit

## print environment info
printf "===> Start time: "; /bin/date
printf "===> Job is running on node: "; /bin/hostname
printf "===> Job running as user: "; /usr/bin/id
voms-proxy-info
printf "===> Job is running in directory: "; /bin/pwd

## execute main.py
echo "`python3 src/main.py -f $1 -m $2 -j $3`"

## transfer files via condor
#xrdcp -f *.parquet root://eosuser.cern.ch/$3
