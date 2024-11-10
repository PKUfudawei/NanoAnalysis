#!/bin/bash
export HOME=`pwd`;
echo "===> Home Directory:"; echo $HOME
source /cvmfs/cms.cern.ch/cmsset_default.sh

## configure python environment
echo "===> Python3 information: "; which python3; which pip3
echo "===> Python installing/upgrading modules"
python3 -m pip install -U coffea --user

## print environment info
printf "===> Start time: "; /bin/date
printf "===> Job is running on node: "; /bin/hostname
printf "===> Job running as user: "; /usr/bin/id; voms-proxy-info
printf "===> Job is running in directory: "; /bin/pwd

## generate customized NanoAOD file
wget https://raw.githubusercontent.com/PKUfudawei/customizedNanoAOD/master/main.py; mv main.py customized_NanoAOD.py
python3 customized_NanoAOD.py -m $1 -f $3

rm -rf *.parq *.pkl
## execute main.py
python3 src/main.py -m $1 -p $2 -f custom_nano.root

## transfer files via condor
#xrdcp -f *.parquet root://eosuser.cern.ch/$3
