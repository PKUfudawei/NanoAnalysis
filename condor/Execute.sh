#!/bin/bash
export HOME=`pwd`;
echo "===> Home Directory:"; echo $HOME
source /cvmfs/cms.cern.ch/cmsset_default.sh


## print environment info
printf "===> Start time: "; /bin/date
printf "===> Job is running on node: "; /bin/hostname
printf "===> Job running as user: "; /usr/bin/id; voms-proxy-info
printf "===> Job is running in directory: "; /bin/pwd
printf "===> Running on file: "; echo $3


## generate customized NanoAOD file
wget https://raw.githubusercontent.com/PKUfudawei/customizedNanoAOD/master/main.py; mv main.py customized_NanoAOD.py
python3 customized_NanoAOD.py -m $1 -f $3


## configure python environment
echo "===> Python3 information: "; which python3; which pip3
echo "===> Python installing/upgrading modules"
pip3 install coffea==2024.8.0 --user


rm -rf *.parq *.pkl
## execute main.py
python3 src/main.py -m $1 -p $2 -f custom_nano.root -n $5


## transfer files via condor
xrdcp -f *.parquet root://eosuser.cern.ch//eos/user/d/dfu/bbgamma_ntuple/condor/$4
xrdcp -r *.yaml root://eosuser.cern.ch//eos/user/d/dfu/bbgamma_ntuple/condor/$4
