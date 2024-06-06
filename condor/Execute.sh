#!/bin/bash
export HOME=`pwd`
echo "===> Home Directory:"; echo $HOME
source /cvmfs/cms.cern.ch/cmsset_default.sh

## configure python environment
export cmssw=CMSSW_12_5_5
echo "===> Initializing $cmssw"; cmsrel $cmssw; cd $cmssw/src; cmsenv; cd -

pip3 install -U pip -i https://pypi.tuna.tsinghua.edu.cn/simple --user
echo "===> Python3 information: "; which python3; which pip3
echo "===> Python installing/upgrading modules"
for i in {1..4}; do
python3 -m pip install -U coffea -i https://pypi.tuna.tsinghua.edu.cn/simple --user
done

## print environment info
printf "===> Start time: "; /bin/date
printf "===> Job is running on node: "; /bin/hostname
printf "===> Job running as user: "; /usr/bin/id; voms-proxy-info
printf "===> Job is running in directory: "; /bin/pwd

## execute main.py
python3 src/main.py -m $1 -p $2 -f $3

## transfer files via condor
#xrdcp -f *.parquet root://eosuser.cern.ch/$3
