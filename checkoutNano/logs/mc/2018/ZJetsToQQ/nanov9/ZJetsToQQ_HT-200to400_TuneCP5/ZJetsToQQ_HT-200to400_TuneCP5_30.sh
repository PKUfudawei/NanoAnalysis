#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=18.89m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=39fe8b6c
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/F47411B5-A65B-3A49-8240-A195321E8429.root"
FILE1="davs://grid-webdav.physik.rwth-aachen.de:2889/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/F47411B5-A65B-3A49-8240-A195321E8429.root"
FILE2="davs://node12.datagrid.cea.fr:26633/dpm/datagrid.cea.fr/home/cms/trivcat/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/F47411B5-A65B-3A49-8240-A195321E8429.root"
FILE3="davs://storage01.lcg.cscs.ch:2880//pnfs/lcg.cscs.ch/cms/trivcat/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/F47411B5-A65B-3A49-8240-A195321E8429.root"
FILE4="davs://xrootd-vanderbilt.sites.opensciencegrid.org:1094/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/F47411B5-A65B-3A49-8240-A195321E8429.root"
FILE5="davs://xrootd.hep.kbfi.ee:1094/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/F47411B5-A65B-3A49-8240-A195321E8429.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZJetsToQQ/nanov9//ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZJetsToQQ/nanov9//ZJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM//F47411B5-A65B-3A49-8240-A195321E8429.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZJetsToQQ/nanov9/ZJetsToQQ_HT-200to400_TuneCP5//ZJetsToQQ_HT-200to400_TuneCP5_30.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZJetsToQQ/nanov9/ZJetsToQQ_HT-200to400_TuneCP5//ZJetsToQQ_HT-200to400_TuneCP5_30.done
