#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=18.49m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=11eb51a4
FILE0="davs://dcache-cms-webdav-wan.desy.de:2880/pnfs/desy.de/cms/tier2/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/43DB8DB0-7F8F-114D-B363-DA38742C4414.root"
FILE1="davs://eoscms.cern.ch:443/eos/cms/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/43DB8DB0-7F8F-114D-B363-DA38742C4414.root"
FILE2="davs://maite.iihe.ac.be:2880/pnfs/iihe/cms/ph/sc4/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/43DB8DB0-7F8F-114D-B363-DA38742C4414.root"
FILE3="davs://osg-se.sprace.org.br:1094/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/43DB8DB0-7F8F-114D-B363-DA38742C4414.root"
FILE4="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/43DB8DB0-7F8F-114D-B363-DA38742C4414.root"
FILE5="davs://xrootd.rcac.purdue.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/43DB8DB0-7F8F-114D-B363-DA38742C4414.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZJetsToQQ/nanov9//ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZJetsToQQ/nanov9//ZJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM//43DB8DB0-7F8F-114D-B363-DA38742C4414.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZJetsToQQ/nanov9/ZJetsToQQ_HT-600to800_TuneCP5//ZJetsToQQ_HT-600to800_TuneCP5_34.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZJetsToQQ/nanov9/ZJetsToQQ_HT-600to800_TuneCP5//ZJetsToQQ_HT-600to800_TuneCP5_34.done
