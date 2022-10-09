#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=2.31m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=7248657d
FILE0="davs://ccdavcms.in2p3.fr:2880/disk/data/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/100000/5285B20A-E3CB-1A46-9692-812874E47731.root"
FILE1="davs://cmsio7.rc.ufl.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/100000/5285B20A-E3CB-1A46-9692-812874E47731.root"
FILE2="davs://cmswebdav-kit.gridka.de:2880/pnfs/gridka.de/cms/disk-only/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/100000/5285B20A-E3CB-1A46-9692-812874E47731.root"
FILE3="davs://madhatter.csc.fi:2880/pnfs/csc.fi/data/cms/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/100000/5285B20A-E3CB-1A46-9692-812874E47731.root"
FILE4="davs://xrootd.cmsaf.mit.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/100000/5285B20A-E3CB-1A46-9692-812874E47731.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/WJetsToQQ/nanov9//WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/WJetsToQQ/nanov9//WJetsToQQ_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM//5285B20A-E3CB-1A46-9692-812874E47731.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/WJetsToQQ/nanov9/WJetsToQQ_HT-600to800_TuneCP5//WJetsToQQ_HT-600to800_TuneCP5_6.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/WJetsToQQ/nanov9/WJetsToQQ_HT-600to800_TuneCP5//WJetsToQQ_HT-600to800_TuneCP5_6.done
