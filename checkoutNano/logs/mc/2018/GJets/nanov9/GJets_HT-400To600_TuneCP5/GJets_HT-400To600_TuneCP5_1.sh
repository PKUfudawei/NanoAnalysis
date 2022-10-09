#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=12.5m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=88118a80
FILE0="davs://cmsxrd.ts.infn.it:1094/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/4DEF4AA9-9736-594E-AB77-3F33D317DDBF.root"
FILE1="davs://gfe02.grid.hep.ph.ic.ac.uk:2880/pnfs/hep.ph.ic.ac.uk/data/cms/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/4DEF4AA9-9736-594E-AB77-3F33D317DDBF.root"
FILE2="davs://lcgsedr08.jinr.ru:2880/pnfs/jinr.ru/data/cms/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/4DEF4AA9-9736-594E-AB77-3F33D317DDBF.root"
FILE3="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/4DEF4AA9-9736-594E-AB77-3F33D317DDBF.root"
FILE4="davs://webdav.echo.stfc.ac.uk:1094/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/4DEF4AA9-9736-594E-AB77-3F33D317DDBF.root"
FILE5="davs://xrootd-redir-stageout.ultralight.org:1095/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/4DEF4AA9-9736-594E-AB77-3F33D317DDBF.root"
FILE6="davs://xrootd.phy.bris.ac.uk:1094/xrootd/cms/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2550000/4DEF4AA9-9736-594E-AB77-3F33D317DDBF.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/GJets/nanov9//GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/GJets/nanov9//GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM//4DEF4AA9-9736-594E-AB77-3F33D317DDBF.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/GJets/nanov9/GJets_HT-400To600_TuneCP5//GJets_HT-400To600_TuneCP5_1.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/GJets/nanov9/GJets_HT-400To600_TuneCP5//GJets_HT-400To600_TuneCP5_1.done
