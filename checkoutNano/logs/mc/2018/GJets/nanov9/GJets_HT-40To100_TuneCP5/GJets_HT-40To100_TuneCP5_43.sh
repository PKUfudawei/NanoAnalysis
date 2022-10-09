#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=8.97m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=0deabe9d
FILE0="davs://ccdavcms.in2p3.fr:2880/disk/data/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/FE0E4BC1-8333-D740-B33D-16E9B4214A98.root"
FILE1="davs://cmswebdav-kit.gridka.de:2880/pnfs/gridka.de/cms/disk-only/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/FE0E4BC1-8333-D740-B33D-16E9B4214A98.root"
FILE2="davs://cmsxrd.ts.infn.it:1094/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/FE0E4BC1-8333-D740-B33D-16E9B4214A98.root"
FILE3="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/FE0E4BC1-8333-D740-B33D-16E9B4214A98.root"
FILE4="davs://se01.indiacms.res.in:443/dpm/indiacms.res.in/home/cms/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/FE0E4BC1-8333-D740-B33D-16E9B4214A98.root"
FILE5="davs://xfer-cms.cr.cnaf.infn.it:8443/cmsdisk/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/FE0E4BC1-8333-D740-B33D-16E9B4214A98.root"
FILE6="davs://xrootd-redir-stageout.ultralight.org:1095/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/FE0E4BC1-8333-D740-B33D-16E9B4214A98.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/GJets/nanov9//GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/GJets/nanov9//GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM//FE0E4BC1-8333-D740-B33D-16E9B4214A98.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/GJets/nanov9/GJets_HT-40To100_TuneCP5//GJets_HT-40To100_TuneCP5_43.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/GJets/nanov9/GJets_HT-40To100_TuneCP5//GJets_HT-40To100_TuneCP5_43.done
