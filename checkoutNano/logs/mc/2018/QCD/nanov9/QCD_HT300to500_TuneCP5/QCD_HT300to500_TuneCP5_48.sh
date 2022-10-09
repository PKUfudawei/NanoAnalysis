#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=19.01m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=2024f5b4
FILE0="davs://ccsrm.ihep.ac.cn:443/dpm/ihep.ac.cn/home/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/70000/BB5AF3A2-BBFE-7546-A338-991923EE9542.root"
FILE1="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/70000/BB5AF3A2-BBFE-7546-A338-991923EE9542.root"
FILE2="davs://dcache-cms-webdav-wan.desy.de:2880/pnfs/desy.de/cms/tier2/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/70000/BB5AF3A2-BBFE-7546-A338-991923EE9542.root"
FILE3="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/70000/BB5AF3A2-BBFE-7546-A338-991923EE9542.root"
FILE4="davs://webdav.echo.stfc.ac.uk:1094/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/70000/BB5AF3A2-BBFE-7546-A338-991923EE9542.root"
FILE5="davs://xfer-cms.cr.cnaf.infn.it:8443/cmsdisk/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/70000/BB5AF3A2-BBFE-7546-A338-991923EE9542.root"
FILE6="davs://xrootd-redir-stageout.ultralight.org:1095/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/70000/BB5AF3A2-BBFE-7546-A338-991923EE9542.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM//BB5AF3A2-BBFE-7546-A338-991923EE9542.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT300to500_TuneCP5//QCD_HT300to500_TuneCP5_48.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT300to500_TuneCP5//QCD_HT300to500_TuneCP5_48.done
