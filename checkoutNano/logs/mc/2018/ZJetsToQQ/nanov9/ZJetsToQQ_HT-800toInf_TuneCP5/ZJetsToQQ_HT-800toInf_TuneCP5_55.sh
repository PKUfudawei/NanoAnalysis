#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=19.09m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=e887d8e0
FILE0="davs://ccdavcms.in2p3.fr:2880/disk/data/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/F11DC48D-05E4-AC4A-AA68-38DE7F7C85A3.root"
FILE1="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/F11DC48D-05E4-AC4A-AA68-38DE7F7C85A3.root"
FILE2="davs://webdav.echo.stfc.ac.uk:1094/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/F11DC48D-05E4-AC4A-AA68-38DE7F7C85A3.root"
FILE3="davs://xfer-cms.cr.cnaf.infn.it:8443/cmsdisk/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/F11DC48D-05E4-AC4A-AA68-38DE7F7C85A3.root"
FILE4="davs://xrootd-local.unl.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/F11DC48D-05E4-AC4A-AA68-38DE7F7C85A3.root"
FILE5="davs://xrootd.rcac.purdue.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/F11DC48D-05E4-AC4A-AA68-38DE7F7C85A3.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZJetsToQQ/nanov9//ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZJetsToQQ/nanov9//ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM//F11DC48D-05E4-AC4A-AA68-38DE7F7C85A3.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZJetsToQQ/nanov9/ZJetsToQQ_HT-800toInf_TuneCP5//ZJetsToQQ_HT-800toInf_TuneCP5_55.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZJetsToQQ/nanov9/ZJetsToQQ_HT-800toInf_TuneCP5//ZJetsToQQ_HT-800toInf_TuneCP5_55.done
