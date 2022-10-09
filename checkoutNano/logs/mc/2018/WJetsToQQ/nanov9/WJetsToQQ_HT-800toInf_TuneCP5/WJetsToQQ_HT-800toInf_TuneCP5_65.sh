#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=18.51m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=92982920
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/E6053139-4C0B-7945-BC7D-5581950CA8F8.root"
FILE1="davs://dcache-cms-webdav-wan.desy.de:2880/pnfs/desy.de/cms/tier2/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/E6053139-4C0B-7945-BC7D-5581950CA8F8.root"
FILE2="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/E6053139-4C0B-7945-BC7D-5581950CA8F8.root"
FILE3="davs://webdav.recas.ba.infn.it:8443/cms/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/E6053139-4C0B-7945-BC7D-5581950CA8F8.root"
FILE4="davs://xfer-cms.cr.cnaf.infn.it:8443/cmsdisk/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/E6053139-4C0B-7945-BC7D-5581950CA8F8.root"
FILE5="davs://xrootd.cmsaf.mit.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/E6053139-4C0B-7945-BC7D-5581950CA8F8.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/WJetsToQQ/nanov9//WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/WJetsToQQ/nanov9//WJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM//E6053139-4C0B-7945-BC7D-5581950CA8F8.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/WJetsToQQ/nanov9/WJetsToQQ_HT-800toInf_TuneCP5//WJetsToQQ_HT-800toInf_TuneCP5_65.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/WJetsToQQ/nanov9/WJetsToQQ_HT-800toInf_TuneCP5//WJetsToQQ_HT-800toInf_TuneCP5_65.done
