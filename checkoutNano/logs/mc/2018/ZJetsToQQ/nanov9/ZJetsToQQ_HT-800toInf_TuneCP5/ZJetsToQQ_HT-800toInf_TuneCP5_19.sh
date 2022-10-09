#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=10.51m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=5888c86b
FILE0="davs://ccdavcms.in2p3.fr:2880/disk/data/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/B04F2FAF-6C01-394E-B407-A4C3D4632FD5.root"
FILE1="davs://ccsrm.ihep.ac.cn:443/dpm/ihep.ac.cn/home/cms/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/B04F2FAF-6C01-394E-B407-A4C3D4632FD5.root"
FILE2="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/B04F2FAF-6C01-394E-B407-A4C3D4632FD5.root"
FILE3="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/B04F2FAF-6C01-394E-B407-A4C3D4632FD5.root"
FILE4="davs://srmcms.pic.es:8459/pnfs/pic.es/data/cms/disk/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/B04F2FAF-6C01-394E-B407-A4C3D4632FD5.root"
FILE5="davs://storage01.lcg.cscs.ch:2880//pnfs/lcg.cscs.ch/cms/trivcat/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/B04F2FAF-6C01-394E-B407-A4C3D4632FD5.root"
FILE6="davs://xrootd.cmsaf.mit.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/30000/B04F2FAF-6C01-394E-B407-A4C3D4632FD5.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZJetsToQQ/nanov9//ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZJetsToQQ/nanov9//ZJetsToQQ_HT-800toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM//B04F2FAF-6C01-394E-B407-A4C3D4632FD5.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZJetsToQQ/nanov9/ZJetsToQQ_HT-800toInf_TuneCP5//ZJetsToQQ_HT-800toInf_TuneCP5_19.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZJetsToQQ/nanov9/ZJetsToQQ_HT-800toInf_TuneCP5//ZJetsToQQ_HT-800toInf_TuneCP5_19.done
