#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=5.93m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=1301dc08
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/C55CCD8A-261B-CE48-BA6F-8A735A0CB965.root"
FILE1="davs://cmswebdav-kit.gridka.de:2880/pnfs/gridka.de/cms/disk-only/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/C55CCD8A-261B-CE48-BA6F-8A735A0CB965.root"
FILE2="davs://sbgse1.in2p3.fr:443/dpm/in2p3.fr/home/cms/phedex/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/C55CCD8A-261B-CE48-BA6F-8A735A0CB965.root"
FILE3="davs://srm.ciemat.es:2880/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/C55CCD8A-261B-CE48-BA6F-8A735A0CB965.root"
FILE4="davs://xrootd.cmsaf.mit.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/C55CCD8A-261B-CE48-BA6F-8A735A0CB965.root"
FILE5="davs://xrootd.hep.kbfi.ee:1094/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/260000/C55CCD8A-261B-CE48-BA6F-8A735A0CB965.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/WJetsToQQ/nanov9//WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/WJetsToQQ/nanov9//WJetsToQQ_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM//C55CCD8A-261B-CE48-BA6F-8A735A0CB965.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/WJetsToQQ/nanov9/WJetsToQQ_HT-200to400_TuneCP5//WJetsToQQ_HT-200to400_TuneCP5_8.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/WJetsToQQ/nanov9/WJetsToQQ_HT-200to400_TuneCP5//WJetsToQQ_HT-200to400_TuneCP5_8.done
