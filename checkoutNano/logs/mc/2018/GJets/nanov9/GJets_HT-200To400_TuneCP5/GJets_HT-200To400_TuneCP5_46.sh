#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=6.02m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=a0bd0335
FILE0="davs://ccdavcms.in2p3.fr:2880/disk/data/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/C829FA59-0C0C-8044-9E98-D53E74EABF1F.root"
FILE1="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/C829FA59-0C0C-8044-9E98-D53E74EABF1F.root"
FILE2="davs://cmsxrd.ts.infn.it:1094/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/C829FA59-0C0C-8044-9E98-D53E74EABF1F.root"
FILE3="davs://node12.datagrid.cea.fr:26633/dpm/datagrid.cea.fr/home/cms/trivcat/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/C829FA59-0C0C-8044-9E98-D53E74EABF1F.root"
FILE4="davs://se01.indiacms.res.in:443/dpm/indiacms.res.in/home/cms/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/C829FA59-0C0C-8044-9E98-D53E74EABF1F.root"
FILE5="davs://srmcms.pic.es:8459/pnfs/pic.es/data/cms/disk/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/C829FA59-0C0C-8044-9E98-D53E74EABF1F.root"
FILE6="davs://xrootd-redir-stageout.ultralight.org:1095/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/C829FA59-0C0C-8044-9E98-D53E74EABF1F.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/GJets/nanov9//GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/GJets/nanov9//GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM//C829FA59-0C0C-8044-9E98-D53E74EABF1F.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/GJets/nanov9/GJets_HT-200To400_TuneCP5//GJets_HT-200To400_TuneCP5_46.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/GJets/nanov9/GJets_HT-200To400_TuneCP5//GJets_HT-200To400_TuneCP5_46.done
