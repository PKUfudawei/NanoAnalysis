#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=9.28m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=76c70212
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/E5D24FB3-8395-6B47-9B94-668A0EEADB68.root"
FILE1="davs://cmswebdav-kit.gridka.de:2880/pnfs/gridka.de/cms/disk-only/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/E5D24FB3-8395-6B47-9B94-668A0EEADB68.root"
FILE2="davs://cmsxrd.ts.infn.it:1094/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/E5D24FB3-8395-6B47-9B94-668A0EEADB68.root"
FILE3="davs://dcache-cms-webdav-wan.desy.de:2880/pnfs/desy.de/cms/tier2/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/E5D24FB3-8395-6B47-9B94-668A0EEADB68.root"
FILE4="davs://ingrid-se08.cism.ucl.ac.be:1094/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/E5D24FB3-8395-6B47-9B94-668A0EEADB68.root"
FILE5="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/E5D24FB3-8395-6B47-9B94-668A0EEADB68.root"
FILE6="davs://se01.indiacms.res.in:443/dpm/indiacms.res.in/home/cms/store/mc/RunIISummer20UL18NanoAODv9/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/40000/E5D24FB3-8395-6B47-9B94-668A0EEADB68.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/GJets/nanov9//GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/GJets/nanov9//GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM//E5D24FB3-8395-6B47-9B94-668A0EEADB68.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/GJets/nanov9/GJets_HT-400To600_TuneCP5//GJets_HT-400To600_TuneCP5_16.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/GJets/nanov9/GJets_HT-400To600_TuneCP5//GJets_HT-400To600_TuneCP5_16.done
