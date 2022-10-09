#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=10.69m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=d7f29781
FILE0="davs://cms-t2-se01.sdfarm.kr:2880/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/DC5CA186-A9D4-1340-9490-C96172A1723A.root"
FILE1="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/DC5CA186-A9D4-1340-9490-C96172A1723A.root"
FILE2="davs://cmswebdav-kit.gridka.de:2880/pnfs/gridka.de/cms/disk-only/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/DC5CA186-A9D4-1340-9490-C96172A1723A.root"
FILE3="davs://gfe02.grid.hep.ph.ic.ac.uk:2880/pnfs/hep.ph.ic.ac.uk/data/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/DC5CA186-A9D4-1340-9490-C96172A1723A.root"
FILE4="davs://se01.indiacms.res.in:443/dpm/indiacms.res.in/home/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/DC5CA186-A9D4-1340-9490-C96172A1723A.root"
FILE5="davs://xfer-cms.cr.cnaf.infn.it:8443/cmsdisk/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/DC5CA186-A9D4-1340-9490-C96172A1723A.root"
FILE6="davs://xrootd.cmsaf.mit.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/DC5CA186-A9D4-1340-9490-C96172A1723A.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM//DC5CA186-A9D4-1340-9490-C96172A1723A.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT200to300_TuneCP5//QCD_HT200to300_TuneCP5_49.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT200to300_TuneCP5//QCD_HT200to300_TuneCP5_49.done
