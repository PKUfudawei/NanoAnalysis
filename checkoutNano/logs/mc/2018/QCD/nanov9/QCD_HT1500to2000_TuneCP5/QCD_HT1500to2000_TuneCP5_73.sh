#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=11.28m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=e06ea3ad
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/80000/A1CBD045-6B7B-3145-9BF9-CE227C2F2B8E.root"
FILE1="davs://grid-webdav.physik.rwth-aachen.de:2889/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/80000/A1CBD045-6B7B-3145-9BF9-CE227C2F2B8E.root"
FILE2="davs://polgrid4.in2p3.fr:443//dpm/in2p3.fr/home/cms/trivcat/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/80000/A1CBD045-6B7B-3145-9BF9-CE227C2F2B8E.root"
FILE3="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/80000/A1CBD045-6B7B-3145-9BF9-CE227C2F2B8E.root"
FILE4="davs://se01.indiacms.res.in:443/dpm/indiacms.res.in/home/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/80000/A1CBD045-6B7B-3145-9BF9-CE227C2F2B8E.root"
FILE5="davs://srmcms.pic.es:8459/pnfs/pic.es/data/cms/disk/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/80000/A1CBD045-6B7B-3145-9BF9-CE227C2F2B8E.root"
FILE6="davs://webdav.recas.ba.infn.it:8443/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/80000/A1CBD045-6B7B-3145-9BF9-CE227C2F2B8E.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM//A1CBD045-6B7B-3145-9BF9-CE227C2F2B8E.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT1500to2000_TuneCP5//QCD_HT1500to2000_TuneCP5_73.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT1500to2000_TuneCP5//QCD_HT1500to2000_TuneCP5_73.done
