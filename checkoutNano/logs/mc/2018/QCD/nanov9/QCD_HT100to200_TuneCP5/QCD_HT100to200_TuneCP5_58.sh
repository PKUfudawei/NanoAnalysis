#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=2.24m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=c343fca4
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/60000/E084BDFA-CC6C-0447-A826-E369D1B983E7.root"
FILE1="davs://cmsrm-xrootd01.roma1.infn.it:2880/mc/RunIISummer20UL18NanoAODv9/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/60000/E084BDFA-CC6C-0447-A826-E369D1B983E7.root"
FILE2="davs://grid143.kfki.hu:443/dpm/kfki.hu/home/cms/phedex/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/60000/E084BDFA-CC6C-0447-A826-E369D1B983E7.root"
FILE3="davs://lcgsedr08.jinr.ru:2880/pnfs/jinr.ru/data/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/60000/E084BDFA-CC6C-0447-A826-E369D1B983E7.root"
FILE4="davs://se01.indiacms.res.in:443/dpm/indiacms.res.in/home/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/60000/E084BDFA-CC6C-0447-A826-E369D1B983E7.root"
FILE5="davs://srmcms.pic.es:8459/pnfs/pic.es/data/cms/disk/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/60000/E084BDFA-CC6C-0447-A826-E369D1B983E7.root"
FILE6="davs://xrootd.rcac.purdue.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/60000/E084BDFA-CC6C-0447-A826-E369D1B983E7.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM//E084BDFA-CC6C-0447-A826-E369D1B983E7.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT100to200_TuneCP5//QCD_HT100to200_TuneCP5_58.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT100to200_TuneCP5//QCD_HT100to200_TuneCP5_58.done
