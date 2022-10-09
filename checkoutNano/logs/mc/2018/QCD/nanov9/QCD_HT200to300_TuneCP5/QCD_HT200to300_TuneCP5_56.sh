#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=17.12m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=4095274c
FILE0="davs://cmsxrootd.hep.wisc.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/EF655A34-A51D-7041-8A41-3081A009E5A2.root"
FILE1="davs://lcgsedr08.jinr.ru:2880/pnfs/jinr.ru/data/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/EF655A34-A51D-7041-8A41-3081A009E5A2.root"
FILE2="davs://mover.pp.rl.ac.uk:2880/pnfs/pp.rl.ac.uk/data/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/EF655A34-A51D-7041-8A41-3081A009E5A2.root"
FILE3="davs://se01.indiacms.res.in:443/dpm/indiacms.res.in/home/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/EF655A34-A51D-7041-8A41-3081A009E5A2.root"
FILE4="davs://t2-xrdcms.lnl.infn.it:2880/pnfs/lnl.infn.it/data/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/EF655A34-A51D-7041-8A41-3081A009E5A2.root"
FILE5="davs://webdav.echo.stfc.ac.uk:1094/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/EF655A34-A51D-7041-8A41-3081A009E5A2.root"
FILE6="davs://xrootd-local.unl.edu:1094/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2830000/EF655A34-A51D-7041-8A41-3081A009E5A2.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM//EF655A34-A51D-7041-8A41-3081A009E5A2.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT200to300_TuneCP5//QCD_HT200to300_TuneCP5_56.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT200to300_TuneCP5//QCD_HT200to300_TuneCP5_56.done
