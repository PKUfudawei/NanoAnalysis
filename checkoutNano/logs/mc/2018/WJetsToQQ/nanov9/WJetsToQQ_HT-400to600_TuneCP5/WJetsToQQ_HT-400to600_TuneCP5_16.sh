#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=10.47m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=2023960b
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2510000/0C42A20A-8ADA-7149-A533-0C9CD957F443.root"
FILE1="davs://madhatter.csc.fi:2880/pnfs/csc.fi/data/cms/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2510000/0C42A20A-8ADA-7149-A533-0C9CD957F443.root"
FILE2="davs://node12.datagrid.cea.fr:26633/dpm/datagrid.cea.fr/home/cms/trivcat/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2510000/0C42A20A-8ADA-7149-A533-0C9CD957F443.root"
FILE3="davs://srm.ciemat.es:2880/pnfs/ciemat.es/data/cms/prod/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2510000/0C42A20A-8ADA-7149-A533-0C9CD957F443.root"
FILE4="davs://webdav.echo.stfc.ac.uk:1094/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2510000/0C42A20A-8ADA-7149-A533-0C9CD957F443.root"
FILE5="davs://xrootd-redir-stageout.ultralight.org:1095/store/mc/RunIISummer20UL18NanoAODv9/WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2510000/0C42A20A-8ADA-7149-A533-0C9CD957F443.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/WJetsToQQ/nanov9//WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/WJetsToQQ/nanov9//WJetsToQQ_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v2/NANOAODSIM//0C42A20A-8ADA-7149-A533-0C9CD957F443.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/WJetsToQQ/nanov9/WJetsToQQ_HT-400to600_TuneCP5//WJetsToQQ_HT-400to600_TuneCP5_16.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/WJetsToQQ/nanov9/WJetsToQQ_HT-400to600_TuneCP5//WJetsToQQ_HT-400to600_TuneCP5_16.done
