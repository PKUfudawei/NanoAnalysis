#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=3.89m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=be59e24b
FILE0="davs://ccdavcms.in2p3.fr:2880/disk/data/store/mc/RunIIAutumn18NanoAODv7/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/6FC7FD72-70AE-B74A-8133-CBCC6F7B6ED3.root"
FILE1="davs://cms-se.sdfarm.kr:1094/xrd/store/mc/RunIIAutumn18NanoAODv7/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/6FC7FD72-70AE-B74A-8133-CBCC6F7B6ED3.root"
FILE2="davs://cms-t2-se01.sdfarm.kr:2880/store/mc/RunIIAutumn18NanoAODv7/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/6FC7FD72-70AE-B74A-8133-CBCC6F7B6ED3.root"
FILE3="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIIAutumn18NanoAODv7/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/6FC7FD72-70AE-B74A-8133-CBCC6F7B6ED3.root"
FILE4="davs://maite.iihe.ac.be:2880/pnfs/iihe/cms/ph/sc4/store/mc/RunIIAutumn18NanoAODv7/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/6FC7FD72-70AE-B74A-8133-CBCC6F7B6ED3.root"
FILE5="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIIAutumn18NanoAODv7/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/6FC7FD72-70AE-B74A-8133-CBCC6F7B6ED3.root"
FILE6="davs://storage01.lcg.cscs.ch:2880//pnfs/lcg.cscs.ch/cms/trivcat/store/mc/RunIIAutumn18NanoAODv7/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/6FC7FD72-70AE-B74A-8133-CBCC6F7B6ED3.root"
FILE7="srm://ccsrm.in2p3.fr:8443/srm/managerv2?SFN=/pnfs/in2p3.fr/data/cms/data/store/mc/RunIIAutumn18NanoAODv7/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/6FC7FD72-70AE-B74A-8133-CBCC6F7B6ED3.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6 $FILE7; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//nanov7//QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//nanov7//QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//6FC7FD72-70AE-B74A-8133-CBCC6F7B6ED3.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT700to1000_TuneCP5//QCD_HT700to1000_TuneCP5_21.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT700to1000_TuneCP5//QCD_HT700to1000_TuneCP5_21.done
