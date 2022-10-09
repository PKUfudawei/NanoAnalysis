#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=15.21m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=063d3ca0
FILE0="davs://ccdavcms.in2p3.fr:2880/disk/data/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/5844EC70-38B4-CF48-87C2-DF74AF87A4C2.root"
FILE1="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/5844EC70-38B4-CF48-87C2-DF74AF87A4C2.root"
FILE2="davs://storage01.lcg.cscs.ch:2880//pnfs/lcg.cscs.ch/cms/trivcat/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/5844EC70-38B4-CF48-87C2-DF74AF87A4C2.root"
FILE3="davs://t2-xrdcms.lnl.infn.it:2880/pnfs/lnl.infn.it/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/5844EC70-38B4-CF48-87C2-DF74AF87A4C2.root"
FILE4="davs://xrootd-local.unl.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/5844EC70-38B4-CF48-87C2-DF74AF87A4C2.root"
FILE5="davs://xrootd.cmsaf.mit.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/5844EC70-38B4-CF48-87C2-DF74AF87A4C2.root"
FILE6="srm://srm-cms-mss.jinr-t1.ru:8443/srm/managerv2?SFN=/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/260000/5844EC70-38B4-CF48-87C2-DF74AF87A4C2.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M1200_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M1200_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//5844EC70-38B4-CF48-87C2-DF74AF87A4C2.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M1200_13TeV_TuneCP5//ZpToHGamma_M1200_13TeV_TuneCP5_1.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M1200_13TeV_TuneCP5//ZpToHGamma_M1200_13TeV_TuneCP5_1.done
