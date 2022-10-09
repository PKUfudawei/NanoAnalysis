#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=1.0m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=ea13e59c
FILE0="davs://cmswebdav-kit.gridka.de:2880/pnfs/gridka.de/cms/disk-only/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/120000/3B3DA099-D02D-5F4F-A678-580F2B744E75.root"
FILE1="davs://gfe02.grid.hep.ph.ic.ac.uk:2880/pnfs/hep.ph.ic.ac.uk/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/120000/3B3DA099-D02D-5F4F-A678-580F2B744E75.root"
FILE2="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/120000/3B3DA099-D02D-5F4F-A678-580F2B744E75.root"
FILE3="davs://srm.ciemat.es:2880/pnfs/ciemat.es/data/cms/prod/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/120000/3B3DA099-D02D-5F4F-A678-580F2B744E75.root"
FILE4="davs://xrootd-local.unl.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/120000/3B3DA099-D02D-5F4F-A678-580F2B744E75.root"
FILE5="davs://xrootd-redir-stageout.ultralight.org:1095/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/120000/3B3DA099-D02D-5F4F-A678-580F2B744E75.root"
FILE6="srm://srm-cms-mss.jinr-t1.ru:8443/srm/managerv2?SFN=/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/120000/3B3DA099-D02D-5F4F-A678-580F2B744E75.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//3B3DA099-D02D-5F4F-A678-580F2B744E75.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M700_13TeV_TuneCP5//ZpToHGamma_M700_13TeV_TuneCP5_1.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M700_13TeV_TuneCP5//ZpToHGamma_M700_13TeV_TuneCP5_1.done
