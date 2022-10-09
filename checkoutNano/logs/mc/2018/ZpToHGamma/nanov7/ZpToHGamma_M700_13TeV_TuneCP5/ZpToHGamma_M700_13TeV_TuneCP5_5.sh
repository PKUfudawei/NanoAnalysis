#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=1.9m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=73e262e7
FILE0="davs://cmswebdav-kit.gridka.de:2880/pnfs/gridka.de/cms/disk-only/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/221BB4D4-565A-1D4C-919B-2148E27682D4.root"
FILE1="davs://dcache-cms-webdav-wan.desy.de:2880/pnfs/desy.de/cms/tier2/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/221BB4D4-565A-1D4C-919B-2148E27682D4.root"
FILE2="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/221BB4D4-565A-1D4C-919B-2148E27682D4.root"
FILE3="davs://se.cis.gov.pl:443/dpm/cis.gov.pl/home/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/221BB4D4-565A-1D4C-919B-2148E27682D4.root"
FILE4="davs://xrootd.cmsaf.mit.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/221BB4D4-565A-1D4C-919B-2148E27682D4.root"
FILE5="srm://srm-cms-mss.jinr-t1.ru:8443/srm/managerv2?SFN=/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/221BB4D4-565A-1D4C-919B-2148E27682D4.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M700_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//221BB4D4-565A-1D4C-919B-2148E27682D4.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M700_13TeV_TuneCP5//ZpToHGamma_M700_13TeV_TuneCP5_5.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M700_13TeV_TuneCP5//ZpToHGamma_M700_13TeV_TuneCP5_5.done
