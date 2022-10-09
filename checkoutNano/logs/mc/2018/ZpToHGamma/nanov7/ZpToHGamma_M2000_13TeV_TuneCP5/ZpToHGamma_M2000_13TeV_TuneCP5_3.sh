#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=3.85m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=7191cb2d
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/ED6F538F-5A06-DD44-B702-8D0BAC9C2321.root"
FILE1="davs://cmswebdav-kit.gridka.de:2880/pnfs/gridka.de/cms/disk-only/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/ED6F538F-5A06-DD44-B702-8D0BAC9C2321.root"
FILE2="davs://ingrid-se08.cism.ucl.ac.be:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/ED6F538F-5A06-DD44-B702-8D0BAC9C2321.root"
FILE3="davs://storage01.lcg.cscs.ch:2880//pnfs/lcg.cscs.ch/cms/trivcat/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/ED6F538F-5A06-DD44-B702-8D0BAC9C2321.root"
FILE4="davs://xrootd.hep.kbfi.ee:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/ED6F538F-5A06-DD44-B702-8D0BAC9C2321.root"
FILE5="davs://xrootd.rcac.purdue.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/ED6F538F-5A06-DD44-B702-8D0BAC9C2321.root"
FILE6="root://antares.stfc.ac.uk:1094//eos/antares/prod/cms//store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/ED6F538F-5A06-DD44-B702-8D0BAC9C2321.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//ED6F538F-5A06-DD44-B702-8D0BAC9C2321.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M2000_13TeV_TuneCP5//ZpToHGamma_M2000_13TeV_TuneCP5_3.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M2000_13TeV_TuneCP5//ZpToHGamma_M2000_13TeV_TuneCP5_3.done
