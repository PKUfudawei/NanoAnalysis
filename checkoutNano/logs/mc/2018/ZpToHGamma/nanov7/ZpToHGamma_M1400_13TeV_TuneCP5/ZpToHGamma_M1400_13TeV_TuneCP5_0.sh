#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=7.65m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=52923dd6
FILE0="davs://cmsrm-xrootd01.roma1.infn.it:2880/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1400_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/C4C8923F-DEF0-7341-9033-AA3B150BC17C.root"
FILE1="davs://cmsxrootd.hep.wisc.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1400_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/C4C8923F-DEF0-7341-9033-AA3B150BC17C.root"
FILE2="davs://webdav.echo.stfc.ac.uk:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1400_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/C4C8923F-DEF0-7341-9033-AA3B150BC17C.root"
FILE3="davs://xfer-cms.cr.cnaf.infn.it:8443/cmsdisk/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1400_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/C4C8923F-DEF0-7341-9033-AA3B150BC17C.root"
FILE4="davs://xrootd.hep.kbfi.ee:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1400_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/C4C8923F-DEF0-7341-9033-AA3B150BC17C.root"
FILE5="davs://xrootd.rcac.purdue.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1400_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/C4C8923F-DEF0-7341-9033-AA3B150BC17C.root"
FILE6="srm://srm-cms-mss.jinr-t1.ru:8443/srm/managerv2?SFN=/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1400_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/C4C8923F-DEF0-7341-9033-AA3B150BC17C.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M1400_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M1400_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//C4C8923F-DEF0-7341-9033-AA3B150BC17C.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M1400_13TeV_TuneCP5//ZpToHGamma_M1400_13TeV_TuneCP5_0.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M1400_13TeV_TuneCP5//ZpToHGamma_M1400_13TeV_TuneCP5_0.done
