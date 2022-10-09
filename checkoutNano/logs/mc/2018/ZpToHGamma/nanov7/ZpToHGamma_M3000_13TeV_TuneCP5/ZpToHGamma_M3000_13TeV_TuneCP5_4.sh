#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=16.93m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=b92f3b2e
FILE0="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M3000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/B5F673E5-827A-4B4A-A677-3DCF8B53FB2E.root"
FILE1="davs://se01.indiacms.res.in:443/dpm/indiacms.res.in/home/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M3000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/B5F673E5-827A-4B4A-A677-3DCF8B53FB2E.root"
FILE2="davs://stwebdav.pi.infn.it:8443/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M3000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/B5F673E5-827A-4B4A-A677-3DCF8B53FB2E.root"
FILE3="davs://webdav.echo.stfc.ac.uk:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M3000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/B5F673E5-827A-4B4A-A677-3DCF8B53FB2E.root"
FILE4="davs://webdav.recas.ba.infn.it:8443/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M3000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/B5F673E5-827A-4B4A-A677-3DCF8B53FB2E.root"
FILE5="davs://xrootd.rcac.purdue.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M3000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/B5F673E5-827A-4B4A-A677-3DCF8B53FB2E.root"
FILE6="srm://storm-fe-cms.cr.cnaf.infn.it:8444/srm/managerv2?SFN=/cmstape/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M3000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/B5F673E5-827A-4B4A-A677-3DCF8B53FB2E.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M3000_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M3000_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//B5F673E5-827A-4B4A-A677-3DCF8B53FB2E.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M3000_13TeV_TuneCP5//ZpToHGamma_M3000_13TeV_TuneCP5_4.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M3000_13TeV_TuneCP5//ZpToHGamma_M3000_13TeV_TuneCP5_4.done
