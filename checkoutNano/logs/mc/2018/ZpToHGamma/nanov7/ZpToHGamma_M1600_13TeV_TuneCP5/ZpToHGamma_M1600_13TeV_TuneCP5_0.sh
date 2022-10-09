#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=19.21m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=f217d200
FILE0="davs://cmsrm-xrootd01.roma1.infn.it:2880/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1600_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/C8A1D997-B70D-FD47-8204-C0434ED27337.root"
FILE1="davs://grse001.inr.troitsk.ru:443/dpm/inr.troitsk.ru/home/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1600_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/C8A1D997-B70D-FD47-8204-C0434ED27337.root"
FILE2="davs://mover.pp.rl.ac.uk:2880/pnfs/pp.rl.ac.uk/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1600_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/C8A1D997-B70D-FD47-8204-C0434ED27337.root"
FILE3="davs://xfer-cms.cr.cnaf.infn.it:8443/cmsdisk/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1600_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/C8A1D997-B70D-FD47-8204-C0434ED27337.root"
FILE4="davs://xrootd.cmsaf.mit.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1600_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/C8A1D997-B70D-FD47-8204-C0434ED27337.root"
FILE5="davs://xrootd.rcac.purdue.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1600_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/C8A1D997-B70D-FD47-8204-C0434ED27337.root"
FILE6="srm://storm-fe-cms.cr.cnaf.infn.it:8444/srm/managerv2?SFN=/cmstape/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M1600_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/280000/C8A1D997-B70D-FD47-8204-C0434ED27337.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M1600_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M1600_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//C8A1D997-B70D-FD47-8204-C0434ED27337.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M1600_13TeV_TuneCP5//ZpToHGamma_M1600_13TeV_TuneCP5_0.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M1600_13TeV_TuneCP5//ZpToHGamma_M1600_13TeV_TuneCP5_0.done
