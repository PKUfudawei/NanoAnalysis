#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=18.99m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=d10575a3
FILE0="davs://cmsxrootd.hep.wisc.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/77682A1C-EC7D-764E-AC57-5B6E39F05130.root"
FILE1="davs://eoscms.cern.ch:443/eos/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/77682A1C-EC7D-764E-AC57-5B6E39F05130.root"
FILE2="davs://gfe02.grid.hep.ph.ic.ac.uk:2880/pnfs/hep.ph.ic.ac.uk/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/77682A1C-EC7D-764E-AC57-5B6E39F05130.root"
FILE3="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/77682A1C-EC7D-764E-AC57-5B6E39F05130.root"
FILE4="davs://t2-xrdcms.lnl.infn.it:2880/pnfs/lnl.infn.it/data/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/77682A1C-EC7D-764E-AC57-5B6E39F05130.root"
FILE5="davs://xrootd-local.unl.edu:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/77682A1C-EC7D-764E-AC57-5B6E39F05130.root"
FILE6="srm://ccsrm.in2p3.fr:8443/srm/managerv2?SFN=/pnfs/in2p3.fr/data/cms/data/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/77682A1C-EC7D-764E-AC57-5B6E39F05130.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M2000_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//77682A1C-EC7D-764E-AC57-5B6E39F05130.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M2000_13TeV_TuneCP5//ZpToHGamma_M2000_13TeV_TuneCP5_0.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M2000_13TeV_TuneCP5//ZpToHGamma_M2000_13TeV_TuneCP5_0.done
