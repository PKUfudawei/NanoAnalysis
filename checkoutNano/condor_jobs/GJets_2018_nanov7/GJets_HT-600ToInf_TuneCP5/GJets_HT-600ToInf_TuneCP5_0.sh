#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=0.18m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=9b9a6cbb
FILE0="davs://ccdavcms.in2p3.fr:2880/disk/data/store/mc/RunIIAutumn18NanoAODv7/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/36EEB84F-2600-D643-83BA-F69417922C05.root"
FILE1="davs://cmsxrd.ts.infn.it:1094/store/mc/RunIIAutumn18NanoAODv7/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/36EEB84F-2600-D643-83BA-F69417922C05.root"
FILE2="davs://dcache-cms-webdav-wan.desy.de:2880/pnfs/desy.de/cms/tier2/store/mc/RunIIAutumn18NanoAODv7/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/36EEB84F-2600-D643-83BA-F69417922C05.root"
FILE3="davs://eos.grid.vbc.ac.at:8443/eos/vbc/experiments/cms/store/mc/RunIIAutumn18NanoAODv7/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/36EEB84F-2600-D643-83BA-F69417922C05.root"
FILE4="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/36EEB84F-2600-D643-83BA-F69417922C05.root"
FILE5="davs://t2-xrdcms.lnl.infn.it:2880/pnfs/lnl.infn.it/data/cms/store/mc/RunIIAutumn18NanoAODv7/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/36EEB84F-2600-D643-83BA-F69417922C05.root"
FILE6="davs://xrootd.cmsaf.mit.edu:1094/store/mc/RunIIAutumn18NanoAODv7/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/36EEB84F-2600-D643-83BA-F69417922C05.root"
FILE7="davs://xrootd.rcac.purdue.edu:1094/store/mc/RunIIAutumn18NanoAODv7/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/36EEB84F-2600-D643-83BA-F69417922C05.root"
FILE8="srm://cmsdcatape.fnal.gov:8443/srm/managerv2?SFN=/11/store/mc/RunIIAutumn18NanoAODv7/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/100000/36EEB84F-2600-D643-83BA-F69417922C05.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6 $FILE7 $FILE8; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//nanov7//GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//nanov7//GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21_ext1-v1/NANOAODSIM//36EEB84F-2600-D643-83BA-F69417922C05.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/GJets_2018_nanov7/GJets_HT-600ToInf_TuneCP5//GJets_HT-600ToInf_TuneCP5_0.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/GJets_2018_nanov7/GJets_HT-600ToInf_TuneCP5//GJets_HT-600ToInf_TuneCP5_0.done
