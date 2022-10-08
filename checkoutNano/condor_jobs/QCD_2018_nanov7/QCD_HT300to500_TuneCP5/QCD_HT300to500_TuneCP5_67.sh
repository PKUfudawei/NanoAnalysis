#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=13.78m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=3068e1fe
FILE0="davs://cms-se.sdfarm.kr:1094/xrd/store/mc/RunIIAutumn18NanoAODv7/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/8AA88BD5-6620-6F4E-90E2-3E3C31DA995A.root"
FILE1="davs://cms-t2-se01.sdfarm.kr:2880/store/mc/RunIIAutumn18NanoAODv7/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/8AA88BD5-6620-6F4E-90E2-3E3C31DA995A.root"
FILE2="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIIAutumn18NanoAODv7/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/8AA88BD5-6620-6F4E-90E2-3E3C31DA995A.root"
FILE3="davs://maite.iihe.ac.be:2880/pnfs/iihe/cms/ph/sc4/store/mc/RunIIAutumn18NanoAODv7/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/8AA88BD5-6620-6F4E-90E2-3E3C31DA995A.root"
FILE4="davs://polgrid4.in2p3.fr:443//dpm/in2p3.fr/home/cms/trivcat/store/mc/RunIIAutumn18NanoAODv7/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/8AA88BD5-6620-6F4E-90E2-3E3C31DA995A.root"
FILE5="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIIAutumn18NanoAODv7/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/8AA88BD5-6620-6F4E-90E2-3E3C31DA995A.root"
FILE6="davs://sbgse1.in2p3.fr:443/dpm/in2p3.fr/home/cms/phedex/store/mc/RunIIAutumn18NanoAODv7/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/8AA88BD5-6620-6F4E-90E2-3E3C31DA995A.root"
FILE7="srm://cmsdcatape.fnal.gov:8443/srm/managerv2?SFN=/11/store/mc/RunIIAutumn18NanoAODv7/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/8AA88BD5-6620-6F4E-90E2-3E3C31DA995A.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6 $FILE7; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//nanov7//QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//nanov7//QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//8AA88BD5-6620-6F4E-90E2-3E3C31DA995A.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT300to500_TuneCP5//QCD_HT300to500_TuneCP5_67.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT300to500_TuneCP5//QCD_HT300to500_TuneCP5_67.done
