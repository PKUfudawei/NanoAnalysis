#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=15.9m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=84409e98
FILE0="davs://cmsio7.rc.ufl.edu:1094/store/mc/RunIIAutumn18NanoAODv7/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/D9B3BA20-1AD4-0F4D-A707-1B5B7A661149.root"
FILE1="davs://cmsrm-xrootd01.roma1.infn.it:2880/mc/RunIIAutumn18NanoAODv7/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/D9B3BA20-1AD4-0F4D-A707-1B5B7A661149.root"
FILE2="davs://eymir.grid.metu.edu.tr:443//dpm/grid.metu.edu.tr/home/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/D9B3BA20-1AD4-0F4D-A707-1B5B7A661149.root"
FILE3="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIIAutumn18NanoAODv7/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/D9B3BA20-1AD4-0F4D-A707-1B5B7A661149.root"
FILE4="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/D9B3BA20-1AD4-0F4D-A707-1B5B7A661149.root"
FILE5="davs://storage01.lcg.cscs.ch:2880//pnfs/lcg.cscs.ch/cms/trivcat/store/mc/RunIIAutumn18NanoAODv7/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/D9B3BA20-1AD4-0F4D-A707-1B5B7A661149.root"
FILE6="srm://cmsdcatape.fnal.gov:8443/srm/managerv2?SFN=/11/store/mc/RunIIAutumn18NanoAODv7/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/100000/D9B3BA20-1AD4-0F4D-A707-1B5B7A661149.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//nanov7//QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//nanov7//QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//D9B3BA20-1AD4-0F4D-A707-1B5B7A661149.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT200to300_TuneCP5//QCD_HT200to300_TuneCP5_13.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT200to300_TuneCP5//QCD_HT200to300_TuneCP5_13.done
