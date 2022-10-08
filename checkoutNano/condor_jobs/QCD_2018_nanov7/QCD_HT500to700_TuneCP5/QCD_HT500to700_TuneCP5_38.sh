#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=10.29m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=59bdf881
FILE0="davs://cms-se.sdfarm.kr:1094/xrd/store/mc/RunIIAutumn18NanoAODv7/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/70000/DA7FF71E-32BF-E84E-B0C4-5ED7F0120A96.root"
FILE1="davs://cms-t2-se01.sdfarm.kr:2880/store/mc/RunIIAutumn18NanoAODv7/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/70000/DA7FF71E-32BF-E84E-B0C4-5ED7F0120A96.root"
FILE2="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIIAutumn18NanoAODv7/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/70000/DA7FF71E-32BF-E84E-B0C4-5ED7F0120A96.root"
FILE3="davs://cmswebdav-kit.gridka.de:2880/pnfs/gridka.de/cms/disk-only/store/mc/RunIIAutumn18NanoAODv7/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/70000/DA7FF71E-32BF-E84E-B0C4-5ED7F0120A96.root"
FILE4="davs://maite.iihe.ac.be:2880/pnfs/iihe/cms/ph/sc4/store/mc/RunIIAutumn18NanoAODv7/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/70000/DA7FF71E-32BF-E84E-B0C4-5ED7F0120A96.root"
FILE5="davs://webdav.echo.stfc.ac.uk:1094/store/mc/RunIIAutumn18NanoAODv7/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/70000/DA7FF71E-32BF-E84E-B0C4-5ED7F0120A96.root"
FILE6="davs://xrootd.cmsaf.mit.edu:1094/store/mc/RunIIAutumn18NanoAODv7/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/70000/DA7FF71E-32BF-E84E-B0C4-5ED7F0120A96.root"
FILE7="root://antares.stfc.ac.uk:1094//eos/antares/prod/cms//store/mc/RunIIAutumn18NanoAODv7/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/70000/DA7FF71E-32BF-E84E-B0C4-5ED7F0120A96.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6 $FILE7; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//nanov7//QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//nanov7//QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//DA7FF71E-32BF-E84E-B0C4-5ED7F0120A96.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT500to700_TuneCP5//QCD_HT500to700_TuneCP5_38.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT500to700_TuneCP5//QCD_HT500to700_TuneCP5_38.done
