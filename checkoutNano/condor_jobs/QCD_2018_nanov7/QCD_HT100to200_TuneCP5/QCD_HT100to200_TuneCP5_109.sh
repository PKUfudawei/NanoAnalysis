#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=4.98m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=7a51e649
FILE0="davs://cmsxrootd.hep.wisc.edu:1094/store/mc/RunIIAutumn18NanoAODv7/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/A50538B0-84FA-504B-A7D0-915B4D956C8F.root"
FILE1="davs://eoscms.cern.ch:443/eos/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/A50538B0-84FA-504B-A7D0-915B4D956C8F.root"
FILE2="davs://gfe02.grid.hep.ph.ic.ac.uk:2880/pnfs/hep.ph.ic.ac.uk/data/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/A50538B0-84FA-504B-A7D0-915B4D956C8F.root"
FILE3="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/A50538B0-84FA-504B-A7D0-915B4D956C8F.root"
FILE4="davs://srmcms.pic.es:8459/pnfs/pic.es/data/cms/disk/store/mc/RunIIAutumn18NanoAODv7/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/A50538B0-84FA-504B-A7D0-915B4D956C8F.root"
FILE5="davs://xrootd.rcac.purdue.edu:1094/store/mc/RunIIAutumn18NanoAODv7/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/A50538B0-84FA-504B-A7D0-915B4D956C8F.root"
FILE6="srm://cmsdcatape.fnal.gov:8443/srm/managerv2?SFN=/11/store/mc/RunIIAutumn18NanoAODv7/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/A50538B0-84FA-504B-A7D0-915B4D956C8F.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//nanov7//QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//nanov7//QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//A50538B0-84FA-504B-A7D0-915B4D956C8F.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT100to200_TuneCP5//QCD_HT100to200_TuneCP5_109.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT100to200_TuneCP5//QCD_HT100to200_TuneCP5_109.done
