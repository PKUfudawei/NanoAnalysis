#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=1.52m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=ee687931
FILE0="davs://cms-se.sdfarm.kr:1094/xrd/store/mc/RunIIAutumn18NanoAODv7/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/64618189-5454-AB4D-BBDB-C43E75BD241D.root"
FILE1="davs://cms-t2-se01.sdfarm.kr:2880/store/mc/RunIIAutumn18NanoAODv7/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/64618189-5454-AB4D-BBDB-C43E75BD241D.root"
FILE2="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIIAutumn18NanoAODv7/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/64618189-5454-AB4D-BBDB-C43E75BD241D.root"
FILE3="davs://node12.datagrid.cea.fr:26633/dpm/datagrid.cea.fr/home/cms/trivcat/store/mc/RunIIAutumn18NanoAODv7/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/64618189-5454-AB4D-BBDB-C43E75BD241D.root"
FILE4="davs://osg-se.sprace.org.br:1094/store/mc/RunIIAutumn18NanoAODv7/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/64618189-5454-AB4D-BBDB-C43E75BD241D.root"
FILE5="davs://polgrid4.in2p3.fr:443//dpm/in2p3.fr/home/cms/trivcat/store/mc/RunIIAutumn18NanoAODv7/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/64618189-5454-AB4D-BBDB-C43E75BD241D.root"
FILE6="davs://stwebdav.pi.infn.it:8443/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/64618189-5454-AB4D-BBDB-C43E75BD241D.root"
FILE7="srm://cmssrm-kit.gridka.de:8443/srm/managerv2?SFN=/pnfs/gridka.de/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/270000/64618189-5454-AB4D-BBDB-C43E75BD241D.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6 $FILE7; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//nanov7//QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//nanov7//QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//64618189-5454-AB4D-BBDB-C43E75BD241D.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT1000to1500_TuneCP5//QCD_HT1000to1500_TuneCP5_42.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT1000to1500_TuneCP5//QCD_HT1000to1500_TuneCP5_42.done
