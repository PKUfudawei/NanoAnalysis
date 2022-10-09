#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=17.28m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=1236ee2a
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2820000/B7A2ACB5-FE18-A649-8493-EEC140B7177A.root"
FILE1="davs://cmsrm-xrootd01.roma1.infn.it:2880/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2820000/B7A2ACB5-FE18-A649-8493-EEC140B7177A.root"
FILE2="davs://node12.datagrid.cea.fr:26633/dpm/datagrid.cea.fr/home/cms/trivcat/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2820000/B7A2ACB5-FE18-A649-8493-EEC140B7177A.root"
FILE3="davs://osg-se.sprace.org.br:1094/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2820000/B7A2ACB5-FE18-A649-8493-EEC140B7177A.root"
FILE4="davs://redirector.t2.ucsd.edu:1095/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2820000/B7A2ACB5-FE18-A649-8493-EEC140B7177A.root"
FILE5="davs://se01.indiacms.res.in:443/dpm/indiacms.res.in/home/cms/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2820000/B7A2ACB5-FE18-A649-8493-EEC140B7177A.root"
FILE6="davs://xfer-cms.cr.cnaf.infn.it:8443/cmsdisk/store/mc/RunIISummer20UL18NanoAODv9/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/2820000/B7A2ACB5-FE18-A649-8493-EEC140B7177A.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/QCD/nanov9//QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM//B7A2ACB5-FE18-A649-8493-EEC140B7177A.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT1500to2000_TuneCP5//QCD_HT1500to2000_TuneCP5_4.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/QCD/nanov9/QCD_HT1500to2000_TuneCP5//QCD_HT1500to2000_TuneCP5_4.done
