#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=17.74m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=4a94cffb
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIIAutumn18NanoAODv7/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/FF1F1377-DDE7-204A-B6AC-347DA2CE77BF.root"
FILE1="davs://cmsxrootd.hep.wisc.edu:1094/store/mc/RunIIAutumn18NanoAODv7/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/FF1F1377-DDE7-204A-B6AC-347DA2CE77BF.root"
FILE2="davs://gfe02.grid.hep.ph.ic.ac.uk:2880/pnfs/hep.ph.ic.ac.uk/data/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/FF1F1377-DDE7-204A-B6AC-347DA2CE77BF.root"
FILE3="davs://mover.pp.rl.ac.uk:2880/pnfs/pp.rl.ac.uk/data/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/FF1F1377-DDE7-204A-B6AC-347DA2CE77BF.root"
FILE4="davs://se-wbdv.jinr-t1.ru:2880/pnfs/jinr-t1.ru/data/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/FF1F1377-DDE7-204A-B6AC-347DA2CE77BF.root"
FILE5="davs://t2-xrdcms.lnl.infn.it:2880/pnfs/lnl.infn.it/data/cms/store/mc/RunIIAutumn18NanoAODv7/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/FF1F1377-DDE7-204A-B6AC-347DA2CE77BF.root"
FILE6="root://antares.stfc.ac.uk:1094//eos/antares/prod/cms//store/mc/RunIIAutumn18NanoAODv7/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/60000/FF1F1377-DDE7-204A-B6AC-347DA2CE77BF.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//nanov7//QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//nanov7//QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//FF1F1377-DDE7-204A-B6AC-347DA2CE77BF.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT50to100_TuneCP5//QCD_HT50to100_TuneCP5_6.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/condor_jobs/QCD_2018_nanov7/QCD_HT50to100_TuneCP5//QCD_HT50to100_TuneCP5_6.done
