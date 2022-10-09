#!/bin/bash
export X509_USER_PROXY=/home/pku/fudawei/.proxy
voms-proxy-info
TO_SLEEP=16.27m
echo "===> Sleep How Long: ${TO_SLEEP}in"
sleep ${TO_SLEEP}
# Load data
echo '===> PWD: '$PWD
export WORK_PATH=$PWD
export FILECHECKSUM=0aa0558b
FILE0="davs://cmsdcadisk.fnal.gov:2880/dcache/uscmsdisk/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/69A54DEC-2FF5-7648-B0FA-8DD60CE49BDE.root"
FILE1="davs://cmsrm-xrootd01.roma1.infn.it:2880/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/69A54DEC-2FF5-7648-B0FA-8DD60CE49BDE.root"
FILE2="davs://grse001.inr.troitsk.ru:443/dpm/inr.troitsk.ru/home/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/69A54DEC-2FF5-7648-B0FA-8DD60CE49BDE.root"
FILE3="davs://node12.datagrid.cea.fr:26633/dpm/datagrid.cea.fr/home/cms/trivcat/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/69A54DEC-2FF5-7648-B0FA-8DD60CE49BDE.root"
FILE4="davs://osg-se.sprace.org.br:1094/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/69A54DEC-2FF5-7648-B0FA-8DD60CE49BDE.root"
FILE5="davs://stwebdav.pi.infn.it:8443/cms/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/69A54DEC-2FF5-7648-B0FA-8DD60CE49BDE.root"
FILE6="srm://storm-fe-cms.cr.cnaf.infn.it:8444/srm/managerv2?SFN=/cmstape/store/mc/RunIIAutumn18NanoAODv7/ZpToHGamma_M2200_13TeV_TuneCP5_madgraph-pythia8/NANOAODSIM/Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/110000/69A54DEC-2FF5-7648-B0FA-8DD60CE49BDE.root"
for ifile in $FILE0 $FILE1 $FILE2 $FILE3 $FILE4 $FILE5 $FILE6; do
    sleep 2s
    export X509_USER_PROXY=/home/pku/fudawei/.proxy
    gfal-copy -t 86400 -T 86400 -f $ifile /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M2200_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM/
    if [ $? -eq 0 ]; then break; fi
done

# Checksum
check_file=`gfal-sum /data/pubfs/fudawei/samples//mc/2018/ZpToHGamma/nanov7//ZpToHGamma_M2200_13TeV_TuneCP5_madgraph-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_102X_upgrade2018_realistic_v21-v1/NANOAODSIM//69A54DEC-2FF5-7648-B0FA-8DD60CE49BDE.root ADLER32`
[ ${check_file#* } = ${FILECHECKSUM} ]
[ $? -eq 0 ] && mv /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M2200_13TeV_TuneCP5//ZpToHGamma_M2200_13TeV_TuneCP5_0.jid /home/pku/fudawei/GitRepo/NanoAnalysis/checkoutNano/logs/mc/2018/ZpToHGamma/nanov7/ZpToHGamma_M2200_13TeV_TuneCP5//ZpToHGamma_M2200_13TeV_TuneCP5_0.done
