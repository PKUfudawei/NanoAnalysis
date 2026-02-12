echo $1
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git -c advice.detachedHead=false clone --depth 1 --branch v10.2.1 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
scramv1 b clean; scramv1 b # always make a clean build

cd ../../../..


combine -M HybridNew $1
combine -M HybridNew $1 --expectedFromGrid=0.025
combine -M HybridNew $1 --expectedFromGrid=0.16
combine -M HybridNew $1 --expectedFromGrid=0.5
combine -M HybridNew $1 --expectedFromGrid=0.84
combine -M HybridNew $1 --expectedFromGrid=0.975
