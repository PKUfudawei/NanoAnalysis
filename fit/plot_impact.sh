text2workspace.py datacard/Run2/1000/SRH.txt

rm -rf *.root;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 --doInitialFit;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 --doFits;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 --output impacts/Run2/impacts_Hbb_1000.json;

rm -rf *.root;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t -1 --expectSignal 0 --robustFit 1 --doInitialFit;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t -1 --expectSignal 0 --robustFit 1 --doFits;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t -1 --expectSignal 0 --robustFit 1 --output impacts/Run2/impacts_Hbb_1000_r0.json;

rm -rf *.root;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t -1 --expectSignal 1 --robustFit 1 --doInitialFit;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t -1 --expectSignal 1 --robustFit 1 --doFits;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t -1 --expectSignal 1 --robustFit 1 --output impacts/Run2/impacts_Hbb_1000_r1.json;

rm -rf *.root;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t 1 --expectSignal 0 --robustFit 1 --doInitialFit;
python3 impact_helper.py
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t 1 --expectSignal 0 --robustFit 1 --doFits;
python3 impact_helper.py
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t 1 --expectSignal 0 --robustFit 1 --output impacts/Run2/impacts_Hbb_1000_t1_r0.json;

rm -rf *.root;
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t 1 --expectSignal 1 --robustFit 1 --doInitialFit;
python3 impact_helper.py
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t 1 --expectSignal 1 --robustFit 1 --doFits;
python3 impact_helper.py
combineTool.py -M Impacts -d datacard/Run2/1000/SRH.root -m 125 -t 1 --expectSignal 1 --robustFit 1 --output impacts/Run2/impacts_Hbb_1000_t1_r1.json;


plotImpacts.py -i impacts/Run2/impacts_Hbb_1000.json -o ../plots/fit/Run2/1000/impacts_Hbb_1000
plotImpacts.py -i impacts/Run2/impacts_Hbb_1000_r0.json -o ../plots/fit/Run2/1000/impacts_Hbb_1000_r0
plotImpacts.py -i impacts/Run2/impacts_Hbb_1000_r1.json -o ../plots/fit/Run2/1000/impacts_Hbb_1000_r1
plotImpacts.py -i impacts/Run2/impacts_Hbb_1000_t1_r0.json -o ../plots/fit/Run2/1000/impacts_Hbb_1000_t1_r0
plotImpacts.py -i impacts/Run2/impacts_Hbb_1000_t1_r1.json -o ../plots/fit/Run2/1000/impacts_Hbb_1000_t1_r1 