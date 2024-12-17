combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --doInitialFit --expectSignal 0 -t -1;
combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --doFits --expectSignal 0 -t -1;
combineTool.py -M Impacts -d  datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --output impacts_Hbb_1000_r0.json --expectSignal 0 -t -1;
plotImpacts.py -i impacts_Hbb_1000_r0.json -o ../plots/fit/Run2/impacts_Hbb_1000_r0