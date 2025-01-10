text2workspace.py datacard/700/TagHbb.txt
combineTool.py -M Impacts -d datacard/700/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --doInitialFit -t -1 --expectSignal 0;
combineTool.py -M Impacts -d datacard/700/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --doFits -t -1 --expectSignal 0;
combineTool.py -M Impacts -d  datacard/700/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --output impacts_Hbb_700_r0.json -t -1 --expectSignal 0;
plotImpacts.py -i impacts_Hbb_700_r0.json -o ../plots/fit/Run2/impacts_Hbb_700_r0