text2workspace.py datacard/1000/TagHbb.txt
combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --doInitialFit;
combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --doFits;
combineTool.py -M Impacts -d  datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --output impacts_Hbb_1000.json;
plotImpacts.py -i impacts_Hbb_1000.json -o ../plots/fit/Run2/impacts_Hbb_1000 --blind