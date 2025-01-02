text2workspace.py datacard/2000/TagHbb.txt
combineTool.py -M Impacts -d datacard/2000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --doInitialFit -t -1;
combineTool.py -M Impacts -d datacard/2000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --doFits -t -1;
combineTool.py -M Impacts -d  datacard/2000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --output impacts_Hbb_2000.json -t -1;
plotImpacts.py -i impacts_Hbb_2000.json -o ../plots/fit/Run2/impacts_Hbb_2000