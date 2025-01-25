text2workspace.py datacard/1000/TagHbb.txt
combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -100 --rMax 100 --robustFit 1 -t 1 --expectSignal 0 --doInitialFit;
combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -100 --rMax 100 --robustFit 1 -t 1 --expectSignal 0 --doFits;
combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -100 --rMax 100 --robustFit 1 -t 1 --expectSignal 0 --output impacts_Hbb_1000.json;
plotImpacts.py -i impacts_Hbb_1000.json -o ../plots/fit/Run2/impacts_Hbb_1000