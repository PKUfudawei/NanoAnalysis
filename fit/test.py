text2workspace.py datacard/3000/TagHbb.txt
combineTool.py -M Impacts -d datacard/3000/TagHbb.root -m 125 --rMin -1 --rMax 10 --robustFit 1  --doInitialFit;
python3 work.py;
combineTool.py -M Impacts -d datacard/3000/TagHbb.root -m 125 --rMin -1 --rMax 10 --robustFit 1  --doFits;
python3 work.py;
combineTool.py -M Impacts -d datacard/3000/TagHbb.root -m 125 --rMin -1 --rMax 10 --robustFit 1  --output impacts_Hbb_3000.json;
python3 work.py;
plotImpacts.py -i impacts_Hbb_3000.json -o ../plots/fit/Run2/impacts_Hbb_3000