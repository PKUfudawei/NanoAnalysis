text2workspace.py datacard/1000/TagHbb.txt
combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 1 --robustFit 1 -t 1 --expectSignal 0 --doInitialFit;
python3 work.py;
combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 1 --robustFit 1 -t 1 --expectSignal 0 --doFits;
python3 work.py;
combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 1 --robustFit 1 -t 1 --expectSignal 0 --output impacts_Hbb_2000_r0_t1.json;
python3 work.py;
plotImpacts.py -i impacts_Hbb_2000_r0_t1.json -o ../plots/fit/Run2/impacts_Hbb_2000_r0_t1