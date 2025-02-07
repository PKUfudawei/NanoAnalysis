text2workspace.py datacard/2016/1000/TagHbb.txt

rm -rf *.root;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 100 --robustFit 1 --doInitialFit;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 100 --robustFit 1 --doFits;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 100 --robustFit 1 --output impacts/2016/impacts_Hbb_1000.json;

rm -rf *.root;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 1 -t -1 --expectSignal 0 --robustFit 1 --doInitialFit;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 1 -t -1 --expectSignal 0 --robustFit 1 --doFits;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 1 -t -1 --expectSignal 0 --robustFit 1 --output impacts/2016/impacts_Hbb_1000_r0.json;

rm -rf *.root;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 -t -1 --expectSignal 1 --robustFit 1 --doInitialFit;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 -t -1 --expectSignal 1 --robustFit 1 --doFits;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 -t -1 --expectSignal 1 --robustFit 1 --output impacts/2016/impacts_Hbb_1000_r1.json;

rm -rf *.root;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 1 -t 1 --expectSignal 0 --robustFit 1 --doInitialFit;
python3 work.py
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 1 -t 1 --expectSignal 0 --robustFit 1 --doFits;
python3 work.py
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 1 -t 1 --expectSignal 0 --robustFit 1 --output impacts/2016/impacts_Hbb_1000_t1_r0.json;

rm -rf *.root;
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 -t 1 --expectSignal 1 --robustFit 1 --doInitialFit;
python3 work.py
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 -t 1 --expectSignal 1 --robustFit 1 --doFits;
python3 work.py
combineTool.py -M Impacts -d datacard/2016/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 -t 1 --expectSignal 1 --robustFit 1 --output impacts/2016/impacts_Hbb_1000_t1_r1.json;


plotImpacts.py -i impacts/2016/impacts_Hbb_1000.json -o ../plots/fit/2016/impacts_Hbb_1000 --blind
plotImpacts.py -i impacts/2016/impacts_Hbb_1000_r0.json -o ../plots/fit/2016/impacts_Hbb_1000_r0
plotImpacts.py -i impacts/2016/impacts_Hbb_1000_r1.json -o ../plots/fit/2016/impacts_Hbb_1000_r1
plotImpacts.py -i impacts/2016/impacts_Hbb_1000_t1_r0.json -o ../plots/fit/2016/impacts_Hbb_1000_t1_r0
plotImpacts.py -i impacts/2016/impacts_Hbb_1000_t1_r1.json -o ../plots/fit/2016/impacts_Hbb_1000_t1_r1 