combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --doInitialFit --expectSignal 0 -t -1;
combineTool.py -M Impacts -d datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --doFits --expectSignal 0 -t -1;
combineTool.py -M Impacts -d  datacard/1000/TagHbb.root -m 125 --rMin -1 --rMax 2 --robustFit 1 --output impacts.json --expectSignal 0 -t -1;
plotImpacts.py -i impacts.json -o impacts