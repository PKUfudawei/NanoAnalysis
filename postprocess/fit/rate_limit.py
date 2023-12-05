import os

for file in set(os.listdir('datacard/Run2/')):
    if file.endswith('.txt') and 'SR1' in file:
        mass = int(file.split('_')[1])
        os.system(f"combineCards.py SR1=datacard/Run2/{file} SR2=datacard/Run2/{file.replace('SR1', 'SR2')} > datacard/Run2/{file.replace('SR1', 'combine')}")

for file in os.listdir('datacard/Run2/'):
        os.system(f"combine -M AsymptoticLimits datacard/Run2/{file} -n f{file.split('.')[0].replace('datacard_', '.')}")
        # --setParameters pdfindex_Tag0=0 --freezeParameters pdfindex_Tag0
