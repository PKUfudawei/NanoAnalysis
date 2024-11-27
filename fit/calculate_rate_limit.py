import os

for mass in os.listdir('datacard/Run2/'):
    for file in set(os.listdir(f'datacard/Run2/{mass}')):
        if file.endswith('.txt') and 'SR1' in file:
            os.system(f"combineCards.py SR1=datacard/Run2/{file} SR2=datacard/Run2/{file.replace('SR1', 'SR2')} > datacard/Run2/{file.replace('_SR1', '')}")

for mass in os.listdir('datacard/Run2/'):
    for file in set(os.listdir(f'datacard/Run2/{mass}')):
        if not file.endswith('.txt'):
            continue
        os.system(f"combine -M AsymptoticLimits datacard/Run2/{mass}/{file} -n {mass}_{file.split('.')[0].replace('Tag', '')} --run blind --setParameters pdfindex_CR1=0,pdfindex_CR2=0 --freezeParameters pdfindex_CR1,pdfindex_CR2")
        #os.system(f"combine -M MultiDimFit datacard/Run2/{file} -n f{file.split('.')[0].replace('datacard_', '.')} -t -1 --expectSignal 1")
        print(f'Above: {file}')
