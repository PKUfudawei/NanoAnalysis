import os

for mass in os.listdir(f'datacard/Run2'):
    for year in ['Run2']:
        prefix = f'datacard/{year}/{mass}'

        for file in set(os.listdir(f'datacard/{year}/{mass}')):
            if file.endswith('.txt') and 'SR1' in file:
                os.system(f"cd {prefix}; combineCards.py SR1={file} SR2={file.replace('SR1', 'SR2')} > {file.replace('_SR1', '')}")
    if '_' in mass:
        continue

    os.makedirs(f'datacard/Run2/{mass}', exist_ok=True)
    os.system(f"combineCards.py year2016=datacard/2016/{mass}/TagHbb.txt year2017=datacard/2017/{mass}/TagHbb.txt year2018=datacard/2018/{mass}/TagHbb.txt > datacard/Run2/{mass}/TagHbb.txt")
"""
for mass in os.listdir('datacard/'):
    for file in set(os.listdir(f'datacard/{mass}')):
        if not file.endswith('.txt'):
            continue
        os.system(f"combine -M AsymptoticLimits datacard/{mass}/{file} -n .{mass}.{file.split('.')[0].replace('Tag', '')} --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2")
        #os.system(f"combine -M MultiDimFit datacard/Run2/{file} -n f{file.split('.')[0].replace('datacard_', '.')} -t -1 --expectSignal 1")
        print(f'Above: {file}')
"""
