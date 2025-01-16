import os

for mass in os.listdir('datacard/'):
    for file in set(os.listdir(f'datacard/{mass}')):
        if file.endswith('.txt') and 'SR1' in file:
            os.system(f"combineCards.py SR1=datacard/{mass}/{file} SR2=datacard/{mass}/{file.replace('SR1', 'SR2')} > datacard/{mass}/{file.replace('_SR1', '')}")

for mass in os.listdir('datacard/'):
    for file in set(os.listdir(f'datacard/{mass}')):
        if not file.endswith('.txt'):
            continue
        os.system(f"combine -M AsymptoticLimits datacard/{mass}/{file} -n .{mass}.{file.split('.')[0].replace('Tag', '')} -t -1 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2")
        #os.system(f"combine -M MultiDimFit datacard/Run2/{file} -n f{file.split('.')[0].replace('datacard_', '.')} -t -1 --expectSignal 1")
        print(f'Above: {file}')

