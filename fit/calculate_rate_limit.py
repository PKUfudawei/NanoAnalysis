import os

for mass in os.listdir('datacard/Run2'):
    for file in set(os.listdir(f'datacard/Run2/{mass}')):
        os.system(f"combine -M AsymptoticLimits datacard/Run2/{mass}/{file} -n .{mass}.{file.split('.')[0]} --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2")
        #os.system(f"combine -M MultiDimFit datacard/Run2/{file} -n f{file.split('.')[0].replace('datacard_', '.')} -t -1 --expectSignal 1")
