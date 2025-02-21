import os

MASS = [700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 3000, 3500]

for year in ['oldRun2']:
    os.makedirs(f'./AsymptoticLimits/{year}', exist_ok=True)
    for mass in os.listdir(f'datacard/{year}'):
        for file in set(os.listdir(f'datacard/{year}/{mass}')):
            print(f'\ndatacard/{year}/{mass}/{file}:')
            os.system(f"combine -M AsymptoticLimits datacard/{year}/{mass}/{file} -n .{mass}.{file.split('.')[0]} --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2")
            os.system(f"mv *AsymptoticLimits*.root ./AsymptoticLimits/{year}/")
            #os.system(f"combine -M MultiDimFit datacard/Run2/{file} -n f{file.split('.')[0].replace('datacard_', '.')} -t -1 --expectSignal 1")
