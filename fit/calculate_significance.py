import os

MASS = [700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 3000, 3500]

os.makedirs(f'./Significance/expected', exist_ok=True)
os.makedirs(f'./Significance/observed', exist_ok=True)

for year in ['Run2']:
    for mass in os.listdir(f'datacard/{year}'):
        for file in set(os.listdir(f'datacard/{year}/{mass}')):
            print(f'\ndatacard/{year}/{mass}/{file}:')
            if 'SRH' in file:
                rMax = 2
            else:
                rMax = 100
            
            os.system(f"combine -M Significance datacard/{year}/{mass}/{file} -n .{mass}.{file.split('.')[0]} --rMin -1 --rMax 2 -t -1 --expectSignal 1 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2")
            os.system(f"mv *Significance*.root ./Significance/expected/")
            
            os.system(f"combine -M Significance datacard/{year}/{mass}/{file} -n .{mass}.{file.split('.')[0]} --rMin -1 --rMax {rMax} --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2")
            os.system(f"mv *Significance*.root ./Significance/observed/")
