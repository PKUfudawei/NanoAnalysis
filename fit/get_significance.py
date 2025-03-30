import os, uproot
import pandas as pd


def run_significance(input_dir='./', year='Run2', out_dir='./Significance', width=None):
    os.makedirs(os.path.join(out_dir, year), exist_ok=True)
    datacard_dir = os.path.join(input_dir, 'datacard', year)
    for mass in os.listdir(datacard_dir):
        for f in set(os.listdir(os.path.join(datacard_dir, mass))):
            if (width is not None and not f.endswith(f'_{width}.txt')) or not f.endswith(f'.txt'):
                continue
            file = os.path.join('datacard', year, mass, f)
            print(file)
            rMax = 3 if 'SRH' in file else 10
            os.system(f"rm -rf {os.path.join(input_dir, '*Significance*.root')}")
            os.system(f"cd {input_dir}; combine -M Significance {file} -n .{mass}.{f.split('.')[0]}.expected --rMin -1 --rMax 3 -t -1 --expectSignal 1 --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2")
            os.system(f"cd {input_dir}; combine -M Significance {file} -n .{mass}.{f.split('.')[0]}.observed --rMin -1 --rMax {rMax} --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2")
            os.system(f"mv {os.path.join(input_dir, '*Significance*.root')} {os.path.join(out_dir, year)}")
            print()


def convert_to_csv(in_dir='./Significance', year='Run2', signal_region='SRH', width='N', out_dir='./significance_results'):
    significance = {m: {'expected': -1, 'observed': -1} for m in range(700, 3501, 50)} 

    root_dir = os.path.join(in_dir, year)
    for f in set(os.listdir(root_dir)):
        if not (f.endswith('.root') and 'Significance' in f and f'{signal_region}_{width}' in f):
            continue
        mass = int(f.split('.')[1])
        mode = f.split('.')[3]
        if mode not in significance[mass].keys():
            raise ValueError()
        print(f)
        stats = uproot.open(os.path.join(root_dir, f))
        significance[mass][mode] = float(stats['limit']['limit'].array()[0])

    table = {
        'mass': [m for m in range(700, 3501, 50)],
        'significance_exp': [significance[m]['expected'] for m in range(700, 3501, 50)],
        'significance_obs': [significance[m]['observed'] for m in range(700, 3501, 50)],
    }

    df = pd.DataFrame(table).round(4)
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(os.path.join(out_dir, f'significance_{signal_region}_{width}.csv'), index=False, float_format="%.4f")


if __name__ == "__main__":
    run_significance()
    convert_to_csv(signal_region='SRH', width='N')
    convert_to_csv(signal_region='SRH1', width='N')
    convert_to_csv(signal_region='SRH2', width='N')
    convert_to_csv(signal_region='SRZ', width='N')
    convert_to_csv(signal_region='SRZ1', width='N')
    convert_to_csv(signal_region='SRZ2', width='N')
    convert_to_csv(signal_region='SRZ', width='W')
    convert_to_csv(signal_region='SRZ1', width='W')
    convert_to_csv(signal_region='SRZ2', width='W')
    convert_to_csv(signal_region='SRZ', width='VW')
    convert_to_csv(signal_region='SRZ1', width='VW')
    convert_to_csv(signal_region='SRZ2', width='VW')
