import os, uproot, yaml
import pandas as pd


def run_rate_limit(input_dir='./', year='Run2', out_dir='./AsymptoticLimits', width='narrow'):
    os.makedirs(os.path.join(out_dir, year), exist_ok=True)
    datacard_dir = os.path.join(input_dir, 'datacard', year)
    for mass in os.listdir(datacard_dir):
        if width == 'narrow' and '_' in mass:
            continue
        elif width == 'wide' and '_5p6' not in mass:
            continue
        elif width == 'very_wide' and '_10p0' not in mass:
            continue
        elif width not in ['narrow', 'wide', 'very_wide']:
            raise ValueError()
        for f in set(os.listdir(os.path.join(datacard_dir, mass))):
            if not f.endswith('.txt'):
                continue
            file = os.path.join('datacard', year, mass, f)
            print(file)
            os.system(f"rm -rf {os.path.join(input_dir, '*AsymptoticLimits*.root')}")
            os.system(f"cd {input_dir}; combine -M AsymptoticLimits {file} -n .{mass}.{f.split('.')[0]} --X-rtd REMOVE_CONSTANT_ZERO_POINT=1 --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2")
            os.system(f"mv {os.path.join(input_dir, '*AsymptoticLimits*.root')} {os.path.join(out_dir, year)}")
            print()
           

def convert_to_csv(in_dir='./AsymptoticLimits', year='Run2', signal_region='SRH', width=None, out_dir='./limit_results', xsec_file='../src/parameters/cross-section.yaml'):
    rate_limit = {} 

    root_dir = os.path.join(in_dir, year)
    for f in os.listdir(root_dir):
        if not (f.endswith('.root') and 'AsymptoticLimits' in f):
            continue
        mass, SR = f.split('.')[1:3]
        if SR != signal_region:
            continue
        if width == 'narrow' and '_' in mass:
            continue
        elif width == 'wide' and '_5p6' not in mass:
            continue
        elif width == 'very_wide' and '_10p0' not in mass:
            continue
        elif width not in ['narrow', 'wide', 'very_wide']:
            raise ValueError()

        mass = int(mass.split('_')[0])
        stats = uproot.open(os.path.join(root_dir, f))
        limits = stats['limit']['limit'].array()
        rate_limit[mass] = {
            'limit_m2': float(limits[0]),
            'limit_m1': float(limits[1]),
            'exp': float(limits[2]),
            'limit_p1': float(limits[3]),
            'limit_p2': float(limits[4]),
            'obs': float(limits[5]),
        }

    if signal_region == 'SRH':
        with open(xsec_file, 'r', encoding='utf-8') as f:
            xsec_info = yaml.safe_load(f)['ZpToHG']
        # convert pb to fb
        xsec = {int(k.split('_M')[1].split('_')[0]): float(v)*1e3 for k, v in xsec_info.items()}
    else:
        xsec = {m: 10 for m in range(700, 3501, 50)}

    limit = {
        'mass': [m for m in range(700, 3501, 50)],
        'limit_p1': [rate_limit[m]['limit_p1'] * xsec[m] for m in range(700, 3501, 50)],
        'limit_p2': [rate_limit[m]['limit_p2'] * xsec[m] for m in range(700, 3501, 50)],
        'limit_m1': [rate_limit[m]['limit_m1'] * xsec[m] for m in range(700, 3501, 50)],
        'limit_m2': [rate_limit[m]['limit_m2'] * xsec[m] for m in range(700, 3501, 50)],
        'exp': [rate_limit[m]['exp'] * xsec[m] for m in range(700, 3501, 50)],
        'obs': [rate_limit[m]['obs'] * xsec[m] for m in range(700, 3501, 50)],
    }

    df = pd.DataFrame(limit).round(4)
    name = f'limit_{signal_region}.csv' if signal_region == 'SRH' else f'limit_{signal_region}_{width}.csv'
    df.to_csv(os.path.join(out_dir, name), index=False, float_format="%.4f")


if __name__ == "__main__":
    run_rate_limit(width='wide')
    #convert_to_csv(signal_region='SRH')
    #convert_to_csv(signal_region='SRZ', width='narrow')
    convert_to_csv(signal_region='SRZ', width='wide')
    #convert_to_csv(signal_region='SRZ', width='very_wide')

