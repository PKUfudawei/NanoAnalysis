import awkward as ak
import numpy as np
import os, uproot, yaml, correctionlib, gc
from glob import glob

source = 'glopart_v2'
base = f'/data/bond/fudawei/public/bbgamma_output/{source}'

filepath = {
    year: {
        r'Top+$\gamma$': sum([glob(f'{base}/mc/{year}/{channel}/*') for channel in ['TTGJets']], []),
        r'Z+$\gamma$':  sum([glob(f'{base}/mc/{year}/{channel}/*') for channel in ['ZGToLLG', 'ZGToJJG']], []),
        r'W+$\gamma$': sum([glob(f'{base}/mc/{year}/{channel}/*') for channel in ['WGToJJG', 'WGToLNuG']], []),
        r'Fake $\gamma$': sum([glob(f'{base}/mc/{year}/{channel}/*') for channel in ['ST', 'TTJets', 'ZJetsToNuNu', 'ZJetsToLL', 'ZJetsToQQ', 'WJetsToLNu', 'WJetsToQQ', 'QCD']], []),
        r'QCD+$\gamma$': sum([glob(f'{base}/mc/{year}/{channel}/*') for channel in ['GJets']], []),
        'data': sum([glob(f'{base}/data/{year}/{channel}/*') for channel in ["EGamma" if year=="2018" else "SinglePhoton"]], []),
    } for year in ['2016pre', '2016post', '2017', '2018']
}

BKG = set([k for k in filepath['2018'].keys()]) - set(['data'])

for year in filepath:
    for job_dir in sum([glob(f'{base}/mc/{year}/{channel}/*') for channel in ['ZpToHG', 'GluGluToZG']], []):
        mass = int(job_dir.split('M')[1].split('_')[0])
        m = mass//1000 if mass%1000 == 0 else mass/1000
        if 'ZpToHG' in job_dir:
            channel = r'$Z^\prime$'+f'({m}TeV,N)'+r'$\to H\gamma$'
        elif '_0p014' in job_dir: 
            channel = f'S({m}TeV,N)'+r'$\to Z\gamma$'
        elif '_5p6' in job_dir: 
            channel = f'S({m}TeV,W)'+r'$\to Z\gamma$'
        elif '_10p0' in job_dir: 
            channel = f'S({m}TeV,VW)'+r'$\to Z\gamma$'
        filepath[year][channel] = [job_dir]

events = {year: {k: [] for k in filepath[year]} for year in filepath}


for year in filepath:
    for channel, folders in filepath[year].items():
        for folder in folders:
            mode = '_'.join(folder.split('/')[-4:-1])
            if len(os.listdir(folder)) != 2:
                continue
            if os.path.exists(os.path.join(folder, f'{mode}.parq')):
                os.rename(os.path.join(folder, f'{mode}.parq'), os.path.join(folder, f'{mode}.parquet'))
            folder_events = ak.from_parquet(os.path.join(folder, f'{mode}.parquet'))
            folder_events = folder_events[folder_events.photon_pt > 225]
            with open(os.path.join(folder, f'{mode}.yaml'), 'r', encoding='utf-8') as f:
                folder_stats = yaml.safe_load(f)
            if 'mc' in mode and ak.mean(folder_events.event_final_weight) > np.sqrt(folder_stats['final']):
                continue

            events[year][channel].append(folder_events)
        if len(events[year][channel]) == 0:
            del events[year][channel]
        else:
            events[year][channel] = ak.concatenate(events[year][channel], axis=0)


def decomposite(signal: ak.Array):
    HWW_decay_mode = ak.fill_none(signal.gen_HWW_decay_mode, 0)

    signal['HWW_4q'] = signal['gen_ZpToH(WW)Gamma'] & (HWW_decay_mode >= 32)
    
    signal['HWW_qqlv'] = signal['gen_ZpToH(WW)Gamma'] & (HWW_decay_mode > 16) & (HWW_decay_mode < 32)
    
    signal['HWW_lvlv'] = signal['gen_ZpToH(WW)Gamma'] & (HWW_decay_mode <= 16)

    return signal

DECOMPOSITE = True
if DECOMPOSITE:
    for y in events:
        for c in list(events[y].keys()):
            if c in BKG or c=='data':
                continue
            if r'H\gamma' in c:
                events[y][c.replace('H', 'H(bb)')] = events[y][c][events[y][c]['gen_ZpToH(bb)Gamma']]
                #events[y][c.replace('H', 'H(cc)')] = events[y][c][events[y][c]['gen_ZpToH(cc)Gamma']]
                del events[y][c]
            if r'Z\gamma' in c:
                events[y][c.replace('Z', 'Z(bb)')] = events[y][c][events[y][c]['gen_GluGluToZ(bb)Gamma']]
                #events[y][c.replace('Z', 'Z(cc)')] = events[y][c][events[y][c]['gen_GluGluToZ(cc)Gamma']]
                del events[y][c]

nlo = correctionlib.CorrectionSet.from_file('src/parameters/nlo_kfactor.json')
list(nlo.keys())

weight = {y: {c: events[y][c].event_final_weight for c in events[y]} for y in events}

for y in events:
    for c in events[y]:
        if c == r'QCD+$\gamma$':
            weight[y][c] = weight[y][c] * nlo['gjets_nlo_kfactor_qcd'].evaluate(events[y][c]['AK8jet_pt'])
        elif c == r'W+$\gamma$':
            weight[y][c] = weight[y][c] * nlo['wjets_nlo_kfactor_qcd'].evaluate(events[y][c]['AK8jet_pt'])
        elif c == r'Z+$\gamma$':
            weight[y][c] = weight[y][c] * nlo['zjets_nlo_kfactor_qcd'].evaluate(events[y][c]['AK8jet_pt'])

gc.collect()
for year in events:
    weight[year]['bkg'] = ak.concatenate([weight[year][c] for c in BKG], axis=0)
    events[year]['bkg'] = ak.concatenate([events[year][c] for c in BKG], axis=0)

gc.collect()
weight['2016'] = {k: ak.concatenate([weight[y][k] for y in ['2016pre', '2016post']], axis=0) for k in weight['2016pre']}
events['2016'] = {k: ak.concatenate([events[y][k] for y in ['2016pre', '2016post']], axis=0) for k in events['2016pre']}

gc.collect()

weight['Run2'] = {k: ak.concatenate([weight[y][k] for y in ['2016', '2017', '2018']], axis=0) for k in weight['2018']}


gc.collect()

events['Run2'] = {k: ak.concatenate([events[y][k] for y in ['2016', '2017', '2018']], axis=0) for k in events['2018']}

mass_cut = {
    'Z': [80, 110],
    'H': [110, 150],
}

for y in ['Run2', '2016', '2017', '2018']:
    for k in events[y]:
        if '(bb)' not in k:
            continue
        mass, fatjet = round(float(k.split('TeV')[0].split('(')[1])*1000), k.split('(')[1][-1]

        width = k.split(',')[1].split(')')[0]
        if width == 'W':
            postfix = '_5p6'
        elif width == 'VW':
            postfix = '_10p0'
        else:
            postfix = ''

        cut = f"""(
        (events[y][k]['photon-jet_deltaR'] > 1.1) &
        (np.abs(events[y][k].photon_eta) < 1.4442) &
        (np.abs(events[y][k].AK8jet_eta) < 2.4) &
        (events[y][k].photon_pt/events[y][k]['photon+jet_mass'] > 0.35) &
        (events[y][k]['AK8jet_particleNet_mass'] > {mass_cut[fatjet][0]}) & 
        (events[y][k]['AK8jet_particleNet_mass'] < {mass_cut[fatjet][1]}) &
        (events[y][k]['AK8jet_Xbb_tagger'] > 0.8)
        )"""

        _events = events[y][k][eval(cut)]
        seen = set()
        flag = []

        
        for i in range(len(_events)):
            event_id = (_events.event_run[i], _events.event_luminosityBlock[i], _events.event_event[i])
            if event_id not in seen:
                flag.append(True)
                seen.add(event_id)
            else:
                flag.append(False)

        if len(flag) - ak.sum(flag) > 0:
            print(y, k, len(flag) - ak.sum(flag))

        events_ = _events[flag]
        out_tree = {
            'fit_mass': np.array(events_['photon+jet_mass']),
            'jet_mass': np.array(events_['AK8jet_particleNet_mass']),
            'weight': np.array(events_['event_final_weight']),
            'tagger': np.array(events_['AK8jet_Xbb_tagger']),
            'event_run': np.array(events_['event_run']),
            'event_luminosityBlock': np.array(events_['event_luminosityBlock']),
            'event_event': np.array(events_['event_event']),
        }
        for k in events_.fields:
            if k.startswith('photon+jet_mass_'):
                out_tree[k] = np.array(events_[k])

        os.makedirs(f'./fit/input/{y}/{mass}{postfix}', exist_ok=True)
        with uproot.recreate(f'./fit/input/{y}/{mass}{postfix}/bbgamma_SR{fatjet}.root') as file:
            file['Events'] = out_tree
