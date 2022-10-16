#!/usr/bin/env python3
import argparse
import os
import glob
import yaml
import json
import subprocess

def parse_commanline():
    parser = argparse.ArgumentParser(description='Script to check if each condor job is done')
    parser.add_argument('-d', '--directory', help='To specify base directory', default=os.path.abspath('./datasets'))
    parser.add_argument('-t', '--type', help='To specify jobs in mc/ or data/', choices=('data', 'mc'), default='*')
    parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
    parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
    parser.add_argument('-v', '--version', help='To specify jobs in which nanoAOD version', default='*')
    args = parser.parse_args()
    return args

def dataset_to_filelist(card_path: str):
    name = os.path.join(*card_path.split('/')[-4:-1])
    with open(card_path, 'r') as f:
        dataset = yaml.load(f, Loader=yaml.FullLoader)

    for (k, v) in dataset.items():
        filelist = []
        query_str = f"\"file dataset={v} system=rucio\""
        output = subprocess.check_output(f"/cvmfs/cms.cern.ch/common/dasgoclient -query={query_str} -json", shell=True, encoding='utf-8')
        output = json.loads(output)
        for file_info in output:
            filelist.append('root://cms-xrd-global.cern.ch/'+file_info['file'][0]['name'])
        
        if not os.path.exists(f'./filelists/{name}'):
            os.makedirs(f'./filelists/{name}')
        with open(f'./filelists/{name}/{k}.txt', 'w') as f:
            f.write('\n'.join(filelist))
        if not os.path.exists(f'./output/{name}/{k}'):
            os.makedirs(f'./output/{name}/{k}')
        if not os.path.exists(f'./log/{name}/{k}'):
            os.makedirs(f'./log/{name}/{k}')
            
    print(f'==> Generated filelists/{name}/{k}.txt from {card_path}')
    return True

def filelist_to_submit(filelist: str, template: str):
    name = os.path.join(*filelist.split('.')[-2].split('/')[-4:])
    path = os.path.join('./submit', *name.split('/')[:-1])
    if not os.path.exists(path):
        os.makedirs(path)
    with open(f'./submit/{name}.submit', 'w') as f:
        f.write(template.replace('$template', name))
        
    print(f'==> Generated submit/{name}.submit from {filelist}')
    return True


def main() -> None:
    args = parse_commanline()
    dataset_cards = os.path.join(args.directory, args.type, args.year, args.channel, args.version+'.yaml')
    dataset_cards = set(glob.glob(dataset_cards))
    print(f'==> Start generating filelist(s) from {len(dataset_cards)} dataset cards')
    
    succeeded = 0
    for dataset_card in dataset_cards:            
        succeeded += dataset_to_filelist(card_path=dataset_card)
    print(f'==> Successfully generated {succeeded} filelist(s) in total')
    
    print('################################################'*3)    
    filelists = os.path.join('filelists', args.type, args.year, args.channel, args.version+'.txt')
    filelists = set(glob.glob(filelists))
    print(f'==> Start generating condor-submit job(s) from   {len(filelists)} filelists')
    
    with open('./template.submit', 'r') as f:
        template = f.read()
        
    succeeded = 0
    for filelist in filelists:         
        succeeded += filelist_to_submit(filelist=filelist, template=template)
    print(f'==> Successfully generated {succeeded} condor-submit file(s) in total')
    
    
if __name__ == "__main__":
    main()