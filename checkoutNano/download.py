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

def download_from_dataset_card(card_path: str, out_basedir: str):
    channel, year, type = card_path.split('/')[-2], card_path.split('/')[-3], card_path.split('/')[-4]
    out_basedir = os.path.join(out_basedir, type, year, channel)
    with open(card_path, 'r') as f:
        dataset = yaml.load(f, Loader=yaml.FullLoader)

    for (k, v) in dataset.items():
        out_dir = out_basedir + v['dataset'] if v['dataset'].startswith('/') else os.path.join(out_basedir, v['dataset'])
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        query_str = f"\"file dataset={v['dataset']} system=rucio\""
        output = subprocess.check_output(f"/cvmfs/cms.cern.ch/common/dasgoclient -query={query_str} -json", shell=True, encoding='utf-8')
        output = json.loads(output)
        for file_info in output:
            file = 'root://cms-xrd-global.cern.ch/'+file_info['file'][0]['name']
            print(f'xrdcp {file} {out_dir}')
            os.system(f'xrdcp {file} {out_dir}')
        print()
        print(f'\t\tDataset of {k} is done')
        print()

    return True
    
def main():
    args = parse_commanline()
    dataset_cards = '/'.join([args.directory, args.type, args.year, args.channel, args.version])+'.yaml'
    dataset_cards = set(glob.glob(dataset_cards))
    total = len(dataset_cards)
    
    with open('./config/condor_cfg.yaml', 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
    finished = 0
    for dataset_card in dataset_cards:
        finished += download_from_dataset_card(card_path=dataset_card, out_basedir=config['out_dir'])
        print('====================================================================')
        print('==> Downloading process:\t%.1f'%(100*finished/total))
        print('====================================================================')
        
if __name__ == "__main__":
    main()