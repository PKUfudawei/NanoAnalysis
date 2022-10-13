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

def download_from_dataset_card(card_path: str, out_dir: str):
    with open(card_path, 'r') as f:
        dataset = yaml.load(f, Loader=yaml.FullLoader)
        
    for (k, v) in dataset.items():
        out_dir = os.path.join(out_dir, k)+v['dataset']
        query_str = f"\"file dataset={v['dataset']} system=rucio\""
        output = subprocess.check_output(f"/cvmfs/cms.cern.ch/common/dasgoclient -query={query_str} -json", shell=True, encoding='utf-8')
        file = 'root://cms-xrd-global.cern.ch/'+json.loads(output)[0]['file'][0]['name']
        os.system(f'xrdcp {file} {out_dir}')

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
        finished += download_from_dataset_card(card_path=dataset_card, out_dir=config['out_dir'])
        print('====================================================================')
        print('==> Downloading process:\t%.1f'%(100*finished/total))
        print('====================================================================')
        
if __name__ == "__main__":
    main()