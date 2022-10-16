import argparse
import os
import glob
import yaml
import json
import subprocess

def parse_commanline():
    parser = argparse.ArgumentParser(description='Script to check if each condor job is done')
    parser.add_argument('-d', '--directory', help='To specify base directory', default=os.path.abspath('../datasets'))
    parser.add_argument('-t', '--type', help='To specify jobs in mc/ or data/', choices=('data', 'mc'), default='*')
    parser.add_argument('-y', '--year', help='To specify jobs in which year', default='*')
    parser.add_argument('-c', '--channel', help='To specify jobs in which channel', default='*')
    parser.add_argument('-v', '--version', help='To specify jobs in which nanoAOD version', default='*')
    args = parser.parse_args()
    return args

def filelist_from_dataset_card(card_path: str):
    channel, year, type = card_path.split('/')[-2], card_path.split('/')[-3], card_path.split('/')[-4]
    with open(card_path, 'r') as f:
        dataset = yaml.load(f, Loader=yaml.FullLoader)

    filelist = []
    for (k, v) in dataset.items():
        query_str = f"\"file dataset={v['dataset']} system=rucio\""
        output = subprocess.check_output(f"/cvmfs/cms.cern.ch/common/dasgoclient -query={query_str} -json", shell=True, encoding='utf-8')
        output = json.loads(output)
        for file_info in output:
            filelist.append('root://cms-xrd-global.cern.ch/'+file_info['file'][0]['name'])
        with open(f'./{channel}.txt', 'w') as f:
            f.write('\n'.join(filelist))
    print(f'==> Generated {channel}.txt from {card_path}')
    
    return True

def main() -> None:
    args = parse_commanline()
    dataset_cards = '/'.join([args.directory, args.type, args.year, args.channel, args.version])+'.yaml'
    dataset_cards = set(glob.glob(dataset_cards))
    print(f'==> Start generating {len(dataset_cards)} filelist(s) from {dataset_cards}')
    
    succeeded = 0
    for dataset_card in dataset_cards:
        succeeded += filelist_from_dataset_card(card_path=dataset_card)
    print(f'==> Successfully generated {succeeded} file list(s) in total')
        
if __name__ == "__main__":
    main()