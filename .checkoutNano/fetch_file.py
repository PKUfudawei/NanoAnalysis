import argparse
import utils.functions as func
import yaml


parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('mode', help='which mode: look, download, status', choices=('look', 'download', 'status'), default='status')
parser.add_argument('-t','--type', help='which type: data, mc', choices=('data', 'mc'), default='mc')
parser.add_argument('-v','--version', help='which nanoAOD version: nanov7, nanov8, nanov9, ..., default: nanov9', default='nanov9')
parser.add_argument('-c', '--channel', help='which channel: vhjj, ssww, your own channel, ...', default='vhjj')
parser.add_argument('-y', '--year', help='which year to run: 2016, 2017, 2018', default='2018')
parser.add_argument('-f', '--flavour', help='job flavour: espresso, longlunch, workday, ...', default='testmatch')
parser.add_argument('-cf', '--checkfile', help='check if the local file is zombie, can be slow', action='store_true', default=False)
parser.add_argument('-ds', '--dataset', help='The das dataset, only for mode: look', default='/DoubleMuon/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD')
parser.add_argument('--dryrun', help='only generate submit files, but not run', action='store_true', default=False)
# avoid downloading at the same time
parser.add_argument('--sleep', help='random sleep time from zero to setting time (min), default 20', default="20")
group = parser.add_mutually_exclusive_group()  # type: argparse._MutuallyExclusiveGroup
group.add_argument('-a', '--all', help='run on all, default: True', action='store_false', default=True)
group.add_argument('-s', '--sample', help='run on which sample, can be more than 1', nargs='+')
# ,default=['WpWpJJ_EWK','WpWpJJ_EWK_powheg','WmWmJJ_EWK_powheg'])
args = parser.parse_args()


def get_cfg():
    cfg = {
        'year': args.year,
        'channel': args.channel,
        'nano_ver': args.version,
        'sample_type': args.type,
        'dryrun': args.dryrun,
        'checkfile': args.checkfile,
        'sleep': args.sleep,
        'all': args.all,
        'sample': args.sample,
        'jflavour': args.flavour,
    }

    with open('./config/condor_cfg.yaml', 'r') as f:
        condor_cfg = yaml.load(f, Loader=yaml.FullLoader)
    job_dir = condor_cfg['job_dir']
    out_dir = condor_cfg['out_dir']  # cms connect output place
    
    with open(
        f"./datasets/{cfg['sample_type']}/{cfg['year']}/{cfg['channel']}/{cfg['nano_ver']}.yaml", mode='r'
    ) as f:
        ds_yml = yaml.load(f, Loader=yaml.FullLoader)
    #with open(f"./datasets/{ds_cfg[cfg['channel']][cfg['nano_ver']][cfg['sample_type']][int(cfg['year'])]}", 'r') as f:
    #    ds_yml = yaml.load(f, Loader=yaml.FullLoader)
    print(ds_yml)
    cfg['job_dir'] = job_dir
    cfg['out_dir'] = out_dir
    cfg['ds_yml'] = ds_yml
    
    return cfg


def to_download(year, sample_type):
    print("===> to_submit:", "year:", year, "sample_type:", sample_type)
    cfg = get_cfg()

    if cfg['sample']:
        for isp in cfg['sample']:
            if isp in cfg['ds_yml']:
                func.submit_command(cfg, isp)
            else:
                print("--- submit info:", "no sample: ", isp, "in", year)
    else:
        for isp in cfg['ds_yml']:
            func.submit_command(cfg, isp)
    print("<=== to_submit END :)")


def to_status(year, sample_type):
    print("===> to_status: year:", year, "sample_type:", sample_type)
    cfg = get_cfg()
    if cfg['sample']:
        for isp in cfg['sample']:
            if isp in cfg['ds_yml']:
                func.status_command(cfg, isp)
            else:
                print("--- status info:", "no sample: ", isp, "in", year)
    else:
        for isp in cfg['ds_yml']:
            func.status_command(cfg, isp)
    print("<=== to_status END :)")


def have_a_look(ds):
    use_gfal, file_dict = func.get_file_info(ds)
    for ifile, ifname in enumerate(file_dict):
        print("===> File index:", ifile, ", file name:", ifname)
        print("<=== Available path:",file_dict[ifname])


if __name__ == '__main__':
    if args.mode == "look":
        have_a_look(args.dataset)
    elif args.mode == 'download':
        to_download(args.year, args.type)
    elif args.mode == 'status':
        to_status(args.year, args.type)
    else:
        print("===>", "no such mode")
