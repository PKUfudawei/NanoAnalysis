# checkoutNano
The repo to check or download **NanoAOD** files from the datasets.

To use this repo, download firstly,
```bash
git clone https://github.com/PKU-Hep-Group/checkoutNano.git
cd checkoutNano
```

# Configuration
Several configurations files are in `config` folder:

1. In `condor_cfg.yaml`:
   - `job_dir`: define the place to put the job submission configurations, logs, and excutable files
   - `out_dir`: the place to put the datasets downloaded from DAS, EOS, or your own directories
2. To view `datasets` structure hierarchy: `tree datasets`

# How to run
The python script to use is `fetch_file.py`.

## The arguements
Firstly, check the arguements in `fetch_file.py`,

```shell
$ python3 fetch_file.py --help
usage: fetch_file.py [-h] [-t {data,mc}] [-v VERSION] [-c CHANNEL] [-y YEAR]
                     [-f FLAVOUR] [-cf] [-ds DATASET] [--dryrun]
                     [--sleep SLEEP] [-a | -s SAMPLE [SAMPLE ...]]
                     {look,submit,status}

manual to this script

positional arguments:
  {look,download,status}  which mode: look, download, status

optional arguments:
  -h, --help            show this help message and exit
  -t {data,mc}, --type {data,mc}
                        which type: data, mc
  -v VERSION, --version VERSION
                        which nanoAOD version: nanov7, nanov8, nanov9, ...,
                        default: nanov9
  -c CHANNEL, --channel CHANNEL
                        which channel: vhjj, ssww, your own channel, ...
  -y YEAR, --year YEAR  which year to run: 2016, 2017, 2018
  -f FLAVOUR, --flavour FLAVOUR
                        job flavour: espresso, longlunch, workday, ...
  -cf, --checkfile      check if the local file is zombie, can be slow
  -ds DATASET, --dataset DATASET
                        The das dataset, only for mode: look
  --dryrun              only generate submit files, but not run
  --sleep SLEEP         random sleep time from zero to setting time (min),
                        default 20
  -a, --all             run on all, default: True
  -s SAMPLE [SAMPLE ...], --sample SAMPLE [SAMPLE ...]
                        run on which sample, can be more than 1
```

## Some examples

Set the proxy firstly,

```bash
voms-proxy-init -voms cms -rfc -valid 192:00 --out ~/.proxy # "--out ~/.proxy" is just needed for once
export X509_USER_PROXY=$HOME/.proxy
```

Just to have a look at the files in the dataset using the query command from [dasgoclient](https://cmsweb.cern.ch/das/cli),
```bash
python3 fetch_file.py look -ds /DoubleMuon/Run2018B-UL2018_MiniAODv2_NanoAODv9-v1/NANOAOD
```

Generate configuration files with option **--dryrun**, will not submit condor jobs to download files, it's always good to check before submitting jobs.
```bash
python3 fetch_file.py download -y 2018 -c vhjj -t data -v nanov9 --dryrun
```

To submit condor jobs, you can run,
```bash
# for data
python3 fetch_file.py download -y 2018 -c vhjj -t data -v nanov9
# for mc
python3 fetch_file.py download -y 2018 -c vhjj -t mc -v nanov9
# only for some samples
python3 fetch_file.py download -y 2018 -c vhjj -t data -v nanov9 -s EGamma_Run2018D EGamma_Run2018A
```

To check the status, the script will only check locally exist files,
```bash
# for data
python3 fetch_file.py status -y 2018 -c vhjj -t data -v nanov9
# for mc
python3 fetch_file.py status -y 2018 -c vhjj -t mc -v nanov9
# only for some samples
python3 fetch_file.py status -y 2018 -c vhjj -t data -v nanov9 -s EGamma_Run2018D EGamma_Run2018A
```

### Resubmit jobs

In case the jobs fail to downlaod sample files, you can simply rerun the `dowanload` mode with `fetch_file.py`, exist files will not be downloaded again. Or just go the job configration folder, use `bash` commands to resubmit the jobs, e.g.,
```bash
for i in *jid; do sed -i "s/testmatch/testmatch/g" ${i/jid/sub}; condor_submit ${i/jid/sub}; done
```

## For future
We need to download files to a common place, then the script can check if the files are already there, if yes, it will not download it again.
