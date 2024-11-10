## Parse datasets to generate filelists and submit files
```bash
./parseDatasets.py (-t $type -y $year -c $channel -v $version -m $machine -rm $bool)
```

## To submit condor jobs
```bash
./condor_submit.py (-t $type -y $year -c $channel -b $bool -j $job)
```

## *References*
[CERN Batch Docs](https://batchdocs.web.cern.ch/index.html)  
[HTCondor Manual](https://htcondor.readthedocs.io/en/latest/index.html)
