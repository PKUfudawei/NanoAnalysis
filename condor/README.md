## To submit condor jobs
```bash
myschedd bump
type=
year=
channel=
for i in submit/$type/$year/$channel/*.submit; do condor_submit $i; done
```

## *References*
[CERN Batch Docs](https://batchdocs.web.cern.ch/index.html)
[HTCondor Manual](https://htcondor.readthedocs.io/en/latest/index.html)