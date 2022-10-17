# #To submit condor jobs
```bash
myschedd bump
type=
year=
channel=
for i in submit/$type/$year/$channel/*.submit; do condor_submit $i; done
```
