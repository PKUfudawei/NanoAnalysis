# NanoAnalysis
[HEP-EX](https://inspirehep.net/) analyzing framework for NanoAOD-like data in [CMS](https://cms.cern/) experiment

## Brief introduction
- ### datasets/
Hierarchically stores dataset info in the structure of  `datasets/$type/$year/$channel/$division/$version.yaml`  
- ### src/
Containing all source code to process NanoAOD into slimmed files (e.g., `*.parq`) and to calculate correction and uncertainty by [coffea](https://coffeateam.github.io/coffea/) (read `*.root` files into vector-like structure, e.g., [NanoEvents](https://coffeateam.github.io/coffea/modules/coffea.nanoevents.html#module-coffea.nanoevents)) and [awkward](https://awkward-array.org/quickstart.html) (a numpy-like library but more capable to process nested and variable-sized data in json-like format)  
- ### condor/
Use `parseDatasets.py` to parse filelists from `*.yaml` and generate corresponding job description files as `*.submit`. Then use integrated `condor_submit.py` script to submit condor jobs for processing NanoAOD files with codes in `src/`  
- ### notebooks/
Jupyter-notebboks to analyze/postprocess `*.parquet` coming from processed NanoAOD data. Produce the final plots here!  
- ### .checkoutNano/
It is deprecated but codes in it may be reused  
