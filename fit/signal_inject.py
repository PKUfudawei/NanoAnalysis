import os

for m in [700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 3000, 3500]:
    os.system(f"combine -M MultiDimFit datacard/Run2/datacard_{m}_combine.txt -t 200 --expectSignal=1 --setParameters pdfindex_SR1=0,pdfindex_SR2=0 -n .{m}")
