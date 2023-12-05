import os

for index, generate_model in enumerate(['dijet2', 'expow1', 'invpow2']):
    os.system(f"combine -M GenerateOnly datacard/Run2/datacard_1000_SR1.txt -t 1000 --fixedSignalStrength=0 --setParameters pdfindex_Tag0={2*index} --freezeParameters pdfindex_Tag0 -n .{generate_model} --saveToys")
    os.system(f"combine -M GoodnessOfFit datacard/Run2/datacard_1000_SR1.txt -t 1000 --fixedSignalStrength=0 --setParameters pdfindex_Tag0={2*index} --freezeParameters pdfindex_Tag0 -n .{generate_model} --toysFile higgsCombine.{generate_model}.GenerateOnly.mH120.123456.root --algo=saturated")
    higher_model = generate_model[:-1] + str(int(generate_model[-1])+1)
    os.system(f"combine -M GoodnessOfFit datacard/Run2/datacard_1000_SR1.txt -t 1000 --fixedSignalStrength=0 --setParameters pdfindex_Tag0={2*index+1} --freezeParameters pdfindex_Tag0 -n .{higher_model} --toysFile higgsCombine.{generate_model}.GenerateOnly.mH120.123456.root --algo=saturated")
    
    os.system(f"combine -M GoodnessOfFit datacard/Run2/datacard_1000_SR1.txt --fixedSignalStrength=0 --setParameters pdfindex_Tag0={2*index} --freezeParameters pdfindex_Tag0 --algo=saturated -t -1 -n .data_{generate_model}")
    os.system(f"combine -M GoodnessOfFit datacard/Run2/datacard_1000_SR1.txt --fixedSignalStrength=0 --setParameters pdfindex_Tag0={2*index+1} --freezeParameters pdfindex_Tag0 --algo=saturated -t -1 -n .data_{higher_model}")
    