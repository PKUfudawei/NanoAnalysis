import ROOT
import numpy as np
import argparse, os, yaml

ROOT.gROOT.SetBatch(True)

N_toys = 1000
year = 'Run2'


def plot(truth_function, signal_mass, SR, fit_function='dijet2', r_truth=0):
    name = f"{truth_function}_{fit_function}"

    # Open file with fits
    f = ROOT.TFile("fitDiagnosticsTest.root")
    t = f.Get("tree_fit_sb")

    hist_pull = ROOT.TH1F("pull_%s"%name, "Pull distribution: generate=%s, fit=%s"%(truth_function, fit_function), 50, -5, 5)
    hist_pull.GetXaxis().SetTitle("Pull = (r_{truth}-r_{fit})/#sigma_{fit}")
    hist_pull.GetYaxis().SetTitle("Entries")

    sigma_values = np.array([])

    for i_toy in range(N_toys):
        # Best-fit value
        t.GetEntry(i_toy*3)
        r_fit = getattr(t, "r")
        t.GetEntry(i_toy*3+1)
        r_lo = getattr(t, "r")
        t.GetEntry(i_toy*3+2)
        r_hi = getattr(t, "r")
        diff = r_fit - r_truth
        # Use uncertainty depending on where mu_truth is relative to mu_fit
        if diff < 0:
            sigma = abs(r_hi - r_fit)
        else:
            sigma = abs(r_fit - r_lo)
        if sigma != 0:
            sigma_values = np.append(sigma_values, sigma)
        else: 
            sigma = sigma_values.mean()
        if sigma != 0:
            hist_pull.Fill(diff/sigma)

    canv = ROOT.TCanvas()
    hist_pull.Draw()

    # Fit Gaussian to pull distribution
    ROOT.gStyle.SetOptFit(111)
    hist_pull.Fit("gaus")
    fit = hist_pull.GetFunction("gaus")

    if not os.path.exists(f"../plots/fit/{year}/"):
        os.makedirs(f"../plots/fit/{year}/")
    canv.SaveAs(f"../plots/fit/{year}/bias_pull_{name}_{signal_mass}_{SR}.pdf")

    return fit.GetParameter(0), fit.GetParameter(1), fit.GetParameter(2)


def main():
    functions = ['dijet3', 'expow1', 'expow2', 'invpow2', 'invpow3']
    r_truth = 0
    A, mean, sigma = {func:{'SR1':{}, 'SR2':{}} for func in functions}, {func:{'SR1':{}, 'SR2':{}} for func in functions}, {func:{'SR1':{}, 'SR2':{}} for func in functions}
    
    for index, func in enumerate(functions):
        for SR in ['SR1', 'SR2']:
            for m in [700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 3000, 3500]:
                os.system(f"combine -M GenerateOnly datacard/Run2/datacard_{m}_{SR}.txt --setParameters pdfindex_{SR}={index+1} --toysFrequentist -t {N_toys} --expectSignal {r_truth} --saveToys -m 125 --freezeParameters pdfindex_{SR} --rMin -10 --rMax 10")
                os.system(f"combine -M FitDiagnostics datacard/Run2/datacard_{m}_{SR}.txt --setParameters pdfindex_{SR}=0 --toysFile higgsCombineTest.GenerateOnly.mH125.123456.root -t {N_toys} -m 125 --rMin -10 --rMax 10 --freezeParameters pdfindex_{SR} --cminDefaultMinimizerStrategy=0")
                A[func][SR][m], mean[func][SR][m], sigma[func][SR][m] = plot(truth_function=func, signal_mass=m, SR=SR, r_truth=r_truth)

    with open('./bias_pull.yaml', 'w', encoding='utf-8') as f:
        yaml.dump({'A': A, 'mean': mean, 'sigma': sigma}, f)


if __name__ == "__main__":
    main()
