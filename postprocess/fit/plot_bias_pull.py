import ROOT
import numpy as np
import argparse, os

ROOT.gROOT.SetBatch(True)


def parse_commandline():
    parser = argparse.ArgumentParser(description='parametric fitting')
    parser.add_argument('-y', '--year', help='To specify which year', choices=('2016pre', '2016post', '2017', '2018', 'Run2'), default='Run2')
    parser.add_argument('-t', '--truth_function', help='To specify which truth function', choices=('expow2', 'invpow2'))
    parser.add_argument('-m', '--signal_mass', help='To specify the mass of signal resonance', type=int)
    parser.add_argument('-R', '--SR', help='To specify which signal region', choices=('SR1', 'SR2'))
    args = parser.parse_args()
    return args


def main(truth_function, signal_mass, SR, fit_function='dijet2'):
    N_toys = 1000
    r_truth = 1

    year = 'Run2'

    name = "truth_%s_fit_%s"%(truth_function, fit_function)

    # Open file with fits
    f = ROOT.TFile("higgsCombine.bias_%s.MultiDimFit.mH120.123456.root"%name)
    t = f.Get("limit")

    hist_pull = ROOT.TH1F("pull_%s"%name, "Pull distribution: truth=%s, fit=%s"%(truth_function, fit_function), 80, -4, 4)
    hist_pull.GetXaxis().SetTitle("Pull = (r_{truth}-r_{fit})/#sigma_{fit}")
    hist_pull.GetYaxis().SetTitle("Entries")

    sigma_values = np.array([])

    for i_toy in range( N_toys ):
        # Best-fit value
        t.GetEntry(i_toy*3)
        r_fit = getattr(t, "r")

        # -1 sigma value
        t.GetEntry(i_toy*3+1)
        r_lo = getattr(t, "r")

        # +1 sigma value
        t.GetEntry(i_toy*3+2)
        r_hi = getattr(t, "r")

        diff = r_truth-r_fit
        # Use uncertainty depending on where mu_truth is relative to mu_fit
        if diff > 0: sigma = abs(r_hi-r_fit)
        else: sigma = abs(r_lo-r_fit)

        if sigma != 0:
            sigma_values = np.append( sigma_values , sigma )
        else: 
            sigma = sigma_values.mean()

        if sigma != 0:
            hist_pull.Fill( diff/sigma )

    canv = ROOT.TCanvas()
    hist_pull.Draw()

    # Fit Gaussian to pull distribution
    ROOT.gStyle.SetOptFit(111)
    hist_pull.Fit("gaus")

    canv.SaveAs(f"../plots/fit/{year}/bias_pull_{name}_{signal_mass}_{SR}.pdf")


if __name__ == "__main__":
    args = parse_commandline()
    m = args.signal_mass
    truth_function = args.truth_function
    SR = args.SR
    os.system(f"combine -M GenerateOnly datacard/Run2/{truth_function}/datacard_{m}_{SR}.root -t 1000 -n .generate_{truth_function} --expectSignal 1 --saveToys")
    os.system(f"combine -M MultiDimFit datacard/Run2/dijet2/datacard_{m}_{SR}.root -t 1000 -n .bias_truth_{truth_function}_fit_dijet2 --expectSignal 1 --toysFile higgsCombine.generate_{truth_function}.GenerateOnly.mH120.123456.root --algo singles;")
    main(truth_function=truth_function, signal_mass=m, SR=SR)
