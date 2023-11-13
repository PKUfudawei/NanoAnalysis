#!/usr/bin/env python3
import ROOT, os, yaml, argparse
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)


def parse_commandline():
    parser = argparse.ArgumentParser(description='parametric fitting')
    parser.add_argument('-y', '--year', help='To specify which year', choices=('2016pre', '2016post', '2017', '2018', 'Run2'))
    parser.add_argument('-m', '--signal_mass', help='To specify the mass of signal resonance', type=int)
    parser.add_argument('-R', '--SR', help='To specify which signal region', choices=('SR1', 'SR2'))
    parser.add_argument('-l', '--fit_range_low', help='To specify the lower bound of fitting range', default=720, type=int)
    parser.add_argument('-u', '--fit_range_up', help='To specify the higher bound of fitting range', default=4000, type=int)
    args = parser.parse_args()
    return args


def fit_signal(year):
    # # Signal modelling
    f = ROOT.TFile(f"input/{year}/mc_signal_Hbb_{signal_mass}.root", "r")
    # Load TTree
    tree = f.Get("Events")

    # Define mass and weight variables
    mass_Zprime = ROOT.RooRealVar("mass_Zprime", "mass_Zprime", signal_mass, fit_range_low, fit_range_up)
    weight = ROOT.RooRealVar("weight", "weight", 0, -10, 10)
    mass_Higgs = ROOT.RooRealVar("mass_Higgs", "mass_Higgs", 125, 0, 999)
    tagger_Hbb = ROOT.RooRealVar("tagger_Hbb", "tagger_Hbb", 0, 0, 2)

    # Convert to RooDataSet

    mc = ROOT.RooDataSet("signal_Hbb", "signal_Hbb", tree, ROOT.RooArgSet(mass_Zprime, weight, mass_Higgs, tagger_Hbb), SR_cut, "weight")

    # Lets plot the signal mass distribution
    can = ROOT.TCanvas()
    plot = mass_Zprime.frame()
    mc.plotOn(plot)
    plot.Draw()
    can.Update()
    if not os.path.exists(f'../plots/fit/{year}'):
        os.makedirs(f'../plots/fit/{year}')
    can.SaveAs(f"../plots/fit/{year}/signal_Hbb_{signal_mass}_mass_Zprime.pdf")

    # Introduce RooRealVars into the workspace for the fitted variable
    x0 = ROOT.RooRealVar("x0", "x0", signal_mass, signal_mass - 200, signal_mass + 200)
    sigmaL = ROOT.RooRealVar("sigmaL", "sigmaL", 100, 0, 600)
    sigmaR = ROOT.RooRealVar("sigmaR", "sigmaR", 100, 0, 600)
    alphaL = ROOT.RooRealVar("alphaL", "alphaL", 2, 0, 7)
    alphaR = ROOT.RooRealVar("alphaR", "alphaR", 2, 0, 7)
    nL = ROOT.RooRealVar("nL", "nL", 0.1, -10, 10)
    nR = ROOT.RooRealVar("nR", "nR", -1, -10, 10)

    # Define the Gaussian with mean=MH and width=sigma
    model_signal = ROOT.RooCrystalBall("model_signal", "model_signal", mass_Zprime, x0, sigmaL, sigmaR, alphaL, nL, sigmaR, nR)

    # Fit Gaussian to MC events and plot
    model_signal.fitTo(mc, ROOT.RooFit.SumW2Error(True))
    model_signal.fitTo(mc, ROOT.RooFit.SumW2Error(True))

    x0.setConstant(True)
    sigmaL.setConstant(True)
    sigmaR.setConstant(True)
    alphaL.setConstant(True)
    alphaR.setConstant(True)
    nL.setConstant(True)
    nR.setConstant(True)

    can = ROOT.TCanvas()
    plot = mass_Zprime.frame()
    mc.plotOn(plot)
    model_signal.plotOn(plot, ROOT.RooFit.LineColor(2))
    plot.Draw()
    can.Update()
    can.Draw()

    mass_Zprime.setBins(160)
    # hist = ROOT.RooDataHist("hist", "hist", mass_Zprime, mc)
    # print("==> chi^2/ndf = ", ROOT.RooChi2Var('chi2/ndf', 'chi2/ndf', model_signal, hist))
    # text.Draw()
    can.SaveAs(f"../plots/fit/{year}/model_signal_{signal_mass}_{signal_region}.pdf")

    sig_model_dir = f'output/{year}/signal'
    if not os.path.exists(sig_model_dir):
        os.makedirs(sig_model_dir)
    f_out = ROOT.TFile(f"{sig_model_dir}/workspace_signal_{signal_mass}_{signal_region}.root", "RECREATE")
    w_sig = ROOT.RooWorkspace("workspace_signal", "workspace_signal")
    getattr(w_sig, "import")(model_signal)
    w_sig.Print()
    w_sig.Write()
    f_out.Close()

    with open(f'{sig_model_dir}/fit_info_signal_{signal_mass}_{signal_region}.yaml', 'w', encoding='utf-8') as f:
        info = {
            'x0': x0.getVal(),
            'sigmaL': sigmaL.getVal(),
            'sigmaR': sigmaR.getVal(),
            'alphaL': alphaL.getVal(),
            'alphaR': alphaR.getVal(),
            'nL': nL.getVal(),
            'nR': nR.getVal(),
            'event_sum': mc.sumEntries()
        }
        yaml.dump(info, f)


def background_fit(year):
    # # Background modelling
    f = ROOT.TFile(f"input/{year}/data_Hbb.root", "r")
    # Load TTree
    tree = f.Get("Events")

    # Define mass and weight variables
    mass_Zprime = ROOT.RooRealVar("mass_Zprime", "mass_Zprime", 1500, fit_range_low, fit_range_up)
    weight = ROOT.RooRealVar("weight", "weight", 1, -10, 10)
    mass_Higgs = ROOT.RooRealVar("mass_Higgs", "mass_Higgs", 125, 0, 999)
    tagger_Hbb = ROOT.RooRealVar("tagger_Hbb", "tagger_Hbb", 0, 0, 2)

    # Convert to RooDataSet
    data_sideband = ROOT.RooDataSet("data_sideband", "data_sideband", tree, ROOT.RooArgSet(mass_Zprime, weight, mass_Higgs, tagger_Hbb), sideband_cut, "weight")

    n_bins = (fit_range_up - fit_range_low) // 20
    binning = ROOT.RooFit.Binning(n_bins, fit_range_low, fit_range_up)

    # mass_Higgs.setBins(n_bins)
    # data_sideband_hist = ROOT.RooDataHist("data_sideband_hist", "data_sideband_hist", mass_Higgs, data_sideband)

    # Lets plot the signal mass distribution
    can = ROOT.TCanvas()
    plot = mass_Zprime.frame()
    data_sideband.plotOn(plot, binning)
    plot.Draw()
    can.Update()
    can.SaveAs(f"../plots/fit/{year}/data_sideband_mass_Zprime.pdf")

   
    p0 = ROOT.RooRealVar("p0", "p0", 1e3, 0, 1e6)
    p1 = ROOT.RooRealVar("p1", "p2", 1, -10, 10)
    p2 = ROOT.RooRealVar("p2", "p2", -1, -10, 10)
    model_background = ROOT.RooGenericPdf("model_background", "model_background", "TMath::Power(@0, @1 + @2 * TMath::Log(@0))", ROOT.RooArgList(mass_Zprime, p1, p2))
    """
    p0 = ROOT.RooRealVar("p0", "p0", 1e3, 0, 1e4)
    p1 = ROOT.RooRealVar("p1", "p1", -0.01, -1, 0)
    p2 = ROOT.RooRealVar("p2", "p2", -0.001, -0.1, 0)
    model_background = ROOT.RooGenericPdf("model_background", "model_background", "@0 * TMath::Power(@3, @1) * TMath::Exp(@2 * @3)", ROOT.RooArgList(p0, p1, p2, mass_Zprime))

    p0 = ROOT.RooRealVar("p0", "p0", 1e3, 0, 1e6)
    p1 = ROOT.RooRealVar("p1", "p1", -0.000001, -0.001, 0)
    p2 = ROOT.RooRealVar("p2", "p2", 10, 0, 2000)
    model_background = ROOT.RooGenericPdf("model_background", "model_background", "@0 * TMath::Power(1 + @1*@3, @2)", ROOT.RooArgList(p0, p1, p2, mass_Zprime))
    """
    # Fit model to data sidebands
    model_background.fitTo(data_sideband, ROOT.RooFit.SumW2Error(True))
    model_background.fitTo(data_sideband, ROOT.RooFit.SumW2Error(True))

    # Let's plot the model fit to the data
    can = ROOT.TCanvas()
    plot = mass_Zprime.frame()
    # We have to be careful with the normalisation as we only fit over sidebands
    # First do an invisible plot of the full data set
    data_sideband.plotOn(plot, binning, ROOT.RooFit.MarkerColor(0), ROOT.RooFit.LineColor(0))
    model_background.plotOn(plot, ROOT.RooFit.LineColor(2))
    data_sideband.plotOn(plot, binning)
    plot.Draw()
    can.Update()
    can.Draw()
    can.SaveAs(f"../plots/fit/{year}/model_background.pdf")

    background_norm = ROOT.RooRealVar("model_background_norm", "Number of background events", data_sideband.numEntries(), 0, 100 * data_sideband.numEntries())
    background_norm.setConstant(False)
    p0.setConstant(True)
    p1.setConstant(True)
    p2.setConstant(True)

    bkg_model_dir = f'output/{year}/background/dijet2'
    if not os.path.exists(bkg_model_dir):
        os.makedirs(bkg_model_dir)
    f_out = ROOT.TFile(f"{bkg_model_dir}/workspace_background_{signal_mass}_{signal_region}.root", "RECREATE")
    w_bkg = ROOT.RooWorkspace("workspace_background", "workspace_background")
    getattr(w_bkg, "import")(data_sideband)
    # getattr(w_bkg, "import")(data_sideband_hist)
    getattr(w_bkg, "import")(background_norm)
    getattr(w_bkg, "import")(model_background)
    w_bkg.Print()
    w_bkg.Write()
    f_out.Close()

    with open(f'{bkg_model_dir}/fit_info_background_{signal_mass}_{signal_region}.yaml', 'w', encoding='utf-8') as f:
        info = {
            'p0': p0.getVal(),
            'p1': p1.getVal(),
            'p2': p2.getVal(),
            'event_sum': data_sideband.sumEntries()
        }
        yaml.dump(info, f)


def get_SR_data(year):
    # # Data in SR
    f = ROOT.TFile(f"input/{year}/data_Hbb.root", "r")
    tree = f.Get("Events")

    # Define mass and weight variables
    mass_Zprime = ROOT.RooRealVar("mass_Zprime", "mass_Zprime", 1500, fit_range_low, fit_range_up)
    weight = ROOT.RooRealVar("weight", "weight", 0, -10, 10)
    mass_Higgs = ROOT.RooRealVar("mass_Higgs", "mass_Higgs", 125, 0, 500)
    tagger_Hbb = ROOT.RooRealVar("tagger_Hbb", "tagger_Hbb", 0, 0, 2)

    # Convert to RooDataSet
    data_SR = ROOT.RooDataSet("data_SR", "data_SR", tree, ROOT.RooArgSet(mass_Zprime, weight, mass_Higgs, tagger_Hbb), SR_cut, "weight")

    n_bins = (fit_range_up - fit_range_low) // 20
    binning = ROOT.RooFit.Binning(n_bins, fit_range_low, fit_range_up)

    # Lets plot the signal mass distribution
    can = ROOT.TCanvas()
    plot = mass_Zprime.frame()
    data_SR.plotOn(plot, binning)
    plot.Draw()
    can.Update()
    can.SaveAs(f"../plots/fit/{year}/data_SR_mass_Zprime.pdf")

    f_out = ROOT.TFile(f"./output/{year}/workspace_data_{signal_mass}.root", "RECREATE")
    w = ROOT.RooWorkspace("workspace_data", "workspace_data")
    getattr(w, "import")(data_SR)
    w.Print()
    w.Write()
    f_out.Close()


def get_SR_bkg_MC(year):
    # # Data in SR
    f = ROOT.TFile(f"input/{year}/mc_background.root", "r")
    tree = f.Get("Events")

    # Define mass and weight variables
    mass_Zprime = ROOT.RooRealVar("mass_Zprime", "mass_Zprime", 1500, fit_range_low, fit_range_up)
    weight = ROOT.RooRealVar("weight", "weight", 0, -999, 999)
    mass_Higgs = ROOT.RooRealVar("mass_Higgs", "mass_Higgs", 125, 0, 500)
    tagger_Hbb = ROOT.RooRealVar("tagger_Hbb", "tagger_Hbb", 0, 0, 2)

    # Convert to RooDataSet
    bkg_mc = ROOT.RooDataSet("bkg_mc", "bkg_mc", tree, ROOT.RooArgSet(mass_Zprime, weight, mass_Higgs, tagger_Hbb), SR_cut, "weight")

    n_bins = (fit_range_up - fit_range_low) // 20
    binning = ROOT.RooFit.Binning(n_bins, fit_range_low, fit_range_up)

    # mass_Zprime.setBins(n_bins)
    # bkg_mc_hist = ROOT.RooDataHist("data_sideband_hist", "data_sideband_hist", mass_Zprime, bkg_mc)

    # Lets plot the signal mass distribution
    can = ROOT.TCanvas()
    plot = mass_Zprime.frame()
    bkg_mc.plotOn(plot, binning)
    plot.Draw()
    can.Update()
    can.SaveAs(f"../plots/fit/{year}/bkg_mc_{signal_mass}_{signal_region}_mass_Zprime.pdf")

    if not os.path.exists(f'output/{year}/test_bkg_mc'):
        os.makedirs(f'output/{year}/test_bkg_mc')
    f_out = ROOT.TFile(f"./output/{year}/test_bkg_mc/workspace_bkg_mc_{signal_mass}_{signal_region}.root", "RECREATE")
    w = ROOT.RooWorkspace("workspace_bkg_mc", "workspace_bkg_mc")
    getattr(w, "import")(bkg_mc)
    w.Print()
    w.Write()
    f_out.Close()


if __name__ == "__main__":
    args = parse_commandline()
    year = args.year
    signal_mass = args.signal_mass
    signal_region = args.SR

    SR_binning = {
        'SR1': (0.7, 0.9),
        'SR2': (0.9, 2)
    }
    tagger_cut_low, tagger_cut_high = SR_binning[signal_region]
    fit_range_low, fit_range_up = args.fit_range_low, args.fit_range_up

    Fit_signal = True
    Fit_background = True
    SR_cut = f"(mass_Higgs>110) & (mass_Higgs<140) & (tagger_Hbb>{tagger_cut_low}) & (tagger_Hbb<{tagger_cut_high})"
    sideband_cut = f"(((mass_Higgs>50) & (mass_Higgs<70)) | ((mass_Higgs>100) & (mass_Higgs<110)) | (mass_Higgs>140)) & (tagger_Hbb>{tagger_cut_low}) & (tagger_Hbb<{tagger_cut_high})"

    if Fit_signal:
        fit_signal(year)
    if Fit_background:
        background_fit(year)

    get_SR_bkg_MC(year)
