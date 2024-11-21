#!/usr/bin/env python3
import ROOT, os, yaml, argparse, uproot
import numpy as np
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)


def parse_commandline():
    parser = argparse.ArgumentParser(description='parametric fitting')
    parser.add_argument('-y', '--year', help='To specify which year', choices=('2016pre', '2016post', '2017', '2018', 'Run2'), default='Run2')
    parser.add_argument('-m', '--signal_mass', help='To specify the mass of signal resonance', type=int, default=None)
    parser.add_argument('-R', '--SR', help='To specify which signal region', choices=('SR1', 'SR2', None), default=None)
    parser.add_argument('-d', '--fit_range_down', help='To specify the lower bound of fitting range', default=650, type=int)
    parser.add_argument('-u', '--fit_range_up', help='To specify the higher bound of fitting range', default=4000, type=int)
    args = parser.parse_args()
    return args
                     

def fit_signal(year, fatjet, signal_mass, SR):
    with open('../src/parameters/uncertainty/shape_uncertainties.yaml', 'r', encoding='utf-8') as f:
        shape_uncertainties = yaml.safe_load(f)
    
    f = uproot.open(f"input/{year}/{signal_mass}/{fatjet}bb_gamma.root")
    sigma = np.std(f['Events']['fit_mass'].array())
    
    # # Signal modelling
    f = ROOT.TFile(f"input/{year}/{signal_mass}/{fatjet}bb_gamma.root", "r")
    # Load TTree
    tree = f.Get("Events")

    # Define mass and weight variables
    fit_mass = ROOT.RooRealVar("fit_mass", "fit_mass", signal_mass, signal_mass-5*sigma, signal_mass+5*sigma)
    weight = ROOT.RooRealVar("weight", "weight", 0.1, 0, 100)
    jet_mass = ROOT.RooRealVar("jet_mass", "jet_mass", 125, 0, 999)
    tagger = ROOT.RooRealVar("tagger", "tagger", 0, 0, 2)

    # Convert to RooDataSet

    mc = ROOT.RooDataSet("signal", "signal", tree, ROOT.RooArgSet(fit_mass, weight, jet_mass, tagger), SR_cut, "weight")

    # Lets plot the signal mass distribution
    can = ROOT.TCanvas()
    plot = fit_mass.frame()
    mc.plotOn(plot)
    plot.Draw()
    can.Update()
    if not os.path.exists(f'../plots/fit/{year}/{signal_mass}'):
        os.makedirs(f'../plots/fit/{year}/{signal_mass}')
    can.SaveAs(f"../plots/fit/{year}/{signal_mass}/fit_variable_{fatjet}bb_{signal_mass}_{SR}.pdf")

    # Introduce RooRealVars into the workspace for the fitted variable
    x0 = ROOT.RooRealVar("x0", "x0", signal_mass, signal_mass - 200, signal_mass + 200)
    sigmaL = ROOT.RooRealVar("sigmaL", "sigmaL", sigma, 0, 5*sigma)
    sigmaR = ROOT.RooRealVar("sigmaR", "sigmaR", sigma, 0, 5*sigma)
    alphaL = ROOT.RooRealVar("alphaL", "alphaL", 1, 0.1, 5)
    alphaR = ROOT.RooRealVar("alphaR", "alphaR", 1, 0.1, 5)
    nL = ROOT.RooRealVar("nL", "nL", 1, 0.2, 3)
    nR = ROOT.RooRealVar("nR", "nR", 1, 0.2, 3)

    JES = ROOT.RooRealVar("JES", "JES", 0, -5, 5)
    JER = ROOT.RooRealVar("JER", "JER", 0, -5, 5)
    PU = ROOT.RooRealVar("PU", "PU", 0, -5, 5)
    JES.setConstant(True); JER.setConstant(True); PU.setConstant(True);
    mean = ROOT.RooFormulaVar("mean", "mean", 
        "@0*(1+%f*@1+%f*@2+%f*@3)"%tuple(shape_uncertainties[SR][signal_mass]['x0'][k] for k in ['JES', 'JER', 'PU']), 
        ROOT.RooArgList(x0, JES, JER, PU))
    widthL = ROOT.RooFormulaVar("widthL", "widthL", 
        "@0*(1+%f*@1+%f*@2+%f*@3)"%tuple(shape_uncertainties[SR][signal_mass]['sigmaL'][k] for k in ['JES', 'JER', 'PU']), 
        ROOT.RooArgList(sigmaL, JES, JER, PU))
    widthR = ROOT.RooFormulaVar("widthR", "widthR", 
        "@0*(1+%f*@1+%f*@2+%f*@3)"%tuple(shape_uncertainties[SR][signal_mass]['sigmaR'][k] for k in ['JES', 'JER', 'PU']), 
        ROOT.RooArgList(sigmaR, JES, JER, PU))

    # Define the Gaussian with mean=MH and width=sigma
    model_signal = ROOT.RooCrystalBall("model_bbgamma", "model_bbgamma", fit_mass, mean, widthL, widthR, alphaL, nL, alphaR, nR)
    signal_norm = ROOT.RooRealVar("model_bbgamma_norm", f"Number of signal events in Tag {fatjet}bb+gamma", mc.sumEntries(), 0, 100*mc.sumEntries())

    # Fit Gaussian to MC events and plot
    model_signal.fitTo(mc, ROOT.RooFit.SumW2Error(True))

    x0.setConstant(True)
    sigmaL.setConstant(True)
    sigmaR.setConstant(True)
    alphaL.setConstant(True)
    alphaR.setConstant(True)
    nL.setConstant(True)
    nR.setConstant(True)
    signal_norm.setConstant(True)

    can = ROOT.TCanvas()
    plot = fit_mass.frame()
    mc.plotOn(plot)
    model_signal.plotOn(plot, ROOT.RooFit.LineColor(2))
    plot.Draw()
    can.Update()
    can.Draw()

    fit_mass.setBins(160)
    # hist = ROOT.RooDataHist("hist", "hist", fit_mass, mc)
    # print("==> chi^2/ndf = ", ROOT.RooChi2Var('chi2/ndf', 'chi2/ndf', model_signal, hist))
    # text.Draw()
    can.SaveAs(f"../plots/fit/{year}/model_{fatjet}bb_{signal_mass}_{SR}.pdf")

    sig_model_dir = f'workspace/{year}/{signal_mass}'
    if not os.path.exists(sig_model_dir):
        os.makedirs(sig_model_dir)
    f_out = ROOT.TFile(f"{sig_model_dir}/{fatjet}bb_{signal_mass}_{SR}.root", "RECREATE")
    w_sig = ROOT.RooWorkspace("workspace_signal", "workspace_signal")
    getattr(w_sig, "import")(model_signal)
    getattr(w_sig, "import")(signal_norm)
    w_sig.Print()
    w_sig.Write()
    f_out.Close()

    with open(f'{sig_model_dir}/fit_info_{fatjet}bb_{SR}.yaml', 'w', encoding='utf-8') as f:
        info = {
            'x0': x0.getVal(),
            'mean': mean.getVal(),
            'sigmaL': sigmaL.getVal(),
            'widthL': widthL.getVal(),
            'sigmaR': sigmaR.getVal(),
            'widthR': widthR.getVal(),
            'alphaL': alphaL.getVal(),
            'alphaR': alphaR.getVal(),
            'nL': nL.getVal(),
            'nR': nR.getVal(),
            'event_sum': mc.sumEntries(),
            'norm': signal_norm.getVal()
        }
        yaml.dump(info, f)


def fit_background(year, CR):
    bkg_model_dir = f'workspace/{year}'
    
    # # Background modelling
    f = ROOT.TFile(f"input/{year}/data.root", "r")
    # Load TTree
    tree = f.Get("Events")

    # Define mass and weight variables
    fit_mass = ROOT.RooRealVar("fit_mass", "fit_mass", 1500, fit_range_down, fit_range_up)
    weight = ROOT.RooRealVar("weight", "weight", 1, -10, 10)
    jet_mass = ROOT.RooRealVar("jet_mass", "jet_mass", 125, 0, 999)
    tagger = ROOT.RooRealVar("tagger", "tagger", 0, 0, 2)

    # Convert to RooDataSet
    data_CR = ROOT.RooDataSet("data_CR", "data_CR", tree, ROOT.RooArgSet(fit_mass, weight, jet_mass, tagger), CR_cut, "weight")


    n_bins = (fit_range_up - fit_range_down) // 20
    binning = ROOT.RooFit.Binning(n_bins, fit_range_down, fit_range_up)

    # Lets plot the signal mass distribution
    can = ROOT.TCanvas()
    plot = fit_mass.frame()
    data_CR.plotOn(plot, binning)
    plot.Draw()
    can.Update()
    can.SaveAs(f"../plots/fit/{year}/data_CR_fit_mass.pdf")

    ## Multiple background models
    model, p1, p2, p3 = {}, {}, {}, {}

    # dijet2 model
    p1['dijet2'] = ROOT.RooRealVar("p1", "p1", 1, -10, 100)
    p2['dijet2'] = ROOT.RooRealVar("p2", "p2", -1, -10, 10)
    model['dijet2'] = ROOT.RooGenericPdf("model_background_dijet2", "model_background_dijet2", "TMath::Power(@0, @1 + @2 * TMath::Log(@0))", ROOT.RooArgList(fit_mass, p1['dijet2'], p2['dijet2']))

    # dijet3 model
    p1['dijet3'] = ROOT.RooRealVar("p1", "p1", 1, -10, 10)
    p2['dijet3'] = ROOT.RooRealVar("p2", "p2", -1, -10, 10)
    p3['dijet3'] = ROOT.RooRealVar("p3", "p3", -0.1, -10, 10)
    model['dijet3'] = ROOT.RooGenericPdf("model_background_dijet3", "model_background_dijet3", "TMath::Power(@0, @1 + @2 * TMath::Log(@0) + @3 * TMath::Power(TMath::Log(@0), 2))", ROOT.RooArgList(fit_mass, p1['dijet3'], p2['dijet3'], p3['dijet3']))

    # expow1 model
    p1['expow1'] = ROOT.RooRealVar("p1", "p1", -0.1, -10, 0)
    model['expow1'] = ROOT.RooGenericPdf("model_background_expow1", "model_background_expow1", "TMath::Power(@0, @1)", ROOT.RooArgList(fit_mass, p1['expow1']))
    
    # expow2 model
    p1['expow2'] = ROOT.RooRealVar("p1", "p1", -0.01, -5, 0)
    p2['expow2'] = ROOT.RooRealVar("p2", "p2", -0.001, -0.1, 0)
    model['expow2'] = ROOT.RooGenericPdf("model_background_expow2", "model_background_expow2", "TMath::Power(@0, @1) * TMath::Exp(@2 * @0)", ROOT.RooArgList(fit_mass, p1['expow2'], p2['expow2']))

    # invpow2 model
    p1['invpow2'] = ROOT.RooRealVar("p1", "p1", -0.000001, -0.001, 0)
    p2['invpow2'] = ROOT.RooRealVar("p2", "p2", 10, 0, 2000)
    model['invpow2'] = ROOT.RooGenericPdf("model_background_invpow2", "model_background_invpow2", "TMath::Power(1 + @1*@0, @2)", ROOT.RooArgList(fit_mass, p1['invpow2'], p2['invpow2']))

    # invpow3 model
    p1['invpow3'] = ROOT.RooRealVar("p1", "p1", -0.000001, -0.001, 0)
    p2['invpow3'] = ROOT.RooRealVar("p2", "p2", 10, 0, 2000)
    p3['invpow3'] = ROOT.RooRealVar("p3", "p3", -0.1, -1, 10)
    model['invpow3'] = ROOT.RooGenericPdf("model_background_invpow3", "model_background_invpow3", "TMath::Power(1 + @1*@0, @2 + @3*@0)", ROOT.RooArgList(fit_mass, p1['invpow3'], p2['invpow3'], p3['invpow3']))

    # Make a RooCategory object: this will control which PDF is "active"
    category = ROOT.RooCategory(f"pdfindex_{CR}", "Index of Pdf which is active")

    # Make a RooArgList of the models
    models = ROOT.RooArgList()

    # Fit model to data sidebands
    for k in model:
        model[k].fitTo(data_CR, ROOT.RooFit.SumW2Error(True))
        p1[k].setConstant(True)
        if k in p2:
            p2[k].setConstant(True)
        if k in p3:
            p3[k].setConstant(True)
        models.add(model[k])

    # Build the RooMultiPdf object
    multipdf = ROOT.RooMultiPdf(f"multipdf", "MultiPdf", category, models)
    background_norm = ROOT.RooRealVar(f"multipdf_norm", "Number of background events", data_CR.numEntries(), 0, 100 * data_CR.numEntries())
    background_norm.setConstant(False)

    # Let's plot the model fit to the data
    can = ROOT.TCanvas()
    fit_mass.setBins(100)
    plot = fit_mass.frame()
    # We have to be careful with the normalisation as we only fit over sidebands
    # First do an invisible plot of the full data set
    data_CR.plotOn(plot, binning, ROOT.RooFit.MarkerColor(0), ROOT.RooFit.LineColor(0))
    model['dijet2'].plotOn(plot, ROOT.RooFit.LineColor(2))
    data_CR.plotOn(plot, binning)
    plot.Draw()
    can.Update()
    can.Draw()
    can.SaveAs(f"../plots/fit/{year}/background_{CR}.pdf")


    if not os.path.exists(bkg_model_dir):
        os.makedirs(bkg_model_dir)
    with open(f'{bkg_model_dir}/fit_info_background_{CR}.yaml', 'w', encoding='utf-8') as f:
        info = {
            'p1': {func: p1[func].getVal() for func in p1},
            'p2': {func: p2[func].getVal() for func in p2},
            'p3': {func: p3[func].getVal() for func in p3},
            'CR_num': data_CR.sumEntries(),
            'norm': background_norm.getVal()
        }
        yaml.dump(info, f)

    if not os.path.exists(bkg_model_dir):
        os.makedirs(bkg_model_dir)
    f_out = ROOT.TFile(f"{bkg_model_dir}/data_{CR}.root", "RECREATE")
    w_bkg = ROOT.RooWorkspace("workspace_CR", "workspace_CR")
    getattr(w_bkg, "import")(data_CR)
    getattr(w_bkg, "import")(category)
    getattr(w_bkg, "import")(background_norm)
    getattr(w_bkg, "import")(multipdf)
    w_bkg.Print()
    w_bkg.Write()
    f_out.Close()


def get_SR_data(year, SR):
    # # Data in SR
    f = ROOT.TFile(f"input/{year}/data.root", "r")
    tree = f.Get("Events")

    # Define mass and weight variables
    fit_mass = ROOT.RooRealVar("fit_mass", "fit_mass", 1500, fit_range_down, fit_range_up)
    weight = ROOT.RooRealVar("weight", "weight", 0, -10, 10)
    jet_mass = ROOT.RooRealVar("jet_mass", "jet_mass", 125, 0, 500)
    tagger = ROOT.RooRealVar("tagger", "tagger", 0, 0, 2)

    # Convert to RooDataSet
    data_SR = ROOT.RooDataSet("data_SR", "data_SR", tree, ROOT.RooArgSet(fit_mass, weight, jet_mass, tagger), SR_cut, "weight")

    n_bins = (fit_range_up - fit_range_down) // 20
    binning = ROOT.RooFit.Binning(n_bins, fit_range_down, fit_range_up)

    # Lets plot the signal mass distribution
    can = ROOT.TCanvas()
    plot = fit_mass.frame()
    data_SR.plotOn(plot, binning)
    plot.Draw()
    can.Update()
    can.SaveAs(f"../plots/fit/{year}/data_{SR}_fit_mass.pdf")

    data_dir = f"./workspace/{year}"
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    f_out = ROOT.TFile(f"{data_dir}/data_{SR}_{fatjet}bb.root", "RECREATE")
    w = ROOT.RooWorkspace("workspace_SR", "workspace_SR")
    getattr(w, "import")(data_SR)
    w.Print()
    w.Write()
    f_out.Close()


if __name__ == "__main__":
    args = parse_commandline()
    year = args.year
    fit_range_down, fit_range_up = args.fit_range_down, args.fit_range_up                    
    signal_mass = [700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 3000, 3500]
    
    if args.signal_mass is not None:
        signal_mass = [args.signal_mass]
    if args.SR is not None:
        signal_region = [args.SR]
    else:
        signal_region = ['SR1', 'SR2']

    Fit_signal = True
    Fit_background = True

    tagger_cut = {
        'SR1': [0.8, 0.98],
        'SR2': [0.98, 2],
    }
    mass_SR = {
        'Z': [80, 110],
        'H': [110, 150],
    }
    for SR in signal_region:
        CR = SR.replace('S', 'C')
        tagger_cut_low, tagger_cut_high = tagger_cut[SR]
        CR_cut = f"""(
            (((jet_mass>50) & (jet_mass<{mass_SR['Z'][0]})) | (jet_mass>{mass_SR['H'][1]})) & 
            (tagger>{tagger_cut_low}) & (tagger<{tagger_cut_high})
        )"""
        if Fit_background:
            fit_background(year, CR)

        for fatjet in ['H', 'Z']:
            mass_low, mass_high = mass_SR[fatjet]
            SR_cut = f"""(
                (jet_mass>{mass_low}) & (jet_mass<{mass_high}) & 
                (tagger>{tagger_cut_low}) & (tagger<{tagger_cut_high})
            )"""
            get_SR_data(year, SR)
            
            for m in signal_mass:
                if Fit_signal:
                    fit_signal(year, fatjet, m, SR)
