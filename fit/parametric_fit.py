#!/usr/bin/env python3
import ROOT, os, yaml, argparse, uproot
import numpy as np
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)


def parse_commandline():
    parser = argparse.ArgumentParser(description='parametric fitting')
    parser.add_argument('-y', '--year', help='To specify which year', choices=('2016pre', '2016post', '2016', '2017', '2018', 'Run2'), default='Run2')
    parser.add_argument('-m', '--signal_mass', help='To specify the mass of signal resonance', type=int, default=None)
    parser.add_argument('-R', '--SR', help='To specify which signal region', choices=('SR1', 'SR2', None), default=None)
    parser.add_argument('-d', '--fit_range_down', help='To specify the lower bound of fitting range', default=650, type=int)
    parser.add_argument('-u', '--fit_range_up', help='To specify the higher bound of fitting range', default=4000, type=int)
    args = parser.parse_args()
    return args


def fit_signal(year, jet, signal_mass, region, cut):
    signal_region = region[:2]+jet+region[2:]
    m = int(str(signal_mass).split('_')[0])

    with open('../src/parameters/uncertainty/systematics.yaml', 'r', encoding='utf-8') as f:
        systematics = yaml.safe_load(f)

    if '_' in str(signal_mass):
        sigma = m * (0.056 if '_5p6' in signal_mass else 0.10)
    else:
        f = uproot.open(f"input/{year}/{signal_mass}/bbgamma_SR{jet}.root")
        sigma = np.std(f['Events']['fit_mass'].array())

    ## Signal modeling
    f = ROOT.TFile(f"input/{year}/{signal_mass}/bbgamma_SR{jet}.root", "r")
    # Load TTree
    tree = f.Get("Events")

    # Define mass and weight variables
    fit_mass = ROOT.RooRealVar("fit_mass", "fit_mass", m, max(fit_range_down, m-5*sigma), min(fit_range_up, m+5*sigma))
    weight = ROOT.RooRealVar("weight", "weight", 0.1, 0, 100)
    jet_mass = ROOT.RooRealVar("jet_mass", "jet_mass", 125, 0, 999)
    tagger = ROOT.RooRealVar("tagger", "tagger", 0.5, 0, 2)

    # Convert to RooDataSet

    mc = ROOT.RooDataSet("signal", "signal", tree, ROOT.RooArgSet(fit_mass, weight, jet_mass, tagger), cut, "weight")

    # Lets plot the signal mass distribution
    can = ROOT.TCanvas()
    plot = fit_mass.frame()
    mc.plotOn(plot)
    plot.Draw()
    can.Update()
    os.makedirs(f'../plots/fit/{year}/{signal_mass}', exist_ok=True)
    can.SaveAs(f"../plots/fit/{year}/{signal_mass}/fit_variable_{signal_region}_{signal_mass}.pdf")

    # Introduce RooRealVars into the workspace for the fitted variable
    x0 = ROOT.RooRealVar("x0", "x0", m, m - 300, m + 300)
    sigmaL = ROOT.RooRealVar("sigmaL", "sigmaL", sigma, 5, 5*sigma)
    sigmaR = ROOT.RooRealVar("sigmaR", "sigmaR", sigma, 5, 5*sigma)
    alphaL = ROOT.RooRealVar("alphaL", "alphaL", 1, 0.01, 5)
    alphaR = ROOT.RooRealVar("alphaR", "alphaR", 1, 0.01, 5)
    nL = ROOT.RooRealVar("nL", "nL", 1, 0.1, 5)
    nR = ROOT.RooRealVar("nR", "nR", 1, 0.1, 5)

    if year != 'Run2':
        JES = ROOT.RooRealVar(f"JES_{year}", f"JES_{year}", 0, -5, 5)
    else:
        JES = ROOT.RooRealVar("JES", f"JES", 0, -5, 5)
    JER = ROOT.RooRealVar("JER", "JER", 0, -5, 5)
    PES = ROOT.RooRealVar("PES", "PES", 0, -5, 5)
    PER = ROOT.RooRealVar("PER", "PER", 0, -5, 5)
    JES.setConstant(True); JER.setConstant(True); PES.setConstant(True); PER.setConstant(True)
    mean = ROOT.RooFormulaVar("mean", "mean", 
        "@0*(1+%f*@1+%f*@2)"%(systematics['JES'][region][m]-1, (systematics['PES'][region][m]-1)), 
        ROOT.RooArgList(x0, JES, PES))
    widthL = ROOT.RooFormulaVar("widthL", "widthL", 
        "@0*(1+%f*@1+%f*@2)"%(systematics['JER'][region][jet][signal_mass]-1, systematics['PER'][region][jet][signal_mass]-1), 
        ROOT.RooArgList(sigmaL, JER, PER))
    widthR = ROOT.RooFormulaVar("widthR", "widthR", 
        "@0*(1+%f*@1+%f*@2)"%(systematics['JER'][region][jet][signal_mass]-1, systematics['PER'][region][jet][signal_mass]-1),
        ROOT.RooArgList(sigmaR, JER, PER))

    # Define the Gaussian with mean=MH and width=sigma
    model_signal = ROOT.RooCrystalBall(f"model_bbgamma_{signal_region}", f"model_bbgamma_{signal_region}", fit_mass, mean, widthL, widthR, alphaL, nL, alphaR, nR)
    signal_norm = ROOT.RooRealVar(f"model_bbgamma_{signal_region}_norm", f"Number of signal events in {signal_region}", mc.sumEntries(), 1e-3, 5*mc.sumEntries())

    # Fit Gaussian to MC events and plot
    model_signal.fitTo(mc, ROOT.RooFit.SumW2Error(True))

    x0.setConstant(True)
    sigmaL.setConstant(True)
    sigmaR.setConstant(True)
    alphaL.setConstant(True)
    alphaR.setConstant(True)
    nL.setConstant(True)
    nR.setConstant(True)
    #signal_norm.setConstant(True)

    sig_model_dir = f'workspace/{year}/{signal_mass}'
    os.makedirs(sig_model_dir, exist_ok=True)
    f_out = ROOT.TFile(f"{sig_model_dir}/signal_{signal_region}.root", "RECREATE")
    w_sig = ROOT.RooWorkspace("workspace_signal", "workspace_signal")
    getattr(w_sig, "import")(model_signal)
    #getattr(w_sig, "import")(signal_norm)
    w_sig.Print()
    w_sig.Write()
    f_out.Close()

    with open(f'{sig_model_dir}/signal_{signal_region}.yaml', 'w', encoding='utf-8') as f:
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
            'norm': signal_norm.getVal(),
            'sigma': float(sigma)
        }
        yaml.dump(info, f)


def fit_background(year, jet, region, cut):
    region = region[:2]+jet+region[2:]
    bkg_model_dir = f'workspace/{year}'
    
    # # Background modelling
    f = ROOT.TFile(f"input/{year}/data.root", "r")
    # Load TTree
    tree = f.Get("Events")

    # Define mass and weight variables
    fit_mass = ROOT.RooRealVar("fit_mass", "fit_mass", 1500, fit_range_down, fit_range_up)
    weight = ROOT.RooRealVar("weight", "weight", 1, -10, 10)
    jet_mass = ROOT.RooRealVar("jet_mass", "jet_mass", 125, 0, 999)
    tagger = ROOT.RooRealVar("tagger", "tagger", 0.5, 0, 2)

    # Convert to RooDataSet
    data_region = ROOT.RooDataSet(f"data_{region}", f"data_{region}", tree, ROOT.RooArgSet(fit_mass, weight, jet_mass, tagger), cut, "weight")

    ## Multiple background models
    model, p1, p2, p3 = {}, {}, {}, {}
    energy = 1e2

    # expow1 model
    p1['expow1'] = ROOT.RooRealVar("p1_expow1", "p1_expow1", -1, -10, 0)
    model['expow1'] = ROOT.RooGenericPdf("model_background_expow1", "model_background_expow1", f"TMath::Power(@0/{energy}, @1)", ROOT.RooArgList(fit_mass, p1['expow1']))

    # expow2 model
    p1['expow2'] = ROOT.RooRealVar("p1_expow2", "p1_expow2", -1, -10, 0)
    p2['expow2'] = ROOT.RooRealVar("p2_expow2", "p2_expow2", -1e-2, -0.5, 0.5)
    model['expow2'] = ROOT.RooGenericPdf("model_background_expow2", "model_background_expow2", f"TMath::Power(@0/{energy}, @1) * TMath::Exp(@2 * @0/{energy})", ROOT.RooArgList(fit_mass, p1['expow2'], p2['expow2']))

    # dijet2 model
    p1['dijet2'] = ROOT.RooRealVar("p1_dijet2", "p1_dijet2", -2, -10, 0)
    p2['dijet2'] = ROOT.RooRealVar("p2_dijet2", "p2_dijet2", -1, -5, 0)
    model['dijet2'] = ROOT.RooGenericPdf("model_background_dijet2", "model_background_dijet2", f"TMath::Power(@0/{energy}, @1 + @2 * TMath::Log(@0/{energy}))", ROOT.RooArgList(fit_mass, p1['dijet2'], p2['dijet2']))

    # dijet3 model
    p1['dijet3'] = ROOT.RooRealVar("p1_dijet3", "p1_dijet3", -1, -10, 0)
    p2['dijet3'] = ROOT.RooRealVar("p2_dijet3", "p2_dijet3", -1, -5, 0)
    p3['dijet3'] = ROOT.RooRealVar("p3_dijet3", "p3_dijet3", -1e-3, -0.1, 0.1)
    model['dijet3'] = ROOT.RooGenericPdf("model_background_dijet3", "model_background_dijet3", f"TMath::Power(@0/{energy}, @1 + @2 * TMath::Log(@0/{energy}) + @3 * TMath::Power(TMath::Log(@0/{energy}), 2))", ROOT.RooArgList(fit_mass, p1['dijet3'], p2['dijet3'], p3['dijet3']))

    # invpow2 model
    p1['invpow2'] = ROOT.RooRealVar("p1_invpow2", "p1_invpow2", 1e-2, 0, 10)
    p2['invpow2'] = ROOT.RooRealVar("p2_invpow2", "p2_invpow2", -2, -30, 0)
    model['invpow2'] = ROOT.RooGenericPdf("model_background_invpow2", "model_background_invpow2", f"TMath::Power(1 + @1*@0/{energy}, @2)", ROOT.RooArgList(fit_mass, p1['invpow2'], p2['invpow2']))

    # invpow3 model
    p1['invpow3'] = ROOT.RooRealVar("p1_invpow3", "p1_invpow3", 1e-2, 0, 10)
    p2['invpow3'] = ROOT.RooRealVar("p2_invpow3", "p2_invpow3", -2, -30, 0)
    p3['invpow3'] = ROOT.RooRealVar("p3_invpow3", "p3_invpow3", -0.5, -1, 1)
    model['invpow3'] = ROOT.RooGenericPdf("model_background_invpow3", "model_background_invpow3", f"TMath::Power(1 + @1*@0/{energy}, @2 + @3*@0/{energy})", ROOT.RooArgList(fit_mass, p1['invpow3'], p2['invpow3'], p3['invpow3']))

    # Make a RooCategory object: this will control which PDF is "active"
    category = ROOT.RooCategory(f"pdfindex_{region}", "Index of Pdf which is active")
    # Make a RooArgList of the models
    models = ROOT.RooArgList()

    # Fit models
    for k in ['expow1', 'expow2', 'dijet2', 'dijet3', 'invpow2', 'invpow3']:
        model[k].fitTo(data_region, ROOT.RooFit.SumW2Error(True))
        p1[k].setConstant(True)
        if k in p2:
            p2[k].setConstant(True)
        if k in p3:
            p3[k].setConstant(True)
        models.add(model[k])
    
    # Build the RooMultiPdf object
    multipdf = ROOT.RooMultiPdf(f"multipdf_{region}", f"multipdf_{region}", category, models)

    background_norm = ROOT.RooRealVar(f"multipdf_{region}_norm", "Number of background events", data_region.sumEntries(), 0, 100 * data_region.sumEntries())
    background_norm.setConstant(False)
    os.makedirs(bkg_model_dir, exist_ok=True)
    with open(f'{bkg_model_dir}/background_{region}.yaml', 'w', encoding='utf-8') as f:
        info = {
            'p1': {k: p1[k].getVal() for k in p1},
            'p2': {k: p2[k].getVal() for k in p2},
            'p3': {k: p3[k].getVal() for k in p3},
            'event_sum': data_region.sumEntries(),
            'norm': background_norm.getVal()
        }
        yaml.dump(info, f)


    #fit_mass.setBins(320)
    #data_hist = ROOT.RooDataHist(f"data_{region}_hist", f"data_{region}_hist", fit_mass, data_region)
    os.makedirs(bkg_model_dir, exist_ok=True)
    f_out = ROOT.TFile(f"{bkg_model_dir}/data_{region}.root", "RECREATE")
    w_bkg = ROOT.RooWorkspace(f"workspace_{region}", f"workspace_{region}")
    getattr(w_bkg, "import")(data_region)
    getattr(w_bkg, "import")(multipdf)
    for k in model:
        getattr(w_bkg, "import")(model[k])
    getattr(w_bkg, "import")(background_norm)
    #getattr(w_bkg, "import")(data_hist)
    w_bkg.Print()
    w_bkg.Write()
    f_out.Close()


def get_SR_data(year, region, cut, jet):
    # # Data in SR
    f = ROOT.TFile(f"input/{year}/data_CR_to_SR{jet}.root", "r")
    #f = ROOT.TFile(f"input/{year}/background_mc.root", "r")
    tree = f.Get("Events")

    # Define mass and weight variables
    fit_mass = ROOT.RooRealVar("fit_mass", "fit_mass", 1500, fit_range_down, fit_range_up)
    weight = ROOT.RooRealVar("weight", "weight", 1, 0, 10)
    jet_mass = ROOT.RooRealVar("jet_mass", "jet_mass", 125, 0, 500)
    tagger = ROOT.RooRealVar("tagger", "tagger", 0, 0, 2)

    # Convert to RooDataSet
    data_region = ROOT.RooDataSet(f"data_{region}", f"data_{region}", tree, ROOT.RooArgSet(fit_mass, weight, jet_mass, tagger), cut, "weight")

    data_dir = f"./workspace/{year}"
    os.makedirs(data_dir, exist_ok=True)
    f_out = ROOT.TFile(f"{data_dir}/data_{jet}bb_{region}.root", "RECREATE")
    w = ROOT.RooWorkspace(f"workspace_{region}", f"workspace_{region}")
    getattr(w, "import")(data_region)
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
    mass_cut = {
        'Z': [80, 110],
        'H': [110, 150],
    }
    for SR in signal_region:
        CR = SR.replace('S', 'C')
        tagger_cut_low, tagger_cut_high = tagger_cut[SR]
        CR_cut = f"""(
            (((jet_mass>50) & (jet_mass<{mass_cut['Z'][0]})) | (jet_mass>{mass_cut['H'][1]})) & 
            (tagger>{tagger_cut_low}) & (tagger<{tagger_cut_high})
        )"""

        for jet in ['H', 'Z']:
            mass_low, mass_high = mass_cut[jet]
            SR_cut = f"""(
                (jet_mass>{mass_low}) & (jet_mass<{mass_high}) & 
                (tagger>{tagger_cut_low}) & (tagger<{tagger_cut_high})
            )"""
            #get_SR_data(year, SR, SR_cut, jet)
            if Fit_background:
                fit_background(year, jet, SR, SR_cut)

            if Fit_signal:
                for m in signal_mass:
                    fit_signal(year, jet, m, SR, SR_cut)
                    if jet == 'Z':
                        fit_signal(year, jet, str(m)+'_5p6', SR, SR_cut)
                        fit_signal(year, jet, str(m)+'_10p0', SR, SR_cut)
