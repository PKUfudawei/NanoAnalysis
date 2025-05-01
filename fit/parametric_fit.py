#!/usr/bin/env python3
import ROOT, os, yaml, argparse, uproot
import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)


def parse_commandline():
    parser = argparse.ArgumentParser(description='parametric fitting')
    parser.add_argument('-i', '--in_dir', help='To specify the input directory', type=str, default='./slimmed_ntuple/')
    parser.add_argument('-o', '--out_dir', help='To specify the input directory', type=str, default='./workspace/')
    parser.add_argument('-y', '--year', help='To specify which year', choices=('2016pre', '2016post', '2016', '2017', '2018', 'Run2'), default='Run2')
    parser.add_argument('-m', '--signal_mass', help='To specify the mass of signal resonance', type=int, default=None)
    parser.add_argument('-R', '--signal_region', help='To specify which signal region', choices=('SRH1_N', 'SRH2_N', 'SRZ1_N', 'SRZ2_N', 'SRZ1_W', 'SRZ2_W', 'SRZ1_VW', 'SRZ2_VW', 'CR1', 'CR2', None), default=None)
    parser.add_argument('-d', '--fit_range_low', help='To specify the lower bound of fitting range', default=650, type=int)
    parser.add_argument('-u', '--fit_range_high', help='To specify the higher bound of fitting range', default=4000, type=int)
    args = parser.parse_args()
    return args


def Garwood_eror(number, direction):
    upper = np.array([1.84, 3.30, 4.64, 5.92, 7.16, 8.38, 9.58, 10.77, 11.95, 13.11, 14.27])
    lower = np.array([0, 0.17, 0.71, 1.37, 2.09, 2.84, 3.62, 4.42, 5.23, 6.06, 6.89])
    center = np.arange(11)
    if number < 0:
        raise ValueError('number < 0 !')
    elif number > 10:
        return np.sqrt(number)
    elif direction=='up':
        return (upper-center)[int(number)]
    elif direction=='down':
        return (center-lower)[int(number)]
    else:
        return None


def get_signal_norm(file_path, mass, cut, x_min, x_max):
    ## Signal modeling
    f = ROOT.TFile(file_path, "r")
    # Load TTree
    tree = f.Get("Events")

    # Define mass and weight variables
    fit_mass = ROOT.RooRealVar("fit_mass", "fit_mass", mass, x_min, x_max)
    weight = ROOT.RooRealVar("weight", "weight", 0.1, 0, 100)
    jet_mass = ROOT.RooRealVar("jet_mass", "jet_mass", 125, 0, 999)
    tagger = ROOT.RooRealVar("tagger", "tagger", 0.5, 0, 2)

    # Convert to RooDataSet
    mc = ROOT.RooDataSet("signal", "signal", tree, ROOT.RooArgSet(fit_mass, weight, jet_mass, tagger), cut, "weight")

    return mc.sumEntries() #mc.sumEntries('', 'fit_range')


def plot_signal_fit(model, result, fit_variable, mc, signal_region, mass, x_max, x_min, bin_width=50):
    bins = int((x_max-x_min)/bin_width)

    # Create a canvas and split it into two pads
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    top_pad = ROOT.TPad("top_pad", "top_pad", 0, 0.3, 1, 1)  # Top pad (main plot)
    bottom_pad = ROOT.TPad("bottom_pad", "bottom_pad", 0, 0, 1, 0.3)  # Bottom pad (pull plot)
    top_pad.Draw()
    bottom_pad.Draw()


    # Draw the main plot in the top pad
    top_pad.cd()
    top_pad.SetBottomMargin(0.02)  # Reduce margin between pads

    legend = ROOT.TLegend(0.7, 0.5, 0.9, 0.9)
    legend.SetBorderSize(0)
    legend.SetNColumns(1)
    #legend.SetTextSize(0.05)
    legend.SetFillColorAlpha(ROOT.kWhite, 0)

    frame = fit_variable.frame(x_min, x_max, bins)
    mc.plotOn(frame, ROOT.RooFit.DrawOption("PZ"))
    
    # plot errorbands
    model.plotOn(frame, LineColor='kBlue')
    chi2_ndf = frame.chiSquare(len(result.floatParsFinal()))
    legend.AddEntry(frame.getObject(0), "signal MC", "ep")
    legend.AddEntry(frame.getObject(1), f"DCB fit", "l")
    hpull = frame.pullHist(frame.getObject(0).GetName(), frame.getObject(1).GetName())
    for i in range(hpull.GetN()):
        hpull.SetPointError(i, 0, 0, 1, 1)  # Set x-error to 0 and y-error to 1

    #frame.SetMinimum(1e-2)
    frame.SetMaximum(1.2*frame.GetMaximum())
    frame.SetTitle("")
    frame.GetXaxis().SetLabelSize(0)  # Hide x-axis labels
    frame.GetXaxis().SetTickLength(0) # Hide x-axis ticks
    #frame.GetYaxis().SetTitleSize(0.1)
    #frame.GetXaxis().SetTitle('m_{j\gamma} (GeV)')
    frame.Draw()
    legend.Draw()
    
    ##########################################
    # Draw the pull plot in the bottom pad

    bottom_pad.cd()
    bottom_pad.SetTopMargin(0.04)  # Reduce margin between pads
    bottom_pad.SetBottomMargin(0.25)  # Increase bottom margin for labels

    # Create a frame for the pull plot
    pull_frame = fit_variable.frame(x_min, x_max, bins)
    # Add a horizontal line at y = 0 for reference
    zero_line = ROOT.TLine(x_min, 0, x_max, 0)
    zero_line.SetLineColor(ROOT.kBlue)
    zero_line.SetLineWidth(3)
    pull_frame.addObject(zero_line)

    # Calculate and plot the pulls for ExPow1
    pull_frame.addPlotable(hpull, "PZ")

    # plot
    pull_frame.SetTitle("")
    pull_frame.GetYaxis().SetLabelSize(0.1)
    pull_frame.GetYaxis().SetTitle("(MC - fit) / #sigma_{STAT}")
    pull_frame.GetYaxis().SetTitleOffset(0.4)
    pull_frame.GetYaxis().SetTitleSize(0.1)

    pull_frame.GetXaxis().SetTitle('m_{j#gamma} (GeV)')
    pull_frame.GetXaxis().SetTitleSize(0.1)
    pull_frame.GetXaxis().SetLabelSize(0.1)
    pull_frame.Draw()
    pull_frame.SetMaximum(+3)
    pull_frame.SetMinimum(-3)
    #bottom_legend.Draw()
    
    os.makedirs(f'../postprocess/plots/fit/{year}/{mass}', exist_ok=True)
    canvas.SaveAs(f"../postprocess/plots/fit/{year}/{mass}/signal_fit_{signal_region}.pdf")


def fit_signal(in_file, out_dir, mass, signal_region, cut, year='Run2', fit_range_low=650, fit_range_high=4000):
    sig_model_dir = os.path.join(out_dir, year, str(mass))
    os.makedirs(sig_model_dir, exist_ok=True)

    SR, width = signal_region.split('_')
    k = 0.2 if width=='N' else 0.3 if width=='W' else 0.35
    x_max = fit_range_high if (1+k)*mass > fit_range_high else round((1+k)*mass/50)*50
    x_min = fit_range_low if (1-k)*mass < fit_range_low else round((1-k)*mass/50)*50

    ## Load file and tree
    f = ROOT.TFile(in_file, "r")
    tree = f.Get("Events")

    # Define mass and weight variables
    fit_mass = ROOT.RooRealVar("fit_mass", "fit_mass", mass, x_min, x_max)
    weight = ROOT.RooRealVar("weight", "weight", 0.1, 0, 100)
    jet_mass = ROOT.RooRealVar("jet_mass", "jet_mass", 125, 0, 999)
    tagger = ROOT.RooRealVar("tagger", "tagger", 0.5, 0, 2)

    # Convert to RooDataSet
    mc = ROOT.RooDataSet("signal", "signal", tree, ROOT.RooArgSet(fit_mass, weight, jet_mass, tagger), cut, "weight")

    if width == 'N':
        sigma = np.std(uproot.open(in_file)['Events']['fit_mass'].array())
    else:
        sigma = mass * (0.056 if width=='W' else 0.10)
    # Introduce RooRealVars into the workspace for the fitted variable
    x0 = ROOT.RooRealVar("x0", "x0", mass, mass - 100, mass + 100)
    sigmaL = ROOT.RooRealVar("sigmaL", "sigmaL", sigma/2, 10, 2*sigma)
    sigmaR = ROOT.RooRealVar("sigmaR", "sigmaR", sigma/2, 10, 2*sigma)
    alphaL = ROOT.RooRealVar("alphaL", "alphaL", 0.5, 0.3, 4)
    alphaR = ROOT.RooRealVar("alphaR", "alphaR", 1.5, 0.3, 4)
    nL = ROOT.RooRealVar("nL", "nL", 2, 1e-1, 2e2)
    nR = ROOT.RooRealVar("nR", "nR", 2, 0.5, 10)

    # shape uncertainties
    JES_2016 = ROOT.RooRealVar("JES_2016", "JES_2016", 0, -5, 5)
    JES_2017 = ROOT.RooRealVar("JES_2017", "JES_2017", 0, -5, 5)
    JES_2018 = ROOT.RooRealVar("JES_2018", "JES_2018", 0, -5, 5)
    JER = ROOT.RooRealVar("JER", "JER", 0, -5, 5)
    PES = ROOT.RooRealVar("PES", "PES", 0, -5, 5)
    PER = ROOT.RooRealVar("PER", "PER", 0, -5, 5)
    JES_2016.setConstant(True); JES_2017.setConstant(True); JES_2018.setConstant(True); JER.setConstant(True); PES.setConstant(True); PER.setConstant(True)

    systematics = pd.read_csv('../src/parameters/uncertainty/systematics.csv', index_col=(0, 1))
    index = (mass, signal_region)
    mean = ROOT.RooFormulaVar("mean", "mean",
        "@0*(1+%f*(0.264*@1+0.301*@2+0.435*@3)+%f*@4)"%(systematics['JES'][index]-1, (systematics['PES'][index]-1)/2), 
        ROOT.RooArgList(x0, JES_2016, JES_2017, JES_2018, PES))
    widthL = ROOT.RooFormulaVar("widthL", "widthL", 
        "@0*(1+%f*@1+%f*@2)"%(systematics['JER'][index]-1, systematics['PER'][index]-1), 
        ROOT.RooArgList(sigmaL, JER, PER))
    widthR = ROOT.RooFormulaVar("widthR", "widthR", 
        "@0*(1+%f*@1+%f*@2)"%(systematics['JER'][index]-1, systematics['PER'][index]-1),
        ROOT.RooArgList(sigmaR, JER, PER))

    # Fit signal model to MC
    model_signal = ROOT.RooCrystalBall(f"model_bbgamma_{SR}", f"model_bbgamma_{SR}", fit_mass, mean, widthL, widthR, alphaL, nL, alphaR, nR)
    model_signal.fitTo(mc,  ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.SumW2Error(True))
    result = model_signal.fitTo(mc,  ROOT.RooFit.PrintLevel(-1), ROOT.RooFit.SumW2Error(True), Save=True)

    norm = {
        year: get_signal_norm(file_path=in_file, mass=mass, cut=cut, x_min=x_min, x_max=x_max) for year in ['2016', '2017', '2018', 'Run2']
    }
    signal_norm = {
        year: ROOT.RooRealVar(f"model_bbgamma_{SR}_norm_{year}", f"Number of signal events in {SR} {year}", norm[year], 0, 5*norm[year])
        for year in norm
    }

    x0.setConstant(True)
    sigmaL.setConstant(True)
    sigmaR.setConstant(True)
    alphaL.setConstant(True)
    alphaR.setConstant(True)
    nL.setConstant(True)
    nR.setConstant(True)
    #signal_norm.setConstant(True)


    f_out = ROOT.TFile(os.path.join(sig_model_dir, f'{signal_region}.root'), "RECREATE")
    w_sig = ROOT.RooWorkspace("workspace_signal", "workspace_signal")
    getattr(w_sig, "import")(model_signal)
    #getattr(w_sig, "import")(signal_norm)
    w_sig.Print()
    w_sig.Write()
    f_out.Close()

    with open(os.path.join(sig_model_dir, f'{signal_region}.yaml'), 'w', encoding='utf-8') as f:
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
            'sigma': float(sigma),
        }
        for year in ['2016', '2017', '2018', 'Run2']:
            info[f'norm_{year}'] = signal_norm[year].getVal()
        yaml.dump(info, f)

    bin_width_list = [5, 10, 25, 50]
    bin_width = bin_width_list[0]
    for bw in bin_width_list:
        if (x_max-x_min)/bw<=50:
            bin_width = bw
            break

    plot_signal_fit(model=model_signal, result=result, fit_variable=fit_mass, mc=mc, signal_region=signal_region, mass=mass, x_max=x_max, x_min=x_min, bin_width=bin_width)


def plot_b_only_fit(candidates, model, result, fit_variable, data, region, x_min=650, x_max=3700, bin_width=50):
    line_color = {'ExPow1': ROOT.kViolet+2, 'ExPow2': ROOT.kBlue, 'DiJet2': ROOT.kAzure+1, 'DiJet3': ROOT.kGreen+2, 'InvPow2': ROOT.kOrange-3, 'InvPow3': ROOT.kRed+1}
    #line_color = {'ExPow1': ROOT.kRed, 'ExPow2': ROOT.kGreen+1, 'DiJet2': ROOT.kYellow+2, 'DiJet3': ROOT.kCyan, 'InvPow2':ROOT.kBlue, 'InvPow3': ROOT.kMagenta}
    #band_color = {'ExPow1': ROOT.kPink, 'ExPow2': ROOT.kYellow, 'DiJet2': ROOT.kGreen, 'DiJet3': ROOT.kCyan, 'InvPow2':ROOT.kAzure, 'InvPow3': ROOT.kMagenta}

    ## plot
    bins = int((x_max-x_min)/bin_width)

    # Create a canvas and split it into two pads
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    top_pad = ROOT.TPad("top_pad", "top_pad", 0, 0.3, 1, 1)  # Top pad (main plot)
    bottom_pad = ROOT.TPad("bottom_pad", "bottom_pad", 0, 0, 1, 0.31)  # Bottom pad (pull plot)
    top_pad.Draw()
    bottom_pad.Draw()


    # Draw the main plot in the top pad
    top_pad.cd()
    top_pad.SetLogy()
    top_pad.SetBottomMargin(0.02)  # Reduce margin between pads

    legend = ROOT.TLegend(0.49, 0.5, 0.89, 0.89)
    legend.SetBorderSize(0)
    legend.SetNColumns(1)
    #legend.SetTextSize(0.03)
    #legend.SetFillColorAlpha(ROOT.kWhite, 0)

    frame = fit_variable.frame(x_min, x_max, bins)
    data.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack), ROOT.RooFit.LineColor(ROOT.kWhite), ROOT.RooFit.DrawOption("PZ"))

    # Create a histogram from the RooDataSet
    hist_data = data.createHistogram(f"hist_data", fit_variable, ROOT.RooFit.Binning(bins, x_min, x_max))
    # Convert the histogram to a RooHist object
    data_hist = ROOT.RooHist(hist_data)
    for i in range(data_hist.GetN()):
        y = data_hist.GetPointY(i)
        data_hist.SetPointError(i, 0, 0, Garwood_eror(y, 'down'), Garwood_eror(y, 'up'))
    
    # plot errorbands, len(candiates)
    for k in candidates:
        model[k].plotOn(frame, VisualizeError=(result[k], 1), FillColor='kGray', LineColor='kWhite', Name=f'error_{k}')

    chi_square = {}
    for i, k in enumerate(candidates):
        model[k].plotOn(frame, LineColor=line_color[k], Name=k)
        chi_square[i] = frame.chiSquare(len(result[k].floatParsFinal()))
    best_fit_index = min(chi_square, key=lambda i: chi_square[i])

    frame.addPlotable(data_hist, "PZ")
    legend.AddEntry(frame.getObject(2*len(candidates)+1), "Data", "ep")
    for i, k in enumerate(candidates):
        legend.AddEntry(frame.getObject(len(candidates)+1+i), f"{k}, #chi^{{2}}/NDF = {chi_square[i]:.3f}", "l")
    legend.AddEntry(frame.getObject(1), '#sigma_{SYS}', "f")

    hpull = frame.pullHist(frame.getObject(2*len(candidates)+1).GetName(), frame.getObject(len(candidates)+1+best_fit_index).GetName())
    for i in range(hpull.GetN()):
        hpull.SetPointError(i, 0, 0, 1, 1)  # Set x-error to 0 and y-error to 1

    canvas.SetLogy()
    frame.SetMinimum(7e-2)
    frame.SetMaximum(3e2)
    frame.SetTitle("")
    frame.GetXaxis().SetLabelSize(0)  # Hide x-axis labels
    frame.GetXaxis().SetTickLength(0.02)
    frame.GetXaxis().SetTitleSize(0)
    frame.GetYaxis().SetTitleSize(0.07)
    frame.GetYaxis().SetTitle("Events / 50 GeV")
    frame.GetYaxis().SetTitleOffset(0.6)
    frame.Draw()
    legend.Draw()

    cms_label = ROOT.TLatex()
    cms_label.SetNDC(True)
    cms_label.SetTextFont(61)
    cms_label.SetTextSize(0.08)
    cms_label.DrawLatex(0.15, 0.82, "CMS")

    preliminary_label = ROOT.TLatex()
    preliminary_label.SetNDC(True)
    preliminary_label.SetTextFont(52)
    preliminary_label.SetTextSize(0.06)
    preliminary_label.DrawLatex(0.28, 0.82, "Preliminary")

    region_label = ROOT.TLatex()
    region_label.SetNDC(True)
    region_label.SetTextFont(42)
    region_label.SetTextSize(0.06)
    region_label.DrawLatex(0.2, 0.75, region)

    lumi_label = ROOT.TLatex()
    lumi_label.SetNDC(True)
    lumi_label.SetTextFont(42)
    lumi_label.SetTextSize(0.06)
    lumi_label.SetTextAlign(31)
    lumi_label.DrawLatex(0.9, 0.92, "138 fb^{-1} (13 TeV)")

    ##########################################
    # Draw the pull plot in the bottom pad

    bottom_pad.cd()
    bottom_pad.SetTopMargin(0.04)  # Reduce margin between pads
    bottom_pad.SetBottomMargin(0.25)  # Increase bottom margin for labels

    bottom_legend = ROOT.TLegend(0.55, 0.68, 0.89, 0.99)
    bottom_legend.SetBorderSize(0)
    bottom_legend.SetFillColorAlpha(ROOT.kWhite, 0)
    bottom_legend.SetNColumns(1)
    #bottom_legend.SetTextSize(0.08)

    # Create a frame for the pull plot
    pull_frame = fit_variable.frame(x_min, x_max, bins)

    # sigma_sys/sigma_stats
    x_centers = np.array(frame.getObject(2*len(candidates)+1).GetX())
    error_stats = np.array(frame.getObject(2*len(candidates)+1).GetEYhigh()) + np.array(frame.getObject(2*len(candidates)+1).GetEYlow())

    cs = {}
    for k in candidates:
        x=np.array(frame.findObject(f'error_{k}').GetX())
        y=np.array(frame.findObject(f'error_{k}').GetY())
        N=int(len(x)/2)
        error_bar = {}
        for i in x[:N]:
            if not x_min<=i<=x_max:
                continue
            y_up = np.max(y[x==i])
            _y_down = np.min(y[x==i])
            y_down = np.where(_y_down>=0, _y_down, 0)
            error_bar[i] = y_up - y_down
        cs[k]=CubicSpline(x=list(error_bar.keys()), y=list(error_bar.values()))

    error_sys = [cs[k](x_centers) for k in candidates]
    error_sys_over_stats = np.max(error_sys, axis=0) / error_stats

    roohist = ROOT.RooHist()
    for i in range(len(x_centers)):
        x = x_centers[i]
        roohist.addBinWithXYError(x, 0, bin_width/2, bin_width/2, error_sys_over_stats[i], error_sys_over_stats[i])
    roohist.SetFillColor(ROOT.kGray)
    #roohist.SetFillStyle(3001)
    #roohist.SetMarkerSize(0)
    pull_frame.addPlotable(roohist, "E2")

    # Add a horizontal line at y = 0 for reference
    zero_line = ROOT.TLine(x_min, 0, x_max, 0)
    zero_line.SetLineColor(ROOT.kBlack)
    zero_line.SetLineWidth(2)
    pull_frame.addObject(zero_line)
    #bottom_legend.AddEntry(pull_frame.getObject(1), candidates[0], "l")

    # Calculate and plot the pulls for best-fit function
    pull_frame.addPlotable(hpull, "PZ")
    bottom_legend.AddEntry(pull_frame.getObject(2), '(Data - '+candidates[best_fit_index]+') / #sigma_{STAT}', "ep")
    bottom_legend.AddEntry(pull_frame.getObject(0), '#sigma_{SYS}/#sigma_{STAT}', "f")

    # plot
    pull_frame.SetTitle("")
    pull_frame.GetYaxis().SetLabelSize(0.1)
    pull_frame.GetYaxis().SetTitle("Pull")
    pull_frame.GetYaxis().SetTitleOffset(0.3)
    pull_frame.GetYaxis().SetTitleSize(0.15)

    pull_frame.GetXaxis().SetTitle("m_{j#gamma} (GeV)")
    pull_frame.GetXaxis().SetTitleSize(0.15)
    pull_frame.GetXaxis().SetLabelSize(0.08)
    pull_frame.GetXaxis().SetTitleOffset(0.68)
    pull_frame.Draw()
    pull_frame.SetMaximum(+3.5)
    pull_frame.SetMinimum(-3.5)
    bottom_legend.Draw()
    
    os.makedirs('../postprocess/plots/fit/Run2', exist_ok=True)
    plot_name = candidates[0] if len(candidates)==1 else len(candidates)
    canvas.SaveAs(f"../postprocess/plots/fit/Run2/b_only_fit_{region}_{plot_name}.pdf")
    canvas.SaveAs(f"../postprocess/plots/fit/Run2/b_only_fit_{region}_{plot_name}.C")


def fit_background(in_file, out_dir, region, cut, year='Run2', fit_range_low=650, fit_range_high=4000):
    bkg_model_dir = os.path.join(out_dir, year)
    os.makedirs(bkg_model_dir, exist_ok=True)
    
    # # Background modelling
    f = ROOT.TFile(in_file, "r")
    # Load TTree
    tree = f.Get("Events")

    # Define mass and weight variables
    fit_mass = ROOT.RooRealVar("fit_mass", "fit_mass", 1500, fit_range_low, fit_range_high)
    weight = ROOT.RooRealVar("weight", "weight", 1, -10, 1e3)
    jet_mass = ROOT.RooRealVar("jet_mass", "jet_mass", 125, 0, 999)
    tagger = ROOT.RooRealVar("tagger", "tagger", 0.5, 0, 2)

    # Convert to RooDataSet
    data_region = ROOT.RooDataSet(f"data_{region}", f"data_{region}", tree, ROOT.RooArgSet(fit_mass, weight, jet_mass, tagger), cut, "weight")

    ## Multiple background models
    model, p1, p2, p3, result = {}, {}, {}, {}, {}
    energy = 1e2

    # ExPow1 model
    p1['ExPow1'] = ROOT.RooRealVar("p1_ExPow1", "p1_ExPow1", -1, -10, 0)
    model['ExPow1'] = ROOT.RooGenericPdf("model_background_ExPow1", "model_background_ExPow1", f"TMath::Power(@0/{energy}, @1)", ROOT.RooArgList(fit_mass, p1['ExPow1']))

    # ExPow2 model
    p1['ExPow2'] = ROOT.RooRealVar("p1_ExPow2", "p1_ExPow2", -1, -10, 0)
    p2['ExPow2'] = ROOT.RooRealVar("p2_ExPow2", "p2_ExPow2", -1e-2, -0.5, 0.5)
    model['ExPow2'] = ROOT.RooGenericPdf("model_background_ExPow2", "model_background_ExPow2", f"TMath::Power(@0/{energy}, @1) * TMath::Exp(@2 * @0/{energy})", ROOT.RooArgList(fit_mass, p1['ExPow2'], p2['ExPow2']))

    # DiJet2 model
    p1['DiJet2'] = ROOT.RooRealVar("p1_DiJet2", "p1_DiJet2", -2, -10, 0)
    p2['DiJet2'] = ROOT.RooRealVar("p2_DiJet2", "p2_DiJet2", -1, -5, 0)
    model['DiJet2'] = ROOT.RooGenericPdf("model_background_DiJet2", "model_background_DiJet2", f"TMath::Power(@0/{energy}, @1 + @2 * TMath::Log(@0/{energy}))", ROOT.RooArgList(fit_mass, p1['DiJet2'], p2['DiJet2']))

    # DiJet3 model
    p1['DiJet3'] = ROOT.RooRealVar("p1_DiJet3", "p1_DiJet3", -1, -10, 0)
    p2['DiJet3'] = ROOT.RooRealVar("p2_DiJet3", "p2_DiJet3", -1, -5, 0)
    p3['DiJet3'] = ROOT.RooRealVar("p3_DiJet3", "p3_DiJet3", -1e-3, -0.1, 0.1)
    model['DiJet3'] = ROOT.RooGenericPdf("model_background_DiJet3", "model_background_DiJet3", f"TMath::Power(@0/{energy}, @1 + @2 * TMath::Log(@0/{energy}) + @3 * TMath::Power(TMath::Log(@0/{energy}), 2))", ROOT.RooArgList(fit_mass, p1['DiJet3'], p2['DiJet3'], p3['DiJet3']))

    # InvPow2 model
    p1['InvPow2'] = ROOT.RooRealVar("p1_InvPow2", "p1_InvPow2", 1e-2, 0, 10)
    p2['InvPow2'] = ROOT.RooRealVar("p2_InvPow2", "p2_InvPow2", -2, -30, 0)
    model['InvPow2'] = ROOT.RooGenericPdf("model_background_InvPow2", "model_background_InvPow2", f"TMath::Power(1 + @1*@0/{energy}, @2)", ROOT.RooArgList(fit_mass, p1['InvPow2'], p2['InvPow2']))

    # InvPow3 model
    p1['InvPow3'] = ROOT.RooRealVar("p1_InvPow3", "p1_InvPow3", 1e-2, 0, 10)
    p2['InvPow3'] = ROOT.RooRealVar("p2_InvPow3", "p2_InvPow3", -2, -30, 0)
    p3['InvPow3'] = ROOT.RooRealVar("p3_InvPow3", "p3_InvPow3", 0.5, -0.05 if region in ['SRH2', 'CR2'] else -1, 1)
    model['InvPow3'] = ROOT.RooGenericPdf("model_background_InvPow3", "model_background_InvPow3", f"TMath::Power(1 + @1*@0/{energy}, @2 + @3*@0/{energy})", ROOT.RooArgList(fit_mass, p1['InvPow3'], p2['InvPow3'], p3['InvPow3']))

    # Make a RooCategory object: this will control which PDF is "active"
    category = ROOT.RooCategory(f"pdfindex_{region}", "Index of Pdf which is active")
    # Make a RooArgList of the models
    models = ROOT.RooArgList()

    # Fit models
    for k in ['ExPow1', 'ExPow2', 'DiJet2', 'DiJet3', 'InvPow2', 'InvPow3']:
        model[k].fitTo(data_region, ROOT.RooFit.SumW2Error(True))
        result[k] = model[k].fitTo(data_region, ROOT.RooFit.SumW2Error(True), Save=True)
        p1[k].setConstant(True)
        if k in p2:
            p2[k].setConstant(True)
        if k in p3:
            p3[k].setConstant(True)
        models.add(model[k])
        plot_b_only_fit(candidates=[k], model=model, result=result, fit_variable=fit_mass, data=data_region, region=region, x_min=fit_range_low, x_max=fit_range_high-300, bin_width=50)

    plot_b_only_fit(candidates=['ExPow1', 'ExPow2', 'DiJet2', 'DiJet3', 'InvPow2', 'InvPow3'], model=model, result=result, fit_variable=fit_mass, data=data_region, region=region, x_min=fit_range_low, x_max=fit_range_high-300, bin_width=50)

    # Build the RooMultiPdf object
    multipdf = ROOT.RooMultiPdf(f"multipdf_{region}", f"multipdf_{region}", category, models)
    background_norm = ROOT.RooRealVar(f"multipdf_{region}_norm", "Number of background events", data_region.sumEntries(), 0, 100 * data_region.sumEntries())
    background_norm.setConstant(False)

    # store workspace
    with open(f'{bkg_model_dir}/background_{region}.yaml', 'w', encoding='utf-8') as f:
        info = {
            'p1': {k: p1[k].getVal() for k in p1},
            'p2': {k: p2[k].getVal() for k in p2},
            'p3': {k: p3[k].getVal() for k in p3},
            'event_sum': data_region.sumEntries(),
            'norm': background_norm.getVal()
        }
        yaml.dump(info, f)

    f_out = ROOT.TFile(os.path.join(bkg_model_dir, f'data_{region}.root'), "RECREATE")
    w_bkg = ROOT.RooWorkspace(f"workspace_{region}", f"workspace_{region}")
    getattr(w_bkg, "import")(data_region)
    getattr(w_bkg, "import")(multipdf)
    for k in model:
        getattr(w_bkg, "import")(model[k])
    getattr(w_bkg, "import")(background_norm)
    w_bkg.Print()
    w_bkg.Write()
    f_out.Close()


if __name__ == "__main__":
    tagger_cut = {
        'SR1': [0.8, 0.98],
        'SR2': [0.98, 2],
    }
    mass_cut = {
        'Z': [80, 110],
        'H': [110, 150],
    }

    args = parse_commandline()
    fit_range_low, fit_range_high = args.fit_range_low, args.fit_range_high
    year = args.year
    if args.year is None:
        YEARS = ['2016', '2017', '2018', 'Run2']
    else:
        YEARS = [args.year]
    if args.signal_mass is not None:
        signal_mass = [args.signal_mass]
    else:
        signal_mass = [700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 3000, 3500]
    if args.signal_region is not None:
        signal_regions = [args.signal_region]
    else:
        signal_regions = ['SRH1_N', 'SRH2_N', 'SRZ1_N', 'SRZ2_N', 'SRZ1_W', 'SRZ2_W', 'SRZ1_VW', 'SRZ2_VW']

    for year in YEARS:
        for signal_region in signal_regions:
            jet = signal_region[2]
            mass_low, mass_high = mass_cut[jet]
            tagger_region = signal_region[:2]+signal_region[3]
            tagger_cut_low, tagger_cut_high = tagger_cut[tagger_region]
            CR = 'CR1' if '1' in signal_region else 'CR2'
            CR_cut = f"""(
                (((jet_mass>50) & (jet_mass<{mass_cut['Z'][0]})) | (jet_mass>{mass_cut['H'][1]})) & 
                (tagger>{tagger_cut_low}) & (tagger<{tagger_cut_high})
            )"""
            SR_cut = f"""(
                (jet_mass>{mass_low}) & (jet_mass<{mass_high}) & 
                (tagger>{tagger_cut_low}) & (tagger<{tagger_cut_high})
            )"""
            if year=='Run2':
                fit_background(in_file=os.path.join(args.in_dir, year, 'data.root'), out_dir=args.out_dir, region=signal_region.split('_')[0], cut=SR_cut)
                fit_background(in_file=os.path.join(args.in_dir, year, 'data.root'), out_dir=args.out_dir, region=CR, cut=CR_cut)

            for mass in signal_mass:
                fit_signal(
                    in_file=os.path.join(args.in_dir, year, str(mass), f'{signal_region[:3]+signal_region[4:]}.root'), 
                    out_dir=args.out_dir, mass=mass, signal_region=signal_region, cut=SR_cut
                )
