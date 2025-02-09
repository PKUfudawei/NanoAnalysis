import ROOT, os
import numpy as np

def Garwood_eror(number, direction):
    upper = np.array([1.84, 3.30, 4.64, 5.92, 7.16, 8.38, 9.58, 10.77, 11.95, 13.11, 14.27])
    lower = np.array([0, 0.17, 0.71, 1.37, 2.09, 2.84, 3.62, 4.42, 5.23, 6.06, 6.89])
    center = np.arange(11)
    if number > 10:
        return np.sqrt(number)
    elif direction=='up':
        return (upper-center)[int(number)]
    elif direction=='down':
        return (center-lower)[int(number)]
    else:
        return None

line_color = {'expow1': ROOT.kRed, 'expow2': ROOT.kGreen+1, 'dijet2': ROOT.kYellow+2, 'dijet3': ROOT.kCyan, 'invpow2':ROOT.kBlue, 'invpow3': ROOT.kMagenta}
band_color = {'expow1': ROOT.kPink, 'expow2': ROOT.kYellow, 'dijet2': ROOT.kGreen, 'dijet3': ROOT.kCyan, 'invpow2':ROOT.kAzure, 'invpow3': ROOT.kMagenta}


def fit_error_band(candidates, model, result, fit_variable, data_region, SR, jet, x_min=650, x_max=3700, bin_width=50, line_color=line_color):
    ## plot
    bins = int((x_max-x_min)/bin_width)
    plot_name = (candidates[0] if len(candidates) == 1 else f"{len(candidates)}")

    # Create a canvas and split it into two pads
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    top_pad = ROOT.TPad("top_pad", "top_pad", 0, 0.3, 1, 1)  # Top pad (main plot)
    bottom_pad = ROOT.TPad("bottom_pad", "bottom_pad", 0, 0, 1, 0.3)  # Bottom pad (pull plot)
    top_pad.Draw()
    bottom_pad.Draw()

    # Draw the main plot in the top pad
    top_pad.cd()
    top_pad.SetLogy()
    top_pad.SetBottomMargin(0.02)  # Reduce margin between pads

    frame = fit_variable.frame(Bins=bins)

    # Create a histogram from the RooDataSet
    hist_data = data_region.createHistogram("hist_data", fit_variable, ROOT.RooFit.Binning(bins, x_min, x_max))

    # Convert the histogram to a RooHist object
    data_hist = ROOT.RooHist(hist_data)

    for i in range(data_hist.GetN()):
        y = data_hist.GetPointY(i)
        data_hist.SetPointError(i, 0, 0, Garwood_eror(y, 'down'), Garwood_eror(y, 'up'))

    # Plot the data with custom errors
    frame.addPlotable(data_hist, "P")
    model[candidates[0]].plotOn(frame)
    hpull = frame.pullHist()
    # Set the error bars of the pulls to 1
    for i in range(hpull.GetN()):
        hpull.SetPointError(i, 0, 0, 1, 1)  # Set x-error to 0 and y-error to 1
    data_region.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kWhite), ROOT.RooFit.LineColor(ROOT.kWhite))

    legend = ROOT.TLegend(0.5, 0.55, 0.79, 0.89)
    #legend.SetTextSize(0.05)
    legend.SetBorderSize(0)
    legend.AddEntry(frame.getObject(0), "data", "lep")

    for k in candidates:
        model[k].plotOn(frame, VisualizeError=(result[k], 1), FillColor='kGray')
    #legend.AddEntry(frame.getObject(1), 'Stats. Unc.', "f")


    for i, k in enumerate(candidates):
        model[k].plotOn(frame, LineColor=line_color[k])
        chi2_ndf = frame.chiSquare(len(result[k].floatParsFinal()))
        legend.AddEntry(frame.getObject(i + len(candidates) + 3), f"{k}, #chi^{{2}}/NDF = {chi2_ndf:.2f}", "l")

    frame.addPlotable(data_hist, "P")

    canvas.SetLogy()
    frame.SetMinimum(1e-2)
    frame.SetTitle("")
    frame.GetXaxis().SetLabelSize(0)  # Hide x-axis labels
    frame.GetXaxis().SetTickLength(0) # Hide x-axis ticks
    #frame.GetXaxis().SetTitle('m_{j\gamma} [GeV]')
    frame.Draw()
    legend.Draw()

    # Draw the pull plot in the bottom pad
    bottom_pad.cd()
    bottom_pad.SetTopMargin(0.04)  # Reduce margin between pads
    bottom_pad.SetBottomMargin(0.25)  # Increase bottom margin for labels

    # Create a frame for the pull plot
    pull_frame = fit_variable.frame(x_min, x_max, bins)
    pull_frame.SetTitle("")
    pull_frame.GetYaxis().SetLabelSize(0.1)
    pull_frame.GetYaxis().SetTitle(f"Pull w.r.t. {candidates[0]}")
    pull_frame.GetYaxis().SetTitleOffset(0.4)
    pull_frame.GetYaxis().SetTitleSize(0.1)

    pull_frame.GetXaxis().SetTitle('m_{j#gamma} [GeV]')
    pull_frame.GetXaxis().SetTitleSize(0.1)
    pull_frame.GetXaxis().SetLabelSize(0.1)

    # Calculate and plot the pulls for expow1
    pull_frame.addPlotable(hpull, "P")
    pull_frame.Draw()
    pull_frame.SetMaximum(+3)
    pull_frame.SetMinimum(-3)


    # Add a horizontal line at y = 0 for reference
    zero_line = ROOT.TLine(x_min, 0, x_max, 0)
    zero_line.SetLineColor(line_color[candidates[0]])
    #zero_line.SetLineStyle(2)
    zero_line.Draw()

    os.makedirs('../plots/fit/Run2', exist_ok=True)
    canvas.SaveAs(f"../plots/fit/Run2/{SR}_{jet}_{plot_name}_errorband.pdf")

def main(SR, jet, x_max=3700, x_min=650, bin_width=50):
    tagger_cut = {
        'down': {'SR1': 0.8, 'SR2': 0.98},
        'up': {'SR1': 0.98, 'SR2': 2}
    }
    mass_cut = {
        'Z': [80, 110],
        'H': [110, 150]
    }

    CR_cut = f"""(
        (((jet_mass>50) & (jet_mass<80)) | (jet_mass>150)) & 
        (tagger>{tagger_cut['down'][SR]}) & (tagger<{tagger_cut['up'][SR]})
    )"""
    SR_cut = f"""(
        (jet_mass>{mass_cut[jet][0]}) & (jet_mass<{mass_cut[jet][1]}) & 
        (tagger>{tagger_cut['down'][SR]}) & (tagger<{tagger_cut['up'][SR]})
    )"""


    # Read files
    f = ROOT.TFile(f"input/Run2/data.root", "r")
    tree = f.Get("Events")
    fit_mass = ROOT.RooRealVar("fit_mass", "fit_mass", 1500, x_min, x_max)
    weight = ROOT.RooRealVar("weight", "weight", 1, -10, 100)
    jet_mass = ROOT.RooRealVar("jet_mass", "jet_mass", 125, 0, 999)
    tagger = ROOT.RooRealVar("tagger", "tagger", 0.5, 0, 2)
    data_region = ROOT.RooDataSet("data_region", "data_region", tree, ROOT.RooArgSet(fit_mass, weight, jet_mass, tagger), SR_cut, "weight")

    ## Multiple background models
    model, p1, p2, p3, result = {}, {}, {}, {}, {}
    energy = 1e2
    # dijet2 model
    p1['dijet2'] = ROOT.RooRealVar("p1_dijet2", "p1_dijet2", -2, -10, 0)
    p2['dijet2'] = ROOT.RooRealVar("p2_dijet2", "p2_dijet2", -1, -5, 0)
    model['dijet2'] = ROOT.RooGenericPdf("model_background_dijet2", "model_background_dijet2", f"TMath::Power(@0/{energy}, @1 + @2 * TMath::Log(@0/{energy}))", ROOT.RooArgList(fit_mass, p1['dijet2'], p2['dijet2']))

    # dijet3 model
    p1['dijet3'] = ROOT.RooRealVar("p1_dijet3", "p1_dijet3", -1, -10, 0)
    p2['dijet3'] = ROOT.RooRealVar("p2_dijet3", "p2_dijet3", -1, -5, 0)
    p3['dijet3'] = ROOT.RooRealVar("p3_dijet3", "p3_dijet3", -1e-3, -0.1, 0.1)
    model['dijet3'] = ROOT.RooGenericPdf("model_background_dijet3", "model_background_dijet3", f"TMath::Power(@0/{energy}, @1 + @2 * TMath::Log(@0/{energy}) + @3 * TMath::Power(TMath::Log(@0/{energy}), 2))", ROOT.RooArgList(fit_mass, p1['dijet3'], p2['dijet3'], p3['dijet3']))

    # expow1 model
    p1['expow1'] = ROOT.RooRealVar("p1_expow1", "p1_expow1", -1, -10, 0)
    model['expow1'] = ROOT.RooGenericPdf("model_background_expow1", "model_background_expow1", f"TMath::Power(@0/{energy}, @1)", ROOT.RooArgList(fit_mass, p1['expow1']))

    # expow2 model
    p1['expow2'] = ROOT.RooRealVar("p1_expow2", "p1_expow2", -1, -10, 0)
    p2['expow2'] = ROOT.RooRealVar("p2_expow2", "p2_expow2", -1e-2, -0.5, 0.5)
    model['expow2'] = ROOT.RooGenericPdf("model_background_expow2", "model_background_expow2", f"TMath::Power(@0/{energy}, @1) * TMath::Exp(@2 * @0/{energy})", ROOT.RooArgList(fit_mass, p1['expow2'], p2['expow2']))

    # invpow2 model
    p1['invpow2'] = ROOT.RooRealVar("p1_invpow2", "p1_invpow2", 1e-2, 0, 10)
    p2['invpow2'] = ROOT.RooRealVar("p2_invpow2", "p2_invpow2", -2, -30, 0)
    model['invpow2'] = ROOT.RooGenericPdf("model_background_invpow2", "model_background_invpow2", f"TMath::Power(1 + @1*@0/{energy}, @2)", ROOT.RooArgList(fit_mass, p1['invpow2'], p2['invpow2']))

    # invpow3 model
    p1['invpow3'] = ROOT.RooRealVar("p1_invpow3", "p1_invpow3", 1e-2, 0, 10)
    p2['invpow3'] = ROOT.RooRealVar("p2_invpow3", "p2_invpow3", -2, -30, 0)
    p3['invpow3'] = ROOT.RooRealVar("p3_invpow3", "p3_invpow3", -0.5, -1, 1)
    model['invpow3'] = ROOT.RooGenericPdf("model_background_invpow3", "model_background_invpow3", f"TMath::Power(1 + @1*@0/{energy}, @2 + @3*@0/{energy})", ROOT.RooArgList(fit_mass, p1['invpow3'], p2['invpow3'], p3['invpow3']))

    for k in model:
        result[k] = model[k].fitTo(data_region, ROOT.RooFit.SumW2Error(True), Save=True)
        p1[k].setConstant(True)
        if k in p2:
            p2[k].setConstant(True)
        if k in p3:
            p3[k].setConstant(True)

    for k in ['expow1', 'expow2', 'dijet2', 'dijet3', 'invpow2', 'invpow3']:
        fit_error_band(candidates=[k], model=model, result=result, fit_variable=fit_mass, data_region=data_region, SR=SR, jet=jet, x_min=x_min, x_max=x_max, bin_width=bin_width)

    fit_error_band(candidates=['expow1', 'expow2', 'dijet2', 'dijet3', 'invpow2', 'invpow3'], model=model, result=result, fit_variable=fit_mass, data_region=data_region, SR=SR, jet=jet, x_min=x_min, x_max=x_max, bin_width=bin_width)

main('SR1', 'Z')
"""
for SR in ['SR1', 'SR2']:
    for jet in ['H', 'Z']:
        main(SR, jet)
"""
