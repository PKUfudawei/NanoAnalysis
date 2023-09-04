#!/usr/bin/env python3
import ROOT

f = ROOT.TFile("higgsCombine.bestfit.MultiDimFit.mH125.root")
w = f.Get("w")
w.Print("v")

x_low, x_high = 720, 2000
n_bins = (x_high - x_low) // 50
binning = ROOT.RooFit.Binning(n_bins, x_low, x_high)

can = ROOT.TCanvas()
plot = w.var("mass_Zprime").frame()
w.data("data_obs").plotOn(plot, binning)

# Load the S+B model
sb_model = w.pdf("model_s").getPdf("Tag0")

# Prefit
sb_model.plotOn(plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("prefit"))

# Postfit
w.loadSnapshot("MultiDimFit")
sb_model.plotOn(plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("postfit"))
r_bestfit = w.var("r").getVal()

plot.Draw()

leg = ROOT.TLegend(0.55, 0.6, 0.85, 0.85)
leg.AddEntry("prefit", "Prefit S+B model (r=1.00)", "L")
leg.AddEntry("postfit", "Postfit S+B model (r=%.2g)" % r_bestfit, "L")
leg.Draw("Same")

y='2018'
can.Update()
can.SaveAs(f"../plots/fit/{y}/model_mc_SR.pdf")
