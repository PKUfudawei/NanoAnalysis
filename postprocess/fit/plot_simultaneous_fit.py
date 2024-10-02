#!/usr/bin/env python3
import ROOT
ROOT.gStyle.SetOptTitle(0)

y = 'Run2'
mass = 1000
sr = 'SR2'

x_low = 720
if mass == 1000:
    x_high = 2020
elif mass == 2000:
    x_high = 3020
elif mass == 3000:
    x_high = 3520
else:
    x_high = None


def main() -> None:
    f = ROOT.TFile(f"higgsCombine{mass}_{sr}.MultiDimFit.mH120.root")
    w = f.Get("w")
    w.Print("v")

    n_bins = (x_high - x_low) // 50
    binning = ROOT.RooFit.Binning(n_bins, x_low, x_high)

    can = ROOT.TCanvas()
    plot = w.var("mass_Zprime").frame()
    w.data("data_obs").plotOn(plot, binning, ROOT.RooFit.Name("bkg MC"))

    # Load the S+B model
    sb_model = w.pdf("model_s").getPdf("Tag0")

    # Prefit
    sb_model.plotOn(plot, ROOT.RooFit.LineColor(4), ROOT.RooFit.Name("prefit"))

    # Postfit
    w.loadSnapshot("MultiDimFit")
    sb_model.plotOn(plot, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name("postfit"))
    r_bestfit = w.var("r").getVal()

    plot.Draw()

    leg = ROOT.TLegend(0.35, 0.6, 0.85, 0.85)
    leg.AddEntry("bkg MC", "pseudo data", "PE")
    leg.AddEntry("prefit", "Prefit S+B model (r=1.00)", "L")
    leg.AddEntry("postfit", "Postfit S+B model (r=%.2g)" % r_bestfit, "L")
    leg.SetBorderSize(0)
    leg.SetTextSize(0.05)
    leg.Draw("Same")

    can.Update()
    can.SaveAs(f"../plots/fit/{y}/fit_test_bkg_mc_{mass}_{sr}.pdf")


if __name__ == "__main__":
    main()
