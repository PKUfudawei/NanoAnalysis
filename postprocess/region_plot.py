#!/usr/bin/env python3

import ROOT, os

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
ROOT.gStyle.SetLineWidth(2)

c = ROOT.TCanvas("c", "", 800, 470)
c.SetFrameFillStyle(0)
c.SetFrameBorderMode(0)
c.SetFrameLineColor(0)
c.SetFrameLineWidth(0)

# Adjust margins
c.SetLeftMargin(0.08)
c.SetRightMargin(0.01)  # Reduced right margin
c.SetTopMargin(0.08)
c.SetBottomMargin(0.18)

# Axis range adjusted to remove space on the right
h = ROOT.TH1F("h", "", 1, 40, 172.5)
h.SetMinimum(0.76)
h.SetMaximum(1.018)
h.GetXaxis().SetTitleSize(0.058)
h.GetXaxis().SetLabelSize(0.058)
h.GetYaxis().SetTitleSize(0.058)
h.GetYaxis().SetLabelSize(0.058)
h.GetYaxis().SetNdivisions(503)
h.GetXaxis().SetTickLength(0.02)
h.GetYaxis().SetTickLength(0.02)

h.Draw("AXIS")

# Axes arrows
a_xmin, a_xmax = 40, 172.5
a_ymin, a_ymax = 0.76, 1.018

x_arrow = ROOT.TArrow(a_xmin, a_ymin, a_xmax, a_ymin, 0.02, ">")
x_arrow.SetLineWidth(2)
x_arrow.Draw()

y_arrow = ROOT.TArrow(a_xmin, a_ymin, a_xmin, a_ymax, 0.02, ">")
y_arrow.SetLineWidth(2)
y_arrow.Draw()

# Axis labels
x_label = ROOT.TLatex()
x_label.SetTextFont(42)
x_label.SetTextSize(0.075)
x_label.SetTextAlign(21)
x_label.DrawLatex(160, 0.7, "m_{j} (GeV)")

y_label = ROOT.TLatex()
y_label.SetTextFont(42)
y_label.SetTextSize(0.075)
y_label.SetTextAlign(30)
y_label.SetTextAngle(90)
y_label.DrawLatex(31, 1.01, "Xbb score")

# Custom tick label at y = 0.98
y_tick = ROOT.TLatex();
y_tick.SetTextFont(42);
y_tick.SetTextSize(0.053);
y_tick.SetTextAlign(32);
y_tick.DrawLatex(39.5, 0.98, "0.98");

# Horizontal dashed lines
lines = []
for yval in [1.0, 0.98, 0.8]:
    l = ROOT.TLine(50, yval, 172.5, yval)
    l.SetLineStyle(2)
    l.SetLineWidth(2)
    l.Draw()
    lines.append(l)

# Vertical dashed lines at mj = 50, 80, 110, 150
for xval in [50, 80, 110, 150]:
    l = ROOT.TLine(xval, 0.76, xval, 1)
    l.SetLineStyle(2)
    l.SetLineWidth(2)
    l.Draw()
    lines.append(l)

c.SetLeftMargin(0.09);

# Boxes for different regions
cr2_box_A = ROOT.TBox(50, 0.98, 80, 1.0);     cr2_box_A.SetFillColorAlpha(17, 0.5); cr2_box_A.SetLineWidth(0); cr2_box_A.Draw("same")
cr2_box_B = ROOT.TBox(150, 0.98, 172.5, 1.0); cr2_box_B.SetFillColorAlpha(17, 0.5); cr2_box_B.SetLineWidth(0); cr2_box_B.Draw("same")
SRZ2_box  = ROOT.TBox(80, 0.98, 110, 1.0);     SRZ2_box.SetFillColorAlpha(64, 0.5); SRZ2_box.SetLineWidth(0);  SRZ2_box.Draw("same")
SRH2_box  = ROOT.TBox(110, 0.98, 150, 1.0);    SRH2_box.SetFillColorAlpha(50, 0.5); SRH2_box.SetLineWidth(0);  SRH2_box.Draw("same")
cr1_box_A = ROOT.TBox(50, 0.8, 80, 0.98);      cr1_box_A.SetFillColorAlpha(18, 0.5); cr1_box_A.SetLineWidth(0); cr1_box_A.Draw("same")
cr1_box_B = ROOT.TBox(150, 0.8, 172.5, 0.98);  cr1_box_B.SetFillColorAlpha(18, 0.5); cr1_box_B.SetLineWidth(0); cr1_box_B.Draw("same")
SRZ1_box  = ROOT.TBox(80, 0.8, 110, 0.98);     SRZ1_box.SetFillColorAlpha(91, 0.5); SRZ1_box.SetLineWidth(0);  SRZ1_box.Draw("same")
SRH1_box  = ROOT.TBox(110, 0.8, 150, 0.98);    SRH1_box.SetFillColorAlpha(8, 0.5);  SRH1_box.SetLineWidth(0);  SRH1_box.Draw("same")

# Region labels
text = ROOT.TLatex()
text.SetTextFont(42)
text.SetTextSize(0.077)
text.SetTextAlign(22)
text.DrawLatex(65, 0.99, "CR2");  text.DrawLatex(162, 0.99, "CR2")
text.DrawLatex(95, 0.99, "SRZ2"); text.DrawLatex(130, 0.99, "SRH2")
text.DrawLatex(65, 0.90, "CR1");  text.DrawLatex(162, 0.90, "CR1")
text.DrawLatex(95, 0.90, "SRZ1"); text.DrawLatex(130, 0.90, "SRH1")

# Save outputs
c.Draw()
c.Update()
#c.Print("Figure_004-a.pdf")
os.makedirs("./plots", exist_ok=True)
c.SaveAs("./plots/Figure_004-a.pdf")
#os.system("display ./Figure_004-a.png &")
