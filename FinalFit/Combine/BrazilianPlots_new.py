import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, TSpline3
import CMS_lumi, tdrstyle
import subprocess  # to execute shell command
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# CMS style
CMS_lumi.cmsText = "CMS"
CMS_lumi.extraText = "Preliminary"
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = True
tdrstyle.setTDRStyle()

# GET limits from root file
def getLimits(file_name):
    file = TFile(file_name)
    tree = file.Get("limit")

    limits = []
    for quantile in tree:
        limits.append(tree.limit)
        print ">>>   %.2f" % limits[-1]

    return limits[:6]

# PLOT upper limits
def plotUpperLimits(labels):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/
    N = len(labels)
    yellow = TGraph(2 * N)    # yellow band
    green = TGraph(2 * N)     # green band
    median = TGraph(N)        # median line
    observed = TGraph(N)      # observed line

    up2s = []
    for i in range(N):
        file_name = "higgsCombine.HToAATo2B2G_m"+str(labels[i])+"_zh_final_v2.AsymptoticLimits.mH"+str(labels[i])+".root"
        limit = getLimits(file_name)
        up2s.append(limit[4])
        yellow.SetPoint(i, labels[i], limit[4])  # + 2 sigma
        green.SetPoint(i, labels[i], limit[3])   # + 1 sigma
        median.SetPoint(i, labels[i], limit[2])  # median
        observed.SetPoint(i, labels[i], limit[5]) # observed
        green.SetPoint(2 * N - 1 - i, labels[i], limit[1])  # - 1 sigma
        yellow.SetPoint(2 * N - 1 - i, labels[i], limit[0])  # - 2 sigma

    W = 800
    H = 600
    T = 0.08 * H
    B = 0.12 * H
    L = 0.12 * W
    R = 0.04 * W
    c = TCanvas("c", "c", 100, 100, W, H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin(L / W)
    c.SetRightMargin(R / W)
    c.SetTopMargin(T / H)
    c.SetBottomMargin(B / H)
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid()
    c.cd()
    frame = c.DrawFrame(min(labels)-1, 0.001, max(labels)+1, max(up2s) * 2)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("95% upper limit on #sigma / #sigma_{SM}")
    frame.GetXaxis().SetTitle("M_{A} (GeV)")
    frame.SetMinimum(0)
#    frame.SetMaximum(max(up2s) * 2)
    frame.SetMaximum(10)
    frame.GetXaxis().SetLimits(min(labels), max(labels))

    # Create TSpline3 for smooth uncertainty bands
    spline_yellow_up = TSpline3("spline_yellow_up", yellow)
    spline_yellow_down = TSpline3("spline_yellow_down", yellow)
    spline_green_up = TSpline3("spline_green_up", green)
    spline_green_down = TSpline3("spline_green_down", green)

    spline_yellow_up.SetLineColor(ROOT.kOrange)
    spline_yellow_up.SetFillColor(ROOT.kOrange)
    spline_yellow_down.SetLineColor(ROOT.kOrange)
    spline_yellow_down.SetFillColor(ROOT.kOrange)

    spline_green_up.SetLineColor(ROOT.kGreen + 1)
    spline_green_up.SetFillColor(ROOT.kGreen + 1)
    spline_green_down.SetLineColor(ROOT.kGreen + 1)
    spline_green_down.SetFillColor(ROOT.kGreen + 1)

    # Draw the bands first (filled areas)
    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
    green.SetFillColor(ROOT.kGreen + 1)
    green.SetLineColor(ROOT.kGreen + 1)
    green.SetFillStyle(1001)

    yellow.Draw('F')
    green.Draw('Fsame')

    # Smooth the median and observed lines
    spline_median = TSpline3("spline_median", median)
    spline_observed = TSpline3("spline_observed", observed)

    spline_median.SetLineColor(1)
    spline_median.SetLineWidth(2)
    spline_median.SetLineStyle(2)
    spline_median.Draw('Lsame')

    spline_observed.SetLineColor(2)  # Red color for observed limits
    spline_observed.SetLineWidth(2)
    spline_observed.SetLineStyle(1)
#    spline_observed.Draw('Lsame')

    CMS_lumi.CMS_lumi(c, 14, 11)
    ROOT.gPad.SetTicks(1, 1)
    frame.Draw('sameaxis')

    x1 = 0.55
    x2 = x1 + 0.24
    y2 = 0.76
    y1 = 0.60
    legend = TLegend(x1, y1, x2, y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)
    legend.AddEntry(spline_median, "Asymptotic CL_{s} expected", 'L')
    legend.AddEntry(green, "#pm 1 std. deviation", 'f')
    legend.AddEntry(yellow, "#pm 2 std. deviation", 'f')
#    legend.AddEntry(spline_observed, "Observed", 'L')
    legend.Draw()

    print(" ")
    c.SaveAs("UpperLimit_zh_final_v2.png")
    c.Close()

# RANGE of floats
def frange(start, stop, step):
    i = start
    while i <= stop:
        yield i
        i += step

# MAIN
def main():
    labels = []
    for label in range(20, 65, 5):
        labels.append(label)

    plotUpperLimits(labels)

if __name__ == '__main__':
    main()
