import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend, TColor
import CMS_lumi, tdrstyle
import subprocess  # to execute shell command
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# CMS style
#CMS_lumi.cmsText = "CMS"
#CMS_lumi.extraText = "Preliminary"
#CMS_lumi.cmsTextSize = 0.65
#CMS_lumi.outOfFrame = True
#CMS_lumi.relPosX = 0.13
#tdrstyle.setTDRStyle()

# GET limits from root file
def getLimits(file_name):
    file = TFile(file_name)
    tree = file.Get("limit")

    limits = []
    for quantile in tree:
        limits.append(tree.limit)
    return limits[:6]

# PLOT upper limits
def plotUpperLimits(labels):
    N = len(labels)
    yellow = TGraph(2 * N)    # ±2σ band
    green  = TGraph(2 * N)    # ±1σ band
    median = TGraph(N)        # expected (median)

    up2s = []

    # Calculate expected limits at 1 GeV intervals
    print("Reading expected limits from files:")
    for i, mass in enumerate(labels):
        file_name = f"higgsCombine.HToAATo2B2G_m{mass}_com.AsymptoticLimits.mH{mass}.root"
        try:
            limit = getLimits(file_name)
            # order assumed: [-2σ, -1σ, median, +1σ, +2σ, observed]
            yellow.SetPoint(i,             mass, limit[4])   # +2σ
            green.SetPoint(i,              mass, limit[3])   # +1σ
            median.SetPoint(i,             mass, limit[2])   # median
            green.SetPoint(2 * N - 1 - i,  mass, limit[1])   # -1σ
            yellow.SetPoint(2 * N - 1 - i, mass, limit[0])   # -2σ
            up2s.append(limit[4])
            print(f"Mass {mass}: Median Exp = {limit[2]:.4f}")
        except Exception as e:
            print(f"Warning: Could not read expected file for mass {mass}: {e}")
            yellow.SetPoint(i,             mass, 0)
            green.SetPoint(i,              mass, 0)
            median.SetPoint(i,             mass, 0)
            green.SetPoint(2 * N - 1 - i,  mass, 0)
            yellow.SetPoint(2 * N - 1 - i, mass, 0)
            up2s.append(0)

    # Read observed limits at 1 GeV intervals
    N_obs = N
    obspts = TGraph(N_obs)

    obs_vals, obs_errors = [], []
    print("\nReading observed limits from files:")

    for i, mass in enumerate(labels):
        file_name = f"obs/higgsCombine.HToAATo2B2G_m{mass}_com.AsymptoticLimits.mH{mass}.root"
        try:
            limit = getLimits(file_name)
            observed_val = limit[5]  # observed limit

            # Calculate 1-sigma error for observed based on expected ±1σ band
            expected_err = (limit[3] - limit[1]) / 2.0  # Half of the ±1σ range

            obspts.SetPoint(i, mass, observed_val)
            #obspts.SetPointError(i, 0, expected_err)  # x_err=0, y_err=expected_err

            obs_vals.append(observed_val)
            obs_errors.append(expected_err)

            print(f"Mass {mass}: Observed = {observed_val:.4f}, Error = {expected_err:.4f}")

        except Exception as e:
            print(f"Warning: Could not read observed file for mass {mass}: {e}")
            obspts.SetPoint(i, mass, 0)
            obs_vals.append(0)
            obs_errors.append(0)

    # canvas & frame
    W, H = 800, 600
    T, B, L, R = 0.08 * H, 0.12 * H, 0.12 * W, 0.04 * W
    c = TCanvas("c", "c", 100, 100, W, H)
    c.SetLeftMargin(L / W); c.SetRightMargin(R / W)
    c.SetTopMargin(T / H);  c.SetBottomMargin(B / H)
    #c.SetGrid(); c.cd()

    # Calculate y_max including both expected and observed values
    valid_obs_vals   = [val for val in obs_vals if val > 0]
    valid_obs_errors = [obs_vals[i] + obs_errors[i]
                        for i in range(len(obs_vals)) if obs_vals[i] > 0]

    if valid_obs_errors and up2s:
        y_max = max(max(up2s), max(valid_obs_errors)) * 2
    elif up2s:
        y_max = max(up2s) * 2
    else:
        y_max = 1.0

    frame = c.DrawFrame(min(labels) - 1, 0.0, max(labels) + 1, y_max)
    frame.GetYaxis().SetTitle(
        "[#sigma(VH) #times B(H #rightarrow AA #rightarrow bb#gamma#gamma)]/#sigma(VH)_{SM}")
    frame.GetXaxis().SetTitle("M_{A} (GeV)")
    frame.GetXaxis().SetLimits(min(labels), max(labels))
    frame.SetMinimum(0)
    frame.SetMaximum(0.5)

    # bands
    #yellow.SetFillColor(ROOT.kOrange); yellow.SetLineColor(ROOT.kOrange); yellow.SetFillStyle(1001)
    #green.SetFillColor(ROOT.kGreen + 1); green.SetLineColor(ROOT.kGreen + 1); green.SetFillStyle(1001)
    #yellow.Draw("F"); green.Draw("F SAME")

    yellow.SetFillColor(TColor.GetColor("#85D1FBff")); yellow.SetLineColor(TColor.GetColor("#85D1FBff")); yellow.SetFillStyle(1001)
    green.SetFillColor(TColor.GetColor("#FFDF7Fff"));  green.SetLineColor(TColor.GetColor("#FFDF7Fff"));  green.SetFillStyle(1001)
    yellow.Draw("F"); green.Draw("F SAME")

    # expected: dashed black line with straight connections (no smooth spline)
    median.SetLineColor(ROOT.kBlack)
    median.SetLineStyle(2)
    median.SetLineWidth(2)
    median.Draw("L SAME")  # straight lines between points

    # observed: draw with straight lines connecting points (no interpolation)
    obspts.SetMarkerStyle(20)
    obspts.SetMarkerSize(0.8)
    obspts.SetMarkerColor(ROOT.kBlack)
    obspts.SetLineColor(ROOT.kBlack)
    obspts.SetLineWidth(2)
    obspts.Draw("LP SAME")  # line + points + error bars (straight lines)

    ROOT.gPad.SetTicks(1, 1)

    # legend
    legend = TLegend(0.55, 0.58, 0.82, 0.82)
    legend.SetFillStyle(0); legend.SetBorderSize(0); legend.SetTextSize(0.041); legend.SetTextFont(42)
    legend.AddEntry(median, "Median expected", 'L')
    legend.AddEntry(green,  "Expected 68% CL",  'f')
    legend.AddEntry(yellow, "Expected 95% CL",  'f')
    legend.AddEntry(obspts, "Observed",          "LP")
    legend.Draw()

    # CMS prelim text
    latex = ROOT.TLatex(); latex.SetNDC(); latex.SetTextFont(42); latex.SetTextSize(0.05)
    latex.DrawLatex(0.145, 0.835, "#bf{CMS} #it{Preliminary}")
#    latex.DrawLatex(0.145, 0.835, "#bf{CMS}")
    latex.SetTextAlign(31); latex.DrawLatex(0.95, 0.935, "137 fb^{-1} (13 TeV)")
    frame.Draw("sameaxis")

    # Print summary statistics
    print(f"\nSummary:")
    print(f"Masses processed (1 GeV interval): {labels}")
    print(f"Total points: {N}")
    print(f"Total observed points with valid data: {len(valid_obs_vals)}")

    c.SaveAs("UpperLimit_vh_obs_br_paper.png")
    c.SaveAs("UpperLimit_vh_obs_br_paper.pdf")
    c.Close()


# MAIN
def main():
    labels = list(range(20, 61))
    plotUpperLimits(labels)

if __name__ == '__main__':
    main()
