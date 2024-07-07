import ROOT

def create_datacard(data_histogram, signal_histogram, background_histogram, output_file):
    """
    Creates a datacard for CMS analysis with one background histogram.
    
    Args:
        data_histogram (str): Path to the data histogram ROOT file.
        signal_histogram (str): Path to the signal histogram ROOT file.
        background_histogram (str): Path to the background histogram ROOT file.
        output_file (str): Path to the output datacard text file.
    """
    data_file = ROOT.TFile.Open(data_histogram)
    signal_file = ROOT.TFile.Open(signal_histogram)
    bkg_file = ROOT.TFile.Open(background_histogram)
    
    data_hist = data_file.Get("CMS_hgg_mass")
    signal_hist = signal_file.Get("CMS_hgg_mass")
    bkg_hist = bkg_file.Get("CMS_hgg_mass")

    data_obs = data_hist.Integral()
    signal_yield = signal_hist.Integral()
    bkg_yield = bkg_hist.Integral()

    with open(output_file, 'w') as f:
        f.write("imax 1  number of channels\n")
        f.write("jmax 1  number of backgrounds\n")
        f.write("kmax *  number of nuisance parameters (sources of systematical uncertainties)\n\n")
        
        f.write("shapes * * FAKE\n\n")
        
        f.write("bin             bin1\n")
        f.write("observation     {}\n\n".format(data_obs))
        
        f.write("bin             bin1    bin1\n")
        f.write("process         signal  bkg\n")
        f.write("process         0       1\n")
        f.write("rate            {}  {}\n\n".format(signal_yield, bkg_yield))
        
        # Add any nuisance parameters here if needed
        # f.write("lumi    lnN    1.025    -    -\n")
    
    data_file.Close()
    signal_file.Close()
    bkg_file.Close()

# Example usage:
data_histogram = "zh_template/background/allData_2018.root"
signal_histogram = "zh_template/signal/output_ZHToAA2B2G_M20_2018.root"
background_histogram = "zh_template/background/allData_2018.root"
output_file = "datacard_18_zh_test.txt"

create_datacard(data_histogram, signal_histogram, background_histogram, output_file)
