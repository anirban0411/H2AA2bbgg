# Input config file for running trees2ws

trees2wsCfg = {

  # Name of RooDirectory storing input tree
  'inputTreeDir':'',

  # Variables to be added to dataframe: use wildcard * for common strings
  'mainVars':["CMS_hgg_mass","weight","dZ"], # Vars to add to nominal RooDatasets
  'dataVars':["CMS_hgg_mass","weight","dZ"], # Vars for data workspace (trees2ws_data.py script)
  'stxsVar':'stage1p2bin', # Var for STXS splitting: if using option doSTXSSplitting
  'notagVars':["weight"], # Vars to add to NOTAG RooDataset
  'systematicsVars':["CMS_hgg_mass","weight","dZ"], # Variables to add to sytematic RooDataHists
  'theoryWeightContainers':{}, # Theory weights to add to nominal + NOTAG RooDatasets, value corresponds to number of weights (0-N)

  # List of systematics: use string YEAR for year-dependent systematics
#  'systematics':['JetBTagReshapeWeight','TriggerWeight','JECAbsoluteYEAR','JER','PUJIDShift'],
  'systematics':['PU','leptonID','photonID','btag','Trig','PUJetID','JEC','JER'],

  # Analysis categories: python list of cats or use 'auto' to extract from input tree
  'cats':'auto'

}
