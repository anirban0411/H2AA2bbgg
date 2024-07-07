# Config file: options for signal fitting

backgroundScriptCfg = {
  
  # Setup
  'inputWSDir':'/afs/cern.ch/work/a/abala/private/bkg_model/ZH_comb/ws', # location of 'allData.root' file
  'cats':'auto', # auto: automatically inferred from input ws
  'catOffset':0, # add offset to category numbers (useful for categories from different allData.root files)  
  'ext':'test_combined', # extension to add to output directory
  'year':'combined', # Use combined when merging all years in category (for plots)
#  'year':'2018',

  # Job submission options
  'batch':'local', # [condor,SGE,IC,local]
  'queue':'hep.q' # for condor e.g. microcentury
  
}
