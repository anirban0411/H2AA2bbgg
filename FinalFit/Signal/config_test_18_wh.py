# Config file: options for signal fitting

#_year = 'combined'
_year = '2018'

signalScriptCfg = {
  
  # Setup
  'inputWSDir':'/afs/cern.ch/work/a/abala/private/sig_model/WH_18/ws_GG2H',
  'procs':'auto', # if auto: inferred automatically from filenames
  'cats':'auto', # if auto: inferred automatically from (0) workspace
  'ext':'test_%s'%_year,
  'analysis':'example', # To specify which replacement dataset mapping (defined in ./python/replacementMap.py)
  'year':'%s'%_year, # Use 'combined' if merging all years: not recommended
  'massPoints':'20,25,30,35,40,45,50,55,60',
#  'massPoints':'20',

  #Photon shape systematics  
#  'scales':'HighR9EB,HighR9EE,LowR9EB,LowR9EE,Gain1EB,Gain6EB', # separate nuisance per year
#  'scalesCorr':'MaterialCentralBarrel,MaterialOuterBarrel,MaterialForward,FNUFEE,FNUFEB,ShowerShapeHighR9EE,ShowerShapeHighR9EB,ShowerShapeLowR9EE,ShowerShapeLowR9EB', # correlated across years
#  'scalesGlobal':'NonLinearity,Geant4', # affect all processes equally, correlated across years
#  'smears':'HighR9EBPhi,HighR9EBRho,HighR9EEPhi,HighR9EERho,LowR9EBPhi,LowR9EBRho,LowR9EEPhi,LowR9EERho', # separate nuisance per year
  
  'scales':'', # separate nuisance per year
  'scalesCorr' :'',
  'scalesGlobal' : '',
  'smears' : '',
  # Job submission options
  'batch':'local', # ['condor','SGE','IC','local']
  'queue':'espresso',
}
