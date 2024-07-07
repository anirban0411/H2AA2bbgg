#! /bin/bash

# higgs combine commands

datacard="Datacard_mA60_com_test.txt"
output_ws="Datacard_mA60_com_test.root"
prof_name="_mA60_comb_v1"

text2workspace.py $datacard -o $output_ws -m 60

combineTool.py -M Impacts -d $output_ws -m 60 -n impacts_$prof_name --doInitialFit --robustFit 1 -t -1 --redefineSignalPOIs r --rMax 20 --rMin -20 --expectSignal 1  --setRobustFitAlgo=Minuit2,Migrad --X-rtd FITTER_NEW_CROSSING_ALGO --setRobustFitTolerance=0.1 --X-rtd FITTER_NEVER_GIVE_UP --cminFallbackAlgo Minuit2,0:1 --cminDefaultMinimizerStrategy 0 --freezeParameters MH

combineTool.py -M Impacts -d $output_ws -m 60 -n impacts_$prof_name --robustFit 1 --cminDefaultMinimizerStrategy 0 --freezeParameters MH -P r --doFits --expectSignal 1 --rMax 20 --rMin -20 --redefineSignalPOIs r

combineTool.py -M Impacts -d $output_ws -m 60 -n impacts_$prof_name --robustFit 1 --cminDefaultMinimizerStrategy 0 --freezeParameters MH -P r -o impacts_$prof_name.json
plotImpacts.py -i impacts_$prof_name.json -o impacts_$prof_name\_r --POI r
