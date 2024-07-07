    cmsrel CMSSW_10_2_3

    cd CMSSW_10_2_3/src

    cmsenv

    git clone https://github.com/lviliani/LHEWeightsReader.git

    scram b

    cd LHEWeightsReader/test/

    cmsRun lheWeightReader_cfg.py inputFiles="MINIAOD_FILE"
