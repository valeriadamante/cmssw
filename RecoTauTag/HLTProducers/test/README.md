## Tests for tau HLT producers

### Setup to tes L2TauNNTag
You need to have this repository in CMSSW_12_1_0_pre3/src. If you don't have it, the command to run is:
```
export SCRAM_ARCH=slc7_amd64_gcc900
cmsrel CMSSW_12_1_0_pre3
cd CMSSW_12_1_0_pre3/src
cmsenv  
git cms-merge-topic valeriadamante:L2TauTagNN
git cms-addpkg RecoTauTag/HLTProducers
git cms-addpkg HLTrigger/Configuration
scram b -j8
```

You also need the RecoTauTag-Training files. If you don't have it, run the command:
```
git clone git@github.com:valeriadamante/RecoTauTag-TrainingFiles.git RecoTauTag/TrainingFiles
```


### L2TauNNTag test
To get a user configuration file with the paths you need run the following command:
```
hltGetConfiguration /dev/CMSSW_12_0_0/GRun/V6 --paths HLTriggerFirstPath,HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4,HLT_VBF_DoubleTightChargedIsoPFTauHPS20_Trk1_eta2p1_v1,HLTriggerFinalPath  --unprescale --cff
```
In this way all the paths should be included correctly.


You need to authenticate to run the test command since it exploits files that are saved on store. So you need to run the command:
```
voms-proxy-init --rfc --voms cms
```

You are now ready to run the test for L2TauTagNN
```
cd CMSSW_12_1_0_pre3/src
cmsenv
scram b -j 10
cmsRun RecoTauTag/HLTProducers/test/testL2TauTagNN.py
```

if you want to run on a limited number of events :
```
cmsRun RecoTauTag/HLTProducers/test/testL2TauTagNN.py maxEvents=20
```
