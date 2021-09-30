Instruction to run the L2TauNNTag test
cd CMSSW_12_1_0_pre3/src
cmsenv
scram b -j 10
cmsRun testL2TauTagNN.py

if you want to run on a limited number of events :
cmsRun testL2TauTagNN.py maxEvents=100
