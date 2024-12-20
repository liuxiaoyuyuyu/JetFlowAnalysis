#!/bin/sh

#source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc700
cd ~xiaoyul/MYDEMOANALYZER/CMSSW_13_3_0/src/

#cmsenv
eval `scramv1 runtime -sh`
cd ~xiaoyul/MYDEMOANALYZER/CMSSW_13_3_0/src/Run2/basicAna/JetTreeAnalyzer/batch
echo PWD: $PWD
echo $HOSTNAME
../bin/run2_default_data_vn.exe $1 0 1
