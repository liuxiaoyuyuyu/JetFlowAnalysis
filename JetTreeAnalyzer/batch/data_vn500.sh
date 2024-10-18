#!/bin/sh

#source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc700
cd ~xiaoyul/MYDEMOANALYZER/CMSSW_13_3_0/src/

#cmsenv
eval `scramv1 runtime -sh`
cd ~xiaoyul/MYDEMOANALYZER/CMSSW_13_3_0/src/JetFlowAnalysis/JetTreeAnalyzer/batch
echo PWD: $PWD
echo $HOSTNAME
#../bin/new_default_data_vn.exe ./Run3_tree_list/2022/list_C_25/list_job$1 0 1
../bin/era_run3_2023_correction.exe ./Run3_tree_list/2023/list_D_25_MC/list_job$1 0 1
