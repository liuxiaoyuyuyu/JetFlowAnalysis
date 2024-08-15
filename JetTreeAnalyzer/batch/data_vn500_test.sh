#!/bin/sh

#source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc700
cd ~xiaoyul/MYDEMOANALYZER/CMSSW_13_3_0/src/

#cmsenv
eval `scramv1 runtime -sh`
#cmssw-el7
cd ~xiaoyul/MYDEMOANALYZER/CMSSW_13_3_0/src/JetFlowAnalysis/JetTreeAnalyzer/batch
echo PWD: $PWD
echo $HOSTNAME
#../bin/new_default_data_1d2d.exe ./run3_data/list_25/list_job$1 0 1
#../bin/data_QA.exe ./run3_data/list_25/list_job$1 0 1
#../bin/data_QA.exe ./run3_data/list_25_hlt260/list_job$1 0 1
#../bin/data_statistics.exe ./run3_data/list_25/list_job$1 0 1
../bin/new_default_data_vn.exe ./Run3_tree_list/2023/list_25/list_job$1 0 1
