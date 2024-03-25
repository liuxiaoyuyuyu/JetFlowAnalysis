#!/bin/sh

#source /cvmfs/cms.cern.ch/cmsset_default.sh
#export SCRAM_ARCH=slc7_amd64_gcc700
cd ~pgardner/CMSSW_10_6_4_patch1/src/

#cmsenv
eval `scramv1 runtime -sh`
cd ~pgardner/local_not_scp/CMSSW_10_6_4_patch1/src/HeavyIonAnalysis/TrackAnalysis/batch
echo PWD: $PWD
../bin/tree500.exe ./split_list_ak8_all_DCA_v2/list_cor_$1 0 1
