# JetFlowAnalysis
Parker Gardner's code for the jet collectivity study

- [Make Jet Trees](#make-jet-trees)
  - [JetTreeMaker](#jettreemaker)
  - [Local test](#local-test)
  - [Crab jobs](#crab-jobs)
- [Analyze Jet Trees](#analyze-jet-trees) 
  - [JetTreeAnalyzer](#jettreeanalyzer)
  - [Batch jobs](#batch-jobs)
   
## Make Jet Trees
### JetTreeMaker
```Linux
ssh -XY xiaoyul@lxplus.cern.ch
```
CMS release for Run 3 data:
```Linux
cmsrel CMSSW_13_3_0
cd CMSSW_13_3_0/src
cmsenv
```  
>[!Note] 
>`cmsenv` in the CMSSW_13_3_0/src directory will set the specific enviroment for the CMS release CMSSW_13_3_0.
Compile JetTreeMaker
```Linux
git clone https://github.com/liuxiaoyuyuyu/JetFlowAnalysis/
cd JetFlowAnalysis/JetTreeMaker
scram b -j8
```  
>[!Note] 
>`scram b` only compile codes two directories down CMSSW_*/src. In this case: CMSSW_13_3_0/src/JetFlowAnalysis/JetTreeMaker/.

### Local test
```Linux
cd JetFlowAnalysis/JetTreeMaker/python
cmsRun pset.py
```  
### Crab jobs 
```Linux
voms-proxy-init -voms cms
crab submit -c crab.py
```
>[!Note] 
>1. "pset.py" is the macro for local testing, it is also the configuration that will be used by the crab jobs. Therefore, "psetName" in "crab.py" should match "pset.py". Output name should also match in "pset.py" and "crab.py".
>2. Modify "config.Data.outLFNDirBase" in "crab.py" if needed. It is currently set to my eos space in flowcorr. 
>3. InputDataset:
>[Run 3 data](https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#Year_2024)
>For the latest data to use, refer to the "Analysis Summary Table" section of each year.
>5. Set tiggers by "process.hlt.HLTPaths" in "pset.py".
>6. Check task status: `crab status -d`.

| Year   | Trigger   |
|------------|------------|
| 2022| HLT_AK8PFJet500*| 
| 2023| HLT_AK8PFJet500*| 
| 2024|HLT_AK8PFJet500_v27, HLT_AK8PFJet275_Nch40_v5, HLT_AK8PFJet260_v27 | 

How to check tiggers:
1. Go to [CMS OMS](https://cmsoms.cern.ch/cms/run_3/index)
2. Select a recent phyiscs run
3. In the upper bar, hover over "Runs", then go to "Triggers->HLT trigger rates"
4. Find triggers in the "HLT Path Rates" section.

## Analyze Jet Trees
### JetTreeAnalyzer
Main analysis macro in: JetTreeAnalyzer/src/

Header files in: JetTreeAnalyzer/include/ 
<!--
Macros in src/
Main macro( list of files, job number 1-N)
    loads root file
    2PC, jet multiplicity, save histograms
    Line 136-137 comment out
        MC corrections
    Line 141-144 comment out
        HLT efficiency

    Line 79 
    Line 149 Main code starts

Header file
    include/
    3 header files: coordinate tools, constants, Tree details
-->    
Once changes are made to macro and header, back up one directory so you can see both \src and \include. 

Modify "Makefile": change the name to reflect the macro you are working on. (change it in five locations).

Remove old executable with `make clean`. 

Compile new changes with `make` (very fast, unlike scram b). It will produce excutable in a \bin directory.

### Batch jobs
Make data file list that contains all paths to all files (these files were produced by the JetTreeMaker in the last step.)
e.g
```
cd /eos/cms/store/group/phys_heavyions/flowcorr/Run3_jet_trees
find $PWD/JetMET* -name '*root'> TreeList_Run2023Dec_MINIAOD.list
```
Split the file list:
```
cd /afs/cern.ch/user/x/xiaoyul/MYDEMOANALYZER/CMSSW_13_3_0/src/JetFlowAnalysis/JetTreeAnalyzer/batch/
cp /eos/cms/store/group/phys_heavyions/flowcorr/Run3_jet_trees/TreeList_Run2023Dec_MINIAOD.list .
mkdir run3_data
mkdir run3_data/list_25
cd run3_data/list_25
split -l25 -d -a 3 ../../TreeList_Run2023Dec_MINIAOD.list list_job
```

Modify "JetTreeAnalyzer/batch/data_vn500.sh" to reflect the excutable and files that will be used. 
```
chmod +x data_vn500.sh
```

Modify "JetTreeAnalyzer/batch/OnOff.py" accordingly (e.g excutable name, file list name, line_count_frac=line_count/X, X needs to be the same as split -lX). This python script will divide jobs, give args etc. 

Submit condor jobs:
```Linux
python OnOff.py
```


