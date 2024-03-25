# JetFlowAnalysis
Parker Gardner's code for the jet collectivity study

- [Make Jet Trees](#make-jet-trees)
  - [JetTreeMaker](#treemaker)
  - [Local test](#local-test)
  - [Crab jobs](#crab-jobs)
- [Analyze Jet Trees](#anslyze-jet-trees) 
  - [JetTreeAnalyzer](#treeanalyzer)
  - [Batch jobs](#batch-jobs)
   
## Make Jet Trees
### JetTreeMaker
Parker's codes only work with lxplus7** machines for now:
```Linux
ssh -XY xiaoyul@lxplus7.cern.ch
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
>`scram b` only compile codes two directories down CMSSW_*/src. In this case: CMSSW_13_3_0/src/JetFlowAnalysis/JetTreeMaker/; on my lxplus (when learning this anslysis from Parker): CMSSW_13_3_0/src/Dir1/Dir2/

### Local test
```Linux
cd JetFlowAnalysis/JetTreeMaker/python
cmsRun pset.py
```  
### Crab jobs 
```Linux
voms-proxy-init -voms cms
crab submit -c crab.py --dryrun
```
>[!Note] 
>1. "pset.py" is the macro for local testing, it is also the configuration that will be used by the crab jobs. Therefore, "psetName" in "crab.py" should match "pset.py".
>2. Modify "config.Data.outLFNDirBase" in "crab.py" if needed. It is currently set to my eos space in flowcorr. 
>3. InputDataset:
> https://cmsweb.cern.ch/das/request?input=dataset%3D%2FJetMET0%2FRun2023C-22Sep2023_v4-v1%2FMINIAOD&instance=prod/global
>4. Monitor crab jobs:
>check task status: `crab status`;
>online monitor: https://monit-grafana.cern.ch/d/cmsTMDetail/cms-task-monitoring-task-view?from=1709757014000&orgId=11&to=now&var-task=240306_213014%3Axiaoyul_crab_20240306_223010&var-user=xiaoyul  

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

Make data file list that contains all paths to all files, e.g 
(Those files were prodeuced by the JetTreeMaker in the last step.)

Split the file list with "splitfiles" in JetTreeAnalyzer/batch/, usually it only needs to be done once.

### Batch jobs
Modify "JetTreeAnalyzer/batch/data_vn500.sh" to reflect the excutable and files that will be used. 

Modify "JetTreeAnalyzer/batch/OnOff.py" accordingly. This python script will divide jobs, give args etc. 

Submit condor jobs:
```Linux
python OnOff.py
```


