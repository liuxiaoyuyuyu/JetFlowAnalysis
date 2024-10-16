from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from http.client import HTTPException  # updated import for Python 3
#from WMCore.Configuration import Configuration
#config = Configuration()

# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.workArea = '/afs/cern.ch/user/x/xiaoyul/MYDEMOANALYZER/CMSSW_13_3_0/src/JetFlowAnalysis/JetTreeMaker/crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles = ['']

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.splitting = 'Automatic'
#config.Data.totalUnits = 5000
#config.Data.lumiMask = 'Cert_Collisions2023HI_374288_375823_Golden.json'
#config.Data.runRange = '374288-375823'
config.Data.publication = False
config.JobType.allowUndistributedCMSSW = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
#config.Data.ignoreLocality = True
#config.Site.whitelist = ['T1_US_*','T2_US_*','T1_FR_*','T2_FR_*','T2_CH_CERN','T2_BE_IIHE']
config.Site.storageSite = 'T2_CH_CERN'

def submit(config):
    try:
        crabCommand('submit', config = config, dryrun=False)
    except HTTPException as hte:
        print("Failed submitting task: %s" % (hte.headers))  # updated for Python 3

    except ClientException as cle:
        print("Failed submitting task: %s" % (cle))  # updated for Python 3

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

dataMap = {
#            "HIDoubleMuon": { "PD": "/HIDoubleMuon/HIRun2018A-04Apr2019-v1/AOD", "Units": 5, "Memory": 2500, "RunTime": 1800, "PSet": "PbPbSkimAndTree2018_DiMuContBoth_ZDC_cfg.py" },
#            "HISingleMuon": { "PD": "/HISingleMuon/HIRun2018A-04Apr2019-v1/AOD", "Units": 5, "Memory": 2500, "RunTime": 1800, "PSet": "PbPbSkimAndTree2018_DiMuContBoth_ZDC_cfg.py" },
#            "HIDoubleMuonPsiPeri": { "PD": "/HIDoubleMuonPsiPeri/HIRun2018A-04Apr2019-v1/AOD", "Units": 5, "Memory": 2500, "RunTime": 1800, "PSet": "PbPbSkimAndTree2018_DiMuContBoth_ZDC_cfg.py" },
#            "HIForward": { "PD": "/HIForward/HIRun2018A-04Apr2019-v1/AOD", "Units": 30, "Memory": 1800, "RunTime": 1400, "PSet": "PbPbSkimAndTree2018_DiMuContBoth_ZDC_ALLDIMU_cfg.py" },
            }

dataset_list = [
    #"/HIPhysicsRawPrime1/HIRun2023A-PromptReco-v2/MINIAOD",
    #"/HIPhysicsRawPrime2/HIRun2023A-PromptReco-v2/MINIAOD",
    #"/QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM",
    #"/QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2/MINIAODSIM",
    #"/QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM",
    #"/QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2/MINIAODSIM",
    #"/QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM",
    #"/QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2/MINIAODSIM",
    #"/QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM",
    #"/QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2/MINIAODSIM",
    #"/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM",
    #"/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2/MINIAODSIM",
    #"/QCD_PT-3200_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM",
    #"/QCD_PT-3200_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2/MINIAODSIM",
    
    #"/QCD_PT-470to600_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM",
    #"/QCD_PT-470to600_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2/MINIAODSIM",
    #"/QCD_PT-600to800_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM",
    #"/QCD_PT-600to800_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2/MINIAODSIM",
    #"/QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM",
    #"/QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2/MINIAODSIM",
    #"/QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM",
    #"/QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2/MINIAODSIM",
    #"/QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM",
    #"/QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2/MINIAODSIM",
    #"/QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM",
    #"/QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2/MINIAODSIM",
    "/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM",
    "/QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2/MINIAODSIM",
    "/QCD_PT-3200_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM",
    "/QCD_PT-3200_TuneCP5_13p6TeV_pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6_ext1-v2/MINIAODSIM"
    # Add more dataset names as needed
]


#for i in range(1,2):
#    dataMap[("JetCollectivityRun3_2022CD"+str(i))] = { "PD": ("/HIPhysicsRawPrime"+str(i)+"/HIRun2023A-PromptReco-v2/MINIAOD"), "Units": 25, "Memory": 4000, "RunTime": 2100, "PSet": "pset_mc.py" } # UCC

# Iterating over the list and dynamically updating the dictionary
for i, dataset_name in enumerate(dataset_list, start=1):
    dataMap[("JetCollectivityRun3_2022EF" + str(i))] = {
        "PD": dataset_name,
        "Units": 25,
        "Memory": 4000,
        "RunTime": 2100,
        "PSet": "pset_mc.py"
    }

## Submit the muon PDs
for key, val in dataMap.items():
    #config.General.requestName = 'SoS_'+key+'_HIRun2023_PromptReco_MBZDCOR_test_20240903v2'
    config.Data.inputDataset = val["PD"]
    #config.Data.unitsPerJob = val["Units"]
    #config.JobType.maxMemoryMB = val["Memory"]
    #config.JobType.maxJobRuntimeMin = val["RunTime"]
    config.JobType.psetName = val["PSet"]
    #config.Data.outputDatasetTag = config.General.requestName
    #config.Data.outLFNDirBase = '/store/group/phys_heavyions/xiaoyul/Run3_2022_jet_trees/MC_CD/' 
    config.Data.outLFNDirBase = '/store/group/phys_heavyions/xiaoyul/Run3_2022_jet_trees/MC_EF/' 

    print("Submitting CRAB job for: " + val["PD"])  # updated for Python 3

    submit(config)
