from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.psetName = 'pset_mc.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['MATCH.root']
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 4000
#config.JobType.inputFiles = ['HeavyIonRPVRcd_PbPb2018_offline.db']
config.section_('Data')
#config.Data.inputDataset = '/JetHT/Run2018D-12Nov2019_UL2018-v4/MINIAOD'
#config.Data.inputDataset = '/JetMET0/Run2023B-19Dec2023-v1/MINIAOD'
#config.Data.inputDataset = '/JetMET1/Run2023D-22Sep2023_v2-v1/MINIAOD'
#config.Data.inputDataset = '/QCD_PT-600to800_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM' 
config.Data.inputDataset = '/QCD_PT-600to800_TuneCP5_13p6TeV_pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5_ext1-v2/MINIAODSIM' 

config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader'

#config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
#config.Data.lumiMask = '/afs/cern.ch/user/x/xiaoyul/Cert_Collisions2023_366442_370790_Golden_JSON.txt'
#.txt file convert from raw data of /eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json
#config.Data.lumiMask = '/afs/cern.ch/user/x/xiaoyul/Cert_Collisions2022_355100_362760_Golden_JSON.txt'

config.Data.splitting = 'Automatic'
#config.Data.unitsPerJob = 150
config.Data.publication = False

#config.Data.outLFNDirBase = '/store/user/pgardner/MINIAOD_2018_UL_D_ak8_new'
#config.Data.outLFNDirBase = '/store/group/phys_heavyions/flowcorr/Run3_jet_trees/'

config.Data.outLFNDirBase = '/store/group/phys_heavyions/xiaoyul/Run3_2022_jet_trees/MC_CD/'

#config.Data.totalUnits        = 1 # root file or lumi section for test only

config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

#blacklist,Vanderbuilt 



#config.Site.storageSite = 'T3_US_Rice'
