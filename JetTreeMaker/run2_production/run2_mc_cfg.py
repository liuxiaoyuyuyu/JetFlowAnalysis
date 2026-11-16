import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag
from FWCore.ParameterSet.VarParsing import VarParsing


options = VarParsing()
options.register(
    "inputFiles",
    [],
    VarParsing.multiplicity.list,
    VarParsing.varType.string,
    "Input MiniAOD files",
)
options.register(
    "maxEvents",
    -1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Maximum number of events",
)
options.register(
    "outputFile",
    "output_UL_ak8.root",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "TFileService output file",
)
options.register(
    "globalTag",
    "106X_upgrade2018_realistic_v16_L1v1",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Run-2 UL MC global tag",
)
options.parseArguments()


process = cms.Process("Analysis")
process.load("Configuration.Geometry.GeometryDB_cff")
process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.GlobalTag = GlobalTag(process.GlobalTag, options.globalTag, "")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
)
process.TFileService = cms.Service(
    "TFileService",
    fileName=cms.string(options.outputFile),
)

# These settings reproduce the 2022 Run-2 custom-tree content.  The analyzer
# records HLT decisions but the path does not filter on them, so all input
# events contribute to the global q-hat normalization.
process.analyzer = cms.EDAnalyzer(
    "TrackAnalyzer",
    minJetPt=cms.untracked.double(350.0),
    maxJetEta=cms.untracked.double(2.4),
    trackPtMin=cms.untracked.double(0.01),
    doTrack=cms.untracked.bool(True),
    doGen=cms.untracked.bool(True),
    vertexSrc=cms.InputTag("offlineSlimmedPrimaryVertices"),
    packedCandSrc=cms.InputTag("packedPFCandidates"),
    lostTracksSrc=cms.InputTag("lostTracks"),
    beamSpotSrc=cms.untracked.InputTag("offlineBeamSpot"),
    jets2=cms.InputTag("slimmedJetsAK8"),
    genEvtInfo=cms.InputTag("generator"),
    packedGen=cms.InputTag("packedGenParticles"),
    genJets=cms.InputTag("slimmedGenJetsAK8"),
    puSummaryInfo=cms.InputTag("slimmedAddPileupInfo"),
)

process.analysis = cms.Path(process.analyzer)

