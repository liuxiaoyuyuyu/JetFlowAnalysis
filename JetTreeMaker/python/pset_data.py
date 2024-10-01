import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi
process = cms.Process('Analysis')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True),
                                   TryToContinue = cms.untracked.vstring('ProductNotFound'))


# input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"file:/eos/cms/store/group/phys_heavyions/davidlw/JetMET0/Jet260_RAWMINIAODSkim/240210_104720/0000/skim_23.root"
        "file:/eos/user/x/xiaoyul/13cdeae8-0c59-401a-a6f9-defef27d998e.root"
        #"file:/afs/cern.ch/user/p/pgardner/StorageArea/Output_all_275_24.root"
        #"file:/afs/cern.ch/user/p/pgardner/CMSSW_10_6_4_patch1/src/HeavyIonsAnalysis/TrackAnalysis/python/puppi/highMultJets.root"
        #"file:highMultJets.root"
        #"file:/eos/cms/store/group/phys_heavyions/flowcorr/JetHT/Ak8Jet500Skim_JetHT_Run2018D/210629_021937/0000/ppRun2UL_MINIAOD_2.root"
        )
    )

# number of events to process, set to -1 to process all events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

# load Global Tag, geometry, etc.
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

# root output
process.TFileService = cms.Service("TFileService",
    #fileName = cms.string("output_UL_ak8.root"))
    fileName = cms.string("MATCH.root"))

process.hlt = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.hlt.HLTPaths = ['HLT_AK8PF*'] # for allphysics
process.hlt.HLTPaths = ['HLT_AK8PFJet500*'] # for allphysics
#500 is the highest jet energy without prescale
process.hlt.andOr = cms.bool(True)
process.hlt.throw = cms.bool(False)
process.hlt500= process.hlt.clone()
process.eventFilterHLT500 = cms.Sequence(process.hlt500)

process.analyzerOffline = cms.EDAnalyzer('TrackAnalyzer',
    minJetPt = cms.untracked.double(100),
    maxJetEta = cms.untracked.double(2.5),
    doTrack = cms.untracked.bool(False),
    trackPtMin = cms.untracked.double(0.01),
    vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
    packedCandSrc = cms.InputTag("packedPFCandidates"),
    jets2 = cms.InputTag('slimmedJetsAK8'),
    lostTracksSrc = cms.InputTag("lostTracks"),
    beamSpotSrc = cms.untracked.InputTag('offlineBeamSpot'),
    #pfjetH = cms.InputTag('hltAK8PFJetsCorrectedMatchedToCaloJets200')
)
#process.analyzer260 = process.analyzer.clone()
#process.analyzerOnline = cms.EDAnalyzer('TrackAnalyzer_jet',
#    vertexSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),
#    packedCandSrc = cms.InputTag("packedPFCandidates"),
#    jets2 = cms.InputTag("hltAK8PFJetsCorrectedMatchedToCaloJets200","HLTX")
#)
#
#main forest sequence
process.runAnalyzer500 = cms.Path(
    #process.analyzerOnline *
    process.eventFilterHLT500 *
    process.analyzerOffline
    )

