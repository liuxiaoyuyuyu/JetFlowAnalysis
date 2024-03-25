#ifndef TRKANALYZER_H
#define TRKANALYZER_H

// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <functional>

// CMSSW user include files
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
//#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmAlgorithm.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "FWCore/Common/interface/Provenance.h"
#include "FWCore/Framework/interface/EventPrincipal.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
//#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/TriggerResults.h" 
#include "FWCore/Common/interface/TriggerNames.h"

// Particle Flow
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//Gen info
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


// Vertex significance
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"

// Root include files
#include "TTree.h"



#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"



class TrackAnalyzer_jet : public edm::one::EDAnalyzer<> {

  public:
    explicit TrackAnalyzer_jet(const edm::ParameterSet&);
    ~TrackAnalyzer_jet() override;

  private:
    void beginJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;

    void fillJets2(const edm::Event& iEvent);
    void clearVectors();

     // ----------member data ---------------------------
     edm::Service<TFileService> fs;

     edm::InputTag packedCandLabel_;
     edm::EDGetTokenT<edm::View<pat::PackedCandidate>> packedCandSrc_;

     edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;


//reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("pfjetH"
     edm::InputTag pfjetH_;
     edm::EDGetTokenT<reco::PFJetCollection> jets1Token_;

     edm::InputTag jets2_;
     edm::EDGetTokenT<pat::JetCollection> jets2Token_;
     
     edm::InputTag vertexSrcLabel_;
     edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

     //edm::InputTag RawSrcLabel_;
     //edm::EDGetTokenT<reco::PFJet> RawSrc_;

     //edm::EDGetTokenT<reco::BeamSpot> beamSpotProducer_;
     //edm::EDGetTokenT< std::vector< PileupSummaryInfo > > puSummary_;

     //double trackPtMin_; 

     //edm::InputTag vertexSrcLabel_;
     //edm::EDGetTokenT<reco::VertexCollection> vertexSrc_;

     //edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
     //edm::EDGetTokenT< std::vector<reco::GenJet> > packedGenJetToken_;
     //edm::EDGetTokenT< GenEventInfoProduct > genEvtInfo_;

     //edm::InputTag lostTracksLabel_;
     //edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostTracksSrc_;

     //edm::InputTag chi2MapLabel_;
     //edm::EDGetTokenT<edm::ValueMap<float>> chi2Map_;

     //edm::InputTag chi2MapLostLabel_;
     //edm::EDGetTokenT<edm::ValueMap<float>> chi2MapLost_;


     // Root object
     TTree* trackTree_;

     //Branch entries
     long int nRun;
     bool MAINpassTrig1;//didHLTFire400;
     bool MAINpassTrig2;//didHLTFire400;
     bool MAINpassTrig3;//didHLTFire400;
     long int nEv;
     long int nLumi;

     std::vector< float > online_eta;
     std::vector< float > online_phi;
     std::vector< int   > online_Nch;

     std::vector< float > jetEta;
     std::vector< float > jetPt;
     std::vector< float > jetPhi;
     std::vector< int   > jetNumDaughters;
     std::vector< int   > chargedMultiplicity;
     int jetN;


     std::vector<std::vector<float>>		dau_PuppiW;
     std::vector<std::vector<int>>		    dau_chg;
     std::vector<std::vector<float>>		dau_pt;
     std::vector<std::vector<float>>        dau_ptError;
     std::vector<std::vector<float>>		dau_eta;
     std::vector<std::vector<float>>		dau_phi;
     std::vector<std::vector<float>>      	dau_theta;
     std::vector<std::vector<float>>        dau_ZDCAsig;
     std::vector<std::vector<float>>        dau_XYDCAsig;

     std::vector<std::vector<float>>		dau_vz;
     std::vector<std::vector<float>>		dau_vy;
     std::vector<std::vector<float>>		dau_vx;

    float minJetPt;
    float maxJetEta;

    //edm::EDGetTokenT<edm::TriggerResults> tok_triggerResults_;
};

void TrackAnalyzer_jet::clearVectors(){

  jetEta.clear();
  jetPt.clear();
  jetPhi.clear();
  jetNumDaughters.clear();
  chargedMultiplicity.clear();

dau_PuppiW.clear();
dau_chg.clear();
dau_pt.clear();
dau_ptError.clear();
dau_eta.clear();
dau_phi.clear();
dau_theta.clear();
dau_XYDCAsig.clear();
dau_ZDCAsig.clear();

dau_vz.clear();
dau_vy.clear();
dau_vx.clear();

online_eta.clear();
online_phi.clear();
online_Nch.clear();


}

#endif 
