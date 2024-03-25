#include "JetFlowAnalysis/JetTreeMaker/interface/reduced.h"
#include "JetFlowAnalysis/JetTreeMaker/interface/coordinateTools.h"
#include "math.h"
#include <iostream>
#include <iomanip>
#include "FWCore/Framework/interface/Event.h"
using TMath::ATan;
using TMath::Exp;

TrackAnalyzer_jet::TrackAnalyzer_jet(const edm::ParameterSet& iConfig)
{
    // minJetPt = iConfig.getUntrackedParameter<double>("minJetPt",400);
    // maxJetEta = iConfig.getUntrackedParameter<double>("maxJetEta",1.6);

    // trackPtMin_ = iConfig.getUntrackedParameter<double>("trackPtMin",0.01);

    packedCandLabel_ = iConfig.getParameter<edm::InputTag>("packedCandSrc");
    packedCandSrc_ = consumes<edm::View<pat::PackedCandidate>>(packedCandLabel_);
    vertexSrcLabel_ = iConfig.getParameter<edm::InputTag>("vertexSrc");
    vertexSrc_ = consumes<reco::VertexCollection>(vertexSrcLabel_);

    // lostTracksLabel_ = iConfig.getParameter<edm::InputTag>("lostTracksSrc");
    // lostTracksSrc_ = consumes<edm::View<pat::PackedCandidate>>(lostTracksLabel_);

    //beamSpotProducer_ = consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc",edm::InputTag("offlineBeamSpot")));


    //tok_triggerResults_ = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults::HLT"));

    jets2Token_      = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets2"));
}
//--------------------------------------------------------------------------------------------------
TrackAnalyzer_jet::~TrackAnalyzer_jet()
{
}
//--------------------------------------------------------------------------------------------------
void TrackAnalyzer_jet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    nEv   = (int)iEvent.id().event();
    nRun  = (int)iEvent.id().run();
    nLumi = (int)iEvent.luminosityBlock();

    clearVectors();

    //***********************************************
    //***********************************************
    //***********************************************
    //maybe delete this section??
    /*
       edm::Handle<edm::TriggerResults> trigResults; //our trigger result object
       edm::InputTag trigResultsTag("TriggerResults","","HLT"); //make sure have correct process on MC
    //data process=HLT, MC depends, Spring11 is REDIGI311X
    iEvent.getByLabel(trigResultsTag,trigResults);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResults);   

    if(trigNames.size() >0){
    for(int i = 1; i<500; i++){
    //HLT_AK8PFJet400
    std::string pathName="HLT_AK8PFJet400_v"+std::to_string(i);
    bool passTrig1=(trigNames.triggerIndex(pathName) < trigNames.size() && trigResults->wasrun(trigNames.triggerIndex(pathName)) && trigResults->accept(trigNames.triggerIndex(pathName))); 
    MAINpassTrig1=passTrig1;
    if(passTrig1) break;
    }
    }
    //if(trigNames.size() >1){
    if(trigNames.size() >0){
    for(int i = 1; i<500; i++){
    std::string pathName="HLT_AK8PFJet500_v"+std::to_string(i);
    bool passTrig2=(trigNames.triggerIndex(pathName) < trigNames.size() && trigResults->wasrun(trigNames.triggerIndex(pathName)) && trigResults->accept(trigNames.triggerIndex(pathName))); 
    MAINpassTrig2=passTrig2;
    if(passTrig2) break;
    }
    }
    //if(trigNames.size() >1){
    if(trigNames.size() >0){
    for(int i = 1; i<500; i++){
    std::string pathName="HLT_AK8PFJet550_v"+std::to_string(i);
    bool passTrig3=(trigNames.triggerIndex(pathName) < trigNames.size() && trigResults->wasrun(trigNames.triggerIndex(pathName)) && trigResults->accept(trigNames.triggerIndex(pathName))); 
    MAINpassTrig3=passTrig3;
    if(passTrig3) break;
    }
    }

*/

    fillJets2(iEvent);

    trackTree_->Fill();
}
//--------------------------------------------------------------------------------------------------
void TrackAnalyzer_jet::fillJets2(const edm::Event& iEvent) {


    const reco::VertexCollection* recoVertices;
    edm::Handle<reco::VertexCollection> vertexCollection;
    iEvent.getByToken(vertexSrc_,vertexCollection);
    recoVertices = vertexCollection.product();

    edm::Handle<pat::JetCollection> jets2;
    iEvent.getByToken(jets2Token_, jets2);
    edm::Handle<edm::View<pat::PackedCandidate>> cands;
    iEvent.getByToken(packedCandSrc_,cands);

    int passer = 0;
    for (const pat::Jet &j :  *jets2) {
        if (!j.isPFJet()) continue;
        if (j.pt() < 20 || fabs(j.eta()) > 2.5) continue;
        if (j.chargedMultiplicity() < 1) continue;
        jetPt.push_back(j.pt());
        jetEta.push_back(j.eta());
        jetPhi.push_back(j.phi());
        jetNumDaughters.push_back(j.numberOfDaughters());
        chargedMultiplicity.push_back(j.chargedMultiplicity());

        std::vector<float>		    vPuppiW;
        std::vector<int>		    vcharge;
        std::vector<unsigned int>	vVertRef;
        std::vector<float>		    vpt;
        std::vector<float>      vptError;
        std::vector<float>		veta;
        std::vector<float>		vphi;
        std::vector<float>      vtheta;

        std::vector<float>      trkxysig;
        std::vector<float>      trkzsig;

        std::vector<float>		vdauVZ;
        std::vector<float>		vdauVY;
        std::vector<float>		vdauVX;

        for( unsigned int dID=0; dID < j.numberOfDaughters();  dID++){
            const pat::PackedCandidate &dau = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(dID));
            vPuppiW.push_back(  dau.puppiWeight());
            vcharge.push_back(	dau.charge());
            vpt.push_back(		dau.pt());

            if(dau.hasTrackDetails()){
                reco::Track const& t = dau.pseudoTrack();
                vptError.push_back( t.ptError() );
                math::XYZPoint v( recoVertices->at(0).position().x(), recoVertices->at(0).position().y(), recoVertices->at(0).position().z() );
                //trkDzFirstVtx.push_back( t.dz( v ) );
                //trkDzErrFirstVtx.push_back( sqrt( t.dzError()*t.dzError() + recoVertices->at(0).zError() * recoVertices->at(0).zError() ) );
                //trkDxyFirstVtx.push_back( t.dxy( v ) );
                //trkDxyErrFirstVtx.push_back( sqrt( t.dxyError()*t.dxyError() + recoVertices->at(0).xError() * recoVertices->at(0).yError() ) );
                trkzsig.push_back(t.dz( v )/(sqrt( t.dzError()*t.dzError() + recoVertices->at(0).zError() * recoVertices->at(0).zError() )));
                trkxysig.push_back(t.dxy( v )/(sqrt( t.dxyError()*t.dxyError() + recoVertices->at(0).xError() * recoVertices->at(0).yError() )));
            } else {
                vptError.push_back(-1);
                //trkDzFirstVtx.push_back(-999);
                //trkDzErrFirstVtx.push_back( 1 );
                //trkDxyFirstVtx.push_back(-999);
                //trkDxyErrFirstVtx.push_back( 1 );
                trkzsig.push_back(-999);
                trkxysig.push_back(-999);
            }

            veta.push_back(		dau.eta());
            vphi.push_back(		dau.phi());
            vtheta.push_back(       dau.theta());

            float dauVZ    = dau.vertex().z();
            float dauVY    = dau.vertex().y();
            float dauVX    = dau.vertex().x();

            vdauVZ.push_back(dauVZ);
            vdauVY.push_back(dauVY);
            vdauVX.push_back(dauVX);
        }
        dau_PuppiW.push_back(	vPuppiW);
        dau_chg.push_back(	    vcharge);
        dau_pt.push_back(       vpt);
        dau_ptError.push_back(  vptError);
        dau_eta.push_back(      veta);
        dau_phi.push_back(      vphi);
        dau_theta.push_back(    vtheta);
        dau_XYDCAsig.push_back( trkxysig   );
        dau_ZDCAsig.push_back(  trkzsig   );

        dau_vz.push_back(       vdauVZ);
        dau_vy.push_back(       vdauVY);
        dau_vx.push_back(       vdauVX);

        passer = passer +1;

    }
    jetN = passer;


}



// ------------ method called once each job just before starting event loop  ------------
void TrackAnalyzer_jet::beginJob()
{
    trackTree_ = fs->make<TTree>("trackTree","v1");
    //jetTree_ = fs->make<TTree>("jetTree","v1");

    // event
    trackTree_->Branch("nRun",&nRun,"nRun/I");
    trackTree_->Branch("nEv",&nEv,"nEv/I");
    trackTree_->Branch("nLumi",&nLumi,"nLumi/I");


    // Jets
    trackTree_->Branch("jetNumDaughters",&jetNumDaughters);
    trackTree_->Branch("jetEta",&jetEta);
    trackTree_->Branch("jetPt",&jetPt);
    trackTree_->Branch("jetPhi",&jetPhi);
    trackTree_->Branch("chargedMultiplicity",&chargedMultiplicity);
    trackTree_->Branch("jetN",&jetN);

    trackTree_->Branch("dau_PuppiW",		&dau_PuppiW); 
    trackTree_->Branch("dau_chg",		&dau_chg); 
    trackTree_->Branch("dau_pt",		&dau_pt);
    trackTree_->Branch("dau_ptError",     &dau_ptError);
    trackTree_->Branch("dau_eta",		&dau_eta);	 
    trackTree_->Branch("dau_phi",		&dau_phi );
    trackTree_->Branch("dau_theta",	&dau_theta);
    trackTree_->Branch("dau_XYDCAsig",    &dau_XYDCAsig);
    trackTree_->Branch("dau_ZDCAsig",     &dau_ZDCAsig);

    trackTree_->Branch("dau_vz",		&dau_vz	 );
    trackTree_->Branch("dau_vy",		&dau_vy	 );
    trackTree_->Branch("dau_vx",		&dau_vx	 );

}

// ------------ method called once each job just after ending the event loop  ------------
void
TrackAnalyzer_jet::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer_jet);
