#include "HeavyIonsAnalysis/TrackAnalysis/interface/TrackAnalyzer.h"
#include "HeavyIonsAnalysis/TrackAnalysis/interface/coordinateTools.h"
#include "math.h"
#include <iostream>
#include <iomanip>
//#include "FWCore/Framework/interface/Event.h"
using TMath::ATan;
using TMath::Exp;

TrackAnalyzer::TrackAnalyzer(const edm::ParameterSet& iConfig)
{
    minJetPt = iConfig.getUntrackedParameter<double>("minJetPt",100);
    maxJetEta = iConfig.getUntrackedParameter<double>("maxJetEta",2.5);

    doTrack_ = iConfig.getUntrackedParameter<bool>("doTrack",true);
    doGen = iConfig.getUntrackedParameter< bool >("doGen",false);

    trackPtMin_ = iConfig.getUntrackedParameter<double>("trackPtMin",0.01);

    vertexSrcLabel_ = iConfig.getParameter<edm::InputTag>("vertexSrc");
    vertexSrc_ = consumes<reco::VertexCollection>(vertexSrcLabel_);

    packedCandLabel_ = iConfig.getParameter<edm::InputTag>("packedCandSrc");
    packedCandSrc_ = consumes<edm::View<pat::PackedCandidate>>(packedCandLabel_);

    lostTracksLabel_ = iConfig.getParameter<edm::InputTag>("lostTracksSrc");
    lostTracksSrc_ = consumes<edm::View<pat::PackedCandidate>>(lostTracksLabel_);

    beamSpotProducer_ = consumes<reco::BeamSpot>(iConfig.getUntrackedParameter<edm::InputTag>("beamSpotSrc",edm::InputTag("offlineBeamSpot")));

    if(doGen){
        genEvtInfo_ = consumes< GenEventInfoProduct >(iConfig.getParameter<edm::InputTag>("genEvtInfo"));
        packedGenToken_ = consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packedGen"));
        packedGenJetToken_ = consumes< std::vector< reco::GenJet > >(iConfig.getParameter<edm::InputTag>("genJets"));
        puSummary_ = consumes< std::vector< PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("puSummaryInfo")); 
    }

    tok_triggerResults_ = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults::HLT"));

    //jets1Token_      = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets1"));
    jets2Token_      = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets2"));
    //jets3Token_      = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets3"));
    //jets4Token_      = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets4"));
}
//--------------------------------------------------------------------------------------------------
TrackAnalyzer::~TrackAnalyzer()
{
}
//--------------------------------------------------------------------------------------------------
void TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    nEv   = iEvent.id().event();
    nRun  = iEvent.id().run();
    nLumi = iEvent.luminosityBlock();

    clearVectors();

    //***********************************************
    //***********************************************
    //***********************************************
    //maybe delete this section??

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


/*    std::string pathName="HLT_AK8PFJet500_v*";
    //bool passTrig1 = (trigNames.triggerIndex(pathName) < trigNames.size());
    //bool passTrig2=0;
    bool passTrig3=0;
    //if(passTrig1){
    //cout << "good 1" << endl;
    //passTrig2 = (trigResults->wasrun(trigNames.triggerIndex(pathName)));
    //cout << "good 2" << endl;
    if(trigNames.size() >1){
    passTrig3 = (trigResults->accept(trigNames.triggerIndex(pathName))); 
    cout << "good 3" << endl;}
    //}
    MAINpassTrig=(passTrig3);
*/
/*
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByToken(tok_triggerResults_, triggerResults);
    const auto& filterNames = iEvent.triggerNames(*triggerResults);


    const auto& index400 = filterNames.triggerIndex("HLT_AK8PFJet400*");
    bool tempDidHLTFire400 =  (index400 < filterNames.size() && triggerResults->wasrun(index400) && triggerResults->accept(index400));
    didHLTFire400 = tempDidHLTFire400;

    const auto& index500 = filterNames.triggerIndex("HLT_AK8PFJet500*");
    bool tempDidHLTFire500 =  (index500 < filterNames.size() && triggerResults->wasrun(index500) && triggerResults->accept(index500));
    didHLTFire500 = tempDidHLTFire500;

    const auto& index550 = filterNames.triggerIndex("HLT_AK8PFJet550iiiiiii*");
    bool tempDidHLTFire550 =  (index550 < filterNames.size() && triggerResults->wasrun(index550) && triggerResults->accept(index550));
    didHLTFire550 =  tempDidHLTFire550;



    const auto& index550NOT = filterNames.triggerIndex("HLT_AK8PFJet550iiiiiii*");
    bool tempDidHLTFire550NOT =  (index550NOT < filterNames.size() && triggerResults->wasrun(index550NOT) && !triggerResults->accept(index550NOT));
    didHLTFire550NOT =  tempDidHLTFire550NOT;
*/

    fillVertices(iEvent);

    fillJets2(iEvent);
    if(doGen) fillGen(iEvent);  

    if(doTrack_) fillTracks(iEvent, iSetup); 
    trackTree_->Fill();
}
//--------------------------------------------------------------------------------------------------
void TrackAnalyzer::fillVertices(const edm::Event& iEvent) {
    // Fill reconstructed vertices.
    const reco::VertexCollection* recoVertices;
    edm::Handle<reco::VertexCollection> vertexCollection;
    iEvent.getByToken(vertexSrc_,vertexCollection);
    recoVertices = vertexCollection.product();
    unsigned int nVertex = recoVertices->size();
    for (unsigned int i = 0; i < nVertex; ++i) {
        xVtx.push_back( recoVertices->at(i).position().x() );
        yVtx.push_back( recoVertices->at(i).position().y() );
        zVtx.push_back( recoVertices->at(i).position().z() );
        xErrVtx.push_back( recoVertices->at(i).xError() );
        yErrVtx.push_back( recoVertices->at(i).yError() );
        zErrVtx.push_back( recoVertices->at(i).zError() );
        ndofVtx.push_back( recoVertices->at(i).ndof() );
        isFakeVtx.push_back( recoVertices->at(i).isFake() );
        //number of tracks having a weight in vtx fit above 0.5
        nTracksVtx.push_back( recoVertices->at(i).nTracks() );
        float ptSum = 0;
        for( auto ref = recoVertices->at(i).tracks_begin(); ref != recoVertices->at(i).tracks_end(); ref++){
            ptSum += (*ref)->pt();
        }
        ptSumVtx.push_back( ptSum );
    }   
}
//--------------------------------------------------------------------------------------------------
void TrackAnalyzer::fillJets2(const edm::Event& iEvent) {
    edm::Handle<pat::JetCollection> jets2;
    iEvent.getByToken(jets2Token_, jets2);

    const reco::VertexCollection* recoVertices;
    edm::Handle<reco::VertexCollection> vertexCollection;
    iEvent.getByToken(vertexSrc_,vertexCollection);
    recoVertices = vertexCollection.product();

    edm::Handle<edm::View<pat::PackedCandidate>> cands;
    iEvent.getByToken(packedCandSrc_,cands);

    int passer = 0;
    // Let's compute the fraction of charged pt from particles with dz < 0.1 cm
    for (const pat::Jet &j :  *jets2) {
        if (j.pt() < 350 || fabs(j.eta()) > 2.4) continue;

        jetPt.push_back(j.pt());
        jetEta.push_back(j.eta());
        jetPhi.push_back(j.phi());
        jetTheta.push_back(j.theta());
        jetMass.push_back(j.mass());
        jetNumDaughters.push_back(j.numberOfDaughters());
        chargedMultiplicity.push_back(j.chargedMultiplicity());
        muonMultiplicity.push_back(j.muonMultiplicity());

        std::vector<int>		vcharge;
        std::vector<int>	        vpid;
        std::vector<unsigned int>	vVertRef;
        std::vector<float>		vpt;

        std::vector<float>      vptError;
        std::vector<float>      trkDzFirstVtx;
        std::vector<float>      trkDzErrFirstVtx;
        std::vector<float>      trkDxyFirstVtx;
        std::vector<float>      trkDxyErrFirstVtx;
        std::vector<float>      trkxysig;
        std::vector<float>      trkzsig;

        std::vector<float>		veta;
        std::vector<float>		vphi;
        std::vector<float>              vtheta;
        std::vector<float>		vdauVZ;
        std::vector<float>		vdauVY;
        std::vector<float>		vdauVX;
        std::vector<float>		vdauVrefZ;
        std::vector<float>		vdauVrefY;
        std::vector<float>		vdauVrefX;
        std::vector<float>		vp_difZ;
        std::vector<float>		vp_difY;
        std::vector<float>		vp_difX;

        std::vector<double>		vptSTAR;
        std::vector<double>		vetaSTAR;
        std::vector<double>		vphiSTAR;
        std::vector<double>             vthetaSTAR;

        float Ndau_pt_sum = 0;

        for( unsigned int dID=0; dID < j.numberOfDaughters();  dID++){
            const pat::PackedCandidate &dau = dynamic_cast<const pat::PackedCandidate &>(*j.daughter(dID));

            vcharge.push_back(	dau.charge());
            vpid.push_back(		dau.pdgId());
            vpt.push_back(		dau.pt());


            if(dau.hasTrackDetails()){
                reco::Track const& t = dau.pseudoTrack();
                vptError.push_back( t.ptError() );

                math::XYZPoint v( recoVertices->at(0).position().x(), recoVertices->at(0).position().y(), recoVertices->at(0).position().z() );
                trkDzFirstVtx.push_back( t.dz( v ) );
                trkDzErrFirstVtx.push_back( sqrt( t.dzError()*t.dzError() + recoVertices->at(0).zError() * recoVertices->at(0).zError() ) );
                trkDxyFirstVtx.push_back( t.dxy( v ) );
                trkDxyErrFirstVtx.push_back( sqrt( t.dxyError()*t.dxyError() + recoVertices->at(0).xError() * recoVertices->at(0).yError() ) );
                trkzsig.push_back(t.dz( v )/(sqrt( t.dzError()*t.dzError() + recoVertices->at(0).zError() * recoVertices->at(0).zError() )));
                trkxysig.push_back(t.dxy( v )/(sqrt( t.dxyError()*t.dxyError() + recoVertices->at(0).xError() * recoVertices->at(0).yError() )));
            } else {
                vptError.push_back(-1);
                trkDzFirstVtx.push_back(-999);
                trkDzErrFirstVtx.push_back( 1 );
                trkDxyFirstVtx.push_back(-999);
                trkDxyErrFirstVtx.push_back( 1 );
                trkzsig.push_back(-999);
                trkxysig.push_back(-999);
            }

            veta.push_back(		dau.eta());
            vphi.push_back(		dau.phi());
            vtheta.push_back(       dau.theta());

            Ndau_pt_sum = Ndau_pt_sum + dau.pt();
            float dauVZ    = dau.vertex().z();
            float dauVrefZ = recoVertices->at(dau.vertexRef().key()).position().z();
            float dauVY    = dau.vertex().y();
            float dauVrefY = recoVertices->at(dau.vertexRef().key()).position().y();
            float dauVX    = dau.vertex().x();
            float dauVrefX = recoVertices->at(dau.vertexRef().key()).position().x();
            float V_percent_difZ = 100 * fabs(fabs(dauVrefZ/dauVZ) - 1);
            float V_percent_difY = 100 * fabs(fabs(dauVrefY/dauVY) - 1);
            float V_percent_difX = 100 * fabs(fabs(dauVrefX/dauVX) - 1);

            vVertRef.push_back(  	dau.vertexRef().key());
            vdauVZ.push_back(	dauVZ);
            vdauVY.push_back(       dauVY);
            vdauVX.push_back(       dauVX);
            vdauVrefZ.push_back(    dauVrefZ);
            vdauVrefY.push_back(    dauVrefY);
            vdauVrefX.push_back(    dauVrefX);
            vp_difZ.push_back(    	V_percent_difZ);
            vp_difY.push_back(      V_percent_difY);
            vp_difX.push_back(      V_percent_difX);

            double jet_dau_pt    =  ptWRTJet((double)(j.pt()), (double)(j.eta()), (double)(j.phi()), (double)(dau.pt()), (double)(dau.eta()), (double)(dau.phi()));
            double jet_dau_eta   =  etaWRTJet((double)(j.pt()), (double)(j.eta()), (double)(j.phi()), (double)(dau.pt()), (double)(dau.eta()), (double)(dau.phi()));
            double jet_dau_phi   =  phiWRTJet((double)(j.pt()), (double)(j.eta()), (double)(j.phi()), (double)(dau.pt()), (double)(dau.eta()), (double)(dau.phi()));
            double jet_dau_theta = 2*ATan(Exp(-(etaWRTJet((double)(j.pt()), (double)(j.eta()), (double)(j.phi()), (double)(dau.pt()), (double)(dau.eta()), (double)(dau.phi())))));

            vptSTAR.push_back(jet_dau_pt);
            vetaSTAR.push_back(jet_dau_eta);
            vphiSTAR.push_back(jet_dau_phi);
            vthetaSTAR.push_back(jet_dau_theta);

        }
        dau_pt_sum.push_back(   Ndau_pt_sum);

        dau_chg.push_back(	vcharge);
        dau_pid.push_back(      vpid);
        dau_vref.push_back(     vVertRef);

        dau_pt.push_back(       vpt);
        dau_ptError.push_back(       vptError);
        dau_XYDCAsig.push_back( trkxysig   );
        dau_ZDCAsig.push_back(  trkzsig   );
        dau_eta.push_back(      veta);
        dau_phi.push_back(      vphi);
        dau_theta.push_back(    vtheta);

        dau_pt_STAR.push_back(       vptSTAR);
        dau_eta_STAR.push_back(      vetaSTAR);
        dau_phi_STAR.push_back(      vphiSTAR);
        dau_theta_STAR.push_back(    vthetaSTAR);

        dau_vz.push_back(       vdauVZ);
        dau_vy.push_back(       vdauVY);
        dau_vx.push_back(       vdauVX);
        dau_vrefz.push_back(    vdauVrefZ);
        dau_vrefy.push_back(    vdauVrefY);
        dau_vrefx.push_back(    vdauVrefX);
        dau_vp_difZ.push_back(  vp_difZ);
        dau_vp_difY.push_back(  vp_difY);
        dau_vp_difX.push_back(  vp_difX);

        passer = passer +1;
        //if(dau.charge()  ==0 ) continue;
        //const edm::Ref<std::vector<reco::Vertex> > JetVRef = dau.vertexRef().key();
        //int JetVRef_int = dau.vertexRef().key();
        //int fromPV = dau.fromPV(dau.vertexRef().key() );
        //std::cout << j.pt() <<"  "<<  JetVRef_int << "  " << fromPV << std::endl;

        //float dauVRef = recoVertices->at(dau.vertexRef().key()).position().z();
        //float dauV    = dau.vertex().z();
        //std::cout << std::setprecision(4)<< std::endl;
        //std::cout << j.pt() << " " << Ndau_pt_sum << std::endl;

    }
    jetN = passer;


}


void TrackAnalyzer::fillGen(const edm::Event& iEvent){

    edm::Handle< GenEventInfoProduct > genEvt;
    iEvent.getByToken(genEvtInfo_,genEvt);
    genWeight = genEvt->weight();
    genQScale = genEvt->qScale();
    genSignalProcessID = genEvt->signalProcessID();

    edm::Handle<edm::View<pat::PackedGenParticle> > packed;
    iEvent.getByToken(packedGenToken_,packed);

    edm::Handle< std::vector< reco::GenJet > > genJets;
    iEvent.getByToken(packedGenJetToken_,genJets);

    edm::Handle< std::vector< PileupSummaryInfo > > puSummary;
    iEvent.getByToken(puSummary_, puSummary);

    //jets and their  constituents
    for(size_t i=0; i<(*genJets).size();i++){
        const reco::GenJet * jt = &(*genJets)[i];
        if( jt->pt()<minJetPt ) continue;
        if( fabs(jt->eta()) > maxJetEta ) continue;

        genJetPt.push_back(jt->pt());
        genJetEta.push_back(jt->eta());
        genJetPhi.push_back(jt->phi());

        //std::cout << "Jet pt " << jt->pt() << " eta: " << jt->eta() << " phi: " << jt->phi() << std::endl;   
        int chargedMult = 0;
        std::vector< float > tempPt;
        std::vector< float > tempEta;
        std::vector< float > tempPhi;
        std::vector< int > tempPiD;
        std::vector< int > tempChg;
        for( size_t j = 0; j < jt->numberOfDaughters(); j++){
            const reco::Candidate * p = jt->daughter(j);

            tempChg.push_back(p->charge());
            tempPiD.push_back(p->pdgId());

            //check photons for mom to see if direct
            //if(p->pdgId() == 22) std::cout << (p->mother())->pdgId() << std::endl;
            tempPt.push_back(p->pt());
            tempEta.push_back(p->eta());
            tempPhi.push_back(p->phi());

            if( p->charge() == 0 ) continue;

            chargedMult++;

            //std::cout << "PdgID: " << p->pdgId() << " pt " << p->pt() << " eta: " << p->eta() << " phi: " << p->phi() << std::endl;    
        } 
        genJetChargedMultiplicity.push_back( chargedMult );
        gendau_pt.push_back(tempPt);
        gendau_eta.push_back(tempEta);   
        gendau_phi.push_back(tempPhi);   
        gendau_pid.push_back(tempPiD);   
        gendau_chg.push_back(tempChg);   

    }

    //pileup info
    for(unsigned int i = 0; i<(*puSummary).size(); i++){
        int bc = (*puSummary).at(i).getBunchCrossing(); 
        if(bc==0){//ignore out of time pu (shouldn't be in mini AOD though)
            pu = (*puSummary).at(i).getPU_NumInteractions();
            puTrue = (*puSummary).at(i).getTrueNumInteractions();
            puZ = (*puSummary).at(i).getPU_zpositions();
            puPthat = (*puSummary).at(i).getPU_pT_hats();
            puSumPt0p5 = (*puSummary).at(i).getPU_sumpT_highpT();
            puSumPt0p1 = (*puSummary).at(i).getPU_sumpT_lowpT();
            puNTrk0p5 = (*puSummary).at(i).getPU_ntrks_highpT();
            puNTrk0p1 = (*puSummary).at(i).getPU_ntrks_lowpT();
        }
    }
}


//--------------------------------------------------------------------------------------------------
void TrackAnalyzer::fillTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    edm::Handle<edm::View<pat::PackedCandidate>> cands;
    // edm::Handle<edm::ValueMap<float>> chi2Map;

    //loop over packed cands, then loop over lost tracks
    for(int i = 0; i<2; i++){

        if(i==0){
            iEvent.getByToken(packedCandSrc_,cands);
            //   iEvent.getByToken(chi2Map_,chi2Map);
        }
        if(i==1){
            iEvent.getByToken(lostTracksSrc_,cands);
            // iEvent.getByToken(chi2MapLost_,chi2Map);
        }

        for(unsigned it = 0; it<cands->size(); ++it){
            const pat::PackedCandidate & c = (*cands)[it];

            if(!c.hasTrackDetails()) continue;

            reco::Track const& t = c.pseudoTrack();

            if(t.pt() < trackPtMin_) continue;

            trkPt.push_back( t.pt() );
            trkPtError.push_back( t.ptError() );
            trkEta.push_back( t.eta() );
            trkPhi.push_back( t.phi() );
            trkCharge.push_back( (char) t.charge() );
            trkPDGId.push_back( c.pdgId() );
            trkNHits.push_back( (char) t.numberOfValidHits() );
            trkNPixHits.push_back( (char) t.hitPattern().numberOfValidPixelHits() );
            trkNLayers.push_back( (char) t.hitPattern().trackerLayersWithMeasurement() );
            highPurity.push_back( t.quality(reco::TrackBase::qualityByName("highPurity")));
            // trkNormChi2.push_back( (*chi2Map)[cands->ptrAt(it)] );

            //DCA info for associated vtx
            trkAssociatedVtxIndx.push_back( c.vertexRef().key() );
            trkAssociatedVtxQuality.push_back( c.fromPV(c.vertexRef().key() ));
            trkDzAssociatedVtx.push_back( c.dz( c.vertexRef()->position() ) );
            trkDzErrAssociatedVtx.push_back( sqrt( c.dzError()*c.dzError() + c.vertexRef()->zError() * c.vertexRef()->zError() ) );
            trkDxyAssociatedVtx.push_back( c.dxy( c.vertexRef()->position() ) );
            trkDxyErrAssociatedVtx.push_back( sqrt( c.dxyError()*c.dxyError() + c.vertexRef()->xError() * c.vertexRef()->yError() ) );

            //DCA info for first (highest pt) vtx
            if( !xVtx.empty() ){
                math::XYZPoint v(xVtx.at(0),yVtx.at(0), zVtx.at(0));   
                trkFirstVtxQuality.push_back( c.fromPV( 0 ));
                trkDzFirstVtx.push_back( c.dz( v ) );
                trkDzErrFirstVtx.push_back( sqrt( c.dzError()*c.dzError() + zErrVtx.at(0) * zErrVtx.at(0) ) );
                trkDxyFirstVtx.push_back( c.dxy( v ) );
                trkDxyErrFirstVtx.push_back( sqrt( c.dxyError()*c.dxyError() + xErrVtx.at(0) * yErrVtx.at(0) ) );
            }
        }
    }
}


// ------------ method called once each job just before starting event loop  ------------
void TrackAnalyzer::beginJob()
{
    trackTree_ = fs->make<TTree>("trackTree","v1");
    //jetTree_ = fs->make<TTree>("jetTree","v1");

    // event
    trackTree_->Branch("nRun",&nRun,"nRun/l");
    trackTree_->Branch("nEv",&nEv,"nRun/i");

    //trackTree_->Branch("didHLTFire400",&didHLTFire400);
    //trackTree_->Branch("didHLTFire500",&didHLTFire500);
    trackTree_->Branch("didHLTFire400",&MAINpassTrig1);
    trackTree_->Branch("didHLTFire500",&MAINpassTrig2);
    trackTree_->Branch("didHLTFire550",&MAINpassTrig3);
    //trackTree_->Branch("didHLTFire550NOT",&didHLTFire550NOT);

    trackTree_->Branch("nLumi",&nLumi,"nLumi/i");

    // vertex
    trackTree_->Branch("xVtx",&xVtx);
    trackTree_->Branch("yVtx",&yVtx);
    trackTree_->Branch("zVtx",&zVtx);
    trackTree_->Branch("nTracksVtx",&nTracksVtx);
    trackTree_->Branch("ptSumVtx",&ptSumVtx);

    // Tracks
    trackTree_->Branch("trkPt",&trkPt);
    trackTree_->Branch("trkEta",&trkEta);
    trackTree_->Branch("trkPhi",&trkPhi);
    trackTree_->Branch("trkCharge",&trkCharge);
    trackTree_->Branch("trkPDFId",&trkPDGId);
    trackTree_->Branch("trkNHits",&trkNHits);
    trackTree_->Branch("highPurity",&highPurity);
    trackTree_->Branch("trkAssociatedVtxIndx",&trkAssociatedVtxIndx);

    // Jets
    trackTree_->Branch("jetNumDaughters",&jetNumDaughters);
    trackTree_->Branch("jetEta",&jetEta);
    trackTree_->Branch("jetPt",&jetPt);
    trackTree_->Branch("jetPhi",&jetPhi);
    trackTree_->Branch("jetTheta",&jetTheta);
    trackTree_->Branch("jetMass",&jetMass);
    trackTree_->Branch("muonMultiplicity",&muonMultiplicity);
    trackTree_->Branch("chargedMultiplicity",&chargedMultiplicity);
    trackTree_->Branch("jetN",&jetN);


    trackTree_->Branch("dau_pt_sum",      &dau_pt_sum);
    trackTree_->Branch("dau_chg",		&dau_chg); 
    trackTree_->Branch("dau_pid",		&dau_pid);	 
    trackTree_->Branch("dau_vref",	&dau_vref);
    trackTree_->Branch("dau_pt",		&dau_pt);
    trackTree_->Branch("dau_ptError",     &dau_ptError);
    trackTree_->Branch("dau_XYDCAsig",    &dau_XYDCAsig);
    trackTree_->Branch("dau_ZDCAsig",     &dau_ZDCAsig);

    trackTree_->Branch("dau_eta",		&dau_eta);	 
    trackTree_->Branch("dau_phi",		&dau_phi );
    trackTree_->Branch("dau_theta",	&dau_theta);

    trackTree_->Branch("dau_pt_STAR",		&dau_pt_STAR);
    trackTree_->Branch("dau_eta_STAR",		&dau_eta_STAR);	 
    trackTree_->Branch("dau_phi_STAR",		&dau_phi_STAR );
    trackTree_->Branch("dau_theta_STAR",	&dau_theta_STAR);

    trackTree_->Branch("dau_vz",		&dau_vz	 );
    trackTree_->Branch("dau_vy",		&dau_vy	 );
    trackTree_->Branch("dau_vx",		&dau_vx	 );
    trackTree_->Branch("dau_vrefz",	&dau_vrefz);
    trackTree_->Branch("dau_vrefy",	&dau_vrefy);
    trackTree_->Branch("dau_vrefx",	&dau_vrefx);
    trackTree_->Branch("dau_vp_difZ",	&dau_vp_difZ);
    trackTree_->Branch("dau_vp_difY",	&dau_vp_difY);
    trackTree_->Branch("dau_vp_difX",	&dau_vp_difX);

    if(doGen){ 
        trackTree_->Branch("genQScale",&genQScale);
        trackTree_->Branch("genWeight",&genWeight);
        trackTree_->Branch("genSignalProcessID",&genSignalProcessID);

        trackTree_->Branch("genJetEta",&genJetEta);
        trackTree_->Branch("genJetPt",&genJetPt);
        trackTree_->Branch("genJetPhi",&genJetPhi);
        trackTree_->Branch("genJetChargedMultiplicity",&genJetChargedMultiplicity);

        trackTree_->Branch("genDau_chg",		&gendau_chg); 
        trackTree_->Branch("genDau_pid",		&gendau_pid);	 
        trackTree_->Branch("genDau_pt",		&gendau_pt);
        trackTree_->Branch("genDau_eta",		&gendau_eta);	 
        trackTree_->Branch("genDau_phi",		&gendau_phi );

        trackTree_->Branch("nPu",&pu);  
        trackTree_->Branch("nTruePu",&puTrue);
        trackTree_->Branch("puZ",&puZ);
        trackTree_->Branch("puPthat",&puPthat);  
        trackTree_->Branch("puSumPt0p1",&puSumPt0p1);  
        trackTree_->Branch("puSumPt0p5",&puSumPt0p5);  
        trackTree_->Branch("puNTrk0p1",&puNTrk0p1);  
        trackTree_->Branch("puNTrk0p5",&puNTrk0p5);  
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void
TrackAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAnalyzer);
