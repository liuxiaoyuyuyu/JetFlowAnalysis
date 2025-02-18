#define MyClass_cxx
// touch

//#define _doGEN
#define _doRECO

#include "include/MyTrimPythia.h"
#include "include/coordinateTools.h"
#include "include/mc_constants.h"

#include <iostream>
#include <iomanip>

#include <vector>
#include "math.h"
#include <numeric>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>

#include <TStyle.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TRandom3.h"

using TMath::ATan;
using TMath::Exp;


//---------------------------------------------------------------------CUTS
const double ptscale = 1000;
const double etascale = 3.2;
const double phiscale = 2*3.14159;
const double multscale = 100;


bool F_eventpass(std::vector< float > *jetPt, int jetnumber){

    if(jetPt->size() < 1)           return false;
    if(jetPt->size() > 4)           return false;
    if(jetnumber > 200)             return false;
    //if((*jetPt)[0] < jetPtCut) return false;

    return true;
}

bool isSubstring(string s1, string s2)
{
    int M = s1.length();
    int N = s2.length();

    /* A loop to slide pat[] one by one */
    for (int i = 0; i <= N - M; i++) {
        int j;

        /* For current index i, check for
         *  pattern match */
        for (j = 0; j < M; j++)
            if (s2[i + j] != s1[j])
                break;

        if (j == M)
            return 1;
    }

    return 0;
}


//MyClass::Loop(int job, std::string fList){
void MyClass::Loop(int job, std::string fList){

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    //Initializing Histograms

    const int nqHats = 8;
    const float qHatBoundaries[nqHats+1] = {470, 600, 800, 1000, 1400, 1800, 2400, 3200, 13000};
    float xs[nqHats] = {552.1, 156.5, 26.28, 7.47, 0.6484, 0.08743, 0.005236, 0.0001357};
    TH1D * qHatHist = new TH1D("qhatHit",";;qHat",nqHats,qHatBoundaries);

    for(int ff = 0; ff < fileList.size(); ff++){
        cout << "genQ file: " << ff << endl;
        fFile = TFile::Open(fileList.at(ff).c_str(),"read");

        TTree *tree = (TTree*)fFile->Get("analyzer/trackTree");
        if(!fFile || fFile->IsZombie()){
                std::cout << "File " << ff+1 << " out of " << fileList.size() <<" is a Zombie, skipped."<< std::endl;
                continue;
        } 
        Init(tree);

        for( int i = 0; i < tree->GetEntries(); i++){
            tree->GetEntry(i);
            qHatHist->Fill(genQScale);
        }

        fFile->Close();
    }

    const int i40 = 40;
    const int i470 = 470;
    const int i3500 = 3500;
    const int i50 = 50;
    const int i500 = 500;
    const int i1500 = 1500;
    const int i150 = 150;
    const int i100 = 100;

    TH1D * hpthat = new TH1D("pthat",";#hat{q};#sigma (pb)",i40,i470,i3500);
    TH1D * hleadingJetPt = new TH1D("leadingJetPt",";Leading p_{T}^{gen};#sigma (pb)",i50,i500,i1500);
    TH1D * hjetPtReco = new TH1D("JetPtReco",";p_{T}^{gen};#sigma (pb)",i150,i100,i3500);
    TH1D * hjetPtGen = new  TH1D("JetPtGen",";p_{T}^{gen};#sigma (pb)",i150,i100,i3500);

    TH1D* hJet_Pass550     = new TH1D("hJet_Pass550"  ,"hJet_Pass550"  , trackbin,bin0,trackbin);
    
    TH1D* hBinDist_cor_single = new TH1D("hBinDist_cor_single","hBinDist_cor_single",bin360,bin0,bin120);
    TH1D* hBinDist_unc_single = new TH1D("hBinDist_unc_single","hBinDist_unc_single",bin120,bin0,bin120);

    TH2D* h_lab_JetMult_pT=new TH2D("lab_JetMult_pT","lab_JetMult_pT",340,100,3500,120,0.5,120.5);
    TH2D* h_lab_JetMult_phi=new TH2D("lab_JetMult_phi","lab_JetMult_phi",30,-TMath::Pi(),TMath::Pi(),120,0.5,120.5);
    TH2D* h_lab_JetMult_eta=new TH2D("lab_JetMult_eta","lab_JetMult_eta",34,-1.7,1.7,120,0.5,120.5);

    TH2D* h_lab_cor_JetMult_pT=new TH2D("lab_cor_JetMult_pT","lab_cor_JetMult_pT",340,100,3500,360,0.5,120.5);
    TH2D* h_lab_cor_JetMult_phi=new TH2D("lab_cor_JetMult_phi","lab_cor_JetMult_phi",30,-TMath::Pi(),TMath::Pi(),360,0.5,120.5);
    TH2D* h_lab_cor_JetMult_eta=new TH2D("lab_cor_JetMult_eta","lab_cor_JetMult_eta",34,-1.7,1.7,360,0.5,120.5);

    TH1D* h_lab_track_eta =new TH1D("lab_track_eta","lab_track_eta", 50, -2.5,2.5);
    TH1D* h_lab_track_pt  =new TH1D("lab_track_pt","lab_track_pt",600,0,600);
    TH1D* h_lab_track_phi =new TH1D("lab_track_phi","lab_track_phi",30,-TMath::Pi(),TMath::Pi());

    TH1D* h_lab_cor_track_eta =new TH1D("lab_cor_track_eta","lab_cor_track_eta", 50, -2.5,2.5);
    TH1D* h_lab_cor_track_pt  =new TH1D("lab_cor_track_pt","lab_cor_track_pt",600,0,600);
    TH1D* h_lab_cor_track_phi =new TH1D("lab_cor_track_phi","lab_cor_track_phi",30,-TMath::Pi(),TMath::Pi());
    
    TH1D* h_jet_jT[trackbin];
    TH1D* h_jet_etastar[trackbin];
    TH1D* h_jet_cor_jT[trackbin];
    TH1D* h_jet_cor_etastar[trackbin];

    TH2D* h_jet_jT_etastar[trackbin];
    TH2D* h_jet_cor_jT_etastar[trackbin];
    
    TH1D* hBinDist_cor[trackbin];
    TH1D* hBinDist_unc[trackbin];
    
    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_cor[wtrk-1]    = new TH1D(Form("hBinDist_cor_%d",wtrk),Form("hBinDist_cor_%d",wtrk), bin360, bin0, bin120);
        hBinDist_unc[wtrk-1]    = new TH1D(Form("hBinDist_unc_%d",wtrk),Form("hBinDist_unc_%d",wtrk), bin360, bin0, bin120);
        
        h_jet_jT[wtrk-1]=new TH1D(Form("jet_jT_%d",wtrk),Form("jet_jT_%d",wtrk),200,0,20);
        h_jet_etastar[wtrk-1]=new TH1D(Form("jet_etastar_%d",wtrk),Form("jet_etastar_%d",wtrk),100,0,10);
        h_jet_cor_jT[wtrk-1]=new TH1D(Form("jet_cor_jT_%d",wtrk),Form("jet_cor_jT_%d",wtrk),200,0,20);
        h_jet_cor_etastar[wtrk-1]=new TH1D(Form("jet_cor_etastar_%d",wtrk),Form("jet_cor_etastar_%d",wtrk),100,0,10);
        
        h_jet_jT_etastar[wtrk-1]=new TH2D(Form("jet_jT_etastar_%d",wtrk),Form("jet_jT_etastar_%d",wtrk),100,0,10,200,0,20);
        h_jet_cor_jT_etastar[wtrk-1]=new TH2D(Form("jet_cor_jT_etastar_%d",wtrk),Form("jet_cor_jT_etastar_%d",wtrk),100,0,10,200,0,20);
    }


    TH2D* hReco2D[fileList.size()];
    TH2D* hGen2D[fileList.size()];

    std::cout << "Starting event loop" << std::endl;
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;
    for(int f = 0; f<fileList.size(); f++){
        int f_from_file = f;
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("analyzer/trackTree");
        if(!fFile || fFile->IsZombie()){
                std::cout << "File " << f+1 << " out of " << fileList.size() <<" is a Zombie, skipped."<< std::endl;
                continue;
        } 
        Init(tree);
        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;

        vector<string> era_vec;
        vector<string> matched_cor_table_vec;
        era_vec.push_back("Run2018_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("Run2017_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("Run2016_"); matched_cor_table_vec.push_back("16_v2");
        era_vec.push_back("Run2016APV_"); matched_cor_table_vec.push_back("16apv_v2");
        int i_keep=999;
        for(int i=0; i<era_vec.size(); i++){
            if(isSubstring(era_vec[i],fileList.at(f).c_str())){
                i_keep=i;
                continue;
            }
        }

        if(i_keep==999) continue;

        string filename = "~/StorageArea/"+matched_cor_table_vec[i_keep]+".root";
        TFile *f_pt_eta_DCA_lookup = new TFile(filename.c_str());
        
        hReco2D[f] = (TH2D*)f_pt_eta_DCA_lookup->Get("h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All")->Clone(Form("h2_DCA_Dau_Reco_Pt_Eta_Lab_All_%d",f));
        hGen2D[f] = (TH2D*)f_pt_eta_DCA_lookup->Get("h2_Dau_Gen_Pt_Eta_Lab_All")->Clone(Form("h2_Dau_Gen_Pt_Eta_Lab_All_%d",f));
        hReco2D[f]->Divide(hGen2D[f]);
        int thisEffTable =f_from_file;
        
        
        //ENTERING EVENT LOOP 
        for (Long64_t ievent=0; ievent <nentries; ievent ++){
            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;

            //float genWeightPy = xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale));
            float genWeightPy =1.0; 
            hpthat->Fill(genQScale, genWeightPy );
            
            #ifdef _doGEN
            if(genJetPt->size()==0) continue;
            if(genJetChargedMultiplicity->size()==0) continue;
            hleadingJetPt->Fill(genJetPt->at(0), genWeightPy);

            //a loose cut requiring the size of jetPt and genJetPt no greater than 4.
            if(!F_eventpass(genJetPt, jetN)) continue;

            //==============ENTERING GEN JET LOOP==================
            for(int ijet=0; ijet < genJetPt->size(); ijet++){
                long int NNtrk = (genDau_pt->at(ijet)).size();
                if( fabs(((*genJetEta)[ijet])) > jetEtaCut ) continue;
                if( (*genJetPt)[ijet] < jetPtCut_Jet   ) continue;
                hjetPtGen->Fill((*genJetPt)[ijet]);

                double jet_HLT_weight = 1.0;
                int tkBool[trackbin] = {0};

                double n_ChargeMult_DCA_labPt_Eta_exclusion_Cor =0;

                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;//charge
                    if(fabs((*genDau_pt)[ijet][A_trk])  < 0.3)     continue;//lab pt
                    if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4)     continue;//lab eta

                    double nUnc_weight = 1.0;//(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    n_ChargeMult_DCA_labPt_Eta_exclusion_Cor += (1.0/nUnc_weight);
                }

                h_lab_cor_JetMult_pT->Fill((*genJetPt)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion_Cor);
                h_lab_cor_JetMult_phi->Fill((*genJetPhi)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion_Cor);
                h_lab_cor_JetMult_eta->Fill((*genJetEta)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion_Cor);
                
                hBinDist_cor_single            ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor, 1.0*jet_HLT_weight);
                for(int i = 0; i < trackbin; i++){
                    if( n_ChargeMult_DCA_labPt_Eta_exclusion_Cor >= trackbinbounds[i] && n_ChargeMult_DCA_labPt_Eta_exclusion_Cor < trackbinboundsUpper[i]){
                        tkBool[i] = 1;
                        hBinDist_cor[i]             ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor,1.0*jet_HLT_weight);
                    }
                }

                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;
                    if((*genDau_pt)[ijet][A_trk] < 0.3) continue;
                    if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4) continue;

                    double jet_dau_pt = ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    double jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet]   , (double)(*genJetPhi)[ijet]  , (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet]  , (double)(*genJetPhi)[ijet]   , (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);

                    //getting et pt efficiency for A track in beam frame
                    //double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    double Atrk_weight=1.0; 

                    h_lab_cor_track_eta->Fill((*genDau_eta)[ijet][A_trk]); 
                    h_lab_cor_track_pt->Fill((*genDau_pt)[ijet][A_trk]); 
                    h_lab_cor_track_phi->Fill((*genDau_phi)[ijet][A_trk]); 

                    //in the case where the total jet mult and the individual daughter pt is acceptable for this track bin and pt bin, we increase the Ntrig count.
                    for(int i = 0; i < trackbin; i++){
                        if(tkBool[i]==1){
                            h_jet_cor_jT[i]->Fill(jet_dau_pt,1.0/(Atrk_weight));
                            h_jet_cor_etastar[i]->Fill(jet_dau_eta,1.0/(Atrk_weight));
                            h_jet_cor_jT_etastar[i]->Fill(jet_dau_eta,jet_dau_pt,1.0*jet_HLT_weight/(Atrk_weight));
                        } 
                    }

                }//here ends the boolean array creations.
            }
            #endif

            #ifdef _doRECO
            if(jetPt->size()==0) continue;
            if(chargedMultiplicity->size()==0) continue;
            hleadingJetPt->Fill(jetPt->at(0), genWeightPy);

            //a loose cut requiring the size of jetPt and genJetPt no greater than 4.
            if(!F_eventpass(jetPt, jetN)) continue;

            //==============ENTERING GEN JET LOOP==================
            for(int ijet=0; ijet < jetPt->size(); ijet++){
                long int NNtrk = (dau_pt->at(ijet)).size();
                if( fabs(((*jetEta)[ijet])) > jetEtaCut ) continue;
                if( (*jetPt)[ijet] < jetPtCut_Jet   ) continue;
                hjetPtReco->Fill((*jetPt)[ijet]);

                double jet_HLT_weight = 1.0;
                int tkBool[trackbin] = {0};

                int n_ChargeMult_DCA_labPt_Eta_exclusion =0;
                double n_ChargeMult_DCA_labPt_Eta_exclusion_Cor =0;

                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    if((*dau_chg)[ijet][A_trk] == 0) continue;//charge
                    if(fabs((*dau_pt)[ijet][A_trk])  < 0.3)     continue;//lab pt
                    if(fabs((*dau_eta)[ijet][A_trk]) > 2.4)     continue;//lab eta

                    double dauptnow = (*dau_pt)[ijet][A_trk];

                    double pterr_R_pt =0;
                    double dcaXY = 0;
                    double dcaZ = 0;
                    pterr_R_pt= ( (*dau_ptError)[ijet][A_trk] ) / ( (*dau_pt)[ijet][A_trk] );
                    dcaXY = (*dau_XYDCAsig)[ijet][A_trk];
                    dcaZ = (*dau_ZDCAsig)[ijet][A_trk];

                    if(pterr_R_pt   > 0.1   && dauptnow > 0.5) continue;
                    if(fabs(dcaXY)  > 3     && dauptnow > 0.5) continue;
                    if(fabs(dcaZ)   > 3     && dauptnow > 0.5) continue;
                    n_ChargeMult_DCA_labPt_Eta_exclusion += 1;
                    double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    n_ChargeMult_DCA_labPt_Eta_exclusion_Cor += (1.0/nUnc_weight);
                }

                h_lab_JetMult_pT->Fill((*jetPt)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion);
                h_lab_JetMult_phi->Fill((*jetPhi)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion);
                h_lab_JetMult_eta->Fill((*jetEta)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion);

                h_lab_cor_JetMult_pT->Fill((*jetPt)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion_Cor);
                h_lab_cor_JetMult_phi->Fill((*jetPhi)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion_Cor);
                h_lab_cor_JetMult_eta->Fill((*jetEta)[ijet],n_ChargeMult_DCA_labPt_Eta_exclusion_Cor);
                
                hBinDist_cor_single            ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor, 1.0*jet_HLT_weight);
                hBinDist_unc_single            ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion, 1.0*jet_HLT_weight);
                for(int i = 0; i < trackbin; i++){
                    if( n_ChargeMult_DCA_labPt_Eta_exclusion >= trackbinbounds[i] && n_ChargeMult_DCA_labPt_Eta_exclusion < trackbinboundsUpper[i]){
                        tkBool[i] = 1;
                        hBinDist_unc[i]             ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion,1.0*jet_HLT_weight);
                        hBinDist_cor[i]             ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor,1.0*jet_HLT_weight);
                    }
                }

                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    if((*dau_chg)[ijet][A_trk] == 0) continue;
                    if((*dau_pt)[ijet][A_trk] < 0.3) continue;
                    if(fabs((*dau_eta)[ijet][A_trk]) > 2.4) continue;

                    double dauptnow = (*dau_pt)[ijet][A_trk];
                    double pterr_R_pt =0;
                    double dcaXY = 0;
                    double dcaZ = 0;
                    pterr_R_pt= ( (*dau_ptError)[ijet][A_trk] ) / ( (*dau_pt)[ijet][A_trk] );
                    dcaXY = (*dau_XYDCAsig)[ijet][A_trk];
                    dcaZ = (*dau_ZDCAsig)[ijet][A_trk];
                    if(pterr_R_pt   > 0.1   && dauptnow > 0.5) continue;
                    if(fabs(dcaXY)  > 3     && dauptnow > 0.5) continue;
                    if(fabs(dcaZ)   > 3     && dauptnow > 0.5) continue;

                    double jet_dau_pt = ptWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet], (double)(*jetPhi)[ijet], (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_eta   = etaWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet]     , (double)(*jetPhi)[ijet]  , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet]    , (double)(*jetPhi)[ijet]  , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);

                    //getting et pt efficiency for A track in beam frame
                    double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    //double Atrk_weight=1.0; 

                    h_lab_track_eta->Fill((*dau_eta)[ijet][A_trk]); 
                    h_lab_track_pt->Fill((*dau_pt)[ijet][A_trk]); 
                    h_lab_track_phi->Fill((*dau_phi)[ijet][A_trk]); 
                    
                    h_lab_cor_track_eta->Fill((*dau_eta)[ijet][A_trk],1.0/(Atrk_weight)); 
                    h_lab_cor_track_pt->Fill((*dau_pt)[ijet][A_trk],1.0/(Atrk_weight)); 
                    h_lab_cor_track_phi->Fill((*dau_phi)[ijet][A_trk],1.0/(Atrk_weight)); 

                    //in the case where the total jet mult and the individual daughter pt is acceptable for this track bin and pt bin, we increase the Ntrig count.
                    for(int i = 0; i < trackbin; i++){
                        if(tkBool[i]==1){
                            h_jet_jT[i]->Fill(jet_dau_pt,1.0*jet_HLT_weight);
                            h_jet_etastar[i]->Fill(jet_dau_eta,1.0*jet_HLT_weight);
                            h_jet_jT_etastar[i]->Fill(jet_dau_eta,jet_dau_pt,1.0*jet_HLT_weight);
                            
                            h_jet_cor_jT[i]->Fill(jet_dau_pt,1.0*jet_HLT_weight/(Atrk_weight));
                            h_jet_cor_etastar[i]->Fill(jet_dau_eta,1.0*jet_HLT_weight/(Atrk_weight));
                            h_jet_cor_jT_etastar[i]->Fill(jet_dau_eta,jet_dau_pt,1.0*jet_HLT_weight/(Atrk_weight));
                        } 
                    }

                }//here ends the boolean array creations.
            }
            
            #endif
        }
        fFile->Close();
    }

    string subList = fList.substr(fList.size() - 14);
    #ifdef _doGEN
    TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/xiaoyul/MC/13TeV/closure/nophatWeight/nojetmatch/gen/job_%s.root",subList.c_str()), "recreate");
    #endif

    #ifdef _doRECO
    TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/xiaoyul/MC/13TeV/closure/nophatWeight/nojetmatch/reco/job_%s.root",subList.c_str()), "recreate");
    #endif

    qHatHist->Write();
    hpthat->Write();
    hleadingJetPt->Write(); 
    hjetPtReco->Write();
    hjetPtGen->Write();

    hJet_Pass550->Write();
    hBinDist_cor_single   ->Write();
    hBinDist_unc_single   ->Write(); 
    
    h_lab_JetMult_pT      ->Write();
    h_lab_JetMult_phi     ->Write();
    h_lab_JetMult_eta     ->Write();

    h_lab_cor_JetMult_pT  ->Write();  
    h_lab_cor_JetMult_phi ->Write();
    h_lab_cor_JetMult_eta ->Write();

    h_lab_track_eta       ->Write();
    h_lab_track_pt        ->Write();
    h_lab_track_phi       ->Write();

    h_lab_cor_track_eta   ->Write();
    h_lab_cor_track_pt    ->Write();
    h_lab_cor_track_phi   ->Write();


    for(int i=0;i<trackbin;i++){
        h_jet_jT[i]->Write();
        h_jet_etastar[i]->Write();
        h_jet_cor_jT[i]->Write();
        h_jet_cor_etastar[i]->Write();

        h_jet_jT_etastar[i]->Write();
        h_jet_cor_jT_etastar[i]->Write();

        hBinDist_cor[i]->Write();
        hBinDist_unc[i]->Write();
    }
    
    fS_tempA->Close();
}
//Code enters execution here
int main(int argc, const char* argv[])
{
    if(argc != 4)
    {
        std::cout << "Usage: Z_mumu_Channel <fileList> <jobNumber> <nJobs>" << std::endl;
        return 1;
    }


    //read input parameters
    std::string fList = argv[1];
    std::string buffer;
    std::vector<std::string> listOfFiles;
    std::ifstream inFile(fList.data());

    int job = (int)std::atoi(argv[2]);
    int nJobs = (int)std::atoi(argv[3]);


    //read the file list and spit it into a vector of strings based on how the parallelization is to be done
    //each vector is a separate subset of the fileList based on the job number
    if(!inFile.is_open())
    {
        std::cout << "Error opening jet file. Exiting." <<std::endl;
        return 1;
    }
    else
    {
        int line = 0;
        while(true)
        {
            inFile >> buffer;
            if(inFile.eof()) break;
            if( line%nJobs == job) listOfFiles.push_back(buffer);
            line++;
        }
    }

    //create the MyClass Object
    MyClass m = MyClass(listOfFiles);
    m.Loop(job, fList);

    return 0;
}

