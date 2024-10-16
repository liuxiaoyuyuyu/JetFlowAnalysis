#define MyClass_cxx
// touch
#include "include/MyTrimPythia.h"
#include "include/coordinateTools.h"

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

const int bin100 = 100;
const int bin200 = 200;
const float bin1 = 1;
const float bin1a = 2;
const float bin1b = 2;
const float etarange = 3.2;
const float pi = 3.14159;
const float bin0= 0;


//---------------------------------------------------------------------CUTS
const float jetEtaCut   = 1.6;
const float jetPtCut    = 550.0;
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

//sort jet phi in descending order 
std::vector<int> F_JetSort(std::vector< float > * vphi, std::vector< float > * vphi2) {
    std::vector<int> result(vphi->size());
    float mR = *std::min_element(vphi->begin(), vphi->end());
    float mG = *std::min_element(vphi2->begin(), vphi2->end());
    float xR = *std::max_element(vphi->begin(), vphi->end());
    float xG = *std::max_element(vphi2->begin(), vphi2->end());
    //Jet phi is in [-pi,pi], if jet phi is around 0 or Pi(-Pi),
    //it is likely that gen jet and reco jet has opposite-sign-similar-fabs phi. 
    //if we don't order jet phi in fabs, reco and gen jet won't get matched. 
    if(fabs(xR)> 3.07 || fabs(mR) < .07 || fabs(xG) > 3.07 || fabs(mG) < .07 ){
        //iota:populates a range with sequentially increasing values
        std::iota(result.begin(), result.end(), 0);
        std::sort(result.begin(), result.end(),
                [&](int A, int B) -> bool {
                return fabs((*vphi)[A]) > fabs((*vphi)[B]);
                });
    }else{
        std::iota(result.begin(), result.end(), 0);
        std::sort(result.begin(), result.end(),
                [&](int A, int B) -> bool {
                return (*vphi)[A] > (*vphi)[B];
                });
    }
    return result;
}


void  keepMatch(std::vector<int> keep, std::vector<int>& indices){
    for(int i = keep.size() - 1; i >= 0; i--){
        if(keep[i] < 2){
            keep.erase(keep.begin() + i);
            indices.erase(indices.begin() + i);
        }
    }
}

//MyClass::Loop(int job, std::string fList){
void MyClass::Loop(int job, std::string fList){

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    //Initializing Histograms

    double ptBins[39] = {
        0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
        1.0, 1.05, 1.1, 1.15, 1.2,  1.3, 1.4,  1.5, 1.6,  1.7, 
        1.8,  1.9, 2.0,  2.5, 3.0,  4.0, 5.0,  7.5, 10.0
    };     
    int len = 38;//sizeof(ptBins)/sizeof(ptBins[0]);

    int bin120=120;int bin500=500;int bin3000=3000;float bin1d6=1.6;float binpi=3.14159;int bin20=20;int bin4=4;int bin10=10;int bin80=80;

    TH1D* h_MAIN_DCA_Dau_Reco_Pt_Lab_All                     = new TH1D("h_MAIN_DCA_Dau_Reco_Pt_Lab_All                 ","h_DCA_Dau_Reco_Pt_Lab_All                 ",len, ptBins);                    
    TH1D* h_MAIN_DCA_Dau_Reco_Phi_Lab_All                    = new TH1D("h_MAIN_DCA_Dau_Reco_Phi_Lab_All                ","h_DCA_Dau_Reco_Phi_Lab_All                ",bin80, -pi,pi);                        
    TH1D* h_MAIN_DCA_Dau_Reco_Eta_Lab_All                    = new TH1D("h_MAIN_DCA_Dau_Reco_Eta_Lab_All                ","h_DCA_Dau_Reco_Eta_Lab_All                ",bin80,-bin4,bin4);                   
    TH2D* h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All                = new TH2D("h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All            ","h2_DCA_Dau_Reco_Pt_Eta_Lab_All            ",len, ptBins,bin80,-bin4,bin4); 
    TH2D* h2_MAIN_DCA_Dau_Reco_Phi_Eta_Lab_All               = new TH2D("h2_MAIN_DCA_Dau_Reco_Phi_Eta_Lab_All           ","h2_DCA_Dau_Reco_Phi_Eta_Lab_All           ",bin80, -pi,pi,bin80,-bin4,bin4);     
    TH2D* h2_MAIN_DCA_Dau_Reco_Phi_Pt_Lab_All                = new TH2D("h2_MAIN_DCA_Dau_Reco_Phi_Pt_Lab_All            ","h2_DCA_Dau_Reco_Phi_Pt_Lab_All            ",bin80, -pi,pi,len, ptBins);     

    TH1D* h_MAIN_ORpt_DCA_Dau_Reco_Pt_Lab_All                     = new TH1D("h_MAIN_ORpt_DCA_Dau_Reco_Pt_Lab_All                 ","h_DCA_ORpt_Dau_Reco_Pt_Lab_All                 ",len, ptBins);                    
    TH1D* h_MAIN_ORpt_DCA_Dau_Reco_Phi_Lab_All                    = new TH1D("h_MAIN_ORpt_DCA_Dau_Reco_Phi_Lab_All                ","h_DCA_ORpt_Dau_Reco_Phi_Lab_All                ",bin80, -pi,pi);                        
    TH1D* h_MAIN_ORpt_DCA_Dau_Reco_Eta_Lab_All                    = new TH1D("h_MAIN_ORpt_DCA_Dau_Reco_Eta_Lab_All                ","h_DCA_ORpt_Dau_Reco_Eta_Lab_All                ",bin80,-bin4,bin4);                   
    TH2D* h2_MAIN_ORpt_DCA_Dau_Reco_Pt_Eta_Lab_All                = new TH2D("h2_MAIN_ORpt_DCA_Dau_Reco_Pt_Eta_Lab_All            ","h2_DCA_ORpt_Dau_Reco_Pt_Eta_Lab_All            ",len, ptBins,bin80,-bin4,bin4); 
    TH2D* h2_MAIN_ORpt_DCA_Dau_Reco_Phi_Eta_Lab_All               = new TH2D("h2_MAIN_ORpt_DCA_Dau_Reco_Phi_Eta_Lab_All           ","h2_DCA_ORpt_Dau_Reco_Phi_Eta_Lab_All           ",bin80, -pi,pi,bin80,-bin4,bin4);     
    TH2D* h2_MAIN_ORpt_DCA_Dau_Reco_Phi_Pt_Lab_All                = new TH2D("h2_MAIN_ORpt_DCA_Dau_Reco_Phi_Pt_Lab_All            ","h2_DCA_ORpt_Dau_Reco_Phi_Pt_Lab_All            ",bin80, -pi,pi,len, ptBins);     

    TH1D* h_Dau_Gen_Pt_Lab_All                     = new TH1D("h_Dau_Gen_Pt_Lab_All                 ","h_Dau_Gen_Pt_Lab_All                 ",len, ptBins);                    
    TH1D* h_Dau_Gen_Phi_Lab_All                    = new TH1D("h_Dau_Gen_Phi_Lab_All                ","h_Dau_Gen_Phi_Lab_All                ",bin80, -pi,pi);                        
    TH1D* h_Dau_Gen_Eta_Lab_All                    = new TH1D("h_Dau_Gen_Eta_Lab_All                ","h_Dau_Gen_Eta_Lab_All                ",bin80,-bin4,bin4);                   
    TH2D* h2_Dau_Gen_Pt_Eta_Lab_All                = new TH2D("h2_Dau_Gen_Pt_Eta_Lab_All            ","h2_Dau_Gen_Pt_Eta_Lab_All            ",len, ptBins,bin80,-bin4,bin4); 
    TH2D* h2_Dau_Gen_Phi_Eta_Lab_All               = new TH2D("h2_Dau_Gen_Phi_Eta_Lab_All           ","h2_Dau_Gen_Phi_Eta_Lab_All           ",bin80, -pi,pi,bin80,-bin4,bin4);     
    TH2D* h2_Dau_Gen_Phi_Pt_Lab_All                = new TH2D("h2_Dau_Gen_Phi_Pt_Lab_All            ","h2_Dau_Gen_Phi_Pt_Lab_All            ",bin80, -pi,pi,len, ptBins);     


    const int nqHats = 8;
    const float qHatBoundaries[nqHats+1] = {470, 600, 800, 1000, 1400, 1800, 2400, 3200, 13000};
    //float xs[nqHats] = {552.1, 156.5, 26.28, 7.47, 0.6484, 0.08743, 0.005236, 0.0001357};
    //Run 3 weights:
    float xs[nqHats] = {628.7, 179.4, 30.8, 8.94, 0.81, 0.1153, 0.007641, 0.0002326};
    TH1D * qHatHist = new TH1D("qhatHit",";;qHat",nqHats,qHatBoundaries);

    for(int ff = 0; ff < fileList.size(); ff++){
        cout << "genQ file: " << ff << endl;
        fFile = TFile::Open(fileList.at(ff).c_str(),"read");

        TTree *tree = (TTree*)fFile->Get("analyzerOffline/trackTree");
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

    std::cout << "Starting event loop" << std::endl;
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;
    for(int f = 0; f<fileList.size(); f++){
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("analyzerOffline/trackTree");
        Init(tree);
        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;

        TFile* jet_veto_file[2];
        jet_veto_file[0]=new TFile("~/StorageArea/Summer22_23Sep2023_RunCD_v1.root","read");
        jet_veto_file[1]=new TFile("~/StorageArea/Summer22EE_23Sep2023_RunEFG_v1.root","read");
        TH2D* jet_veto_map=(TH2D*)jet_veto_file[1]->Get("jetvetomap");

        for (Long64_t ievent=0; ievent <nentries; ievent ++){
            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;

            hpthat->Fill(genQScale, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
            if(genJetPt->size()==0) continue;
            if(genJetChargedMultiplicity->size()==0) continue;
            hleadingJetPt->Fill(genJetPt->at(0), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
            float genWeightPy = xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale));

            if(!F_eventpass(jetPt, jetN)) continue;
            if(!F_eventpass(genJetPt, jetN)) continue;

            //apply jet veto map
            //jet veto, reference: https://cms-jerc.web.cern.ch/Recommendations/#run-3
            //"The safest procedure would be to veto events if ANY jet with a loose selection lies in the veto regions."
            bool jetvetoBool=0;
            for(int ijet=0; ijet < jetPt->size(); ijet++){
                //if( (*jetPt)[ijet] > 15){
                    if(jet_veto_map->GetBinContent(jet_veto_map->FindBin((*jetEta)[ijet],(*jetPhi)[ijet]))>0){
                        jetvetoBool=1;
                        break;
                    } 
                //}

            }

            if(jetvetoBool) continue;

            int rjN = jetPhi->size();
            int gjN = genJetPhi->size();
            std::vector<int> indicesR;
            if(rjN>1){
                indicesR = F_JetSort(jetPhi,genJetPhi);
            }else{
                indicesR.push_back(0);
            }     
            std::vector<int> indicesG;
            if(gjN>1){
                indicesG = F_JetSort(genJetPhi,jetPhi);
            }else{
                indicesG.push_back(0);
            }

            if(indicesG.size() != indicesR.size()){
                // Option One
                // Option One
                if(indicesG.size() > indicesR.size()){
                    //initiate a vector "keep" of size indicesG.size() and fill it with zeros. 
                    std::vector<int> keep(indicesG.size(),0);
                    for(int i = 0; i< indicesG.size(); i++){
                        for(int j = 0; j< indicesR.size(); j++){ 
                            if(fabs((*genJetPhi)[indicesG[i]]) > 3.07 || fabs((*genJetPhi)[indicesG[i]]) < .07 || fabs((*jetPhi)[indicesR[j]]) > 3.07 || fabs((*jetPhi)[indicesR[j]]) < .07 ){
                                if(fabs(fabs((*jetPhi)[indicesR[j]]) - fabs((*genJetPhi)[indicesG[i]])) < 0.1){
                                    keep[i] += 1;
                                }
                            }else{
                                if(fabs((*jetPhi)[indicesR[j]] - (*genJetPhi)[indicesG[i]]) < 0.1){
                                    keep[i] += 1;
                                }
                            }
                            if(fabs((*jetEta)[indicesR[j]] - (*genJetEta)[indicesG[i]]) < 0.1){
                                keep[i] += 1;}
                        }
                    }
                    keepMatch(keep, indicesG);
                    keep.clear();
                }

                // Option Two
                // Option Two
                if(indicesR.size() > indicesG.size()){
                    std::vector<int> keep(indicesR.size(),0);
                    for(int i = 0; i< indicesR.size(); i++){
                        for(int j = 0; j< indicesG.size(); j++){ 
                            if(fabs((*genJetPhi)[indicesG[j]]) > 3.07 || fabs((*genJetPhi)[indicesG[j]]) < .07 || fabs((*jetPhi)[indicesR[i]]) > 3.07 || fabs((*jetPhi)[indicesR[i]]) < .07 ){
                                if(fabs(fabs((*jetPhi)[indicesR[i]]) - fabs((*genJetPhi)[indicesG[j]])) < 0.1){
                                    keep[i] += 1;
                                }
                            }else{
                                if(fabs((*jetPhi)[indicesR[i]] - (*genJetPhi)[indicesG[j]]) < 0.1){
                                    keep[i] += 1;
                                }
                            }
                            if(fabs((*jetEta)[indicesR[i]] - (*genJetEta)[indicesG[j]]) < 0.1){
                                keep[i] += 1;}
                        }
                    }
                    keepMatch(keep, indicesR);
                    keep.clear();
                }
            }
            //End Non-Matching Length Options
            //End Non-Matching Length Options

            std::vector<int> keepRi(indicesR.size());
            std::vector<int> keepGi(indicesG.size());

            for(int i = 0; i< indicesR.size(); i++){
                for(int j = 0; j< indicesG.size(); j++){ 
                    if(fabs((*genJetPhi)[indicesG[j]]) > 3.07 || fabs((*genJetPhi)[indicesG[j]]) < .07 || fabs((*jetPhi)[indicesR[i]]) > 3.07 || fabs((*jetPhi)[indicesR[i]]) < .07 ){
                        if(fabs(fabs((*jetPhi)[indicesR[i]]) - fabs((*genJetPhi)[indicesG[j]])) < 0.1){
                            keepRi[i] += 1;
                            if(fabs((*jetEta)[indicesR[i]] - (*genJetEta)[indicesG[j]]) < 0.1){
                                keepRi[i] += 1;
                            }else{ keepRi[i] -= 1;}
                        }
                    }else{
                        if(fabs((*jetPhi)[indicesR[i]] - (*genJetPhi)[indicesG[j]]) < 0.1){
                            keepRi[i] += 1;
                            if(fabs((*jetEta)[indicesR[i]] - (*genJetEta)[indicesG[j]]) < 0.1){
                                keepRi[i] += 1;
                            }else{ keepRi[i] -= 1;}
                        }
                    }
                }
            }
            for(int i = 0; i< indicesG.size(); i++){
                for(int j = 0; j< indicesR.size(); j++){ 
                    if(fabs((*genJetPhi)[indicesG[i]]) > 3.07 || fabs((*genJetPhi)[indicesG[i]]) < .07 || fabs((*jetPhi)[indicesR[j]]) > 3.07 || fabs((*jetPhi)[indicesR[j]]) < .07 ){
                        if(fabs(fabs((*jetPhi)[indicesR[j]]) - fabs((*genJetPhi)[indicesG[i]])) < 0.1){
                            keepGi[i] += 1;
                            if(fabs((*jetEta)[indicesR[j]] - (*genJetEta)[indicesG[i]]) < 0.1){
                                keepGi[i] += 1;
                            }else{ keepGi[i] -= 1;}
                        }
                    }else{
                        if(fabs((*jetPhi)[indicesR[j]] - (*genJetPhi)[indicesG[i]]) < 0.1){
                            keepGi[i] += 1;
                            if(fabs((*jetEta)[indicesR[j]] - (*genJetEta)[indicesG[i]]) < 0.1){
                                keepGi[i] += 1;
                            }else{ keepGi[i] -= 1;}
                        }
                    }
                }
            }
            keepMatch(keepRi, indicesR);
            keepMatch(keepGi, indicesG);

            int jetCounter1 = jetPt->size();
            int jetCounter2 = genJetPt->size();
            int jetCounter = min(jetCounter1,jetCounter2);
            if(jetCounter == 0) continue;
            if(jetCounter>4) continue;

            for(int ijet=0; ijet < indicesG.size(); ijet++){

                hjetPtReco->Fill(jetPt->at(indicesR[ijet]), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
                hjetPtGen->Fill(genJetPt->at(indicesG[ijet]), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );

                double rphi = (*jetPhi)[indicesR[ijet]];
                double rpt = (*jetPt)[indicesR[ijet]];
                double reta = (*jetEta)[indicesR[ijet]];
                double rchM = (*chargedMultiplicity)[indicesR[ijet]];
                double gphi = (*genJetPhi)[indicesG[ijet]];
                double gpt = (*genJetPt)[indicesG[ijet]];
                double geta = (*genJetEta)[indicesG[ijet]];
                double gchM = (*genJetChargedMultiplicity)[indicesG[ijet]];
                if( rphi == 0 ||  gphi == 0){continue;}
                if( rpt  == 0 ||  gpt  == 0){continue;}
                if( reta == 0 ||  geta == 0){continue;}
                if( rchM == 0 ||  gchM == 0){continue;}


                if( fabs((*jetEta)[indicesR[ijet]]) > jetEtaCut) continue;
                if( fabs((*jetPt)[indicesR[ijet]])  < jetPtCut)  continue;
                //daughter loop:		
                //reco
                long int NNtrk_R = (dau_pt->at( indicesR[ijet] )).size();
                for(int  A_trk=0; A_trk < NNtrk_R; A_trk++ ){
                    int kjet = indicesR[ijet];
                    if(fabs((*dau_eta)[kjet][A_trk]) > 2.4 ) continue;
                    if((*dau_chg)[kjet][A_trk] == 0) continue;

                    double dauptnow = (*dau_pt)[kjet][A_trk];

                    double pterr_R_pt =0;
                    double dcaXY = 0;
                    double dcaZ = 0;
                    pterr_R_pt= ( (*dau_ptError)[kjet][A_trk] ) / ( (*dau_pt)[kjet][A_trk] );
                    dcaXY = (*dau_XYDCAsig)[kjet][A_trk];
                    dcaZ = (*dau_ZDCAsig)[kjet][A_trk];

                    if(pterr_R_pt   > 0.1   && dauptnow > 0.5) continue;
                    if(fabs(dcaXY)  > 3     && dauptnow > 0.5) continue;
                    if(fabs(dcaZ)   > 3     && dauptnow > 0.5) continue;

                    h_MAIN_DCA_Dau_Reco_Pt_Lab_All	->Fill((*dau_pt)[kjet][A_trk],genWeightPy);
                    h_MAIN_DCA_Dau_Reco_Phi_Lab_All	->Fill((*dau_phi)[kjet][A_trk],genWeightPy);
                    h_MAIN_DCA_Dau_Reco_Eta_Lab_All	->Fill((*dau_eta)[kjet][A_trk],genWeightPy);

                    h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All	->Fill((*dau_pt)[kjet][A_trk]	,(*dau_eta)[kjet][A_trk],genWeightPy);
                    h2_MAIN_DCA_Dau_Reco_Phi_Eta_Lab_All->Fill((*dau_phi)[kjet][A_trk]	,(*dau_eta)[kjet][A_trk],genWeightPy);
                    h2_MAIN_DCA_Dau_Reco_Phi_Pt_Lab_All	->Fill((*dau_phi)[kjet][A_trk]	,(*dau_pt)[kjet][A_trk],genWeightPy);

                }
                for(int  A_trk=0; A_trk < NNtrk_R; A_trk++ ){
                    int kjet = indicesR[ijet];
                    if(fabs((*dau_eta)[kjet][A_trk]) > 2.4 ) continue;
                    if((*dau_chg)[kjet][A_trk] == 0) continue;

                    double dauptnow = (*dau_pt)[kjet][A_trk];

                    double pterr_R_pt =0;
                    double dcaXY = 0;
                    double dcaZ = 0;
                    pterr_R_pt= ( (*dau_ptError)[kjet][A_trk] ) / ( (*dau_pt)[kjet][A_trk] );
                    dcaXY = (*dau_XYDCAsig)[kjet][A_trk];
                    dcaZ = (*dau_ZDCAsig)[kjet][A_trk];

                    if(pterr_R_pt   > 0.1   || dauptnow < 0.5) continue;
                    if(fabs(dcaXY)  > 3     || dauptnow < 0.5) continue;
                    if(fabs(dcaZ)   > 3     || dauptnow < 0.5) continue;

                    h_MAIN_ORpt_DCA_Dau_Reco_Pt_Lab_All	->Fill((*dau_pt)[kjet][A_trk],genWeightPy);
                    h_MAIN_ORpt_DCA_Dau_Reco_Phi_Lab_All	->Fill((*dau_phi)[kjet][A_trk],genWeightPy);
                    h_MAIN_ORpt_DCA_Dau_Reco_Eta_Lab_All	->Fill((*dau_eta)[kjet][A_trk],genWeightPy);

                    h2_MAIN_ORpt_DCA_Dau_Reco_Pt_Eta_Lab_All	->Fill((*dau_pt)[kjet][A_trk]	,(*dau_eta)[kjet][A_trk],genWeightPy);
                    h2_MAIN_ORpt_DCA_Dau_Reco_Phi_Eta_Lab_All   ->Fill((*dau_phi)[kjet][A_trk]	,(*dau_eta)[kjet][A_trk],genWeightPy);
                    h2_MAIN_ORpt_DCA_Dau_Reco_Phi_Pt_Lab_All	->Fill((*dau_phi)[kjet][A_trk]	,(*dau_pt)[kjet][A_trk],genWeightPy);


                }
                //gen
                long int NNtrk_G = (genDau_pt->at( indicesG[ijet] )).size();
                for(int  A_trk=0; A_trk < NNtrk_G; A_trk++ ){
                    int kjet = indicesG[ijet];
                    if((*genDau_chg)[kjet][A_trk] == 0) continue;
                    if(fabs((*genDau_eta)[kjet][A_trk]) > 2.4 ) continue;

                    h_Dau_Gen_Pt_Lab_All	->Fill((*genDau_pt)[kjet][A_trk],genWeightPy);
                    h_Dau_Gen_Phi_Lab_All	->Fill((*genDau_phi)[kjet][A_trk],genWeightPy);
                    h_Dau_Gen_Eta_Lab_All	->Fill((*genDau_eta)[kjet][A_trk],genWeightPy);

                    h2_Dau_Gen_Pt_Eta_Lab_All	->Fill((*genDau_pt)[kjet][A_trk]	,(*genDau_eta)[kjet][A_trk],genWeightPy);
                    h2_Dau_Gen_Phi_Eta_Lab_All->Fill((*genDau_phi)[kjet][A_trk]	,(*genDau_eta)[kjet][A_trk],genWeightPy);
                    h2_Dau_Gen_Phi_Pt_Lab_All	->Fill((*genDau_phi)[kjet][A_trk]	,(*genDau_pt)[kjet][A_trk],genWeightPy);

                }
            }
            keepRi.clear();;
            keepGi.clear();
            indicesG.clear();
            indicesR.clear();
        }
        fFile->Close();
    }

    string subList = fList.substr(fList.size() - 3);
    //TFile* fS_tempA = new TFile(Form("out_correction_18/job_%s.root",subList.c_str()), "recreate");
    //TFile* fS_tempA = new TFile(Form("out_correction_17/job_%s.root",subList.c_str()), "recreate");
    //TFile* fS_tempA = new TFile(Form("out_correction_16/job_%s.root",subList.c_str()), "recreate");
    //TFile* fS_tempA = new TFile(Form("correct_by_era/out_correction_16APV/job_%s.root",subList.c_str()), "recreate");
    TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/xiaoyul/Run3_2022_root_out/correction_EFG/job_%s.root",subList.c_str()), "recreate");
   
    hpthat->Write(); 
    hleadingJetPt->Write(); 
    hjetPtReco->Write();
    hjetPtGen->Write();

    h_MAIN_DCA_Dau_Reco_Pt_Lab_All                   ->Write(); 
    h_MAIN_DCA_Dau_Reco_Phi_Lab_All                  ->Write(); 
    h_MAIN_DCA_Dau_Reco_Eta_Lab_All                  ->Write(); 
    h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All              ->Write(); 
    h2_MAIN_DCA_Dau_Reco_Phi_Eta_Lab_All             ->Write(); 
    h2_MAIN_DCA_Dau_Reco_Phi_Pt_Lab_All              ->Write(); 
    h_MAIN_ORpt_DCA_Dau_Reco_Pt_Lab_All                   ->Write(); 
    h_MAIN_ORpt_DCA_Dau_Reco_Phi_Lab_All                  ->Write(); 
    h_MAIN_ORpt_DCA_Dau_Reco_Eta_Lab_All                  ->Write(); 
    h2_MAIN_ORpt_DCA_Dau_Reco_Pt_Eta_Lab_All              ->Write(); 
    h2_MAIN_ORpt_DCA_Dau_Reco_Phi_Eta_Lab_All             ->Write(); 
    h2_MAIN_ORpt_DCA_Dau_Reco_Phi_Pt_Lab_All              ->Write(); 

    h_Dau_Gen_Pt_Lab_All                    ->Write(); 
    h_Dau_Gen_Phi_Lab_All                   ->Write(); 
    h_Dau_Gen_Eta_Lab_All                   ->Write(); 
    h2_Dau_Gen_Pt_Eta_Lab_All               ->Write(); 
    h2_Dau_Gen_Phi_Eta_Lab_All              ->Write(); 
    h2_Dau_Gen_Phi_Pt_Lab_All               ->Write(); 

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

