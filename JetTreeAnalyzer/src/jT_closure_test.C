#define MyClass_cxx
// touch
#include "include/MyTrimPythia.h"
#include "include/coordinateTools.h"
#include "include/1d2d_constants.h"

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

std::vector<int> F_JetSort(std::vector< float > * vphi, std::vector< float > * vphi2) {
    std::vector<int> result(vphi->size());
    float mR = *std::min_element(vphi->begin(), vphi->end());
    float mG = *std::min_element(vphi2->begin(), vphi2->end());
    float xR = *std::max_element(vphi->begin(), vphi->end());
    float xG = *std::max_element(vphi2->begin(), vphi2->end());
    if(fabs(xR)> 3.07 || fabs(mR) < .07 || fabs(xG) > 3.07 || fabs(mG) < .07 ){
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

    const int nqHats = 8;
    const float qHatBoundaries[nqHats+1] = {470, 600, 800, 1000, 1400, 1800, 2400, 3200, 13000};
    float xs[nqHats] = {552.1, 156.5, 26.28, 7.47, 0.6484, 0.08743, 0.005236, 0.0001357};
    TH1D * qHatHist = new TH1D("qhatHit",";;qHat",nqHats,qHatBoundaries);

    for(int ff = 0; ff < fileList.size(); ff++){
        cout << "genQ file: " << ff << endl;
        fFile = TFile::Open(fileList.at(ff).c_str(),"read");

        TTree *tree = (TTree*)fFile->Get("analyzer/trackTree");
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

    TH1D* reco_hBinDist_cor_single = new TH1D("reco_hBinDist_single","reco_hBinDist_single",bin360,bin0,bin120);
    TH1D* reco_hBinDist_unc_single = new TH1D("reco_hBinDist_unc_single","reco_hBinDist_unc_single",bin120,bin0,bin120);
    TH1D* reco_h_jet_unc_jT[trackbin];
    TH1D* reco_h_jet_unc_etastar[trackbin];
    TH1D* reco_h_jet_cor_jT[trackbin];
    TH1D* reco_h_jet_cor_etastar[trackbin];
    TH1D* reco_hBinDist_cor[trackbin];
    TH1D* reco_hBinDist_unc[trackbin];
    
    TH1D* gen_hBinDist_unc_single = new TH1D("gen_hBinDist_single","gen_hBinDist_single",bin120,bin0,bin120);
    TH1D* gen_h_jet_unc_jT[trackbin];
    TH1D* gen_h_jet_unc_etastar[trackbin];
    TH1D* gen_hBinDist_unc[trackbin];

    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        reco_hBinDist_unc[wtrk-1]    = new TH1D(Form("reco_hBinDist_unc_%d",wtrk),Form("reco_hBinDist_unc_%d",wtrk), bin360, bin0, bin120);
        reco_h_jet_unc_jT[wtrk-1]=new TH1D(Form("reco_jet_jT_unc_%d",wtrk),Form("reco_jet_jT_unc_%d",wtrk),200,0,20);
        reco_h_jet_unc_etastar[wtrk-1]=new TH1D(Form("reco_jet_etastar_unc_%d",wtrk),Form("reco_jet_etastar_unc_%d",wtrk),100,0,10);
        reco_hBinDist_cor[wtrk-1]    = new TH1D(Form("reco_hBinDist_%d",wtrk),Form("reco_hBinDist_%d",wtrk), bin360, bin0, bin120);
        reco_h_jet_cor_jT[wtrk-1]=new TH1D(Form("reco_jet_jT_%d",wtrk),Form("reco_jet_jT_%d",wtrk),200,0,20);
        reco_h_jet_cor_etastar[wtrk-1]=new TH1D(Form("reco_jet_etastar_%d",wtrk),Form("reco_jet_etastar_%d",wtrk),100,0,10);
        gen_hBinDist_unc[wtrk-1]    = new TH1D(Form("gen_hBinDist_%d",wtrk),Form("gen_hBinDist_%d",wtrk), bin360, bin0, bin120);
        gen_h_jet_unc_jT[wtrk-1]=new TH1D(Form("gen_jet_jT_%d",wtrk),Form("gen_jet_jT_%d",wtrk),200,0,20);
        gen_h_jet_unc_etastar[wtrk-1]=new TH1D(Form("gen_jet_etastar_%d",wtrk),Form("gen_jet_etastar_%d",wtrk),100,0,10);
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

            hpthat->Fill(genQScale, xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
            if(genJetPt->size()==0) continue;
            if(genJetChargedMultiplicity->size()==0) continue;
            hleadingJetPt->Fill(genJetPt->at(0), xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale)) );
            float genWeightPy = xs[qHatHist->FindBin(genQScale) - 1 ] / qHatHist->GetBinContent(qHatHist->FindBin(genQScale));

            //a loose cut requiring the size of jetPt and genJetPt no greater than 4.
            if(!F_eventpass(jetPt, jetN)) continue;
            if(!F_eventpass(genJetPt, jetN)) continue;

            //**************************MATCH JETS***********************
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
                if(indicesG.size() > indicesR.size()){
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
            //**************************END MATCH JETS***********************
            
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
                if( fabs((*jetPt)[indicesR[ijet]])  < jetPtCut_Jet)  continue;
        //daughter loop:		
            //gen: calculate Nch
                int gen_n_ChargeMult_DCA_labPt_Eta_exclusion =0;
                
                long int NNtrk_G = (genDau_pt->at( indicesG[ijet] )).size();
                for(int  A_trk=0; A_trk < NNtrk_G; A_trk++ ){
                    int kjet = indicesG[ijet];
                    if((*genDau_chg)[kjet][A_trk] == 0) continue;
                    if(fabs((*genDau_pt)[kjet][A_trk])  < 0.3)     continue;//lab pt
                    if(fabs((*genDau_eta)[kjet][A_trk]) > 2.4 ) continue;

                    gen_n_ChargeMult_DCA_labPt_Eta_exclusion += 1;
                }

                gen_hBinDist_unc_single->Fill(gen_n_ChargeMult_DCA_labPt_Eta_exclusion,1.0*genWeightPy);          
            
            //reco: calculate Nch and corrected Nch.
                int tkBool[trackbin] = {0};
                int n_ChargeMult_DCA_labPt_Eta_exclusion =0;
                double n_ChargeMult_DCA_labPt_Eta_exclusion_Cor =0;
                
                long int NNtrk_R = (dau_pt->at( indicesR[ijet] )).size();
                for(int  A_trk=0; A_trk < NNtrk_R; A_trk++ ){
                    int kjet = indicesR[ijet];
                    if((*dau_chg)[kjet][A_trk] == 0) continue;
                    if(fabs((*dau_pt)[kjet][A_trk])  < 0.3)     continue;//lab pt
                    if(fabs((*dau_eta)[kjet][A_trk]) > 2.4 ) continue;

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
                    n_ChargeMult_DCA_labPt_Eta_exclusion += 1;
                    double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[kjet][A_trk] , (*dau_eta)[kjet][A_trk] )));
                    n_ChargeMult_DCA_labPt_Eta_exclusion_Cor += (1.0/nUnc_weight);
                }

                for(int i = 0; i < trackbin; i++){
                    if( n_ChargeMult_DCA_labPt_Eta_exclusion >= trackbinbounds[i] && n_ChargeMult_DCA_labPt_Eta_exclusion < trackbinboundsUpper[i]){
                        tkBool[i] = 1;
                        reco_hBinDist_cor[i]->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor,1.0*genWeightPy);
                        reco_hBinDist_unc[i]->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion,1.0*genWeightPy);
                        gen_hBinDist_unc[i]->Fill(gen_n_ChargeMult_DCA_labPt_Eta_exclusion,1.0*genWeightPy);
                    }
                }

                reco_hBinDist_cor_single->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor,1.0*genWeightPy);          
                reco_hBinDist_unc_single->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion,1.0*genWeightPy);          

            //gen:fill jT dis.

                for(int  A_trk=0; A_trk < NNtrk_G; A_trk++ ){
                    int kjet = indicesG[ijet];
                    if((*genDau_chg)[kjet][A_trk] == 0) continue;
                    if(fabs((*genDau_pt)[kjet][A_trk])  < 0.3)     continue;//lab pt
                    if(fabs((*genDau_eta)[kjet][A_trk]) > 2.4 ) continue;

                    double jet_dau_pt = ptWRTJet((double)(*genJetPt)[kjet], (double)(*genJetEta)[kjet], (double)(*genJetPhi)[kjet], (double)(*genDau_pt)[kjet][A_trk], (double)(*genDau_eta)[kjet][A_trk], (double)(*genDau_phi)[kjet][A_trk]);
                    double jet_dau_eta   = etaWRTJet((double)(*genJetPt)[kjet], (double)(*genJetEta)[kjet]  , (double)(*genJetPhi)[kjet]  , (double)(*genDau_pt)[kjet][A_trk], (double)(*genDau_eta)[kjet][A_trk], (double)(*genDau_phi)[kjet][A_trk]);
                
                    for(int i = 0; i < trackbin; i++){
                        if(tkBool[i]==1){
                            gen_h_jet_unc_jT[i]->Fill(jet_dau_pt,1.0*genWeightPy);
                            gen_h_jet_unc_etastar[i]->Fill(jet_dau_eta,1.0*genWeightPy);
                        } 
                    }
                }

            //reco: fill jT dis.
                for(int  A_trk=0; A_trk < NNtrk_R; A_trk++ ){
                    int kjet = indicesR[ijet];
                    if((*dau_chg)[kjet][A_trk] == 0) continue;
                    if(fabs((*dau_pt)[kjet][A_trk])  < 0.3)     continue;//lab pt
                    if(fabs((*dau_eta)[kjet][A_trk]) > 2.4 ) continue;

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
                    
                    double jet_dau_pt = ptWRTJet((double)(*jetPt)[kjet], (double)(*jetEta)[kjet], (double)(*jetPhi)[kjet], (double)(*dau_pt)[kjet][A_trk], (double)(*dau_eta)[kjet][A_trk], (double)(*dau_phi)[kjet][A_trk]);
                    double jet_dau_eta   = etaWRTJet((double)(*jetPt)[kjet], (double)(*jetEta)[kjet] , (double)(*jetPhi)[kjet]  , (double)(*dau_pt)[kjet][A_trk], (double)(*dau_eta)[kjet][A_trk], (double)(*dau_phi)[kjet][A_trk]);
                    
                    double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[kjet][A_trk] , (*dau_eta)[kjet][A_trk] )));
                    for(int i = 0; i < trackbin; i++){
                        if(tkBool[i]==1){
                            reco_h_jet_cor_jT[i]->Fill(jet_dau_pt,1.0*genWeightPy/(Atrk_weight));
                            reco_h_jet_cor_etastar[i]->Fill(jet_dau_eta,1.0*genWeightPy/(Atrk_weight));
                            reco_h_jet_unc_jT[i]->Fill(jet_dau_pt,1.0*genWeightPy);
                            reco_h_jet_unc_etastar[i]->Fill(jet_dau_eta,1.0*genWeightPy);
                        } 
                    }
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
    TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/xiaoyul/MC/13TeV/closure/job_%s.root",subList.c_str()), "recreate");
   
    hpthat->Write();
    hleadingJetPt->Write(); 
    hjetPtReco->Write();
    hjetPtGen->Write();

    reco_hBinDist_cor_single->Write();
    reco_hBinDist_unc_single->Write();
    gen_hBinDist_unc_single->Write();

    for(int i=0;i<trackbin;i++){
        reco_h_jet_unc_jT[i]->Write();
        reco_h_jet_unc_etastar[i]->Write();
        reco_h_jet_cor_jT[i]->Write();
        reco_h_jet_cor_etastar[i]->Write();
        reco_hBinDist_cor[i]->Write();
        reco_hBinDist_unc[i]->Write();

        gen_h_jet_unc_jT[i]->Write();
        gen_h_jet_unc_etastar[i]->Write();
        gen_hBinDist_unc[i]->Write();
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

