#define MyClass_cxx

#include "include/MyTrim.h"
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

//Bool folds

bool F_eventpass(std::vector< float > *jetPt, int jetnumber, float jetPtCut){

    if(jetPt->size() < 1)           return false;
    if(jetnumber > 200)              return false;
    if((*jetPt)[0] < jetPtCut) return false;
    return true;
}

bool F_jetpass(std::vector< float > * jetEta, std::vector< float > * jetPt, int     ijet, float   jetPtCut){
    if(fabs( (*jetEta)[ijet])   >jetEtaCut)   return false;
    if((*jetPt)[ijet]           <jetPtCut)    return false;
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


void MyClass::Loop(int job, std::string fList){

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    const int i40 = 40;
    const int i470 = 470;
    const int i3500 = 3500;
    const int i50 = 50;
    const int i500 = 500;
    const int i1500 = 1500;
    const int i150 = 150;
    const int i100 = 100;
    const int PUbin = 1;
    double track_eta_lim = 5.0;
    //Initializing Histograms
    TH2D* hPairs        = new TH2D("hPairs","hPairs",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);
    TH2D* hNtrig        = new TH2D("hNtrig","hNtrig",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);
    TH2D* hNtrigCorrected        = new TH2D("hNtrigCorrected","hNtrigCorrected",  trackbin, bin0,trackbin, ptbin, bin0, ptbin);


    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,bin0,trackbin);
    TH1D* hJet_Pass550     = new TH1D("hJet_Pass550"  ,"hJet_Pass550"  , trackbin,bin0,trackbin);
    TH1D* hJet_Pass550_hltCor     = new TH1D("hJet_Pass550_hltCor"  ,"hJet_Pass550_hltCor"  , trackbin,bin0,trackbin);
    TH1D* hJetPt = new TH1D("hJetPt"  ,"hJetPt" ,i150,i100,i3500);
    TH1D* hBinDist_cor_single = new TH1D("hBinDist_cor_single","hBinDist_cor_single",bin360,bin0,bin120);
    TH1D* hBinDist_unc_single = new TH1D("hBinDist_unc_single","hBinDist_unc_single",bin120,bin0,bin120);


    TH2D* hEPDrawCor[trackbin][ptbin][PUbin];
    TH2D* hSignalShiftedCor[trackbin][ptbin][PUbin];
    TH2D* hBckrndShiftedCor[trackbin][ptbin][PUbin];
    TH1D* hTotalWeight                = new TH1D("hTotalWeight","hTotalWeight"                          , bin500,-bin2,bin3);
    TH1D* hAvg_Atrk_Weight            = new TH1D("hAvg_Atrk_Weight","hAvg_Atrk_Weight"                  , bin500,-bin2,bin3);
    TH1D* hAvg_NtrigCorrected_Bump    = new TH1D("hAvg_NtrigCorrected_Bump","hAvg_NtrigCorrected_Bump"  , bin22 ,-bin2,bin20);

    TH1D* hBinDist_cor[trackbin];
    TH1D* hBinDist_unc[trackbin];


    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_cor[wtrk-1]    = new TH1D(Form("hBinDist_cor_%d",wtrk),Form("hBinDist_cor_%d",wtrk), bin360, bin0, bin120);
        hBinDist_unc[wtrk-1]    = new TH1D(Form("hBinDist_unc_%d",wtrk),Form("hBinDist_unc_%d",wtrk), bin120, bin0, bin120);
        for(int wppt = 1; wppt<ptbin+1; wppt++){
            for(int wpPU = 1; wpPU<PUbin+1; wpPU++){
                hBckrndShiftedCor[wtrk-1][wppt-1][wpPU-1]   = new TH2D(Form("hBckrndS_Cor_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form("hBckrndS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                hSignalShiftedCor[wtrk-1][wppt-1][wpPU-1]   = new TH2D(Form("hSignalS_Cor_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form("hSignalS_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                hEPDrawCor[wtrk-1][wppt-1][wpPU-1]          = new TH2D(Form( "hEPDraw_Cor_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form( "hEPDraw_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
            }
        }
    }

    //NEW THING                        
    //NEW THING                        
    //NEW THING                   
    /*      
    TH2D* hReco2D[fileList.size()];
    TH2D* hGen2D[fileList.size()];
    TH1D* hdid500;
    TH1D* hdid400;

    TFile *f_jet_HLT_lookup = new TFile("did400500_v2_all.root");
    hdid500 = (TH1D*)f_jet_HLT_lookup->Get("d500")->Clone("did500");
    hdid400 = (TH1D*)f_jet_HLT_lookup->Get("d400")->Clone("did400");
    hdid500->Divide(hdid400);
    */
    // MAIN CODE BEGINS
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

        /* Below is all about matching MC to Data era
        vector<string> era_vec;
        vector<string> matched_cor_table_vec;
        era_vec.push_back("2018D_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018C_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018B_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2018A_"); matched_cor_table_vec.push_back("18_v2");
        era_vec.push_back("2017B_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017C_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017D_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017E_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2017F_"); matched_cor_table_vec.push_back("17_v2");
        era_vec.push_back("2016Bv2_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016C_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016D_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016E_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016F-HIPM_"); matched_cor_table_vec.push_back("16apv_v2");
        era_vec.push_back("2016F_"); matched_cor_table_vec.push_back("16_v2");
        era_vec.push_back("2016G_"); matched_cor_table_vec.push_back("16_v2");
        era_vec.push_back("2016H_"); matched_cor_table_vec.push_back("16_v2");
        int i_keep=999;
        for(int i=0; i<era_vec.size(); i++){
            if(isSubstring(era_vec[i],fileList.at(f).c_str())){
                i_keep=i;
                continue;
            }
        }

        string filename = "~/StorageArea/"+matched_cor_table_vec[i_keep]+".root";
        TFile *f_pt_eta_DCA_lookup = new TFile(filename.c_str());
        
        hReco2D[f] = (TH2D*)f_pt_eta_DCA_lookup->Get("h2_MAIN_DCA_Dau_Reco_Pt_Eta_Lab_All")->Clone(Form("h2_DCA_Dau_Reco_Pt_Eta_Lab_All_%d",f));
        hGen2D[f] = (TH2D*)f_pt_eta_DCA_lookup->Get("h2_Dau_Gen_Pt_Eta_Lab_All")->Clone(Form("h2_Dau_Gen_Pt_Eta_Lab_All_%d",f));
        hReco2D[f]->Divide(hGen2D[f]);
        int thisEffTable =f_from_file;

        //Above is all about matching MC to Data era */

        //========ENTERING EVENT LOOP========
        for (Long64_t ievent=0; ievent <nentries; ievent ++){
            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;

            //cut on jetPt and jetN (number of [HLT-passed] jets in an evnet)
            if(!F_eventpass(jetPt, jetN, jetPtCut_Event)){
                continue;
            }

            hEvent_Pass->Fill(1);

            int jetCounter = jetPt->size();
            if(jetCounter == 0) continue;

            //========ENTERING JET LOOP========
            for(int kjet=0; kjet < jetCounter; kjet++){

                int ijet = kjet; //indicesG[kjet];
                long int NNtrk = (dau_pt->at(ijet)).size();
                //long int NNtrk = (dau_pt_STAR->at(ijet)).size();
                gRandom->SetSeed(0);
                double eta_smear;

                if( fabs(((*jetEta)[kjet])+eta_smear) > jetEtaCut ) continue;
                if( (*jetPt)[kjet] < jetPtCut_Jet) continue;
                hJetPt->Fill((*jetPt)[kjet]);
                // filling distrivutions within track bins
                // ALSO VERY IMPORTANLTY changing the tkBool to 1 for this particular jet. This will be usefull later wen I create conditons for filling other historgams.

                double jet_HLT_weight = 1.0;
                 /*
                 //COMMENT OUT NEXT 3 LINES
                 //HIGH LEVEL TRIGGER "is there a jet with cone size 8 and pt above X"
                 if((*jetPt)[ijet] < 880 && (*jetPt)[ijet] > 550){
                     jet_HLT_weight = 1.0/(hdid500->GetBinContent(hdid500->FindBin((*jetPt)[ijet]) ) );
                 }
                */

                int tkBool[trackbin] = {0};

                int n_ChargeMult_DCA_labPt_Eta_exclusion =0;
                double n_ChargeMult_DCA_labPt_Eta_exclusion_Cor =0;
                
                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    if((*dau_chg)[kjet][A_trk] == 0) continue;//charge
                    if(fabs((*dau_pt)[kjet][A_trk])  < 0.3)     continue;//lab pt
                    if(fabs((*dau_eta)[kjet][A_trk]) > 2.4)     continue;//lab eta

                    double dauptnow = (*dau_pt)[kjet][A_trk]; double pterr_R_pt =0; double dcaXY = 0; double dcaZ = 0;
                    pterr_R_pt= ( (*dau_ptError)[kjet][A_trk] ) / ( (*dau_pt)[kjet][A_trk] );
                    dcaXY = (*dau_XYDCAsig)[kjet][A_trk];
                    dcaZ = (*dau_ZDCAsig)[kjet][A_trk];

                    if(pterr_R_pt   > 0.1   && dauptnow > 0.5)  continue;
                    if(fabs(dcaXY)  > 3     && dauptnow > 0.5)  continue;
                    if(fabs(dcaZ)   > 3     && dauptnow > 0.5)  continue;
                    n_ChargeMult_DCA_labPt_Eta_exclusion += 1;


                    //ABOVE is normal mult after cuts
                    // Comment below and replace 
                    //          "n_ChargeMult_DCA_labPt_Eta_exclusion_Cor"
                    //          with
                    //          "n_ChargeMult_DCA_labPt_Eta_exclusion"
                    /*
                    double nUnc_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    n_ChargeMult_DCA_labPt_Eta_exclusion_Cor += (1.0/nUnc_weight);
                    */
                }
                /*
                hBinDist_cor_single->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor, 1.0*jet_HLT_weight);
                */
                hBinDist_unc_single->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion, 1.0*jet_HLT_weight);
                for(int i = 0; i < trackbin; i++){
                    if( n_ChargeMult_DCA_labPt_Eta_exclusion >= trackbinbounds[i] && n_ChargeMult_DCA_labPt_Eta_exclusion < trackbinboundsUpper[i]){
                        tkBool[i] = 1;
                        hJet_Pass->Fill(i);
                        if((*jetPt)[kjet] >= jetPtCut_Event) hJet_Pass550->Fill(i);
                        if((*jetPt)[kjet] >= jetPtCut_Event) hJet_Pass550_hltCor->Fill(i,1.0*jet_HLT_weight);
                        /*
                        hBinDist_cor[i]->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor,1.0*jet_HLT_weight);
                        */
                        hBinDist_unc[i]->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion,1.0*jet_HLT_weight);
                    }
                }

                int Ntrig[trackbin][ptbin] = {0};
                double NtrigCorrected[trackbin][ptbin] = {0};
                int A_ptBool[NNtrk][ptbin] = {0};

                //****** TWO PARTICLE MULT CODE *******
                // VERY IMPORTANT calculating daughter pt wrt to jet axis.
                // So this needs to be 2d vector, for pt bin and for daughter index.
                // In each case I create a true falsse for that daughter falling in to the specific pt bin.
                // Note this is NOT jet x daughter. It's pt bin x daughtr

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

                    for(int i = 0; i < ptbin; i++){
                        if(jet_dau_pt >= ptbinbounds_lo[i] && jet_dau_pt < ptbinbounds_hi[i]){
                            A_ptBool[A_trk][i] = 1;
                        }
                    }
                    // 0.0 - 3, 0.3 - 3, 2.0 - 3
                    // some particle has pt 1
                    // A_ptBool = [1,1,0]
                    //
                    // 0-1000, 60+, 80+
                    // Jet with mult 90
                    // tkBool [1,1,1]

                    //NEW THING                        
                    //NEW THING                        
                    //NEW THING                        

                    //getting et pt efficiency for A track in beam frame
                    //COMMENT OUT IF YOU DONT HAVE MC TABLES
                    /*
                    double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    */

                    double Atrk_weight = 1.0;

                    //in the case where the total jet mult and the individual daughter pt is acceptable for this track bin and pt bin, we increase the Ntrig count.
                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                NtrigCorrected[i][j] += (1.0/Atrk_weight);
                                Ntrig[i][j] += 1;
                            }
                        }
                    }

                }//here ends the boolean array creations.

                //continuation of main loops. Here is where the 2D Corr plots are created using the above booleans 

                for(int i = 0; i < trackbin; i++){
                    for(int j = 0; j < ptbin; j++){
                        hNtrig->Fill(i,j,Ntrig[i][j]);
                        //hNtrigCorrected->Fill(i,j, NtrigCorrected[i][j]);
                        //hAvg_NtrigCorrected_Bump->Fill(NtrigCorrected[i][j] - Ntrig[i][j]);
                    }
                }

                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    //this should be redundant if it passes the bools above? i guess it helps skip daughters faster. maybe i can reindex and run through the daughters quickly by aranging all the charged dauhghter sat the front.
                    if((*dau_chg)[ijet][A_trk] == 0) continue;
                    if((*dau_pt)[ijet][A_trk] < 0.3) continue;
                    if(fabs((*dau_eta)[ijet][A_trk]) > 2.4) continue;
                    //if((*highPurity)[ijet][A_trk] == 0) continue;

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

                    gRandom->SetSeed(0);
                    double phi_smear=0;

                    double jet_dau_pt    =  ptWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_eta   = etaWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);

                    double jet_dau_theta = 2*ATan(Exp(-(jet_dau_eta)));
                    if(jet_dau_eta > track_eta_lim) continue;
                    
                    //getting eta pt efficiency for A track in beam frame
                    double Atrk_weight = 1.0; //(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    hAvg_Atrk_Weight->Fill(1.0/(Atrk_weight));

                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                int k_PU=0;
                                hEPDrawCor[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight/(Atrk_weight * NtrigCorrected[i][j] ));
                            }
                        }
                    }
                    
                    if(A_trk == NNtrk - 1) continue;

                    for(long int T_trk=A_trk+1; T_trk< NNtrk; T_trk++ ){

                        if((*dau_chg)[ijet][T_trk] == 0) continue;
                        if((*dau_pt)[ijet][T_trk] < 0.3) continue;
                        if(fabs((*dau_eta)[ijet][T_trk]) > 2.4) continue;

                        double T_dauptnow = (*dau_pt)[ijet][T_trk];
                        double T_pterr_R_pt =0;
                        double T_dcaXY = 0;
                        double T_dcaZ = 0;
                        T_pterr_R_pt= ( (*dau_ptError)[ijet][T_trk] ) / ( (*dau_pt)[ijet][T_trk] );
                        T_dcaXY = (*dau_XYDCAsig)[ijet][T_trk];
                        T_dcaZ = (*dau_ZDCAsig)[ijet][T_trk];
                        if(T_pterr_R_pt   > 0.1   && T_dauptnow > 0.5) continue;
                        if(fabs(T_dcaXY)  > 3     && T_dauptnow > 0.5) continue;
                        if(fabs(T_dcaZ)   > 3     && T_dauptnow > 0.5) continue;


                        //We use E Scheme jet axis, others include "winner take all"
                        double T_jet_dau_pt    =  ptWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);
                        double T_jet_dau_eta   = etaWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);
                        double T_jet_dau_phi   = phiWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);

                        if(T_jet_dau_eta > track_eta_lim) continue;

                        //ETA AND PHI DIFFERNCES
                        double deltaEta = (jet_dau_eta - T_jet_dau_eta);
                        double deltaPhi = (TMath::ACos(TMath::Cos(jet_dau_phi - T_jet_dau_phi)));
                        //NEW THING                        
                        //NEW THING                        
                        //NEW THING                        
                        //getting et pt efficiency for T track in beam frame
                        double Ttrk_weight = 1.0;//(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][T_trk] , (*dau_eta)[ijet][T_trk] )));

                        // jet mults of 0-30,30-69,90+
                        // pt 0-1,1-2,2-3
                        // BELOW, i=3, j=3

                        for(int i = 0; i < trackbin; i++){
                            for(int j = 0; j < ptbin; j++){
                                if(tkBool[i] + A_ptBool[A_trk][j] + A_ptBool[T_trk][j] == 3){

                                    //cout << "1.0/(Atrk_weight*Ttrk_weight*NtrigCorrected[i][j]) " << 1.0/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ) << " and components: " << Atrk_weight << " " << Ttrk_weight << " " << NtrigCorrected[i][j] << endl;

                                    hPairs->Fill(i,j);

                                    if(Atrk_weight == 0){/* cout  << "Atrk_weight = 0!" << endl;
                                    */
                                        continue;}
                                    if(Ttrk_weight == 0){/* cout << "Ttrk_weight = 0!" << endl;
                                    */
                                        continue;}
                                    if(NtrigCorrected[i][j] == 0){ cout << "NtrigCorrected[i][j] = 0!" << endl; continue;}

                                    hTotalWeight->Fill(1.0 /(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));

                                    int k_PU=0;
                                    hSignalShiftedCor[i][j][k_PU]->Fill(deltaEta, deltaPhi,                  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                    hSignalShiftedCor[i][j][k_PU]->Fill(-deltaEta, deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                    hSignalShiftedCor[i][j][k_PU]->Fill(deltaEta, -deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                    hSignalShiftedCor[i][j][k_PU]->Fill(-deltaEta, -deltaPhi,                1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                    hSignalShiftedCor[i][j][k_PU]->Fill( deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                    hSignalShiftedCor[i][j][k_PU]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                }
                            }
                        }
                    }//T_trk end
                }//A_trk end
            }
        }
        fFile->Close();
        //ABOVE loops every single root file in fList and produce the entire signal dist
        //it also makes EPDraw, single particle 2D eta phi dist for ALL jets in ALL events, within a jet mult class

        //EPDraw takes all particles from every event that meets our selection, jet mult, daughter pT, etc
        //  combines them all into one 2D distribution of Eta-Phi (with respect to jet axis)
        //  later for background we will pull out particles from the EPDraw distribution
        //  when we pull particles from EPDraw we are selecting partciles from ALL events and so there should be no single-jet correlations

        int backMult =10;
        for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
            std::cout << wtrk << "/" << trackbin << std::endl;
            for(int wppt = 1; wppt < ptbin+1; wppt++){ 
                std::cout << wppt << "/" << ptbin << std::endl;
                for(int wpPU = 1; wpPU < PUbin+1; wpPU++){
                    std::cout << wpPU << "/" << PUbin << std::endl;
                    //Nent is the number of pairs in the signal which we will try to 10x
                    //Xent is the number of pseudoparticles requried such that when we build the pairs nCp = Xent CHOOSE 2 will give 
                    //us 10 times as many pairs as we have in the signal histogrm.

                    //long int NENT =  hPairsPU[wpPU-1]->GetBinContent(wtrk, wppt);
                    long int NENT =  hPairs->GetBinContent(wtrk, wppt);
                    long int XENT =  ((1+floor(sqrt(1+(4*2*backMult*NENT))))/2) ;
                    //float A_ETA[XENT] = {0};
                    //float A_PHI[XENT] = {0};
                    float A_ETA_Cor[XENT] = {0};
                    float A_PHI_Cor[XENT] = {0};
                    for(int x = 0; x<XENT; x++){
                        gRandom->SetSeed(0);

                        double WEta1_Cor, WPhi1_Cor;//making the pseudoparticles
                        hEPDrawCor[wtrk-1][wppt-1][wpPU-1]->GetRandom2(WEta1_Cor, WPhi1_Cor);
                        A_ETA_Cor[x] = WEta1_Cor;
                        A_PHI_Cor[x] = WPhi1_Cor;
                    }
                    for(long int i = 0; i < (XENT-1); i++){
                        for(long int j = (i+1); j < XENT; j++){

                            double WdeltaEta_Cor = (A_ETA_Cor[i]-A_ETA_Cor[j]);
                            double WdeltaPhi_Cor = (TMath::ACos(TMath::Cos(A_PHI_Cor[i]-A_PHI_Cor[j])));
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaEta_Cor, WdeltaPhi_Cor,   1);//./XENT);
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpPU-1]->Fill(-WdeltaEta_Cor, WdeltaPhi_Cor,  1);//./XENT);
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaEta_Cor, -WdeltaPhi_Cor,  1);//./XENT);
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpPU-1]->Fill(-WdeltaEta_Cor, -WdeltaPhi_Cor, 1);//./XENT);
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpPU-1]->Fill(WdeltaEta_Cor, 2*TMath::Pi() - WdeltaPhi_Cor, 1);//./XENT);
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpPU-1]->Fill(-WdeltaEta_Cor,2*TMath::Pi() - WdeltaPhi_Cor, 1);//./XENT);

                        }
                    }
                }
            }
        }

        //We have theBck and the Sig for all jets
        //closed all the root files

        string subList = fList.substr(fList.size() - 3);
        TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/flowcorr/root_out_1d2d/round3/job_%s_%d.root",subList.c_str(),f), "recreate");
        for(int wtrk =1; wtrk <trackbin+1; wtrk++){
            hBinDist_cor[wtrk-1]->Write();
            hBinDist_unc[wtrk-1]->Write();
            for(int wppt =1; wppt <ptbin+1; wppt++){
                for(int wpPU =1; wpPU<PUbin+1; wpPU++){
                    // In "1d2d_constants.h" we 3 jet mults, 2 pT bins
                    // 0-30,30-60,60+
                    // 0-3 GeV, 0.3 - 3 Gev
                    // hSignalShiftedCor_0_to_30_and_0_3
                    // hSignalShiftedCor_30_to_60_and_0_3

                    hSignalShiftedCor[wtrk-1][wppt-1][wpPU-1]->Write(Form("hSigS_Cor_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ));
                    hBckrndShiftedCor[wtrk-1][wppt-1][wpPU-1]->Write(Form("hBckS_Cor_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU ));
                    hEPDrawCor[wtrk-1][wppt-1][wpPU-1]->Write(Form("hEPD_Cor_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1]),(int)(10*ptbinbounds_hi[wppt-1]),wpPU  ));
                }
            }
        }
        hPairs->Write();
        hNtrig->Write();
        hNtrigCorrected->Write();
        hTotalWeight->Write();
        hAvg_Atrk_Weight->Write();
        hAvg_NtrigCorrected_Bump->Write();
        hJetPt->Write();
        hEvent_Pass->Write();
        hJet_Pass->Write();
        hJet_Pass550->Write();
        hJet_Pass550_hltCor->Write();
        hBinDist_cor_single->Write();
        hBinDist_unc_single->Write();
        
        fS_tempA->Close();

        for(int wtrk =1; wtrk <trackbin+1; wtrk++){
            hBinDist_cor[wtrk-1]->Reset();
            hBinDist_unc[wtrk-1]->Reset();
            for(int wppt =1; wppt <ptbin+1; wppt++){
                for(int wpPU =1; wpPU<PUbin+1; wpPU++){                 
                    // In "1d2d_constants.h" we 3 jet mults, 2 pT bins
                    // 0-30,30-60,60+
                    // 0-3 GeV, 0.3 - 3 Gev
                    // hSignalShiftedCor_0_to_30_and_0_3
                    // hSignalShiftedCor_30_to_60_and_0_3

                    hSignalShiftedCor[wtrk-1][wppt-1][wpPU-1]->Reset();
                    hBckrndShiftedCor[wtrk-1][wppt-1][wpPU-1]->Reset();
                    hEPDrawCor[wtrk-1][wppt-1][wpPU-1]->Reset();
                }
            }
        }
        hPairs->Reset();
        hNtrig->Reset();
        hNtrigCorrected->Reset();
        hTotalWeight->Reset();
        hAvg_Atrk_Weight->Reset();
        hAvg_NtrigCorrected_Bump->Reset();
        hJetPt->Reset();
        hEvent_Pass->Reset();
        hJet_Pass->Reset();
        hJet_Pass550->Reset();
        hJet_Pass550_hltCor->Reset();
        hBinDist_cor_single->Reset();
        hBinDist_unc_single->Reset(); 
    }//end looping over files
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
    std::ifstream inFile(fList.data());//open the file 

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

