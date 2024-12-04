#define MyClass_cxx

#include "include/MyTrim.h"
#include "include/coordinateTools.h"
#include "include/1d2d_constants_jT.h"

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

    bool smear_bool = 0;

    const int i40 = 40;
    const int i470 = 470;
    const int i3500 = 3500;
    const int i50 = 50;
    const int i500 = 500;
    const int i1500 = 1500;
    const int i150 = 150;
    const int i100 = 100;
    const int i340 =340;
    const int PUbin = 1;
    double track_eta_lim = 5.0;
    //Initializing Histograms
    TH2D* hPairs[ptbin_A];
    for(int i=0;i<ptbin_A;i++){
        hPairs[i]= new TH2D(Form("hPairs_A_%d",i),Form("hPairs_A_%d",i),  trackbin, bin0,trackbin, ptbin_T, bin0, ptbin_T);
    }
    TH2D* hNtrig        = new TH2D("hNtrig","hNtrig",  trackbin, bin0,trackbin, ptbin_T, bin0, ptbin_T);
    TH2D* hNtrigCorrected        = new TH2D("hNtrigCorrected","hNtrigCorrected",  trackbin, bin0,trackbin, ptbin_T, bin0, ptbin_T);

    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,bin0,trackbin);
    TH1D* hJet_Pass550     = new TH1D("hJet_Pass550"  ,"hJet_Pass550"  , trackbin,bin0,trackbin);
    TH1D* hJet_Pass550_hltCor     = new TH1D("hJet_Pass550_hltCor"  ,"hJet_Pass550_hltCor"  , trackbin,bin0,trackbin);
    TH1D* hJetPt = new TH1D("hJetPt"  ,"hJetPt" ,i340,i100,i3500);
    TH1D* hJetPt_wo_ptcut = new TH1D("hJetPt_wo_ptcut"  ,"hJetPt_wo_ptcut" ,i340,i100,i3500);
    TH1D* hBinDist_cor_single = new TH1D("hBinDist_cor_single","hBinDist_cor_single",bin360,bin0,bin120);
    TH1D* hBinDist_unc_single = new TH1D("hBinDist_unc_single","hBinDist_unc_single",bin120,bin0,bin120);

    TH2D* h_lab_JetMult_pT=new TH2D("lab_JetMult_pT","lab_JetMult_pT",340,100,3500,120,0.5,120.5);
    TH2D* h_lab_JetMult_phi=new TH2D("lab_JetMult_phi","lab_JetMult_phi",30,-TMath::Pi(),TMath::Pi(),120,0.5,120.5);
    TH2D* h_lab_JetMult_eta=new TH2D("lab_JetMult_eta","lab_JetMult_eta",34,-1.7,1.7,120,0.5,120.5);

    TH2D* h_lab_cor_JetMult_pT=new TH2D("lab_cor_JetMult_pT","lab_cor_JetMult_pT",340,100,3500,360,0.5,120.5);
    TH2D* h_lab_cor_JetMult_phi=new TH2D("lab_cor_JetMult_phi","lab_cor_JetMult_phi",30,-TMath::Pi(),TMath::Pi(),360,0.5,120.5);
    TH2D* h_lab_cor_JetMult_eta=new TH2D("lab_cor_JetMult_eta","lab_cor_JetMult_eta",34,-1.7,1.7,360,0.5,120.5);
    
    TH1D* h_jet_jT=new TH1D("jet_jT","jet_jT",200,0,20);
    TH1D* h_jet_etastar=new TH1D("jet_etastar","jet_etastar",100,0,10);
    TH1D* h_jet_cor_jT=new TH1D("jet_cor_jT","jet_cor_jT",200,0,20);
    TH1D* h_jet_cor_etastar=new TH1D("jet_cor_etastar","jet_cor_etastar",100,0,10);

    TH2D* hEPDrawCor_T[trackbin][ptbin_T][PUbin];
    TH2D* hEPDrawCor_A[trackbin][ptbin_A][PUbin];
    TH2D* hSignalShiftedCor[trackbin][ptbin_T][ptbin_A][PUbin];
    TH2D* hBckrndShiftedCor[trackbin][ptbin_T][ptbin_A][PUbin];
    TH2D* hEPDrawUnc_T[trackbin][ptbin_T][PUbin];
    TH2D* hEPDrawUnc_A[trackbin][ptbin_A][PUbin];
    TH2D* hSignalShiftedUnc[trackbin][ptbin_T][ptbin_A][PUbin];
    TH2D* hBckrndShiftedUnc[trackbin][ptbin_T][ptbin_A][PUbin];
    
    TH1D* hTotalWeight                = new TH1D("hTotalWeight","hTotalWeight"                          , bin500,-bin2,bin3);
    TH1D* hAvg_Atrk_Weight            = new TH1D("hAvg_Atrk_Weight","hAvg_Atrk_Weight"                  , bin500,-bin2,bin3);
    TH1D* hAvg_NtrigCorrected_Bump    = new TH1D("hAvg_NtrigCorrected_Bump","hAvg_NtrigCorrected_Bump"  , bin22 ,-bin2,bin20);

    //For recording Nch dis. and calculating <Nch> per Nch bin.
    TH1D* hBinDist_cor[trackbin];
    TH1D* hBinDist_unc[trackbin];

    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_cor[wtrk-1]    = new TH1D(Form("hBinDist_cor_%d",wtrk),Form("hBinDist_cor_%d",wtrk), bin360, bin0, bin120);
        hBinDist_unc[wtrk-1]    = new TH1D(Form("hBinDist_unc_%d",wtrk),Form("hBinDist_unc_%d",wtrk), bin360, bin0, bin120);
        for(int wpPU = 1; wpPU<PUbin+1; wpPU++){
            for(int wppt = 1; wppt<ptbin_T+1; wppt++){
                for(int wpptA = 1; wpptA<ptbin_A+1; wpptA++){
                    hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]   = new TH2D(Form("hBkg_Cor_trk_%d_ptT_%d_ptA_%d_PU_%d",wtrk,wppt,wpptA,wpPU) ,
                    Form("hBkg_Cor_%d_to_%d_and_T_%d_to_%d_A_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_T_lo[wppt-1]),(int)(10*pt_T_hi[wppt-1]),(int)(10*pt_A_lo[wpptA-1]),(int)(10*pt_A_hi[wpptA-1]),wpPU),
                    41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                    
                    hSignalShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]   = new TH2D(Form("hSig_Cor_trk_%d_ptT_%d_ptA_%d_PU_%d",wtrk,wppt,wpptA,wpPU) ,
                    Form("hSig_Cor_%d_to_%d_and_T_%d_to_%d_A_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_T_lo[wppt-1]),(int)(10*pt_T_hi[wppt-1]),(int)(10*pt_A_lo[wpptA-1]),(int)(10*pt_A_hi[wpptA-1]),wpPU) ,
                    41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);


                    hBckrndShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]   = new TH2D(Form("hBkg_Unc_trk_%d_ptT_%d_ptA_%d_PU_%d",wtrk,wppt,wpptA,wpPU) ,
                    Form("hBkg_Unc_%d_to_%d_and_T_%d_to_%d_A_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_T_lo[wppt-1]),(int)(10*pt_T_hi[wppt-1]),(int)(10*pt_A_lo[wpptA-1]),(int)(10*pt_A_hi[wpptA-1]),wpPU),
                    41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                    
                    hSignalShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]   = new TH2D(Form("hSig_Unc_trk_%d_ptT_%d_ptA_%d_PU_%d",wtrk,wppt,wpptA,wpPU) ,
                    Form("hSig_Unc_%d_to_%d_and_T_%d_to_%d_A_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_T_lo[wppt-1]),(int)(10*pt_T_hi[wppt-1]),(int)(10*pt_A_lo[wpptA-1]),(int)(10*pt_A_hi[wpptA-1]),wpPU) ,
                    41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                }
            }
            
            for(int wppt=1;wppt<ptbin_T+1;wppt++){
                hEPDrawCor_T[wtrk-1][wppt-1][wpPU-1]          = new TH2D(Form( "hEPDraw_Cor_T_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form( "hEPDraw_Cor_T_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_T_lo[wppt-1]),(int)(10*pt_T_hi[wppt-1]),wpPU) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
                hEPDrawUnc_T[wtrk-1][wppt-1][wpPU-1]          = new TH2D(Form( "hEPDraw_Unc_T_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form( "hEPDraw_Unc_T_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_T_lo[wppt-1]),(int)(10*pt_T_hi[wppt-1]),wpPU) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
            }
            
            for(int wpptA=1;wpptA<ptbin_A+1;wpptA++){
                hEPDrawCor_A[wtrk-1][wpptA-1][wpPU-1]          = new TH2D(Form( "hEPDraw_Cor_A_trk_%d_ppt_%d_PU_%d",wtrk,wpptA,wpPU) ,Form( "hEPDraw_Cor_A_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_A_lo[wpptA-1]),(int)(10*pt_A_hi[wpptA-1]),wpPU) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
                hEPDrawUnc_A[wtrk-1][wpptA-1][wpPU-1]          = new TH2D(Form( "hEPDraw_Unc_A_trk_%d_ppt_%d_PU_%d",wtrk,wpptA,wpPU) ,Form( "hEPDraw_Unc_A_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_A_lo[wpptA-1]),(int)(10*pt_A_hi[wpptA-1]),wpPU) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
            }
        }
    }

    TH2D* hjet_eta_phi_before_veto=new TH2D(Form("jet_eta_phi_before_veto"),Form("jet_eta_phi_before_veto"),82,-5.191,5.191,72,-TMath::Pi(),TMath::Pi());
    TH2D* hjet_eta_phi_after_veto=new TH2D(Form("jet_eta_phi_after_veto"),Form("jet_eta_phi_after_veto"),82,-5.191,5.191,72,-TMath::Pi(),TMath::Pi());


    TH2D* hReco2D[fileList.size()];
    TH2D* hGen2D[fileList.size()];
    TH2D* jet_veto_map[fileList.size()];
    /*
    TH1D* hdid500;
    TH1D* hdid400;

    TFile *f_jet_HLT_lookup = new TFile("did400500_v2_all.root");
    hdid500 = (TH1D*)f_jet_HLT_lookup->Get("d500")->Clone("did500");
    hdid400 = (TH1D*)f_jet_HLT_lookup->Get("d400")->Clone("did400");
    hdid500->Divide(hdid400);
    */

    std::cout << "Starting event loop" << std::endl;
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;
    for(int f = 0; f<fileList.size(); f++){
        int f_from_file = f;
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("analyzerOffline/trackTree");
        Init(tree);

        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;

        vector<string> era_vec;
        vector<string> matched_cor_table_vec;
        vector<string> matched_jet_veto_vec;
        era_vec.push_back("2022C"); matched_cor_table_vec.push_back("correction_2022CD"); matched_jet_veto_vec.push_back("Summer22_23Sep2023_RunCD_v1");
        era_vec.push_back("2022D"); matched_cor_table_vec.push_back("correction_2022CD"); matched_jet_veto_vec.push_back("Summer22_23Sep2023_RunCD_v1");
        era_vec.push_back("2022E"); matched_cor_table_vec.push_back("correction_2022EFG"); matched_jet_veto_vec.push_back("Summer22EE_23Sep2023_RunEFG_v1");
        era_vec.push_back("2022F"); matched_cor_table_vec.push_back("correction_2022EFG"); matched_jet_veto_vec.push_back("Summer22EE_23Sep2023_RunEFG_v1");
        era_vec.push_back("2022G"); matched_cor_table_vec.push_back("correction_2022EFG"); matched_jet_veto_vec.push_back("Summer22EE_23Sep2023_RunEFG_v1");
        era_vec.push_back("2023C"); matched_cor_table_vec.push_back("correction_2023C");
        era_vec.push_back("2023D"); matched_cor_table_vec.push_back("correction_2023D");
        era_vec.push_back("2024B"); matched_cor_table_vec.push_back("correction_2024");
        era_vec.push_back("2024C"); matched_cor_table_vec.push_back("correction_2024");
        era_vec.push_back("2024D"); matched_cor_table_vec.push_back("correction_2024");
        era_vec.push_back("2024E"); matched_cor_table_vec.push_back("correction_2024");
        era_vec.push_back("2024F"); matched_cor_table_vec.push_back("correction_2024");
        era_vec.push_back("2024G"); matched_cor_table_vec.push_back("correction_2024");
        era_vec.push_back("2024H"); matched_cor_table_vec.push_back("correction_2024");
        era_vec.push_back("2024I"); matched_cor_table_vec.push_back("correction_2024");
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

        bool applyjetveto=0;
        if(i_keep<5){
            applyjetveto=1;
            string jetvetofilename = "~/StorageArea/"+matched_jet_veto_vec[i_keep]+".root";
            TFile *f_jet_veto_lookup = new TFile(jetvetofilename.c_str());
            jet_veto_map[f]=(TH2D*)f_jet_veto_lookup->Get("jetvetomap");
        }

        
        // ENTERING EVENT LOOP
        for (Long64_t ievent=0; ievent <nentries; ievent ++){
            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;

            if(!F_eventpass(jetPt, jetN, jetPtCut_Event)){
                continue;
            }
            
            hEvent_Pass->Fill(1);

            int jetCounter = jetPt->size();
            if(jetCounter == 0) continue;

            if(applyjetveto){
                //apply jet veto map
                //jet veto, reference: https://cms-jerc.web.cern.ch/Recommendations/#run-3
                //"The safest procedure would be to veto events if ANY jet with a loose selection lies in the veto regions."
                bool jetvetoBool=0;
                for(int ijet=0; ijet < jetCounter; ijet++){
                    hjet_eta_phi_before_veto->Fill((*jetEta)[ijet],(*jetPhi)[ijet]);
                    
                    //if( (*jetPt)[ijet] > 15){
                        if(jet_veto_map[thisEffTable]->GetBinContent(jet_veto_map[thisEffTable]->FindBin((*jetEta)[ijet],(*jetPhi)[ijet]))>0){
                            jetvetoBool=1;
                            break;
                        } 
                    //}

                }

                if(jetvetoBool) continue;

            }

            //=================ENTERING JET LOOP==================
            for(int ijet=0; ijet < jetCounter; ijet++){
                hjet_eta_phi_after_veto->Fill((*jetEta)[ijet],(*jetPhi)[ijet]);
                
                long int NNtrk = (dau_pt->at(ijet)).size();
                gRandom->SetSeed(0);
                double eta_smear;
                eta_smear=0;
                if( fabs(((*jetEta)[ijet])+eta_smear) > jetEtaCut ) continue;
                if( (*jetPt)[ijet] < jetPtCut_Jet   ) continue;
                hJetPt->Fill((*jetPt)[ijet]);
                // filling distrivutions within track bins
                // ALSO VERY IMPORTANLTY changing the tkBool to 1 for this particular jet. This will be usefull later wen I create conditons for filling other historgams.

                double jet_HLT_weight = 1.0;
                /*
                if((*jetPt)[ijet] < 880 && (*jetPt)[ijet] > 550){
                    jet_HLT_weight = 1.0/(hdid500->GetBinContent(hdid500->FindBin((*jetPt)[ijet]) ) );
                }
                */
                int tkBool[trackbin] = {0};

                int n_ChargeMult_DCA_labPt_Eta_exclusion =0;
                double n_ChargeMult_DCA_labPt_Eta_exclusion_Cor =0;

                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    if((*dau_chg)[ijet][A_trk] == 0) continue;//charge
                    if(fabs((*dau_pt)[ijet][A_trk])  < 0.3)     continue;//lab pt
                    if(fabs((*dau_eta)[ijet][A_trk]) > 2.4)     continue;//lab eta

                    double dauptnow = (*dau_pt)[ijet][A_trk]; double pterr_R_pt =0; double dcaXY = 0; double dcaZ = 0;
                    pterr_R_pt= ( (*dau_ptError)[ijet][A_trk] ) / ( (*dau_pt)[ijet][A_trk] );
                    dcaXY = (*dau_XYDCAsig)[ijet][A_trk];
                    dcaZ = (*dau_ZDCAsig)[ijet][A_trk];

                    if(pterr_R_pt   > 0.1   && dauptnow > 0.5)  continue;
                    if(fabs(dcaXY)  > 3     && dauptnow > 0.5)  continue;
                    if(fabs(dcaZ)   > 3     && dauptnow > 0.5)  continue;
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
                        hJet_Pass                   ->Fill(i);
                        if((*jetPt)[ijet] >= jetPtCut_Event) hJet_Pass550                           ->Fill(i);
                        if((*jetPt)[ijet] >= jetPtCut_Event) hJet_Pass550_hltCor                   ->Fill(i,1.0*jet_HLT_weight);

                        hBinDist_cor[i]             ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor,1.0*jet_HLT_weight);
                        hBinDist_unc[i]             ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion,1.0*jet_HLT_weight);
                    }
                }

                int Ntrig[trackbin][ptbin_T] = {0};
                double NtrigCorrected[trackbin][ptbin_T] = {0};
                int A_ptBool[NNtrk][ptbin_A] = {0};
                int T_ptBool[NNtrk][ptbin_T] = {0};

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

                    for(int i = 0; i < ptbin_T; i++){
                        if(jet_dau_pt >= pt_T_lo[i] && jet_dau_pt < pt_T_hi[i]){
                            T_ptBool[A_trk][i] = 1;
                        }
                    }
                    for(int i = 0; i < ptbin_A; i++){
                        if(jet_dau_pt >= pt_A_lo[i] && jet_dau_pt < pt_A_hi[i]){
                            A_ptBool[A_trk][i] = 1;
                        }
                    }

                    //getting et pt efficiency for A track in beam frame
                    double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    //double Atrk_weight=1.0; 

                    //in the case where the total jet mult and the individual daughter pt is acceptable for this track bin and pt bin, we increase the Ntrig count.
                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin_T; j++){
                            if(tkBool[i] + T_ptBool[A_trk][j] == 2){
                                NtrigCorrected[i][j] += (1.0/Atrk_weight);
                                Ntrig[i][j] += 1;
                                //cout<<"itrack="<<i<<" ipT="<<j<<" NtrigCorrected="<<NtrigCorrected[i][j]<<endl;
                                //cout<<"itrack="<<i<<" ipT="<<j<<" Ntrig="<<Ntrig[i][j]<<endl;
                            }
                        }
                    }

                }//here ends the boolean array creations.
                //continuation of main loops. Here is where the 2D Corr plots are created using the above booleans and 
                for(int i = 0; i < trackbin; i++){
                    for(int j = 0; j < ptbin_T; j++){
                        hNtrig->Fill(i,j,Ntrig[i][j]);
                        hNtrigCorrected->Fill(i,j, NtrigCorrected[i][j]);
                        hAvg_NtrigCorrected_Bump->Fill(NtrigCorrected[i][j] - Ntrig[i][j]);
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
                    double phi_smear;
                        phi_smear=0;

                    double jet_dau_pt    =  ptWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_eta   = etaWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_theta = 2*ATan(Exp(-(jet_dau_eta)));
                    
                    //getting et pt efficiency for A track in beam frame
                    //double Atrk_weight = 1.0;//(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    double Atrk_weight = (hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    hAvg_Atrk_Weight->Fill(1.0/(Atrk_weight));
                    
                    if(jet_dau_eta > track_eta_lim) continue;
                    
                    h_jet_jT->Fill(jet_dau_pt);
                    h_jet_etastar->Fill(jet_dau_eta);
                    h_jet_cor_jT->Fill(jet_dau_pt,1.0/(Atrk_weight));
                    h_jet_cor_etastar->Fill(jet_dau_eta,1.0/(Atrk_weight));

                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin_A; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                int k_PU=0;
                                hEPDrawCor_A[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight/Atrk_weight);
                                hEPDrawUnc_A[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight);
                            }
                        }
                        for(int j = 0; j < ptbin_T; j++){
                            if(tkBool[i] + T_ptBool[A_trk][j] == 2){
                                int k_PU=0;
                                hEPDrawCor_T[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight/Atrk_weight);
                                hEPDrawUnc_T[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight);
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

                        double T_jet_dau_pt    =  ptWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);
                        double T_jet_dau_eta   = etaWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);
                        double T_jet_dau_phi   = phiWRTJet( (double)(*jetPt)[ijet], (double)(*jetEta)[ijet] +eta_smear    , (double)(*jetPhi)[ijet] +phi_smear      , (double)(*dau_pt)[ijet][T_trk], (double)(*dau_eta)[ijet][T_trk], (double)(*dau_phi)[ijet][T_trk]);


                        if(T_jet_dau_eta > track_eta_lim) continue;

                        double deltaEta = (jet_dau_eta - T_jet_dau_eta);
                        double deltaPhi = (TMath::ACos(TMath::Cos(jet_dau_phi - T_jet_dau_phi)));
                        
                        //getting et pt efficiency for T track in beam frame
                        //double Ttrk_weight =1.0; //(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][T_trk] , (*dau_eta)[ijet][T_trk] )));
                        double Ttrk_weight =(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][T_trk] , (*dau_eta)[ijet][T_trk] )));

                        for(        int i = 0; i < trackbin; i++){
                            for(    int j = 0; j < ptbin_T;    j++){
                                for( int k=0;k<ptbin_A;k++){
                                    if((tkBool[i] + A_ptBool[A_trk][k] + T_ptBool[T_trk][j] == 3)
                                    ||(tkBool[i] + A_ptBool[T_trk][k] + T_ptBool[A_trk][j] == 3)){

                                        //cout << "1.0/(Atrk_weight*Ttrk_weight*NtrigCorrected[i][j]) " << 1.0/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ) << " and components: " << Atrk_weight << " " << Ttrk_weight << " " << NtrigCorrected[i][j] << endl;

                                        hPairs[k]->Fill(i,j);

                                        if(Atrk_weight == 0){/* cout  << "Atrk_weight = 0!" << endl;
                                        */
                                            continue;}
                                        if(Ttrk_weight == 0){/* cout << "Ttrk_weight = 0!" << endl;
                                        */
                                            continue;}
                                        if(NtrigCorrected[i][j] == 0){ cout << "NtrigCorrected[i][j] = 0!" << endl; continue;}

                                        hTotalWeight->Fill(1.0 /(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));

                                        int k_PU=0;
                                        hSignalShiftedUnc[i][j][k][k_PU]->Fill(deltaEta, deltaPhi,                  1.0*jet_HLT_weight/(Ntrig[i][j] ));
                                        hSignalShiftedUnc[i][j][k][k_PU]->Fill(-deltaEta, deltaPhi,                 1.0*jet_HLT_weight/(Ntrig[i][j] ));
                                        hSignalShiftedUnc[i][j][k][k_PU]->Fill(deltaEta, -deltaPhi,                 1.0*jet_HLT_weight/(Ntrig[i][j] ));
                                        hSignalShiftedUnc[i][j][k][k_PU]->Fill(-deltaEta, -deltaPhi,                1.0*jet_HLT_weight/(Ntrig[i][j] ));
                                        hSignalShiftedUnc[i][j][k][k_PU]->Fill( deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Ntrig[i][j] ));
                                        hSignalShiftedUnc[i][j][k][k_PU]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Ntrig[i][j] ));
                                        
                                        hSignalShiftedCor[i][j][k][k_PU]->Fill(deltaEta, deltaPhi,                  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                        hSignalShiftedCor[i][j][k][k_PU]->Fill(-deltaEta, deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                        hSignalShiftedCor[i][j][k][k_PU]->Fill(deltaEta, -deltaPhi,                 1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                        hSignalShiftedCor[i][j][k][k_PU]->Fill(-deltaEta, -deltaPhi,                1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                        hSignalShiftedCor[i][j][k][k_PU]->Fill( deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                        hSignalShiftedCor[i][j][k][k_PU]->Fill(-deltaEta,2*TMath::Pi() - deltaPhi,  1.0*jet_HLT_weight/(Atrk_weight * Ttrk_weight * NtrigCorrected[i][j] ));
                                    }
                                }
                            }
                        }
                    }//T_trk end
                }//A_trk end
            }//ijet end
        }//ievent end
        
        fFile->Close();
    }//f in filelist end 

    int backMult =10;
    for(int wtrk = 1; wtrk < trackbin+1; wtrk++){
        for(int wppt = 1; wppt < ptbin_T+1; wppt++){ 
            for(int wpptA = 1; wpptA < ptbin_A+1; wpptA++){ 
                for(int wpPU = 1; wpPU < PUbin+1; wpPU++){
                    //Nent is the number of pairs in the signal which we will try to 10x
                    //Xent is the number of pseudoparticles requried such that when we build the pairs nCp = Xent CHOOSE 2 will give 
                    //us 10 times as many pairs as we have in the signal histogrm.

                    //long int NENT =  hPairsPU[wpPU-1]->GetBinContent(wtrk, wppt);
                    long int NENT =  hPairs[wpptA-1]->GetBinContent(wtrk, wppt);
                    //long int XENT =  ((1+floor(sqrt(1+(4*2*backMult*NENT))))/2);
                    long int XENT =  floor(sqrt(backMult*NENT));
                    float A_ETA[XENT] = {0};
                    float A_PHI[XENT] = {0};
                    float A_ETA_Cor[XENT] = {0};
                    float A_PHI_Cor[XENT] = {0};
                    float T_ETA[XENT] = {0};
                    float T_PHI[XENT] = {0};
                    float T_ETA_Cor[XENT] = {0};
                    float T_PHI_Cor[XENT] = {0};
                    
                    //cout<<"wtrk="<<wtrk<<" wppt="<<wppt<<" wpptA="<<wpptA<<endl;

                    for(int x = 0; x<XENT; x++){
                        gRandom->SetSeed(0);

                        double WEta1, WPhi1;//making the pseudoparticles
                        hEPDrawUnc_A[wtrk-1][wpptA-1][wpPU-1]->GetRandom2(WEta1, WPhi1);
                        A_ETA[x] = WEta1;
                        A_PHI[x] = WPhi1;
                        //cout<<"A_eta="<<A_ETA[x]<<" A_phi"<<A_PHI[x]<<endl;

                        double WEta1_Cor, WPhi1_Cor;//making the pseudoparticles
                        hEPDrawCor_A[wtrk-1][wpptA-1][wpPU-1]->GetRandom2(WEta1_Cor, WPhi1_Cor);
                        A_ETA_Cor[x] = WEta1_Cor;
                        A_PHI_Cor[x] = WPhi1_Cor;
                        //cout<<"A_eta_cor="<<A_ETA_Cor[x]<<" A_phi_cor"<<A_PHI_Cor[x]<<endl;

                        double WEta2, WPhi2;//making the pseudoparticles
                        hEPDrawUnc_T[wtrk-1][wppt-1][wpPU-1]->GetRandom2(WEta2, WPhi2);
                        T_ETA[x] = WEta2;
                        T_PHI[x] = WPhi2;
                        //cout<<"T_eta="<<T_ETA[x]<<" T_phi"<<T_PHI[x]<<endl;
                        
                        double WEta2_Cor, WPhi2_Cor;//making the pseudoparticles
                        hEPDrawCor_T[wtrk-1][wppt-1][wpPU-1]->GetRandom2(WEta2_Cor, WPhi2_Cor);
                        T_ETA_Cor[x] = WEta2_Cor;
                        T_PHI_Cor[x] = WPhi2_Cor;
                        //cout<<"T_eta_cor="<<T_ETA_Cor[x]<<" T_phi_cor"<<T_PHI_Cor[x]<<endl;
                    }
                    
                    for(long int i = 0; i < XENT; i++){
                        for(long int j = 0; j < XENT; j++){
                            double WdeltaEta = (A_ETA[i]-T_ETA[j]);
                            double WdeltaPhi = (TMath::ACos(TMath::Cos(A_PHI[i]-T_PHI[j])));
                            hBckrndShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(WdeltaEta, WdeltaPhi,   1);//./XENT);
                            hBckrndShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(-WdeltaEta, WdeltaPhi,  1);//./XENT);
                            hBckrndShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(WdeltaEta, -WdeltaPhi,  1);//./XENT);
                            hBckrndShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(-WdeltaEta, -WdeltaPhi, 1);//./XENT);
                            hBckrndShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(WdeltaEta, 2*TMath::Pi() - WdeltaPhi, 1);//./XENT);
                            hBckrndShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(-WdeltaEta,2*TMath::Pi() - WdeltaPhi, 1);//./XENT);

                            double WdeltaEta_Cor = (A_ETA_Cor[i]-T_ETA_Cor[j]);
                            double WdeltaPhi_Cor = (TMath::ACos(TMath::Cos(A_PHI_Cor[i]-T_PHI_Cor[j])));
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(WdeltaEta_Cor, WdeltaPhi_Cor,   1);//./XENT);
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(-WdeltaEta_Cor, WdeltaPhi_Cor,  1);//./XENT);
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(WdeltaEta_Cor, -WdeltaPhi_Cor,  1);//./XENT);
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(-WdeltaEta_Cor, -WdeltaPhi_Cor, 1);//./XENT);
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(WdeltaEta_Cor, 2*TMath::Pi() - WdeltaPhi_Cor, 1);//./XENT);
                            hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Fill(-WdeltaEta_Cor,2*TMath::Pi() - WdeltaPhi_Cor, 1);//./XENT);

                        }
                    }
                }
            }
        }
    }

    string subList = fList.substr(fList.size() - 3);
    TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/xiaoyul/Run3_root_out/corrected/jTDiff/job_%s.root",subList.c_str()), "recreate");
    for(int wtrk =1; wtrk <trackbin+1; wtrk++){
        hBinDist_cor[wtrk-1]->Write();
        hBinDist_unc[wtrk-1]->Write();
        for(int wpPU =1; wpPU<PUbin+1; wpPU++){
            for(int wppt =1; wppt <ptbin_T+1; wppt++){
                for(int wpptA =1; wpptA <ptbin_A+1; wpptA++){
                    /*
                    hSignalShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Write(Form("hSigS_Unc_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1][wpptA-1]),(int)(10*ptbinbounds_hi[wppt-1][wpptA-1]),wpPU ));
                    hBckrndShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Write(Form("hBckS_Unc_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1][wpptA-1]),(int)(10*ptbinbounds_hi[wppt-1][wpptA-1]),wpPU ));
                    hSignalShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Write(Form("hSigS_Cor_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1][wpptA-1]),(int)(10*ptbinbounds_hi[wppt-1][wpptA-1]),wpPU ));
                    hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Write(Form("hBckS_Cor_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1] ,(int)(10*ptbinbounds_lo[wppt-1][wpptA-1]),(int)(10*ptbinbounds_hi[wppt-1][wpptA-1]),wpPU ));
                    */
                    hSignalShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Write();
                    hBckrndShiftedUnc[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Write();
                    hSignalShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Write();
                    hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Write();
                }
            }
            for(int wppt=1;wppt<ptbin_T+1;wppt++){
                hEPDrawUnc_T[wtrk-1][wppt-1][wpPU-1]->Write();
                hEPDrawCor_T[wtrk-1][wppt-1][wpPU-1]->Write();
            }
            for(int wpptA=1;wpptA<ptbin_A+1;wpptA++){
                hEPDrawUnc_A[wtrk-1][wpptA-1][wpPU-1]->Write();
                hEPDrawCor_A[wtrk-1][wpptA-1][wpPU-1]->Write();
            }
        }
    }
    for(int i=0;i<ptbin_A;i++){
        hPairs[i]->Write();
    }
    hNtrig->Write();
    hNtrigCorrected->Write();
    hTotalWeight->Write();
    hAvg_Atrk_Weight->Write();
    hAvg_NtrigCorrected_Bump->Write();
    hJetPt->Write();
    hJetPt_wo_ptcut->Write(); 
    hEvent_Pass->Write();
    hJet_Pass->Write();
    hJet_Pass550->Write();
    hJet_Pass550_hltCor->Write();
    hBinDist_cor_single->Write();
    hBinDist_unc_single->Write();
    
    h_jet_jT->Write();
    h_jet_etastar->Write();
    h_jet_cor_jT->Write();
    h_jet_cor_etastar->Write();

    h_lab_JetMult_pT->Write();
    h_lab_JetMult_phi->Write();
    h_lab_JetMult_eta->Write();
    h_lab_cor_JetMult_pT->Write();
    h_lab_cor_JetMult_phi->Write();
    h_lab_cor_JetMult_eta->Write();

    hjet_eta_phi_before_veto->Write();
    hjet_eta_phi_after_veto->Write();

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

