#define MyClass_cxx

#include "include/MyTrimMC.h"
#include "include/coordinateTools.h"
#include "include/mc_constants_jT.h"

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
    TH2D* hNtrigCorrected        = new TH2D("hNtrigCorrected","hNtrigCorrected",  trackbin, bin0,trackbin, ptbin_T, bin0, ptbin_T);

    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);
    TH1D* hJet_Pass     = new TH1D("hJet_Pass"  ,"hJet_Pass"  , trackbin,bin0,trackbin);
    TH1D* hJet_Pass550     = new TH1D("hJet_Pass550"  ,"hJet_Pass550"  , trackbin,bin0,trackbin);
    TH1D* hJetPt = new TH1D("hJetPt"  ,"hJetPt" ,i340,i100,i3500);
    TH1D* hJetPt_wo_ptcut = new TH1D("hJetPt_wo_ptcut"  ,"hJetPt_wo_ptcut" ,i340,i100,i3500);
    TH1D* hBinDist_cor_single = new TH1D("hBinDist_cor_single","hBinDist_cor_single",bin360,bin0,bin120);

    TH2D* h_lab_cor_JetMult_pT=new TH2D("lab_cor_JetMult_pT","lab_cor_JetMult_pT",340,100,3500,360,0.5,120.5);
    TH2D* h_lab_cor_JetMult_phi=new TH2D("lab_cor_JetMult_phi","lab_cor_JetMult_phi",30,-TMath::Pi(),TMath::Pi(),360,0.5,120.5);
    TH2D* h_lab_cor_JetMult_eta=new TH2D("lab_cor_JetMult_eta","lab_cor_JetMult_eta",34,-1.7,1.7,360,0.5,120.5);
    
    TH1D* h_jet_cor_jT=new TH1D("jet_cor_jT","jet_cor_jT",200,0,20);
    TH1D* h_jet_cor_etastar=new TH1D("jet_cor_etastar","jet_cor_etastar",100,0,10);

    TH2D* hEPDrawCor_T[trackbin][ptbin_T][PUbin];
    TH2D* hEPDrawCor_A[trackbin][ptbin_A][PUbin];
    TH2D* hSignalShiftedCor[trackbin][ptbin_T][ptbin_A][PUbin];
    TH2D* hBckrndShiftedCor[trackbin][ptbin_T][ptbin_A][PUbin];
    
    //For recording Nch dis. and calculating <Nch> per Nch bin.
    TH1D* hBinDist_cor[trackbin];

    for(int wtrk = 1; wtrk<trackbin+1; wtrk++){
        hBinDist_cor[wtrk-1]    = new TH1D(Form("hBinDist_cor_%d",wtrk),Form("hBinDist_cor_%d",wtrk), bin360, bin0, bin120);
        for(int wpPU = 1; wpPU<PUbin+1; wpPU++){
            for(int wppt = 1; wppt<ptbin_T+1; wppt++){
                for(int wpptA = 1; wpptA<ptbin_A+1; wpptA++){
                    hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]   = new TH2D(Form("hBkg_Cor_trk_%d_ptT_%d_ptA_%d_PU_%d",wtrk,wppt,wpptA,wpPU) ,
                    Form("hBkg_Cor_%d_to_%d_and_T_%d_to_%d_A_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_T_lo[wppt-1]),(int)(10*pt_T_hi[wppt-1]),(int)(10*pt_A_lo[wpptA-1]),(int)(10*pt_A_hi[wpptA-1]),wpPU),
                    41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);
                    
                    hSignalShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]   = new TH2D(Form("hSig_Cor_trk_%d_ptT_%d_ptA_%d_PU_%d",wtrk,wppt,wpptA,wpPU) ,
                    Form("hSig_Cor_%d_to_%d_and_T_%d_to_%d_A_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_T_lo[wppt-1]),(int)(10*pt_T_hi[wppt-1]),(int)(10*pt_A_lo[wpptA-1]),(int)(10*pt_A_hi[wpptA-1]),wpPU) ,
                    41,-(20*EtaBW)-(0.5*EtaBW),(20*EtaBW)+(0.5*EtaBW),33,-(8*PhiBW)-0.5*PhiBW,(24*PhiBW)+0.5*PhiBW);

                }
            }
            
            for(int wppt=1;wppt<ptbin_T+1;wppt++){
                hEPDrawCor_T[wtrk-1][wppt-1][wpPU-1]          = new TH2D(Form( "hEPDraw_Cor_T_trk_%d_ppt_%d_PU_%d",wtrk,wppt,wpPU) ,Form( "hEPDraw_Cor_T_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_T_lo[wppt-1]),(int)(10*pt_T_hi[wppt-1]),wpPU) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
            }
            
            for(int wpptA=1;wpptA<ptbin_A+1;wpptA++){
                hEPDrawCor_A[wtrk-1][wpptA-1][wpPU-1]          = new TH2D(Form( "hEPDraw_Cor_A_trk_%d_ppt_%d_PU_%d",wtrk,wpptA,wpPU) ,Form( "hEPDraw_Cor_A_%d_to_%d_and_%d_to_%d_w_PU_%d",trackbinbounds[wtrk-1],trackbinboundsUpper[wtrk-1],(int)(10*pt_A_lo[wpptA-1]),(int)(10*pt_A_hi[wpptA-1]),wpPU) , EPD_xb   , EPD_xlo, EPD_xhi , EPD_yb      , EPD_ylo    , EPD_yhi);
            }
        }
    }

    std::cout << "Starting event loop" << std::endl;
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;
    for(int f = 0; f<fileList.size(); f++){
        int f_from_file = f;
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("trackTree");
        Init(tree);

        std::cout << "File " << f+1 << " out of " << fileList.size() << std::endl;
        Long64_t nbytes = 0, nb = 0;
        Long64_t nentries = fChain->GetEntriesFast();
        cout<<"Total Entries is:"<<endl;
        cout<< nentries <<endl;

        // ENTERING EVENT LOOP
        for (Long64_t ievent=0; ievent <nentries; ievent ++){
            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;

            if(!F_eventpass(genJetPt, genJetPt->size(), jetPtCut_Event)){
                continue;
            }
            
            hEvent_Pass->Fill(1);

            int jetCounter = genJetPt->size();
            if(jetCounter == 0) continue;

            //=================ENTERING JET LOOP==================
            for(int ijet=0; ijet < jetCounter; ijet++){
                
                long int NNtrk = (genDau_pt->at(ijet)).size();
                gRandom->SetSeed(0);
                double eta_smear;
                eta_smear=0;
                if( fabs(((*genJetEta)[ijet])+eta_smear) > jetEtaCut ) continue;
                if( (*genJetPt)[ijet] < jetPtCut_Jet   ) continue;
                hJetPt->Fill((*genJetPt)[ijet]);
                // filling distrivutions within track bins
                // ALSO VERY IMPORTANLTY changing the tkBool to 1 for this particular jet. This will be usefull later wen I create conditons for filling other historgams.

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
                        hJet_Pass                   ->Fill(i);
                        if((*genJetPt)[ijet] >= jetPtCut_Event) hJet_Pass550                           ->Fill(i);

                        hBinDist_cor[i]             ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion_Cor,1.0*jet_HLT_weight);
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
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;
                    if((*genDau_pt)[ijet][A_trk] < 0.3) continue;
                    if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4) continue;

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
                    double Atrk_weight=1.0; 

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
                    }
                }

                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    //this should be redundant if it passes the bools above? i guess it helps skip daughters faster. maybe i can reindex and run through the daughters quickly by aranging all the charged dauhghter sat the front.
                    if((*genDau_chg)[ijet][A_trk] == 0) continue;
                    if((*genDau_pt)[ijet][A_trk] < 0.3) continue;
                    if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4) continue;
                    //if((*highPurity)[ijet][A_trk] == 0) continue;

                    gRandom->SetSeed(0);
                    double phi_smear;
                    phi_smear=0;

                    double jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet] +eta_smear    , (double)(*genJetPhi)[ijet] +phi_smear      , (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    double jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet] +eta_smear    , (double)(*genJetPhi)[ijet] +phi_smear      , (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet] +eta_smear    , (double)(*genJetPhi)[ijet] +phi_smear      , (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    double jet_dau_theta = 2*ATan(Exp(-(jet_dau_eta)));
                    
                    //getting et pt efficiency for A track in beam frame
                    double Atrk_weight = 1.0;//(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][A_trk] , (*dau_eta)[ijet][A_trk] )));
                    
                    if(jet_dau_eta > track_eta_lim) continue;
                    
                    h_jet_cor_jT->Fill(jet_dau_pt,1.0/(Atrk_weight));
                    h_jet_cor_etastar->Fill(jet_dau_eta,1.0/(Atrk_weight));

                    for(int i = 0; i < trackbin; i++){
                        for(int j = 0; j < ptbin_A; j++){
                            if(tkBool[i] + A_ptBool[A_trk][j] == 2){
                                int k_PU=0;
                                hEPDrawCor_A[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight/Atrk_weight);
                            }
                        }
                        for(int j = 0; j < ptbin_T; j++){
                            if(tkBool[i] + T_ptBool[A_trk][j] == 2){
                                int k_PU=0;
                                hEPDrawCor_T[i][j][k_PU]->Fill(jet_dau_eta, jet_dau_phi, 1.0*jet_HLT_weight/Atrk_weight);
                            }
                        }
                    }
                    if(A_trk == NNtrk - 1) continue;
                    for(long int T_trk=A_trk+1; T_trk< NNtrk; T_trk++ ){

                        if((*genDau_chg)[ijet][T_trk] == 0) continue;
                        if((*genDau_pt)[ijet][T_trk] < 0.3) continue;
                        if(fabs((*genDau_eta)[ijet][T_trk]) > 2.4) continue;

                        double T_jet_dau_pt    =  ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet] +eta_smear    , (double)(*genJetPhi)[ijet] +phi_smear      , (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                        double T_jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet] +eta_smear    , (double)(*genJetPhi)[ijet] +phi_smear      , (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                        double T_jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet] +eta_smear    , (double)(*genJetPhi)[ijet] +phi_smear      , (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);


                        if(T_jet_dau_eta > track_eta_lim) continue;

                        double deltaEta = (jet_dau_eta - T_jet_dau_eta);
                        double deltaPhi = (TMath::ACos(TMath::Cos(jet_dau_phi - T_jet_dau_phi)));
                        
                        //getting et pt efficiency for T track in beam frame
                        double Ttrk_weight =1.0; //(hReco2D[thisEffTable]->GetBinContent(hReco2D[thisEffTable]->FindBin( (*dau_pt)[ijet][T_trk] , (*dau_eta)[ijet][T_trk] )));

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

                        double WEta1_Cor, WPhi1_Cor;//making the pseudoparticles
                        hEPDrawCor_A[wtrk-1][wpptA-1][wpPU-1]->GetRandom2(WEta1_Cor, WPhi1_Cor);
                        A_ETA_Cor[x] = WEta1_Cor;
                        A_PHI_Cor[x] = WPhi1_Cor;
                        //cout<<"A_eta_cor="<<A_ETA_Cor[x]<<" A_phi_cor"<<A_PHI_Cor[x]<<endl;

                        double WEta2_Cor, WPhi2_Cor;//making the pseudoparticles
                        hEPDrawCor_T[wtrk-1][wppt-1][wpPU-1]->GetRandom2(WEta2_Cor, WPhi2_Cor);
                        T_ETA_Cor[x] = WEta2_Cor;
                        T_PHI_Cor[x] = WPhi2_Cor;
                        //cout<<"T_eta_cor="<<T_ETA_Cor[x]<<" T_phi_cor"<<T_PHI_Cor[x]<<endl;
                    }
                    
                    for(long int i = 0; i < XENT; i++){
                        for(long int j = 0; j < XENT; j++){
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
    TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/xiaoyul/MC/13TeV/pythia8/jTDiff/job_%s.root",subList.c_str()), "recreate");
    for(int wtrk =1; wtrk <trackbin+1; wtrk++){
        hBinDist_cor[wtrk-1]->Write();
        for(int wpPU =1; wpPU<PUbin+1; wpPU++){
            for(int wppt =1; wppt <ptbin_T+1; wppt++){
                for(int wpptA =1; wpptA <ptbin_A+1; wpptA++){
                    hSignalShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Write();
                    hBckrndShiftedCor[wtrk-1][wppt-1][wpptA-1][wpPU-1]->Write();
                }
            }
            for(int wppt=1;wppt<ptbin_T+1;wppt++){
                hEPDrawCor_T[wtrk-1][wppt-1][wpPU-1]->Write();
            }
            for(int wpptA=1;wpptA<ptbin_A+1;wpptA++){
                hEPDrawCor_A[wtrk-1][wpptA-1][wpPU-1]->Write();
            }
        }
    }
    for(int i=0;i<ptbin_A;i++){
        hPairs[i]->Write();
    }
    hNtrig->Write();
    hNtrigCorrected->Write();
    hJetPt->Write();
    hJetPt_wo_ptcut->Write(); 
    hEvent_Pass->Write();
    hJet_Pass->Write();
    hBinDist_cor_single->Write();
    
    h_jet_cor_jT->Write();
    h_jet_cor_etastar->Write();

    h_lab_cor_JetMult_pT->Write();
    h_lab_cor_JetMult_phi->Write();
    h_lab_cor_JetMult_eta->Write();

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

