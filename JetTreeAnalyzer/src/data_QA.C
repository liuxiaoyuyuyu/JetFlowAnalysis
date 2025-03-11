#define MyClass_cxx

//#include "include/MyTrim.h"
#include "include/MyTrimMC.h"
#include "include/coordinateTools.h"
//#include "include/1d2d_constants.h"
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


void MyClass::Loop(int job, std::string fList){

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    //Initializing Histograms
    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);

    TH1D* h_jet_jT=new TH1D("jet_jT","jet_jT",200,0,20);
    TH1D* h_jet_etastar=new TH1D("jet_etastar","jet_etastar",100,0,10);
    TH1D* hBinDist_cor_single = new TH1D("hBinDist_cor_single","hBinDist_cor_single",bin120,bin0,bin120);


    // MAIN CODE BEGINS
    std::cout << "Starting event loop" << std::endl;
    std::cout << "Total Number of Files in this Job: " << fileList.size() << std::endl;
    for(int f = 0; f<fileList.size(); f++){
        fFile = TFile::Open(fileList.at(f).c_str(),"read");
        TTree *tree = (TTree*)fFile->Get("analyzerOffline/trackTree");
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

        //========ENTERING EVENT LOOP========
        for (Long64_t ievent=0; ievent <nentries; ievent ++){
            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;

            //cut on jetPt and jetN (number of [HLT-passed] jets in an evnet)
            //if(!F_eventpass(jetPt, jetN, jetPtCut_Event)){
            //    continue;
            //}

            hEvent_Pass->Fill(1);

            int jetCounter = genJetPt->size();
            if(jetCounter == 0) continue;

            //========ENTERING JET LOOP========
            for(int ijet=0; ijet < jetCounter; ijet++){

                //long int NNtrk = (dau_pt->at(ijet)).size();
                long int NNtrk = (genDau_pt->at(ijet)).size();
                if( fabs(((*genJetEta)[ijet])) > jetEtaCut ) continue;
                if( (*genJetPt)[ijet] < jetPtCut_Jet   ) continue;

                int n_ChargeMult_DCA_labPt_Eta_exclusion =0;
                
                for(int  A_trk=0; A_trk < NNtrk; A_trk++ ){
                    /*
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

                    if(pterr_R_pt   > 0.1   && dauptnow > 0.5)  continue;
                    if(fabs(dcaXY)  > 3     && dauptnow > 0.5)  continue;
                    if(fabs(dcaZ)   > 3     && dauptnow > 0.5)  continue;
                    n_ChargeMult_DCA_labPt_Eta_exclusion += 1;

                    double jet_dau_pt    =  ptWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] , (double)(*jetPhi)[ijet] , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_eta   = etaWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] , (double)(*jetPhi)[ijet] , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet((double)(*jetPt)[ijet], (double)(*jetEta)[ijet] , (double)(*jetPhi)[ijet] , (double)(*dau_pt)[ijet][A_trk], (double)(*dau_eta)[ijet][A_trk], (double)(*dau_phi)[ijet][A_trk]);
                    */

                    if((*genDau_chg)[ijet][A_trk] == 0) continue;//charge
                    if(fabs((*genDau_pt)[ijet][A_trk])  < 0.3)     continue;//lab pt
                    if(fabs((*genDau_eta)[ijet][A_trk]) > 2.4)     continue;//lab eta
                    n_ChargeMult_DCA_labPt_Eta_exclusion += 1;
                    
                    double jet_dau_pt = ptWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet], (double)(*genJetPhi)[ijet], (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    double jet_dau_eta   = etaWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet]   , (double)(*genJetPhi)[ijet]  , (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    double jet_dau_phi   = phiWRTJet((double)(*genJetPt)[ijet], (double)(*genJetEta)[ijet]  , (double)(*genJetPhi)[ijet]   , (double)(*genDau_pt)[ijet][A_trk], (double)(*genDau_eta)[ijet][A_trk], (double)(*genDau_phi)[ijet][A_trk]);
                    
                    h_jet_jT->Fill(jet_dau_pt);
                    h_jet_etastar->Fill(jet_dau_eta);
                }
                hBinDist_cor_single            ->Fill(n_ChargeMult_DCA_labPt_Eta_exclusion);
            
            }
        }
        fFile->Close();
    
    }//end looping over files
    
    string subList = fList.substr(fList.size() - 3);
    //TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/flowcorr/root_out_qa/hlt260/job_%s_%d.root",subList.c_str(),f), "recreate");
    TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/flowcorr/root_out_qa/hlt500/job_%s_%d.root",subList.c_str(),f), "recreate");
    hEvent_Pass->Write();
    hBinDist_cor_single->Write();
    h_jet_jT->Write();
    h_jet_etastar->Write();
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

