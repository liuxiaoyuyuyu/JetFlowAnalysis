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


void MyClass::Loop(int job, std::string fList){

    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);

    //Initializing Histograms
    TH1D* hEvent_Pass   = new TH1D("hEvent_Pass","hEvent_Pass", trackbin,bin0,trackbin);
    TH2D* hRun_Lumi = new TH2D("hRun_Lumi","hRun_Lumi",6200,366300.5,372500.5,10000,-0.5,9999.5);

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

        //========ENTERING EVENT LOOP========
        for (Long64_t ievent=0; ievent <nentries; ievent ++){
            Long64_t jevent = LoadTree(ievent);
            nb = fChain->GetEntry(ievent);   nbytes += nb;

            hRun_Lumi->Fill(nRun,nLumi); 
            //cut on jetPt and jetN (number of [HLT-passed] jets in an evnet)
            if(!F_eventpass(jetPt, jetN, 0)){
                continue;
            }

            hEvent_Pass->Fill(1);

        }
        fFile->Close();

        string subList = fList.substr(fList.size() - 3);
        TFile* fS_tempA = new TFile(Form("/eos/cms/store/group/phys_heavyions/flowcorr/root_out_qa/statistics/job_stat_qa_%s_%d.root",subList.c_str(),f), "recreate");
        hEvent_Pass->Write();
        hRun_Lumi->Write(); 
        fS_tempA->Close();
        
        hEvent_Pass->Close();
        hRun_Lumi->Close(); 
    
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

