//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 17 21:43:29 2020 by ROOT version 6.12/07
// from TTree trackTree/v1
// found on file: outputMiniAOD.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <string>

using namespace std;

class MyClass {
public :
   TFile * fFile;
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   std::vector<std::string> fileList;

   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t           nRun;
   ULong64_t           nEv;
   UInt_t           nLumi;
   vector<int>     *jetNumDaughters;
   vector<float>   *jetEta;
   vector<float>   *jetPt;
   vector<float>   *jetPhi;
   vector<int>     *chargedMultiplicity;
   Int_t           jetN;
   //vector<vector<float> > *dau_PuppiW;
   vector<vector<int> > *dau_chg;
   vector<vector<float> > *dau_pt;
   vector<vector<float> > *dau_ptError;
   vector<vector<float> > *dau_eta;
   vector<vector<float> > *dau_phi;
   vector<vector<float> > *dau_theta;
   vector<vector<float> > *dau_XYDCAsig;
   vector<vector<float> > *dau_ZDCAsig;
   //vector<vector<float> > *dau_vz;
   //vector<vector<float> > *dau_vy;
   //vector<vector<float> > *dau_vx;

   // List of branches
   TBranch        *b_nRun;   //!
   TBranch        *b_nEv;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_jetNumDaughters;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_chargedMultiplicity;   //!
   TBranch        *b_jetN;   //!
   //TBranch        *b_dau_PuppiW;   //!
   TBranch        *b_dau_chg;   //!
   TBranch        *b_dau_pt;   //!
   TBranch        *b_dau_ptError;   //!
   TBranch        *b_dau_eta;   //!
   TBranch        *b_dau_phi;   //!
   TBranch        *b_dau_theta;   //!
   TBranch        *b_dau_XYDCAsig;   //!
   TBranch        *b_dau_ZDCAsig;   //!
   //TBranch        *b_dau_vz;   //!
   //TBranch        *b_dau_vy;   //!
   //TBranch        *b_dau_vx;   //!

   MyClass(std::vector<std::string> _fileList);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int job, std::string fList);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

MyClass::MyClass(std::vector<std::string> _fileList) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("total_output_ak8.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("total_output_ak8.root");
//      }
//      TDirectory * dir = (TDirectory*)f->Get("total_output_ak8.root:/analyzer");
//      dir->GetObject("trackTree",tree);
//
//   }
    fileList = _fileList;
}

MyClass::~MyClass()
{
   //if (!fChain) return;
   //delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jetNumDaughters = 0;
   jetEta = 0;
   jetPt = 0;
   jetPhi = 0;
   chargedMultiplicity = 0;
   //dau_PuppiW = 0;
   dau_chg = 0;
   dau_pt = 0;
   dau_ptError = 0;
   dau_eta = 0;
   dau_phi = 0;
   dau_theta = 0;
   dau_XYDCAsig = 0;
   dau_ZDCAsig = 0;
   //dau_vz = 0;
   //dau_vy = 0;
   //dau_vx = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("nEv", &nEv, &b_nEv);
   fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   fChain->SetBranchAddress("jetNumDaughters", &jetNumDaughters, &b_jetNumDaughters);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("chargedMultiplicity", &chargedMultiplicity, &b_chargedMultiplicity);
   fChain->SetBranchAddress("jetN", &jetN, &b_jetN);
   //fChain->SetBranchAddress("dau_PuppiW", &dau_PuppiW, &b_dau_PuppiW);
   fChain->SetBranchAddress("dau_chg", &dau_chg, &b_dau_chg);
   fChain->SetBranchAddress("dau_pt", &dau_pt, &b_dau_pt);
   fChain->SetBranchAddress("dau_ptError", &dau_ptError, &b_dau_ptError);
   fChain->SetBranchAddress("dau_eta", &dau_eta, &b_dau_eta);
   fChain->SetBranchAddress("dau_phi", &dau_phi, &b_dau_phi);
   fChain->SetBranchAddress("dau_theta", &dau_theta, &b_dau_theta);
   fChain->SetBranchAddress("dau_XYDCAsig", &dau_XYDCAsig, &b_dau_XYDCAsig);
   fChain->SetBranchAddress("dau_ZDCAsig", &dau_ZDCAsig, &b_dau_ZDCAsig);
   //fChain->SetBranchAddress("dau_vz", &dau_vz, &b_dau_vz);
   //fChain->SetBranchAddress("dau_vy", &dau_vy, &b_dau_vy);
   //fChain->SetBranchAddress("dau_vx", &dau_vx, &b_dau_vx);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
