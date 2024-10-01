//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 30 23:40:55 2024 by ROOT version 6.26/11
// from TTree trackTree/v1
// found on file: MATCH_mc_test.root
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
   ULong64_t       nRun;
   UInt_t          nEv;
   UInt_t          nLumi;
   vector<float>   *xVtx;
   vector<float>   *yVtx;
   vector<float>   *zVtx;
   vector<int>     *nTracksVtx;
   vector<float>   *ptSumVtx;
   vector<int>     *jetNumDaughters;
   vector<float>   *jetEta;
   vector<float>   *jetPt;
   vector<float>   *jetPhi;
   vector<float>   *jetTheta;
   vector<float>   *jetMass;
   vector<int>     *muonMultiplicity;
   vector<int>     *chargedMultiplicity;
   Int_t           jetN;
   vector<float>   *dau_pt_sum;
   vector<vector<int> > *dau_chg;
   vector<vector<int> > *dau_pid;
   vector<vector<unsigned int> > *dau_vref;
   vector<vector<float> > *dau_pt;
   vector<vector<float> > *dau_ptError;
   vector<vector<float> > *dau_XYDCAsig;
   vector<vector<float> > *dau_ZDCAsig;
   vector<vector<float> > *dau_eta;
   vector<vector<float> > *dau_phi;
   vector<vector<float> > *dau_theta;
   vector<vector<float> > *dau_vz;
   vector<vector<float> > *dau_vy;
   vector<vector<float> > *dau_vx;
   vector<vector<float> > *dau_vrefz;
   vector<vector<float> > *dau_vrefy;
   vector<vector<float> > *dau_vrefx;
   vector<vector<float> > *dau_vp_difZ;
   vector<vector<float> > *dau_vp_difY;
   vector<vector<float> > *dau_vp_difX;
   Float_t         genQScale;
   Float_t         genWeight;
   Int_t           genSignalProcessID;
   vector<float>   *genJetEta;
   vector<float>   *genJetPt;
   vector<float>   *genJetPhi;
   vector<int>     *genJetChargedMultiplicity;
   vector<vector<int> > *genDau_chg;
   vector<vector<int> > *genDau_pid;
   vector<vector<float> > *genDau_pt;
   vector<vector<float> > *genDau_eta;
   vector<vector<float> > *genDau_phi;

   // List of branches
   TBranch        *b_nRun;   //!
   TBranch        *b_nRun;   //!
   TBranch        *b_nLumi;   //!
   TBranch        *b_xVtx;   //!
   TBranch        *b_yVtx;   //!
   TBranch        *b_zVtx;   //!
   TBranch        *b_nTracksVtx;   //!
   TBranch        *b_ptSumVtx;   //!
   TBranch        *b_jetNumDaughters;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetTheta;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_muonMultiplicity;   //!
   TBranch        *b_chargedMultiplicity;   //!
   TBranch        *b_jetN;   //!
   TBranch        *b_dau_pt_sum;   //!
   TBranch        *b_dau_chg;   //!
   TBranch        *b_dau_pid;   //!
   TBranch        *b_dau_vref;   //!
   TBranch        *b_dau_pt;   //!
   TBranch        *b_dau_ptError;   //!
   TBranch        *b_dau_XYDCAsig;   //!
   TBranch        *b_dau_ZDCAsig;   //!
   TBranch        *b_dau_eta;   //!
   TBranch        *b_dau_phi;   //!
   TBranch        *b_dau_theta;   //!
   TBranch        *b_dau_vz;   //!
   TBranch        *b_dau_vy;   //!
   TBranch        *b_dau_vx;   //!
   TBranch        *b_dau_vrefz;   //!
   TBranch        *b_dau_vrefy;   //!
   TBranch        *b_dau_vrefx;   //!
   TBranch        *b_dau_vp_difZ;   //!
   TBranch        *b_dau_vp_difY;   //!
   TBranch        *b_dau_vp_difX;   //!
   TBranch        *b_genQScale;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genSignalProcessID;   //!
   TBranch        *b_genJetEta;   //!
   TBranch        *b_genJetPt;   //!
   TBranch        *b_genJetPhi;   //!
   TBranch        *b_genJetChargedMultiplicity;   //!
   TBranch        *b_genDau_chg;   //!
   TBranch        *b_genDau_pid;   //!
   TBranch        *b_genDau_pt;   //!
   TBranch        *b_genDau_eta;   //!
   TBranch        *b_genDau_phi;   //!


   //MyClass(TTree *tree=0);
   MyClass(std::vector<std::string> _fileList);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int job, std::string fList);
   //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef MyClass_cxx
MyClass::MyClass(std::vector<std::string> _fileList) : fChain(0)
//MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MATCH_mc_test.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("MATCH_mc_test.root");
//      }
//      TDirectory * dir = (TDirectory*)f->Get("MATCH_mc_test.root:/analyzerOffline");
//      dir->GetObject("trackTree",tree);

//   }
   //Init(tree);
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
   xVtx = 0;
   yVtx = 0;
   zVtx = 0;
   nTracksVtx = 0;
   ptSumVtx = 0;
   jetNumDaughters = 0;
   jetEta = 0;
   jetPt = 0;
   jetPhi = 0;
   jetTheta = 0;
   jetMass = 0;
   muonMultiplicity = 0;
   chargedMultiplicity = 0;
   dau_pt_sum = 0;
   dau_chg = 0;
   dau_pid = 0;
   dau_vref = 0;
   dau_pt = 0;
   dau_ptError = 0;
   dau_XYDCAsig = 0;
   dau_ZDCAsig = 0;
   dau_eta = 0;
   dau_phi = 0;
   dau_theta = 0;
   dau_vz = 0;
   dau_vy = 0;
   dau_vx = 0;
   dau_vrefz = 0;
   dau_vrefy = 0;
   dau_vrefx = 0;
   dau_vp_difZ = 0;
   dau_vp_difY = 0;
   dau_vp_difX = 0;
   genJetEta = 0;
   genJetPt = 0;
   genJetPhi = 0;
   genJetChargedMultiplicity = 0;
   genDau_chg = 0;
   genDau_pid = 0;
   genDau_pt = 0;
   genDau_eta = 0;
   genDau_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nRun", &nRun, &b_nRun);
   fChain->SetBranchAddress("nEv", &nEv, &b_nRun);
   fChain->SetBranchAddress("nLumi", &nLumi, &b_nLumi);
   fChain->SetBranchAddress("xVtx", &xVtx, &b_xVtx);
   fChain->SetBranchAddress("yVtx", &yVtx, &b_yVtx);
   fChain->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
   fChain->SetBranchAddress("nTracksVtx", &nTracksVtx, &b_nTracksVtx);
   fChain->SetBranchAddress("ptSumVtx", &ptSumVtx, &b_ptSumVtx);
   fChain->SetBranchAddress("jetNumDaughters", &jetNumDaughters, &b_jetNumDaughters);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetTheta", &jetTheta, &b_jetTheta);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("muonMultiplicity", &muonMultiplicity, &b_muonMultiplicity);
   fChain->SetBranchAddress("chargedMultiplicity", &chargedMultiplicity, &b_chargedMultiplicity);
   fChain->SetBranchAddress("jetN", &jetN, &b_jetN);
   fChain->SetBranchAddress("dau_pt_sum", &dau_pt_sum, &b_dau_pt_sum);
   fChain->SetBranchAddress("dau_chg", &dau_chg, &b_dau_chg);
   fChain->SetBranchAddress("dau_pid", &dau_pid, &b_dau_pid);
   fChain->SetBranchAddress("dau_vref", &dau_vref, &b_dau_vref);
   fChain->SetBranchAddress("dau_pt", &dau_pt, &b_dau_pt);
   fChain->SetBranchAddress("dau_ptError", &dau_ptError, &b_dau_ptError);
   fChain->SetBranchAddress("dau_XYDCAsig", &dau_XYDCAsig, &b_dau_XYDCAsig);
   fChain->SetBranchAddress("dau_ZDCAsig", &dau_ZDCAsig, &b_dau_ZDCAsig);
   fChain->SetBranchAddress("dau_eta", &dau_eta, &b_dau_eta);
   fChain->SetBranchAddress("dau_phi", &dau_phi, &b_dau_phi);
   fChain->SetBranchAddress("dau_theta", &dau_theta, &b_dau_theta);
   fChain->SetBranchAddress("dau_vz", &dau_vz, &b_dau_vz);
   fChain->SetBranchAddress("dau_vy", &dau_vy, &b_dau_vy);
   fChain->SetBranchAddress("dau_vx", &dau_vx, &b_dau_vx);
   fChain->SetBranchAddress("dau_vrefz", &dau_vrefz, &b_dau_vrefz);
   fChain->SetBranchAddress("dau_vrefy", &dau_vrefy, &b_dau_vrefy);
   fChain->SetBranchAddress("dau_vrefx", &dau_vrefx, &b_dau_vrefx);
   fChain->SetBranchAddress("dau_vp_difZ", &dau_vp_difZ, &b_dau_vp_difZ);
   fChain->SetBranchAddress("dau_vp_difY", &dau_vp_difY, &b_dau_vp_difY);
   fChain->SetBranchAddress("dau_vp_difX", &dau_vp_difX, &b_dau_vp_difX);
   fChain->SetBranchAddress("genQScale", &genQScale, &b_genQScale);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genSignalProcessID", &genSignalProcessID, &b_genSignalProcessID);
   fChain->SetBranchAddress("genJetEta", &genJetEta, &b_genJetEta);
   fChain->SetBranchAddress("genJetPt", &genJetPt, &b_genJetPt);
   fChain->SetBranchAddress("genJetPhi", &genJetPhi, &b_genJetPhi);
   fChain->SetBranchAddress("genJetChargedMultiplicity", &genJetChargedMultiplicity, &b_genJetChargedMultiplicity);
   fChain->SetBranchAddress("genDau_chg", &genDau_chg, &b_genDau_chg);
   fChain->SetBranchAddress("genDau_pid", &genDau_pid, &b_genDau_pid);
   fChain->SetBranchAddress("genDau_pt", &genDau_pt, &b_genDau_pt);
   fChain->SetBranchAddress("genDau_eta", &genDau_eta, &b_genDau_eta);
   fChain->SetBranchAddress("genDau_phi", &genDau_phi, &b_genDau_phi);
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
