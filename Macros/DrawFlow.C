void DrawFlow(){
    TH2D* hSignal[10][2];//Ntrk,pt
    TH2D* hBkg[10][2];
    //TH2D* h2DCorr[5][2];
    TH1D* h1DFlow[10][2];
    //TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/ana_run3.root","READ");
    //TFile* f_bkg= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/Bkg_highMult_jets_run3.root","READ");
    //TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/new_default_complete_vn.root","READ");
    //TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/ana_run3_allNch.root","READ");
    TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/ana_run3_2022_jetveto_C_newBkg.root","READ");
    //TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/output_parkersbin.root","READ");
    //TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/output_parkersbin_newBkg_for_high_Nch.root","READ");
    
    //int   trackbinbounds[5]= {76,78,80,81,82};
    int ptbinbounds[2]={3,5};
    const int   trackbinbounds[10]         = { 0,20,30,40,50,59,66,76,83,78};
    const int   trackbinboundsUpper[10]    = {20,30,40,50,59,66,76,83,1000,1000};
    
    float ptname[2]={0.3,0.5};
    int YPlo=28;
    int YPhi=34;
    TH1D* hJetPass = (TH1D*)f->Get("hJet_Pass550");
    //TH1D* hJetPass = (TH1D*)f->Get("hJet_Pass550_hltCor");
    
    //for(int i=0;i<5;i++){
    for(int i=0;i<10;i++){
        for(int j=0;j<2;j++){
            //hSignal[i][j]=(TH2D*)f->Get(Form("hSigS_Cor_%d_to_1000_and_%d_to_30_w_PU_1",trackbinbounds[i],ptbinbounds[j]));
            //hBkg[i][j]=(TH2D*)f_bkg->Get(Form("hBckS_Cor_%d_to_1000_and_%d_to_30_w_PU_1",trackbinbounds[i],ptbinbounds[j]));
            
            hSignal[i][j]=(TH2D*)f->Get(Form("hSigS_Cor_%d_to_%d_and_%d_to_30_w_PU_1",trackbinbounds[i],trackbinboundsUpper[i],ptbinbounds[j]));
            hBkg[i][j]=(TH2D*)f->Get(Form("hBckS_Cor_%d_to_%d_and_%d_to_30_w_PU_1",trackbinbounds[i],trackbinboundsUpper[i],ptbinbounds[j]));
            
            //h2DCorr[i][j]=(TH2D*)hSignal[i][j]->Clone(Form("h2DCorr_Cor_%d_to_1000_and_%d_to_30_w_PU_1",trackbinbounds[i],ptbinbounds[j]));
            //h2DCorr[i][j]->Reset();
            double eta_bw=hSignal[i][j]->GetXaxis()->GetBinWidth(1);//x:delta eta; y:delta phi
            double phi_bw=hSignal[i][j]->GetYaxis()->GetBinWidth(1);

            hSignal[i][j]->Scale(1.0/(hJetPass->GetBinContent(i+1,j+1)));
            hSignal[i][j]->Scale(1./(phi_bw));

            TH1D *histfit1 = (TH1D*) hSignal[i][j]->ProjectionY("",YPlo,YPhi)->Clone();
            TH1D *histfit2 = (TH1D*) hBkg[i][j]->ProjectionY("",YPlo,YPhi)->Clone();

            histfit1->Divide(histfit2);
            histfit1->Scale(hBkg[i][j]->GetMaximum());

            h1DFlow[i][j]=(TH1D*)histfit1->Clone(); 
            std::string function = "[0]/(TMath::Pi()*2)*(1+2*([1]*TMath::Cos(x)+[2]*TMath::Cos(2*x)+[3]*TMath::Cos(3*x)+[4]*TMath::Cos(4*x)+[5]*TMath::Cos(5*x)))";
            TF1 func1("deltaPhi1", function.c_str(), -0.5*TMath::Pi(), 1.5*TMath::Pi());
            func1.SetParameter(0, histfit1->GetMaximum());
            func1.SetParameter(1, 0.1);
            func1.SetParameter(2, 0.1);
            func1.SetParameter(3, 0.1);
            func1.SetParameter(4, 0.1);
            func1.SetParameter(5, 0.1);
            h1DFlow[i][j]->Fit(&func1, "m E q");
            

            // Set the histogram title
            //h1DFlow[i][j]->SetTitle(Form("Nch>%d, %.2f<jT<3 GeV/c",trackbinbounds[i],ptname[j]));
            h1DFlow[i][j]->SetTitle(Form("%d<Nch<%d, %.2f<jT<3 GeV/c",trackbinbounds[i],trackbinboundsUpper[i],ptname[j]));
             
        }
    }
    
    TCanvas *c1 = new TCanvas("canvas", "Fourier Series Fits", 800, 1200);
    c1->Divide(2, 5); // Divide canvas into 2 columns and 5 rows
    for(int i=0;i<5;i++){
        for(int j=0;j<2;j++){
            c1->cd(i*2+j+1);
            h1DFlow[i][j]->Draw();
        }
    } 
    c1->SaveAs("/Users/xl155/Documents/JetFlow_Run3_data/Flow_run3_allNch_2022_C_jetveto.pdf(");
    //c1->SaveAs("/Users/xl155/Documents/JetFlow_Run3_data/Flow_run2_unc_new.pdf(");

    TCanvas *c2 = new TCanvas("canvas", "Fourier Series Fits", 800, 1200);
    c2->Divide(2, 5); // Divide canvas into 2 columns and 5 rows
    for(int i=0;i<5;i++){
        for(int j=0;j<2;j++){
            c2->cd(i*2+j+1);
            h1DFlow[i+5][j]->Draw();
        }
    } 
    c2->SaveAs("/Users/xl155/Documents/JetFlow_Run3_data/Flow_run3_allNch_2022_C_jetveto.pdf)");
    //c2->SaveAs("/Users/xl155/Documents/JetFlow_Run3_data/Flow_run2_unc_new.pdf)");
    
    /* 
    //TFile* fout= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/2DCorr_run2_PG_ana_eta_phi_bw.root","RECREATE");
    //TFile* fout= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/2DCorr_run3.root","RECREATE");
    TFile* fout= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/1DFlow_run3_allNch_2023_2024.root","RECREATE");
    //TFile* fout= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/1DFlow_run2_unc_new.root","RECREATE");
    
    //for(int i=0;i<5;i++){
    for(int i=0;i<10;i++){
        for(int j=0;j<2;j++){
            //h2DCorr[i][j]->Write();
            h1DFlow[i][j]->Write();
        }
    }
    */

    f->Close();
    delete f;
    //fout->Close();
    //delete fout;

}