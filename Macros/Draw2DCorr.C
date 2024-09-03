void Draw2DCorr(){
    TH2D* hSignal[5][2];//Ntrk,pt
    TH2D* hBkg[5][2];
    TH2D* h2DCorr[5][2];
    //TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/ana_run3_allNch_newBkg_for_high_Nch.root","READ");
    TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/output_parkersbin_newBkg_for_high_Nch.root","READ");
    //TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/new_default_complete_vn.root","READ");
    //int   trackbinbounds[5]= {76,78,80,81,82};
    //int ptbinbounds[2]={3,5};
    //int   trackbinbounds[5]= {83,78,80,81,82};
    int ptbinbounds[2]={3,5};
    const int   trackbinbounds[10]         = { 0,20,30,40,50,59,66,76,83,78};
    const int   trackbinboundsUpper[10]    = {20,30,40,50,59,66,76,83,1000,1000};
    TH1D* hJetPass = (TH1D*)f->Get("hJet_Pass550");
    //TH1D* hJetPass = (TH1D*)f->Get("hJet_Pass550_hltCor");
    for(int i=0;i<10;i++){
        for(int j=0;j<2;j++){
    //for(int i=0;i<1;i++){
        //for(int j=0;j<1;j++){
            hSignal[i][j]=(TH2D*)f->Get(Form("hSigS_Cor_%d_to_%d_and_%d_to_30_w_PU_1",trackbinbounds[i],trackbinboundsUpper[i],ptbinbounds[j]));
            hBkg[i][j]=(TH2D*)f->Get(Form("hBckS_Cor_%d_to_%d_and_%d_to_30_w_PU_1",trackbinbounds[i],trackbinboundsUpper[i],ptbinbounds[j]));
            h2DCorr[i][j]=(TH2D*)hSignal[i][j]->Clone(Form("h2DCorr_Cor_%d_to_%d_and_%d_to_30_w_PU_1",trackbinbounds[i],trackbinboundsUpper[i],ptbinbounds[j]));
            h2DCorr[i][j]->Reset();
            double x_bw=hSignal[i][j]->GetXaxis()->GetBinWidth(1);
            double y_bw=hSignal[i][j]->GetYaxis()->GetBinWidth(1);
            for(int ieta=0;ieta<h2DCorr[i][j]->GetNbinsX();ieta++){
                for(int iphi=0;iphi<h2DCorr[i][j]->GetNbinsY();iphi++){
                    double tmp_njet=hJetPass->GetBinContent(i+1);
                    //double tmp_njet=hJetPass->GetBinContent(9);
                    double tmp_sig=hSignal[i][j]->GetBinContent(ieta+1,iphi+1);
                    double tmp_bkg=hBkg[i][j]->GetBinContent(ieta+1,iphi+1);
                    double tmp_b00=hBkg[i][j]->GetBinContent(hBkg[i][j]->FindBin(0,0));
                    //double tmp_bin=tmp_b00*tmp_sig/(tmp_bkg*tmp_njet);
                    double tmp_bin=tmp_b00*tmp_sig/(tmp_bkg*tmp_njet*x_bw*y_bw);
                    //double tmp_bin=tmp_b00*tmp_sig/(tmp_bkg*tmp_njet*x_bw);
                    //cout<<tmp_bin<<endl;
                    //if(std::isnan(tmp_bin)) continue;
                    if(!std::isfinite(tmp_bin)) continue;
                    h2DCorr[i][j]->SetBinContent(ieta+1,iphi+1,tmp_bin); 
                }
            }
        
        }
    }
    
    //TFile* fout= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/2DCorr_run2_PG_ana_eta_phi_bw.root","RECREATE");
    TFile* fout= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/2DCorr_run2_allNch.root","RECREATE");
    
    for(int i=0;i<10;i++){
        for(int j=0;j<2;j++){
    //for(int i=0;i<1;i++){
        //for(int j=0;j<1;j++){
            h2DCorr[i][j]->Write();
        }
    }

    f->Close();
    delete f;
    fout->Close();
    delete fout;
}