void DrawFlow(){
    TH2D* hSignal[5][2];//Ntrk,pt
    TH2D* hBkg[5][2];
    TH2D* hBkg_ratio[5][2];
    TH1D* hSignal_proj_to_deltaphi[5][3];
    TH1D* hBkg_ratio_proj_to_deltaphi[5][3];
    //TH2D* h2DCorr[5][2];
    TH1D* h1DFlow[5][2];
    TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/ana_run3.root","READ");
    //TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/new_default_complete_vn.root","READ");
    int   trackbinbounds[5]= {76,78,80,81,82};
    int ptbinbounds[2]={3,5};
    TH1D* hJetPass = (TH1D*)f->Get("hJet_Pass550");
    //TH1D* hJetPass = (TH1D*)f->Get("hJet_Pass550_hltCor");
    for(int i=0;i<5;i++){
        for(int j=0;j<2;j++){
    //for(int i=0;i<1;i++){
        //for(int j=0;j<1;j++){
            hSignal[i][j]=(TH2D*)f->Get(Form("hSigS_Cor_%d_to_1000_and_%d_to_30_w_PU_1",trackbinbounds[i],ptbinbounds[j]));
            hBkg[i][j]=(TH2D*)f->Get(Form("hBckS_Cor_%d_to_1000_and_%d_to_30_w_PU_1",trackbinbounds[i],ptbinbounds[j]));
            /*
            hBkg_ratio[i][j]=(TH2D*)hBkg[i][j]->Clone(Form("hBck_ratio_Cor_%d_to_1000_and_%d_to_30_w_PU_1",trackbinbounds[i],ptbinbounds[j]));
            hBkg_ratio[i][j]->Reset(); 
            for(int ieta=0;ieta<hBkg_ratio[i][j]->GetNbinsX();ieta++){
                for(int iphi=0;iphi<hBkg_ratio[i][j]->GetNbinsY();iphi++){
                    double tmp_ratio=hBkg[i][j]->GetBinContent(hBkg[i][j]->FindBin(0,0))/hBkg[i][j]->GetBinContent(ieta+1,iphi+1);
                    double tmp_sgmrt00=hBkg[i][j]->GetBinError(hBkg[i][j]->FindBin(0,0))/hBkg[i][j]->GetBinContent(hBkg[i][j]->FindBin(0,0));
                    double tmp_sgmrt=hBkg[i][j]->GetBinError(ieta+1,iphi+1)/hBkg[i][j]->GetBinContent(ieta+1,iphi+1);
                    hBkg_ratio[i][j]->SetBinContent(ieta+1,iphi+1,tmp_ratio);
                    hBkg_ratio[i][j]->SetBinError(ieta+1,iphi+1,tmp_ratio*sqrt(tmp_sgmrt00*tmp_sgmrt00+tmp_sgmrt*tmp_sgmrt)); 
                }
            }
            */
            hSignal[i][j]->GetXaxis()->SetRange(hSignal[i][j]->GetXaxis()->FindBin(2),hSignal[i][j]->GetNbinsX());
            //hBkg_ratio[i][j]->GetXaxis()->SetRange(hSignal[i][j]->GetXaxis()->FindBin(2),hSignal[i][j]->GetNbinsX());
            hBkg[i][j]->GetXaxis()->SetRange(hBkg[i][j]->GetXaxis()->FindBin(2),hBkg[i][j]->GetNbinsX());
        

            hSignal_proj_to_deltaphi[i][j]=(TH1D*)hSignal[i][j]->ProjectionY(); 
            hBkg_ratio_proj_to_deltaphi[i][j]=(TH1D*)hBkg[i][j]->ProjectionY(); 

            h1DFlow[i][j]=(TH1D*)hSignal_proj_to_deltaphi[i][j]->Clone(Form("h1DFlow_uncor_%d_to_1000_and_%d_to_30_w_PU_1",trackbinbounds[i],ptbinbounds[j]));
            h1DFlow[i][j]->Reset();

            //h2DCorr[i][j]=(TH2D*)hSignal[i][j]->Clone(Form("h2DCorr_Cor_%d_to_1000_and_%d_to_30_w_PU_1",trackbinbounds[i],ptbinbounds[j]));
            //h2DCorr[i][j]->Reset();
            double eta_bw=hSignal[i][j]->GetXaxis()->GetBinWidth(1);//x:delta eta; y:delta phi
            double phi_bw=hSignal[i][j]->GetYaxis()->GetBinWidth(1);
            
            for(int iphi=0;iphi<h1DFlow[i][j]->GetNbinsX();iphi++){
                double tmp_njet=hJetPass->GetBinContent(i+1);
                double tmp_bin=hSignal_proj_to_deltaphi[i][j]->GetBinContent(iphi+1)
                /(hBkg_ratio_proj_to_deltaphi[i][j]->GetBinContent(iphi+1)*tmp_njet);
                double tmp_sgmrt_S=hSignal_proj_to_deltaphi[i][j]->GetBinError(iphi+1)/hSignal_proj_to_deltaphi[i][j]->GetBinContent(iphi+1);
                double tmp_sgmrt_B=hBkg_ratio_proj_to_deltaphi[i][j]->GetBinError(iphi+1)/hBkg_ratio_proj_to_deltaphi[i][j]->GetBinContent(iphi+1);
                if(!std::isfinite(tmp_bin)) continue;
                h1DFlow[i][j]->SetBinContent(iphi+1,tmp_bin);
                h1DFlow[i][j]->SetBinError(iphi+1,tmp_bin*sqrt(tmp_sgmrt_S*tmp_sgmrt_S+tmp_sgmrt_B*tmp_sgmrt_B)); 

            }
        
        }
    }
    
    //TFile* fout= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/2DCorr_run2_PG_ana_eta_phi_bw.root","RECREATE");
    //TFile* fout= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/2DCorr_run3.root","RECREATE");
    TFile* fout= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/1DFlow_run3.root","RECREATE");
    
    for(int i=0;i<5;i++){
        for(int j=0;j<2;j++){
    //for(int i=0;i<1;i++){
        //for(int j=0;j<1;j++){
            //h2DCorr[i][j]->Write();
            h1DFlow[i][j]->Write();
        }
    }

    f->Close();
    delete f;
    fout->Close();
    delete fout;

}