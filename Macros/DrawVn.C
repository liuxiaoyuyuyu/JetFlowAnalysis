void DrawVn(){
    TH2D* hSignal[10][2];//Ntrk,pt
    TH2D* hBkg[10][2];
    TH1D* h1DFlow[10][2];
    TH1D* hBinDist[10];
    
    TFile* f= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/ana_run3_2022_jetveto_C_newBkg.root","READ");
    TFile* fstat= new TFile("/Users/xl155/Documents/JetFlow_Run3_data/ana_run3_2022_jetveto_C.root","READ");
    
    //int   trackbinbounds[5]= {76,78,80,81,82};
    int ptbinbounds[2]={3,5};
    const int   trackbinbounds[10]         = { 0,20,30,40,50,59,66,76,83,78};
    const int   trackbinboundsUpper[10]    = {20,30,40,50,59,66,76,83,1000,1000};
    
    float ptname[2]={0.3,0.5};
    int YPlo=28;
    int YPhi=34;
    TH1D* hJetPass = (TH1D*)f->Get("hJet_Pass550");
    //TH1D* hJetPass = (TH1D*)f->Get("hJet_Pass550_hltCor");
    
    double mAve_Nch[10];
    double mAve_Nch_err[10];
    double mFit_par[2][5][10];
    double mFit_par_err[2][5][10];

    //for(int i=0;i<5;i++){
    for(int i=0;i<10;i++){
        hBinDist[i]=(TH1D*)fstat->Get(Form("hBinDist_unc_%d",i+1)); 
        mAve_Nch[i]=hBinDist[i]->GetMean();
        mAve_Nch_err[i]=0.0;

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

            for(int ip=0;ip<5;ip++){
                mFit_par[j][ip][i]=func1.GetParameter(ip+1);
                mFit_par_err[j][ip][i]=func1.GetParError(ip+1);
                mFit_par_err[j][ip][i]*=sqrt(2); 
            } 

            // Set the histogram title
            //h1DFlow[i][j]->SetTitle(Form("Nch>%d, %.2f<jT<3 GeV/c",trackbinbounds[i],ptname[j]));
            h1DFlow[i][j]->SetTitle(Form("%d<Nch<%d, %.2f<jT<3 GeV/c",trackbinbounds[i],trackbinboundsUpper[i],ptname[j]));
             
        }
    }

    int color[5]={632,419,600,800,616};

    TGraphErrors* gVn[2][5];
    for(int j=0;j<2;j++){
        for(int ip=0;ip<5;ip++){
            gVn[j][ip]=new TGraphErrors(10,mAve_Nch,mFit_par[j][ip],mAve_Nch_err,mFit_par_err[j][ip]);
            gVn[j][ip]->SetMarkerStyle(20);
            gVn[j][ip]->SetMarkerSize(0.8);
            gVn[j][ip]->SetLineColor(color[ip]);
            gVn[j][ip]->SetMarkerColor(color[ip]);
        }
    }

    TLine* l1=new TLine(9.0,0.0,90.0,0.0);
    TCanvas* cVn=new TCanvas("cVn","cVn",600,500);
    gVn[0][1]->Draw("AP");
    gVn[0][1]->SetTitle("Run 3 2022 C, jet veto applied");
    gVn[0][1]->GetXaxis()->SetTitle("N_{ch}^{j}");
    gVn[0][1]->GetXaxis()->SetRange(5,95); 
    gVn[0][1]->GetYaxis()->SetTitle("V^{*}_{n#Delta}{2,|#Delta#eta^{*}|>2}"); 
    gVn[0][1]->GetYaxis()->SetTitleOffset(1.2); 
    gVn[0][1]->SetMaximum(0.1);
    gVn[0][1]->SetMinimum(-0.27);
    gVn[0][0]->Draw("P");
    gVn[0][2]->Draw("P");
    l1->Draw("same");

    /* 
    TCanvas *c1 = new TCanvas("canvas", "Fourier Series Fits", 800, 1200);
    c1->Divide(2, 5); // Divide canvas into 2 columns and 5 rows
    for(int i=0;i<5;i++){
        for(int j=0;j<2;j++){
            c1->cd(i*2+j+1);
            h1DFlow[i][j]->Draw();
        }
    } 
    c1->SaveAs("/Users/xl155/Documents/JetFlow_Run3_data/Flow_run3_allNch_2023_2024.pdf(");
    //c1->SaveAs("/Users/xl155/Documents/JetFlow_Run3_data/Flow_run2_unc_new.pdf(");

    TCanvas *c2 = new TCanvas("canvas", "Fourier Series Fits", 800, 1200);
    c2->Divide(2, 5); // Divide canvas into 2 columns and 5 rows
    for(int i=0;i<5;i++){
        for(int j=0;j<2;j++){
            c2->cd(i*2+j+1);
            h1DFlow[i+5][j]->Draw();
        }
    } 
    c2->SaveAs("/Users/xl155/Documents/JetFlow_Run3_data/Flow_run3_allNch_2023_2024.pdf)");
    //c2->SaveAs("/Users/xl155/Documents/JetFlow_Run3_data/Flow_run2_unc_new.pdf)");
    */
    f->Close();
    delete f;
}