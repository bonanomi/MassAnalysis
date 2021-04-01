void fitMass_125(int year=2018){

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	TString chanName[3]={"4mu","4e","2e2mu"};

	RooRealVar *ZZMass =new RooRealVar("ZZMass","",105,140);
	RooRealVar *chan=new RooRealVar("chan","",0,1000);
	RooRealVar *weight=new RooRealVar("weight","",1.e-9,1.e9);

	for (int j =0;j<3;j++){
		// ofstream outmass ;
		// outmass.open(Form("shape/sim_massParam_htxs_stage1_reco_cat%s%d.txt",chanName[j].Data(),year));
		TChain *t = new TChain ("SelectedTree");
		t->Add("ggH125_redTree_2016.root");
		RooDataSet *data = new RooDataSet("data","",t,RooArgSet(*ZZMass,*weight,*chan),Form("chan==%d",j+1),"weight");
		RooRealVar *a1_1=new RooRealVar("a1_1","",0.00111,-0.2,0.2);
		RooRealVar *a1_2=new RooRealVar("a1_2","",0.00146,-0.1,0.1);
		RooRealVar *a1_0=new RooRealVar("a1_0","",1.2917,1.,5.);
		RooRealVar *a2_1=new RooRealVar("a2_1","",0.00130,-0.2,0.2);
		RooRealVar *a2_2=new RooRealVar("a2_2","",0.00127,-0.1,0.1);
		RooRealVar *a2_0=new RooRealVar("a2_0","",1.8125,1.,5.);
		RooRealVar *n1_1=new RooRealVar("n1_1","",0.00138,-0.2,0.2);
		RooRealVar *n1_2=new RooRealVar("n1_2","",0.00185,-0.1,0.1);
		RooRealVar *n1_0=new RooRealVar("n1_0","",2.0312,1.,7.);
		RooRealVar *n2_1=new RooRealVar("n2_1","",0.001,-0.2,0.2);
		RooRealVar *n2_2=new RooRealVar("n2_2","",-0.001,-0.1,0.1);
		RooRealVar *n2_0=new RooRealVar("n2_0","",3.0369,1.,5.);

		RooRealVar *mean_1=new RooRealVar("mean_1","",1,0,1.5);
		RooRealVar *mean_2=new RooRealVar("mean_2","",0.99945,-1.,1.);
		RooRealVar *mean_0=new RooRealVar("mean_0","",124.8424,120.,130);

		RooRealVar *sigma_1=new RooRealVar("sigma_1","",0.00860,-0.5,0.5);
		RooRealVar *sigma_2=new RooRealVar("sigma_2","",0.00026,-1.,1.);
		RooRealVar *sigma_0=new RooRealVar("sigma_0","",1.5,0.5,3.);

		RooRealVar* mean_l_0 = new RooRealVar("landau_mean_0","",130,110,140);
		RooRealVar* mean_l_1 = new RooRealVar("landau_mean_1","",0,-1.5,1.5);
		RooRealVar* mean_l_2 = new RooRealVar("landau_mean_2","",0,-1,1);
		RooRealVar* sigma_l_0 = new RooRealVar("landau_sigma_0","",15,2,20);
		RooRealVar* sigma_l_1 = new RooRealVar("landau_sigma_1","",0.,-1,1);
		RooRealVar* sigma_l_2 = new RooRealVar("landau_sigma_2","",0,-1,1);
		RooRealVar* frac_0 = new RooRealVar("frac_0","",0.65,0,1);
		RooRealVar* frac_1 = new RooRealVar("frac_1","",-0.1,0.1);
		RooRealVar* frac_2 = new RooRealVar("frac_2","",-0.1,0.1);
		
		RooConstVar *MH=new RooConstVar("MH","",125);

		RooFormulaVar* a1=new RooFormulaVar(Form("a1_125"),"","@0+@1*(MH-125)",RooArgList(*a1_0,*a1_1,*MH));
		RooFormulaVar* a2=new RooFormulaVar(Form("a2_125"),"","@0+@1*(MH-125)",RooArgList(*a2_0,*a2_1,*MH));
		RooFormulaVar* n1=new RooFormulaVar(Form("n1_125"),"","@0+@1*(MH-125)",RooArgList(*n1_0,*n1_1,*MH));
		RooFormulaVar* n2=new RooFormulaVar(Form("n2_125"),"","@0+@1*(MH-125)",RooArgList(*n2_0,*n2_1,*MH));
		RooFormulaVar* mean=new RooFormulaVar(Form("mean_125"),"","@0+@1*(MH-125)",RooArgList(*mean_0,*mean_1,*MH));
		RooFormulaVar* sigma=new RooFormulaVar(Form("sigma_125"),"","@0+@1*(MH-125)",RooArgList(*sigma_0,*sigma_1,*MH));
		RooFormulaVar* sigma_l=new RooFormulaVar(Form("sigma_l_125"),"","@0+@1*(MH-125)",RooArgList(*sigma_l_0,*sigma_l_1,*MH));
		RooFormulaVar* mean_l=new RooFormulaVar(Form("mean_l_125"),"","@0+@1*(MH-125)",RooArgList(*mean_l_0,*mean_l_1,*MH));

		RooDoubleCBFast *cpdftmp = new RooDoubleCBFast(Form("DCBall_125"),"",*ZZMass,*mean,*sigma,*a1,*n1,*a2,*n2);

		// RooDoubleCBFast *cpdftmp = new RooDoubleCBFast("DCBall","",*ZZMass,*mean_0,*sigma_0,*a1_0,*n1_0,*a2_0,*n2_0);
		cpdftmp->fitTo(*data,InitialHesse(true),Strategy(2));

		RooPlot *frame = ZZMass->frame();
		TString cat_type =Form("mh125");
		data->plotOn(frame);
		cpdftmp->plotOn(frame); //,Slice(massrc,cat_type),ProjWData(massrc,*data_all)) ;
		frame->Draw();
		gPad->Print(Form("simFit_%s.pdf",chanName[j].Data()));
	}
}
