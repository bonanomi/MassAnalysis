void fitMass(int year=2018, TString pmode = "ggH"){

	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	TString chanName[3]={"4mu","4e","2e2mu"};
	int mass [5]={125,124,120,126,130};

	RooRealVar *ZZMass =new RooRealVar("ZZMass","",105,140);
	RooRealVar *chan=new RooRealVar("chan","",0,1000);
	RooRealVar *weight=new RooRealVar("weight","",1.e-9,1.e9);

	for (int j =0;j<3;j++){
		ofstream outmass ;
		outmass.open(Form("sim_massParam_%s_%s_%d.txt",pmode.Data(), chanName[j].Data(),year));
		RooDataSet *data[5];
		for (int i =0;i<5;i++){
			TChain *t = new TChain ("SelectedTree");
			t->Add(Form("%s%d_redTree_%d.root",pmode.Data(),mass[i],year));
			data[i]=new RooDataSet(Form("data%d_%d",i,mass[i]),"",t,RooArgSet(*ZZMass,*weight,*chan),Form("chan==%d",j+1),"weight");
		}

		RooCategory massrc("massrc","");
		RooSimultaneous simPdf("simPdf","simultaneous pdf",massrc) ;

		RooDataSet* data_stage_cat[5]; 
		RooDataSet* data_all;
		RooDoubleCBFast* cpdf[5]; 
		RooAddPdf* addpdf[5]; 

		for (int i=0;i<5;i++){
			TString cat_type =Form("mh%d",mass[i]);
			massrc.defineType(cat_type,mass[i]);
		}

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

		for (int i=0;i<5;i++){
			RooConstVar *MH=new RooConstVar("MH","",mass[i]);
			RooFormulaVar* a1=new RooFormulaVar(Form("a1_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*a1_0,*a1_1,*MH));
			RooFormulaVar* a2=new RooFormulaVar(Form("a2_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*a2_0,*a2_1,*MH));
			RooFormulaVar* n1=new RooFormulaVar(Form("n1_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*n1_0,*n1_1,*MH));
			RooFormulaVar* n2=new RooFormulaVar(Form("n2_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*n2_0,*n2_1,*MH));
			RooFormulaVar* mean=new RooFormulaVar(Form("mean_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*mean_0,*mean_1,*MH));
			RooFormulaVar* sigma=new RooFormulaVar(Form("sigma_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*sigma_0,*sigma_1,*MH));
			RooFormulaVar* sigma_l=new RooFormulaVar(Form("sigma_l_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*sigma_l_0,*sigma_l_1,*MH));
			RooFormulaVar* mean_l=new RooFormulaVar(Form("mean_l_%d",mass[i]),"","@0+@1*(MH-125)",RooArgList(*mean_l_0,*mean_l_1,*MH));
			RooFormulaVar* frac=new RooFormulaVar(Form("frac_%d",mass[i]),"","@0",RooArgList(*frac_0));

			TString cat_type =Form("mh%d",mass[i]);
			data_stage_cat[i] = new RooDataSet(Form("data_%d",i),"",RooArgSet(*ZZMass,*weight),Index(massrc),Import(cat_type,*data[i]),WeightVar("weight")); 
			cpdf[i] = new RooDoubleCBFast(Form("DCBall_%d",i),"",*ZZMass,*mean,*sigma,*a1,*n1,*a2,*n2);
			RooLandau *pdf_landau = new RooLandau("landau","",*ZZMass, *mean_l, *sigma_l);
			addpdf[i] = new RooAddPdf(Form("apdf_%d",i),"",RooArgList(*cpdf[i],*pdf_landau),*frac);
			cout<< data_stage_cat[i]->sumEntries()<<endl;
			
			// Additional Landau contribution for WH and ZH
			if (pmode.Contains("Wplus")|| pmode.Contains("Wminus") || pmode.Contains("ZH"))
				simPdf.addPdf(*addpdf[i],cat_type);
			else
				simPdf.addPdf(*cpdf[i],cat_type);

			if(i==0)
				data_all = data_stage_cat[i];
			else
				data_all->append(*data_stage_cat[i]);

			if(i==0){
				RooDoubleCBFast *cpdftmp = new RooDoubleCBFast("DCBall","",*ZZMass,*mean_0,*sigma_0,*a1_0,*n1_0,*a2_0,*n2_0);
				RooLandau *pdf_landautmp = new RooLandau("landautmp","",*ZZMass, *mean_l_0, *sigma_l_0);
				RooAddPdf *addpdftmp = new RooAddPdf(Form("apdf_%d",i),"",RooArgList(*cpdftmp,*pdf_landautmp),*frac);
				if (pmode.Contains("Wplus")|| pmode.Contains("Wminus") || pmode.Contains("ZH"))
					addpdftmp->fitTo(*data_stage_cat[i],InitialHesse(true),Strategy(2));
				else
					cpdftmp->fitTo(*data_stage_cat[i],InitialHesse(true),Strategy(2));

				sigma_0->setConstant(1);
				mean_0->setConstant(1);
				mean_l_0->setConstant(1);
			}
		}
		simPdf.Print("v");
		data_all->Print("v");

		simPdf.fitTo(*data_all,InitialHesse(true),Strategy(2)) ;

		TString a1outs= Form ("a1_%s_%d \t %.4f+%.5f*(MH-125)",chanName[j].Data(),year,a1_0->getVal(),a1_1->getVal());
		TString a2outs= Form ("a2_%s_%d \t %.4f+%.5f*(MH-125)",chanName[j].Data(),year,a2_0->getVal(),a2_1->getVal());
		TString n1outs= Form ("n1_%s_%d \t %.4f+%.5f*(MH-125)",chanName[j].Data(),year,n1_0->getVal(),n1_1->getVal());
		TString n2outs= Form ("n2_%s_%d \t %.4f+%.5f*(MH-125)",chanName[j].Data(),year,n2_0->getVal(),n2_1->getVal());
		TString meanouts= Form ("mean_%s_%d \t %.4f+%.5f*(MH-125)",chanName[j].Data(),year,mean_0->getVal(),mean_1->getVal());
		TString sigmaouts= Form ("sigma_%s_%d \t %.4f+%.5f*(MH-125)",chanName[j].Data(),year,sigma_0->getVal(),sigma_1->getVal());

		outmass<< a1outs<<endl;
		outmass<< a2outs<<endl;
		outmass<< meanouts<<endl;
		outmass<< sigmaouts<<endl;
		outmass<< n1outs<<endl;
		outmass<< n2outs<<endl;
		
		// Additional Landau contribution for WH and ZH
		if (pmode.Contains("Wplus")|| pmode.Contains("Wminus") || pmode.Contains("ZH")){
			outmass<<Form ("mean_l_%s_%d \t %.4f+%.5f*(MH-125)",chanName[j].Data(),year,mean_l_0->getVal(),mean_l_1->getVal())<<endl;
			outmass<<Form ("sigma_l_%s_%d \t %.4f+%.5f*(MH-125)",chanName[j].Data(),year,sigma_l_0->getVal(),sigma_l_1->getVal())<<endl;
			outmass<<Form ("frac_%s_%d \t %.4f",chanName[j].Data(),year,frac_0->getVal())<<endl;
		}
		cout<<"where1"<<endl;
		for (int i =0;i<5;i++){
			RooPlot *frame = ZZMass->frame();
			TString cat_type =Form("mh%d",mass[i]);
			data_all->plotOn(frame,Cut(Form("massrc==massrc::mh%d",mass[i]))) ;
			simPdf.plotOn(frame,Slice(massrc,cat_type),ProjWData(massrc,*data_all)) ;
			frame->Draw();
			gPad->Print(Form("simFit_%d_%s_%d.png",mass[i],chanName[j].Data(),year));
		}
		outmass.close();
	}
}
