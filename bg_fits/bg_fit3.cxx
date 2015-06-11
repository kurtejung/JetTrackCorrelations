#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLatex.h"

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include "../HIN_14_016_functions.h"


using namespace std;

Int_t findbin(double x);
Float_t find_v2(int in_pt_i, int in_Cent_j);
TH1D *Rebin_dPhi(TH1D* hold);
TH1D *Rebin_dEta(TH1D* hold);
TString make_name(TString stem, int g, int i, int j, int l, TString &centlabel, TString &pTlabel);

Int_t bg_fit3(int gstart = 0, int gend = 6,bool skip_pp = kTRUE, bool Is_NonLeading = kFALSE)
{

#include "../HIN_14_016_universals.h"

  gROOT->ForceStyle();
  gStyle->SetOptDate(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(1);

  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);


  TFile *fdata[12];
  TFile *fout[12];

  double phimin = -1.5;
  double phimax = 4.5;

  float x, y;
  TF1 *fit0 = new TF1("fit0","[0]",-3.,3.);

  Float_t v2_assoc[12][5][4][2];
  Float_t v2_jet[12][4][2];
  Float_t v2_jet_test[12][5][4][2];
  Float_t v2_jet_err[12][4][2];

  TH1D *raw_eta[12][5][4][2];

  TH2D *yield[12][5][4][2];
  TH1D *yield_proj[12][5][4][2];
 
  TH1D *yield_rebin[12][5][4][2];
  TH1D *in_hist_test[12][5][4][2];
  TH1D *in_hist[12][5][4][2];
  TH1D *in_hist_shift[12][5][4][2];
  TH1D *in_hist_rebin[12][5][4][2];
  TH1D *in_hist_zoomed[12][5][4][2];
  TH1D *in_hist_zoomed_rebin[12][5][4][2];
  TH1D *in_hist_ref[12][5][4][2];
  TH2D *in_hist2[12][5][4][2];
  TH2D *out_hist[12][5][4][2];
  TH1D *out_hist1[12][5][4][2];
  TH1D *out_hist1_rebin[12][5][4][2];
  TH1D *out_hist_rebin[12][5][4][2];
  TH2D *test_diff[12][5][4][2];
  TH1D *test_diff2[12][5][4][2];
  TH1D *test_diff3[12][5][4][2]; 
  TH1D *test_diff4[12][5][4][2];
  TH1D *test_diff5[12][5][4][2];
  
  TH1D *sideband_hist[12][5][4][2];

  TH2D *result[12][5][4][2];
  TH2D *result2[12][5][4][2];

  TF1 *fitfunc[12][5][4][2];
  TF1 *fitmin[12][5][4][2];
  TF1 *fitmax[12][5][4][2];
  TF1 *funcmax[12][5][4][2];
  TF1 *funcmin[12][5][4][2];
  
  double error[12][5][4][2];

  TCanvas *cfit[12][2];
  TCanvas *czoomed[12][2];
  TCanvas *cdiff_phi[12][5];
  TCanvas *cdiff_eta[12][5];
  TCanvas *craweta[12][5];
  
  double alpha_array_pT0[5];
  double alpha_array_pT1[5];
  double alpha_array_pT2[5];
  double alpha_array_pT3[5];
  
  double beta_array_pT0[5];
  double beta_array_pT1[5];
  double beta_array_pT2[5];
  double beta_array_pT3[5];


  double V2_array_pT0[5];
  double V2_array_pT1[5];
  double V2_array_pT2[5];
  double V2_array_pT3[5];

  double alpha_array_Cent0[4];
  double alpha_array_Cent1[4];
  double alpha_array_Cent2[4];
  double alpha_array_Cent3[4];
  double alpha_array_Cent4[4];

  double beta_array_Cent0[4];
  double beta_array_Cent1[4];
  double beta_array_Cent2[4];
  double beta_array_Cent3[4];
  double beta_array_Cent4[4];

  double V2_array_Cent0[4];
  double V2_array_Cent1[4];
  double V2_array_Cent2[4];
  double V2_array_Cent3[4];
  double V2_array_Cent4[4];


  double alpha_error_pT0[5];
  double alpha_error_pT1[5];
  double alpha_error_pT2[5];
  double alpha_error_pT3[5];


  double beta_error_pT0[5];
  double beta_error_pT1[5];
  double beta_error_pT2[5];
  double beta_error_pT3[5];


  double V2_error_pT0[5];
  double V2_error_pT1[5];
  double V2_error_pT2[5];
  double V2_error_pT3[5];

  double alpha_error_Cent0[4];
  double alpha_error_Cent1[4];
  double alpha_error_Cent2[4];
  double alpha_error_Cent3[4];
  double alpha_error_Cent4[4];

  double beta_error_Cent0[4];
  double beta_error_Cent1[4];
  double beta_error_Cent2[4];
  double beta_error_Cent3[4];
  double beta_error_Cent4[4];

  double V2_error_Cent0[4];
  double V2_error_Cent1[4];
  double V2_error_Cent2[4];
  double V2_error_Cent3[4];
  double V2_error_Cent4[4];

  double pT_bins0[5] = {1.3,2.3,3.3,5.8};
  double pT_bins1[5] = {1.4,2.4,3.4,5.9};
  double pT_bins2[5] = {1.5,2.5,3.5,6.};
  double pT_bins3[5] = {1.6,2.6,3.6,6.1};
  double pT_error[5] = {0.5,0.5,0.5,2.};
  double Cent_bins0[4] = {73.,38.,18.,3.};
  double Cent_bins1[4] = {74.,39.,19.,4.};
  double Cent_bins2[4] = {75.,41.,20.,5.};
  double Cent_bins3[4] = {76.,42.,21.,6.};
  double Cent_error[4] = {25.,10.,10.,5.};

	
  float yleft[12][5][4][2];
  float yright[12][5][4][2];
	

  int bin;

  //----------------------------------------
  //   v2_jet[g][i][l] (for all centralities)
  //-------------------------------------------

  //Inclusive
  v2_jet[0][0][0] = 0.054; 
  v2_jet[0][1][0] = 0.054;
  v2_jet[0][2][0] = 0.061;
  v2_jet[0][3][0] = 0.084;
 
  //SubLeading
  v2_jet[2][0][0] = 0.054;
  v2_jet[2][1][0] = 0.054;
  v2_jet[2][2][0] = 0.037;
  v2_jet[2][3][0] = 0.022;

  //Leading
  v2_jet[4][0][0] = 0.054;
  v2_jet[4][1][0] = 0.054;
  v2_jet[4][2][0] = 0.061;
  v2_jet[4][3][0] = 0.084;



  // v2_jet_error
 

  //Inclusive
  v2_jet_err[0][0][0] = 0.025; 
  v2_jet_err[0][1][0] = 0.025;
  v2_jet_err[0][2][0] = 0.03;
  v2_jet_err[0][3][0] = 0.035;
 
  //SubLeading
  v2_jet_err[2][0][0] = 0.025;
  v2_jet_err[2][1][0] = 0.025;
  v2_jet_err[2][2][0] = 0.03;
  v2_jet_err[2][3][0] = 0.035;

  //Leading
  v2_jet_err[4][0][0] = 0.025;
  v2_jet_err[4][1][0] = 0.025;
  v2_jet_err[4][2][0] = 0.03;
  v2_jet_err[4][3][0] = 0.035;
 


  for(int j = 0; j<4; j++){
    v2_jet[6][j][0] = v2_jet[0][j][0];
    v2_jet[8][j][0] = v2_jet[2][j][0];
    v2_jet[10][j][0] = v2_jet[4][j][0];
    v2_jet[6][j][1]  =  v2_jet[6][j][0];
    v2_jet[8][j][1]  =  v2_jet[8][j][0];
    v2_jet[10][j][1] =  v2_jet[10][j][0];
  }

  double bglevel, A_AS, V_2, V_1, V_3,alpha, beta, bglevelmax, bglevelmin, A_ASmax, A_ASmin, V_2max, V_1max, V_1min, V_3max, V_3min, V_2min, fitfunc0,fitmin0,fitmax0,histvalue, evalpt, error2, bglevel_err, V_1_err, V_2_err, V_3_err, A_AS_err, alpha_err, beta_err,bc, temp_bc, temp,temp_err, ymax, ymin, etamin1_val, etamax1_val, etamin2_val, etamax2_val, refvalue, referror, projmax,projmin;

  int etamin1, etamax1, etamin2, etamax2;

  TString in_name, plotname, outname, funcname, centlabel, datalabel, pTlabel,checkcanvasnameEta,checkcanvasnamePhi,fitcanvasname, diff_etacanvasname, diff_phicanvasname,zoomedcanvasname, alphaCentcanvasname,alphapTcanvasname, betapTcanvasname, betaCentcanvasname;

  //------------------------------------------
  //Fit function and generic parameter sets
  //------------------------------------------

 TF1 *gen_gaus = new TF1("gen_gaus",
			  " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))  \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))");
 
  gen_gaus->SetParName(0,"Bkg level");
  gen_gaus->SetParName(1,"V_{1}");
  gen_gaus->SetParName(2,"V_{2}");
  gen_gaus->SetParName(3,"V_{3}");
  gen_gaus->SetParName(4,"A_{AS}");
  gen_gaus->SetParName(5,"#alpha_{AS}");
  gen_gaus->SetParName(6,"#beta_{AS}");


  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
  
  for(int g = gstart; g<gend; g++){

    if(skip_pp==kTRUE&&(g==1||g==3||g==5||g==7||g==9||g==11)){continue;}

    
    switch(g){
    case 0: 
    
      if(Is_NonLeading == kFALSE){
	fdata[g] = new TFile("../me_correct/PbPb_Inclusive_Correlations.root","READ");
	fout[g] = new TFile("PbPb_Inclusive_Yield_and_Bkg.root", "RECREATE");
	datalabel = "Inclusive";
      }
      
      if(Is_NonLeading == kTRUE){
	fdata[g] = new TFile("../me_correct/PbPb_NonLeading_Correlations.root","READ");
	fout[g] = new TFile("PbPb_NonLeading_Yield_and_Bkg.root", "RECREATE");
	datalabel = "NonLeading";
      }
      break;
       
    case 1:
      //	fdata[g] = new TFile("pp_Inclusive120_2Dyield_and_SummedBkg_files.root","READ");
	fdata[g] = new TFile("../me_correct/pp_Inclusive_Correlations.root","READ");
	fout[g] = new TFile("pp_Inclusive_Yield_and_Bkg.root", "RECREATE");
	datalabel = "Inclusive";
           
      break;
     
    case 2:
      fdata[g] = new TFile("../me_correct/PbPb_SubLeading_Correlations.root","READ");
      fout[g] = new TFile("PbPb_SubLeading_Yield_and_Bkg.root", "RECREATE");
      datalabel = "SubLeading";
      break;
      
    case 3:
      //  fdata[g] = new TFile("pp_SubLeading50_2Dyield_and_SummedBkg_files.root","READ");
      fdata[g] = new TFile("../me_correct/pp_SubLeading_Correlations.root","READ");
      fout[g] = new TFile("pp_SubLeading_Yield_and_Bkg.root", "RECREATE");
      datalabel = "SubLeading";
      break;
    
    case 4:
      fdata[g] = new TFile("../me_correct/PbPb_Leading_Correlations.root","READ");
      fout[g] = new TFile("PbPb_Leading_Yield_and_Bkg.root", "RECREATE");
      datalabel = "Leading";
      break;

    case 5:
      //  fdata[g] = new TFile("pp_Leading120_2Dyield_and_SummedBkg_files.root","READ");
      fdata[g] = new TFile("../me_correct/pp_Leading_Correlations.root","READ");
      fout[g] = new TFile("pp_Leading_Yield_and_Bkg.root", "RECREATE");
      datalabel = "Leading";
      break;

    case 6:
      fdata[g] = new TFile("../Pelin_New_Output/PYTHIAHYDJET_Inclusive_Correlations.root","READ");
      fout[g] = new TFile("HYDJET_Inc120_2Dyield_and_NewBkg_files.root", "RECREATE");
      datalabel = "Inclusive";
      break;
       
    case 7:
      fdata[g] = new TFile("../Pelin_New_Output/PYTHIA_Inclusive_Correlations.root","READ");
      fout[g] = new TFile("PYTHIA_Inc120_2Dyield_and_NewBkg_files.root", "RECREATE");
      datalabel = "Inclusive";
      break;
     
    case 8:
      fdata[g] = new TFile("../Pelin_New_Output/PYTHIAHYDJET_SubLeading_Correlations.root","READ");
      fout[g] = new TFile("HYDJET_Subleading50_2Dyield_and_NewBkg_files.root", "RECREATE");
      datalabel = "SubLeading";
      break;
      
    case 9:
      fdata[g] = new TFile("../Pelin_New_Output/PYTHIA_SubLeading_Correlations.root","READ");
      fout[g] = new TFile("PYTHIA_Subleading50_2Dyield_and_NewBkg_files.root", "RECREATE");
      datalabel = "SubLeading";
      break;
    
    case 10:
      fdata[g] = new TFile("../Pelin_New_Output/PYTHIAHYDJET_Leading_Correlations.root","READ");
      fout[g] = new TFile("HYDJET_Leading120_2Dyield_and_NewBkg_files.root", "RECREATE");
      datalabel = "Leading";
      break;

    case 11:
      fdata[g] = new TFile("../Pelin_New_Output/PYTHIA_Leading_Correlations.root","READ");
      fout[g] = new TFile("PYTHIA_Leading120_2Dyield_and_NewBkg_files.root", "RECREATE");
      datalabel = "Leading";
      break;
   
    }

    
    //----------------------------------------------
    //     Start i and l loops -- set up canvases etc.
    //------------------------------------------------

    for(int l = 0; l<1; l++){


      if(l>0){ continue; } //no Gen/Reco for data
   

      fitcanvasname = "FitCanvas";
      fitcanvasname+= g;
      fitcanvasname+= l;
    
      cfit[g][l] = new TCanvas(fitcanvasname,"",0,0,1500,1600);
      cfit[g][l]->Divide(4,4,0.000001,0.000001);
    
     

      zoomedcanvasname = "ZoomedCanvas";
      zoomedcanvasname+= g;
      zoomedcanvasname+= l;


      czoomed[g][l] = new TCanvas(zoomedcanvasname,"",0,0,1500,1600);
      czoomed[g][l]->Divide(4,4,0.000001,0.000001);
    
      
      for(int i=0; i<4; i++){

	//----------------------------------------------------
	//  Start of j loop -- individual histo analysis begins
	//-----------------------------------------------------
	
	for (int j=0; j<4; j++){

	  cfit[g][l]->cd(4*i+j+1);

	  
	
	  in_name = make_name("Yield_",g,i,j,l,centlabel,pTlabel);


	  //******************************************//

	  //  SET BACKGROUND SAMPLING REGION HERE//

	  //******************************************//

	  projmax = 2.5;
	  projmin = 1.5;

	  
	  etamin1_val = -projmax+0.001;
	  etamax1_val = -projmin-0.001;
	  etamin2_val =  projmin+0.001;
	  etamax2_val =  projmax -0.001;
	    
	    

	  //******************************************//
	


	//-------------------------------------
	//  Get and assign all histograms 
	//-------------------------------------


	cout<<"                   "<<endl;
	cout<<"------------------------------------------------------------------------------------"<<endl;
	cout<<datalabel<<"    "<<in_name<<endl;
	cout<<"------------------------------------------------------------------------------------"<<endl;

	if(g==6||g==8||g==10){ in_name.ReplaceAll("Gen","PYTHIAHYDJET"); }
	if(g==7||g==9||g==11){ in_name.ReplaceAll("Gen","PYTHIA"); }

	cout<<in_name<<endl;

	
	
	yield[g][i][j][l] = (TH2D*)fdata[g]->Get(in_name)->Clone(in_name);

	if((g==2||g==4)&&i==3&&j==0){
	  TString Temp_name = in_name;
	  Temp_name.ReplaceAll("Yield_PbPb","Raw_Yield");
	  
	  yield[g][i][j][l] = (TH2D*)fdata[g]->Get(Temp_name)->Clone(in_name);
	  
	  Temp_name.ReplaceAll("Raw_Yield","Mixed_Event");
	  Temp_name.ReplaceAll("Cent50_Cent100","Cent30_Cent50");

	  TH2D *temp_me = (TH2D*)fdata[g]->Get(Temp_name)->Clone(Temp_name);
	  yield[g][i][j][l]->Divide(temp_me);
	  
	}

	llimiteta = yield[g][i][j][l]->GetXaxis()->FindBin(-etalim+0.001);
	rlimiteta = yield[g][i][j][l]->GetXaxis()->FindBin(etalim-0.001);
	
	nbins = (rlimiteta-llimiteta+1);


	TString yieldprojname = in_name;
	yieldprojname.ReplaceAll("Yield","Yield_proj");
	yield_proj[g][i][j][l] = (TH1D*)yield[g][i][j][l]->ProjectionY(yieldprojname,llimiteta,rlimiteta);
	yield_proj[g][i][j][l]->Scale(1./nbins);
	


	TString bkgname = in_name;
	bkgname.ReplaceAll("Yield","Summed_bkg");

	TString tempname = "temp_name";
	tempname+=g;
	tempname+=i;
	tempname+=j;
	tempname+=l;


	etamin1 = findbin(etamin1_val);
	etamax1 = findbin(etamax1_val);
	etamin2 = findbin(etamin2_val);
	etamax2 = findbin(etamax2_val);
 
	nbins = etamax1-etamin1+etamax2-etamin2+2;
	
	in_hist[g][i][j][l] = yield[g][i][j][l]->ProjectionY(bkgname,etamin1,etamax1);
	in_hist[g][i][j][l]->Add(yield[g][i][j][l]->ProjectionY(tempname,etamin2,etamax2));
	in_hist[g][i][j][l]->Scale(1./nbins);

	outname= in_name;
	outname.ReplaceAll("Yield","Fit_bkg");    
	out_hist[g][i][j][l] = (TH2D*)yield[g][i][j][l]->Clone(outname);
	
	
	TString bkgname2 = in_name;
	bkgname2.ReplaceAll("Yield","Summed2D_bkg");
	in_hist2[g][i][j][l] = (TH2D*)yield[g][i][j][l]->Clone(bkgname2);
	
	cout<<"got histos"<<endl;
	
	// Deal properly with bad bins:

	for(int k = 0; k<100; k++){

	  if(g==4&&i==2&&j==0){
	    cout<<k<<" "<<in_hist[g][i][j][l]->GetBinCenter(k)<<" "<<in_hist[g][i][j][l]->GetBinContent(k)<<" "<<temp_bc<<" "<<in_hist[g][i][j][l]->GetBinError(k)<<" "<<temp_err<<endl;
	  }
	 
	    temp_bc = in_hist[g][i][j][l]->GetBinContent(k-1);
	    temp_err = in_hist[g][i][j][l]->GetBinError(k-1);
	   
	 
	    if(in_hist[g][i][j][l]->GetBinContent(k)<0.2*temp_bc&&in_hist[g][i][j][l]->GetBinError(k)<0.5*temp_err){
	      cout<<"******Replacing bin content**********"<<endl;
	      cout<<in_hist[g][i][j][l]->GetBinContent(k)<<" "<<temp_bc<<" "<<in_hist[g][i][j][l]->GetBinError(k)<<" "<<temp_err<<endl;
	      in_hist[g][i][j][l]->SetBinContent(k,temp_bc);
	      in_hist[g][i][j][l]->SetBinError(k,temp_err);
	    
	    }
		
	}

	//-------------------------------------
	//      PbPb Fits
	//-------------------------------------


	//	cfit[g][l]->cd(4*i+j+1);
		

	bglevel= in_hist[g][i][j][l]->GetBinContent(in_hist[g][i][j][l]->FindBin(TMath::Pi()/2));
	A_AS =in_hist[g][i][j][l]->GetMaximum();
	gen_gaus->SetParameter(0,bglevel);
	gen_gaus->SetParLimits(0,0.,1.0);


	if(g==0||g==2||g==4){
	  v2_assoc[g][i][j][l] = find_v2(i,j);
	  V_2 = v2_assoc[g][i][j][l]*v2_jet[g][j][l];
	    
	  if(V_2<0){cout<<"v2_jet = "<<v2_jet[g][j][l]<<" v2_assoc = "<<v2_assoc[g][i][j][l]<<" V2 = "<<V_2<<" WHAT?????"<<endl; return -1;}
	}


	V_3= 0.;
	gen_gaus->FixParameter(3,0.);

	if(g==0||g==2||g==4){

	  switch(j){
	    case 0:
	      gen_gaus->FixParameter(1,-0.002);
	      break;
	    case 1:
	      gen_gaus->FixParameter(1,-0.004);
	      break;
	    case 2:
	      gen_gaus->FixParameter(1,-0.005);
	      break;
	    case 3:
	      gen_gaus->FixParameter(1,-0.005);
	      break;
	    }  
	  
	  gen_gaus->FixParameter(2,V_2);
	  gen_gaus->FixParameter(3,0.);
	  gen_gaus->SetParameter(5,0.4);
	  gen_gaus->SetParLimits(5,0.2,0.8);
	  gen_gaus->SetParameter(6,1.7);
	  gen_gaus->SetParLimits(6,1.1,2.);
	
	  if(g==2){
	    gen_gaus->SetParLimits(5,0.2,0.6);
  
	    switch(j){
	    case 0:
	      gen_gaus->FixParameter(1,0.002);
	      break;
	    case 1:
	      gen_gaus->FixParameter(1,0.004);
	      break;
	    case 2:
	      gen_gaus->FixParameter(1,0.005);
	      break;
	    case 3:
	      gen_gaus->FixParameter(1,0.005);
	      break;
	    }  

	   
	  }

	 
	 
	}

	if(g==0||g==2||g==4){

	
	  in_hist[g][i][j][l]->Fit("gen_gaus","","",1.51,4.);

	  
	  A_AS = gen_gaus->GetParameter(4);
	  alpha = gen_gaus->GetParameter(5);
	  beta  = gen_gaus->GetParameter(6);

	  A_AS_err = gen_gaus->GetParError(4);
	  alpha_err = gen_gaus->GetParError(5);
	  beta_err = gen_gaus->GetParError(6);

	  // gen_gaus->FixParameter(4,A_AS);
	  gen_gaus->FixParameter(5, alpha);  
	  gen_gaus->FixParameter(6, beta);

	  in_hist[g][i][j][l]->Fit("gen_gaus","","",-1.57,1.57);

	  bglevel = gen_gaus->GetParameter(0);
	  V_1     = gen_gaus->GetParameter(1);
	  V_2     = gen_gaus->GetParameter(2);

	  bglevel_err = gen_gaus->GetParError(0);
	  V_1_err = gen_gaus->GetParError(1);
	 
	  V_2_err = v2_jet_err[g][j][l]*v2_assoc[g][i][j][l];

	  //	  V_2_err = gen_gaus->GetParError(2);
	  
	  v2_jet_test[g][i][j][l] = V_2/v2_assoc[g][i][j][l];

	}

	//--------------------------------
	//    Fits for pp data/ PYTHIA
	//-------------------------------	  


	if(g==1||g==3||g==5||g==7||g==9||g==11){
	  gen_gaus->SetParameter(0,bglevel);
	  gen_gaus->ReleaseParameter(1);
	  gen_gaus->SetParameter(1,0.01);
	  gen_gaus->SetParLimits(1,-0.6,0.6);
	  gen_gaus->FixParameter(2,0.);	
	  gen_gaus->FixParameter(3,0.);
	  gen_gaus->SetParameter(4,A_AS);
	  gen_gaus->SetParameter(5,0.4);
	  gen_gaus->SetParLimits(5,0.1,2.0);
	  gen_gaus->ReleaseParameter(6);
	  gen_gaus->SetParameter(6,1.1);
	  gen_gaus->SetParLimits(6,1.,2.);


	  in_hist[g][i][j][l]->Fit("gen_gaus","","",-1.,4.5);
	  
	  bglevel_err = gen_gaus->GetParError(0);
	  V_1_err   = gen_gaus->GetParError(1);
	  V_2_err   = gen_gaus->GetParError(2);
	  V_3_err   = gen_gaus->GetParError(3);
	  A_AS_err  = gen_gaus->GetParError(4);
	  alpha_err = gen_gaus->GetParError(5);
	  beta_err  = gen_gaus->GetParError(6);



	  bglevel = gen_gaus->GetParameter(0);
	  V_1   = gen_gaus->GetParameter(1);
	  V_2   = gen_gaus->GetParameter(2);
	  V_3   = gen_gaus->GetParameter(3);
	  A_AS  = gen_gaus->GetParameter(4);
	  alpha = gen_gaus->GetParameter(5);
	  beta  = gen_gaus->GetParameter(6);

	}

	//------------------------------
	// Save parameters in arrays
	//---------------------------------



	switch(j){
	case 0: alpha_array_pT0[i] = alpha;
	  alpha_error_pT0[i] = alpha_err;
	  beta_array_pT0[i] = beta;
	  beta_error_pT0[i] = beta_err;
	  V2_array_pT0[i] = V_2;
	  V2_error_pT0[i] = V_2_err;
	  break;

	case 1: alpha_array_pT1[i] = alpha;
	  alpha_error_pT1[i] = alpha_err;
	  beta_array_pT1[i] = beta;
	  beta_error_pT1[i] = beta_err;
	  V2_array_pT1[i] = V_2;
	  V2_error_pT1[i] = V_2_err; 
	  break;

	case 2: alpha_array_pT2[i] = alpha;
	  alpha_error_pT2[i] = alpha_err;
	  beta_array_pT2[i] = beta;
	  beta_error_pT2[i] = beta_err;
	  V2_array_pT2[i] = V_2;
	  V2_error_pT2[i] = V_2_err;
	  break;

	case 3: alpha_array_pT3[i] = alpha;
	  alpha_error_pT3[i] = alpha_err;
	  beta_array_pT3[i] = beta;
	  beta_error_pT3[i] = beta_err;
	  V2_array_pT3[i] = V_2;
	  V2_error_pT3[i] = V_2_err;
	  break;
	}


	switch(i){
	case 0:alpha_array_Cent0[j] = alpha;
	  alpha_error_Cent0[j] = alpha_err;
	  beta_array_Cent0[j] = beta;
	  beta_error_Cent0[j] = beta_err;
	  V2_array_Cent0[j] = V_2;
	  V2_error_Cent0[j] = V_2_err;
	  break; 

	case 1:alpha_array_Cent1[j] = alpha;
	  alpha_error_Cent1[j] = alpha_err;
	  beta_array_Cent1[j] = beta;
	  beta_error_Cent1[j] = beta_err;
	  V2_array_Cent1[j] = V_2;
	  V2_error_Cent1[j] = V_2_err;
	  break;	  

	case 2:alpha_array_Cent2[j] = alpha;
	  alpha_error_Cent2[j] = alpha_err;
	  beta_array_Cent2[j] = beta;
	  beta_error_Cent2[j] = beta_err;
	  V2_array_Cent2[j] = V_2;
	  V2_error_Cent2[j] = V_2_err;
	  break;	  

	case 3:alpha_array_Cent3[j] = alpha;
	  alpha_error_Cent3[j] = alpha_err;
	  beta_array_Cent3[j] = beta;
	  alpha_error_Cent3[j] = beta_err;
	  V2_array_Cent3[j] = V_2;
	  V2_error_Cent3[j] = V_2_err;
	  break;	  

	case 4:alpha_array_Cent4[j] = alpha;
	  alpha_error_Cent4[j] = alpha_err;
	  beta_array_Cent4[j] = beta;
	  beta_error_Cent4[j] = beta_err;
	  V2_array_Cent4[j] = V_2;
	  V2_error_Cent4[j] = V_2_err;
	  break;
	}
	

	bglevelmax = bglevel + bglevel_err;
	bglevelmin = bglevel - bglevel_err;
      
	//	A_ASmax = A_AS + A_AS_err;
	//A_ASmin = A_AS - A_AS_err;
     

	A_ASmax = A_AS;
	A_ASmin = A_AS;
     
	V_1max = V_1 + V_1_err;
	V_1min = V_1 - V_1_err;
 
	V_2max = V_2 + V_2_err;
	V_2min = V_2 - V_2_err;

	V_3max =0.;
	V_3min = 0.;
	
	//----------------------------------------------------
	//This whole section is just about drawing fit curves
	//------------------------------------------------------

	funcname = "fitfunc";	funcname+=g;   funcname+=i;	funcname+=j;
	fitfunc[g][i][j][l] = new TF1( funcname,"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))+[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))+ TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",phimin,phimax);

	funcname = "fitmax";	funcname+=g;  funcname+=i;	funcname+=j;
	fitmax[g][i][j][l] = new TF1(funcname,"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))+[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))+ TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",phimin,phimax);

	funcname = "fitmin";	funcname+=g;  funcname+=i;	funcname+=j;
	fitmin[g][i][j][l] = new TF1( funcname,"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))+[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))+ TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",phimin,phimax);
      
	funcname = "funcmax";	funcname+=g;	funcname+=i;	funcname+=j;
	funcmax[g][i][j][l] = new TF1( funcname,"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))+[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))+ TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))+[7]",phimin,phimax);

	funcname = "funcmin";	funcname+=g; 	funcname+=i;	funcname+=j;
	funcmin[g][i][j][l] = new TF1( funcname,"[0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))+[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))+ TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))+[7]",phimin,phimax);

	/*
	cout<<bglevel<<" "<<bglevelmax<<" "<<bglevelmin<<endl;
	cout<<V_1<<" "<<V_1max<<" "<<V_1min<<endl;
	cout<<V_2<<" "<<V_2max<<" "<<V_2min<<endl;
	cout<<V_3<<" "<<V_3max<<" "<<V_3min<<endl;
	cout<<A_AS<<" "<<A_ASmax<<" "<<A_ASmin<<endl;
	*/

	fitfunc[g][i][j][l]->SetParameter(0,bglevel);
	fitfunc[g][i][j][l]->SetParameter(1,V_1);
	fitfunc[g][i][j][l]->SetParameter(2,V_2);
	fitfunc[g][i][j][l]->SetParameter(3,V_3);
	fitfunc[g][i][j][l]->SetParameter(4,A_AS);
	fitfunc[g][i][j][l]->SetParameter(5,alpha);
	fitfunc[g][i][j][l]->SetParameter(6,beta);

	fitmax[g][i][j][l]->SetParameter(0,bglevelmax);
	fitmax[g][i][j][l]->SetParameter(1,V_1max);
	fitmax[g][i][j][l]->SetParameter(2,V_2max);
	fitmax[g][i][j][l]->SetParameter(3,V_3max);
	fitmax[g][i][j][l]->SetParameter(4,A_ASmax);
	fitmax[g][i][j][l]->SetParameter(5,alpha);
	fitmax[g][i][j][l]->SetParameter(6,beta);

	fitmin[g][i][j][l]->SetParameter(0,bglevelmin);
	fitmin[g][i][j][l]->SetParameter(1,V_1min);
	fitmin[g][i][j][l]->SetParameter(2,V_2min);
	fitmin[g][i][j][l]->SetParameter(3,V_3min);
	fitmin[g][i][j][l]->SetParameter(4,A_ASmin);
	fitmin[g][i][j][l]->SetParameter(5,alpha);
	fitmin[g][i][j][l]->SetParameter(6,beta);

	fitmax[g][i][j][l]->SetLineColor(kBlack);
	fitmin[g][i][j][l]->SetLineColor(kBlue);  

     	in_hist[g][i][j][l]->SetLineColor(2);
	in_hist[g][i][j][l]->SetMarkerColor(2);
	in_hist[g][i][j][l]->SetMarkerSize(1);
	in_hist[g][i][j][l]->SetMarkerStyle(10);
	in_hist[g][i][j][l]->GetXaxis()->SetLabelSize(0.05);
	in_hist[g][i][j][l]->GetYaxis()->SetLabelSize(0.05);
	in_hist[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	in_hist[g][i][j][l]->GetXaxis()->SetTitleSize(0.05);
      

	fitfunc[g][i][j][l]->SetLineColor(kBlack);

	fitfunc0 = fitfunc[g][i][j][l]->Eval(0,0,0);
	fitmax0 = fitmax[g][i][j][l]->Eval(0,0,0);
	fitmin0 = fitmin[g][i][j][l]->Eval(0,0,0);

	error[g][i][j][l] = TMath::Max(TMath::Abs(fitmax0 - fitfunc0),TMath::Abs(fitfunc0-fitmin0));
    
	error2 = error[g][i][j][l];


	//	cout<<"********Error is: "<< error2<<" ******************"<<endl;

	funcmax[g][i][j][l]->SetParameter(0,bglevel);
	funcmax[g][i][j][l]->SetParameter(1,V_1);
	funcmax[g][i][j][l]->SetParameter(2,V_2);
	funcmax[g][i][j][l]->SetParameter(3,V_3);
	funcmax[g][i][j][l]->SetParameter(4,A_AS);
	funcmax[g][i][j][l]->SetParameter(5,alpha);
	funcmax[g][i][j][l]->SetParameter(6,beta);
	funcmax[g][i][j][l]->SetParameter(7,error2);
	funcmax[g][i][j][l]->SetLineColor(kBlue);

	funcmin[g][i][j][l]->SetParameter(0,bglevel);
	funcmin[g][i][j][l]->SetParameter(1,V_1);
	funcmin[g][i][j][l]->SetParameter(2,V_2);	
	funcmin[g][i][j][l]->SetParameter(3,V_3);
	funcmin[g][i][j][l]->SetParameter(4,A_AS);
	funcmin[g][i][j][l]->SetParameter(5,alpha);
	funcmin[g][i][j][l]->SetParameter(6,beta);     
	funcmin[g][i][j][l]->SetParameter(7,-error2);      
	funcmin[g][i][j][l]->SetLineColor(kBlue);



	//---------------------------------------------
	//    Fill "New Bkg"
	//---------------------------------------------

	//

	for(int k = 0; k<101;k++){
   
	  evalpt = in_hist[g][i][j][l]->GetBinCenter(k);
	  histvalue = (fitfunc[g][i][j][l]->Eval(evalpt));

	  refvalue = in_hist[g][i][j][l]->GetBinContent(k);
	  referror = in_hist[g][i][j][l]->GetBinError(k);
	

	  for(int m = 0; m<101; m++){

   
	    in_hist2[g][i][j][l]->SetBinContent(m,k,refvalue);
	    in_hist2[g][i][j][l]->SetBinError(m,k,referror*sqrt(nbins));
	 
	    out_hist[g][i][j][l]->SetBinContent(m,k,histvalue);
	    out_hist[g][i][j][l]->SetBinError(m,k,error2*sqrt(nbins));

	  }
	}

	//--------------------------------------------
	//   Draw fit plots
	//---------------------------------------------

	cfit[g][l]->cd(4*i+j+1);
	
	TString outhistproj_name = in_name;
	outhistproj_name.ReplaceAll("Yield","New1d_bkg");

	out_hist1[g][i][j][l] = out_hist[g][i][j][l]->ProjectionY(outhistproj_name,50,50);  //okay to be single bin, because we made it that way!!
	

	out_hist1[g][i][j][l]->SetLineColor(1);
	out_hist1[g][i][j][l]->SetMarkerColor(1);
	out_hist1[g][i][j][l]->SetMarkerSize(1);
	out_hist1[g][i][j][l]->SetMarkerStyle(10);

	switch(i){
	case 0: ymax =yield_proj[g][i][j][l]->GetMaximum()+0.008;  break;
	case 1: ymax =yield_proj[g][i][j][l]->GetMaximum()+0.006;  break;
	case 2: ymax =yield_proj[g][i][j][l]->GetMaximum()+0.005;  break;
	case 3: ymax =yield_proj[g][i][j][l]->GetMaximum()+0.012;  break;

	}
	
	ymin =yield_proj[g][i][j][l]->GetMinimum()-0.002;

	if(g<6){in_hist_rebin[g][i][j][l] = Rebin_dPhi_full(in_hist[g][i][j][l]);
	}else{
	  in_hist_rebin[g][i][j][l] = Rebin_dPhi_full_old(in_hist[g][i][j][l]);
	}

	in_hist_rebin[g][i][j][l]->Rebin(2);
	in_hist_rebin[g][i][j][l]->Scale(1/2.);



	if((in_hist_rebin[g][i][j][l]->GetMaximum()+.004)>ymax){ymax = in_hist_rebin[g][i][j][l]->GetMaximum()+.004;}
	//	in_hist_rebin[g][i][j][l]->SetAxisRange(-1.5,4.8);

	cfit[g][l]->cd(4*i+j+1);

	in_hist_rebin[g][i][j][l]->SetMaximum(ymax);
	in_hist_rebin[g][i][j][l]->SetMinimum(ymin);
	in_hist_rebin[g][i][j][l]->Draw();
	
	
	if(g<6){	yield_rebin[g][i][j][l]= Rebin_dPhi_full(yield_proj[g][i][j][l]);
	} else {
	  yield_rebin[g][i][j][l]= Rebin_dPhi_full_old(yield_proj[g][i][j][l]);
	}


	TString yield_rebin_name = in_name;
	yield_rebin_name += "_1D_rebin";
	yield_rebin_name +=i;
	yield_rebin_name+=j;
	yield_rebin[g][i][j][l]->SetName(yield_rebin_name);
	


	yield_rebin[g][i][j][l]->SetMarkerColor(kMagenta+2);

	yield_rebin[g][i][j][l]->SetLineColor(kMagenta+2);
	yield_rebin[g][i][j][l]->SetMarkerStyle(10);

	yield_rebin[g][i][j][l]->SetMaximum(ymax);
	yield_rebin[g][i][j][l]->SetMinimum(ymin);

	funcmax[g][i][j][l]->Draw("same");
	funcmin[g][i][j][l]->Draw("same");

	yield_rebin[g][i][j][l]->Draw("same");
	in_hist_rebin[g][i][j][l]->Draw("same");
	
	fitfunc[g][i][j][l]->Draw("same");

	if(g<6){	out_hist_rebin[g][i][j][l] = (TH1D*)Rebin_dPhi_full(out_hist1[g][i][j][l]);
	} else{
	  out_hist_rebin[g][i][j][l] = (TH1D*)Rebin_dPhi_full_old(out_hist1[g][i][j][l]);
	}

	TString out_hist_rebin_name = outname;
	out_hist_rebin_name += "rebin";
	out_hist_rebin_name += g;
	out_hist_rebin_name += i;
	out_hist_rebin_name+=j;

	//	out_hist_rebin[g][i][j][l]->Draw("same");

	
	TLegend *lfit = new TLegend(0.65,0.7,0.88,0.94);
	lfit->SetFillColor(kWhite);
	lfit->SetLineColor(kWhite);
	lfit->AddEntry(yield_rebin[g][i][j][l],"Yield");
	lfit->AddEntry(in_hist_rebin[g][i][j][l],"Sum. Bkg");
	lfit->AddEntry(fitfunc[g][i][j][l],"Fit Bkg");
	lfit->AddEntry(funcmax[g][i][j][l],"Fit Error","l");
	lfit->SetTextSize(ts2-0.01);
	lfit->Draw("same");

	TPaveText *pave;
	if(j==0){
	  pave = new TPaveText(0.18,0.8,0.45,0.9,"NDC");
	}else{
	  pave = new TPaveText(0.05,0.8,0.45,0.9,"NDC");
	}  
	pave->SetName("pave");
	pave->SetFillColor(0);
	pave->SetTextAlign(11);
	pave->AddText(centlabel);
	pave->AddText(pTlabel);
	pave->SetTextSize(ts2-0.01);
	pave->Draw("same");
	

	cout<<"******"<<centlabel<<"********"<<endl;
	czoomed[g][l]->cd(4*i+j+1);

      	
	TString in_name_zoomed = in_name;
	in_name_zoomed+= "_zoomed";
	in_name_zoomed+=g;
	in_name_zoomed+=i;
	in_name_zoomed+=j;
	in_name_zoomed+=l;


	in_hist_zoomed[g][i][j][l]= (TH1D*)in_hist_rebin[g][i][j][l]->Clone(in_name_zoomed);
	in_hist_zoomed[g][i][j][l]->SetAxisRange(-1.5,1.5,"x");

	in_hist_zoomed[g][i][j][l]->Draw();



	funcmax[g][i][j][l]->Draw("same");
	funcmin[g][i][j][l]->Draw("same");

	yield_rebin[g][i][j][l]->Draw("same");
	
	in_hist_zoomed[g][i][j][l]->Draw("same");
	fitfunc[g][i][j][l]->Draw("same");
	//	out_hist_rebin[g][i][j][l]->Draw("same");


	lfit->Draw("same");
	pave->Draw("same");

	//-------------------------------
	//And finally produce subtracted yield save it all!
	//------------------------------ 
      
	TString resultname = in_name;
	resultname.ReplaceAll("Yield","Result");

	result[g][i][j][l] = (TH2D*)yield[g][i][j][l]->Clone(resultname);
	result[g][i][j][l]->Add(out_hist[g][i][j][l],-1.);


     
	TString resultname2 = in_name;
	resultname2.ReplaceAll("Yield","SummedResult");
      
     
	result2[g][i][j][l] = (TH2D*)yield[g][i][j][l]->Clone(resultname2);

	result2[g][i][j][l]->Add(in_hist2[g][i][j][l],-1.);

	in_hist[g][i][j][l]->Write();
	out_hist[g][i][j][l]->Write();
	result[g][i][j][l]->Write();
	in_hist2[g][i][j][l]->Write();
	result2[g][i][j][l]->Write();
	fitfunc[g][i][j][l]->Write();



	
	//	drawlabels(g,i,j);

      
	
	} //closes the [j] (centrality) loop

      } //closes the [i] (pT) loop


	
      TString fitsavename = in_name;

      fitsavename.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      fitsavename.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt8_TrkPt999",datalabel);
      fitsavename.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      if(g==6||g==8||g==10){fitsavename.ReplaceAll("Yield","Fits_HYDJET");}
      if(g==7||g==9||g==11){fitsavename.ReplaceAll("Yield","Fits_PYTHIA");}
      fitsavename.ReplaceAll("Yield","Fits");
	
      fitsavename +=".pdf";
      cfit[g][l]->SaveAs(fitsavename);
     
      fitsavename.ReplaceAll(".pdf",".png");
      cfit[g][l]->SaveAs(fitsavename);
    
      TString zoomedsavename = in_name;
      zoomedsavename.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      zoomedsavename.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt8_TrkPt999",datalabel);
      zoomedsavename.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      if(g==6||g==8||g==10){zoomedsavename.ReplaceAll("Yield","Zoomed_HYDJET");}
      if(g==7||g==9||g==11){zoomedsavename.ReplaceAll("Yield","Zoomed_PYTHIA");}
      zoomedsavename.ReplaceAll("Yield","Zoomed");

      zoomedsavename +=".pdf";
      czoomed[g][l]->SaveAs(zoomedsavename);
      zoomedsavename.ReplaceAll(".pdf",".png");
      czoomed[g][l]->SaveAs(zoomedsavename);
            
    

    }
  }

  for(int g = gstart; g<gend; g++){
    
    if(skip_pp==kTRUE&&(g==1||g==3||g==5||g==7||g==9||g==11)){continue;}
    

    if(g<2){
      datalabel = "Inclusive";
    }else if(g>1&&g<4){
      datalabel = "SubLeading";
    }else if(g>3){
      datalabel = "Leading";
    }

    for(int l = 0; l<1; l++){
   

 
      TString  rawcanvasnameEta = "RawCanvasEta";
      rawcanvasnameEta+= g;
      rawcanvasnameEta+= l;
      craweta[g][l] = new TCanvas(rawcanvasnameEta," ",10,10,1500,1600);
      craweta[g][l]->Divide(4,4,0.001,0.001);

         
      for(Int_t i = 0; i<4; i++){

	for(Int_t j = 0; j<4; j++){

	  craweta[g][l]->cd(4*i+j+1);

	  in_name = make_name("Yield_",g,i,j,l,centlabel,pTlabel);

	  int phileft = yield[g][i][j][l]->GetYaxis()->FindBin(-.999);
	  int phiright = yield[g][i][j][l]->GetYaxis()->FindBin(.999);

	  
	  TString rawetaname = in_name;
	  rawetaname.ReplaceAll("Yield_","RawEta_");

	  cout<<rawetaname<<endl;

	  raw_eta[g][i][j][l] = (TH1D*)yield[g][i][j][l]->ProjectionX(rawetaname,phileft,phiright);

	  raw_eta[g][i][j][l]->SetAxisRange(-2.999,2.999);


	
	  float min = raw_eta[g][i][j][l]->GetMinimum()-0.1;
	  float max = raw_eta[g][i][j][l]->GetMaximum()+0.3;

	  raw_eta[g][i][j][l]->SetMinimum(min);
	  raw_eta[g][i][j][l]->SetMaximum(max);

	  raw_eta[g][i][j][l]->SetMarkerColor(kBlack);
	  raw_eta[g][i][j][l]->SetLineColor(kBlack);
	  raw_eta[g][i][j][l]->SetMarkerStyle(10);

	  
	  Int_t loc = 4*i+j+1;

	  craweta[g][l]->cd(4*i+j+1);

	  cout<<g<<" "<<i<<" "<<j<<" "<<4*i+j+1<<" "<<loc<<" "<<raw_eta[g][i][j][l]->GetBinContent(50)<<endl;
	  
	  raw_eta[g][i][j][l]->Draw();

	  TLine *l1 = new TLine(-2.5,min,-2.5,min+(max-min)/3);
	  l1->SetLineWidth(3);
	  l1->SetLineColor(kRed);
	  l1->Draw("same");

	  TLine *l2 = new TLine(2.5,min,2.5,min+(max-min)/3);
	  l2->SetLineWidth(3);
	  l2->SetLineColor(kRed);
	  l2->Draw("same");

	  TLine *l3 = new TLine(-2.0,min,-2.0,min+(max-min)/3);
	  l3->SetLineWidth(3);
	  l3->SetLineColor(kViolet);
	  l3->Draw("same");

	  TLine *l4 = new TLine(2.0,min,2.0,min+(max-min)/3);
	  l4->SetLineWidth(3);
	  l4->SetLineColor(kViolet);
	  l4->Draw("same");

	  TLine *l5 = new TLine(-1.5,min,-1.5,min+(max-min)/3);
	  l5->SetLineWidth(3);
	  l5->SetLineColor(kGreen);
	  l5->Draw("same");

	  TLine *l6 = new TLine(1.5,min,1.5,min+(max-min)/3);
	  l6->SetLineWidth(3);
	  l6->SetLineColor(kGreen);
	  l6->Draw("same");


	  if(i==0&&j==3){
	    TLegend *legend = new TLegend(0.55,0.6,0.95,0.95,"","NDC");
	    legend->AddEntry(l1,"|#Delta#eta|=2.5","l");
	    legend->AddEntry(l3,"|#Delta#eta|=2.0","l");
	    legend->AddEntry(l5,"|#Delta#eta|=1.5","l");
	    legend->SetTextSize(0.08);
	    legend->SetLineWidth(0);
	    legend->SetLineColor(0);
	    legend->SetFillColor(0);
	    legend->Draw("same");

	  }


	  raw_eta[g][i][j][l]->Draw("same");
	  TPaveText *pave;
	  if(j==0){
	    pave = new TPaveText(0.18,0.8,0.45,0.9,"NDC");
	  }else{
	    pave = new TPaveText(0.05,0.8,0.45,0.9,"NDC");
	  }  
	  pave->SetName("pave");
	  pave->SetFillColor(0);
	  pave->SetTextAlign(11);
	  pave->AddText(centlabel);
	  pave->AddText(pTlabel);
	  pave->SetTextSize(ts2-0.01);
	  pave->Draw("same");


	  craweta[g][l]->Update();
	}//j
      
      }//i
     
    
      TString rawetasavename = in_name;
      rawetasavename.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      rawetasavename.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      rawetasavename.ReplaceAll("Yield","Raweta");

      rawetasavename +=".pdf";
      cout<<rawetasavename<<endl;
      craweta[g][l]->Update();
      craweta[g][l]->SaveAs(rawetasavename);

      rawetasavename.ReplaceAll(".pdf",".png");
      craweta[g][l]->SaveAs(rawetasavename);



    }

  }

  for(int g = gstart; g<gend; g++){
    
    if(skip_pp==kTRUE&&(g==1||g==3||g==5||g==7||g==9||g==11)){continue;}
    

    if(g<2){
      datalabel = "Inclusive";
    }else if(g>1&&g<4){
      datalabel = "SubLeading";
    }else if(g>3){
      datalabel = "Leading";
    }

    for(int l = 0; l<1; l++){
   
 

      diff_etacanvasname = "Diff_EtaCanvas";
      diff_etacanvasname+= g;
      diff_etacanvasname+= l;
    
      cdiff_eta[g][l] = new TCanvas(diff_etacanvasname,"",0,0,1500,1600);
      cdiff_eta[g][l]->Divide(4,4,0.000,0.000);
    
      diff_phicanvasname = "Diff_PhiCanvas";
      diff_phicanvasname+= g;
      diff_phicanvasname+= l;
    
      cdiff_phi[g][l] = new TCanvas(diff_phicanvasname,"",0,0,1500,1600);
      cdiff_phi[g][l]->Divide(4,4,0.000001,0.000001);
  
            
      for(int i = 0; i<4; i++){
	for(int j = 0; j<4; j++){

	  in_name = make_name("Yield_",g,i,j,l,centlabel,pTlabel);


	  //  Difference check plots
	  cdiff_phi[g][l]->cd(4*i+j+1);

	  TString testname = "test_fit_deviation";
	  testname+=g;
	  testname+=i;
	  testname+=j;
	  testname+=l;


	  test_diff[g][i][j][l] =(TH2D*)result[g][i][j][l]->Clone(testname);
	  test_diff[g][i][j][l]->GetXaxis()->SetRangeUser(-3.,3.);
	  test_diff[g][i][j][l]->GetYaxis()->SetRangeUser(-1.5,1.5);

	
	  llimitphi = test_diff[g][i][j][l]->GetYaxis()->FindBin(-1.5+0.001);
	  rlimitphi = test_diff[g][i][j][l]->GetYaxis()->FindBin(1.5-0.001);


	  llimiteta = test_diff[g][i][j][l]->GetXaxis()->FindBin(-1.5+0.001);
	  rlimiteta = test_diff[g][i][j][l]->GetXaxis()->FindBin(1.5-0.001);

	  
	  for(int k = llimiteta; k<rlimiteta+1; k++){
	    for(int m = 0; m<101; m++){    
	    
	      test_diff[g][i][j][l]->SetBinContent(k,m,0.);
	      test_diff[g][i][j][l]->SetBinError(k,m,0.);
	    }
	  }



	  float dx_eta = test_diff[g][i][j][l]->GetXaxis()->GetBinWidth(1);
	  float dx_phi = test_diff[g][i][j][l]->GetYaxis()->GetBinWidth(1);
	  test_diff[g][i][j][l]->Scale(1./dx_eta/dx_phi);

	  TString test_diff2name = testname+="phi";
	  TString test_diff3name = testname+="phi2";
	  TString test_diff4name = testname+="in_hist_clone";
	  TString test_diff5name = testname+="eta";

	  test_diff2[g][i][j][l] = (TH1D*)test_diff[g][i][j][l]->ProjectionY(test_diff2name,etamin1,etamax1);
	  test_diff3[g][i][j][l] = (TH1D*)test_diff[g][i][j][l]->ProjectionY(test_diff3name,etamin2,etamax2);

	  test_diff2[g][i][j][l]->Add(test_diff3[g][i][j][l]);

	  test_diff2[g][i][j][l]->Scale(1./nbins);
	
	  test_diff2[g][i][j][l]->SetAxisRange(-1.5,1.5,"x");

	  test_diff2[g][i][j][l]->SetMinimum(-0.5);
	  test_diff2[g][i][j][l]->SetMaximum(0.5);

	  test_diff2[g][i][j][l]->Draw();

	  test_diff2[g][i][j][l]->Fit("fit0","","",-1.5,1.5);
	


	  test_diff4[g][i][j][l] =(TH1D*) in_hist[g][i][j][l]->Clone(test_diff4name);
	
	  test_diff4[g][i][j][l]->Add(out_hist1[g][i][j][l],-1.);
	
	  test_diff4[g][i][j][l]->Scale(1./dx_phi/dx_eta);

	  test_diff4[g][i][j][l]->Fit("fit0","","",-1.5,1.5);
	
	  test_diff2[g][i][j][l]->Draw();
	


	  TLine *linePhi = new TLine(-1.5,0.,1.5,0.);
	  linePhi->Draw("same");


	  test_diff4[g][i][j][l]->Draw("same");





	  cdiff_eta[g][l]->cd(4*i+j+1);
	
	  test_diff5[g][i][j][l] = (TH1D*)test_diff[g][i][j][l]->ProjectionX(test_diff5name,llimitphi,rlimitphi);

	  nbins = rlimitphi-llimitphi+1;

	  test_diff5[g][i][j][l]->Scale(1./nbins);

	
	  test_diff5[g][i][j][l]->SetAxisRange(-3.,3.,"x");

	  test_diff5[g][i][j][l]->SetMinimum(-0.5);
	  test_diff5[g][i][j][l]->SetMaximum(0.5);

	  //	test_diff5[g][i][j][l]->Draw();

	  sideband_hist[g][i][j][l] = (TH1D*)Rebin_dEta3(test_diff5[g][i][j][l]);

	  sideband_hist[g][i][j][l]->Fit("fit0","","q0",-2.99,-1.51);

	  yleft[g][i][j][l] = fit0->GetParameter(0);
	
	  sideband_hist[g][i][j][l]->Fit("fit0","","q0",1.51,2.99);

	  yright[g][i][j][l] = fit0->GetParameter(0);





	  sideband_hist[g][i][j][l]->SetMarkerColor(kBlack);
	  sideband_hist[g][i][j][l]->SetMarkerStyle(10);
	  sideband_hist[g][i][j][l]->SetLineColor(kBlack);

  
      
	  sideband_hist[g][i][j][l]->GetYaxis()->SetLabelSize(ts);
	  if(i==3){	sideband_hist[g][i][j][l]->GetYaxis()->SetLabelSize(ts2);}
	  sideband_hist[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	  sideband_hist[g][i][j][l]->GetXaxis()->SetTitle("#Delta#eta");
	  sideband_hist[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	  sideband_hist[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	  sideband_hist[g][i][j][l]->GetYaxis()->SetTitle("1/N_{jet} dN/(d#Delta#eta dp_{T})");
	  sideband_hist[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	  sideband_hist[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);

 
	  sideband_hist[g][i][j][l]->GetXaxis()->CenterTitle();
	  sideband_hist[g][i][j][l]->GetYaxis()->CenterTitle();
	  if(i<3){
	    sideband_hist[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	  }
	  if(j>0){
	    sideband_hist[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	    sideband_hist[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	  }



	  sideband_hist[g][i][j][l]->SetMinimum(-0.3);

	  sideband_hist[g][i][j][l]->Draw();
	
	  TLine *lineleft = new TLine(-3.0,yleft[g][i][j][l],-1.5,yleft[g][i][j][l]);
	  lineleft->SetLineColor(kRed);
	  lineleft->SetLineWidth(3);
	  lineleft->Draw();
	
	  TLine *lineright = new TLine(1.5,yright[g][i][j][l],3.0,yright[g][i][j][l]);
	  lineright->SetLineColor(kViolet);
	  lineright->SetLineWidth(3);
	  lineright->Draw();

	  drawlabels(g,i,j);

	  TLine *lineEta = new TLine(-3.,0.,3.,0.);
	  lineEta->SetLineStyle(2);
	  lineEta->Draw("same");



	}//j
      
      }//i
     
    
      TString diffsavename = in_name;
      diffsavename.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      diffsavename.ReplaceAll("Yield","Phi_EdgeFit");
	
      diffsavename +=".pdf";
      cdiff_phi[g][l]->SaveAs(diffsavename);
      diffsavename.ReplaceAll(".pdf",".png");
      cdiff_phi[g][l]->SaveAs(diffsavename);
            
    

      diffsavename.ReplaceAll("Phi","Eta");
      cdiff_eta[g][l]->SaveAs(diffsavename);

      diffsavename.ReplaceAll(".png",".pdf");
      cdiff_eta[g][l]->SaveAs(diffsavename);
    
      cdiff_eta[g][l]->Close();
      cdiff_phi[g][l]->Close();
     	 
      
    } // l   
    fout[g]->Close();  
    cout<<"done file "<<g<<endl;
  } //Closes the [g] (files) loop

  /*

    for(int g = 0; g<6; g++){
    for(int i = 0; i<4; i++){
    for(int j = 0; j<4; j++){

    cout<<g<<" "<<i<<" "<<j<<" "<<yleft[g][i][j][0]<<" "<<yright[g][i][j][0]<<endl;
    }
    }
    }
  */

  return 0;
}



Int_t findbin(double x) {
  int out_val = (x+5)*10.+1.;
    return out_val;
}


Float_t find_v2(int in_pT_i, int in_Cent_j){
 
  TGraphErrors *v2_dist[8];
  float v2_assoc_temp[8];
  
  float v2_assoc;
  int first_file,last_file;
  
  double pT_bins[5] = {1.5,2.5,3.5,6.,12.};
  double Cent_bins[4] = {75.,40.,20.,5.};
  

  float eval_pT = pT_bins[in_pT_i];
  
  vector<TString>* v2_by_Cent = new vector<TString>;
 
  v2_by_Cent->push_back("v2_pt_ep_cen0_10_eta08.txt");
  v2_by_Cent->push_back("v2_pt_ep_cen10_20_eta08.txt");
  v2_by_Cent->push_back("v2_pt_ep_cen20_30_eta08.txt");
  v2_by_Cent->push_back("v2_pt_ep_cen30_40_eta08.txt");
  v2_by_Cent->push_back("v2_pt_ep_cen40_50_eta08.txt");
  v2_by_Cent->push_back("v2_pt_ep_cen50_60_eta08.txt");
  v2_by_Cent->push_back("v2_pt_ep_cen60_70_eta08.txt");
  v2_by_Cent->push_back("v2_pt_ep_cen70_80_eta08.txt");

  
  vector<float> pt;
  vector<float> v2;
  vector<float> stat;
  vector<float> sys;
  vector<float> errtot;

  string item;
  string line;

  if(in_Cent_j==0){
    first_file = 5;
    last_file = 7;
  }
  if(in_Cent_j==1){
    first_file = 3;
    last_file = 4;
  }
  if(in_Cent_j==2){
    first_file = 1;
    last_file = 2;
 }
  if(in_Cent_j==3){
    first_file = 0;
    last_file = 0;
  }

  for(int this_file = first_file; this_file<=last_file; this_file++){
 

  ifstream inFile ( v2_by_Cent->at(this_file) );   
    
    pt.clear();
    v2.clear();
    stat.clear();
    sys.clear();
    errtot.clear();

    int linenum = 0;
  
    // read data
    // ---------
    while (getline (inFile, line) )  {
      // cout << "\nLine #" << linenum << ":" << endl;
      istringstream linestream(line);
      // cout << TString(line[0]) << endl;
      // if ( TString(line[0]).IsDigit() ) {
      // 	cout << "begins with a digit" << endl;
      // } else{
      // 	cout << "begins with something else" << endl;
      // }
      if ( !(TString(line[0]).IsDigit()) ){ continue; }

      if ( !getline (linestream, item, '\t') ) { cerr << "PROBLEM." << endl; return -1; }
      pt.push_back( atof(item.data()) ) ;

      if ( !getline (linestream, item, '\t') ) { cerr << "PROBLEM." << endl; return -1; }
      v2.push_back( atof(item.data()) ) ;

      if ( !getline (linestream, item, '\t') ) { cerr << "PROBLEM." << endl; return -1; }
      stat.push_back( atof(item.data()) ) ;

      if ( !getline (linestream, item, '\t') ) { cerr << "PROBLEM." << endl; return -1; }
      sys.push_back( atof(item.data()) ) ;

      errtot.push_back( TMath::Sqrt( stat.back()*stat.back() + sys.back()*sys.back() ) );


    
      /*
      cout << " pt = " << pt.at(linenum)
       	   << " v2 = " << v2.at(linenum)
       	   << " stat = " << stat.at(linenum)
       	   << " sys = " << sys.at(linenum)
       	   << " tot = " << errtot.at(linenum)
       	   << endl;
      linenum++;
      */
    }
    
    inFile.close();


    float * fpt = &(pt.at(0));
    float * fv2 = &(v2.at(0));
  

    v2_dist[this_file] = new TGraphErrors( pt.size(),
					fpt, fv2,
					0, (float *) &(errtot.at(0)) ); 
    /*
    TString v2test_name = "v2_Test_Plot";
    v2test_name += in_pT_i;
    v2test_name += in_Cent_j;
    v2test_name += this_file;
    v2test_name+= ".pdf";

    TCanvas *c1 = new TCanvas("c1");
    v2_dist[this_file]->Draw("ALPE");
    c1->SaveAs(v2test_name);
*/

    v2_assoc_temp[this_file] = v2_dist[this_file]->Eval(eval_pT);
    

  }



  if(in_Cent_j==0){
    v2_assoc = (v2_assoc_temp[5]+v2_assoc_temp[6]+v2_assoc_temp[7])/3.;
    //  cout<<eval_pT<<"     "<<v2_assoc_temp[5]<<"      "<<v2_assoc_temp[6]<<"     "<<v2_assoc_temp[7]<<endl;
  }
 if(in_Cent_j==1){
    v2_assoc = (v2_assoc_temp[3]+v2_assoc_temp[4])/2.;
    //  cout<<eval_pT<<"     "<<v2_assoc_temp[3]<<"     "<<v2_assoc_temp[4]<<endl;
  }
  if(in_Cent_j==2){
    v2_assoc = (v2_assoc_temp[1]+v2_assoc_temp[2])/2.;
    //  cout<<eval_pT<<"     "<<v2_assoc_temp[1]<<v2_assoc_temp[2]<<endl;
  }
  if(in_Cent_j==3){
    v2_assoc = v2_assoc_temp[0];
    // cout<<eval_pT<<"     "<<v2_assoc_temp[0]<<endl;
  }

  // cout<<in_Cent_j<<"    v2_assoc = "<<v2_assoc<<endl;


  return v2_assoc;
}





 


