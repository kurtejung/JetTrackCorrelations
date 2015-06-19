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


Float_t find_v2(int in_pt_i, int in_Cent_j);
TH1D *Rebin_dPhi(TH1D* hold);
TH1D *Rebin_dEta(TH1D* hold);
TString make_name(TString stem, int g, int i, int j, int l, TString &centlabel, TString &pTlabel);

Int_t bg_fit_free(int gstart = 0, int gend = 6,bool skip_pp = kTRUE, bool Is_NonLeading = kFALSE)
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
  Float_t v1_tot[12][5][4][2];
  Float_t v3_tot[12][5][4][2];
  Float_t v2_jet_test[12][5][4][2];
  Float_t v2_jet_err[12][4][2];

  TH1D *raw_eta[12][5][4][2];

  TH2D *yield[12][5][4][2];
  TH1D *yield_proj[12][5][4][2];
 
  TH1D *yield_rebin[12][5][4][2];
  TH1D *in_hist_test[12][5][4][2];
  TH1D *in_hist[12][5][4][2];
  TH1D *in_hist3[12][5][4][2];
  TH1D *in_hist4[12][5][4][2];

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

  double error[12][5][4][2];

  TCanvas *cfit[12][2];
  TCanvas *czoomed[12][2];
  TCanvas *cdiff_phi[12][5];
  TCanvas *cdiff_eta[12][5];
  TCanvas *craweta[12][5];
 
	
  float yleft[12][5][4][2];
  float yright[12][5][4][2];
	

  int bin;

  //----------------------------------------
  //   v2_jet[g][i][l] (for all centralities)
  //-------------------------------------------
  

  double bglevel, A_AS, V_2, V_1, V_3,alpha, beta, bglevelmax, bglevelmin, A_ASmax, A_ASmin, V_2max, V_1max, V_1min, V_3max, V_3min, V_2min, fitfunc0,fitmin0,fitmax0,histvalue, evalpt, error2, bglevel_err, V_1_err, V_2_err, V_3_err, A_AS_err, alpha_err, beta_err,bc, temp_bc, temp,temp_err, ymax, ymin, etamin1_val, etamax1_val, etamin2_val, etamax2_val, refvalue, referror, projmax,projmin;

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
      //    if(g>5&&l==0){ continue;}  

      fitcanvasname = "FitCanvas";
      fitcanvasname+= g;
      fitcanvasname+= l;
    
      cfit[g][l] = new TCanvas(fitcanvasname,"",0,0,1500,1600);
      cfit[g][l]->Divide(4,4,0.000001,0.000001);
    


      for(int i=0; i<4; i++){

	if(g>5&&i==4){ continue;}  //  No MC for pT>8

	//----------------------------------------------------
	//  Start of j loop -- individual histo analysis begins
	//-----------------------------------------------------
	
	for (int j=0; j<4; j++){

	  //	  if(j>0&&(g==7||g==9||g==11)){continue;}
      
	  cfit[g][l]->cd(4*i+j+1);
	  
	
	  in_name = make_name("Yield_",g,i,j,l,centlabel,pTlabel);


	  //******************************************//

	  //  SET BACKGROUND SAMPLING REGION HERE//

	  //******************************************//

	  projmax = 3.0;
	  projmin = 1.5;
		 
	  if(g<6){
	    etamin1_val = -projmax+0.001;
	    etamax1_val = -projmin-0.001;
	    etamin2_val =  projmin+0.001;
	    etamax2_val =  projmax -0.001;
	    
	  }


	  if(g>5){
	    if(g%2==0){
	      etamin1_val = -2.19999;
	      etamax1_val = -1.201;
	      etamin2_val =  1.201;
	      etamax2_val =  2.19999;
	    }else{
	      etamin1_val = -1.6999;
	      etamax1_val = -1.01;
	      etamin2_val =  1.01;
	      etamax2_val =  1.6999;
	    }
	  } 



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

	  if((g==2||g==4)&&i>1&&j==0){
	    TString Temp_name = in_name;
	    Temp_name.ReplaceAll("Yield_PbPb","Raw_Yield");
	  
	    yield[g][i][j][l] = (TH2D*)fdata[g]->Get(Temp_name)->Clone(in_name);
	  
	    Temp_name.ReplaceAll("Raw_Yield","Mixed_Event");
	    Temp_name.ReplaceAll("Cent50_Cent100","Cent30_Cent50");

	    TH2D *temp_me = (TH2D*)fdata[g]->Get(Temp_name)->Clone(Temp_name);
	    yield[g][i][j][l]->Divide(temp_me);
	  
	  }

	  if(g==2&&i==3&&j==3){
	    TString Temp_name = in_name;
	    Temp_name.ReplaceAll("Yield_PbPb","Raw_Yield");
	  
	    yield[g][i][j][l] = (TH2D*)fdata[g]->Get(Temp_name)->Clone(in_name);
	  
	    Temp_name.ReplaceAll("Raw_Yield","Mixed_Event");
	    Temp_name.ReplaceAll("Cent0_Cent10","Cent10_Cent30");

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


	  Int_t etamin1 = yield[g][i][j][l]->GetXaxis()->FindBin(etamin1_val);
	  Int_t etamax1 = yield[g][i][j][l]->GetXaxis()->FindBin(etamax1_val);
	  Int_t etamin2 = yield[g][i][j][l]->GetXaxis()->FindBin(etamin2_val);
	  Int_t etamax2 = yield[g][i][j][l]->GetXaxis()->FindBin(etamax2_val);
 
	  nbins = etamax1-etamin1+etamax2-etamin2+2;
	
	  in_hist[g][i][j][l] = yield[g][i][j][l]->ProjectionY(bkgname,etamin1,etamax1);
	  in_hist[g][i][j][l]->Add(yield[g][i][j][l]->ProjectionY(tempname,etamin2,etamax2));
	  in_hist[g][i][j][l]->Scale(1./nbins);

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

	}
      }
    }
  }
  
  cout<<"ready to start fitting"<<endl;


  int g = 4;
  //-------------------------------------
  //      PbPb Fits
  //-------------------------------------
  for(int l = 0; l<1; l++){
      
    for(int j = 0; j<4; j++){

      for(int i = 0; i<4; i++){

	
	in_name = make_name("Yield_",g,i,j,l,centlabel,pTlabel);


	for(int k = 1; k<51; k++){

	

	  in_hist[4][i][j][l]->SetBinContent(k+50,in_hist[2][i][j][l]->GetBinContent(k));
	  in_hist[4][i][j][l]->SetBinError(k+50,in_hist[2][i][j][l]->GetBinError(k));
	 
	}

      
	cfit[4][l]->cd(4*i+j+1);
		

	bglevel= in_hist[4][i][j][l]->GetBinContent(in_hist[4][i][j][l]->FindBin(TMath::Pi()/2));
	A_AS =in_hist[4][i][j][l]->GetMaximum();
	gen_gaus->SetParameter(0,bglevel);
	gen_gaus->SetParLimits(0,0.,1.0);
	gen_gaus->SetParLimits(2,0.,1.0);
	gen_gaus->SetParLimits(1,-1.,1.0);
	gen_gaus->SetParLimits(3,-1.,1.0);

	if(g==0||g==2||g==4){
	  v2_assoc[4][i][j][l] = find_v2(i,j);
	  V_2 = v2_assoc[4][i][j][l]*v2_jet[4][j][l];
	    
	  //	  if(V_2<0){cout<<"v2_jet = "<<v2_jet[4][j][l]<<" v2_assoc = "<<v2_assoc[4][i][j][l]<<" V2 = "<<V_2<<" WHAT?????"<<endl; return -1;}
	}


	V_3= 0.;
	V_1 = 0.;
	//	gen_gaus->FixParameter(1,0.);
	//	gen_gaus->FixParameter(3,0.);
	gen_gaus->FixParameter(4,0.);
	gen_gaus->FixParameter(5,0.);
	gen_gaus->FixParameter(6,0.);
	
      if(g==0||g==2||g==4){

	in_hist[4][i][j][l]->Fit("gen_gaus","","",-TMath::Pi()/2.,3.*TMath::Pi()/2.);

	bglevel = gen_gaus->GetParameter(0);
	V_1     = gen_gaus->GetParameter(1);
	V_2     = gen_gaus->GetParameter(2);
	V_3     = gen_gaus->GetParameter(3);

	bglevel_err = gen_gaus->GetParError(0);
	V_1_err = gen_gaus->GetParError(1);
	V_2_err = gen_gaus->GetParError(2);
	v2_jet_test[4][i][j][l] = V_2/v2_assoc[4][i][j][l];

	v1_tot[4][i][j][l] = V_1;
	v3_tot[4][i][j][l] = V_3;

  
      }

      //--------------------------------
      //    Fits for pp data/ PYTHIA
      //-------------------------------	  


      if(g==1||g==3||g==5||g==7||g==9||g==11){
	gen_gaus->SetParameter(0,bglevel);
	gen_gaus->ReleaseParameter(1);
	gen_gaus->SetParameter(1,0.01);
	gen_gaus->SetParLimits(1,0.,0.6);
	gen_gaus->FixParameter(2,0.);	
	gen_gaus->FixParameter(3,0.);
	gen_gaus->SetParameter(4,A_AS);
	gen_gaus->SetParameter(5,0.4);
	gen_gaus->SetParLimits(5,0.1,2.0);
	gen_gaus->ReleaseParameter(6);
	gen_gaus->SetParameter(6,1.1);
	gen_gaus->SetParLimits(6,1.,2.);


	in_hist[4][i][j][l]->Fit("gen_gaus","","",-1.,4.5);
	  
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
      
      //--------------------------------------------
      //   Draw fit plots
      //---------------------------------------------



      TString outhistproj_name = in_name;
      outhistproj_name.ReplaceAll("Yield","New1d_bkg");

      switch(i){
      case 0: ymax =yield_proj[4][i][j][l]->GetMaximum()+0.008;  break;
      case 1: ymax =yield_proj[4][i][j][l]->GetMaximum()+0.006;  break;
      case 2: ymax =yield_proj[4][i][j][l]->GetMaximum()+0.005;  break;
      case 3: ymax =yield_proj[4][i][j][l]->GetMaximum()+0.012;  break;

      }
	
      ymin =yield_proj[4][i][j][l]->GetMinimum()-0.002;
      in_hist_rebin[4][i][j][l] = Rebin_dPhi_full(in_hist[4][i][j][l]);
          
      in_hist_rebin[4][i][j][l]->Rebin(2);
      in_hist_rebin[4][i][j][l]->Scale(1/2.);


      if((in_hist_rebin[4][i][j][l]->GetMaximum()+.004)>ymax){ymax = in_hist_rebin[4][i][j][l]->GetMaximum()+.004;}
      //	in_hist_rebin[4][i][j][l]->SetAxisRange(-1.5,4.8);

      in_hist_rebin[4][i][j][l]->SetMaximum(ymax);
      in_hist_rebin[4][i][j][l]->SetMinimum(ymin);
      in_hist_rebin[4][i][j][l]->Draw();
	
	
      yield_rebin[4][i][j][l]= Rebin_dPhi_full(yield_proj[4][i][j][l]);
   

      TString yield_rebin_name = in_name;
      yield_rebin_name += "_1D_rebin";
      yield_rebin_name +=i;
      yield_rebin_name+=j;
      yield_rebin[4][i][j][l]->SetName(yield_rebin_name);
	


      yield_rebin[4][i][j][l]->SetMarkerColor(kMagenta+2);

      yield_rebin[4][i][j][l]->SetLineColor(kMagenta+2);
      yield_rebin[4][i][j][l]->SetMarkerStyle(10);

      yield_rebin[4][i][j][l]->SetMaximum(ymax);
      yield_rebin[4][i][j][l]->SetMinimum(ymin);

      yield_rebin[4][i][j][l]->Draw("same");
      in_hist_rebin[4][i][j][l]->Draw("same");
	
      TLegend *lfit = new TLegend(0.65,0.7,0.88,0.94);
      lfit->SetFillColor(kWhite);
      lfit->SetLineColor(kWhite);
      lfit->AddEntry(yield_rebin[4][i][j][l],"Yield");
      lfit->AddEntry(in_hist_rebin[4][i][j][l],"Sum. Bkg");
      lfit->SetTextSize(ts2-0.01);
      lfit->Draw("same");


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
	
      cout<<i<<" "<<j<<endl;

      } //closes the [i] (centrality) loop

    

    } //closes the [j] (pT) loop
    cout<<"about to save"<<endl;

  cfit[4][l]->SaveAs("DijetFreeFitsGluedBackgroundV1V2V3.pdf");
  cfit[4][l]->SaveAs("DijetFreeFitsGluedBackgroundV1V2V3.png");

  }

  for(int j = 0; j<4; j++){
    for(int i = 0; i<4; i++){

      if(i==0&&j==0)   cout<<"CentIndex TrkPtIndex JetV2 V1 V3"<<endl;
      cout<<j<<" "<<i<<" "<<v2_jet_test[4][i][j][0]<<" "<<v1_tot[4][i][j][0]<<" "<<v3_tot[4][i][j][0]<<endl;

      //  cout<<"Jet v2 err for"<<g<<" "<<j<<" "<<i<<" is:"<<v2_jet_test_err[4][i][j][l]<<endl;

      //	cout<<g<<" "<<i<<" "<<j<<" "<<yleft[4][i][j][0]<<" "<<yright[4][i][j][0]<<endl;
    }
  }
  
  return 0;
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





 


