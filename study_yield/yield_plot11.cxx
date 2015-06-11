#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TF2.h"
#include "TTree.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TExec.h"
#include "TLatex.h"

#include <iostream>
#include <vector>
#include <fstream>

#include "../HIN_14_016_functions.h"


using namespace std;

Int_t yield_plot11(int gstart= 0, int gend= 6, bool Closure_and_pp_subtracted = kTRUE, bool Is_NonLeading = kFALSE){
  // Default values run over data only.  gstart = 6, gend = 12 are values for full MC run

  if(Closure_and_pp_subtracted == kTRUE && Is_NonLeading == kTRUE){
    cerr<<"We have no closures for NonLeading data."<<endl;
      return -1;
  }


#include "../HIN_14_016_universals.h"

  //Set Style:

  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
    
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);

  
  

  TF1 *fit0 = new TF1("fit0","[0]",-3.,3.);
 
  float x;
  TF1 *do_offset = new TF1("do_offset","-1.*[0]+x-x",-3.,3.);
  float offset;
   
  TFile *fnewbg[14];
  TFile *fout[14];
  TFile *fclosures[12];

  TH2D *result[12][5][4][2];
  TH2D *result2[12][5][4][2];
  TH2D *fitbg[12][5][4][2];

  TH1D *check_old_phi[12][5][4][2];
  TH1D *check_old_phi_rebin[12][5][4][2];
  TH1D *check_old_eta[12][5][4][2];
  TH1D *check_old_eta_rebin[12][5][4][2];
 
  TH1D *check_new_eta[12][5][4][2];
  TH1D *check_new_eta_rebin[12][5][4][2];
  TH1D *check_new_eta_rebin2[12][5][4][2];
  TH1D *check_new_eta_syst[12][5][4][2];
  TH1D *check_new_eta_min[12][5][4][2];
  TH1D *check_new_eta_max[12][5][4][2];
 
  TH1D *check_new_phi[12][5][4][2];
  TH1D *check_new_phi_rebin[12][5][4][2];
  TH1D *check_new_phi_syst[12][5][4][2];
  TH1D *check_new_phi_max[12][5][4][2];
  TH1D *check_new_phi_min[12][5][4][2];

  
  TH1D *HYDJET_PYTHIA_eta_new[12][5][4][2];
  TH1D *HYDJET_PYTHIA_phi_new[12][5][4][2];
 
  TH1D *closure_phi[12][5][4][2];
  TH1D *closure_eta[12][5][4][2];


  
  TH1D *PbPb_pp_phi[12][5][4][2];
  TH1D *PbPb_pp_phi_syst[12][5][4][2];
  TH1D *PbPb_pp_eta[12][5][4][2];
  TH1D *PbPb_pp_eta_syst[12][5][4][2];

  TH1D *Error_BG_phi[12][5][4][2];
  TH1D *Error_Closure_phi[12][5][4][2];
  TH1D *Error_Relative_phi[12][5][4][2];
  TH1D *Error_BG_Closure_phi[12][5][4][2];
  TH1D *Error_BG_Closure_Relative_phi[12][5][4][2];
  

  TH1D *Error_BG_eta[12][5][4][2];
  TH1D *Error_Closure_eta[12][5][4][2];
  TH1D *Error_Relative_eta[12][5][4][2];
  TH1D *Error_BG_Closure_eta[12][5][4][2];
  TH1D *Error_BG_Closure_Relative_eta[12][5][4][2];

  TH1D *Integral_phi_Pt[12][4][2];
  TH1D *Integral_eta_Pt[12][4][2];

  TH1D *incl_dist_ref_eta[12][5][4][2];
  TH1D *incl_dist_ref_phi[12][5][4][2];

  TGraphAsymmErrors *Integral_eta_syst[12][4][2];
  TGraphAsymmErrors *Integral_phi_syst[12][4][2];

  TCanvas *ccheckEta_wide[12][2];
  TCanvas *ccheckEta[12][2];
  TCanvas *ccheckPhi[12][2];
  
  TCanvas *cHminusP_phi[12][2];
  TCanvas *cHminusP_eta[12][2];
 
  double Integral_eta[12][5][4][2];
  double Integral_eta_max[12][5][4][2];
  double Integral_eta_min[12][5][4][2];
  double Integral_eta_error_max[12][5][4][2];
  double Integral_eta_error_min[12][5][4][2];
  double Integral_eta_error[12][5][4][2];
 
  Double_t error_A[12][5][4][2];
  Double_t error_B[12][5][4][2];
  Double_t error_C[12][5][4][2];
  Double_t error_D[12][5][4][2];
  Double_t error_E[12][5][4][2];
  Double_t error_F[12][5][4][2];

  double syst_err[12][5][4][2];
  double Integral_phi[12][5][4][2];
  double Integral_phi_max[12][5][4][2];
  double Integral_phi_min[12][5][4][2];
  double Integral_phi_error_max[12][5][4][2];
  double Integral_phi_error_min[12][5][4][2];
  double Integral_phi_error[12][5][4][2];
 
  TString in_name, plotname, outname, funcname, centlabel, datalabel, pTlabel,checkcanvasnameEta,checkcanvasnamePhi,rawetacanvasname,HminusPphicanvasname, HminusPetacanvasname,PbPbminuscanvas_phi,PbPbminuscanvas_eta, checkcanvasnameEta_wide, EtaClosureName;

  double eta_inttest, phi_inttest;

  vector<double> pTbin_centers;
  pTbin_centers.push_back(1.5);
  pTbin_centers.push_back(2.5);
  pTbin_centers.push_back(3.5);
  pTbin_centers.push_back(6.0);
  vector<double> pTbin_errors;
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(2.);


  vector<double> Integral_eta_value;
  vector<double> Integral_eta_error_up;
  vector<double> Integral_eta_error_down;
 

  

  vector<double> Integral_phi_value;
  vector<double> Integral_phi_error_up;
  vector<double> Integral_phi_error_down;

  vector<double> Closure_integral_eta0;
  vector<double> Closure_integral_phi0;
  vector<double> Closure_integral_eta1;
  vector<double> Closure_integral_phi1;
  vector<double> Closure_integral_eta2;
  vector<double> Closure_integral_phi2;
  vector<double> Closure_integral_eta3;
  vector<double> Closure_integral_phi3;
 

  TGraph *Error_up_eta_pT[12][4][2];
    
  TGraph *Error_down_eta_pT[12][4][2];
   
  TGraph *Closure_integral_eta_pT[12][4][2];
  TGraph *Closure_integral_phi_pT[12][4][2];
   
 
 
  Double_t syst_error_tot,lphiedge1,lphiedge2,rphiedge1,rphiedge2,letaedge1,letaedge2,retaedge1,retaedge2,most_err,second_err,evalpt, syst_error_bg_c,re,ce,bc,PbPb_err, pp_err, tot_err, check_ymax, check_ymin, dx_eta, dx_phi, err_up, err_down,pp_error,PbPb_error;

  TF1 *gaus1d = new TF1("gaus1d","[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))");

  TF1 *gaus_phi[12][5][4];
  TF1 *gaus_eta[12][5][4];
 
  TLegend *lcheck, *leta, *lHminusP;
  
  TLine *linePhi, *lineEta;
  /*
  TFile *f2 = new TFile("pbpb_pp_data_inc120_dEta.root");
  TFile *f3 = new TFile("pbpb_pp_data_inc120_dPhi.root");
  */

  
  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
 
  for(int g=gstart; g<gend; g++){

    //  There will only be one big "g-loop".

    switch(g){
    case 0:
      if(Is_NonLeading == kFALSE){
	fnewbg[g] = new TFile("../bg_fits/PbPb_Inclusive_Yield_and_Bkg.root", "READ");

	fclosures[g] = new TFile("../spill_over/Inclusive_SpillOvers.root","READ");
	datalabel = "Inclusive";
      }
      
      if(Is_NonLeading == kTRUE){
	fnewbg[g] = new TFile("../bg_fits/PbPb_NonLeading_Yield_and_Bkgs.root", "READ");
	datalabel = "NonLeading";
	
      }
      break;

    case 1:
      if(Is_NonLeading == kFALSE){
	fnewbg[g] = new TFile("../bg_fits/pp_Inclusive_Yield_and_Bkg.root", "READ");
	if(Closure_and_pp_subtracted == kTRUE){	
	  fout[g] = new TFile("Inclusive_Data_AllPlots.root", "RECREATE");}
	if(Closure_and_pp_subtracted == kFALSE){	fout[g] = new TFile("Inclusive_Data_NoSpillOver_AllPlots.root", "RECREATE");}
	datalabel = "Inclusive";
      }	
    
      if(Is_NonLeading == kTRUE){
	fnewbg[g] = new TFile("../bg_fits/pp_Inclusive_Yield_and_Bkg.root", "READ");
	fout[g] = new TFile("NonLeading_Data_NoSpillOver_AllPlots.root", "RECREATE");
	datalabel = "NonLeading";
      }	
      break;

    case 2:
      fnewbg[g] = new TFile("../bg_fits/PbPb_SubLeading_Yield_and_Bkg.root", "READ");
      fclosures[g] = new TFile("../spill_over/SubLeading_SpillOvers.root","READ");
      datalabel = "SubLeading";
      break;

    case 3:
      fnewbg[g] = new TFile("../bg_fits/pp_SubLeading_Yield_and_Bkg.root", "READ");
      if(Closure_and_pp_subtracted ==kTRUE){    fout[g] = new TFile("SubLeading_Data_AllPlots.root", "RECREATE");}
      if(Closure_and_pp_subtracted ==kFALSE){    fout[g] = new TFile("SubLeading_Data_NoSpillOver_AllPlots.root", "RECREATE");}
      datalabel = "SubLeading";
      break;

    case 4:
      fnewbg[g] = new TFile("../bg_fits/PbPb_Leading_Yield_and_Bkg.root", "READ");
      fclosures[g] = new TFile("../spill_over/Leading_SpillOvers.root","READ");
      datalabel = "Leading";
      break;
     
    case 5:
      fnewbg[g] = new TFile("../bg_fits/pp_Leading_Yield_and_Bkg.root", "READ");
      if(Closure_and_pp_subtracted == kTRUE){  fout[g] = new TFile("Leading_Data_AllPlots.root", "RECREATE");}
      if(Closure_and_pp_subtracted == kFALSE){  fout[g] = new TFile("Leading_Data_NoSpillOver_AllPlots.root", "RECREATE");}
      datalabel = "Leading";
      break;
     
    case 6:
      fnewbg[g] = new TFile("../bg_fits/HYDJET_Inc120_2Dyield_and_NewBkg_files.root", "READ");
      datalabel = "Inclusive";
      break;

    case 7:
      fnewbg[g] = new TFile("../bg_fits/PYTHIA_Inc120_2Dyield_and_NewBkg_files.root", "READ");
      fout[g] = new TFile("Inclusive_Closures.root", "RECREATE");
      datalabel = "Inclusive";
      break;
    
    case 8:
      fnewbg[g] = new TFile("../bg_fits/HYDJET_Subleading50_2Dyield_and_NewBkg_files.root", "READ");
      datalabel = "SubLeading";
      break;
     
    case 9:
      fnewbg[g] = new TFile("../bg_fits/PYTHIA_Subleading50_2Dyield_and_NewBkg_files.root", "READ");
      fout[g] = new TFile("SubLeading_Closures.root", "RECREATE");
      datalabel = "SubLeading";
      break;
      
    case 10:
      fnewbg[g] = new TFile("../bg_fits/HYDJET_Leading120_2Dyield_and_NewBkg_files.root", "READ");
      datalabel = "Leading";
      break;

    case 11:
      fnewbg[g] = new TFile("../bg_fits/PYTHIA_Leading120_2Dyield_and_NewBkg_files.root", "READ");
      fout[g] = new TFile("Leading_Closures.root", "RECREATE");
      datalabel = "Leading";
      break;
    }



    //----------------------------------------------
    //     Start i and l loops -- set up canvases
    //------------------------------------------------

    for(int l = 0; l<1; l++){
    
      // cout << g << endl; continue;
      if(l>0&&g<6){ continue; } //no Gen/Reco for data

      //  if(l<1&&g>5){ continue; } //not interested in Gen
     
   
      checkcanvasnameEta = "CheckCanvasEta";
      checkcanvasnameEta+= g;
      checkcanvasnameEta+= l;
      ccheckEta[g][l] = new TCanvas(checkcanvasnameEta," ",10,10,1500,1600);
      ccheckEta[g][l]->Divide(4,4,0.,0.);
     
      checkcanvasnameEta_wide = "CheckCanvasEtaWide";
      checkcanvasnameEta_wide+= g;
      checkcanvasnameEta_wide+= l;
      ccheckEta_wide[g][l] = new TCanvas(checkcanvasnameEta_wide," ",10,10,1500,1600);
      ccheckEta_wide[g][l]->Divide(4,4,0.,0.);

      checkcanvasnamePhi = "CheckCanvasPhi";
      checkcanvasnamePhi+= g;
      checkcanvasnamePhi+= l;
      ccheckPhi[g][l] = new TCanvas(checkcanvasnamePhi," ",10,10,1500,1600);
      ccheckPhi[g][l]->Divide(4,4,0.,0.);



     
      // if(g>5&&l==1){
      if(g>5){

	HminusPetacanvasname = "HminusPcanvas_eta";
	HminusPetacanvasname+=g;
	HminusPetacanvasname+=l;
	cHminusP_eta[g][l] = new TCanvas(HminusPetacanvasname," ",10,10,1500,1600);
	cHminusP_eta[g][l]->Divide(4,4,0.,0.);

	HminusPphicanvasname = "HminusPcanvas_phi";
	HminusPphicanvasname+=g;
	HminusPphicanvasname+=l;
	cHminusP_phi[g][l] = new TCanvas(HminusPphicanvasname," ",10,10,1500,1600);
	cHminusP_phi[g][l]->Divide(4,4,0.,0.);

      }
      //----------------------------------------------------
      //  Start of main i & j loops 
      //-----------------------------------------------------
	
        
    

      for(int i=0; i<4; i++){
	
	for (int j=0; j<4; j++){

	  //  if((g==7||g==9||g==11)&&j>0){continue;}  //only one centrality class for PYTHIA
	  
	  in_name = make_name("Result_",g,i,j,l,centlabel,pTlabel);

	  if(g==6||g==8||g==10){in_name.ReplaceAll("Gen","PYTHIAHYDJET");}
	  if(g==7||g==9||g==11){in_name.ReplaceAll("Gen","PYTHIA");}

	  cout<<in_name<<endl;

	  result[g][i][j][l] = (TH2D*)fnewbg[g]->Get(in_name)->Clone(in_name);

	  TString in_name2 = in_name;
	  in_name2.ReplaceAll("Result","SummedResult");
	
	  result2[g][i][j][l] = (TH2D*)fnewbg[g]->Get(in_name2)->Clone(in_name2);

	  TString fitbg_name = in_name;
	  fitbg_name.ReplaceAll("Result","Fit_bkg");
	  fitbg[g][i][j][l] = (TH2D*)fnewbg[g]->Get(fitbg_name)->Clone(fitbg_name);
	
	 
	  if(g==6||g==8||g==10){in_name.ReplaceAll("PYTHIAHYDJET","Gen");}
	  if(g==7||g==9||g==11){in_name.ReplaceAll("PYTHIA","Gen");}



	  //-------------------------------
	  //dEta projection
	  //------------------------

	  TString summedchecknameEta = in_name;
	  summedchecknameEta.ReplaceAll("Result","Summed_check_Eta");
	    
	  llimiteta = result[g][i][j][l]->GetXaxis()->FindBin(-etalim+.001);
	  rlimiteta = result[g][i][j][l]->GetXaxis()->FindBin(etalim-.001);

	  llimitphi = result[g][i][j][l]->GetYaxis()->FindBin(-philim+.001);
	  rlimitphi = result[g][i][j][l]->GetYaxis()->FindBin(philim-.001);
	    

	  check_old_eta[g][i][j][l] = (TH1D*)result2[g][i][j][l]->ProjectionX(summedchecknameEta,llimitphi,rlimitphi);
	  dx_eta = check_old_eta[g][i][j][l]->GetBinWidth(1);
	  check_old_eta[g][i][j][l]->Scale(1/dx_eta);

	  TString summedchecknameEta_rebin = summedchecknameEta;
	  summedchecknameEta_rebin += "rebin";
	 
	  check_old_eta_rebin[g][i][j][l] = (TH1D*)Rebin_dEta(check_old_eta[g][i][j][l]);
	  check_old_eta_rebin[g][i][j][l]->SetName(summedchecknameEta_rebin);
	 

	  TString newchecknameEta = in_name;
	  newchecknameEta.ReplaceAll("Result", "New_check_Eta");

	  if(g==6||g==8||g==10){newchecknameEta.ReplaceAll("Gen","HYD"); }
	  if(g==7||g==9||g==11){newchecknameEta.ReplaceAll("Gen","PYTH"); }
	  check_new_eta[g][i][j][l] =result[g][i][j][l]->ProjectionX(newchecknameEta,llimitphi,rlimitphi);
	  check_new_eta[g][i][j][l]->Scale(1/dx_eta);
	 
	  TString newchecknameEta_rebin = newchecknameEta;
	  newchecknameEta_rebin += "rebin";

	  if(g<6){
	    if(Is_NonLeading == kFALSE){
	      check_new_eta_rebin[g][i][j][l] = (TH1D*)Rebin_dEta(check_new_eta[g][i][j][l]);
	      check_new_eta_rebin[g][i][j][l]->SetName(newchecknameEta_rebin);
	    }

	    if(Is_NonLeading == kTRUE){
	      check_new_eta_rebin[g][i][j][l] = (TH1D*) check_old_eta_rebin[g][i][j][l]->Clone(newchecknameEta_rebin);
	    }
	    
	  }

	  if(g>5){
 
	    check_new_eta_rebin[g][i][j][l] = (TH1D*)check_old_eta_rebin[g][i][j][l]->Clone(newchecknameEta_rebin);
	    
	  }


	  //-------------------------------
	  //   dPhi projection
	  //------------------------------ 
      
	  TString summedchecknamePhi = in_name;
	  summedchecknamePhi.ReplaceAll("Result","Summed_check_Phi");

	  //ccheckPhi[g][l]->cd(4*i+j+1);
     
	  check_old_phi[g][i][j][l] = result2[g][i][j][l]->ProjectionY(summedchecknamePhi,llimiteta,rlimiteta);
	  dx_phi = check_old_phi[g][i][j][l]->GetBinWidth(1);
	  check_old_phi[g][i][j][l]->Scale(1/dx_phi);

	 
	  TString summedchecknamePhi_rebin = summedchecknamePhi;
	  summedchecknamePhi_rebin += "rebin";

	  check_old_phi_rebin[g][i][j][l] =(TH1D*)Rebin_dPhi(check_old_phi[g][i][j][l]);
	  check_old_phi_rebin[g][i][j][l]->SetName(summedchecknamePhi_rebin);
	 
	  TString newchecknamePhi = in_name;
	  if(g==6||g==8||g==10){newchecknamePhi.ReplaceAll("Gen","HYD"); }
	  if(g==7||g==9||g==11){newchecknamePhi.ReplaceAll("Gen","PYTH"); }

	  newchecknamePhi.ReplaceAll("Result", "New_check_Phi");
	  check_new_phi[g][i][j][l] = result[g][i][j][l]->ProjectionY(newchecknamePhi,llimiteta,rlimiteta);
	  check_new_phi[g][i][j][l]->Scale(1/dx_phi);
	  

	  TString newchecknamePhi_rebin = newchecknamePhi;
	  newchecknamePhi_rebin += "rebin";


	  if(g<6){
	    check_new_phi_rebin[g][i][j][l] = (TH1D*)Rebin_dPhi(check_new_phi[g][i][j][l]);
	    check_new_phi_rebin[g][i][j][l]->SetName(newchecknamePhi_rebin);
	  }

	
	  //SCALING FOR pT 4-8.
	    
	  if(i==3){
	    check_old_eta_rebin[g][i][j][l]->Scale(1/4.);
	    check_new_eta_rebin[g][i][j][l]->Scale(1/4.);
	    check_old_phi_rebin[g][i][j][l]->Scale(1/4.);
	    check_new_phi_rebin[g][i][j][l]->Scale(1/4.);
	  }


	  //----------------------------------------------------
	  //   Calculate Systematic "BG Error" for ALL Check Plots
	  //-----------------------------------------------------



	  if(g<6){
	  
	    //error_A[g][i][j][l]: calculated as bin error at dPhi = 0, propogated over all 2D
	    
	    error_A[g][i][j][l] = fitbg[g][i][j][l]->GetBinError(26,20)*20/dx_phi;

	    //  Systematic source B:  calculated as the average of the 2 greatest out of the 4 most peripheral points in dPhi
	     
	    lphiedge1 = check_new_phi_rebin[g][i][j][l]->GetBinContent(1);
	    lphiedge2 = check_new_phi_rebin[g][i][j][l]->GetBinContent(2);
	    rphiedge1 = check_new_phi_rebin[g][i][j][l]->GetBinContent(16);
	    rphiedge2 = check_new_phi_rebin[g][i][j][l]->GetBinContent(17);

	   

	    if(lphiedge1>0)lphiedge1=0.;
	    if(lphiedge2>0)lphiedge2=0.;
	    if(rphiedge1>0)rphiedge1=0.;
	    if(rphiedge2>0)rphiedge2=0.;
	    most_err = TMath::Min(lphiedge1,lphiedge2);
	    if(rphiedge1<most_err){
	      second_err=most_err;
	      most_err=rphiedge1;
	    }
	    if(rphiedge2<most_err){
	      second_err = most_err;
	      most_err = rphiedge2;
	    }
	    if(rphiedge1>most_err&&rphiedge2>most_err){
	      second_err = TMath::Max(lphiedge1,lphiedge2);
	    }
	    error_B[g][i][j][l] = (most_err+second_err)/2.;


	    //   Systematic source C: calculated as the average of the 2 greatest out of the 4 most peripheral points in dEta

	    letaedge1 = check_new_eta_rebin[g][i][j][l]->GetBinContent(3);
	    letaedge2 = check_new_eta_rebin[g][i][j][l]->GetBinContent(4);
	    retaedge1 = check_new_eta_rebin[g][i][j][l]->GetBinContent(nbounds_eta-4);
	    retaedge2 = check_new_eta_rebin[g][i][j][l]->GetBinContent(nbounds_eta-3);
	    if(letaedge1>0)letaedge1=0.;
	    if(letaedge2>0)letaedge2=0.;
	    if(retaedge1>0)retaedge1=0.;
	    if(retaedge2>0)retaedge2=0.;
	    most_err = TMath::Min(letaedge1,letaedge2);
	    if(retaedge1<most_err){
	      second_err=most_err;
	      most_err=retaedge1;
	    }
	    if(retaedge2<most_err){
	      second_err = most_err;
	      most_err = retaedge2;
	    }
	    if(retaedge1>most_err&&retaedge2>most_err){
	      second_err = TMath::Max(letaedge1,letaedge2);
	    }
	    error_C[g][i][j][l] = (most_err+second_err)/2.;

	    
	    //  Systematic source D: calculated as the average of the absolute values of the deviation of the contents of the 1.5<|dEta|<2.0

	    error_D[g][i][j][l] = (TMath::Abs(check_new_eta_rebin[g][i][j][l]->GetBinContent(2))+TMath::Abs(check_new_eta_rebin[g][i][j][l]->GetBinContent(18)))/2.;
	 

	    error_E[g][i][j][l] = (TMath::Abs(check_new_eta_rebin[g][i][j][l]->GetBinContent(1)-check_new_eta_rebin[g][i][j][l]->GetBinContent(2))+TMath::Abs(check_new_eta_rebin[g][i][j][l]->GetBinContent(18)-check_new_eta_rebin[g][i][j][l]->GetBinContent(19)))/2.;
	      
	    //  With these corrections done, add all errors in quadrature:
	    //   syst_err[g][i][j][l]= TMath::Sqrt(error_B[g][i][j][l]*syst_checkB+error_C[g][i][j][l]*error_C[g][i][j][l]+error_D[g][i][j][l]*error_D[g][i][j][l]+error_E[g][i][j][l]*error_E[g][i][j][l]+syst_checkF*syst_checkF);

	    //D, E, & F only
	    syst_err[g][i][j][l]= TMath::Sqrt(error_D[g][i][j][l]*error_D[g][i][j][l]+error_F[g][i][j][l]*error_F[g][i][j][l]);

	      
	    // ***note: syst_err[g][i][j][l] accounts for all "BG error"
	   

	    if(syst_err[g][i][j][l]<0.01){syst_err[g][i][j][l]=0.01;}
	      
	
	    //----------------
	    // over-rides
	    //----------------
	  
	    if(g==0){
	      switch(i){
	      case 0:
		if(j==0){syst_err[g][i][j][l]=0.15;}
		if(j==1){syst_err[g][i][j][l]=0.16;}
		if(j==2){syst_err[g][i][j][l]=0.20;}
		if(j==3){syst_err[g][i][j][l]=0.22;}
		break;
	      case 1:
		if(j==0){syst_err[g][i][j][l]=0.09;}
		if(j==1){syst_err[g][i][j][l]=0.09;}
		if(j==2){syst_err[g][i][j][l]=0.11;}
		if(j==3){syst_err[g][i][j][l]=0.11;}
		break;
	      case 2:
		if(j==0){syst_err[g][i][j][l]=0.04;}
		if(j==1){syst_err[g][i][j][l]=0.05;}
		if(j==2){syst_err[g][i][j][l]=0.06;}
		if(j==3){syst_err[g][i][j][l]=0.07;}
		break;
	      case 3:
		syst_err[g][i][j][l] = 0.02;
		break;
	      }
	    }
	    if(g==2){
	      switch(i){
	      case 0:
		if(j==0){syst_err[g][i][j][l]=0.2;}
		if(j==1){syst_err[g][i][j][l]=0.23;}
		if(j==2){syst_err[g][i][j][l]=0.27;}
		if(j==3){syst_err[g][i][j][l]=0.3;}
		break;
	      case 1:
		if(j==0){syst_err[g][i][j][l]=0.11;}
		if(j==1){syst_err[g][i][j][l]=0.12;}
		if(j==2){syst_err[g][i][j][l]=0.12;}
		if(j==3){syst_err[g][i][j][l]=0.14;}
		break;
	      case 2:
		if(j==0){syst_err[g][i][j][l]=0.04;}
		if(j==1){syst_err[g][i][j][l]=0.05;}
		if(j==2){syst_err[g][i][j][l]=0.07;}
		if(j==3){syst_err[g][i][j][l]=0.09;}
		break;
	      case 3:
		syst_err[g][i][j][l] = 0.03;
		break;
	      }
	    }
	    if(g==4){
	      switch(i){
	      case 0:
		if(j==0){syst_err[g][i][j][l]=0.15;}
		if(j==1){syst_err[g][i][j][l]=0.16;}
		if(j==2){syst_err[g][i][j][l]=0.20;}
		if(j==3){syst_err[g][i][j][l]=0.22;}
		break;
	      case 1:
		if(j==0){syst_err[g][i][j][l]=0.09;}
		if(j==1){syst_err[g][i][j][l]=0.09;}
		if(j==2){syst_err[g][i][j][l]=0.11;}
		if(j==3){syst_err[g][i][j][l]=0.11;}
		break;
	      case 2:
		if(j==0){syst_err[g][i][j][l]=0.04;}
		if(j==1){syst_err[g][i][j][l]=0.05;}
		if(j==2){syst_err[g][i][j][l]=0.07;}
		if(j==3){syst_err[g][i][j][l]=0.09;}
		break;
	      case 3:
		syst_err[g][i][j][l] = 0.02;
		break;
	      }
	    }	
	 

	    //=------------------------------------------------------------------
	    //  Save individual histos for  different sources of error
	    //  ***Section above defines BG error.  Closure error is calculated from fits in the loop below.  Relative error is Sqrt(3*05*.05)*bincontent
	    //
	    //   Also create systematic error histograms for use in plotting
	    //-------------------------------------------------------------------
	  
	    if(Closure_and_pp_subtracted == kTRUE&&(g==0||g==2||g==4)){
	      EtaClosureName = in_name;
	      EtaClosureName.ReplaceAll("Result_","SpillOvers_Eta_");
	      if(g==2){	  EtaClosureName.ReplaceAll("Eta_","Eta_SubLeading_");        }
	      if(g==4){	  EtaClosureName.ReplaceAll("Eta_","Eta_Leading_");        }
	      EtaClosureName.ReplaceAll("PbPb_","RecoJet_GenTrack_");  //actually recojet gentrack
	      EtaClosureName.ReplaceAll("Pt100_Pt300_","");

	      //   cout<<EtaClosureName<<endl;

	      closure_eta[g][i][j][l] = (TH1D*)fclosures[g]->Get(EtaClosureName)->Clone(EtaClosureName);
	      
	      TString PhiClosureName = EtaClosureName;
	      PhiClosureName.ReplaceAll("Eta","Phi");

	      closure_phi[g][i][j][l] = (TH1D*) fclosures[g]->Get(PhiClosureName)->Clone(PhiClosureName);
	    
	    }
	    
	    TString Error_BG_eta_name = in_name;
	    Error_BG_eta_name.ReplaceAll("Result","Error_BG_eta");
	    Error_BG_eta_name.ReplaceAll("pp","PbPb");

	    TString Error_Closure_eta_name = in_name;
	    Error_Closure_eta_name.ReplaceAll("Result","Error_Closure_eta");
	    Error_Closure_eta_name.ReplaceAll("pp","PbPb");
	      	   
	    TString Error_Relative_eta_name = in_name;
	    Error_Relative_eta_name.ReplaceAll("Result","Error_Relative_eta");
	    Error_Relative_eta_name.ReplaceAll("pp","PbPb");

	    TString Error_BG_Closure_eta_name = in_name;
	    Error_BG_Closure_eta_name.ReplaceAll("Result","Error_BG_Closure_eta");
	    Error_BG_Closure_eta_name.ReplaceAll("pp","PbPb");

	    TString Error_BG_Closure_Relative_eta_name = in_name;
	    Error_BG_Closure_Relative_eta_name.ReplaceAll("Result","Error_BG_Closure_Relative_eta");
	    Error_BG_Closure_Relative_eta_name.ReplaceAll("pp","PbPb");
 
	

	    Error_BG_eta[g][i][j][l] = new TH1D (Error_BG_eta_name,"",nbounds_eta-1,bin_bounds_eta);
	    Error_Closure_eta[g][i][j][l]= new TH1D (Error_Closure_eta_name,"",nbounds_eta-1,bin_bounds_eta);
	    Error_Relative_eta[g][i][j][l]= new TH1D (Error_Relative_eta_name,"",nbounds_eta-1,bin_bounds_eta);
	    Error_BG_Closure_eta[g][i][j][l]= new TH1D (Error_BG_Closure_eta_name,"",nbounds_eta-1,bin_bounds_eta);
	    Error_BG_Closure_Relative_eta[g][i][j][l]= new TH1D (Error_BG_Closure_Relative_eta_name,"",nbounds_eta-1,bin_bounds_eta);
  
	    
	    TString newchecknameEtasyst = newchecknameEta;
	    newchecknameEtasyst.ReplaceAll("New_check","Syst_error");
	    
	    check_new_eta_syst[g][i][j][l] = new TH1D(newchecknameEtasyst,"",nbounds_eta-1,bin_bounds_eta);

	    if(Closure_and_pp_subtracted==kFALSE||(g==1||g==3||g==5)){ce = 0.;}

	    for(int m=1; m<nbounds_eta;m++){

	      
	      re = TMath::Sqrt(3*0.05*0.05)*check_new_eta_rebin[g][i][j][l]->GetBinContent(m);
	      if((g==0||g==2||g==4)&&Closure_and_pp_subtracted == kTRUE){ce = closure_eta[g][i][j][l]->GetBinContent(m)/2; }
	      syst_error_bg_c =TMath::Sqrt( syst_err[g][i][j][l]*syst_err[g][i][j][l]+ce*ce);
		
	      syst_error_tot =TMath::Sqrt(syst_err[g][i][j][l]*syst_err[g][i][j][l]+ce*ce+re*re);
	
	      //	syst_error_tot = syst_err[g][i][j][l];
	
	
	      bc = check_new_eta_rebin[g][i][j][l]->GetBinContent(m);
	      check_new_eta_syst[g][i][j][l]->SetBinContent(m,bc);
	      check_new_eta_syst[g][i][j][l]->SetBinError(m,syst_error_tot);
	      //check_new_eta_syst[g][i][j][l]->SetBinError(m,ce);
	
	      Error_BG_eta[g][i][j][l]->SetBinContent(m,syst_err[g][i][j][l]);
	      Error_Closure_eta[g][i][j][l]->SetBinContent(m,ce);
	      Error_Relative_eta[g][i][j][l]->SetBinContent(m,re);
	      Error_BG_Closure_eta[g][i][j][l]->SetBinContent(m,syst_error_bg_c);
	      Error_BG_Closure_Relative_eta[g][i][j][l]->SetBinContent(m,syst_error_tot);
  	    }
	  	    	  
	    TString Error_BG_phi_name = in_name;
	    Error_BG_phi_name.ReplaceAll("Result","Error_BG_phi");
	    Error_BG_phi_name.ReplaceAll("pp","PbPb");

	    TString Error_Closure_phi_name = in_name;
	    Error_Closure_phi_name.ReplaceAll("Result","Error_Closure_phi");
	    Error_Closure_phi_name.ReplaceAll("pp","PbPb");

	    TString Error_Relative_phi_name = in_name;
	    Error_Relative_phi_name.ReplaceAll("Result","Error_Relative_phi");
	    Error_Relative_phi_name.ReplaceAll("pp","PbPb");

	    TString Error_BG_Closure_phi_name = in_name;
	    Error_BG_Closure_phi_name.ReplaceAll("Result","Error_BG_Closure_phi");
	    Error_BG_Closure_phi_name.ReplaceAll("pp","PbPb");

	    TString Error_BG_Closure_Relative_phi_name = in_name;
	    Error_BG_Closure_Relative_phi_name.ReplaceAll("Result","Error_BG_Closure_Relative_phi");
	    Error_BG_Closure_Relative_phi_name.ReplaceAll("pp","PbPb");
 

	    Error_BG_phi[g][i][j][l] = new TH1D (Error_BG_phi_name,"",nbounds_phi-1,bin_bounds_phi);
	    Error_Closure_phi[g][i][j][l]= new TH1D (Error_Closure_phi_name,"",nbounds_phi-1,bin_bounds_phi);
	    Error_Relative_phi[g][i][j][l]= new TH1D (Error_Relative_phi_name,"",nbounds_phi-1,bin_bounds_phi);
	    Error_BG_Closure_phi[g][i][j][l]= new TH1D (Error_BG_Closure_phi_name,"",nbounds_phi-1,bin_bounds_phi);
	    Error_BG_Closure_Relative_phi[g][i][j][l]= new TH1D (Error_BG_Closure_Relative_phi_name,"",nbounds_phi-1,bin_bounds_phi);
  


	    TString newchecknamePhisyst = newchecknamePhi;
	    newchecknamePhisyst.ReplaceAll("New_check","Syst_error");
	     
	    check_new_phi_syst[g][i][j][l] = new TH1D(newchecknamePhisyst,"",nbounds_phi-1,bin_bounds_phi);

	      
	    //   Fill histograms

	   

	    for(int m=1; m<nbounds_phi;m++){

	      re = TMath::Sqrt(3*0.05*.05)*check_new_phi_rebin[g][i][j][l]->GetBinContent(m);

	      if((g==0||g==2||g==4)&&Closure_and_pp_subtracted==kTRUE){ ce = closure_phi[g][i][j][l]->GetBinContent(m)/2.;}

	      syst_error_bg_c =TMath::Sqrt(syst_err[g][i][j][l]*syst_err[g][i][j][l]+ce*ce);
		
	      syst_error_tot =TMath::Sqrt(syst_err[g][i][j][l]*syst_err[g][i][j][l]+ce*ce+re*re);
		
	      //	syst_error_tot = syst_err[g][i][j][l];
		
	      bc = check_new_phi_rebin[g][i][j][l]->GetBinContent(m);

	      check_new_phi_syst[g][i][j][l]->SetBinContent(m,bc);
	      check_new_phi_syst[g][i][j][l]->SetBinError(m,syst_error_tot);
	      //	cout<<g<<" "<<i<<" "<<j<<" "<<m<<" "<<syst_error_tot<<endl;
		
	      Error_BG_phi[g][i][j][l]->SetBinContent(m,syst_err[g][i][j][l]);
	      Error_Closure_phi[g][i][j][l]->SetBinContent(m,ce);
	      Error_Relative_phi[g][i][j][l]->SetBinContent(m,re);
	      Error_BG_Closure_phi[g][i][j][l]->SetBinContent(m,syst_error_bg_c);
	      Error_BG_Closure_Relative_phi[g][i][j][l]->SetBinContent(m,syst_error_tot);
	    }
   
	    float val_l, val_r, err_l, err_r;

	    float nbins2 = check_new_eta_rebin[g][i][j][l]->GetNbinsX()+1;
	    /*
	    for(int k = 1; k<nbins2/2+1;k++){
	    
	      val_l = check_new_eta_rebin[g][i][j][l]->GetBinContent(k);
	      err_l  = check_new_eta_rebin[g][i][j][l]->GetBinError(k);
	    
	      val_r = check_new_eta_rebin[g][i][j][l]->GetBinContent(nbins2-k);
	      err_r = check_new_eta_rebin[g][i][j][l]->GetBinError(nbins2-k);
	    
	      cout<<k<<" "<<nbins2-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<TMath::Sqrt(err_l*err_l+err_r*err_r)/2<<endl;

	      check_new_eta_rebin[g][i][j][l]->SetBinContent(k,(val_l+val_r)/2);
	      check_new_eta_rebin[g][i][j][l]->SetBinContent(nbins2-k,(val_l+val_r)/2);
	   
	      if(g==0||g==2||g==4){
		check_new_eta_syst[g][i][j][l]->SetBinContent(k,(val_l+val_r)/2);
		check_new_eta_syst[g][i][j][l]->SetBinContent(nbins2-k,(val_l+val_r)/2);
	      }
	      if(k<nbins2/2){
		check_new_eta_rebin[g][i][j][l]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_eta_rebin[g][i][j][l]->SetBinError(nbins2-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      }
	    }

*/

	  }  // All of that was a GIANT Data-Only loop.



	  //------------------------------------------------------------------
	  //      YIELD PLOTS:  First, settings in both dEta and dPhi together
	  //-------------------------------------------------------------------
	 
	  ccheckEta_wide[g][l]->cd(4*i+j+1);

	  switch(i){
	  case 0: 
	    check_ymax = 11.;
	    check_ymin = -1.; 
	    break;
	  case 1: 
	    check_ymax = 6.5; 
	    check_ymin = -1.; 
	    break;
	  case 2: 
	    check_ymax = 4.2; 
	    check_ymin = -.45; 
	    break;
	  case 3: 
	    check_ymax = 2.2; 
	    check_ymin = -0.25; 
	   
	    break;
	  default: 
	    check_ymax = 10.; break;
	    check_ymin = -1.; break;
	  }
	  
	  

	  switch(g){
	  case 0: 
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(10);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(10);
	    check_new_eta_syst[g][i][j][l]->SetFillColor(90);
	    check_new_phi_syst[g][i][j][l]->SetFillColor(90);
	    break;
	  case 1:
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(4);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(4);
	    break;
	  case 2: 
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(34);
	    check_new_eta_rebin[g][i][j][l]->SetMarkerSize(2);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(34);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerSize(2);
	    check_new_eta_syst[g][i][j][l]->SetFillColor(30);
	    check_new_phi_syst[g][i][j][l]->SetFillColor(30);
	    break;
	  case 3:
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(28);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(28);
	    break;  
	  case 4: 
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(21);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(21);
	    check_new_eta_syst[g][i][j][l]->SetFillColor(kOrange-2);
	    check_new_phi_syst[g][i][j][l]->SetFillColor(kOrange-2);
	    break;
	  case 5:
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(28);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(28);
	    break;
	  case 6: 
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(10);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(10);
	    break;
	  case 7:
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(4);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(4);
	    break;
	  case 8: 
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(34);
	    check_new_eta_rebin[g][i][j][l]->SetMarkerSize(2);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(34);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerSize(2);
	    break;
	  case 9:
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(28);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(28);
	    break;  
	  case 10: 
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(21);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(21);
	    break;
	  case 11:
	    check_new_eta_rebin[g][i][j][l]->SetMarkerStyle(25);
	    check_new_phi_rebin[g][i][j][l]->SetMarkerStyle(25);
	    break;
	  default:
	    break;
	  }
	 
	  check_new_eta_rebin[g][i][j][l]->SetMarkerColor(kBlack);
	  check_new_eta_rebin[g][i][j][l]->SetLineColor(kBlack);
	  
	  check_new_phi_rebin[g][i][j][l]->SetMarkerColor(kBlack);
	  check_new_phi_rebin[g][i][j][l]->SetLineColor(kBlack);
	
	
	  if((g==0||g==2||g==4)&&Closure_and_pp_subtracted==kTRUE&&i<3){
	        
	    check_new_eta_rebin[g][i][j][l]->Add(closure_eta[g][i][j][l],-1.);
	    check_old_eta_rebin[g][i][j][l]->Add(closure_eta[g][i][j][l],-1.);
	    check_new_eta_syst[g][i][j][l]->Add(closure_eta[g][i][j][l],-1.);
   
	    check_new_phi_rebin[g][i][j][l]->Add(closure_phi[g][i][j][l],-1.);
	    check_old_phi_rebin[g][i][j][l]->Add(closure_phi[g][i][j][l],-1.);
	    check_new_phi_syst[g][i][j][l]->Add(closure_phi[g][i][j][l],-1.);
	     
	  }

	  //------------------
	  // dEta plotting
	  //-------------------

	  check_new_eta_rebin[g][i][j][l]->SetMaximum(check_ymax);
	  check_new_eta_rebin[g][i][j][l]->SetMinimum(check_ymin);

   
	  check_new_eta_rebin[g][i][j][l]->SetMarkerColor(1);
	  check_new_eta_rebin[g][i][j][l]->SetAxisRange(-2.5+.2,2.5-.2,"x");

	  check_new_eta_rebin[g][i][j][l]->Draw();	 
	  
      
	  check_new_eta_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(tstitle);
	  check_new_eta_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	  check_new_eta_rebin[g][i][j][l]->GetXaxis()->SetTitle("#Delta#eta");
	  check_new_eta_rebin[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	  check_new_eta_rebin[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	  check_new_eta_rebin[g][i][j][l]->GetYaxis()->SetTitle("Y#equiv1/N_{jet} d^{2}N/(d#Delta#etadp_{T}) (GeV/c)^{-1}");
	  check_new_eta_rebin[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	  check_new_eta_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);

 
	  check_new_eta_rebin[g][i][j][l]->GetXaxis()->CenterTitle();
	  check_new_eta_rebin[g][i][j][l]->GetYaxis()->CenterTitle();
	  if(i<3){
	    check_new_eta_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	  }
	  if(j>0){
	    check_new_eta_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	    check_new_eta_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	  }

	  if(g==0||g==2||g==4){
	    check_new_eta_syst[g][i][j][l]->SetMarkerSize(0);
	    check_new_eta_syst[g][i][j][l]->Draw("same e2");
	   
	  }
	 
	  
	    if(g==0||g==2||g==4){
	    check_old_eta_rebin[g][i][j][l]->SetMarkerStyle(10);
	    check_old_eta_rebin[g][i][j][l]->SetMarkerSize(1);
	    check_old_eta_rebin[g][i][j][l]->SetMarkerColor(kRed);
	    //   check_old_eta_rebin[g][i][j][l]->Draw("same");
	    }
	  
	    
	  if(g==1||g==3||g==5){
	    check_new_eta_syst[g-1][i][j][l]->Draw("same e2");
	    check_new_eta_rebin[g-1][i][j][l]->Draw("same");
	  }

	  check_new_eta_rebin[g][i][j][l]->Draw("same");


	  drawlabels_4by4(g,i,j,datalabel);
	      
	  lineEta = new TLine(-2.5,0,2.5,0);
	  lineEta->SetLineStyle(2);
	  lineEta->Draw("same");


	  //---------------
	  // Narrower axis range for presentation
	  //----------------

	  
	  ccheckEta[g][l]->cd(4*i+j+1);

	  TString narrowetaname = "NarrowEta";
	  narrowetaname+=g;
	  narrowetaname+=i;
	  narrowetaname+=j;
	  narrowetaname+=l;

	  check_new_eta_rebin2[g][i][j][l] = (TH1D*)check_new_eta_rebin[g][i][j][l]->Clone(narrowetaname);
	      
	  check_new_eta_rebin2[g][i][j][l]->SetAxisRange(-1.5+.1,1.5-.1,"x");

	  check_new_eta_rebin2[g][i][j][l]->Draw();	 
	  
	  if(g==0||g==2||g==4){
	    check_new_eta_syst[g][i][j][l]->Draw("same e2");
	  }

	  if(g==1||g==3||g==5){
	    check_new_eta_syst[g-1][i][j][l]->Draw("same e2");
	    check_new_eta_rebin[g-1][i][j][l]->Draw("same");
	  }
	  check_new_eta_rebin[g][i][j][l]->Draw("same");

	  drawlabels_4by4(g,i,j,datalabel);
	 
	  lineEta = new TLine(-1.5,0,1.5,0);
	  lineEta->SetLineStyle(2);
	  lineEta->Draw("same");

	  //----------------------------
	  //   Now dPhi
	  //-------------------------------


	  ccheckPhi[g][l]->cd(4*i+j+1);

	
	  check_new_phi_rebin[g][i][j][l]->SetMaximum(check_ymax);
	  check_new_phi_rebin[g][i][j][l]->SetMinimum(check_ymin);
	    
	  check_new_phi_rebin[g][i][j][l]->Draw("same");   
      
 
	  check_new_phi_rebin[g][i][j][l]->GetXaxis()->CenterTitle();
	  check_new_phi_rebin[g][i][j][l]->GetYaxis()->CenterTitle();
	  check_new_phi_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(tstitle);
	  check_new_phi_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	  check_new_phi_rebin[g][i][j][l]->GetXaxis()->SetTitle("#Delta#phi");
	  check_new_phi_rebin[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	  check_new_phi_rebin[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	  check_new_phi_rebin[g][i][j][l]->GetYaxis()->SetTitle("Y#equiv1/N_{jet} d^{2}N/(d#Delta#phidp_{T}) (GeV/c)^{-1}");
	  check_new_phi_rebin[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	  check_new_phi_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);
	  if(i<3){
	    check_new_phi_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	  }
	  if(j>0){
	
	    check_new_phi_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	    check_new_phi_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	  }
	    
	  if(g==0||g==2||g==4){
	    check_new_phi_syst[g][i][j][l]->SetMarkerSize(0);
	    check_new_phi_syst[g][i][j][l]->Draw("same e2");
	  }

	  if(g==1||g==3||g==5){
	    check_new_phi_syst[g-1][i][j][l]->Draw("same e2");
	    check_new_phi_rebin[g-1][i][j][l]->Draw("same");
	  }

	  
	  check_new_phi_rebin[g][i][j][l]->Draw("same");


	  /*

	    check_old_phi_rebin[g][i][j][l]->SetMarkerStyle(10);
	    check_old_phi_rebin[g][i][j][l]->SetMarkerSize(1);
	    check_old_phi_rebin[g][i][j][l]->SetMarkerColor(kRed);
	    check_old_phi_rebin[g][i][j][l]->Draw("same");
	
	  */

	  /*
	    if(g==0){
	    incl_dist_ref_phi[g][i][j][l]->SetMarkerColor(kRed);
	    incl_dist_ref_phi[g][i][j][l]->SetLineColor(kRed);
	    incl_dist_ref_phi[g][i][j][l]->Draw("same");
	    }
	  */


   
	  drawlabels_4by4(g,i,j,datalabel);

	  linePhi = new TLine(-1.51,0,1.51,0);
	  linePhi->SetLineStyle(2);
	  linePhi->Draw("same");

	 
	
	  //---------------------------------------------------------------------------
	  // WRITE ALL PLOTS (associated with correlations & correlations syst. error)
	  //---------------------------------------------------------------------------

	

	  if(g==1||g==3||g==5){
	   	     
	    Error_BG_phi[g-1][i][j][l]->Write();
	    Error_Closure_phi[g-1][i][j][l]->Write();
	    Error_Relative_phi[g-1][i][j][l]->Write();
	    Error_BG_Closure_phi[g-1][i][j][l]->Write();
	    Error_BG_Closure_Relative_phi[g-1][i][j][l]->Write();
	      	      
	    Error_BG_eta[g-1][i][j][l]->Write();
	    Error_Closure_eta[g-1][i][j][l]->Write();
	    Error_Relative_eta[g-1][i][j][l]->Write();
	    Error_BG_Closure_eta[g-1][i][j][l]->Write();
	    Error_BG_Closure_Relative_eta[g-1][i][j][l]->Write();
	    
	    check_new_eta_rebin[g-1][i][j][l]->Write();
	    check_new_phi_rebin[g-1][i][j][l]->Write();
	    check_old_eta_rebin[g-1][i][j][l]->Write();
	    check_old_phi_rebin[g-1][i][j][l]->Write();
	    	    
 
	    check_new_eta_syst[g-1][i][j][l]->Write();
	    check_new_phi_syst[g-1][i][j][l]->Write();
	   

	    check_new_eta_rebin[g][i][j][l]->Write();
	    check_new_phi_rebin[g][i][j][l]->Write();
	    check_old_eta_rebin[g][i][j][l]->Write();
	    check_old_phi_rebin[g][i][j][l]->Write();
      	  
	  }

	  //----------------------------------------------------
	  //              Integral Yield Calculations
	  //   *** Done in the pp loop so that pp plots exist***
	  //     Must be drawn after the [i] (pT) loop is closed
	  //----------------------------------------------------
	     	  

	  

	  if(g==1||g==3||g==5){

	    
	    TString PbPbppname_eta =in_name;
	    PbPbppname_eta.ReplaceAll("Result_pp","PbPb_minus_pp_eta");

	    TString PbPbppname_eta_syst =in_name;
	    PbPbppname_eta_syst.ReplaceAll("Result_pp","PbPb_minus_pp_eta_syst");

	    TString PbPbppname_eta_max =in_name;
	    PbPbppname_eta_max.ReplaceAll("Result_pp","PbPb_check_eta_max");

	    TString PbPbppname_eta_min =in_name;
	    PbPbppname_eta_min.ReplaceAll("Result_pp","PbPb_check_eta_min");
	      
	      

	    PbPb_pp_eta[g][i][j][l] = (TH1D*)check_new_eta_rebin[g-1][i][j][l]->Clone(PbPbppname_eta);
	    PbPb_pp_eta[g][i][j][l]->Add(check_new_eta_rebin[g][i][j][l],-1.);
	      
	    PbPb_pp_eta_syst[g][i][j][l] = (TH1D*)PbPb_pp_eta[g][i][j][l]->Clone(PbPbppname_eta_syst);
	      

	    
	    for(int k=1;k<check_new_eta_rebin[g][i][j][l]->GetNbinsX()+1;k++){
		
	      PbPb_err = check_new_eta_syst[g-1][i][j][l]->GetBinError(k);
	      pp_err = check_new_eta_syst[g][i][j][l]->GetBinError(k);

	      if(isnan(pp_err) || isinf(pp_err) ){pp_err=0.;}
	
	      tot_err = TMath::Sqrt(PbPb_err*PbPb_err+pp_err*pp_err);

	      PbPb_pp_eta_syst[g][i][j][l]->SetBinError(k,tot_err);
		
	    }
	      
	      
	    check_new_eta_max[g][i][j][l] = (TH1D*)check_new_eta_rebin[g-1][i][j][l]->Clone(PbPbppname_eta_max);
	    check_new_eta_max[g][i][j][l]->Add(Error_BG_Closure_eta[g-1][i][j][l]);

	    check_new_eta_min[g][i][j][l] = (TH1D*)check_new_eta_rebin[g-1][i][j][l]->Clone(PbPbppname_eta_min);
	    check_new_eta_min[g][i][j][l]->Add(Error_BG_Closure_eta[g-1][i][j][l],-1.);
	

	    if(i==0){
	      TString IntegratedYieldname_eta = in_name;
	      IntegratedYieldname_eta.ReplaceAll("Result_pp","Integrated_Yield_Eta");
	      IntegratedYieldname_eta.ReplaceAll("Pt100_Pt300_","");
	      IntegratedYieldname_eta.ReplaceAll("_TrkPt1_TrkPt2","");
	      Integral_eta_Pt[g][j][l] = new TH1D(IntegratedYieldname_eta,"",4,pTbins);
	    }
	     	      

	    llimiteta = PbPb_pp_eta[g][i][j][l]->FindBin(-etalim+.001);
	    rlimiteta = PbPb_pp_eta[g][i][j][l]->FindBin(etalim-.001);

	     
	  
	    
	    Integral_eta[g][i][j][l] = check_new_eta_rebin[g-1][i][j][l]->IntegralAndError(llimiteta,rlimiteta,PbPb_error,"width");
	      
	     

	    Integral_eta_error[g][i][j][l] = TMath::Sqrt(PbPb_error*PbPb_error+pp_error*pp_error);  

	    for(int k=llimiteta;k<rlimiteta+1;k++){
	      if((check_new_eta_min[g][i][j][l]->GetBinContent(k)<0)&&(check_new_eta_rebin[g-1][i][j][l]->GetBinContent(k)>=0)){
		check_new_eta_min[g][i][j][l]->SetBinContent(k,0.);
	      }
	    }
	     	     	      
	    Integral_eta_max[g][i][j][l] = check_new_eta_max[g][i][j][l]->Integral(llimiteta,rlimiteta,"width")-Integral_eta[g][i][j][l];
	    
	    Integral_eta_min[g][i][j][l] = check_new_eta_min[g][i][j][l]->Integral(llimiteta,rlimiteta,"width")-Integral_eta[g][i][j][l];

	
	    Integral_eta_error_max[g][i][j][l] = TMath::Sqrt(Integral_eta_max[g][i][j][l]*Integral_eta_max[g][i][j][l]+Integral_eta[g][i][j][l]*Integral_eta[g][i][j][l]*3*0.05*0.05);
	    Integral_eta_error_min[g][i][j][l] = TMath::Sqrt(Integral_eta_min[g][i][j][l]*Integral_eta_min[g][i][j][l]+Integral_eta[g][i][j][l]*Integral_eta[g][i][j][l]*3*0.05*0.05);
	    
	    if(Closure_and_pp_subtracted==kTRUE){Integral_eta[g][i][j][l]+= -1.*check_new_eta_rebin[g][i][j][l]->IntegralAndError(llimiteta,rlimiteta,pp_error,"width");	  }    
 
	       //if(Closure_and_pp_subtracted==kFALSE){Integral_eta[g][i][j][l]+= -1.*check_new_eta_rebin[g][i][j][l]->IntegralAndError(llimiteta,rlimiteta,pp_error,"width");	  }    


	    Integral_eta_Pt[g][j][l]->SetBinContent(i+1,Integral_eta[g][i][j][l]);
	    Integral_eta_Pt[g][j][l]->SetBinError(i+1,Integral_eta_error[g][i][j][l]);
	   

	      
	 
	    //------------------------
	    //    Phi integrals	  
	    //------------------------



	    TString PbPbppname_phi =in_name;
	    PbPbppname_phi.ReplaceAll("Result_pp","PbPb_minus_pp_phi");

	    TString PbPbppname_phi_max =in_name;
	    PbPbppname_phi_max.ReplaceAll("Result_pp","PbPb_minus_phi_max");

	    TString PbPbppname_phi_min =in_name;
	    PbPbppname_phi_min.ReplaceAll("Result_pp","PbPb_minus_phi_min");
	      
	    TString PbPbppname_phi_syst =in_name;
	    PbPbppname_phi_syst.ReplaceAll("Result_pp","PbPb_minus_pp_phi_syst");

	  
	    PbPb_pp_phi[g][i][j][l] = (TH1D*)check_new_phi_rebin[g-1][i][j][l]->Clone(PbPbppname_phi);
	    PbPb_pp_phi[g][i][j][l]->Add(check_new_phi_rebin[g][i][j][l],-1.);
	  
	 
	    
	    PbPb_pp_phi_syst[g][i][j][l] = (TH1D*)PbPb_pp_phi[g][i][j][l]->Clone(PbPbppname_phi_syst);
	      

	      
	    check_new_phi_max[g][i][j][l] = (TH1D*)check_new_phi_rebin[g-1][i][j][l]->Clone(PbPbppname_phi_max);
	    check_new_phi_max[g][i][j][l]->Add(Error_BG_Closure_phi[g-1][i][j][l]);

	    check_new_phi_min[g][i][j][l] = (TH1D*)check_new_phi_rebin[g-1][i][j][l]->Clone(PbPbppname_phi_max);
	    check_new_phi_min[g][i][j][l]->Add(Error_BG_Closure_phi[g-1][i][j][l],-1.);

	      
	    for(int k=1;k<check_new_phi_rebin[g][i][j][l]->GetNbinsX()+1;k++){
		
	      PbPb_err = check_new_phi_syst[g-1][i][j][l]->GetBinError(k);
	      pp_err = check_new_phi_syst[g][i][j][l]->GetBinError(k);
		
	      if(isnan(pp_err) || isinf(pp_err) ){pp_err=0.;}
	      if(isnan(PbPb_err) || isinf(PbPb_err) ){PbPb_err=0.;}

	      tot_err = TMath::Sqrt(PbPb_err*PbPb_err+pp_err*pp_err);

	      if(isnan(tot_err)){
		tot_err = 0.;

		cerr<<"I am so confused by not-a-number :( "<<endl;
		//return -1;
	      }
		
	      PbPb_pp_phi_syst[g][i][j][l]->SetBinError(k,tot_err);

	      //		cout<<g<<" "<<i<<" "<<j<<" "<<l<<k<<"   "<<PbPb_err<<" "<<pp_err<<" "<<tot_err<<endl;

	    }
	 



	    if(i==0){
	      TString IntegratedYieldname_phi = in_name;
	      IntegratedYieldname_phi.ReplaceAll("Result_pp","Integrated_Yield_Phi");
	      IntegratedYieldname_phi.ReplaceAll("Pt100_Pt300_","");
	      IntegratedYieldname_phi.ReplaceAll("_TrkPt1_TrkPt2","");
	      Integral_phi_Pt[g][j][l] = new TH1D(IntegratedYieldname_phi,"",4,pTbins);
	    }
	 
	    llimitphi = PbPb_pp_phi[g][i][j][l]->FindBin(-philim+.0001);
	    rlimitphi = PbPb_pp_phi[g][i][j][l]->FindBin(philim-.0001);
	 

	    Integral_phi[g][i][j][l] = check_new_phi_rebin[g-1][i][j][l]->IntegralAndError(llimitphi,rlimitphi,PbPb_error,"width");
	      
	    Integral_phi_error[g][i][j][l] = TMath::Sqrt(PbPb_error*PbPb_error+pp_error*pp_error);  
	      
	    for(int k=llimitphi;k<rlimitphi+1;k++){
	      if((check_new_phi_min[g][i][j][l]->GetBinContent(k)<0)&&(check_new_phi_rebin[g-1][i][j][l]->GetBinContent(k)>=0)){
		check_new_phi_min[g][i][j][l]->SetBinContent(k,0.);
	      }
	    }
	     	     	      
	    Integral_phi_max[g][i][j][l] = check_new_phi_max[g][i][j][l]->Integral(llimitphi,rlimitphi,"width")-Integral_phi[g][i][j][l];
	    
	    Integral_phi_min[g][i][j][l] = check_new_phi_min[g][i][j][l]->Integral(llimitphi,rlimitphi,"width")-Integral_phi[g][i][j][l];

	
	    Integral_phi_error_max[g][i][j][l] = TMath::Sqrt(Integral_phi_max[g][i][j][l]*Integral_phi_max[g][i][j][l]+Integral_phi[g][i][j][l]*Integral_phi[g][i][j][l]*3*0.05*0.05);
	    Integral_phi_error_min[g][i][j][l] = TMath::Sqrt(Integral_phi_min[g][i][j][l]*Integral_phi_min[g][i][j][l]+Integral_phi[g][i][j][l]*Integral_phi[g][i][j][l]*3*0.05*0.05);

	     if(Closure_and_pp_subtracted == kTRUE){Integral_phi[g][i][j][l]+= -1.*check_new_phi_rebin[g][i][j][l]->IntegralAndError(llimitphi,rlimitphi,pp_error,"width");	  }    


	     //if(Closure_and_pp_subtracted == kFALSE){Integral_phi[g][i][j][l]+= -1.*check_new_phi_rebin[g][i][j][l]->IntegralAndError(llimitphi,rlimitphi,pp_error,"width");	  }    



	    Integral_phi_Pt[g][j][l]->SetBinContent(i+1,Integral_phi[g][i][j][l]);
	    Integral_phi_Pt[g][j][l]->SetBinError(i+1,Integral_phi_error[g][i][j][l]);
	  	   
	    PbPb_pp_eta[g][i][j][l]->Write();
	    PbPb_pp_eta_syst[g][i][j][l]->Write();
	    PbPb_pp_phi[g][i][j][l]->Write();
	    PbPb_pp_phi_syst[g][i][j][l]->Write();
	    
	    
	  } //close Integral loop (for g==1, 3, 5)
	  
	  
	 
	} //close j  because otherwise we ignore all centralities after j=0!

	  //-------------------------------------------------------
	  //        CLOSURE CORRECTION PLOTS & CALCULATION
	  //           Draw and fit HYDJET minus PYTHIA
	  //--------------------------------------------------
	
      
	for(int j = 0; j<4; j++){

	  // if((g==7||g==9||g==11)&&l==1){
	  if((g==7||g==9||g==11)){

	    cout<<"                   "<<endl;
	    cout<<"------------------------------------------------------------------------------------"<<endl;
	    cout<<datalabel<<"    "<<g<<" "<<i<<" "<<j<<" "<<in_name<<endl;
	    cout<<"------------------------------------------------------------------------------------"<<endl;
	     

	    //-------------------
	    // All names at once: 
	    //---------------------

	    in_name = make_name("Result_",g,i,j,l,centlabel,pTlabel);

	    TString HminusPnameEta = in_name;
	    HminusPnameEta.ReplaceAll("Pt100_Pt300_","");
	    HminusPnameEta.ReplaceAll("Result","Raw_HYD_minus_PYTH_Eta");
	     
	    TString gaus_eta_name = in_name;
	    gaus_eta_name.ReplaceAll("Pt100_Pt300_","");
	    gaus_eta_name.ReplaceAll("Result","GausFit_Eta");


	    EtaClosureName = in_name;
	    EtaClosureName.ReplaceAll("Pt100_Pt300_","");
	    EtaClosureName.ReplaceAll("Result","Closures_Eta");
	    
	    TString HminusPnamePhi = HminusPnameEta;
	    HminusPnamePhi.ReplaceAll("Eta","Phi");
	   	     
	    TString gaus_phi_name = gaus_eta_name;
	    gaus_phi_name.ReplaceAll("Eta","Phi");


	    TString PhiClosureName = EtaClosureName;
	    PhiClosureName.ReplaceAll("Eta","Phi");
	   
	    //HYDJET minus PYTHIA eta

	    cHminusP_eta[g][l]->cd(4*i+j+1);

	   

	    HYDJET_PYTHIA_eta_new[g][i][j][l] = (TH1D*)check_new_eta_rebin[g-1][i][j][l]->Clone(HminusPnameEta);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->SetName(HminusPnameEta);
	    
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->Add(check_new_eta_rebin[g][i][0][l],-1.);
	      
	  
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->SetLineColor(1);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->SetMarkerStyle(10);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->SetMarkerColor(1);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->SetMarkerSize(1);


	    HYDJET_PYTHIA_eta_new[g][i][j][l]->SetMaximum(check_ymax);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->SetMinimum(check_ymin);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->Draw();     


	    
 
      
 
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->GetXaxis()->CenterTitle();
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->GetYaxis()->CenterTitle();
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->GetYaxis()->SetLabelSize(tstitle);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->GetXaxis()->SetTitle("#Delta#eta");
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->GetYaxis()->SetTitle("1/N_{jet} d^{2}N/(d#Delta#phi dp_{T})");
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);
	    if(i<3){
	      HYDJET_PYTHIA_eta_new[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	    }
	    if(j>0){
	      HYDJET_PYTHIA_eta_new[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	      HYDJET_PYTHIA_eta_new[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	    }
	      
	    llimiteta = HYDJET_PYTHIA_eta_new[g][i][j][l]->GetXaxis()->FindBin(-1.0+.0001);
	    rlimiteta = HYDJET_PYTHIA_eta_new[g][i][j][l]->GetXaxis()->FindBin(1.0-.0001);
	     
	    //   if(g==11||(g==7&&(!((i==0&&j>0)||(i==1&&j>1))))){
	 
	    /*
	    if(g==11||(g==7)){
		
	      for(int k = llimiteta; k<rlimiteta+1; k++){
		if(HYDJET_PYTHIA_eta_new[g][i][j][l]->GetBinContent(k)<0){
		  HYDJET_PYTHIA_eta_new[g][i][j][l]->SetBinContent(k,0.);
		}
	      }
		
	    }
	    */

	    //  double Yield_eta = check_new_eta_rebin[g][i][j][l]->Integral(llimiteta,rlimiteta,"width")-check_new_eta_rebin[g+1][i][0]][l]->Integral(llimiteta,rlimiteta,"width");

	    double Yield_eta = HYDJET_PYTHIA_eta_new[g][i][j][l]->Integral(llimiteta,rlimiteta,"width");	      
	      
	     	     
	      
	    if(Yield_eta<0.){Yield_eta=0.;}

	      
	    gaus1d->FixParameter(0,0.);
	    gaus1d->FixParameter(1,Yield_eta);
	      
	    gaus1d->ReleaseParameter(2);
	    gaus1d->SetParameter(2,0.2);
	    gaus1d->SetParLimits(2,0.1,0.4);
	    if(i>1){    gaus1d->SetParLimits(2,0.05,0.2);}
	   	    	      
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->Fit("gaus1d","","",-1.,1.);
	      
	  
	    gaus_eta[g][i][j] = new TF1(gaus_eta_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

	    for(int k= 0; k<3; k++){
	
	      double temp = gaus1d->GetParameter(k);
	      gaus_eta[g][i][j]->SetParameter(k, temp);

	      //	cout<<g<<"  "<<i<<"  "<<j<<"  "<<l<<"  "<<"  "<<k<<"  "<<temp<<endl;	
	    }

	    double err_temp= gaus1d->GetParError(2);

	    //cout<<"INTEGRAL CHECK "<<Yield_eta<<" "<<gaus_eta[g][i][j]->Integral(-1.,1.)<<endl;

	    HYDJET_PYTHIA_eta_new[g][i][j][l]->SetMaximum(check_ymax);
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->SetMinimum(check_ymin);
	      
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->Draw();     

	 

	    closure_eta[g][i][j][l] = (TH1D*)HYDJET_PYTHIA_eta_new[g][i][j][l]->Clone(EtaClosureName);
	      

	    for(int k=0; k<20; k++){
	      evalpt = HYDJET_PYTHIA_eta_new[g][i][j][l]->GetBinCenter(k);
	      bc = gaus1d->Eval(evalpt);
	      closure_eta[g][i][j][l]->SetBinContent(k,bc);
	      // closure_eta[g][i][j][l]->SetBinError(k,err_temp);
	      closure_eta[g][i][j][l]->SetBinError(k,0.);
	    }

	    gaus_eta[g][i][j]->SetLineColor(kBlue);
	    gaus_eta[g][i][j]->Draw("same");

	    closure_eta[g][i][j][l]->SetMarkerColor(kBlue);
	    closure_eta[g][i][j][l]->SetLineColor(kBlue);
	    closure_eta[g][i][j][l]->SetMarkerStyle(10);
	    closure_eta[g][i][j][l]->Draw("same");
	    

	    drawlabels_4by4(g,i,j,datalabel);
	      
	    lineEta = new TLine(-2.5,0,2.5,0);
	    lineEta->SetLineStyle(2);
	    lineEta->Draw("same");


	  
	    if(j>0){ lHminusP= new TLegend(legendoffset,texty3,texty1,0.95);
	      lHminusP->SetTextSize(ts);
	    }
	    if(j==0){ 
	      lHminusP = new TLegend(legendoffset2-.04,texty3,texty1,0.95);
	      lHminusP->SetTextSize(ts2);
	    }

	    lHminusP->SetFillColor(kWhite);
	    lHminusP->SetLineColor(kWhite);

	    lHminusP->AddEntry(HYDJET_PYTHIA_eta_new[g][i][j][l],"HYD.-PYTH.","lpfe");
	     	
	    if(i==3){lHminusP->SetTextSize(ts2);}
	    lHminusP->Draw("same");



	    //HYDJET minus PYTHIA phi


	    cHminusP_phi[g][l]->cd(4*i+j+1);

	    HYDJET_PYTHIA_phi_new[g][i][j][l] = (TH1D*)check_new_phi_rebin[g-1][i][j][l]->Clone(HminusPnamePhi);


	    HYDJET_PYTHIA_phi_new[g][i][j][l]->SetName(HminusPnamePhi);

	    HYDJET_PYTHIA_phi_new[g][i][j][l]->Add(check_new_phi_rebin[g][i][0][l],-1.);
	  
	 
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->SetLineColor(1);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->SetMarkerStyle(10);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->SetMarkerColor(1);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->SetMarkerSize(1);


	    HYDJET_PYTHIA_phi_new[g][i][j][l]->SetMaximum(check_ymax);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->SetMinimum(check_ymin);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->Draw();     
      
 
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->GetXaxis()->CenterTitle();
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->GetYaxis()->CenterTitle();
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->GetYaxis()->SetLabelSize(tstitle);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->GetXaxis()->SetTitle("#Delta#phi");
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->GetYaxis()->SetTitle("1/N_{jet} d^{2}N/(d#Delta#phi dp_{T})");
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);
	    if(i<3){
	      HYDJET_PYTHIA_phi_new[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	    }
	    if(j>0){
	      HYDJET_PYTHIA_phi_new[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	      HYDJET_PYTHIA_phi_new[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	    }

	    llimitphi = HYDJET_PYTHIA_phi_new[g][i][j][l]->GetXaxis()->FindBin(-1.0+.01);
	    rlimitphi = HYDJET_PYTHIA_phi_new[g][i][j][l]->GetXaxis()->FindBin(1.0-.001);
	   	      
	    //if(g==11||(g==7&&(!(i==0&&j>0)))){
	  
	    /*
	    if(g==11||(g==7)){
		
	      for(int k = llimitphi; k<rlimitphi+1; k++){
		if(HYDJET_PYTHIA_phi_new[g][i][j][l]->GetBinContent(k)<0){
		  HYDJET_PYTHIA_phi_new[g][i][j][l]->SetBinContent(k,0.);
		}
	      }
		
	    }
	    */

	    // double Yield_phi = check_new_phi_rebin[g][i][j][l]->Integral(llimitphi,rlimitphi,"width")- check_new_phi_rebin[g+1][i][0][l]->Integral(llimitphi,rlimitphi,"width");

	    double Yield_phi = HYDJET_PYTHIA_phi_new[g][i][j][l]->Integral(llimitphi,rlimitphi,"width");	      
	      
	     	

	    //  if(Yield_phi2>Yield_phi){Yield_phi=Yield_phi2;}
	      
	    if(Yield_phi<0){Yield_phi=0;}
	    

	    gaus1d->SetParLimits(2,0.1,0.5);

	    if(g==6&&i==0&&j==3){
	      double width_temp = gaus_phi[g][0][2]->GetParameter(2);
	      gaus1d->FixParameter(2,width_temp);
	    }
	  	     	      
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->Fit("gaus1d");
	    
	    gaus_phi[g][i][j] = new TF1(gaus_phi_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

	    for(int k= 0; k<3; k++){
	      double temp = gaus1d->GetParameter(k);
	      gaus_phi[g][i][j]->SetParameter(k, temp);
	    }

	    err_temp = gaus1d->GetParError(2);
	 
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->SetMaximum(check_ymax);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->SetMinimum(check_ymin);
	    HYDJET_PYTHIA_phi_new[g][i][j][l]->Draw();    

	    gaus_phi[g][i][j]->SetLineColor(kBlue);
	    gaus_phi[g][i][j]->Draw("same"); 


	    closure_phi[g][i][j][l] = (TH1D*)HYDJET_PYTHIA_phi_new[g][i][j][l]->Clone(PhiClosureName);
      
	    for(int k=0; k<20; k++){
	      evalpt = HYDJET_PYTHIA_phi_new[g][i][j][l]->GetBinCenter(k);
	      bc = gaus1d->Eval(evalpt);
	      closure_phi[g][i][j][l]->SetBinContent(k,bc);
	      //  closure_phi[g][i][j][l]->SetBinError(k,err_temp);
	      closure_phi[g][i][j][l]->SetBinError(k,0.);
	    }

	    gaus_phi[g][i][j]->SetLineColor(kBlue);
	    gaus_phi[g][i][j]->Draw("same");

	    closure_phi[g][i][j][l]->SetMarkerColor(kBlue);
	    closure_phi[g][i][j][l]->SetLineColor(kBlue);
	    closure_phi[g][i][j][l]->SetMarkerStyle(10);
	    closure_phi[g][i][j][l]->Draw("same");
	    

	    drawlabels_4by4(g,i,j,datalabel);

	    linePhi->Draw("same");
	    lHminusP->Draw("same");

	 
	      
	    if(i==0&&j==0){
	      Closure_integral_eta0.clear();
	      Closure_integral_phi0.clear();
	      Closure_integral_eta1.clear();
	      Closure_integral_phi1.clear();
	      Closure_integral_eta2.clear();
	      Closure_integral_phi2.clear();
	      Closure_integral_eta3.clear();
	      Closure_integral_phi3.clear();
	    }

	      
	    switch(j){
	    case 0:
	      Closure_integral_eta0.push_back(Yield_eta);
	      Closure_integral_phi0.push_back(Yield_eta);
	      break;
	    case 1:
	      Closure_integral_eta1.push_back(Yield_eta);
	      Closure_integral_phi1.push_back(Yield_eta);
	      break;
	    case 2:
	      Closure_integral_eta2.push_back(Yield_eta);
	      Closure_integral_phi2.push_back(Yield_eta);
	      break;
	    case 3:
	      Closure_integral_eta3.push_back(Yield_eta);
	      Closure_integral_phi3.push_back(Yield_eta);
	      break;
	    }
	    
	    

	    check_new_eta_rebin[g-1][i][j][l]->Write();
	    check_new_phi_rebin[g-1][i][j][l]->Write();
	   
	    check_new_eta_rebin[g][i][j][l]->Write();
	    check_new_phi_rebin[g][i][j][l]->Write();
	    
	    
	    HYDJET_PYTHIA_eta_new[g][i][j][l]->Write();
	    closure_eta[g][i][j][l]->Write();
	    gaus_eta[g][i][j]->Write();

	    HYDJET_PYTHIA_phi_new[g][i][j][l]->Write();
	    closure_phi[g][i][j][l]->Write();
	    gaus_phi[g][i][j]->Write();
	  


	  } //close closure correction calcualtion loop	    

	} //close j
      } // close i

  

      //  if((g==7||g==9||g==11)&&l==1){  
      if((g==7||g==9||g==11)){  

	for(int j = 0; j<4; j++){   //we needed to be in a "j loop" not in an "i loop"
	  
	  switch(j){
	  case 0:
	    Closure_integral_eta_pT[g][j][l] = new TGraph(pTbin_centers.size(), &pTbin_centers[0],&Closure_integral_eta0[0]);
	    Closure_integral_phi_pT[g][j][l] = new TGraph(pTbin_centers.size(), &pTbin_centers[0],&Closure_integral_phi0[0]);
	    break; 
	  case 1:
	    Closure_integral_eta_pT[g][j][l] = new TGraph(pTbin_centers.size(), &pTbin_centers[0],&Closure_integral_eta1[0]);
	    Closure_integral_phi_pT[g][j][l] = new TGraph(pTbin_centers.size(), &pTbin_centers[0],&Closure_integral_phi1[0]);
	    break; 
	  case 2:
	    Closure_integral_eta_pT[g][j][l] = new TGraph(pTbin_centers.size(), &pTbin_centers[0],&Closure_integral_eta2[0]);
	    Closure_integral_phi_pT[g][j][l] = new TGraph(pTbin_centers.size(), &pTbin_centers[0],&Closure_integral_phi2[0]);
	    break; 
	  case 3:
	    Closure_integral_eta_pT[g][j][l] = new TGraph(pTbin_centers.size(), &pTbin_centers[0],&Closure_integral_eta3[0]);
	    Closure_integral_phi_pT[g][j][l] = new TGraph(pTbin_centers.size(), &pTbin_centers[0],&Closure_integral_phi3[0]);
	    break;
	  }


	  TString ClosureIntegralEtaPt_name = in_name;
	  ClosureIntegralEtaPt_name.ReplaceAll("Result","Closure_Integral_Eta");
	  ClosureIntegralEtaPt_name.ReplaceAll("TrkPt4_TrkPt8","");
	
	 
	  switch(j){
	  case 0:
	    ClosureIntegralEtaPt_name.ReplaceAll("Cent0_Cent10_Pt100_Pt300_","Cent50_Cent100");
	  break;
	  case 1:
	    ClosureIntegralEtaPt_name.ReplaceAll("Cent0_Cent10_Pt100_Pt300_","Cent30_Cent50");
	  break;
	  case 2:
	    ClosureIntegralEtaPt_name.ReplaceAll("Cent0_Cent10_Pt100_Pt300_","Cent10_Cent30");
	  break;
	  case 3:
	    ClosureIntegralEtaPt_name.ReplaceAll("Cent0_Cent10_Pt100_Pt300_","Cent0_Cent10");
	  break;
	  }

	  Closure_integral_eta_pT[g][j][l]->SetName(ClosureIntegralEtaPt_name);
	  Closure_integral_eta_pT[g][j][l]->Write();
	  

	  TString ClosureIntegralPhiPt_name = ClosureIntegralEtaPt_name;
	  ClosureIntegralPhiPt_name.ReplaceAll("Eta","Phi");
	  Closure_integral_phi_pT[g][j][l]->SetName(ClosureIntegralPhiPt_name);
	  Closure_integral_phi_pT[g][j][l]->Write();

	}
      }
    

      //-------------------------------------------

      //  Save all canvases

      //--------------------------------------------

	
      TString checksavenamePhi = in_name;
      checksavenamePhi.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenamePhi.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenamePhi.ReplaceAll("Result","Yield_Phi");
      checksavenamePhi +=".png";
      ccheckPhi[g][l]->SaveAs(checksavenamePhi);
   
	

      TString checksavenameEta = in_name;
      checksavenameEta.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenameEta.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenameEta.ReplaceAll("Result","Yield_Eta");
      checksavenameEta +=".png";
      ccheckEta[g][l]->SaveAs(checksavenameEta);


      TString checksavenameEta_wide = in_name;
      checksavenameEta_wide.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenameEta_wide.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      if(g==6||g==8||g==10){checksavenameEta_wide.ReplaceAll("Result","Yield_Eta_wide_HYDJET");}
      if(g==7||g==9||g==11){checksavenameEta_wide.ReplaceAll("Result","Yield_Eta_wide_PYTHIA");}
      checksavenameEta_wide.ReplaceAll("Result","Yield_Eta_wide");
      checksavenameEta_wide +=".png";
      ccheckEta_wide[g][l]->SaveAs(checksavenameEta_wide);
      checksavenameEta_wide.ReplaceAll(".png",".pdf");
      ccheckEta_wide[g][l]->SaveAs(checksavenameEta_wide);



      // if((g==7||g==9||g==11)&&l==1){
      if((g==7||g==9||g==11)){

	TString HminusPsavenameEta = in_name;
	HminusPsavenameEta.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
	HminusPsavenameEta.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
	HminusPsavenameEta.ReplaceAll("Result","HYDJET_minus_PYTHIA_Eta");
	HminusPsavenameEta +=".pdf";
	cHminusP_eta[g][l]->SaveAs(HminusPsavenameEta);

	TString HminusPsavenamePhi = in_name;
	HminusPsavenamePhi.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
	HminusPsavenamePhi.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
	HminusPsavenamePhi.ReplaceAll("Result","HYDJET_minus_PYTHIA_Phi");
	HminusPsavenamePhi +=".pdf";
	cHminusP_phi[g][l]->SaveAs(HminusPsavenamePhi);

      }

      //---------------------------------------------------
      //     Make Integral Histograms (but don't draw them)
      //--------------------------------------------------
	


      if((g==1||g==3||g==5)&&l==0){


	for(int j=0;j<4;j++){  //Start a new j-loop not in an i-loop -- we're plotting integrated yield by pT for each centrality class

	  //   Here make the asymmetric systematic errors:
	  
	  in_name =  make_name("Result_",g,3,j,l,centlabel,pTlabel);
	  

	  Integral_eta_value.clear();
	  Integral_eta_error_up.clear();
	  Integral_eta_error_down.clear();


	  for(int k = 0; k<4; k++){
	    bc = Integral_eta_Pt[g][j][l]->GetBinContent(k+1);
	    Integral_eta_value.push_back(bc);

	    err_up = Integral_eta_error_max[g][k][j][l];
	    Integral_eta_error_up.push_back(err_up);
					   
	    err_down = Integral_eta_error_min[g][k][j][l];
	    Integral_eta_error_down.push_back(err_down);


	     
	  }

	  // return 0; 
	  TString IntegratedYieldname_eta_syst = in_name;
	  IntegratedYieldname_eta_syst.ReplaceAll("Result_pp","Integrated_Yield_Eta_Syst");
	  IntegratedYieldname_eta_syst.ReplaceAll("Pt100_Pt300_","");
	  IntegratedYieldname_eta_syst.ReplaceAll("_TrkPt4_TrkPt8","");

	  Integral_eta_syst[g][j][l] = new TGraphAsymmErrors(pTbin_centers.size(),&pTbin_centers[0],&Integral_eta_value[0],&pTbin_errors[0],&pTbin_errors[0],&Integral_eta_error_down[0],&Integral_eta_error_up[0]);
	    
	  Integral_eta_syst[g][j][l]->SetName(IntegratedYieldname_eta_syst);

	 
	  //-------------
	  //  Now Phi
	  //--------------
	

	   	 

	  Integral_phi_value.clear();
	  Integral_phi_error_up.clear();
	  Integral_phi_error_down.clear();


	  for(int k = 0; k<4; k++){
	    bc = Integral_phi_Pt[g][j][l]->GetBinContent(k+1);
	    Integral_phi_value.push_back(bc);

	    err_up = Integral_phi_error_max[g][k][j][l];
	    Integral_phi_error_up.push_back(err_up);
					   
	    err_down = Integral_phi_error_min[g][k][j][l];
	    Integral_phi_error_down.push_back(err_down);

	    //	      cout<<Integral_phi_error_max[g][k][j][l]<<" "<<Integral_phi_error_up.at(k)<<endl;

	  }

	 
	  TString IntegratedYieldname_phi_syst = in_name;
	  IntegratedYieldname_phi_syst.ReplaceAll("Result_pp","Integrated_Yield_Phi_Syst");
	  IntegratedYieldname_phi_syst.ReplaceAll("Pt100_Pt300_","");
	  IntegratedYieldname_phi_syst.ReplaceAll("_TrkPt4_TrkPt8","");
		  	    
	  Integral_phi_syst[g][j][l] = new TGraphAsymmErrors(pTbin_centers.size(),&pTbin_centers[0],&Integral_phi_value[0],&pTbin_errors[0],&pTbin_errors[0],&Integral_phi_error_down[0],&Integral_phi_error_up[0]);

	  Integral_phi_syst[g][j][l]->SetName(IntegratedYieldname_phi_syst);
	  
	 

	  Integral_eta_Pt[g][j][l]->Write();
	  Integral_phi_Pt[g][j][l]->Write();
	  Integral_eta_syst[g][j][l]->Write();
	  Integral_phi_syst[g][j][l]->Write();

	   
	  //Auxillary TGraphs
   
	  Error_up_eta_pT[g][j][l] = new TGraph(pTbin_centers.size(),&pTbin_centers[0],&Integral_eta_error_up[0]);
	  
	  TString ErrorUpEtaPt_name = in_name;
	  ErrorUpEtaPt_name.ReplaceAll("Result_pp","Integral_Error_Up");
	  ErrorUpEtaPt_name.ReplaceAll("Pt100_Pt300_","");
	  ErrorUpEtaPt_name.ReplaceAll("_TrkPt4_TrkPt8","");
	 
	  Error_up_eta_pT[g][j][l]->SetName(ErrorUpEtaPt_name);
	  Error_up_eta_pT[g][j][l]->Write();
	    

	  Error_down_eta_pT[g][j][l] = new TGraph(pTbin_centers.size(),&pTbin_centers[0],&Integral_eta_error_down[0]);
	  TString ErrorDownEtaPt_name = ErrorUpEtaPt_name;
	  ErrorDownEtaPt_name.ReplaceAll("Up","Down");
	  Error_down_eta_pT[g][j][l]->SetName(ErrorDownEtaPt_name);
	  Error_down_eta_pT[g][j][l]->Write();
	    	    
	    
	}  //  end j loop for integral plotting	


      }  //end only do integrals for pp
       

    }  // Closes the [l] (gen vs. reco) loop
    
    /*

    //to cout errors
    // for(int i = 0; i<4; i++){
      for(int j = 0; j<4; j++){
	//	cout<<g<<" "<<i<<" "<<j<<" "<<error_A[g][i][j][0]<<" "<<error_B[g][i][j][0]<<" "<<error_C[g][i][j][0]<<" "<<error_D[g][i][j][0]<<" "<<error_E[g][i][j][0]<<" "<<error_F[g][i][j][0]<<" "<<syst_err[g][i][j][0]<<endl;
	cout<<g<<" "<<i<<" "<<j<<" "<< Closure_integral_eta_pT[g][j][0]->Eval(pTbin_centers.at(i))<<endl;
	

      }
      //  }
 


      */

  }// closes g loop
 
  return 0;


}  //and we're done.
 
