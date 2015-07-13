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

Int_t study_yield(int gstart= 0, int gend= 6, bool Spillover_and_pp_subtracted = kTRUE, bool Is_NonLeading = kFALSE){
  // Default values run over data only.  gstart = 6, gend = 12 are values for full MC run

  if(Spillover_and_pp_subtracted == kTRUE && Is_NonLeading == kTRUE){
    cerr<<"We have no spillovers for NonLeading data."<<endl;
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
  TFile *f_spillovers[12];
  TFile *f_fragmentation[12];

  TH2D *result[12][5][4];
  TH2D *result2[12][5][4];
  TH2D *fitbg[12][5][4];

  TH1D *check_old_phi[12][5][4];
  TH1D *check_old_phi_rebin[12][5][4];
  TH1D *check_old_eta[12][5][4];
  TH1D *check_old_eta_rebin[12][5][4];
 
  TH1D *check_new_eta[12][5][4];
  TH1D *check_new_eta_rebin[12][5][4];
  TH1D *check_new_eta_rebin2[12][5][4];
  TH1D *check_new_eta_syst[12][5][4];
  TH1D *check_new_eta_min[12][5][4];
  TH1D *check_new_eta_max[12][5][4];
 
  TH1D *check_new_phi[12][5][4];
  TH1D *check_new_phi_rebin[12][5][4];
  TH1D *check_new_phi_syst[12][5][4];
  TH1D *check_new_phi_max[12][5][4];
  TH1D *check_new_phi_min[12][5][4];

  TH1D *spillover_phi[12][5][4];
  TH1D *spillover_eta[12][5][4];
  

  TH1D *jff_residual_phi[12][5][4];
  TH1D *jff_residual_eta[12][5][4];
  TH1D *jff_residual_eta_combined[12][5][4];
  TH1D *jff_residual_phi_combined[12][5][4];
   

  TH1D *PbPb_pp_phi[12][5][4];
  TH1D *PbPb_pp_phi_syst[12][5][4];
  TH1D *PbPb_pp_eta[12][5][4];
  TH1D *PbPb_pp_eta_syst[12][5][4];

  TH1D *Error_BG_phi[12][5][4];
  TH1D *Error_Spillover_phi[12][5][4];
  TH1D *Error_Relative_phi[12][5][4];
  TH1D *Error_BG_Relative_phi[12][5][4];
  TH1D *Error_BG_Spillover_phi[12][5][4];
  TH1D *Error_BG_Spillover_Relative_phi[12][5][4];
  

  TH1D *Error_BG_eta[12][5][4];
  TH1D *Error_Spillover_eta[12][5][4];
  TH1D *Error_Relative_eta[12][5][4];
  TH1D *Error_BG_Relative_eta[12][5][4];
  TH1D *Error_BG_Spillover_eta[12][5][4];
  TH1D *Error_BG_Spillover_Relative_eta[12][5][4];

  TH1D *Integral_phi_Pt[12][4];
  TH1D *Integral_eta_Pt[12][4];

  TH1D *incl_dist_ref_eta[12][5][4];
  TH1D *incl_dist_ref_phi[12][5][4];

  TGraphAsymmErrors *Integral_eta_syst[12][4];
  TGraphAsymmErrors *Integral_phi_syst[12][4];

  TCanvas *ccheckEta_wide[12];
  TCanvas *ccheckEta[12];
  TCanvas *ccheckPhi[12];
  
  float jff_percent[12][5];

  double Integral_eta[12][5][4];
  double Integral_eta_max[12][5][4];
  double Integral_eta_min[12][5][4];
  double Integral_eta_error_max[12][5][4];
  double Integral_eta_error_min[12][5][4];
  double Integral_eta_error[12][5][4];
 
  Double_t error_A[12][5][4];
  Double_t error_B[12][5][4];
  Double_t error_C[12][5][4];
  Double_t error_D[12][5][4];
  Double_t error_E[12][5][4];
  Double_t error_F[12][5][4];

  double syst_err[12][5][4];
  double Integral_phi[12][5][4];
  double Integral_phi_max[12][5][4];
  double Integral_phi_min[12][5][4];
  double Integral_phi_error_max[12][5][4];
  double Integral_phi_error_min[12][5][4];
  double Integral_phi_error[12][5][4];
 
  TString in_name, plotname, outname, funcname, centlabel, datalabel, pTlabel,checkcanvasnameEta,checkcanvasnamePhi,rawetacanvasname,PbPbminuscanvas_phi,PbPbminuscanvas_eta, checkcanvasnameEta_wide, EtaSpilloverName;

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

  TGraph *Error_up_eta_pT[12][4];
  TGraph *Error_down_eta_pT[12][4];
 
  Double_t syst_error_tot,lphiedge1,lphiedge2,rphiedge1,rphiedge2,letaedge1,letaedge2,retaedge1,retaedge2,most_err,second_err,evalpt, syst_error_bg_c,re,ce,fe, bc,PbPb_err, pp_err, tot_err, check_ymax, check_ymin, dx_eta, dx_phi, err_up, err_down,pp_error,PbPb_error;

  TF1 *gaus1d = new TF1("gaus1d","[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))");

  TF1 *gaus_phi[12][5][4];
  TF1 *gaus_eta[12][5][4];
 
  TLegend *lcheck, *leta;
  
  TLine *linePhi, *lineEta;
 
  float relative_error_percent = TMath::Sqrt(0.05*0.05+0.04*0.04+0.03*0.03);
  
  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
 
  for(int g=gstart; g<gend; g++){

    //  There will only be one big "g-loop".

    switch(g){
    case 0:
      if(Is_NonLeading == kFALSE){
	fnewbg[g] = new TFile("../bg_fits/PbPb_Inclusive_Yield_and_Bkg.root", "READ");
	f_spillovers[g] = new TFile("../spill_over/Inclusive_SpillOvers.root","READ");
	f_fragmentation[g] = new TFile("../jff_residual/Inclusive_PHSube0_JFFResiduals.root","READ");
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
	f_fragmentation[g] = new TFile("../jff_residual/Inclusive_Pythia_JFFResiduals.root","READ");
	if(Spillover_and_pp_subtracted == kTRUE){	
	  fout[g] = new TFile("Inclusive_Data_AllPlots.root", "RECREATE");}
	if(Spillover_and_pp_subtracted == kFALSE){	fout[g] = new TFile("Inclusive_Data_NoSpillOver_AllPlots.root", "RECREATE");}
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
      f_spillovers[g] = new TFile("../spill_over/SubLeading_SpillOvers.root","READ");
      f_fragmentation[g] = new TFile("../jff_residual/SubLeading_PHSube0_JFFResiduals.root","READ");
      datalabel = "SubLeading";
      break;

    case 3:
      fnewbg[g] = new TFile("../bg_fits/pp_SubLeading_Yield_and_Bkg.root", "READ");
      f_fragmentation[g] = new TFile("../jff_residual/SubLeading_Pythia_JFFResiduals.root","READ");
      if(Spillover_and_pp_subtracted ==kTRUE){    fout[g] = new TFile("SubLeading_Data_AllPlots.root", "RECREATE");}
      if(Spillover_and_pp_subtracted ==kFALSE){    fout[g] = new TFile("SubLeading_Data_NoSpillOver_AllPlots.root", "RECREATE");}
      datalabel = "SubLeading";
      break;

    case 4:
      fnewbg[g] = new TFile("../bg_fits/PbPb_Leading_Yield_and_Bkg.root", "READ");
      f_spillovers[g] = new TFile("../spill_over/Leading_SpillOvers.root","READ");
      f_fragmentation[g] = new TFile("../jff_residual/Leading_PHSube0_JFFResiduals.root","READ");
      datalabel = "Leading";
      break;
     
    case 5:
      fnewbg[g] = new TFile("../bg_fits/pp_Leading_Yield_and_Bkg.root", "READ");
      f_fragmentation[g] = new TFile("../jff_residual/Leading_Pythia_JFFResiduals.root","READ");
      if(Spillover_and_pp_subtracted == kTRUE){  fout[g] = new TFile("Leading_Data_AllPlots.root", "RECREATE");}
      if(Spillover_and_pp_subtracted == kFALSE){  fout[g] = new TFile("Leading_Data_NoSpillOver_AllPlots.root", "RECREATE");}
      datalabel = "Leading";
      break;
    default:
      cerr<<"MC is now studied separately"<<endl;
      return -1;
      break;
    }


    

    //----------------------------------------------
    //     Start i and l loops -- set up canvases
    //------------------------------------------------
   
      checkcanvasnameEta = "CheckCanvasEta";
      checkcanvasnameEta+= g;
      ccheckEta[g] = new TCanvas(checkcanvasnameEta," ",10,10,1500,1600);
      ccheckEta[g]->Divide(4,4,0.,0.);
     
      checkcanvasnameEta_wide = "CheckCanvasEtaWide";
      checkcanvasnameEta_wide+= g;
      ccheckEta_wide[g] = new TCanvas(checkcanvasnameEta_wide," ",10,10,1500,1600);
      ccheckEta_wide[g]->Divide(4,4,0.,0.);

      checkcanvasnamePhi = "CheckCanvasPhi";
      checkcanvasnamePhi+= g;      ccheckPhi[g] = new TCanvas(checkcanvasnamePhi," ",10,10,1500,1600);
      ccheckPhi[g]->Divide(4,4,0.,0.);

      //----------------------------------------------------
      //  Start of main i & j loops 
      //-----------------------------------------------------
	
        
    

      for(int i=0; i<4; i++){
	
	for (int j=0; j<4; j++){


	  in_name = make_name("Result_",g,i,j,0,centlabel,pTlabel);


	  result[g][i][j] = (TH2D*)fnewbg[g]->Get(in_name)->Clone(in_name);

	  TString in_name2 = in_name;
	  in_name2.ReplaceAll("Result","SummedResult");
	
	  result2[g][i][j] = (TH2D*)fnewbg[g]->Get(in_name2)->Clone(in_name2);

	  TString fitbg_name = in_name;
	  fitbg_name.ReplaceAll("Result","Fit_bkg");
	  fitbg[g][i][j] = (TH2D*)fnewbg[g]->Get(fitbg_name)->Clone(fitbg_name);
	
	  //-------------------------------
	  //dEta projection
	  //------------------------

	  TString summedchecknameEta = in_name;
	  summedchecknameEta.ReplaceAll("Result","Summed_check_Eta");
	    
	  llimiteta = result[g][i][j]->GetXaxis()->FindBin(-etalim+.001);
	  rlimiteta = result[g][i][j]->GetXaxis()->FindBin(etalim-.001);

	  llimitphi = result[g][i][j]->GetYaxis()->FindBin(-philim+.001);
	  rlimitphi = result[g][i][j]->GetYaxis()->FindBin(philim-.001);
	    

	  check_old_eta[g][i][j] = (TH1D*)result2[g][i][j]->ProjectionX(summedchecknameEta,llimitphi,rlimitphi);
	  dx_eta = check_old_eta[g][i][j]->GetBinWidth(1);
	  check_old_eta[g][i][j]->Scale(1/dx_eta);

	  TString summedchecknameEta_rebin = summedchecknameEta;
	  summedchecknameEta_rebin += "rebin";
	 
	  check_old_eta_rebin[g][i][j] = (TH1D*)Rebin_dEta(check_old_eta[g][i][j]);
	  check_old_eta_rebin[g][i][j]->SetName(summedchecknameEta_rebin);
	 

	  TString newchecknameEta = in_name;
	  newchecknameEta.ReplaceAll("Result", "New_check_Eta");

	  check_new_eta[g][i][j] =result[g][i][j]->ProjectionX(newchecknameEta,llimitphi,rlimitphi);
	  check_new_eta[g][i][j]->Scale(1/dx_eta);
	 
	  TString newchecknameEta_rebin = newchecknameEta;
	  newchecknameEta_rebin += "rebin";

	  if(g<6){
	    if(Is_NonLeading == kFALSE){
	      check_new_eta_rebin[g][i][j] = (TH1D*)Rebin_dEta(check_new_eta[g][i][j]);
	      check_new_eta_rebin[g][i][j]->SetName(newchecknameEta_rebin);
	    }

	    if(Is_NonLeading == kTRUE){
	      check_new_eta_rebin[g][i][j] = (TH1D*) check_old_eta_rebin[g][i][j]->Clone(newchecknameEta_rebin);
	    }
	    
	  }

	  //-------------------------------
	  //   dPhi projection
	  //------------------------------ 
      
	  TString summedchecknamePhi = in_name;
	  summedchecknamePhi.ReplaceAll("Result","Summed_check_Phi");

	  //ccheckPhi[g]->cd(4*i+j+1);
     
	  check_old_phi[g][i][j] = result2[g][i][j]->ProjectionY(summedchecknamePhi,llimiteta,rlimiteta);
	  dx_phi = check_old_phi[g][i][j]->GetBinWidth(1);
	  check_old_phi[g][i][j]->Scale(1/dx_phi);

	 
	  TString summedchecknamePhi_rebin = summedchecknamePhi;
	  summedchecknamePhi_rebin += "rebin";

	  check_old_phi_rebin[g][i][j] =(TH1D*)Rebin_dPhi(check_old_phi[g][i][j]);
	  check_old_phi_rebin[g][i][j]->SetName(summedchecknamePhi_rebin);
	 
	  TString newchecknamePhi = in_name;

	  newchecknamePhi.ReplaceAll("Result", "New_check_Phi");
	  check_new_phi[g][i][j] = result[g][i][j]->ProjectionY(newchecknamePhi,llimiteta,rlimiteta);
	  check_new_phi[g][i][j]->Scale(1/dx_phi);
	  

	  TString newchecknamePhi_rebin = newchecknamePhi;
	  newchecknamePhi_rebin += "rebin";


	  if(g<6){
	    check_new_phi_rebin[g][i][j] = (TH1D*)Rebin_dPhi(check_new_phi[g][i][j]);
	    check_new_phi_rebin[g][i][j]->SetName(newchecknamePhi_rebin);
	  }

	  if(j==0){

	    for(int k = 0; k<4; k++){
	      if(g%2!=0&&k>0)continue;

	      TString jff_residual_name_eta = make_name("GenJet_GenTrack_Sube0_JFF_Residual_Eta_Merged_",g,i,k,0,centlabel,pTlabel);
	      jff_residual_name_eta.ReplaceAll("PbPb_","");
	      jff_residual_name_eta.ReplaceAll("pp_","");
	      jff_residual_name_eta.ReplaceAll("Pt100_Pt300_","");

	      if(g%2!=0){
		jff_residual_name_eta = make_name("GenJet_GenTrack_JFF_Residual_Eta_Merged_",g,i,3,0,centlabel,pTlabel);
		jff_residual_name_eta.ReplaceAll("PbPb_","");
		jff_residual_name_eta.ReplaceAll("pp_","");
		jff_residual_name_eta.ReplaceAll("Pt100_Pt300_","");
	      }
	  
	      if(g==2||g==3)jff_residual_name_eta.ReplaceAll("Eta","EtaSubLeading");
	      if(g==4||g==5)jff_residual_name_eta.ReplaceAll("Eta","EtaLeading");

	      cout<<g<<" "<<jff_residual_name_eta<<endl;


	      jff_residual_eta[g][i][k] = (TH1D*)f_fragmentation[g]->Get(jff_residual_name_eta)->Clone(jff_residual_name_eta);

	      TString jff_residual_name_phi = jff_residual_name_eta;
	      jff_residual_name_phi.ReplaceAll("Eta","Phi");
	 
	      jff_residual_phi[g][i][k] = (TH1D*)f_fragmentation[g]->Get(jff_residual_name_phi)->Clone(jff_residual_name_phi);
	    }
	 


	    TString jff_residual_name_eta_combined = "JFF_Residual_Eta_Combined"; jff_residual_name_eta_combined+=g;  jff_residual_name_eta_combined+=i;
	    TString jff_residual_name_phi_combined = "JFF_Residual_Phi_Combined"; jff_residual_name_phi_combined+=g;  jff_residual_name_phi_combined+=i;
					


	    jff_residual_eta_combined[g][i][0] = (TH1D*)jff_residual_eta[g][i][0]->Clone(jff_residual_name_eta_combined);
	    jff_residual_phi_combined[g][i][0] = (TH1D*)jff_residual_phi[g][i][0]->Clone(jff_residual_name_phi_combined);

	    if(g%2==0){
	      jff_residual_eta_combined[g][i][0]->Add( jff_residual_eta[g][i][1]);
	      jff_residual_eta_combined[g][i][0]->Add( jff_residual_eta[g][i][2]);
	      jff_residual_eta_combined[g][i][0]->Add( jff_residual_eta[g][i][3]);

	      jff_residual_eta_combined[g][i][0]->Scale(1./4);
	  
	      jff_residual_phi_combined[g][i][0]->Add( jff_residual_phi[g][i][1]);
	      jff_residual_phi_combined[g][i][0]->Add( jff_residual_phi[g][i][2]);
	      jff_residual_phi_combined[g][i][0]->Add( jff_residual_phi[g][i][3]);

	      jff_residual_phi_combined[g][i][0]->Scale(1./4);
	    }
	  }

	  check_new_eta_rebin[g][i][j]->Add( jff_residual_eta_combined[g][i][0],-1.);
	  check_new_phi_rebin[g][i][j]->Add( jff_residual_phi_combined[g][i][0],-1.);

	  jff_percent[g][i] = TMath::Abs(jff_residual_eta[g][i][0] ->GetBinContent(jff_residual_eta[g][i][0] ->FindBin(0)) - jff_residual_eta_combined[g][i][0]->GetBinContent(jff_residual_eta[g][i][0] ->FindBin(0)))/TMath::Abs(jff_residual_eta_combined[g][i][0]->GetBinContent(jff_residual_eta[g][i][0] ->FindBin(0)));

	  for(int k = 1; k<4; k++){
	    if(g%2==0&&(TMath::Abs(jff_residual_eta[g][i][k] ->GetBinContent(jff_residual_eta[g][i][0] ->FindBin(0))- jff_residual_eta_combined[g][i][0]->GetBinContent(jff_residual_eta[g][i][0] ->FindBin(0)))/TMath::Abs(jff_residual_eta_combined[g][i][0]->GetBinContent(jff_residual_eta[g][i][0] ->FindBin(0)))>jff_percent[g][i])){
	      jff_percent[g][i] = TMath::Abs(jff_residual_eta[g][i][k] ->GetBinContent(jff_residual_eta[g][i][0] ->FindBin(0))- jff_residual_eta_combined[g][i][0]->GetBinContent(jff_residual_eta[g][i][0] ->FindBin(0)))/TMath::Abs(jff_residual_eta_combined[g][i][0]->GetBinContent(jff_residual_eta[g][i][0] ->FindBin(0)));
	    }
	  }

	 	  	
	  //SCALING FOR pT 4-8.
	    
	  if(i==3){
	    check_old_eta_rebin[g][i][j]->Scale(1/4.);
	    check_new_eta_rebin[g][i][j]->Scale(1/4.);
	    check_old_phi_rebin[g][i][j]->Scale(1/4.);
	    check_new_phi_rebin[g][i][j]->Scale(1/4.);
	  }


	  //----------------------------------------------------
	  //   Calculate Systematic "BG Error" for ALL Check Plots
	  //-----------------------------------------------------



	  if(g<6){
	  
	    //error_A[g][i][j]: calculated as bin error at dPhi = 0, propogated over all 2D
	    
	    error_A[g][i][j] = fitbg[g][i][j]->GetBinError(26,20)*20/dx_phi;

	    //  Systematic source B:  calculated as the average of the 2 greatest out of the 4 most peripheral points in dPhi
	     
	    lphiedge1 = check_new_phi_rebin[g][i][j]->GetBinContent(1);
	    lphiedge2 = check_new_phi_rebin[g][i][j]->GetBinContent(2);
	    rphiedge1 = check_new_phi_rebin[g][i][j]->GetBinContent(16);
	    rphiedge2 = check_new_phi_rebin[g][i][j]->GetBinContent(17);

	    most_err = 0.;
	    second_err = 0.;
	   

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

	    error_B[g][i][j] = (most_err+second_err)/2.;


	    //   Systematic source C: calculated as the average of the 2 greatest out of the 4 most peripheral points in dEta

	    letaedge1 = check_new_eta_rebin[g][i][j]->GetBinContent(3);
	    letaedge2 = check_new_eta_rebin[g][i][j]->GetBinContent(4);
	    retaedge1 = check_new_eta_rebin[g][i][j]->GetBinContent(nbounds_eta-4);
	    retaedge2 = check_new_eta_rebin[g][i][j]->GetBinContent(nbounds_eta-3);
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
	    error_C[g][i][j] = (most_err+second_err)/2.;
	    
	    //  Systematic source D: calculated as the average of the absolute values of the deviation of the contents of the 1.5<|dEta|<2.0

	    error_D[g][i][j] = (TMath::Abs(check_new_eta_rebin[g][i][j]->GetBinContent(2))+TMath::Abs(check_new_eta_rebin[g][i][j]->GetBinContent(18)))/2.;
	 
	    //  Systematic source E: calculated as the the average bin-to-bin deviation between the 1.5<|dEta|<2.0 and 2.0<|dEta|<2.5 bins

	    error_E[g][i][j] = (TMath::Abs(check_new_eta_rebin[g][i][j]->GetBinContent(1)-check_new_eta_rebin[g][i][j]->GetBinContent(2))+TMath::Abs(check_new_eta_rebin[g][i][j]->GetBinContent(18)-check_new_eta_rebin[g][i][j]->GetBinContent(19)))/2.;

	    //  error_E[g][i][j] = (check_new_eta_rebin[g][i][j]->GetBinContent(1)-check_new_eta_rebin[g][i][j]->GetBinContent(2)+check_new_eta_rebin[g][i][j]->GetBinContent(19)-check_new_eta_rebin[g][i][j]->GetBinContent(18))/2.;
	      
	    //  With these corrections done, add all errors in quadrature:
	    //   syst_err[g][i][j]= TMath::Sqrt(error_B[g][i][j]*syst_checkB+error_C[g][i][j]*error_C[g][i][j]+error_D[g][i][j]*error_D[g][i][j]+error_E[g][i][j]*error_E[g][i][j]+syst_checkF*syst_checkF);

	    //D, E, & F only
	    syst_err[g][i][j]= TMath::Sqrt(error_D[g][i][j]*error_D[g][i][j]+error_F[g][i][j]*error_F[g][i][j]);

	      
	    // ***note: syst_err[g][i][j] accounts for all "BG error"
	   
	    if(syst_err[g][i][j]<0.01){syst_err[g][i][j]=0.01;}
	      
	
	    //----------------
	    // over-rides
	    //----------------

	    if(g==0){
	      switch(i){
	      case 0:
		if(j==0){syst_err[g][i][j]=0.06;}
		if(j==1){syst_err[g][i][j]=0.14;}
		if(j==2){syst_err[g][i][j]=0.19;}
		if(j==3){syst_err[g][i][j]=0.21;}
		break;
	      case 1:
		if(j==0){syst_err[g][i][j]=0.03;}
		if(j==1){syst_err[g][i][j]=0.04;}
		if(j==2){syst_err[g][i][j]=0.05;}
		if(j==3){syst_err[g][i][j]=0.05;}
		break;
	      case 2:
		if(j==0){syst_err[g][i][j]=0.02;}
		if(j==1){syst_err[g][i][j]=0.02;}
		if(j==2){syst_err[g][i][j]=0.03;}
		if(j==3){syst_err[g][i][j]=0.03;}
		break;
	      case 3:
		syst_err[g][i][j] = 0.02;
		break;
	      }
	    }
	    if(g==2){
	      switch(i){
	      case 0:
		if(j==0){syst_err[g][i][j]=0.14;}
		if(j==1){syst_err[g][i][j]=0.16;}
		if(j==2){syst_err[g][i][j]=0.22;}
		if(j==3){syst_err[g][i][j]=0.27;}
		break;
	      case 1:
		if(j==0){syst_err[g][i][j]=0.05;}
		if(j==1){syst_err[g][i][j]=0.06;}
		if(j==2){syst_err[g][i][j]=0.09;}
		if(j==3){syst_err[g][i][j]=0.10;}
		break;
	      case 2:
		if(j==0){syst_err[g][i][j]=0.03;}
		if(j==1){syst_err[g][i][j]=0.03;}
		if(j==2){syst_err[g][i][j]=0.04;}
		if(j==3){syst_err[g][i][j]=0.05;}
		break;
	      case 3:
		syst_err[g][i][j] = 0.02;
		break;
	      }
	    }
	    if(g==4){
	      switch(i){
	      case 0:
		if(j==0){syst_err[g][i][j]=0.13;}
		if(j==1){syst_err[g][i][j]=0.16;}
		if(j==2){syst_err[g][i][j]=0.20;}
		if(j==3){syst_err[g][i][j]=0.24;}
		break;
	      case 1:
		if(j==0){syst_err[g][i][j]=0.05;}
		if(j==1){syst_err[g][i][j]=0.06;}
		if(j==2){syst_err[g][i][j]=0.07;}
		if(j==3){syst_err[g][i][j]=0.08;}
		break;
	      case 2:
		if(j==0){syst_err[g][i][j]=0.02;}
		if(j==1){syst_err[g][i][j]=0.03;}
		if(j==2){syst_err[g][i][j]=0.04;}
		if(j==3){syst_err[g][i][j]=0.05;}
		break;
	      case 3:
		syst_err[g][i][j] = 0.02;
		break;
	      }
	    }	
	
	    //=------------------------------------------------------------------
	    //  Save individual histos for  different sources of error
	    //  ***Section above defines BG error.  Spillover error is calculated from fits in the loop below.  Relative error is Sqrt(3*05*.05)*bincontent
	    //
	    //   Also create systematic error histograms for use in plotting
	    //-------------------------------------------------------------------
	  
	    if(Spillover_and_pp_subtracted == kTRUE&&(g==0||g==2||g==4)){
	      EtaSpilloverName = in_name;
	      EtaSpilloverName.ReplaceAll("Result_","SpillOvers_Eta_");
	      if(g==2){	  EtaSpilloverName.ReplaceAll("Eta_","Eta_SubLeading_");        }
	      if(g==4){	  EtaSpilloverName.ReplaceAll("Eta_","Eta_Leading_");        }
	      EtaSpilloverName.ReplaceAll("PbPb_","RecoJet_GenTrack_NoSube0_"); 
	      EtaSpilloverName.ReplaceAll("Pt100_Pt300_","");

	      spillover_eta[g][i][j] = (TH1D*)f_spillovers[g]->Get(EtaSpilloverName)->Clone(EtaSpilloverName);
	      
	      TString PhiSpilloverName = EtaSpilloverName;
	      PhiSpilloverName.ReplaceAll("Eta","Phi");

	      spillover_phi[g][i][j] = (TH1D*)f_spillovers[g]->Get(PhiSpilloverName)->Clone(PhiSpilloverName);
	    
	    }
	    
	    TString Error_BG_eta_name = in_name;
	    Error_BG_eta_name.ReplaceAll("Result","Error_BG_eta");
	    Error_BG_eta_name.ReplaceAll("pp","PbPb");

	    TString Error_Spillover_eta_name = in_name;
	    Error_Spillover_eta_name.ReplaceAll("Result","Error_Spillover_eta");
	    Error_Spillover_eta_name.ReplaceAll("pp","PbPb");
	      	   
	    TString Error_Relative_eta_name = in_name;
	    Error_Relative_eta_name.ReplaceAll("Result","Error_Relative_eta");
	    Error_Relative_eta_name.ReplaceAll("pp","PbPb");

	    TString Error_BG_Relative_eta_name = in_name;
	    Error_BG_Relative_eta_name.ReplaceAll("Result","Error_BG_Relative_eta");
	    Error_BG_Relative_eta_name.ReplaceAll("pp","PbPb");


	    TString Error_BG_Spillover_eta_name = in_name;
	    Error_BG_Spillover_eta_name.ReplaceAll("Result","Error_BG_Spillover_eta");
	    Error_BG_Spillover_eta_name.ReplaceAll("pp","PbPb");

	    TString Error_BG_Spillover_Relative_eta_name = in_name;
	    Error_BG_Spillover_Relative_eta_name.ReplaceAll("Result","Error_BG_Spillover_Relative_eta");
	    Error_BG_Spillover_Relative_eta_name.ReplaceAll("pp","PbPb");
 
	

	    Error_BG_eta[g][i][j] = new TH1D (Error_BG_eta_name,"",nbounds_eta-1,bin_bounds_eta);
	    Error_Spillover_eta[g][i][j]= new TH1D (Error_Spillover_eta_name,"",nbounds_eta-1,bin_bounds_eta);
	    Error_Relative_eta[g][i][j]= new TH1D (Error_Relative_eta_name,"",nbounds_eta-1,bin_bounds_eta);
	    Error_BG_Relative_eta[g][i][j]= new TH1D (Error_BG_Relative_eta_name,"",nbounds_eta-1,bin_bounds_eta);
	    Error_BG_Spillover_eta[g][i][j]= new TH1D (Error_BG_Spillover_eta_name,"",nbounds_eta-1,bin_bounds_eta);
	    Error_BG_Spillover_Relative_eta[g][i][j]= new TH1D (Error_BG_Spillover_Relative_eta_name,"",nbounds_eta-1,bin_bounds_eta);
  
	    
	    TString newchecknameEtasyst = newchecknameEta;
	    newchecknameEtasyst.ReplaceAll("New_check","Syst_error");
	    
	    check_new_eta_syst[g][i][j] = new TH1D(newchecknameEtasyst,"",nbounds_eta-1,bin_bounds_eta);

	    if(Spillover_and_pp_subtracted==kFALSE||(g==1||g==3||g==5)){ce = 0.;}

	    for(int m=1; m<nbounds_eta;m++){
	  
	      re = relative_error_percent*check_new_eta_rebin[g][i][j]->GetBinContent(m);

	      fe = jff_percent[g][i]*jff_residual_eta_combined[g][i][0]->GetBinContent(m);
		  
	 
	      if((g==0||g==2||g==4)&&Spillover_and_pp_subtracted == kTRUE){ce = spillover_eta[g][i][j]->GetBinContent(m)/2; }

	      syst_error_bg_c =TMath::Sqrt( syst_err[g][i][j]*syst_err[g][i][j]+ce*ce);
		
	      syst_error_tot =TMath::Sqrt(syst_err[g][i][j]*syst_err[g][i][j]+ce*ce+re*re+fe*fe);
	
	      //	syst_error_tot = syst_err[g][i][j];
	
	
	      bc = check_new_eta_rebin[g][i][j]->GetBinContent(m);
	      check_new_eta_syst[g][i][j]->SetBinContent(m,bc);
	      check_new_eta_syst[g][i][j]->SetBinError(m,syst_error_tot);
	      //check_new_eta_syst[g][i][j]->SetBinError(m,ce);
	
	      Error_BG_eta[g][i][j]->SetBinContent(m,syst_err[g][i][j]);
	      Error_Spillover_eta[g][i][j]->SetBinContent(m,ce);
	      Error_Relative_eta[g][i][j]->SetBinContent(m,re);
	      Error_BG_Relative_eta[g][i][j]->SetBinError(m,sqrt(re*re+fe*fe+syst_err[g][i][j]*syst_err[g][i][j]));
	      Error_BG_Relative_eta[g][i][j]->SetBinContent(m,0.);
	      Error_BG_Spillover_eta[g][i][j]->SetBinContent(m,syst_error_bg_c);
	      Error_BG_Spillover_Relative_eta[g][i][j]->SetBinContent(m,syst_error_tot);
  	    }

	  	    	  
	    TString Error_BG_phi_name = in_name;
	    Error_BG_phi_name.ReplaceAll("Result","Error_BG_phi");
	    Error_BG_phi_name.ReplaceAll("pp","PbPb");

	    TString Error_Spillover_phi_name = in_name;
	    Error_Spillover_phi_name.ReplaceAll("Result","Error_Spillover_phi");
	    Error_Spillover_phi_name.ReplaceAll("pp","PbPb");

	    TString Error_Relative_phi_name = in_name;
	    Error_Relative_phi_name.ReplaceAll("Result","Error_Relative_phi");
	    Error_Relative_phi_name.ReplaceAll("pp","PbPb");

	    TString Error_BG_Relative_phi_name = in_name;
	    Error_BG_Relative_phi_name.ReplaceAll("Result","Error_BG_Relative_phi");
	    Error_BG_Relative_phi_name.ReplaceAll("pp","PbPb");

	    TString Error_BG_Spillover_phi_name = in_name;
	    Error_BG_Spillover_phi_name.ReplaceAll("Result","Error_BG_Spillover_phi");
	    Error_BG_Spillover_phi_name.ReplaceAll("pp","PbPb");

	    TString Error_BG_Spillover_Relative_phi_name = in_name;
	    Error_BG_Spillover_Relative_phi_name.ReplaceAll("Result","Error_BG_Spillover_Relative_phi");
	    Error_BG_Spillover_Relative_phi_name.ReplaceAll("pp","PbPb");
 

	    Error_BG_phi[g][i][j] = new TH1D (Error_BG_phi_name,"",nbounds_phi-1,bin_bounds_phi);
	    Error_Spillover_phi[g][i][j]= new TH1D (Error_Spillover_phi_name,"",nbounds_phi-1,bin_bounds_phi);
	    Error_Relative_phi[g][i][j]= new TH1D (Error_Relative_phi_name,"",nbounds_phi-1,bin_bounds_phi);
	    Error_BG_Relative_phi[g][i][j]= new TH1D (Error_BG_Relative_phi_name,"",nbounds_phi-1,bin_bounds_phi);
	    Error_BG_Spillover_phi[g][i][j]= new TH1D (Error_BG_Spillover_phi_name,"",nbounds_phi-1,bin_bounds_phi);
	    Error_BG_Spillover_Relative_phi[g][i][j]= new TH1D (Error_BG_Spillover_Relative_phi_name,"",nbounds_phi-1,bin_bounds_phi);
  


	    TString newchecknamePhisyst = newchecknamePhi;
	    newchecknamePhisyst.ReplaceAll("New_check","Syst_error");
	     
	    check_new_phi_syst[g][i][j] = new TH1D(newchecknamePhisyst,"",nbounds_phi-1,bin_bounds_phi);

	      
	    //   Fill histograms


	    for(int m=1; m<nbounds_phi;m++){

	      re = relative_error_percent*check_new_phi_rebin[g][i][j]->GetBinContent(m);

	      fe = jff_percent[g][i]*jff_residual_phi_combined[g][i][0]->GetBinContent(m);

	      if((g==0||g==2||g==4)&&Spillover_and_pp_subtracted==kTRUE){ ce = spillover_phi[g][i][j]->GetBinContent(m)/2.;}

	      syst_error_bg_c =TMath::Sqrt(syst_err[g][i][j]*syst_err[g][i][j]+ce*ce);
		
	      syst_error_tot =TMath::Sqrt(syst_err[g][i][j]*syst_err[g][i][j]+ce*ce+re*re+fe*fe);
		
	      //	syst_error_tot = syst_err[g][i][j];
		
	      bc = check_new_phi_rebin[g][i][j]->GetBinContent(m);

	      check_new_phi_syst[g][i][j]->SetBinContent(m,bc);
	      check_new_phi_syst[g][i][j]->SetBinError(m,syst_error_tot);
	 		
	      Error_BG_phi[g][i][j]->SetBinContent(m,syst_err[g][i][j]);
	      Error_Spillover_phi[g][i][j]->SetBinContent(m,ce);
	      Error_Relative_phi[g][i][j]->SetBinContent(m,re);
	      Error_BG_Relative_phi[g][i][j]->SetBinError(m,sqrt(re*re+fe*fe+syst_err[g][i][j]*syst_err[g][i][j]));
	      Error_BG_Relative_phi[g][i][j]->SetBinContent(m,0.);
	      Error_BG_Spillover_phi[g][i][j]->SetBinContent(m,syst_error_bg_c);
	      Error_BG_Spillover_Relative_phi[g][i][j]->SetBinContent(m,syst_error_tot);
	    }
   
	    float val_l, val_r, err_l, err_r;

	    float nbins2 = check_new_eta_rebin[g][i][j]->GetNbinsX()+1;
	

	  }  // All of that was a GIANT Data-Only loop.



	  //------------------------------------------------------------------
	  //      YIELD PLOTS:  First, settings in both dEta and dPhi together
	  //-------------------------------------------------------------------
	 
	  ccheckEta_wide[g]->cd(4*i+j+1);

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
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(10);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(10);
	    check_new_eta_syst[g][i][j]->SetFillColor(90);
	    check_new_phi_syst[g][i][j]->SetFillColor(90);
	    break;
	  case 1:
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(4);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(4);
	    break;
	  case 2: 
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(34);
	    check_new_eta_rebin[g][i][j]->SetMarkerSize(2);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(34);
	    check_new_phi_rebin[g][i][j]->SetMarkerSize(2);
	    check_new_eta_syst[g][i][j]->SetFillColor(30);
	    check_new_phi_syst[g][i][j]->SetFillColor(30);
	    break;
	  case 3:
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(28);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(28);
	    break;  
	  case 4: 
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(21);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(21);
	    check_new_eta_syst[g][i][j]->SetFillColor(kOrange-2);
	    check_new_phi_syst[g][i][j]->SetFillColor(kOrange-2);
	    break;
	  case 5:
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(28);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(28);
	    break;
	  case 6: 
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(10);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(10);
	    break;
	  case 7:
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(4);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(4);
	    break;
	  case 8: 
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(34);
	    check_new_eta_rebin[g][i][j]->SetMarkerSize(2);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(34);
	    check_new_phi_rebin[g][i][j]->SetMarkerSize(2);
	    break;
	  case 9:
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(28);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(28);
	    break;  
	  case 10: 
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(21);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(21);
	    break;
	  case 11:
	    check_new_eta_rebin[g][i][j]->SetMarkerStyle(25);
	    check_new_phi_rebin[g][i][j]->SetMarkerStyle(25);
	    break;
	  default:
	    break;
	  }
	 
	  check_new_eta_rebin[g][i][j]->SetMarkerColor(kBlack);
	  check_new_eta_rebin[g][i][j]->SetLineColor(kBlack);
	  
	  check_new_phi_rebin[g][i][j]->SetMarkerColor(kBlack);
	  check_new_phi_rebin[g][i][j]->SetLineColor(kBlack);
	
	
	  if((g==0||g==2||g==4)&&Spillover_and_pp_subtracted==kTRUE&&i<3){
	        
	    check_new_eta_rebin[g][i][j]->Add(spillover_eta[g][i][j],-1.);
	    check_old_eta_rebin[g][i][j]->Add(spillover_eta[g][i][j],-1.);
	    check_new_eta_syst[g][i][j]->Add(spillover_eta[g][i][j],-1.);
   
	    check_new_phi_rebin[g][i][j]->Add(spillover_phi[g][i][j],-1.);
	    check_old_phi_rebin[g][i][j]->Add(spillover_phi[g][i][j],-1.);
	    check_new_phi_syst[g][i][j]->Add(spillover_phi[g][i][j],-1.);
	     
	  }

	  //------------------
	  // dEta plotting
	  //-------------------

	  check_new_eta_rebin[g][i][j]->SetMaximum(check_ymax);
	  check_new_eta_rebin[g][i][j]->SetMinimum(check_ymin);

   
	  check_new_eta_rebin[g][i][j]->SetMarkerColor(1);
	  check_new_eta_rebin[g][i][j]->SetAxisRange(-2.5+.2,2.5-.2,"x");

	  check_new_eta_rebin[g][i][j]->Draw();	 
	  
      
	  check_new_eta_rebin[g][i][j]->GetYaxis()->SetLabelSize(tstitle);
	  check_new_eta_rebin[g][i][j]->GetXaxis()->SetLabelSize(tstitle);
	  check_new_eta_rebin[g][i][j]->GetXaxis()->SetTitle("#Delta#eta");
	  check_new_eta_rebin[g][i][j]->GetXaxis()->SetTitleSize(tstitle);
	  check_new_eta_rebin[g][i][j]->GetXaxis()->SetTitleOffset(xoffset);
	  check_new_eta_rebin[g][i][j]->GetYaxis()->SetTitle("Y#equiv1/N_{jet} d^{2}N/(d#Delta#etadp_{T}) (GeV/c)^{-1}");
	  check_new_eta_rebin[g][i][j]->GetYaxis()->SetTitleOffset(yoffset);
	  check_new_eta_rebin[g][i][j]->GetYaxis()->SetTitleSize(tstitle);

 
	  check_new_eta_rebin[g][i][j]->GetXaxis()->CenterTitle();
	  check_new_eta_rebin[g][i][j]->GetYaxis()->CenterTitle();
	  if(i<3){
	    check_new_eta_rebin[g][i][j]->GetXaxis()->SetLabelSize(0.);
	  }
	  if(j>0){
	    check_new_eta_rebin[g][i][j]->GetYaxis()->SetTitleSize(0.0);
	    check_new_eta_rebin[g][i][j]->GetYaxis()->SetLabelSize(0.0);
	  }

	  if(g==0||g==2||g==4){
	    check_new_eta_syst[g][i][j]->SetMarkerSize(0);
	    check_new_eta_syst[g][i][j]->Draw("same e2");
	   
	  }
	 
	  
	    if(g==0||g==2||g==4){
	    check_old_eta_rebin[g][i][j]->SetMarkerStyle(10);
	    check_old_eta_rebin[g][i][j]->SetMarkerSize(1);
	    check_old_eta_rebin[g][i][j]->SetMarkerColor(kRed);
	    //   check_old_eta_rebin[g][i][j]->Draw("same");
	    }
	  
	    
	  if(g==1||g==3||g==5){
	    check_new_eta_syst[g-1][i][j]->Draw("same e2");
	    check_new_eta_rebin[g-1][i][j]->Draw("same");
	  }

	  check_new_eta_rebin[g][i][j]->Draw("same");


	  drawlabels_4by4(g,i,j,datalabel);
	      
	  lineEta = new TLine(-2.5,0,2.5,0);
	  lineEta->SetLineStyle(2);
	  lineEta->Draw("same");


	  //---------------
	  // Narrower axis range for presentation
	  //----------------

	  
	  ccheckEta[g]->cd(4*i+j+1);

	  TString narrowetaname = "NarrowEta";
	  narrowetaname+=g;
	  narrowetaname+=i;
	  narrowetaname+=j;

	  check_new_eta_rebin2[g][i][j] = (TH1D*)check_new_eta_rebin[g][i][j]->Clone(narrowetaname);
	      
	  check_new_eta_rebin2[g][i][j]->SetAxisRange(-1.5+.1,1.5-.1,"x");

	  check_new_eta_rebin2[g][i][j]->Draw();	 
	  
	  if(g==0||g==2||g==4){
	    check_new_eta_syst[g][i][j]->Draw("same e2");
	  }

	  if(g==1||g==3||g==5){
	    check_new_eta_syst[g-1][i][j]->Draw("same e2");
	    check_new_eta_rebin[g-1][i][j]->Draw("same");
	  }
	  check_new_eta_rebin[g][i][j]->Draw("same");

	  drawlabels_4by4(g,i,j,datalabel);
	 
	  lineEta = new TLine(-1.5,0,1.5,0);
	  lineEta->SetLineStyle(2);
	  lineEta->Draw("same");

	  //----------------------------
	  //   Now dPhi
	  //-------------------------------


	  ccheckPhi[g]->cd(4*i+j+1);

	
	  check_new_phi_rebin[g][i][j]->SetMaximum(check_ymax);
	  check_new_phi_rebin[g][i][j]->SetMinimum(check_ymin);
	    
	  check_new_phi_rebin[g][i][j]->Draw("same");   
      
 
	  check_new_phi_rebin[g][i][j]->GetXaxis()->CenterTitle();
	  check_new_phi_rebin[g][i][j]->GetYaxis()->CenterTitle();
	  check_new_phi_rebin[g][i][j]->GetYaxis()->SetLabelSize(tstitle);
	  check_new_phi_rebin[g][i][j]->GetXaxis()->SetLabelSize(tstitle);
	  check_new_phi_rebin[g][i][j]->GetXaxis()->SetTitle("#Delta#phi");
	  check_new_phi_rebin[g][i][j]->GetXaxis()->SetTitleSize(tstitle);
	  check_new_phi_rebin[g][i][j]->GetXaxis()->SetTitleOffset(xoffset);
	  check_new_phi_rebin[g][i][j]->GetYaxis()->SetTitle("Y#equiv1/N_{jet} d^{2}N/(d#Delta#phidp_{T}) (GeV/c)^{-1}");
	  check_new_phi_rebin[g][i][j]->GetYaxis()->SetTitleOffset(yoffset);
	  check_new_phi_rebin[g][i][j]->GetYaxis()->SetTitleSize(tstitle);
	  if(i<3){
	    check_new_phi_rebin[g][i][j]->GetXaxis()->SetLabelSize(0.);
	  }
	  if(j>0){
	
	    check_new_phi_rebin[g][i][j]->GetYaxis()->SetTitleSize(0.0);
	    check_new_phi_rebin[g][i][j]->GetYaxis()->SetLabelSize(0.0);
	  }
	    
	  if(g==0||g==2||g==4){
	    check_new_phi_syst[g][i][j]->SetMarkerSize(0);
	    check_new_phi_syst[g][i][j]->Draw("same e2");
	  }

	  if(g==1||g==3||g==5){
	    check_new_phi_syst[g-1][i][j]->Draw("same e2");
	    check_new_phi_rebin[g-1][i][j]->Draw("same");
	  }

	  
	  check_new_phi_rebin[g][i][j]->Draw("same");


	  drawlabels_4by4(g,i,j,datalabel);

	  linePhi = new TLine(-1.51,0,1.51,0);
	  linePhi->SetLineStyle(2);
	  linePhi->Draw("same");

	 
	
	  //---------------------------------------------------------------------------
	  // WRITE ALL PLOTS (associated with correlations & correlations syst. error)
	  //---------------------------------------------------------------------------

	

	  if(g==1||g==3||g==5){
	   	     
	    Error_BG_phi[g-1][i][j]->Write();
	    Error_Spillover_phi[g-1][i][j]->Write();
	    Error_Relative_phi[g-1][i][j]->Write();
	    Error_BG_Relative_phi[g-1][i][j]->Write();
	    Error_BG_Spillover_phi[g-1][i][j]->Write();
	    Error_BG_Spillover_Relative_phi[g-1][i][j]->Write();
	      	      
	    Error_BG_eta[g-1][i][j]->Write();
	    Error_Spillover_eta[g-1][i][j]->Write();
	    Error_Relative_eta[g-1][i][j]->Write();
	    Error_BG_Relative_eta[g-1][i][j]->Write();
	    Error_BG_Spillover_eta[g-1][i][j]->Write();
	    Error_BG_Spillover_Relative_eta[g-1][i][j]->Write();
	    
	    check_new_eta_rebin[g-1][i][j]->Write();
	    check_new_phi_rebin[g-1][i][j]->Write();
	    check_old_eta_rebin[g-1][i][j]->Write();
	    check_old_phi_rebin[g-1][i][j]->Write();
	    	    
 
	    check_new_eta_syst[g-1][i][j]->Write();
	    check_new_phi_syst[g-1][i][j]->Write();
	   

	    check_new_eta_rebin[g][i][j]->Write();
	    check_new_phi_rebin[g][i][j]->Write();
	    check_old_eta_rebin[g][i][j]->Write();
	    check_old_phi_rebin[g][i][j]->Write();
      	  
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
	      
	      

	    PbPb_pp_eta[g][i][j] = (TH1D*)check_new_eta_rebin[g-1][i][j]->Clone(PbPbppname_eta);
	    PbPb_pp_eta[g][i][j]->Add(check_new_eta_rebin[g][i][j],-1.);
	      
	    PbPb_pp_eta_syst[g][i][j] = (TH1D*)PbPb_pp_eta[g][i][j]->Clone(PbPbppname_eta_syst);
	      

	    
	    for(int k=1;k<check_new_eta_rebin[g][i][j]->GetNbinsX()+1;k++){
		
	      PbPb_err = check_new_eta_syst[g-1][i][j]->GetBinError(k);
	      pp_err = check_new_eta_syst[g][i][j]->GetBinError(k);

	      if(isnan(pp_err) || isinf(pp_err) ){pp_err=0.;}
	
	      tot_err = TMath::Sqrt(PbPb_err*PbPb_err+pp_err*pp_err);

	      PbPb_pp_eta_syst[g][i][j]->SetBinError(k,tot_err);
		
	    }
	      
	      
	    check_new_eta_max[g][i][j] = (TH1D*)check_new_eta_rebin[g-1][i][j]->Clone(PbPbppname_eta_max);
	    check_new_eta_max[g][i][j]->Add(Error_BG_Spillover_eta[g-1][i][j]);

	    check_new_eta_min[g][i][j] = (TH1D*)check_new_eta_rebin[g-1][i][j]->Clone(PbPbppname_eta_min);
	    check_new_eta_min[g][i][j]->Add(Error_BG_Spillover_eta[g-1][i][j],-1.);
	

	    if(i==0){
	      TString IntegratedYieldname_eta = in_name;
	      IntegratedYieldname_eta.ReplaceAll("Result_pp","Integrated_Yield_Eta");
	      IntegratedYieldname_eta.ReplaceAll("Pt100_Pt300_","");
	      IntegratedYieldname_eta.ReplaceAll("_TrkPt1_TrkPt2","");
	      Integral_eta_Pt[g][j] = new TH1D(IntegratedYieldname_eta,"",4,pTbins);
	    }
	     	      

	    llimiteta = PbPb_pp_eta[g][i][j]->FindBin(-etalim+.001);
	    rlimiteta = PbPb_pp_eta[g][i][j]->FindBin(etalim-.001);

	     
	  
	    
	    Integral_eta[g][i][j] = check_new_eta_rebin[g-1][i][j]->IntegralAndError(llimiteta,rlimiteta,PbPb_error,"width");
	      
	     

	    Integral_eta_error[g][i][j] = TMath::Sqrt(PbPb_error*PbPb_error+pp_error*pp_error);  

	    for(int k=llimiteta;k<rlimiteta+1;k++){
	      if((check_new_eta_min[g][i][j]->GetBinContent(k)<0)&&(check_new_eta_rebin[g-1][i][j]->GetBinContent(k)>=0)){
		check_new_eta_min[g][i][j]->SetBinContent(k,0.);
	      }
	    }
	     	     	      
	    Integral_eta_max[g][i][j] = check_new_eta_max[g][i][j]->Integral(llimiteta,rlimiteta,"width")-Integral_eta[g][i][j];
	    
	    Integral_eta_min[g][i][j] = check_new_eta_min[g][i][j]->Integral(llimiteta,rlimiteta,"width")-Integral_eta[g][i][j];

	
	    Integral_eta_error_max[g][i][j] = TMath::Sqrt(Integral_eta_max[g][i][j]*Integral_eta_max[g][i][j]+Integral_eta[g][i][j]*Integral_eta[g][i][j]*relative_error_percent*relative_error_percent);
	    Integral_eta_error_min[g][i][j] = TMath::Sqrt(Integral_eta_min[g][i][j]*Integral_eta_min[g][i][j]+Integral_eta[g][i][j]*Integral_eta[g][i][j]*relative_error_percent*relative_error_percent);
	    
	    if(Spillover_and_pp_subtracted==kTRUE){Integral_eta[g][i][j]+= -1.*check_new_eta_rebin[g][i][j]->IntegralAndError(llimiteta,rlimiteta,pp_error,"width");	  }    
 
	       //if(Spillover_and_pp_subtracted==kFALSE){Integral_eta[g][i][j]+= -1.*check_new_eta_rebin[g][i][j]->IntegralAndError(llimiteta,rlimiteta,pp_error,"width");	  }    


	    Integral_eta_Pt[g][j]->SetBinContent(i+1,Integral_eta[g][i][j]);
	    Integral_eta_Pt[g][j]->SetBinError(i+1,Integral_eta_error[g][i][j]);
	   

	      
	 
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

	  
	    PbPb_pp_phi[g][i][j] = (TH1D*)check_new_phi_rebin[g-1][i][j]->Clone(PbPbppname_phi);
	    PbPb_pp_phi[g][i][j]->Add(check_new_phi_rebin[g][i][j],-1.);
	  
	 
	    
	    PbPb_pp_phi_syst[g][i][j] = (TH1D*)PbPb_pp_phi[g][i][j]->Clone(PbPbppname_phi_syst);
	      

	      
	    check_new_phi_max[g][i][j] = (TH1D*)check_new_phi_rebin[g-1][i][j]->Clone(PbPbppname_phi_max);
	    check_new_phi_max[g][i][j]->Add(Error_BG_Spillover_phi[g-1][i][j]);

	    check_new_phi_min[g][i][j] = (TH1D*)check_new_phi_rebin[g-1][i][j]->Clone(PbPbppname_phi_max);
	    check_new_phi_min[g][i][j]->Add(Error_BG_Spillover_phi[g-1][i][j],-1.);

	      
	    for(int k=1;k<check_new_phi_rebin[g][i][j]->GetNbinsX()+1;k++){
		
	      PbPb_err = check_new_phi_syst[g-1][i][j]->GetBinError(k);
	      pp_err = check_new_phi_syst[g][i][j]->GetBinError(k);
		
	      if(isnan(pp_err) || isinf(pp_err) ){pp_err=0.;}
	      if(isnan(PbPb_err) || isinf(PbPb_err) ){PbPb_err=0.;}

	      tot_err = TMath::Sqrt(PbPb_err*PbPb_err+pp_err*pp_err);

	      if(isnan(tot_err)){
		tot_err = 0.;

		cerr<<"I am so confused by not-a-number :( "<<endl;
		//return -1;
	      }
		
	      PbPb_pp_phi_syst[g][i][j]->SetBinError(k,tot_err);

	      //		cout<<g<<" "<<i<<" "<<j<<" "<<l<<k<<"   "<<PbPb_err<<" "<<pp_err<<" "<<tot_err<<endl;

	    }
	 



	    if(i==0){
	      TString IntegratedYieldname_phi = in_name;
	      IntegratedYieldname_phi.ReplaceAll("Result_pp","Integrated_Yield_Phi");
	      IntegratedYieldname_phi.ReplaceAll("Pt100_Pt300_","");
	      IntegratedYieldname_phi.ReplaceAll("_TrkPt1_TrkPt2","");
	      Integral_phi_Pt[g][j] = new TH1D(IntegratedYieldname_phi,"",4,pTbins);
	    }
	 
	    llimitphi = PbPb_pp_phi[g][i][j]->FindBin(-philim+.0001);
	    rlimitphi = PbPb_pp_phi[g][i][j]->FindBin(philim-.0001);
	 

	    Integral_phi[g][i][j] = check_new_phi_rebin[g-1][i][j]->IntegralAndError(llimitphi,rlimitphi,PbPb_error,"width");
	      
	    Integral_phi_error[g][i][j] = TMath::Sqrt(PbPb_error*PbPb_error+pp_error*pp_error);  
	      
	    for(int k=llimitphi;k<rlimitphi+1;k++){
	      if((check_new_phi_min[g][i][j]->GetBinContent(k)<0)&&(check_new_phi_rebin[g-1][i][j]->GetBinContent(k)>=0)){
		check_new_phi_min[g][i][j]->SetBinContent(k,0.);
	      }
	    }
	     	     	      
	    Integral_phi_max[g][i][j] = check_new_phi_max[g][i][j]->Integral(llimitphi,rlimitphi,"width")-Integral_phi[g][i][j];
	    
	    Integral_phi_min[g][i][j] = check_new_phi_min[g][i][j]->Integral(llimitphi,rlimitphi,"width")-Integral_phi[g][i][j];

	
	    Integral_phi_error_max[g][i][j] = TMath::Sqrt(Integral_phi_max[g][i][j]*Integral_phi_max[g][i][j]+Integral_phi[g][i][j]*Integral_phi[g][i][j]*relative_error_percent*relative_error_percent);
	    Integral_phi_error_min[g][i][j] = TMath::Sqrt(Integral_phi_min[g][i][j]*Integral_phi_min[g][i][j]+Integral_phi[g][i][j]*Integral_phi[g][i][j]*relative_error_percent*relative_error_percent);

	     if(Spillover_and_pp_subtracted == kTRUE){Integral_phi[g][i][j]+= -1.*check_new_phi_rebin[g][i][j]->IntegralAndError(llimitphi,rlimitphi,pp_error,"width");	  }    


	     //if(Spillover_and_pp_subtracted == kFALSE){Integral_phi[g][i][j]+= -1.*check_new_phi_rebin[g][i][j]->IntegralAndError(llimitphi,rlimitphi,pp_error,"width");	  }    



	    Integral_phi_Pt[g][j]->SetBinContent(i+1,Integral_phi[g][i][j]);
	    Integral_phi_Pt[g][j]->SetBinError(i+1,Integral_phi_error[g][i][j]);
	  	   
	    PbPb_pp_eta[g][i][j]->Write();
	    PbPb_pp_eta_syst[g][i][j]->Write();
	    PbPb_pp_phi[g][i][j]->Write();
	    PbPb_pp_phi_syst[g][i][j]->Write();
	    
	    
	  } //close Integral loop (for g==1, 3, 5)
	 
	} //close j  because otherwise we ignore all centralities after j=0!


      } // close i

  
      //-------------------------------------------

      //  Save all canvases

      //--------------------------------------------

	
      TString checksavenamePhi = in_name;
      checksavenamePhi.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenamePhi.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenamePhi.ReplaceAll("Result","Yield_Phi");
      checksavenamePhi +=".png";
      ccheckPhi[g]->SaveAs(checksavenamePhi);
   
	

      TString checksavenameEta = in_name;
      checksavenameEta.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenameEta.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenameEta.ReplaceAll("Result","Yield_Eta");
      checksavenameEta +=".png";
      ccheckEta[g]->SaveAs(checksavenameEta);


      TString checksavenameEta_wide = in_name;
      checksavenameEta_wide.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenameEta_wide.ReplaceAll("Cent0_Cent100_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      checksavenameEta_wide.ReplaceAll("Result","Yield_Eta_wide");
      checksavenameEta_wide +=".png";
      ccheckEta_wide[g]->SaveAs(checksavenameEta_wide);
      checksavenameEta_wide.ReplaceAll(".png",".pdf");
      ccheckEta_wide[g]->SaveAs(checksavenameEta_wide);

      //---------------------------------------------------
      //     Make Integral Histograms (but don't draw them)
      //--------------------------------------------------
	


      if(g==1||g==3||g==5){


	for(int j=0;j<4;j++){  //Start a new j-loop not in an i-loop -- we're plotting integrated yield by pT for each centrality class

	  //   Here make the asymmetric systematic errors:
	  
	  in_name =  make_name("Result_",g,3,j,0,centlabel,pTlabel);
	  

	  Integral_eta_value.clear();
	  Integral_eta_error_up.clear();
	  Integral_eta_error_down.clear();


	  for(int k = 0; k<4; k++){
	    bc = Integral_eta_Pt[g][j]->GetBinContent(k+1);
	    Integral_eta_value.push_back(bc);

	    err_up = Integral_eta_error_max[g][k][j];
	    Integral_eta_error_up.push_back(err_up);
					   
	    err_down = Integral_eta_error_min[g][k][j];
	    Integral_eta_error_down.push_back(err_down);


	     
	  }

	  // return 0; 
	  TString IntegratedYieldname_eta_syst = in_name;
	  IntegratedYieldname_eta_syst.ReplaceAll("Result_pp","Integrated_Yield_Eta_Syst");
	  IntegratedYieldname_eta_syst.ReplaceAll("Pt100_Pt300_","");
	  IntegratedYieldname_eta_syst.ReplaceAll("_TrkPt4_TrkPt8","");

	  Integral_eta_syst[g][j] = new TGraphAsymmErrors(pTbin_centers.size(),&pTbin_centers[0],&Integral_eta_value[0],&pTbin_errors[0],&pTbin_errors[0],&Integral_eta_error_down[0],&Integral_eta_error_up[0]);
	    
	  Integral_eta_syst[g][j]->SetName(IntegratedYieldname_eta_syst);

	 
	  //-------------
	  //  Now Phi
	  //--------------
	

	   	 

	  Integral_phi_value.clear();
	  Integral_phi_error_up.clear();
	  Integral_phi_error_down.clear();


	  for(int k = 0; k<4; k++){
	    bc = Integral_phi_Pt[g][j]->GetBinContent(k+1);
	    Integral_phi_value.push_back(bc);

	    err_up = Integral_phi_error_max[g][k][j];
	    Integral_phi_error_up.push_back(err_up);
					   
	    err_down = Integral_phi_error_min[g][k][j];
	    Integral_phi_error_down.push_back(err_down);

	    //	      cout<<Integral_phi_error_max[g][k][j]<<" "<<Integral_phi_error_up.at(k)<<endl;

	  }

	 
	  TString IntegratedYieldname_phi_syst = in_name;
	  IntegratedYieldname_phi_syst.ReplaceAll("Result_pp","Integrated_Yield_Phi_Syst");
	  IntegratedYieldname_phi_syst.ReplaceAll("Pt100_Pt300_","");
	  IntegratedYieldname_phi_syst.ReplaceAll("_TrkPt4_TrkPt8","");
		  	    
	  Integral_phi_syst[g][j] = new TGraphAsymmErrors(pTbin_centers.size(),&pTbin_centers[0],&Integral_phi_value[0],&pTbin_errors[0],&pTbin_errors[0],&Integral_phi_error_down[0],&Integral_phi_error_up[0]);

	  Integral_phi_syst[g][j]->SetName(IntegratedYieldname_phi_syst);
	  
	 

	  Integral_eta_Pt[g][j]->Write();
	  Integral_phi_Pt[g][j]->Write();
	  Integral_eta_syst[g][j]->Write();
	  Integral_phi_syst[g][j]->Write();

	   
	  //Auxillary TGraphs
   
	  Error_up_eta_pT[g][j] = new TGraph(pTbin_centers.size(),&pTbin_centers[0],&Integral_eta_error_up[0]);
	  
	  TString ErrorUpEtaPt_name = in_name;
	  ErrorUpEtaPt_name.ReplaceAll("Result_pp","Integral_Error_Up");
	  ErrorUpEtaPt_name.ReplaceAll("Pt100_Pt300_","");
	  ErrorUpEtaPt_name.ReplaceAll("_TrkPt4_TrkPt8","");
	 
	  Error_up_eta_pT[g][j]->SetName(ErrorUpEtaPt_name);
	  Error_up_eta_pT[g][j]->Write();
	    

	  Error_down_eta_pT[g][j] = new TGraph(pTbin_centers.size(),&pTbin_centers[0],&Integral_eta_error_down[0]);
	  TString ErrorDownEtaPt_name = ErrorUpEtaPt_name;
	  ErrorDownEtaPt_name.ReplaceAll("Up","Down");
	  Error_down_eta_pT[g][j]->SetName(ErrorDownEtaPt_name);
	  Error_down_eta_pT[g][j]->Write();
	    	    
	    
	}  //  end j loop for integral plotting	


      }  //end only do integrals for pp
         
   

  

  }// closes g loop


  //to cout errors
  for(int g = 0; g<6; g++){
    for(int i = 0; i<4; i++){
      for(int j = 0; j<4; j++){
	//	cout<<g<<" "<<i<<" "<<j<<" "<<error_A[g][i][j]<<" "<<error_B[g][i][j]<<" "<<error_C[g][i][j]<<" "<<error_D[g][i][j]<<" "<<error_E[g][i][j]<<" "<<error_F[g][i][j]<<" "<<syst_err[g][i][j]<<endl;

	cout<<g<<" "<<i<<" "<<j<<" "<<jff_residual_eta_combined[g][i][0]->Integral("width")<<" "<<check_new_eta_rebin[g][i][j]->Integral("width")<<" "<<jff_residual_eta_combined[g][i][0]->Integral()/check_new_eta_rebin[g][i][j]->Integral()<<" "<<jff_residual_eta_combined[g][i][0]->Integral()/check_new_eta_rebin[g][i][j]->Integral()*jff_percent[g][i]<<endl;

      }
    }
 

  }
 
  return 0;


}  //and we're done.
 
