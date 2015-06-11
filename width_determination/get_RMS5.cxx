
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
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TExec.h"
#include "TLatex.h"

#include "../HIN_14_016_functions.h"

#include <iostream>
#include <vector>
#include <fstream>


using namespace std;

Int_t get_RMS5(){

  
#include "../HIN_14_016_universals.h"

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.15);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);
  gStyle->SetTextFont(43);


  
  TFile *fin[12];
  TFile *fin_raw[12];
  TFile *fin_ref[12];
  TFile *fmc[6];
  Double_t xAxis[5] = {-100,-50,-30,-10,0}; 
  
  TH1D* int_cent[12][4];
  TH1D* blank[4];
  TH1D* blank2[4];
  TH1D* blank3[4];
  
  

  TH1D *check_new_phi_rebin[12][5][4];
  TH1D *check_new_phi_syst[12][5][4];
  
  TH1D *check_new_eta_rebin[12][5][4];
  TH1D *check_new_eta_syst[12][5][4];

  TH2D *raw_corr[12][5][4];
  TH1D *raw_corr_proj[12][5][4];
  TH1D *raw_corr_rebin[12][5][4];
  
  TH1D *PbPb_pp_eta[12][5][4];
  TH1D *PbPb_pp_eta_syst[12][5][4];

  TH1D *check_new_eta_up[12][5][4];
  TH1D *check_new_eta_down[12][5][4];


  TH1D *PbPb_pp_eta_up[12][5][4];
  TH1D *PbPb_pp_eta_down[12][5][4];



  TH1D *check_new_phi_up[12][5][4];
  TH1D *check_new_phi_down[12][5][4];

  TH1D *PbPb_pp_phi_up[12][5][4];
  TH1D *PbPb_pp_phi_down[12][5][4];



  float ref_bc, ref_err, new_bc, new_err, RMS_val_eta_temp, RMS_err_eta_temp, RMS_val_phi_temp, RMS_err_phi_temp, RMS_diff_val_eta_temp, RMS_diff_err_eta_temp, RMS_diff_val_phi_temp, RMS_diff_err_phi_temp, val_l,val_r,err_l,err_r, spill;

  TGraphAsymmErrors *RMS_eta[12][4];
  TGraphAsymmErrors *RMS_diff_eta[12][4];
  TGraphAsymmErrors *RMS_phi[12][4];
  TGraphAsymmErrors *RMS_diff_phi[12][4];

  TGraphAsymmErrors *RMS_eta_demo[12][4];
  
  TGraphAsymmErrors *RMS_eta_cent[12][4];
  TGraphAsymmErrors *RMS_diff_eta_cent[12][4];
  
  TGraphAsymmErrors *RMS_eta_ref[12][4];
  TGraphAsymmErrors *RMS_phi_ref[12][4];
  TGraphAsymmErrors *RMS_eta_cent_ref[12][4];


  TGraphAsymmErrors *RMS_PbPb_minus_pp_eta[12][4];
  TGraphAsymmErrors *RMS_PbPb_minus_pp_phi[12][4];

  TCanvas *RMS_eta_canvas[12];
  TCanvas *RMS_phi_canvas[12];
  TCanvas *cRMS_eta_cent[12];
  TCanvas *check_fits_diff[12]; 


  vector<float> RMS_val_eta;
  vector<float> RMS_err_eta_up;
  vector<float> RMS_err_eta_down;


  vector<float> RMS_diff_val_eta;
  vector<float> RMS_diff_err_eta_up;
  vector<float> RMS_diff_err_eta_down;

 
  vector<float> RMS_val_phi;
  vector<float> RMS_err_phi_up;
  vector<float> RMS_err_phi_down;


  vector<float> RMS_diff_val_phi;
  vector<float> RMS_diff_err_phi_up;
  vector<float> RMS_diff_err_phi_down;

   
  vector<float> RMS_val_cent_eta;
  vector<float> RMS_err_cent_eta_up;
  vector<float> RMS_err_cent_eta_down;

   
  vector<float> RMS_val_cent_phi;
  vector<float> RMS_err_cent_phi_up;
  vector<float> RMS_err_cent_phi_down;

    
  vector<float> RMS_diff_val_cent_eta;
  vector<float> RMS_diff_err_cent_eta_up;
  vector<float> RMS_diff_err_cent_eta_down;

   
  vector<float> RMS_diff_val_cent_phi;
  vector<float> RMS_diff_err_cent_phi_up;
  vector<float> RMS_diff_err_cent_phi_down;

 
  
  vector<float> RMS_eta_demo_val_0;  
  vector<float> RMS_eta_demo_val_1;  
  vector<float> RMS_eta_demo_val_2;  
  vector<float> RMS_eta_demo_val_3;  

  
  float evalpt, value, exh,exl,eyh,eyl;
  


  TLegend *l40,*l41,*l42;

  TString datalabel, in_name, centlabel, pTlabel;
  TString   leadingrefetaname[5][4];
  float raw_min,raw_max,mixed_min,mixed_max,yield_min,yield_max,result_min,result_max, temp_int, temp_sigma, full_int;
  int counter, lbin, rbin;

  TLine *linePhi,*lineEta,*linePt, *lineCent;

  TLatex *tex24phi,*tex24eta;

  vector<float> pTbin_centers;
  pTbin_centers.push_back(1.5);
  pTbin_centers.push_back(2.5);
  pTbin_centers.push_back(3.5);
  pTbin_centers.push_back(6.0);
  vector<float> pTbin_errors;
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(.5);
  pTbin_errors.push_back(2.);


  vector<float> cent_centers;
  cent_centers.push_back(-75);
  cent_centers.push_back(-40);
  cent_centers.push_back(-20);
  cent_centers.push_back(-5);
  vector<float> cent_errors;
  cent_errors.push_back(25);
  cent_errors.push_back(10);
  cent_errors.push_back(10);
  cent_errors.push_back(5);


  TH1D *HYDJET_PYTHIA_eta[6][5][4];
  TH1D *HYDJET_PYTHIA_phi[6][5][4];
 


  TFile *fout[6];
  
  TCanvas *check_fits_eta[6];
  TCanvas *check_fits_phi[6];
 
  
  float par0, par1, par2, par3, par4, special_err, error_up, error_down,   val_2, err_syst, dx_eta;

  //-------------------------
  //   PARAMETER SET

  float range = 1.5;
  int i_max = 2;
  //----------------------------



  TF1 *double_gaus = new TF1("double_gaus","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))+[2]*TMath::Exp(-pow(TMath::Abs(x)/[3],2))+[4]",-range,range);

  TF1 *gaus_1[12][5][4];
  TF1 *gaus_2[12][5][4];
  TF1 *gaus_tot[12][5][4];

  TF1 *flat_line = new TF1("flat_line","-[0]+x-x",-3.5,3.5);


  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------

  for(int g = 0; g<6; g++){
      
    cout<<"starting getting files for g = "<<g<<endl;
    //Open files
    
    
    
    switch(g){
    case 0:
      //   fin[g] = new TFile("../analysis/Inclusive_Data_NoSpillOver_AllPlots.root", "READ");
      fin[g] = new TFile("../analysis/Inclusive_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/final_plots/RMS_Output_Inclusive.root", "READ");
      fin_raw[g] = new TFile("../me_correct/PbPb_Inclusive_Correlations.root", "READ");
      fmc[g] = new TFile("../spill_over/Inclusive_SpillOvers.root","READ");
      datalabel = "Inclusive";     break;
    case 1:
      //   fin[g] = new TFile("../analysis/Inclusive_Data_NoSpillOver_AllPlots.root", "READ");
      fin[g] = new TFile("../analysis/Inclusive_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/final_plots/RMS_Output_Inclusive.root", "READ");
      fin_raw[g] = new TFile("../me_correct/pp_Inclusive_Correlations.root", "READ");
      fmc[g] = new TFile("../spill_over/Inclusive_SpillOvers.root","READ");
      //   fmc[g] = new TFile("../analysis/Inclusive_Closures.root","READ");
      datalabel = "Inclusive";     break;
    case 2:
      fin[g] = new TFile("../analysis/SubLeading_Data_AllPlots.root", "READ");
      //  fin[g] = new TFile("../analysis/SubLeading_Data_NoSpillOver_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/final_plots/RMS_Output_SubLeading.root", "READ");
      fin_raw[g] = new TFile("../me_correct/PbPb_SubLeading_Correlations.root", "READ");
      fmc[g] = new TFile("../spill_over/SubLeading_SpillOvers.root","READ");
      datalabel = "SubLeading";    break;
    case 3:
      fin[g] = new TFile("../analysis/SubLeading_Data_AllPlots.root", "READ");
      // fin[g] = new TFile("../analysis/SubLeading_Data_NoSpillOver_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/final_plots/RMS_Output_SubLeading.root", "READ");
      fin_raw[g] = new TFile("../me_correct/pp_SubLeading_Correlations.root", "READ");
      //    fmc[g] = new TFile("../analysis/SubLeading_Closures.root","READ");
      fmc[g] = new TFile("../spill_over/SubLeading_SpillOvers.root","READ");
      datalabel = "SubLeading";    break;
    case 4:
      fin[g] = new TFile("../analysis/Leading_Data_AllPlots.root", "READ");
      fin_raw[g] = new TFile("../me_correct/PbPb_Leading_Correlations.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/final_plots/RMS_Output_Leading.root", "READ");
      fmc[g] = new TFile("../spill_over/Leading_SpillOvers.root", "READ");      
      //  fin[g] = new TFile("../analysis/Leading_Data_NoSpillOver_AllPlots.root", "READ");
      datalabel = "Leading";       break;
    case 5:
      //  fin[g] = new TFile("../analysis/Leading_Data_NoSpillOver_AllPlots.root", "READ");
      fin[g] = new TFile("../analysis/Leading_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/final_plots/RMS_Output_Leading.root", "READ");
      fin_raw[g] = new TFile("../me_correct/pp_Leading_Correlations.root", "READ");
      fmc[g] = new TFile("../spill_over/Leading_SpillOvers.root", "READ");      
      //   fmc[g] = new TFile("../analysis/Leading_Closures.root", "READ");      
      datalabel = "Leading";       break;
    }
    

    if(g==1||g==3||g==5){
      fout[g] = new TFile((TString)("RMS_Output_"+datalabel+".root"),"RECREATE");
    }

    cout<<"got files"<<endl;

    TString RMS_eta_canvas_name = "RMS_eta";  RMS_eta_canvas_name+=g;
    RMS_eta_canvas[g] = new TCanvas(RMS_eta_canvas_name,"",10,10,1500,850);
    RMS_eta_canvas[g]->Divide(4,2,0.,0.);


    TString RMS_phi_canvas_name = "RMS_phi";  RMS_phi_canvas_name+=g;
    RMS_phi_canvas[g] = new TCanvas(RMS_phi_canvas_name,"",10,10,1500,850);
    RMS_phi_canvas[g]->Divide(4,2,0.,0.);


    TString cRMS_eta_cent_name = "RMS_eta_cent";  cRMS_eta_cent_name+=g;
    if(g==5){  cRMS_eta_cent[g] = new TCanvas(cRMS_eta_cent_name,"",10,10,1500,800);
      cRMS_eta_cent[g]->Divide(2,2,0.,0.);
    }else{
      cRMS_eta_cent[g] = new TCanvas(cRMS_eta_cent_name,"",10,10,1500,400);
      cRMS_eta_cent[g]->Divide(2,1,0.,0.);
    }

    TString checkfitname_eta = "CheckFitEta";
    checkfitname_eta+=g;
    check_fits_eta[g] = new TCanvas(checkfitname_eta," ",10,10,1500,1600);
    check_fits_eta[g]->Divide(4,4,0.,0.);
    
    TString checkfitname_phi = "CheckFitPhi";
    checkfitname_phi+=g;
    check_fits_phi[g] = new TCanvas(checkfitname_phi," ",10,10,1500,1600);
    check_fits_phi[g]->Divide(4,4,0.,0.);
    


    TString checkfitdiffname = "CheckFitDiff";
    checkfitdiffname+=g;
    check_fits_diff[g] = new TCanvas(checkfitdiffname," ",10,10,1500,800);
    check_fits_diff[g]->Divide(4,2,0.,0.);

    cout<<"did canvases"<<endl;

    

    TCanvas *dummy = new TCanvas("dummy");
  
    for(int j = 0; j<4; j++){

      RMS_val_eta.clear();
      RMS_err_eta_up.clear();
      RMS_err_eta_down.clear();
 

      RMS_diff_val_eta.clear();
      RMS_diff_err_eta_up.clear();
      RMS_diff_err_eta_down.clear();
 
      RMS_val_phi.clear();
      RMS_err_phi_up.clear();
      RMS_err_phi_down.clear();
 

      RMS_diff_val_phi.clear();
      RMS_diff_err_phi_up.clear();
      RMS_diff_err_phi_down.clear();
 



      for(int i=0; i<4; i++){

	dummy->cd();


	in_name = make_name("Yield_",g,i,j,0,centlabel,pTlabel);


	raw_corr[g][i][j] = (TH2D*)fin_raw[g]->Get(in_name)->Clone(in_name);

	TString proj_name = in_name; 
	proj_name.ReplaceAll("Yield","Raw_Eta_Proj");

	lbin = raw_corr[g][i][j]->GetYaxis()->FindBin(-1.0+.00001);
	rbin = raw_corr[g][i][j]->GetYaxis()->FindBin(1.0-.00001);

	raw_corr_proj[g][i][j] = (TH1D*)raw_corr[g][i][j]->ProjectionX(proj_name,lbin,rbin);

	dx_eta = raw_corr_proj[g][i][j]->GetBinWidth(1);

	raw_corr_proj[g][i][j]->Scale(1./dx_eta);
	
	//	raw_corr_rebin[g][i][j] = (TH1D*)Rebin_dEta(raw_corr_proj[g][i][j]);

	//-------------------------------------
	//  Get and assign all histograms 
	//-------------------------------------

	TString newchecknameEta_rebin = in_name;
	  	  
	newchecknameEta_rebin.ReplaceAll("Yield", "New_check_Eta");
	newchecknameEta_rebin += "rebin";
	//	check_new_eta_rebin[g][i][j] = (TH1D*)fin[g]->Get(newchecknameEta_rebin)->Clone(newchecknameEta_rebin);

	check_new_eta_rebin[g][i][j] = (TH1D*)raw_corr_proj[g][i][j]->Clone(newchecknameEta_rebin);
	

	TString newchecknamePhi_rebin = in_name;
	newchecknamePhi_rebin.ReplaceAll("Yield", "New_check_Phi");
	newchecknamePhi_rebin += "rebin";
	check_new_phi_rebin[g][i][j] = (TH1D*)fin[g]->Get(newchecknamePhi_rebin)->Clone(newchecknamePhi_rebin);

       	cout<<"got histos"<<endl;
	
	if(g==0||g==2||g==4){
 
	  TString HminusPnameEta = in_name;
	  HminusPnameEta.ReplaceAll("Pt100_Pt300_","");
	  HminusPnameEta.ReplaceAll("Yield_PbPb","SpillOvers_Eta_RecoJet_GenTrack");
	 
	  TString HminusPnamePhi = in_name;
	  HminusPnamePhi.ReplaceAll("Pt100_Pt300_","");
	  HminusPnamePhi.ReplaceAll("Yield_PbPb","SpillOvers_Phi_RecoJet_GenTrack");

	  if(g==2){

	    HminusPnameEta.ReplaceAll("RecoJet","SubLeading_RecoJet");
	    HminusPnamePhi.ReplaceAll("RecoJet","SubLeading_RecoJet");
	  }else if(g==4){

	    HminusPnameEta.ReplaceAll("RecoJet","Leading_RecoJet");
	    HminusPnamePhi.ReplaceAll("RecoJet","Leading_RecoJet");
	  }

	  cout<<HminusPnameEta<<endl;

	  HYDJET_PYTHIA_eta[g][i][j] = (TH1D*)fmc[g]->Get(HminusPnameEta)->Clone(HminusPnameEta);
	  HYDJET_PYTHIA_phi[g][i][j] = (TH1D*)fmc[g]->Get(HminusPnamePhi)->Clone(HminusPnamePhi);


	  TString newchecknameEtasyst = newchecknameEta_rebin;
	  newchecknameEtasyst.ReplaceAll("New_check","Syst_error");
	  newchecknameEtasyst.ReplaceAll("rebin","");
	  check_new_eta_syst[g][i][j] = (TH1D*)fin[g]->Get(newchecknameEtasyst)->Clone(newchecknameEtasyst);
	

	  TString newchecknamePhisyst = newchecknamePhi_rebin;
	  newchecknamePhisyst.ReplaceAll("New_check","Syst_error");
	  newchecknamePhisyst.ReplaceAll("rebin","");
	  check_new_phi_syst[g][i][j] = (TH1D*)fin[g]->Get(newchecknamePhisyst)->Clone(newchecknamePhisyst);



	  newchecknameEta_rebin.ReplaceAll("rebin","up");
	  check_new_eta_up[g][i][j] = (TH1D*)check_new_eta_rebin[g][i][j]->Clone(newchecknameEta_rebin);
	
	  newchecknameEta_rebin.ReplaceAll("up","down");
	  check_new_eta_down[g][i][j] = (TH1D*)check_new_eta_rebin[g][i][j]->Clone(newchecknameEta_rebin);

	  nbins = check_new_eta_rebin[g][i][j]->GetNbinsX()+1;
	
	  for(int k = 1; k<nbins/2+1;k++){
	    
	    val_l = check_new_eta_rebin[g][i][j]->GetBinContent(k);
	    err_l  = check_new_eta_rebin[g][i][j]->GetBinError(k);
	    
	    val_r = check_new_eta_rebin[g][i][j]->GetBinContent(nbins-k);
	    err_r = check_new_eta_rebin[g][i][j]->GetBinError(nbins-k);
	    
	    cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<TMath::Sqrt(err_l*err_l+err_r*err_r)/2<<endl;

	    check_new_eta_rebin[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_eta_rebin[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	    check_new_eta_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_eta_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	  
	    if(k<nbins/2){
	      check_new_eta_rebin[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      check_new_eta_rebin[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    }
	    if(g==0||g==2||g==4){	    
	      val_2 = check_new_eta_rebin[g][i][j]->GetBinContent(k);

	      err_syst  = check_new_eta_syst[g][i][j]->GetBinError(check_new_eta_syst[g][i][j]->FindBin(check_new_eta_rebin[g][i][j]->GetBinCenter(k)))-check_new_eta_syst[g][i][j]->GetBinError(1);
	      //err_syst  = check_new_eta_syst[g][i][j]->GetBinError(k);
	      //  err_syst  = HYDJET_PYTHIA_eta[g][i][j]->GetBinContent(k)/2;

	      spill = HYDJET_PYTHIA_eta[g][i][j]->GetBinContent(HYDJET_PYTHIA_eta[g][i][j]->FindBin(check_new_eta_rebin[g][i][j]->GetBinCenter(k)));
	      
	      check_new_eta_rebin[g][i][j]->SetBinContent(k,val_2-spill);
	      check_new_eta_rebin[g][i][j]->SetBinContent(nbins-k,val_2-spill);

	      check_new_eta_up[g][i][j]->SetBinContent(k,val_2+err_syst-spill);
	      check_new_eta_up[g][i][j]->SetBinContent(nbins-k,val_2+err_syst-spill);
	      check_new_eta_down[g][i][j]->SetBinContent(k,val_2-err_syst-spill);
	      check_new_eta_down[g][i][j]->SetBinContent(nbins-k,val_2-err_syst-spill);
	   
	      if(k<nbins/2){
		check_new_eta_down[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_eta_down[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_eta_up[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_eta_up[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);	 
	      }

	    }
	  }

	  //	  return 0;

	  newchecknamePhi_rebin.ReplaceAll("rebin","up");
	  check_new_phi_up[g][i][j] = (TH1D*)check_new_phi_rebin[g][i][j]->Clone(newchecknamePhi_rebin);
	
	  newchecknamePhi_rebin.ReplaceAll("up","down");
	  check_new_phi_down[g][i][j] = (TH1D*)check_new_phi_rebin[g][i][j]->Clone(newchecknamePhi_rebin);


	  
	  nbins = check_new_phi_rebin[g][i][j]->GetNbinsX()+1;
	  
	  for(int k = 1; k<nbins/2+1;k++){
	    val_l = check_new_phi_rebin[g][i][j]->GetBinContent(k);
	 
	    err_l  = check_new_phi_rebin[g][i][j]->GetBinError(k);
	    
	    val_r = check_new_phi_rebin[g][i][j]->GetBinContent(nbins-k);
	    err_r = check_new_phi_rebin[g][i][j]->GetBinError(nbins-k);
	    
	    cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<endl;

	    check_new_phi_rebin[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_phi_rebin[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	    check_new_phi_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_phi_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	   
	    if(k<nbins/2){
	      check_new_phi_rebin[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      check_new_phi_rebin[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    }


	    if(g==0||g==2||g==4){
	      val_2 = check_new_phi_rebin[g][i][j]->GetBinContent(k);
	      err_syst  = check_new_phi_syst[g][i][j]->GetBinError(k);
	      //err_syst  = HYDJET_PYTHIA_phi[g][i][j]->GetBinContent(k)/2;
	    
	      check_new_phi_up[g][i][j]->SetBinContent(k,val_2+err_syst);
	      check_new_phi_up[g][i][j]->SetBinContent(nbins-k,val_2+err_syst);
	      check_new_phi_down[g][i][j]->SetBinContent(k,val_2-err_syst);
	      check_new_phi_down[g][i][j]->SetBinContent(nbins-k,val_2-err_syst);

	      if(k<nbins/2){
		check_new_phi_up[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_phi_up[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_phi_down[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
		check_new_phi_down[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      }
	    }
	
	  }
	}

	
	//Difference histos

	if(g==1||g==3||g==5){
	  TString Difference_name = in_name;

	  /*	  	  
	  Difference_name.ReplaceAll("Yield_pp", "PbPb_minus_pp_eta");
	  PbPb_pp_eta[g][i][j] = (TH1D*)fin[g]->Get( Difference_name)->Clone( Difference_name);

	  Difference_name.ReplaceAll("eta","eta_syst");
	  PbPb_pp_eta_syst[g][i][j] = (TH1D*)fin[g]->Get(Difference_name)->Clone(Difference_name);

	  Difference_name += "up";
	  PbPb_pp_eta_up[g][i][j] = (TH1D*)PbPb_pp_eta[g][i][j]->Clone(	Difference_name);
	
	  Difference_name.ReplaceAll("up","down");
	  PbPb_pp_eta_down[g][i][j] = (TH1D*)PbPb_pp_eta[g][i][j]->Clone(Difference_name);
	  */
	 
	  cout<<"starting"<<endl;

	  Difference_name.ReplaceAll("Yield_pp", "PbPb_minus_pp_eta");
	  PbPb_pp_eta[g][i][j] = (TH1D*)check_new_eta_rebin[g-1][i][j]->Clone(Difference_name);

	  Difference_name += "up";
	  PbPb_pp_eta_up[g][i][j] = (TH1D*)check_new_eta_up[g-1][i][j]->Clone(Difference_name);
		
	  Difference_name.ReplaceAll("up","down");
	  PbPb_pp_eta_down[g][i][j] = (TH1D*)check_new_eta_down[g-1][i][j]->Clone(Difference_name);


	  PbPb_pp_eta[g][i][j]->Add(check_new_eta_rebin[g][i][j],-1.);
	  PbPb_pp_eta_up[g][i][j]->Add(check_new_eta_rebin[g][i][j],-1.);
	  PbPb_pp_eta_down[g][i][j]->Add(check_new_eta_rebin[g][i][j],-1.);

	  cout<<"done"<<endl;

	  /*
	  nbins = PbPb_pp_eta[g][i][j]->GetNbinsX()+1;
	  
	  for(int k = 1; k<nbins/2+1;k++){
	    val_l = PbPb_pp_eta[g][i][j]->GetBinContent(k);
	 
	    err_l  = PbPb_pp_eta[g][i][j]->GetBinError(k);
	    
	    val_r = PbPb_pp_eta[g][i][j]->GetBinContent(nbins-k);
	    err_r = PbPb_pp_eta[g][i][j]->GetBinError(nbins-k);
	    
	    cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<TMath::Sqrt(err_l*err_l+err_r*err_r)/2<<endl;

	    PbPb_pp_eta[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    PbPb_pp_eta[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	    PbPb_pp_eta_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    PbPb_pp_eta_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);

	    if(k<nbins/2){
	      PbPb_pp_eta[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      PbPb_pp_eta[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    }
	    
	    val_2 = PbPb_pp_eta[g][i][j]->GetBinContent(k);
	    err_syst  = PbPb_pp_eta_syst[g][i][j]->GetBinError(k)-PbPb_pp_eta_syst[g][i][j]->GetBinError(1);
	    //   err_syst  = PbPb_pp_eta_syst[g][i][j]->GetBinError(k);
	    //err_syst  = HYDJET_PYTHIA_eta[g-1][i][j]->GetBinContent(k)/2;
	    
	    PbPb_pp_eta_up[g][i][j]->SetBinContent(k,val_2+err_syst);
	    PbPb_pp_eta_up[g][i][j]->SetBinContent(nbins-k,val_2+err_syst);
	    PbPb_pp_eta_down[g][i][j]->SetBinContent(k,val_2-err_syst);
	    PbPb_pp_eta_down[g][i][j]->SetBinContent(nbins-k,val_2-err_syst);

	    if(k<nbins/2){
	      PbPb_pp_eta_up[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      PbPb_pp_eta_up[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      PbPb_pp_eta_down[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      PbPb_pp_eta_down[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    }
	  }
	  */
	
	}

	//--------------------------------------
	//    Get RMS and error
	//--------------------------------------

	cout<<"*********************************************"<<endl;

	cout<<in_name<<endl;

	cout<<g<<" "<<i<<" "<<j<<" "<<endl;
	
	double_gaus->ReleaseParameter(0);
	double_gaus->ReleaseParameter(2);
	double_gaus->ReleaseParameter(1);
	double_gaus->ReleaseParameter(3);
	double_gaus->ReleaseParameter(4);

	//	double_gaus->FixParameter(4,0.);

	double_gaus->SetParLimits(0,0.0,10.);  
	double_gaus->SetParLimits(1,0.2,range);
	double_gaus->SetParLimits(2,0.0,10.);  
	

	double_gaus->SetParameter(1,0.3);
	double_gaus->SetParLimits(1,0.1,0.67*range);
	double_gaus->SetParameter(3,0.6);
	double_gaus->SetParLimits(3,0.2,0.67*range);
	double_gaus->SetParameter(0,1.);
	double_gaus->SetParameter(2,1.);

	check_new_eta_rebin[g][i][j]->Fit(double_gaus,"","",-range-1.,range+1.);

	flat_line->SetParameter(0,double_gaus->GetParameter(4));

	check_new_eta_rebin[g][i][j]->Add(flat_line);
	if(g%2==0){
	  check_new_eta_up[g][i][j]->Add(flat_line);
	  check_new_eta_down[g][i][j]->Add(flat_line);
	}

	check_new_eta_rebin[g][i][j]->Fit(double_gaus,"","",-range-1.,range+1.);

	
	gaus_1[g][i][j] = new TF1("gaus_1","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))",-range,range);
	gaus_2[g][i][j] = new TF1("gaus_1","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))",-range,range);
	gaus_tot[g][i][j] =  new TF1("double_gaus","[0]*TMath::Exp(-pow(TMath::Abs(x)/[1],2))+[2]*TMath::Exp(-pow(TMath::Abs(x)/[3],2))",-range,range);


	par0 = double_gaus->GetParameter(0);
	par1 = double_gaus->GetParameter(1);
	par2 = double_gaus->GetParameter(2);
	par3 = double_gaus->GetParameter(3);


	//	temp_sigma = TMath::Sqrt((par0*par1*par1+par2*par3*par3)/(par0+par2));

	gaus_1[g][i][j]->SetParameter(0,par0);
	gaus_1[g][i][j]->SetParameter(1,par1);
	gaus_2[g][i][j]->SetParameter(0,par2);
	gaus_2[g][i][j]->SetParameter(1,par3);

	gaus_tot[g][i][j]->SetParameter(0,par0);
	gaus_tot[g][i][j]->SetParameter(1,par1);
	gaus_tot[g][i][j]->SetParameter(2,par2);
	gaus_tot[g][i][j]->SetParameter(3,par3);


	full_int = gaus_tot[g][i][j]->Integral(0,range);
	temp_int = 0;
	temp_sigma = 0;

	counter = 0;
	while(temp_int <0.67*full_int){
	  temp_sigma = range*counter/50000;
	  temp_int = gaus_tot[g][i][j]->Integral(0, temp_sigma);
	  counter++;
	}

	RMS_val_eta.push_back(temp_sigma);

	if(g==0||g==2||g==4){
	  check_new_eta_up[g][i][j]->Fit(double_gaus,"","",-range-1.,range+1.);

	  double_gaus->SetParameter(4,0.);

	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	
	  error_up = temp_sigma-RMS_val_eta.at(i);
	  /*
	    if(g==4&&i==3&&j==3){
	    error_up = special_err;
	    }else{
	    error_up = temp_sigma-RMS_val_eta.at(i);
	    }
	  */
	 
	
	  check_new_eta_down[g][i][j]->Fit(double_gaus,"","",-range-1.,range+1.);

	  double_gaus->SetParameter(4,0.);

	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	  
		  
	  error_down = temp_sigma-RMS_val_eta.at(i);

	  //	  error_down = error_up;
	
	  if(error_up>0&&error_down<0){
	    RMS_err_eta_up.push_back(TMath::Abs(error_up));
	    RMS_err_eta_down.push_back(TMath::Abs(error_down));
	  }else if(error_up<0&&error_down>0){
	    RMS_err_eta_down.push_back(TMath::Abs(error_up));
	    RMS_err_eta_up.push_back(TMath::Abs(error_down));
	  } else {
	    RMS_err_eta_up.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    RMS_err_eta_down.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    cerr<<"Both errors go same way "<<error_up<<" "<<error_down<<endl;
	    //return -1;
	  }

	  
	  

	} //close error for PbPb only
      

	///------------------
	//   Phi

	//-----------

	double_gaus->FixParameter(4,0.);

	check_new_phi_rebin[g][i][j]->Fit(double_gaus,"","",-range,range);

	full_int = double_gaus->Integral(0,range);
	temp_int = 0;
	temp_sigma = 0;

	counter = 0;
	while(temp_int <0.67*full_int){
	  temp_sigma = range*counter/50000;
	  temp_int = double_gaus->Integral(0, temp_sigma);
	  counter++;
	}

	RMS_val_phi.push_back(temp_sigma);
	

	if(g==0||g==2||g==4){
	  check_new_phi_up[g][i][j]->Fit(double_gaus,"","",-range,range);

	
	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	
	  error_up = temp_sigma-RMS_val_phi.at(i);

	  check_new_phi_down[g][i][j]->Fit(double_gaus,"","",-range,range);

	
	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }

	  //  RMS_err_phi_down.push_back(TMath::Abs(temp_sigma-RMS_val_phi.at(i)));
	

	  error_down = temp_sigma-RMS_val_phi.at(i);

	  if(error_up>0&&error_down<0){
	    RMS_err_phi_up.push_back(TMath::Abs(error_up));
	    RMS_err_phi_down.push_back(TMath::Abs(error_down));
	  }else if(error_up<0&&error_down>0){
	    RMS_err_phi_down.push_back(TMath::Abs(error_up));
	    RMS_err_phi_up.push_back(TMath::Abs(error_down));
	  } else {
	    cerr<<"Both errors go same way (dphi) "<<" "<<error_up<<" "<<error_down<<endl;
	    RMS_err_phi_up.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    RMS_err_phi_down.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    //  return -1;
	  }
	  


	} //close error for PbPb only
	

      }// close i

    
      RMS_eta[g][j] = new TGraphAsymmErrors(pTbin_centers.size(),&pTbin_centers[0],&RMS_val_eta[0],&pTbin_errors[0],&pTbin_errors[0],&RMS_err_eta_down[0],&RMS_err_eta_up[0]);
      TString RMS_eta_name = "RMS_eta"; RMS_eta_name+= g; RMS_eta_name+= j;
      RMS_eta[g][j]->SetName(RMS_eta_name);

      RMS_eta_ref[g][j] = (TGraphAsymmErrors*)fin_ref[g]->Get(RMS_eta_name)->Clone((TString)(RMS_eta_name+="_ref"));
      

      RMS_phi[g][j] = new TGraphAsymmErrors(pTbin_centers.size(),&pTbin_centers[0],&RMS_val_phi[0],&pTbin_errors[0],&pTbin_errors[0],&RMS_err_phi_down[0],&RMS_err_phi_up[0]);
      TString RMS_phi_name = "RMS_phi"; RMS_phi_name+= g; RMS_phi_name+= j;
      RMS_phi[g][j]->SetName(RMS_phi_name);
    
      RMS_phi_ref[g][j] = (TGraphAsymmErrors*)fin_ref[g]->Get(RMS_phi_name)->Clone((TString)(RMS_phi_name+="_ref"));

      RMS_eta_ref[g][j]->SetMarkerColor(kRed);      
      RMS_eta_ref[g][j]->SetLineColor(kRed);      

      RMS_phi_ref[g][j]->SetMarkerColor(kRed);      
      RMS_phi_ref[g][j]->SetLineColor(kRed);      


      if(g==1||g==3||g==5){
	TString RMS_PbPb_minus_pp_eta_name = "RMS_PbPb_minus_pp_eta"; RMS_PbPb_minus_pp_eta_name+= g; RMS_PbPb_minus_pp_eta_name+= j;
	RMS_PbPb_minus_pp_eta[g][j]  = (TGraphAsymmErrors*) RMS_eta[g-1][j]->Clone(RMS_PbPb_minus_pp_eta_name);
    
	for(int k = 0; k<4; k++){
	  
	  float temp_x = pTbin_centers.at(k);
	  
	  float PbPb_val = RMS_eta[g-1][j]->Eval(temp_x);
	  float pp_val = RMS_eta[g][j]->Eval(temp_x);

	  RMS_PbPb_minus_pp_eta[g][j]->SetPoint(k,temp_x,PbPb_val-pp_val);

	}

     
	RMS_eta_canvas[g]->cd(j+1);
	RMS_eta[g-1][j]->Draw("AP");

	RMS_eta[g-1][j]->GetXaxis()->SetRangeUser(.8,8.2);


	switch(g){
	case 1: 
	  RMS_eta[g-1][j]->SetMarkerStyle(10);
	  RMS_eta[g][j]->SetMarkerStyle(4);
	  RMS_eta_ref[g-1][j]->SetMarkerStyle(10);
	  RMS_eta_ref[g][j]->SetMarkerStyle(4);
	  RMS_eta[g-1][j]->SetFillColor(90);
	  RMS_eta[g-1][j]->SetLineColor(90);

	  break;
	case 3: 
	  RMS_eta[g-1][j]->SetMarkerStyle(34);
	  RMS_eta[g-1][j]->SetMarkerSize(2);
	  RMS_eta_ref[g-1][j]->SetMarkerStyle(34);
	  RMS_eta_ref[g-1][j]->SetMarkerSize(2);
	  RMS_eta[g-1][j]->SetFillColor(30);
	  RMS_eta[g-1][j]->SetLineColor(30);
	  RMS_eta[g][j]->SetMarkerStyle(28);
	  break;  
	case 5: 
	  RMS_eta[g-1][j]->SetMarkerStyle(21);
	  RMS_eta_ref[g-1][j]->SetMarkerStyle(21);
	  RMS_eta[g-1][j]->SetFillColor(kOrange-2);
	  RMS_eta[g-1][j]->SetLineColor(kOrange-2);
	  RMS_eta[g][j]->SetMarkerStyle(25);
	  break;
	}

	RMS_eta[g-1][j]->SetTitle("");

	RMS_eta[g][j]->SetMarkerSize(2);
	RMS_eta[g-1][j]->SetMarkerSize(2);

	RMS_eta_ref[g-1][j]->SetTitle("");

	RMS_eta_ref[g][j]->SetMarkerSize(2);
	RMS_eta_ref[g-1][j]->SetMarkerSize(2);

	if(j==0){
	  RMS_eta[g-1][j]->GetYaxis()->SetTitle("#Delta#eta Width");
	  RMS_eta[g-1][j]->GetYaxis()->SetTitleSize(tstitle);
	  RMS_eta[g-1][j]->GetYaxis()->SetLabelSize(ts);
	  RMS_eta[g-1][j]->GetYaxis()->CenterTitle();

	}else{
	  RMS_eta[g-1][j]->GetYaxis()->SetLabelSize(0);
	}

	RMS_eta[g-1][j]->SetMarkerColor(kBlack);
	RMS_eta[g-1][j]->SetMinimum(-0.02);
	RMS_eta[g-1][j]->SetMaximum(0.75);

	RMS_eta[g-1][j]->SetLineWidth(0);


	RMS_eta[g-1][j]->GetXaxis()->SetNdivisions(8);
	RMS_eta[g-1][j]->GetYaxis()->SetNdivisions(306);

	RMS_eta[g-1][j]->Draw("p Z e2");

	RMS_eta[g][j]->SetMarkerSize(2);

	RMS_eta[g-1][j]->Draw("p Z same");
	RMS_eta[g][j]->Draw("p X same");

	RMS_eta_ref[g-1][j]->Draw("same P");
	RMS_eta_ref[g][j]->Draw("same P X");
	//	RMS_eta[g-1][j]->Print();


	if(j==0){

	  l40 = new TLegend(textalign2-.05,texty3-0.05,0.6,texty1-0.05);
	  if(g==1){  
	    l40->AddEntry(RMS_eta[g-1][j],"PbPb Inclusive Jets","lpfe"); 
	    l40->AddEntry(RMS_eta[g][j],"pp Inclusive Jets","p"); 
	  }
	  if(g==3){  
	    l40->AddEntry(RMS_eta[g-1][j],"PbPb Subleading Jets","lpfe"); 
	    l40->AddEntry(RMS_eta[g][j],"pp Subleading Jets","p"); 
	  }
	  if(g==5){  
	    l40->AddEntry(RMS_eta[g-1][j],"PbPb Leading Jets","lpfe"); 
	    l40->AddEntry(RMS_eta[g][j],"pp Leading Jets","p"); 
	  }
	  l40->SetTextFont(43);
	  l40->SetTextSizePixels(tspixels);
	  l40->SetLineColor(kWhite);
	  l40->Draw();
	}

	TLine *line = new TLine(1.,0.,8.,0.);
	line->SetLineStyle(2);
	//	line->Draw("same");

	drawlabels_int_pt2(g,j);

	  
	if(j==3){
	  tex24eta = new TLatex(textalign,texty2,phirangelabel);
	  tex24eta->SetName("tex24phi");
	  tex24eta->SetNDC();
	  tex24eta->SetTextSizePixels(tspixels);
	  tex24eta->Draw();
	}


	RMS_eta[g-1][j]->Write();
	RMS_eta[g][j]->Write();

	RMS_phi[g-1][j]->Write();
	RMS_phi[g][j]->Write();


	//-------------------
	//  Draw  WidthPbPb - Widthpp
	//----------------
      

	RMS_eta_canvas[g]->cd(j+5);
	RMS_PbPb_minus_pp_eta[g][j]->Draw("AP");
	RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetRangeUser(0.8,8.2);
	RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetNdivisions(8);
	RMS_PbPb_minus_pp_eta[g][j]->GetYaxis()->SetNdivisions(204);

	switch(g){
	case 1: 
	  RMS_PbPb_minus_pp_eta[g][j]->SetMarkerStyle(10);
	  RMS_PbPb_minus_pp_eta[g][j]->SetFillColor(90);
	  RMS_PbPb_minus_pp_eta[g][j]->SetLineColor(90);
	  break;
	case 3: 
	  RMS_PbPb_minus_pp_eta[g][j]->SetMarkerStyle(34);
	  RMS_PbPb_minus_pp_eta[g][j]->SetMarkerSize(2);
	  RMS_PbPb_minus_pp_eta[g][j]->SetFillColor(30);
	  RMS_PbPb_minus_pp_eta[g][j]->SetLineColor(30);
	  break;  
	case 5: 
	  RMS_PbPb_minus_pp_eta[g][j]->SetMarkerStyle(21);
	  RMS_PbPb_minus_pp_eta[g][j]->SetFillColor(kOrange-2);
	  RMS_PbPb_minus_pp_eta[g][j]->SetLineColor(kOrange-2);
	  break;
	}


	RMS_PbPb_minus_pp_eta[g][j]->SetMarkerSize(2);
	RMS_PbPb_minus_pp_eta[g][j]->SetTitle("");
	RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetLabelSize(ts);
	RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetTitleSize(tstitle);
	RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->CenterTitle();
	RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetTitleOffset(0.9);
	RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetTitle("Track p_{T} (GeV/c)");



	if(j==0){
	  RMS_PbPb_minus_pp_eta[g][j]->GetYaxis()->SetTitle("#Delta#eta Width_{PbPb} - Width_{pp}");
	  RMS_PbPb_minus_pp_eta[g][j]->GetYaxis()->SetTitleSize(ts);
	  RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetTitleSize(ts);
	  RMS_PbPb_minus_pp_eta[g][j]->GetYaxis()->SetLabelSize(ts2);
	  RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetLabelSize(ts3);
	  RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetLabelOffset(RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->GetLabelOffset()*3.);

	  RMS_PbPb_minus_pp_eta[g][j]->GetXaxis()->SetTitleOffset(1.);
	  RMS_PbPb_minus_pp_eta[g][j]->GetYaxis()->CenterTitle();

	}else{
	  RMS_PbPb_minus_pp_eta[g][j]->GetYaxis()->SetLabelSize(0);
	}


	RMS_PbPb_minus_pp_eta[g][j]->SetMarkerColor(kBlack);
	RMS_PbPb_minus_pp_eta[g][j]->SetMinimum(-0.05);
	RMS_PbPb_minus_pp_eta[g][j]->SetMaximum(0.26);

	RMS_PbPb_minus_pp_eta[g][j]->SetLineWidth(0);

	RMS_PbPb_minus_pp_eta[g][j]->Draw("p Z e2");

	RMS_PbPb_minus_pp_eta[g][j]->Draw("p Z same");




	line->Draw("same");
	
	//-------------------------------------
	//    Start dPhi
	//----------------------------


	TString RMS_PbPb_minus_pp_phi_name = "RMS_PbPb_minus_pp_phi"; RMS_PbPb_minus_pp_phi_name+= g; RMS_PbPb_minus_pp_phi_name+= j;

	RMS_PbPb_minus_pp_phi[g][j]  = (TGraphAsymmErrors*) RMS_phi[g-1][j]->Clone(RMS_PbPb_minus_pp_phi_name);

	for(int k = 0; k<4; k++){
	  
	  float temp_x = pTbin_centers.at(k);
	  
	  float PbPb_val = RMS_phi[g-1][j]->Eval(temp_x);
	  float pp_val = RMS_phi[g][j]->Eval(temp_x);

	  RMS_PbPb_minus_pp_phi[g][j]->SetPoint(k,temp_x,PbPb_val-pp_val);


	}



	//-------------------
	//  Draw  WidthPbPb - Widthpp
	//----------------
   
	RMS_phi_canvas[g]->cd(j+5);
	RMS_PbPb_minus_pp_phi[g][j]->Draw("AP");
	RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetRangeUser(0.8,8.2);
	RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetNdivisions(8);
	RMS_PbPb_minus_pp_phi[g][j]->GetYaxis()->SetNdivisions(204);

	switch(g){
	case 1: 
	  RMS_PbPb_minus_pp_phi[g][j]->SetMarkerStyle(10);
	  RMS_PbPb_minus_pp_phi[g][j]->SetFillColor(90);
	  RMS_PbPb_minus_pp_phi[g][j]->SetLineColor(90);
	  break;
	case 3: 
	  RMS_PbPb_minus_pp_phi[g][j]->SetMarkerStyle(34);
	  RMS_PbPb_minus_pp_phi[g][j]->SetMarkerSize(2);
	  RMS_PbPb_minus_pp_phi[g][j]->SetFillColor(30);
	  RMS_PbPb_minus_pp_phi[g][j]->SetLineColor(30);
	  break;  
	case 5: 
	  RMS_PbPb_minus_pp_phi[g][j]->SetMarkerStyle(21);
	  RMS_PbPb_minus_pp_phi[g][j]->SetFillColor(kOrange-2);
	  RMS_PbPb_minus_pp_phi[g][j]->SetLineColor(kOrange-2);
	  break;
	}


	RMS_PbPb_minus_pp_phi[g][j]->SetMarkerSize(2);
	RMS_PbPb_minus_pp_phi[g][j]->SetTitle("");
	RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetLabelSize(ts);
	RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetTitleSize(tstitle);
	RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->CenterTitle();
	RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetTitleOffset(0.9);
	RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetTitle("Track p_{T} (GeV/c)");


	if(j==0){
	  RMS_PbPb_minus_pp_phi[g][j]->GetYaxis()->SetTitle("#Delta#phi Width_{PbPb} - Width_{pp}");
	  RMS_PbPb_minus_pp_phi[g][j]->GetYaxis()->SetTitleSize(ts);
	  RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetTitleSize(ts);
	  RMS_PbPb_minus_pp_phi[g][j]->GetYaxis()->SetLabelSize(ts2);
	  RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetLabelSize(ts3);
	  RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetTitleOffset(1.);
	  RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->SetLabelOffset( RMS_PbPb_minus_pp_phi[g][j]->GetXaxis()->GetLabelOffset()*3.);
	  RMS_PbPb_minus_pp_phi[g][j]->GetYaxis()->CenterTitle();

	}else{
	  RMS_PbPb_minus_pp_phi[g][j]->GetYaxis()->SetLabelSize(0);
	}

	RMS_PbPb_minus_pp_phi[g][j]->SetMarkerColor(kBlack);
	RMS_PbPb_minus_pp_phi[g][j]->SetMinimum(-0.05);
	RMS_PbPb_minus_pp_phi[g][j]->SetMaximum(0.26);

	RMS_PbPb_minus_pp_phi[g][j]->SetLineWidth(0);

	RMS_PbPb_minus_pp_phi[g][j]->Draw("p Z e2");

	RMS_PbPb_minus_pp_phi[g][j]->Draw("p Z same");






	line->Draw("same");
      

	RMS_phi_canvas[g]->cd(j+1);
	RMS_phi[g-1][j]->Draw("AP");
	RMS_phi[g-1][j]->GetXaxis()->SetRangeUser(0.8,8.2);
	RMS_phi[g-1][j]->GetXaxis()->SetNdivisions(8);
	RMS_phi[g-1][j]->GetYaxis()->SetNdivisions(306);

	switch(g){
	case 1: 
	  RMS_phi[g-1][j]->SetMarkerStyle(10);
	  RMS_phi[g][j]->SetMarkerStyle(4);
	  RMS_phi_ref[g-1][j]->SetMarkerStyle(10);
	  RMS_phi_ref[g][j]->SetMarkerStyle(4);
	  RMS_phi[g-1][j]->SetFillColor(90);
	  RMS_phi[g-1][j]->SetLineColor(90);

	  break;
	case 3: 
	  RMS_phi[g-1][j]->SetMarkerStyle(34);
	  RMS_phi[g-1][j]->SetMarkerSize(2);
	  RMS_phi_ref[g-1][j]->SetMarkerStyle(34);
	  RMS_phi_ref[g-1][j]->SetMarkerSize(2);
	  RMS_phi[g-1][j]->SetFillColor(30);
	  RMS_phi[g-1][j]->SetLineColor(30);
	  RMS_phi[g][j]->SetMarkerStyle(28);
	  break;  
	case 5: 
	  RMS_phi[g-1][j]->SetMarkerStyle(21);
	  RMS_phi_ref[g-1][j]->SetMarkerStyle(21);
	  RMS_phi[g-1][j]->SetFillColor(kOrange-2);
	  RMS_phi[g-1][j]->SetLineColor(kOrange-2);
	  RMS_phi[g][j]->SetMarkerStyle(25);
	  break;
	}

	RMS_phi[g-1][j]->SetTitle("");

	RMS_phi[g][j]->SetMarkerSize(2);
	RMS_phi[g-1][j]->SetMarkerSize(2);

	RMS_phi_ref[g][j]->SetMarkerSize(2);
	RMS_phi_ref[g-1][j]->SetMarkerSize(2);

	if(j==0){
	  RMS_phi[g-1][j]->GetYaxis()->SetTitle("#Delta#phi Width");
	  RMS_phi[g-1][j]->GetYaxis()->SetTitleSize(tstitle);
	  RMS_phi[g-1][j]->GetYaxis()->SetLabelSize(ts);
	  RMS_phi[g-1][j]->GetYaxis()->CenterTitle();

	}else{
	  RMS_phi[g-1][j]->GetYaxis()->SetLabelSize(0);
	}

	RMS_phi[g-1][j]->SetMarkerColor(kBlack);
	RMS_phi[g-1][j]->SetMinimum(-0.02);
	RMS_phi[g-1][j]->SetMaximum(0.75);

	RMS_phi[g-1][j]->SetLineWidth(0);

	RMS_phi[g-1][j]->Draw("p Z e2");
	RMS_phi[g][j]->SetMarkerSize(2);

	RMS_phi[g-1][j]->Draw("p Z same");
	RMS_phi[g][j]->Draw("p X same");


	RMS_phi_ref[g-1][j]->Draw("same P");
	RMS_phi_ref[g][j]->Draw("same P X");
	//	RMS_phi[g-1][j]->Print();

	if(j==0){

	  l40 = new TLegend(textalign2-.05,texty3-0.05,0.6,texty1-0.05);
	  if(g==1){  
	    l40->AddEntry(RMS_phi[g-1][j],"PbPb Inclusive Jets","lpfe"); 
	    l40->AddEntry(RMS_phi[g][j],"pp Inclusive Jets","p"); 
	  }
	  if(g==3){  
	    l40->AddEntry(RMS_phi[g-1][j],"PbPb Subleading Jets","lpfe"); 
	    l40->AddEntry(RMS_phi[g][j],"pp Subleading Jets","p"); 
	  }
	  if(g==5){  
	    l40->AddEntry(RMS_phi[g-1][j],"PbPb Leading Jets","lpfe"); 
	    l40->AddEntry(RMS_phi[g][j],"pp Leading Jets","p"); 

	  }
	  l40->SetTextFont(43);	    
	  l40->SetTextSizePixels(tspixels);
	  l40->SetLineColor(kWhite);
	  l40->Draw();
	}

	//	line->Draw("same");

	drawlabels_int_pt2(g,j);

	  
	if(j==3){
	  tex24phi = new TLatex(textalign,texty2,etarangelabel);
	  tex24phi->SetName("tex24phi");
	  tex24phi->SetNDC();
	  tex24phi->SetTextSizePixels(tspixels);
	  tex24phi->Draw();
	}

    

	//*****************************************
	//*****************************************
	
	// Width of difference as function of cent.
	
	//*****************************************
	//*****************************************
	


	//	RMS_PbPb_minus_pp_phi[g][j]->Print();
	cout<<" starting diff fits"<<endl;
         
	//------------------------
	//   RMS of difference
	//-----------------------

	dummy->cd();

   	for(int i = 0; i<2; i++){
	  ///--------------
	  //   Difference
	  //---------------
	
	  double_gaus->ReleaseParameter(0);
	  double_gaus->ReleaseParameter(2);
	  double_gaus->ReleaseParameter(1);
	  double_gaus->ReleaseParameter(3);
	  double_gaus->ReleaseParameter(4);
	  //double_gaus->SetParameter(4,0.);

	  if(i<2){	double_gaus->SetParLimits(0,0.0,10.);  }
	  //	  double_gaus->SetParLimits(1,0.1,range);
	  if(i<2){	double_gaus->SetParLimits(2,0.0,10.);  }

	  double_gaus->SetParLimits(1,0.2,0.67*range);
	  double_gaus->SetParLimits(3,0.2,0.67*range);
	  
	  double_gaus->SetParameter(1,0.1);
	  double_gaus->SetParameter(3,0.6);
	  double_gaus->SetParameter(0,.5);
	  double_gaus->SetParameter(2,.5);
	  	  
	  if(i>1){
	    double_gaus->FixParameter(2,0.);
	    double_gaus->SetParLimits(1,0.2,0.5*range);
	  }
	  switch(i){
	  case 0:
	    RMS_eta_demo_val_0.push_back(PbPb_pp_eta[g][i][j]->GetRMS());
	    break;
	  case 1:
	    RMS_eta_demo_val_1.push_back(PbPb_pp_eta[g][i][j]->GetRMS());
	    break;
	  case 2:
	    RMS_eta_demo_val_2.push_back(PbPb_pp_eta[g][i][j]->GetRMS());
	    break;
	  case 3:
	    RMS_eta_demo_val_3.push_back(PbPb_pp_eta[g][i][j]->GetRMS());
	    break;
	  }


	  //	  double_gaus->FixParameter(2,0.);


	  PbPb_pp_eta[g][i][j]->Fit(double_gaus,"","",-range-1.,range+1.);

	  flat_line->SetParameter(0,double_gaus->GetParameter(4));

	  PbPb_pp_eta[g][i][j]->Add(flat_line);
	  PbPb_pp_eta_up[g][i][j]->Add(flat_line);
	  PbPb_pp_eta_down[g][i][j]->Add(flat_line);

	  PbPb_pp_eta[g][i][j]->Fit(double_gaus,"","",-range-1.,range+1.);

	  double_gaus->SetParameter(4,0.);
	  
	  /*
	    par0 = double_gaus->GetParameter(0);
	    par1 = double_gaus->GetParameter(1);
	    par2 = double_gaus->GetParameter(2);
	    par3 = double_gaus->GetParameter(3);
	    temp_sigma = TMath::Sqrt((par0*par1*par1+par2*par3*par3)/(par0+par2));
	  */
	  

	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }

	
	  
	  RMS_diff_val_eta.push_back(temp_sigma);
	

	  PbPb_pp_eta_up[g][i][j]->Fit(double_gaus,"","",-range-1.,range+1.);

	  double_gaus->SetParameter(4,0.);
	  /*
	    par0 = double_gaus->GetParameter(0);
	    par1 = double_gaus->GetParameter(1);
	    par2 = double_gaus->GetParameter(2);
	    par3 = double_gaus->GetParameter(3);
	    temp_sigma = TMath::Sqrt((par0*par1*par1+par2*par3*par3)/(par0+par2));
	  */
	 
	
	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	   

	
	  error_up = temp_sigma-RMS_diff_val_eta.at(i);

	 
	  
	  PbPb_pp_eta_down[g][i][j]->Fit(double_gaus,"","",-range,range);

	  double_gaus->SetParameter(4,0.);
	  /*
	    par0 = double_gaus->GetParameter(0);
	    par1 = double_gaus->GetParameter(1);
	    par2 = double_gaus->GetParameter(2);
	    par3 = double_gaus->GetParameter(3);
	    temp_sigma = TMath::Sqrt((par0*par1*par1+par2*par3*par3)/(par0+par2));
	  */
	  

	  full_int = double_gaus->Integral(0,range);
	  temp_int = 0;
	  temp_sigma = 0;

	  counter = 0;
	  while(temp_int <0.67*full_int){
	    temp_sigma = range*counter/50000;
	    temp_int = double_gaus->Integral(0, temp_sigma);
	    counter++;
	  }
	 
	  error_down = temp_sigma-RMS_val_eta.at(i);
 
	  if(error_up>0&&error_down<0){
	    RMS_diff_err_eta_up.push_back(TMath::Abs(error_up));
	    RMS_diff_err_eta_down.push_back(TMath::Abs(error_down));
	  }else if(error_up<0&&error_down>0){
	    RMS_diff_err_eta_down.push_back(TMath::Abs(error_up));
	    RMS_diff_err_eta_up.push_back(TMath::Abs(error_down));
	  } else {
	    cerr<<"Both errors go same way"<<endl;
	    RMS_diff_err_eta_up.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    RMS_diff_err_eta_down.push_back(max(TMath::Abs(error_up),TMath::Abs(error_down)));
	    //return -1;
	  }
	 
	}

      
	RMS_diff_eta[g][j] = new TGraphAsymmErrors(pTbin_centers.size(),&pTbin_centers[0],&RMS_diff_val_eta[0],&pTbin_errors[0],&pTbin_errors[0],&RMS_diff_err_eta_down[0],&RMS_diff_err_eta_up[0]);
      

	TString RMS_diff_eta_name = "RMS_diff_eta"; RMS_diff_eta_name+= g; RMS_diff_eta_name+= j;

	RMS_diff_eta[g][j]->SetName(RMS_diff_eta_name);


      }  //close only pp
      
      cout<<"made it here "<<g<<" "<<" "<<j<<endl;

      for(int i = 0; i<4; i++){
	
	check_fits_eta[g]->cd(4*i+j+1);

	check_new_eta_rebin[g][i][j]->SetAxisRange(-2.499,2.4999);



	check_new_eta_rebin[g][i][j]->SetMarkerColor(kBlack);
	check_new_eta_rebin[g][i][j]->SetLineColor(kBlack);
	check_new_eta_rebin[g][i][j]->SetMarkerStyle(20);
	check_new_eta_rebin[g][i][j]->SetMarkerSize(1);
	check_new_eta_rebin[g][i][j]->Draw();

	if(g==0||g==2||g==4){	
	  //	  check_new_eta_syst[g][i][j]->Draw("same e2");
	  check_new_eta_up[g][i][j]->SetMarkerSize(0.);
	  check_new_eta_up[g][i][j]->Draw("same");
	  check_new_eta_down[g][i][j]->SetMarkerSize(0.);
	  check_new_eta_down[g][i][j]->Draw("same");
	}
	
	gaus_tot[g][i][j]->SetLineWidth(2);
	gaus_tot[g][i][j]->SetLineColor(kBlue);


	gaus_1[g][i][j]->SetLineWidth(2);
	gaus_1[g][i][j]->SetLineStyle(2);
	gaus_1[g][i][j]->SetLineColor(kBlack);
	gaus_1[g][i][j]->Draw("same");
	
	gaus_2[g][i][j]->SetLineWidth(2);
	gaus_2[g][i][j]->SetLineStyle(2);
	gaus_2[g][i][j]->SetLineColor(kBlack);
	gaus_2[g][i][j]->Draw("same");
	
	
	
	check_new_eta_rebin[g][i][j]->Draw("same");
	gaus_tot[g][i][j]->Draw("same");

      
	if(i==0){	drawlabels(g,i,j); 
	}else if(j==0){

	  TLatex *tex01;
	  switch(i){
	  case 0:
	    tex01 = new TLatex(textalign2,texty1,"1<p_{T}^{assoc}<2 GeV/c");
	    break;
	  case 1:
	    tex01 = new TLatex(textalign2,texty1,"2<p_{T}^{assoc}<3 GeV/c");
	    break;
	  case 2:
	    tex01 = new TLatex(textalign2,texty1,"3<p_{T}^{assoc}<4 GeV/c");
	    break;
	  case 3:
	    tex01 = new TLatex(textalign2,texty1,"4<p_{T}^{assoc}<8 GeV/c");
	    break;
	  }
	    
	  tex01->SetName("tex01");
	  tex01->SetNDC();
	  tex01->SetTextSize(ts2);
	  tex01->SetLineWidth(2);
	  tex01->Draw("same");
	}



	check_fits_phi[g]->cd(4*i+j+1);

	check_new_phi_rebin[g][i][j]->SetAxisRange(-1.499,1.4999);



	check_new_phi_rebin[g][i][j]->Draw();

	if(g==0||g==2||g==4){	
	  check_new_phi_syst[g][i][j]->Draw("same e2");
	  check_new_phi_up[g][i][j]->SetMarkerSize(0.);
	  check_new_phi_up[g][i][j]->Draw("same");
	  check_new_phi_down[g][i][j]->SetMarkerSize(0.);
	  check_new_phi_down[g][i][j]->Draw("same");
	}
	/*
	  gaus_tot[g][i][j]->SetLineWidth(2);
	  gaus_tot[g][i][j]->SetLineColor(kBlue);


	  gaus_1[g][i][j]->SetLineWidth(2);
	  gaus_1[g][i][j]->SetLineStyle(2);
	  gaus_1[g][i][j]->SetLineColor(kBlack);
	  gaus_1[g][i][j]->Draw("same");
	
	  gaus_2[g][i][j]->SetLineWidth(2);
	  gaus_2[g][i][j]->SetLineStyle(2);
	  gaus_2[g][i][j]->SetLineColor(kBlack);
	  gaus_2[g][i][j]->Draw("same");
	
	*/
	check_new_phi_rebin[g][i][j]->GetFunction("double_gaus")->SetLineColor(kBlue);	

	check_new_phi_rebin[g][i][j]->Draw("same");
	//	gaus_tot[g][i][j]->Draw("same");

      
	if(i==0){	drawlabels(g,i,j); 
	}else if(j==0){

	  TLatex *tex01;
	  switch(i){
	  case 0:
	    tex01 = new TLatex(textalign2,texty1,"1<p_{T}^{assoc}<2 GeV/c");
	    break;
	  case 1:
	    tex01 = new TLatex(textalign2,texty1,"2<p_{T}^{assoc}<3 GeV/c");
	    break;
	  case 2:
	    tex01 = new TLatex(textalign2,texty1,"3<p_{T}^{assoc}<4 GeV/c");
	    break;
	  case 3:
	    tex01 = new TLatex(textalign2,texty1,"4<p_{T}^{assoc}<8 GeV/c");
	    break;
	  }
	    
	  tex01->SetName("tex01");
	  tex01->SetNDC();
	  tex01->SetTextSize(ts2);
	  tex01->SetLineWidth(2);
	  tex01->Draw("same");
	}

	dummy->cd();

      	
	if((g==1||g==3||g==5)&&i<2){

	  check_fits_diff[g]->cd(4*i+j+1);
	  
	  
	  PbPb_pp_eta[g][i][j]->SetMinimum(-1.);
	  if(i==0){
	    PbPb_pp_eta[g][i][j]->SetMaximum(7.);
	  }  else{
	    PbPb_pp_eta[g][i][j]->SetMaximum(4.);
	  }
	  
	  PbPb_pp_eta[g][i][j]->SetMarkerColor(kBlack);
	  PbPb_pp_eta[g][i][j]->SetLineColor(kBlack);
	  PbPb_pp_eta[g][i][j]->SetMarkerStyle(20);
	  PbPb_pp_eta[g][i][j]->SetMarkerSize(1);

	  PbPb_pp_eta[g][i][j]->Draw();
	  //	  PbPb_pp_eta_syst[g][i][j]->Draw("same e2");


	  PbPb_pp_eta_up[g][i][j]->SetMarkerSize(0);
	  PbPb_pp_eta_down[g][i][j]->SetMarkerSize(0);
	  PbPb_pp_eta_up[g][i][j]->Draw("same");
	  PbPb_pp_eta_down[g][i][j]->Draw("same");

	  PbPb_pp_eta[g][i][j]->GetFunction("double_gaus")->SetLineColor(kBlue);
	  PbPb_pp_eta[g][i][j]->Draw("same");

	  if(i==0){  drawlabels(g,i,j); 
	  }else if(j==0){

	    TLatex *tex01;
	    switch(i){
	    case 0:
	      tex01 = new TLatex(textalign2,texty1,"1<p_{T}^{assoc}<2 GeV/c");
	      break;
	    case 1:
	      tex01 = new TLatex(textalign2,texty1,"2<p_{T}^{assoc}<3 GeV/c");
	      break;
	    case 2:
	      tex01 = new TLatex(textalign2,texty1,"3<p_{T}^{assoc}<4 GeV/c");
	      break;
	    case 3:
	      tex01 = new TLatex(textalign2,texty1,"4<p_{T}^{assoc}<8 GeV/c");
	      break;
	    }
	    
	    tex01->SetName("tex01");
	    tex01->SetNDC();
	    tex01->SetTextSize(ts2);
	    tex01->SetLineWidth(2);
	    tex01->Draw("same");
	  }

	  lineEta = new TLine(-2.,0.,2.,0.);
	  lineEta->SetLineStyle(2.);
	  lineEta->Draw();

	}

	dummy->cd();
      

      }  //end i

    } // end j

    RMS_eta_canvas[g]->cd(0);
    TLatex *canvas_title = new TLatex(0.06,0.95,"CMS Preliminary");
    canvas_title->SetTextSizePixels(tspixels);
    canvas_title->SetTextFont(63);
    canvas_title->Draw();

    TLatex *canvas_title2 = new TLatex(0.292,0.95,"PbPb 166 #mub^{-1} (2.76 TeV)               pp 5.3 pb^{-1} (2.76 TeV)");
    canvas_title2->SetTextSizePixels(tspixels);
    canvas_title2->Draw();

    RMS_phi_canvas[g]->cd(0);
    canvas_title->Draw();
    canvas_title2->Draw();

    if(g%2==0){
      check_fits_eta[g]->SaveAs((TString)"Width_Check_Fits_Eta_PbPb_"+datalabel+".png");
      check_fits_eta[g]->SaveAs((TString)"RMS_Check_Fits_PbPb_"+datalabel+".pdf");
      check_fits_phi[g]->SaveAs((TString)"Width_Check_Fits_Phi_PbPb_"+datalabel+".png");
      check_fits_phi[g]->SaveAs((TString)"RMS_Check_Fits_Phi_PbPb_"+datalabel+".pdf");
    }else{
      check_fits_eta[g]->SaveAs((TString)"RMS_Check_Fits_Eta_pp"+datalabel+".pdf");
      check_fits_diff[g]->SaveAs((TString)"RMS_Check_Fits_Diff_"+datalabel+".pdf");
      check_fits_eta[g]->SaveAs((TString)"Width_Check_Fits_Eta_pp_"+datalabel+".png");
      check_fits_diff[g]->SaveAs((TString)"Width_Check_Fits_Diff_"+datalabel+".png");
      check_fits_phi[g]->SaveAs((TString)"RMS_Check_Fits_Phi_pp"+datalabel+".pdf");
      check_fits_diff[g]->SaveAs((TString)"RMS_Check_Fits_Diff_"+datalabel+".pdf");
      check_fits_phi[g]->SaveAs((TString)"Width_Check_Fits_Phi_pp_"+datalabel+".png");
      check_fits_diff[g]->SaveAs((TString)"Width_Check_Fits_Diff_"+datalabel+".png");
      RMS_eta_canvas[g]->SaveAs((TString)("RMS_Eta_"+datalabel+".pdf"));
      RMS_phi_canvas[g]->SaveAs((TString)("RMS_Phi_"+datalabel+".pdf"));
      RMS_eta_canvas[g]->SaveAs((TString)("Width_Eta_"+datalabel+".png"));
      RMS_phi_canvas[g]->SaveAs((TString)("Width_Phi_"+datalabel+".png"));
    }
   
  
  
    //******************************************************

      
    //---------------------------------------
    //Now for width by centrality plots
    //--------------------------------------
      
    cout<<"ready to start widths by centrality"<<endl;

    for(int i = 0; i<2; i++){
      
      RMS_val_cent_eta.clear();
      RMS_err_cent_eta_up.clear();
      RMS_err_cent_eta_down.clear();


      for(int k=0; k<4; k++){
	  
	evalpt = pTbin_centers.at(i);
	RMS_val_cent_eta.push_back( RMS_eta[g][k]->Eval(evalpt));

	
	//	eyl = RMS_eta[g][k]->GetErrorYlow(i);
	//eyh = RMS_eta[g][k]->GetErrorYhigh(i);
	
	eyl = 0.;
	eyh = 0.;
		
	RMS_err_cent_eta_up.push_back(eyh);
	RMS_err_cent_eta_down.push_back(eyl);

      }
    
      RMS_eta_cent[g][i] = new TGraphAsymmErrors(4,&cent_centers[0], &RMS_val_cent_eta[0], &cent_errors[0], &cent_errors[0],  &RMS_err_cent_eta_down[0], &RMS_err_cent_eta_up[0]);

      TString	RMS_eta_cent_name = "RMS_eta_cent"; 	RMS_eta_cent_name+= g; RMS_eta_cent_name+= i; 
      RMS_eta_cent[g][i]->SetName(RMS_eta_cent_name );
    
    
      RMS_diff_val_cent_eta.clear();
      RMS_diff_err_cent_eta_up.clear();
      RMS_diff_err_cent_eta_down.clear();

      cout<<"made it here starting rms eta cent"<<endl;

      if(g%2!=0){
	

	cout<<"Here"<<endl;

	for(int k=0; k<4; k++){
	  
	  evalpt = pTbin_centers.at(i);
	  RMS_diff_val_cent_eta.push_back(RMS_diff_eta[g][k]->Eval(evalpt));

	  eyl = RMS_diff_eta[g][k]->GetErrorYlow(i);
	  eyh = RMS_diff_eta[g][k]->GetErrorYhigh(i);
	
		
	  RMS_diff_err_cent_eta_up.push_back(eyh);
	  RMS_diff_err_cent_eta_down.push_back(eyl);

	  cout<<g<<" "<<i<<" "<<k<<" "<<evalpt<<" "<<RMS_diff_eta[g][k]->Eval(evalpt)<<" "<<eyl<<" "<<eyh<<endl;

	}
	RMS_diff_eta_cent[g][i] = new TGraphAsymmErrors(4,&cent_centers[0], &RMS_diff_val_cent_eta[0], &cent_errors[0], &cent_errors[0],  &RMS_diff_err_cent_eta_down[0], &RMS_diff_err_cent_eta_up[0]);

	cout<<"made RMS_diff_eta_cent"<<endl;

	switch(i){
	case 0:
	  RMS_eta_demo[g][i] = new TGraphAsymmErrors(4,&cent_centers[0], &RMS_eta_demo_val_0[0], &cent_errors[0], &cent_errors[0],  &RMS_diff_err_cent_eta_down[0], &RMS_diff_err_cent_eta_up[0]);
	  break;
	case 1:
	  RMS_eta_demo[g][i] = new TGraphAsymmErrors(4,&cent_centers[0], &RMS_eta_demo_val_1[0], &cent_errors[0], &cent_errors[0],  &RMS_diff_err_cent_eta_down[0], &RMS_diff_err_cent_eta_up[0]);
	  break;
	case 2:
	  RMS_eta_demo[g][i] = new TGraphAsymmErrors(4,&cent_centers[0], &RMS_eta_demo_val_2[0], &cent_errors[0], &cent_errors[0],  &RMS_diff_err_cent_eta_down[0], &RMS_diff_err_cent_eta_up[0]);
	  break;
	case 3:
	  RMS_eta_demo[g][i] = new TGraphAsymmErrors(4,&cent_centers[0], &RMS_eta_demo_val_3[0], &cent_errors[0], &cent_errors[0],  &RMS_diff_err_cent_eta_down[0], &RMS_diff_err_cent_eta_up[0]);
	  break;
	}


	TString	RMS_diff_eta_cent_name = "RMS_diff_eta_cent"; 	RMS_diff_eta_cent_name+= g; RMS_diff_eta_cent_name+= i; 
	RMS_diff_eta_cent[g][i]->SetName(RMS_diff_eta_cent_name );


    
	cRMS_eta_cent[g]->cd(i+1);

	
	TString histnameblank = "blank_hist";
	histnameblank+=g;
	histnameblank+=i;

	blank[i] = new TH1D(histnameblank,"",4,xAxis);
	
	//Markers etc. for all histograms at once
	//----------------------------------------

     

	RMS_eta_cent[g-1][i]->SetMarkerSize(2);
	RMS_eta_cent[g][i]->SetMarkerSize(2);
	RMS_eta_cent[g][i]->SetMarkerColor(1);
	RMS_eta_cent[g][i]->SetLineColor(1);

	switch(g){
	case 1:
	  RMS_eta_cent[g-1][i]->SetFillColor(90);
	  RMS_eta_cent[g-1][i]->SetMarkerStyle(10);
	  RMS_eta_cent[g][i]->SetMarkerStyle(4);
	  break;     
	case 3:
	  RMS_eta_cent[g-1][i]->SetFillColor(30);
	  RMS_eta_cent[g-1][i]->SetMarkerStyle(34);
	  RMS_eta_cent[g][i]->SetMarkerStyle(28);
	  break;  
	case 5:
	  RMS_eta_cent[g-1][i]->SetFillColor(kOrange-2);
	  RMS_eta_cent[g-1][i]->SetMarkerStyle(21);
	  RMS_eta_cent[g][i]->SetMarkerStyle(25);

	  break;     
	}

	//Plot aesthetics for every canvas.
	//----------------------------------
	
	blank[i]->SetMinimum(0.);
	blank[i]->SetMaximum(0.9);
	blank[i]->GetXaxis()->SetTitle("Centrality (%)");
	blank[i]->GetXaxis()->SetTitleOffset(1.1);
	blank[i]->GetXaxis()->CenterTitle(true);
	blank[i]->GetXaxis()->SetTitleSize(ts);

	blank[i]->GetYaxis()->SetTitle("#Delta#eta Width");
	blank[i]->GetYaxis()->SetTitleSize(0.);
	blank[i]->GetYaxis()->CenterTitle(true);
	//	blank[i]->GetYaxis()->SetLabelOffset(yoffset);
	blank[i]->GetYaxis()->SetLabelSize(0.);
   
	blank[i]->GetYaxis()->SetTickLength(0.025);

	blank[i]->GetXaxis()->SetBinLabel(1,"50-100");
	blank[i]->GetXaxis()->SetBinLabel(2,"30-50");
	blank[i]->GetXaxis()->SetBinLabel(3,"10-30");
	blank[i]->GetXaxis()->SetBinLabel(4," 0-10");
    
	blank[i]->GetXaxis()->SetLabelSize(0.09);
	blank[i]->GetXaxis()->SetLabelOffset(0.015);
	//	blank[i]->GetXaxis()->SetTicks("+-");
	blank[i]->GetXaxis()->LabelsOption("h");
	blank[i]->GetXaxis()->SetTickLength(0.025);

	switch(i){
	case 0: 
	  gPad->SetLeftMargin(0.1);
	  blank[i]->GetYaxis()->SetTitleSize(0.09);
	  blank[i]->GetXaxis()->SetTitleOffset(1.2);
	  blank[i]->GetYaxis()->SetTitleOffset(0.5);
	  blank[i]->GetXaxis()->SetTitleSize(0.06);
	  blank[i]->SetLabelSize(0.95*blank[i]->GetXaxis()->GetLabelSize());
	  blank[i]->GetYaxis()->SetLabelSize(ts);
	  break;
	case 3:
	  gPad->SetRightMargin(0.02);
	  break;
	default:
	  break;
	}
	blank[i]->Draw();
	RMS_eta_cent[g-1][i]->Draw("same p e2");
	RMS_eta_cent[g][i]->SetLineWidth(0);
	RMS_eta_cent[g][i]->Draw("same p X");
	
	if(i==0){ 
	  l40 = new TLegend(textalign2,texty2-.05,0.6,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){
	    l40->AddEntry(RMS_eta_cent[g-1][i],"Inclusive PbPb","lpfe");
	    l40->AddEntry(RMS_eta_cent[g][i],"Inclusive pp","lpe");
	  }
	  if(g==3){
	    l40->AddEntry(RMS_eta_cent[g-1][i],"Subleading PbPb","lpfe");
	    l40->AddEntry(RMS_eta_cent[g][i],"Subleading pp","lpe");
	  }
	  if(g==5){
	    l40->AddEntry(RMS_eta_cent[g-1][i],"Leading PbPb","lpfe");
	    l40->AddEntry(RMS_eta_cent[g][i],"Leading pp","lpe");
	  }
	  l40->Draw("same");
	}
	drawlabels_int_cent(g,i);
	 
	lineCent = new TLine(-100.,0.,0.,0.);
	lineCent->SetLineStyle(2);
	lineCent->SetLineWidth(2);
	lineCent->Draw("same");
	
	//---------------------
	//  RMS of diff
	//--------------------
      

	cout<<"here starting rms of diff"<<endl;      

	  cRMS_eta_cent[g]->cd(i+1);


	  TString histnameblank2 = "blank2_hist";
	  histnameblank2+=g;
	  histnameblank2+=i;

	  blank2[i] = new TH1D(histnameblank2,"",4,xAxis);
	
	  //Markers etc. for all histograms at once
	  //----------------------------------------


	  RMS_diff_eta_cent[g][i]->SetMarkerSize(2);

	  switch(g){
	  case 1:
	    RMS_diff_eta_cent[g][i]->SetFillColor(90);
	    RMS_diff_eta_cent[g][i]->SetMarkerStyle(10);
	    break;     
	  case 3:
	    RMS_diff_eta_cent[g][i]->SetFillColor(30);
	    RMS_diff_eta_cent[g][i]->SetMarkerStyle(34);
	    break;  
	  case 5:
	    RMS_diff_eta_cent[g][i]->SetFillColor(kOrange-2);
	    RMS_diff_eta_cent[g][i]->SetMarkerStyle(21);
	    break;     
	  }

	  //Plot aesthetics for every canvas.
	  //----------------------------------
	
	  blank2[i]->SetMinimum(-0.01);
	  blank2[i]->SetMaximum(1.35);
	  blank2[i]->GetXaxis()->SetTitle("Centrality (%)");
	  blank2[i]->GetXaxis()->CenterTitle(true);
	  blank2[i]->GetXaxis()->SetTitleSize(ts);

	  blank2[i]->GetYaxis()->SetTitle("#Delta#eta Width");
	  blank2[i]->GetYaxis()->SetTitleSize(0.);
	  blank2[i]->GetYaxis()->CenterTitle(true);
	  //	blank2[i]->GetYaxis()->SetLabelOffset(yoffset);
	  blank2[i]->GetYaxis()->SetLabelSize(0.);
   
	  blank2[i]->GetYaxis()->SetTickLength(0.025);

	  blank2[i]->GetXaxis()->SetBinLabel(1,"50-100");
	  blank2[i]->GetXaxis()->SetBinLabel(2,"30-50");
	  blank2[i]->GetXaxis()->SetBinLabel(3,"10-30");
	  blank2[i]->GetXaxis()->SetBinLabel(4," 0-10");
    
	  blank2[i]->GetXaxis()->SetLabelSize(0.09);
	  blank2[i]->GetXaxis()->SetLabelOffset(0.015);
	  //	blank2[i]->GetXaxis()->SetTicks("+-");
	  blank2[i]->GetXaxis()->LabelsOption("h");

	  blank2[i]->GetXaxis()->SetTickLength(0.0);
	  blank2[i]->GetYaxis()->SetNdivisions(8);
	

	  if(i==0){
	    gPad->SetLeftMargin(0.1);
	    blank2[i]->GetYaxis()->SetTitleSize(tstitle+0.01);
	    blank2[i]->GetXaxis()->SetTitleOffset(1.);
	    blank2[i]->GetYaxis()->SetTitleOffset(0.5);
	    blank2[i]->GetXaxis()->SetTitleSize(ts);
	    blank2[i]->GetYaxis()->SetLabelSize(tstitle);
	  
	  }

	  blank2[i]->Draw();

	  RMS_eta_cent[g-1][i]->SetLineWidth(2);
	  RMS_eta_cent[g-1][i]->SetMarkerStyle(1);
	  RMS_eta_cent[g-1][i]->SetFillColor(kBlack);
	  RMS_diff_eta_cent[g][i]->Draw("same p e2");
	  RMS_eta_cent[g-1][i]->Draw("same  ");

	  
	  if(i==0){ 
	    l40 = new TLegend(0.4,texty1+0.05,0.9,texty2);
	    l40->SetName("l40");
	    l40->SetTextSize(ts+0.01);
	    l40->SetFillColor(kWhite);
	    l40->SetLineColor(kWhite);
	    l40->AddEntry(RMS_diff_eta_cent[g][i],"Width of PbPb - pp","pfe");
	    l40->AddEntry(RMS_eta_cent[g-1][i],"Width of PbPb correlation","l");
	    l40->Draw("same");



	    TLatex *pTtex = new TLatex(0.15,texty2,"1<p_{T}^{assoc}<2 GeV/c");
	    pTtex->SetTextSizePixels(tspixels);
	    pTtex->SetNDC();
	    pTtex->SetTextFont(43);
	    pTtex->Draw();

	    TLatex *datatex = new TLatex(0.15,texty1,(TString)(datalabel+" Jets"));
	    datatex->SetTextSizePixels(tspixels);
	    datatex->SetTextFont(43);
	    datatex->SetNDC();
	    datatex->Draw();

	  }else{
	  
	    TLatex *pTtex = new TLatex(0.05,texty3,"2<p_{T}^{assoc}<3");
	    pTtex->SetTextSizePixels(tspixels);
	    pTtex->SetNDC();
	    pTtex->SetTextFont(43);
	    pTtex->Draw();

	    TLatex *tex2 = new TLatex(0.05,texty2,"anti-k_{T} jets,  R=0.3,  |#eta_{jet}|<1.6,  Projected |#Delta#phi|<1");
	    tex2->SetTextSizePixels(tspixels);
	    tex2->SetNDC();
	    tex2->SetTextFont(43);
	    tex2->Draw();

	    if(g<2){
	      TLatex *jettex = new TLatex(0.05,texty1,"p_{T,jet}>120 GeV/c");
	      jettex->SetTextSizePixels(tspixels);
	      jettex->SetNDC();
	      jettex->SetTextFont(43);
	      jettex->Draw();
	    }else{
	      TLatex *jettex = new TLatex(0.05,texty1,"p_{T,jet1}>120 GeV/c,  p_{T,jet2}>50 GeV/c,  #Delta#phi_{1,2}>5#pi/6");
	      jettex->SetTextSizePixels(tspixels);
	      jettex->SetNDC();
	      jettex->SetTextFont(43);
	      jettex->Draw();
	    }



	  }

       
		

	  TLine *line1 = new TLine(-50,-0.05,-50.,0.05);
	  TLine *line2 = new TLine(-30,-0.05,-30.,0.05);
	  TLine *line3 = new TLine(-10,0.-0.05,-10.,0.05);

	  if(g!=5){
	    line1->Draw();
	    line2->Draw();
	    line3->Draw();
	  }


	  gPad->RedrawAxis();	


	  if(g==5){
	    cRMS_eta_cent[g]->cd(i+3);

	    TString histnameblank3 = "blank3_hist";
	    histnameblank3+=g;
	    histnameblank3+=i;
	  
	    blank3[i] = (TH1D*)blank2[i]->Clone(histnameblank3);

	    if(i==0){ 
	      gPad->SetLeftMargin(0.1);
	      blank3[i]->GetYaxis()->SetTitleSize(tstitle);
	      blank3[i]->GetYaxis()->SetTitleOffset(0.55);
	      blank3[i]->GetYaxis()->SetLabelSize(ts);
	      blank3[i]->Draw();
	    }
	    blank3[i]->Draw();
	    RMS_diff_eta_cent[3][i]->Draw("same p e2");
	    RMS_eta_cent[2][i]->Draw("same  ");

	    line1->Draw();
	    line2->Draw();
	    line3->Draw();

	    // lineCent->Draw("same");

	  
	    if(i==0){ 
	      l40 = new TLegend(0.4,texty1+0.05,0.9,texty2);
	      l40->SetName("l40");
	      l40->SetTextSize(ts2);
	      l40->SetFillColor(kWhite);
	      l40->SetLineColor(kWhite);
	      l40->AddEntry(RMS_diff_eta_cent[3][i],"Width of PbPb - pp","pfe");
	      l40->AddEntry(RMS_eta_cent[2][i],"Width of PbPb correlation","l");
	      l40->Draw("same");



	      TLatex *datatex2 = new TLatex(0.15,texty1,"Subleading Jets");
	      datatex2->SetTextSizePixels(tspixels);
	      datatex2->SetNDC();
	      datatex2->SetTextFont(43);
	      datatex2->Draw();


	      TLatex *zero = new TLatex(0.0805,.98,"0");
	      zero->SetNDC();
	      zero->SetTextFont(43);
	      zero->SetTextSizePixels(tspixels);
	      zero->Draw();



	    }
	  

	    gPad->RedrawAxis();
	  
	    cRMS_eta_cent[g]->cd(i+1);
	    //	  blank2[i]->Draw("same");

	  gPad->RedrawAxis();

	}

      }
           		
    } //close i


    if(g==1||g==3||g==5){
      cRMS_eta_cent[g]->cd(0);
      
      if(g==5){
	TLatex *canvas_title = new TLatex(0.06,0.95,"CMS Preliminary");
	canvas_title->SetTextSizePixels(tspixels);
	canvas_title->SetTextFont(63);
	canvas_title->Draw();

	TLatex *canvas_title2 = new TLatex(0.28,0.95,"PbPb 166 #mub^{-1} (2.76 TeV)                       pp 5.3 pb^{-1} (2.76 TeV)");
	canvas_title2->SetTextSizePixels(tspixels);
	canvas_title2->SetTextFont(43);
	canvas_title2->Draw();
	
      }else{
	TLatex *canvas_title = new TLatex(0.06,0.9,"CMS Preliminary");
	canvas_title->SetTextSizePixels(tspixels);
	canvas_title->SetTextFont(63);
	canvas_title->Draw();

	TLatex *canvas_title2 = new TLatex(0.28,0.9,"PbPb 166 #mub^{-1} (2.76 TeV)                       pp 5.3 pb^{-1} (2.76 TeV)");
	canvas_title2->SetTextSizePixels(tspixels);
	canvas_title2->SetTextFont(43);
	canvas_title2->Draw();


      }

      TString RMSSaveName_cent = in_name;
      RMSSaveName_cent.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      RMSSaveName_cent.ReplaceAll("Yield_pp","RMS_Eta_Cent");
      RMSSaveName_cent +=".pdf";
      cRMS_eta_cent[g]->SaveAs(RMSSaveName_cent);
    
      RMSSaveName_cent.ReplaceAll("pdf","png");
      RMSSaveName_cent.ReplaceAll("RMS","WidthOfDifference");
      cRMS_eta_cent[g]->SaveAs(RMSSaveName_cent);
    
      fout[g]->cd();
      fout[g]->Write();
    
    } 

  } //close g 

  cout<<" this should be the last thing...."<<endl;

  return 0;

} //Close main loop

