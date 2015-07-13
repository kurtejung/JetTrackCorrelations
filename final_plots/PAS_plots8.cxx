
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

Int_t PAS_plots8(bool draw_ref = kFALSE){

#include "../HIN_14_016_universals.h"

  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.15);
  gStyle->SetPadLeftMargin  (0.16);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);
  gStyle->SetTextFont(43);
  gStyle->SetCanvasBorderMode(0);
 
 

  TFile *fin[12];
  TFile *fin_ref[12];
  TFile *fin2[12];

  TFile *fbkgfit = new TFile("../bg_fits/PbPb_Leading_Yield_and_Bkg.root","READ");
  TFile *fbkgsummed = new TFile("../me_correct/PbPb_Leading_Correlations.root","READ");
 
  TFile *fbkgfit_sub = new TFile("../bg_fits/PbPb_SubLeading_Yield_and_Bkg.root","READ");
  TFile *fbkgsummed_sub = new TFile("../me_correct/PbPb_SubLeading_Correlations.root","READ");
  


  Double_t xAxis[5] = {-100,-50,-30,-10,0}; 
  TH1D* int_cent[12][4];
  TH1D* blank[4];
  TH1D* blank2[4];

  TH1D *check_new_phi_rebin[12][5][4];
  TH1D *check_new_phi_ref[12][5][4];
  TH1D *check_new_phi_syst[12][5][4];
  
  TH1D *check_new_eta_rebin[12][5][4];
  TH1D *check_new_eta_ref[12][5][4];
  TH1D *check_new_eta_syst[12][5][4];

  TH1D *PbPb_pp_phi[12][5][4];
  TH1D *PbPb_pp_phi_syst[12][5][4];

  TH1D *PbPb_pp_eta[12][5][4];
  TH1D *PbPb_pp_eta_syst[12][5][4];

  TH1D *PbPb_pp_phi_ref[12][5][4];
  TH1D *PbPb_pp_eta_ref[12][5][4];

  TH1D *HYDJET_PYTHIA_eta[12][5][4];
  TH1D *HYDJET_PYTHIA_phi[12][5][4];

 
  TF1 *gaus_phi[12][5][4];
  TF1 *gaus_eta[12][5][4];

  
  TH1D *Integral_phi_syst[12][4];
  TH1D *Integral_phi_Pt[12][4];
 
  TH1D *Integral_eta_Pt[12][4];
  TH1D *Integral_eta_ref_Pt[12][4];
  TGraphAsymmErrors *Integral_eta_syst[12][4];
  TGraphAsymmErrors *Integral_eta_cent[12][4];

  TGraph *Error_up_eta_pT[12][4];
  TGraph *Error_down_eta_pT[12][4];
  TH1D *Error_up_eta_cent[12][4];
  TH1D *Error_down_eta_cent[12][4];

  TLegend *l40,*l41,*l42;

  TString datalabel, in_name, centlabel, pTlabel;
  TString   leadingrefetaname[5][4];

  TCanvas *cFigure1[4],*cFigure2[4],*cFigure3[4],*cFigure4[4],*cFigure5[4],*cFigure6[4],*cFigure7;

  float raw_min,raw_max,mixed_min,mixed_max,yield_min,yield_max,result_min,result_max,  val_l,val_r,err_l,err_r;

 
  TCanvas *cintegral_eta_pT[12];
  TCanvas *cintegral_phi_pT[12];
  TCanvas *cError_up_eta_pT[12];
   
  TCanvas *cintegral_eta_cent[12];
  TCanvas *cintegral_phi_cent[12];
  TCanvas *cError_up_eta_cent[12];
  
  TLine *linePhi,*lineEta,*linePt, *lineCent;

  TLatex *tex24phi,*tex24eta;


  float value, error, evalpt, value2, dx_eta, dx_phi;

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


  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
 
  for(int g = 0; g<6; g++){
    
   
    TString integral_eta_pT_name = "integral_eta_pT";
    integral_eta_pT_name+=g;

    cintegral_eta_pT[g] = new TCanvas(integral_eta_pT_name,"",10,10,1500,500);
    cintegral_eta_pT[g]->Divide(4,1,0.,0.);


    TString integral_phi_pT_name = "integral_phi_pT";
    integral_phi_pT_name+=g;

    cintegral_phi_pT[g] = new TCanvas(integral_phi_pT_name,"",10,10,1500,500);
    cintegral_phi_pT[g]->Divide(4,1,0.,0.);

    TString integral_eta_cent_name = "integral_eta_cent";
    integral_eta_cent_name+=g;

    cintegral_eta_cent[g] = new TCanvas(integral_eta_cent_name,"",10,10,1600,450);
    cintegral_eta_cent[g]->Divide(4,1,0.0001,0.001);

    TString integral_phi_cent_name = "integral_phi_cent";
    integral_phi_cent_name+=g;

    cintegral_phi_cent[g] = new TCanvas(integral_phi_cent_name,"",10,10,1500,500);
    cintegral_phi_cent[g]->Divide(4,1,0.,0.);



    if(g<6){
      TString Error_up_eta_pT_name = "Error_up_eta_pT";
      Error_up_eta_pT_name+=g;

      cError_up_eta_pT[g] = new TCanvas(Error_up_eta_pT_name,"",10,10,1500,400);
      cError_up_eta_pT[g]->Divide(4,1,0.,0.);

     
      TString Error_up_eta_cent_name = "Error_up_eta_cent";
      Error_up_eta_cent_name+=g;

      cError_up_eta_cent[g] = new TCanvas(Error_up_eta_cent_name,"",10,10,1500,400);
      cError_up_eta_cent[g]->Divide(4,1,0.,0.);

          
    }


    //Open files
  
    switch(g){
    case 0:
      fin[g] = new TFile("../study_yield/Inclusive_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/Inclusive_Data_AllPlots.root", "READ");
      //  fin[g] = new TFile("../study_yield/Inclusive_Data_NoSpillOver_AllPlots.root", "READ");
      datalabel = "Inclusive";     break;
    case 1:
      fin[g] = new TFile("../study_yield/Inclusive_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/Inclusive_Data_AllPlots.root", "READ");
      datalabel = "Inclusive";     break;
    case 2:
      fin[g] = new TFile("../study_yield/SubLeading_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/SubLeading_Data_AllPlots.root", "READ");
      //     fin[g] = new TFile("../study_yield/SubLeading_Data_NoSpillOver_AllPlots.root", "READ");
      datalabel = "SubLeading";    break;
    case 3:
      fin[g] = new TFile("../study_yield/SubLeading_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/SubLeading_Data_AllPlots.root", "READ");
      datalabel = "SubLeading";    break;
    case 4:
      fin[g] = new TFile("../study_yield/Leading_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/Leading_Data_AllPlots.root", "READ");
      //    fin[g] = new TFile("../study_yield/Leading_Data_NoSpillOver_AllPlots.root", "READ");
      datalabel = "Leading";       break;
    case 5:
      fin[g] = new TFile("../study_yield/Leading_Data_AllPlots.root", "READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/Leading_Data_AllPlots.root", "READ");
      datalabel = "Leading";       break;
   
    default:
      break;
    }

    cout<<g<<endl;

    for(int i=0; i<4; i++){
  
      for (int j=0; j<4; j++){

	if(g<6){	in_name =	make_name("Result_",g,i,j,0,centlabel,pTlabel);}
	if(g>5){
	  in_name =	make_name("Result_",g,i,j,0,centlabel,pTlabel);
	}

	//-------------------------------------
	//  Get and assign all histograms 
	//-------------------------------------


	TString newchecknamePhi_rebin = in_name;
	TString newchecknameEta_rebin = in_name;


	newchecknameEta_rebin.ReplaceAll("Result", "New_check_Eta");
	newchecknameEta_rebin += "rebin";
	newchecknamePhi_rebin.ReplaceAll("Result", "New_check_Phi");
	newchecknamePhi_rebin += "rebin";
	
	check_new_eta_rebin[g][i][j] = (TH1D*)fin[g]->Get(newchecknameEta_rebin)->Clone((TString)(newchecknameEta_rebin+"new"));
	check_new_eta_rebin[g][i][j]->SetAxisRange(-1.35,1.35,"x");

	
	check_new_phi_rebin[g][i][j] = (TH1D*)fin[g]->Get(newchecknamePhi_rebin)->Clone((TString)(newchecknamePhi_rebin+"new"));
	check_new_phi_rebin[g][i][j]->SetAxisRange(-1.35,1.35,"x");


	if(g==0||g==2||g==4){
	  TString newchecknameEtasyst = newchecknameEta_rebin;
	  newchecknameEtasyst.ReplaceAll("New_check","Syst_error");
	  newchecknameEtasyst.ReplaceAll("rebin","");

	  cout<<newchecknameEtasyst<<endl;

	  check_new_eta_syst[g][i][j] = (TH1D*)fin[g]->Get(newchecknameEtasyst)->Clone((TString)(newchecknameEtasyst+"new"));
	  check_new_eta_syst[g][i][j]->SetAxisRange(-1.4,1.41,"x");

	  TString newchecknamePhisyst = newchecknamePhi_rebin;
	  newchecknamePhisyst.ReplaceAll("New_check","Syst_error");
	  newchecknamePhisyst.ReplaceAll("rebin","");
	
	  check_new_phi_syst[g][i][j] = (TH1D*)fin[g]->Get(newchecknamePhisyst)->Clone((TString)(newchecknamePhisyst+"new"));
	  check_new_phi_syst[g][i][j]->SetAxisRange(-1.4,1.41,"x");

	  TString newchecknameEta_ref = in_name;
	  TString newchecknamePhi_ref = in_name;
	  newchecknameEta_ref.ReplaceAll("Result", "New_check_Eta");

	  newchecknameEta_ref += "ref";

	  newchecknamePhi_ref.ReplaceAll("Result", "New_check_Phi");
	  newchecknamePhi_ref += "ref";


	  check_new_eta_ref[g][i][j] = (TH1D*)fin_ref[g]->Get(newchecknameEtasyst)->Clone(newchecknameEta_ref);
	  check_new_phi_ref[g][i][j] = (TH1D*)fin_ref[g]->Get(newchecknamePhisyst)->Clone(newchecknamePhi_ref);

	  check_new_eta_ref[g][i][j]->SetMarkerSize(1);
	  check_new_phi_ref[g][i][j]->SetMarkerSize(1);
	  check_new_phi_ref[g][i][j]->SetMarkerStyle(20);
	  check_new_eta_ref[g][i][j]->SetMarkerStyle(20);


	}else{

	  TString newchecknameEta_ref = in_name;
	  TString newchecknamePhi_ref = in_name;
	  newchecknameEta_ref.ReplaceAll("Result", "New_check_Eta");

	  newchecknameEta_ref += "ref";

	  newchecknamePhi_ref.ReplaceAll("Result", "New_check_Phi");
	  newchecknamePhi_ref += "ref";

	  check_new_eta_ref[g][i][j] = (TH1D*)fin_ref[g]->Get(newchecknameEta_rebin)->Clone(newchecknameEta_ref);
	  check_new_phi_ref[g][i][j] = (TH1D*)fin_ref[g]->Get(newchecknamePhi_rebin)->Clone(newchecknamePhi_ref);
	}

	check_new_eta_ref[g][i][j]->SetMarkerColor(kRed);
	check_new_eta_ref[g][i][j]->SetLineColor(kRed);





	check_new_phi_ref[g][i][j]->SetMarkerColor(kRed);
	check_new_phi_ref[g][i][j]->SetLineColor(kRed);


	cout<<"here"<<endl;


	nbins = check_new_eta_rebin[g][i][j]->GetNbinsX()+1;
	  
	for(int k = 1; k<nbins/2+1;k++){
	    
	  val_l = check_new_eta_rebin[g][i][j]->GetBinContent(k);
	  err_l  = check_new_eta_rebin[g][i][j]->GetBinError(k);
	    
	  val_r = check_new_eta_rebin[g][i][j]->GetBinContent(nbins-k);
	  err_r = check_new_eta_rebin[g][i][j]->GetBinError(nbins-k);
	    
	  cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<TMath::Sqrt(err_l*err_l+err_r*err_r)/2<<endl;

	  check_new_eta_rebin[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	  check_new_eta_rebin[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	  if(g==0||g==2||g==4){
	    check_new_eta_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_eta_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	  }
	  if(k<nbins/2){
	    check_new_eta_rebin[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    check_new_eta_rebin[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	  }
	}




	nbins = check_new_phi_rebin[g][i][j]->GetNbinsX()+1;
	  
	for(int k = 1; k<nbins/2+1;k++){
	    
	  val_l = check_new_phi_rebin[g][i][j]->GetBinContent(k);
	  err_l  = check_new_phi_rebin[g][i][j]->GetBinError(k);
	    
	  val_r = check_new_phi_rebin[g][i][j]->GetBinContent(nbins-k);
	  err_r = check_new_phi_rebin[g][i][j]->GetBinError(nbins-k);
	    
	  cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<TMath::Sqrt(err_l*err_l+err_r*err_r)/2<<endl;

	  check_new_phi_rebin[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	  check_new_phi_rebin[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	  if(g==0||g==2||g==4){
	    check_new_phi_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    check_new_phi_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	  }
	  if(k<nbins/2){
	    check_new_phi_rebin[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    check_new_phi_rebin[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	  }
	}


	check_new_phi_rebin[g][i][j]->GetYaxis()->SetTitle("Y #equiv #frac{1}{N_{jet}} #frac{d^{2}N}{d#Delta#phidp_{T}}   (GeV/c)^{-1}");
	check_new_eta_rebin[g][i][j]->GetYaxis()->SetTitle("Y #equiv #frac{1}{N_{jet}} #frac{d^{2}N}{d#Delta#etadp_{T}}   (GeV/c)^{-1}");


	
	if(g==1||g==3||g==5){
	  
	  TString PbPbppname_eta =in_name;
	  PbPbppname_eta.ReplaceAll("Result_pp","PbPb_minus_pp_eta");

	  TString PbPbppname_eta_syst =in_name;
	  PbPbppname_eta_syst.ReplaceAll("Result_pp","PbPb_minus_pp_eta_syst");



	  PbPb_pp_eta[g][i][j] = (TH1D*)fin[g-1]->Get(PbPbppname_eta)->Clone(PbPbppname_eta);
	  PbPb_pp_eta_syst[g][i][j] = (TH1D*)fin[g-1]->Get(PbPbppname_eta_syst)->Clone(PbPbppname_eta_syst);
	


	 


	  TString PbPbppname_phi =in_name;
	  PbPbppname_phi.ReplaceAll("Result_pp","PbPb_minus_pp_phi");

	  TString PbPbppname_phi_syst =in_name;
	  PbPbppname_phi_syst.ReplaceAll("Result_pp","PbPb_minus_pp_phi_syst");

	  PbPb_pp_phi[g][i][j] = (TH1D*)fin[g-1]->Get(PbPbppname_phi)->Clone(PbPbppname_phi);
	  PbPb_pp_phi_syst[g][i][j] = (TH1D*)fin[g-1]->Get(PbPbppname_phi_syst)->Clone(PbPbppname_phi_syst);


	  PbPb_pp_eta_ref[g][i][j] = (TH1D*)fin_ref[g-1]->Get(PbPbppname_eta)->Clone((TString)(PbPbppname_eta+"Ref"));

	  PbPb_pp_phi_ref[g][i][j] = (TH1D*)fin_ref[g-1]->Get(PbPbppname_phi)->Clone((TString)(PbPbppname_phi+"Ref"));
	  PbPb_pp_eta_ref[g][i][j]->SetLineColor(kRed);
	  PbPb_pp_eta_ref[g][i][j]->SetMarkerColor(kRed);

	  PbPb_pp_phi_ref[g][i][j]->SetLineColor(kRed);
	  PbPb_pp_phi_ref[g][i][j]->SetMarkerColor(kRed);

	 
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
	   
	  }

	  nbins = PbPb_pp_phi[g][i][j]->GetNbinsX()+1;
	
	  for(int k = 1; k<nbins/2+1;k++){
	    val_l = PbPb_pp_phi[g][i][j]->GetBinContent(k);
	 
	    err_l  = PbPb_pp_phi[g][i][j]->GetBinError(k);
	    
	    val_r = PbPb_pp_phi[g][i][j]->GetBinContent(nbins-k);
	    err_r = PbPb_pp_phi[g][i][j]->GetBinError(nbins-k);
	    
	    cout<<k<<" "<<nbins-k<<" "<<val_l<<" "<<val_r<<" "<<err_l<<" "<<err_r<<" "<<TMath::Sqrt(err_l*err_l+err_r*err_r)/2<<endl;

	    PbPb_pp_phi[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    PbPb_pp_phi[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);
	    PbPb_pp_phi_syst[g][i][j]->SetBinContent(k,(val_l+val_r)/2);
	    PbPb_pp_phi_syst[g][i][j]->SetBinContent(nbins-k,(val_l+val_r)/2);

	    if(k<nbins/2){
	      PbPb_pp_phi[g][i][j]->SetBinError(k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	      PbPb_pp_phi[g][i][j]->SetBinError(nbins-k,TMath::Sqrt(err_l*err_l+err_r*err_r)/2);
	    }
	   
	  }	

	  PbPb_pp_eta[g][i][j]->SetAxisRange(-1.4,1.41,"x");
	  PbPb_pp_eta_syst[g][i][j]->SetAxisRange(-1.4,1.41,"x");
	  PbPb_pp_phi[g][i][j]->SetAxisRange(-1.4,1.41,"x");
	  PbPb_pp_phi_syst[g][i][j]->SetAxisRange(-1.4,1.41,"x");

	
	  if(i==0){ //OBVIOUSLY only one set of integrals per centrality class!
	
	
	    TString IntegratedYieldname_eta = in_name;
	    IntegratedYieldname_eta.ReplaceAll("Result_pp","Integrated_Yield_Eta");
	    IntegratedYieldname_eta.ReplaceAll("Pt100_Pt300_","");
	    IntegratedYieldname_eta.ReplaceAll("_TrkPt1_TrkPt2","");
 

	    TString IntegratedYieldname_eta_ref = in_name;
	    IntegratedYieldname_eta_ref.ReplaceAll("Result_pp","Integrated_Yield_Eta_Ref");
	    IntegratedYieldname_eta_ref.ReplaceAll("Pt100_Pt300_","");
	    IntegratedYieldname_eta_ref.ReplaceAll("_TrkPt1_TrkPt2","");
 
	  
	 
	    Integral_eta_ref_Pt[g][j] = (TH1D*)fin_ref[g-1]->Get(IntegratedYieldname_eta)->Clone(IntegratedYieldname_eta_ref);
	    Integral_eta_ref_Pt[g][j]->SetName(IntegratedYieldname_eta_ref);
	    Integral_eta_ref_Pt[g][j]->SetAxisRange(1.01,7.8,"x");



	  
	 
	    Integral_eta_Pt[g][j] = (TH1D*)fin[g-1]->Get(IntegratedYieldname_eta)->Clone(IntegratedYieldname_eta);
	    Integral_eta_Pt[g][j]->SetAxisRange(1.01,7.8,"x");

 
	  
	    TString IntegratedYieldname_eta_syst = IntegratedYieldname_eta;
	    IntegratedYieldname_eta_syst.ReplaceAll("Eta","Eta_Syst");
	    Integral_eta_syst[g][j] = (TGraphAsymmErrors*)fin[g-1]->Get(IntegratedYieldname_eta_syst)->Clone(IntegratedYieldname_eta_syst);
	 
	    TString IntegratedYieldname_phi = IntegratedYieldname_eta;
	    IntegratedYieldname_phi.ReplaceAll("Eta","Phi");
	    Integral_phi_Pt[g][j] = (TH1D*)fin[g-1]->Get(IntegratedYieldname_phi)->Clone(IntegratedYieldname_phi);
	    Integral_phi_Pt[g][j]->SetAxisRange(1.01,7.8,"x");	  
 


	    TString IntegratedYieldname_phi_syst = IntegratedYieldname_eta_syst;
	    IntegratedYieldname_phi_syst.ReplaceAll("Eta","Phi");
	    Integral_phi_syst[g][j] = (TH1D*)fin[g-1]->Get(IntegratedYieldname_phi_syst)->Clone(IntegratedYieldname_phi_syst);

	
	
	
	    TString ErrorUpEtaPt_name = IntegratedYieldname_eta;
	    ErrorUpEtaPt_name.ReplaceAll("Integrated_Yield_Eta","Integral_Error_Up");
	    Error_up_eta_pT[g][j] = (TGraph*)fin[g-1]->Get(ErrorUpEtaPt_name)->Clone(ErrorUpEtaPt_name);

	    TString ErrorDownEtaPt_name = ErrorUpEtaPt_name;
	    ErrorDownEtaPt_name.ReplaceAll("Up","Down");
	    Error_down_eta_pT[g][j] = (TGraph*)fin[g-1]->Get(ErrorDownEtaPt_name)->Clone(ErrorDownEtaPt_name);
		
	  } //end no pT classes for integrals
	 
	} //end integrals only for  g==1,3,5

	if(g==6||g==8||g==10){
	  
	  TString HminusPnameEta = in_name;
	  HminusPnameEta.ReplaceAll("Pt100_Pt300_","");
	  HminusPnameEta.ReplaceAll("Result","Raw_HYD_minus_PYTH_Eta");
	     
	  TString gaus_eta_name = in_name;
	  gaus_eta_name.ReplaceAll("Pt100_Pt300_","");
	  gaus_eta_name.ReplaceAll("Result","GausFit_Eta");

	  TString HminusPnamePhi = in_name;
	  HminusPnamePhi.ReplaceAll("Pt100_Pt300_","");
	  HminusPnamePhi.ReplaceAll("Result","Raw_HYD_minus_PYTH_Phi");
	     
	  TString gaus_phi_name = in_name;
	  gaus_phi_name.ReplaceAll("Pt100_Pt300_","");
	  gaus_phi_name.ReplaceAll("Result","GausFit_Phi");
	  
	  HYDJET_PYTHIA_eta[g][i][j] = (TH1D*)fin[g]->Get(HminusPnameEta)->Clone(HminusPnameEta);
	  HYDJET_PYTHIA_phi[g][i][j] = (TH1D*)fin[g]->Get(HminusPnamePhi)->Clone(HminusPnamePhi);
	  gaus_eta[g][i][j] = (TF1*)fin[g]->Get(gaus_eta_name)->Clone(gaus_eta_name);
	  gaus_phi[g][i][j] = (TF1*)fin[g]->Get(gaus_phi_name)->Clone(gaus_phi_name);

	} //end 


      } //end j
 
    } //end i
  
 
    //*****************************
    //   Plot integrals
    //****************************
     
    if(g==1||g==3||g==5){
    
      //-----------------------------------
      //    pT integral plotting loop
      //----------------------------------


     

      for(int j=0;j<4;j++){

	cintegral_eta_pT[g]->cd(j+1);
      

	Integral_eta_Pt[g][j]->SetMarkerColor(1);
	Integral_eta_Pt[g][j]->SetLineColor(1);
	Integral_eta_syst[g][j]->SetMarkerSize(0);
	Integral_eta_Pt[g][j]->SetMarkerSize(2);
	Integral_eta_ref_Pt[g][j]->SetMarkerSize(2);

	switch(g){
	case 1: 
	  Integral_eta_Pt[g][j]->SetMarkerStyle(10);
	  Integral_eta_ref_Pt[g][j]->SetMarkerStyle(10);
	  Integral_eta_syst[g][j]->SetFillColor(90);
	  Error_up_eta_pT[g][j]->SetMarkerStyle(10);
	  Error_up_eta_pT[g][j]->SetFillColor(90);
	  Error_down_eta_pT[g][j]->SetMarkerStyle(10);
	  Error_down_eta_pT[g][j]->SetFillColor(90);
	  break;
	case 3:
	  Integral_eta_Pt[g][j]->SetMarkerStyle(34);
	  Integral_eta_ref_Pt[g][j]->SetMarkerStyle(34);
	  Integral_eta_syst[g][j]->SetFillColor(30);
	  Error_up_eta_pT[g][j]->SetMarkerStyle(34);
	  Error_up_eta_pT[g][j]->SetFillColor(30);
	  Error_down_eta_pT[g][j]->SetMarkerStyle(34);
	  Error_down_eta_pT[g][j]->SetFillColor(30);
	  break;
	case 5: 
	  Integral_eta_Pt[g][j]->SetMarkerStyle(21);
	  Integral_eta_ref_Pt[g][j]->SetMarkerStyle(21);
	  Integral_eta_syst[g][j]->SetFillColor(kOrange-2);
	  Error_up_eta_pT[g][j]->SetMarkerStyle(21);
	  Error_up_eta_pT[g][j]->SetFillColor(kOrange-2);
	  Error_down_eta_pT[g][j]->SetMarkerStyle(21);
	  Error_down_eta_pT[g][j]->SetFillColor(kOrange-2);
	  break;
	default:
	  Integral_eta_Pt[g][j]->SetMarkerStyle(10);
	  break;
	}

	Error_down_eta_pT[g][j]->SetMarkerColor(kRed);


	if(g==1){Integral_eta_Pt[g][j]->SetMaximum(4.2);
	  Integral_eta_Pt[g][j]->SetMinimum(-.45);}
	if(g==3){Integral_eta_Pt[g][j]->SetMaximum(10.2);
	  Integral_eta_Pt[g][j]->SetMinimum(-1.);}
	if(g==5){Integral_eta_Pt[g][j]->SetMaximum(4.2);
	  Integral_eta_Pt[g][j]->SetMinimum(-.45);}


	    
	Integral_eta_Pt[g][j]->Draw("p");
	 
	
	Integral_eta_Pt[g][j]->GetYaxis()->SetLabelSize(tstitle);
	   
	Integral_eta_Pt[g][j]->GetXaxis()->SetLabelSize(tstitle);
	Integral_eta_Pt[g][j]->GetXaxis()->SetTitle("Track p_{T} (GeV/c)");
	Integral_eta_Pt[g][j]->GetXaxis()->SetTitleSize(tstitle);
	Integral_eta_Pt[g][j]->GetXaxis()->SetTitleOffset(xoffset);
	Integral_eta_Pt[g][j]->GetYaxis()->SetTitle("(dN/dp_{T})_{PbPb}- (dN/dp_{T})_{pp}   (GeV/c)^{-1}");
	Integral_eta_Pt[g][j]->GetYaxis()->SetTitleOffset(yoffset);
	Integral_eta_Pt[g][j]->GetYaxis()->SetTitleSize(tstitle2);

 
	Integral_eta_Pt[g][j]->GetXaxis()->SetRangeUser(1.000001,7.7);
	Integral_eta_Pt[g][j]->GetXaxis()->CenterTitle();
	Integral_eta_Pt[g][j]->GetYaxis()->CenterTitle();
	   
	if(j>0){
	  Integral_eta_Pt[g][j]->GetYaxis()->SetTitleSize(0.0);
	  Integral_eta_Pt[g][j]->GetYaxis()->SetLabelSize(0.0);
	}



	if(g==5){
	   
	  Integral_eta_syst[g][j]->Draw("same 2");
	  Integral_eta_Pt[g][j]->Draw("same p");
	}
	    

	if(g==3){

	  Integral_eta_syst[g][j]->Draw("same 2");
	  Integral_eta_Pt[g][j]->Draw("same p");
	   
	}

	if(g==1){
	  Integral_eta_syst[g][j]->Draw("same 2");
	  Integral_eta_Pt[g][j]->Draw("same p");

	}
	 
	if(g>5){
	  Integral_eta_syst[g][j]->Draw("same2");

	}
	if(g==1||g==3||g==5){
	  Integral_eta_ref_Pt[g][j]->SetLineColor(kRed);
	  Integral_eta_ref_Pt[g][j]->SetMarkerColor(kRed);
	  Integral_eta_ref_Pt[g][j]->Draw("same p");
	}
	if(j==0){ 
	  l40 = new TLegend(textalign2-0.05,texty2-.1,0.6,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextFont(63);
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){l40->AddEntry(Integral_eta_Pt[g][j],"Inclusive Jets","lp");
	  }
	  if(g==3){l40->AddEntry(Integral_eta_Pt[g][j],"SubLeading Jets","lp");
	  }
	  if(g==5){l40->AddEntry(Integral_eta_Pt[g][j],"Leading Jets","lp");
	  }
	  //	  l40->AddEntry( Integral_eta_ref_Pt[g][j],"PAS Reference", "lp");
	  l40->Draw("same");
	}
	
	
	    
	linePt = new TLine(1.,0,8.,0);
	linePt->SetLineStyle(2);
	linePt->SetLineWidth(2);
	linePt->Draw("same");

	drawlabels_int_pt(g,j);



	cintegral_phi_pT[g]->cd(j+1);
    
	Integral_phi_Pt[g][j]->SetMarkerColor(1);
	Integral_phi_Pt[g][j]->SetLineColor(1);
	Integral_phi_syst[g][j]->SetMarkerSize(0);
	Integral_phi_Pt[g][j]->SetMarkerSize(2);
  
       
	  
	    
	Integral_phi_Pt[g][j]->SetMarkerSize(2);
	Integral_phi_Pt[g][j]->SetMarkerColor(1);
	Integral_phi_Pt[g][j]->SetLineColor(1);
	   

	switch(g){
	case 1: 
	  Integral_phi_Pt[g][j]->SetMarkerStyle(10);
	  Integral_phi_syst[g][j]->SetFillColor(90);
	  break;
	case 3:
	  Integral_phi_Pt[g][j]->SetMarkerStyle(34);
	  Integral_phi_syst[g][j]->SetFillColor(30);
	  break;
	case 5: 
	  Integral_phi_Pt[g][j]->SetMarkerStyle(21);
	  Integral_phi_syst[g][j]->SetFillColor(kOrange-2);
	  break;
	default:
	  Integral_phi_Pt[g][j]->SetMarkerStyle(10);
	  break;
	}


	Integral_phi_Pt[g][j]->SetMinimum(-1.);
	Integral_phi_Pt[g][j]->SetMaximum(10.2);
	 

	Integral_phi_Pt[g][j]->SetMarkerColor(kRed);
	Integral_phi_Pt[g][j]->SetLineColor(kRed);
	Integral_phi_Pt[g][j]->SetMarkerStyle(29);
	Integral_phi_Pt[g][j]->SetMarkerSize(5);
	Integral_phi_Pt[g][j]->Draw("p");
	  


	Integral_phi_Pt[g][j]->GetYaxis()->SetLabelSize(tstitle2);
	   
	Integral_phi_Pt[g][j]->GetXaxis()->SetLabelSize(tstitle);
	Integral_phi_Pt[g][j]->GetXaxis()->SetTitle("Track p_{T} (GeV/c)");
	Integral_phi_Pt[g][j]->GetXaxis()->SetTitleSize(tstitle);
	Integral_phi_Pt[g][j]->GetXaxis()->SetTitleOffset(xoffset);
	Integral_phi_Pt[g][j]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV/c)^{-1}");
	Integral_phi_Pt[g][j]->GetYaxis()->SetTitleOffset(yoffset);
	Integral_phi_Pt[g][j]->GetYaxis()->SetTitleSize(tstitle2);

 
	Integral_phi_Pt[g][j]->GetXaxis()->SetRangeUser(1.000001,7.7);
	Integral_phi_Pt[g][j]->GetXaxis()->CenterTitle();
	Integral_phi_Pt[g][j]->GetYaxis()->CenterTitle();
	   
	if(j>0){
	  Integral_phi_Pt[g][j]->GetYaxis()->SetTitleSize(0.0);
	  Integral_phi_Pt[g][j]->GetYaxis()->SetLabelSize(0.0);
	}
	  
	Integral_eta_Pt[g][j]->Draw("same");


	if(j==0){
	  l40 = new TLegend(textalign2-0.05,texty2-.1,0.6,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextSizePixels(tspixels);
	  l40->SetTextFont(43);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){l40->AddEntry(Integral_eta_Pt[g][j],"Inclusive Jets Integrated #Delta#eta","lp");
	    l40->AddEntry(Integral_phi_Pt[g][j],"Inclusive Jets Integrated #Delta#phi","lp");}
	  if(g==3){l40->AddEntry(Integral_eta_Pt[g][j],"SubLeading Jets Int. #Delta#eta","lp");
	    l40->AddEntry(Integral_phi_Pt[g][j],"SubLeading Jets Int. #Delta#phi","lp");}
	  if(g==5){l40->AddEntry(Integral_eta_Pt[g][j],"Leading Jets Integrated #Delta#eta","lp");
	    l40->AddEntry(Integral_phi_Pt[g][j],"Leading Jets Integrated #Delta#phi","lp");}
	  l40->Draw("same");
	}
	



	//	if(j==0){l40->Draw("same");}
	    
	linePt->Draw("same");

	  
	drawlabels_int_pt(g,j);



	// Draw error up evolution plot

	cError_up_eta_pT[g]->cd(j+1);
	
	gStyle->SetOptTitle(0);

	Error_up_eta_pT[g][j]->SetMinimum(-0.1);
	Error_up_eta_pT[g][j]->SetMaximum(4.);
	
	
	    
	Error_up_eta_pT[g][j]->Draw("p");
	


	Error_up_eta_pT[g][j]->GetYaxis()->SetLabelSize(tstitle);
	   
	Error_up_eta_pT[g][j]->GetXaxis()->SetLabelSize(tstitle);
	Error_up_eta_pT[g][j]->GetXaxis()->SetTitle("Track p_{T} (GeV/c)");
	Error_up_eta_pT[g][j]->GetXaxis()->SetTitleSize(tstitle);
	Error_up_eta_pT[g][j]->GetXaxis()->SetTitleOffset(xoffset);
	Error_up_eta_pT[g][j]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV/c)^{-1}");
	Error_up_eta_pT[g][j]->GetYaxis()->SetTitleOffset(yoffset);
	Error_up_eta_pT[g][j]->GetYaxis()->SetTitleSize(tstitle2);

 
	Error_up_eta_pT[g][j]->GetXaxis()->SetRangeUser(1.000001,7.7);
	Error_up_eta_pT[g][j]->GetXaxis()->CenterTitle();
	Error_up_eta_pT[g][j]->GetYaxis()->CenterTitle();
	   
	if(j>0){
	  Error_up_eta_pT[g][j]->GetYaxis()->SetTitleSize(0.0);
	  Error_up_eta_pT[g][j]->GetYaxis()->SetLabelSize(0.0);
	}
	
	Error_up_eta_pT[g][j]->SetMarkerSize(2);
	Error_up_eta_pT[g][j]->Draw();
		
	drawlabels_int_pt(g,j);

	if(j==0){ 
	  l40 = new TLegend(textalign2-0.05,texty2-.1,0.6,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){l40->AddEntry(Error_up_eta_pT[g][j],"Inclusive ErrorUp","lp");
	    l40->AddEntry(Error_down_eta_pT[g][j],"Inclusive ErrorDown","lp");}
	  if(g==3){l40->AddEntry(Error_up_eta_pT[g][j],"SubLeading ErrorDown","lp");
	    l40->AddEntry(Error_down_eta_pT[g][j],"SubLeading ErrorDown","lp");}
	  if(g==5){l40->AddEntry(Error_up_eta_pT[g][j],"Leading ErrorUp","lp");
	    l40->AddEntry(Error_down_eta_pT[g][j],"Leading ErrorDown","lp");}
	  l40->Draw("same");
	}
	
	linePt->Draw("same");
	
	Error_down_eta_pT[g][j]->SetLineColor(kRed);
	Error_down_eta_pT[g][j]->SetMarkerSize(2);
	Error_down_eta_pT[g][j]->Draw("same pl");
	


      } //close j


      TString IntegralSaveName_eta = in_name;
      IntegralSaveName_eta.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      IntegralSaveName_eta.ReplaceAll("Result","Integral_Yield");
      IntegralSaveName_eta.ReplaceAll("pp","PbPb_minus_pp_eta");
      IntegralSaveName_eta +=".pdf";
      cintegral_eta_pT[g]->SaveAs(IntegralSaveName_eta);


      TString IntegralSaveName_phi = in_name;
      IntegralSaveName_phi.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      IntegralSaveName_phi.ReplaceAll("Result","Integral_Yield");
      IntegralSaveName_phi.ReplaceAll("pp","PbPb_minus_pp_phi");
      IntegralSaveName_phi +=".pdf";
      cintegral_phi_pT[g]->SaveAs(IntegralSaveName_phi);
    

      TString ErrorUpSaveNamePt = in_name;
      ErrorUpSaveNamePt.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      ErrorUpSaveNamePt.ReplaceAll("Result","Error_pT");
      ErrorUpSaveNamePt.ReplaceAll("pp","PbPb_minus_pp_phi");
      ErrorUpSaveNamePt +=".pdf";
      cError_up_eta_pT[g]->SaveAs(ErrorUpSaveNamePt);

      
    

      //---------------------------------------
      //Now for integral by centrality plots
      //--------------------------------------

      for(int i = 0; i<4; i++){

	cintegral_eta_cent[g]->cd(i+1);
	

	TString histnamecent = "integral_cent";
	histnamecent+=g;
	histnamecent+=i;

	int_cent[g][i] = new TH1D(histnamecent,"",4,xAxis);

	TString histnameblank = "blank_hist";
	histnameblank+=g;
	histnameblank+=i;

	blank[i] = new TH1D(histnameblank,"",4,xAxis);

	
	TString histnameblank2 = "blank_hist2";
	histnameblank2+=g;
	histnameblank2+=i;

	blank2[i] = new TH1D(histnameblank2,"",4,xAxis);

	
	for(int k=0; k<4; k++){

	  value = Integral_eta_Pt[g][k]->GetBinContent(i+1);
	  error = Integral_eta_syst[g][k]->GetErrorYhigh(i);

	  int_cent[g][i]->SetBinContent(k+1,value);
	  int_cent[g][i]->SetBinError(k+1,error);

	}

	

	TString ErrorUpEtaCent_name = "ErrorUpEtaCent";
	ErrorUpEtaCent_name+=g;
	ErrorUpEtaCent_name+=i;
	Error_up_eta_cent[g][i] = new TH1D(ErrorUpEtaCent_name,"",4,xAxis);

	TString ErrorDownEtaCent_name = "ErrorDownEtaCent";
	ErrorDownEtaCent_name+=g;
	ErrorDownEtaCent_name+=i;
	Error_down_eta_cent[g][i] = new TH1D(ErrorDownEtaCent_name,"",4,xAxis);

		
	for(int k=0; k<4; k++){
	  
	  evalpt = pTbin_centers.at(i);

	  value = Error_up_eta_pT[g][k]->Eval(evalpt);
	 
	  Error_up_eta_cent[g][i]->SetBinContent(k+1,value);
	  
	  value = Error_down_eta_pT[g][k]->Eval(evalpt);
	  
	  Error_down_eta_cent[g][i]->SetBinContent(k+1,value);
	  
	}


	//Markers etc. for all histograms at once
	//----------------------------------------


	int_cent[g][i]->SetMarkerSize(2);
	int_cent[g][i]->SetMarkerColor(1);
	int_cent[g][i]->SetLineColor(1);

	switch(g){
	case 1:
	  int_cent[g][i]->SetFillColor(90);
	  int_cent[g][i]->SetMarkerStyle(10);
	  Error_up_eta_cent[g][i]->SetMarkerStyle(10);
	  Error_down_eta_cent[g][i]->SetMarkerStyle(10);
	  break;     
	case 3:
	  int_cent[g][i]->SetFillColor(30);
	  int_cent[g][i]->SetMarkerStyle(34);
	  Error_up_eta_cent[g][i]->SetMarkerStyle(34);
	  Error_down_eta_cent[g][i]->SetMarkerStyle(34);
	  break;  
   	case 5:
	  int_cent[g][i]->SetFillColor(kOrange-2);
	  int_cent[g][i]->SetMarkerStyle(21);
	  Error_up_eta_cent[g][i]->SetMarkerStyle(21);
	  Error_down_eta_cent[g][i]->SetMarkerStyle(21);
	  break;     
	}

	Error_down_eta_cent[g][i]->SetMarkerColor(kRed);
	
	//Plot aesthetics for every canvas.
	//----------------------------------
	
	blank[i]->SetMinimum(-1.);
	blank[i]->SetMaximum(10.2);
	blank[i]->GetXaxis()->SetTitle("Centrality (%)");
	blank[i]->GetXaxis()->SetTitleOffset(1.1);
	blank[i]->GetXaxis()->CenterTitle(true);
	blank[i]->GetXaxis()->SetTitleSize(ts);

	blank[i]->GetYaxis()->SetTitle("(dN/dp_{T})_{PbPb}- (dN/dp_{T})_{pp} (GeV/c)^{-1}");
	blank[i]->GetYaxis()->SetTitleSize(0.);
	blank[i]->GetYaxis()->CenterTitle(true);
	//	blank[i]->GetYaxis()->SetLabelOffset(yoffset);
	//	blank[i]->GetYaxis()->SetLabelSize(0.);
   
	blank[i]->GetYaxis()->SetTickLength(0.025);

	blank[i]->GetXaxis()->SetBinLabel(1,"50-100");
	blank[i]->GetXaxis()->SetBinLabel(2,"30-50");
	blank[i]->GetXaxis()->SetBinLabel(3,"10-30");
	blank[i]->GetXaxis()->SetBinLabel(4," 0-10");
    
	blank[i]->GetXaxis()->SetLabelSize(0.08);
	blank[i]->GetXaxis()->SetLabelOffset(0.015);
	
	blank[i]->GetXaxis()->LabelsOption("h");
	blank[i]->GetXaxis()->SetTickLength(0.0);

	switch(i){
	case 0: 
	  //  gPad->SetLeftMargin(0.2);
	  blank[i]->GetYaxis()->SetTitleSize(ts);
	  blank[i]->GetXaxis()->SetTitleOffset(1.1);
	  blank[i]->GetXaxis()->SetTitleSize(0.07);
	  blank[i]->SetLabelSize(0.95*blank[i]->GetXaxis()->GetLabelSize());
	  blank[i]->GetYaxis()->SetLabelSize(ts2);
	  break;
	case 3:
	  // gPad->SetRightMargin(0.02);
	  break;
	default:
	  break;
	}


	TLine *line1, *line2, *line3;


	switch(i){
	case 0:
	  blank[i]->SetMaximum(6.1);
	  line1 = new TLine(-50,-0.5,-50.,-0.1);
	  line2 = new TLine(-30,-0.5,-30.,-0.1);
	  line3 = new TLine(-10,-0.5,-10.,-0.1);
	  break;
	case 1:	 
	  blank[i]->SetMaximum(1.8);
	  line1 = new TLine(-50,-0.38,-50.,-0.22);
	  line2 = new TLine(-30,-0.38,-30.,-0.22);
	  line3 = new TLine(-10,-0.38,-10.,-0.22);	 
	  break;
	case 2:
	  blank[i]->SetMaximum(0.7);
	  line1 = new TLine(-50,-0.34,-50.,-0.26);
	  line2 = new TLine(-30,-0.34,-30.,-0.26);
	  line3 = new TLine(-10,-0.34,-10.,-0.26);
	  break;
	case 3:	
	  blank[i]->SetMaximum(0.4);
	  line1 = new TLine(-50,-0.33,-50.,-0.27);
	  line2 = new TLine(-30,-0.33,-30.,-0.27);
	  line3 = new TLine(-10,-0.33,-10.,-0.27);
	  break;
	}

	blank[i]->GetYaxis()->SetLabelSize(ts);

	//	blank[i]->SetMaximum(1.1);
	blank[i]->SetMinimum(-.3);

	blank[i]->GetYaxis()->SetNdivisions(408);
	blank[i]->GetYaxis()->SetTitleOffset(0.8);
	blank[i]->Draw();


	
	line1->Draw();
	line2->Draw();
	line3->Draw();
	


	
	if(i==0){ 
	  l40 = new TLegend(.17,texty2-.05,.8,texty3-0.05);
	  l40->SetName("l40");
	  l40->SetTextFont(tspixels);
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){l40->AddEntry(int_cent[g][i],"Inclusive Jets","lpfe");}
	  if(g==3){l40->AddEntry(int_cent[g][i],"Subleading Jets","lpfe");}
	  if(g==5){l40->AddEntry(int_cent[g][i],"Leading Jets","lpfe");}
	  l40->Draw("same");
	}
	drawlabels_int_cent2(g,i);
	 
	lineCent = new TLine(-100.,0.,0.,0.);
	lineCent->SetLineStyle(2);
	lineCent->SetLineWidth(2);
	lineCent->Draw("same");
	
	int_cent[g][i]->Draw("same bp0 e2");

	lineCent->Draw("same");

	gPad->RedrawAxis();

	//----------------------------
	// Draw error up evolution plot
	//----------------------------
	cError_up_eta_cent[g]->cd(i+1);
	

	gStyle->SetOptTitle(0);


	blank2[i]->GetXaxis()->SetTitle("Centrality (%)");
	blank2[i]->GetXaxis()->SetTitleOffset(1.1);
	blank2[i]->GetXaxis()->CenterTitle(true);
	blank2[i]->GetXaxis()->SetTitleSize(ts);

	blank2[i]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV/c)^{-1}");
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
	blank2[i]->GetXaxis()->SetTicks("+-");
	blank2[i]->GetXaxis()->LabelsOption("h");
	blank2[i]->GetXaxis()->SetTickLength(0.025);

	switch(i){
	case 0: 
	  gPad->SetLeftMargin(0.2);
	  blank2[i]->GetYaxis()->SetTitleSize(0.09);
	  blank2[i]->GetXaxis()->SetTitleOffset(1.2);
	  blank2[i]->GetXaxis()->SetTitleSize(0.06);
	  blank2[i]->SetLabelSize(0.9*blank2[i]->GetXaxis()->GetLabelSize());
	  blank2[i]->GetYaxis()->SetLabelSize(ts);
	  break;
	case 3:
	  gPad->SetRightMargin(0.02);
	  break;
	default:
	  break;
	}
	
	
	blank2[i]->SetMinimum(-0.1);
	blank2[i]->SetMaximum(4.);


	blank2[i]->Draw();





	Error_up_eta_cent[g][i]->SetMarkerSize(2);
	Error_up_eta_cent[g][i]->Draw("same p");
	Error_up_eta_cent[g][i]->Draw("same pl");
		
	
	Error_down_eta_cent[g][i]->SetLineColor(kRed);
	Error_down_eta_cent[g][i]->SetMarkerSize(2);
	Error_down_eta_cent[g][i]->Draw("same pl");
	


	drawlabels_int_cent2(g,i);

	if(i==0){ 
	  l40 = new TLegend(textalign2,texty2-.1,0.6,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextSizePixels(tspixels);
	  l40->SetFillColor(kWhite);
	  l40->SetLineColor(kWhite);
	  if(g==1){l40->AddEntry(Error_up_eta_pT[g][i],"Inclusive ErrorUp","lp");
	    l40->AddEntry(Error_down_eta_pT[g][i],"Inclusive ErrorDown","lp");}
	  if(g==3){l40->AddEntry(Error_up_eta_pT[g][i],"SubLeading ErrorDown","lp");
	    l40->AddEntry(Error_down_eta_pT[g][i],"SubLeading ErrorDown","lp");}
	  if(g==5){l40->AddEntry(Error_up_eta_pT[g][i],"Leading ErrorUp","lp");
	    l40->AddEntry(Error_down_eta_pT[g][i],"Leading ErrorDown","lp");}
	  l40->Draw("same");
	}	
	
	lineCent->Draw("same");
	
	
	
      } //close i



      cintegral_eta_cent[g]->cd(0);

      TLatex *canvas_title = new TLatex(0.05,0.9,"CMS Preliminary");
      canvas_title->SetTextSizePixels(tspixels);
      canvas_title->SetTextFont(63);
      canvas_title->Draw();

      TLatex *canvas_title2 = new TLatex(0.295,0.9,"PbPb 166 #mub^{-1} (2.76 TeV)                     pp 5.3 pb^{-1} (2.76 TeV)");
      canvas_title2->SetTextSizePixels(tspixels);
      canvas_title2->Draw();


      TString IntegralSaveName_cent = in_name;
      IntegralSaveName_cent.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      IntegralSaveName_cent.ReplaceAll("Result","Integral_Yield");
      IntegralSaveName_cent.ReplaceAll("pp","PbPb_minus_pp_cent");
      IntegralSaveName_cent +=".pdf";
      cintegral_eta_cent[g]->SaveAs(IntegralSaveName_cent);
      IntegralSaveName_cent.ReplaceAll(".pdf",".png");
      cintegral_eta_cent[g]->SaveAs(IntegralSaveName_cent);



      TString ErrorUpSaveName_cent = in_name;
      ErrorUpSaveName_cent.ReplaceAll("Cent0_Cent10_Pt100_Pt300_TrkPt4_TrkPt8",datalabel);
      ErrorUpSaveName_cent.ReplaceAll("Result","Error_Cent");
      ErrorUpSaveName_cent.ReplaceAll("pp","PbPb_minus_pp_cent");
      ErrorUpSaveName_cent +=".pdf";
      cError_up_eta_cent[g]->SaveAs(ErrorUpSaveName_cent);


    }// close "only integrate for pp"

  } // close g

  cout<<"Now starting on PAS Plots..."<<endl;




  //*************************************
  //   Now draw PAS plots
  //***********************************

  for(int i = 0;i <4; i++){
   
    
    TString pTrange;

    switch(i){
    case 0: 
      pTrange = "TrkPt1_TrkPt2"; 
      raw_min = 10.;
      raw_max = 50.;
      mixed_min = 10./34.;
      mixed_max = 50/34.;
      yield_min = 35.;
      yield_max = 43.;
      result_min = -.5;
      result_max = 7.5;
      break;
    case 1: 
      pTrange = "TrkPt2_TrkPt3"; 
      raw_min = 1.;
      raw_max = 10.;
      mixed_min = 1./4.5;
      mixed_max = 10./4.5;
      yield_min = 4.;
      yield_max = 10.;
      result_min = -.45;
      result_max = 5.55;
      break;
    case 2: 
      pTrange = "TrkPt3_TrkPt4"; 
      raw_min = 0.;
      raw_max = 5.;
      mixed_min = 0.;
      mixed_max = 6.;
      yield_min = 0.;
      yield_max = 5.;
      result_min = -.45;
      result_max = 4.55;
      break;
    case 3: 
      pTrange = "TrkPt4_TrkPt8"; 
      raw_min = 0.;
      raw_max = 2.;
      mixed_min = 0.;
      mixed_max = 10.;
      yield_min = -.1;
      yield_max = 1.9;
      result_min = -.1;
      result_max = 1.9;
      break;
    }
       

    TF1 *gen_gaus = new TF1("gen_gaus",
			  " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))  \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);
   
    TF1 *gen_gaus_up = new TF1("gen_gaus_up",
			       " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6])) \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6]))) ",-TMath::Pi()/2,3*TMath::Pi()/2);
   
  

    TF1 *gen_gaus_down = new TF1("gen_gaus_down",
				 " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6])) \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);
   

    TF1 *gen_gaus_level = new TF1("gen_gaus_level",
			  " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))  \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);



    TF1 *gen_gaus_bg = new TF1("gen_gaus_bg",
			  " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))  \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);

     TF1 *sub_gen_gaus = new TF1("sub_gen_gaus",
			    " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6])) \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);
   
    TF1 *sub_gen_gaus_up = new TF1("sub_gen_gaus_up",
			       " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6])) \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);
   
  

    TF1 *sub_gen_gaus_down = new TF1("sub_gen_gaus_down",
				 " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6])) \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);
   
    

    TF1 *sub_gen_gaus_bg = new TF1("sub_gen_gaus_bg",
			  " [0]*(1+2.0*[1]*TMath::Cos(1.*x)+2.0*[2]*TMath::Cos(2.*x)+2.0*[3]*TMath::Cos(3.*x))\
                           +[4]*(TMath::Exp(-pow(TMath::Abs(x+TMath::Pi())/[5],[6]))  \
                               + TMath::Exp(-pow(TMath::Abs(x-TMath::Pi())/[5],[6])))",-TMath::Pi()/2,3*TMath::Pi()/2);

    TString figure1_name = "PAS_Figure_1_";
    figure1_name+=pTrange;
    cFigure1[i] = new TCanvas(figure1_name," ",10,10,1500,1200);
    cFigure1[i]->Divide(3,2,0.001,0.001);
  
    TString figure2_name = "PAS_Figure_2_";
    figure2_name+=pTrange;
    cFigure2[i] = new TCanvas(figure2_name," ",10,10,1500,1200);
    cFigure2[i]->Divide(3,2,0.001,0.001);
  
    TString figure3_name = "PAS_Figure_3_";
    figure3_name+=pTrange;
    cFigure3[i] = new TCanvas(figure3_name," ",10,10,1550,850);
    cFigure3[i]->Divide(4,2,0.,0.);

    TString figure4_name = "PAS_Figure_4_";
    figure4_name+=pTrange;
    cFigure4[i] = new TCanvas(figure4_name," ",10,10,1550,850);
    cFigure4[i]->Divide(4,2,0.,0.);

    TString figure5_name = "PAS_Figure_5_";
    figure5_name+=pTrange;
    cFigure5[i] = new TCanvas(figure5_name," ",10,10,1550,1250);
    cFigure5[i]->Divide(4,3,0.,0.);

    TString figure6_name = "PAS_Figure_6_";
    figure6_name+=pTrange;
    cFigure6[i] = new TCanvas(figure6_name," ",10,10,1550,1250);
    cFigure6[i]->Divide(4,3,0.,0.);
 
    
  
    //--------------------
    //   PAS FIGURE 1
    //--------------------------
    cFigure1[i]->cd(1);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);


    TString temp_name = make_name("Raw_Yield_",0,i,3,0,centlabel,pTlabel);

    temp_name.ReplaceAll("PbPb_","");
       

    cout<<temp_name<<endl;

    TH2D *raw_corr = (TH2D*)fbkgsummed->Get(temp_name)->Clone(temp_name);
    
    if(i==3){raw_corr->Scale(1/4.);  }
  
     
    dx_eta = raw_corr->GetXaxis()->GetBinWidth(1);
    dx_phi = raw_corr->GetYaxis()->GetBinWidth(1);
     
    raw_corr->Rebin2D(2,4);
    raw_corr->Scale(1/dx_eta/dx_phi/8.);
    raw_corr->SetMinimum(raw_min);
    raw_corr->SetMaximum(raw_max);
    raw_corr->SetLineColor(kBlack);
    raw_corr->GetXaxis()->SetRangeUser(-3.,3.);
    raw_corr->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    raw_corr->GetXaxis()->SetTitle("#Delta#eta");
    raw_corr->GetYaxis()->SetTitle("#Delta#phi");
    raw_corr->GetXaxis()->CenterTitle();
    raw_corr->GetYaxis()->CenterTitle();
    raw_corr->GetZaxis()->CenterTitle();
    raw_corr->GetYaxis()->SetTitleSize(0.06);
    raw_corr->GetXaxis()->SetTitleSize(0.06);
    raw_corr->GetXaxis()->SetTitleOffset(1.);
    raw_corr->GetYaxis()->SetTitleOffset(1.);
    raw_corr->GetZaxis()->SetTitle("S(#Delta#eta,#Delta#phi)");
    raw_corr->GetZaxis()->SetTitleSize(0.06);
    raw_corr->GetZaxis()->SetTitleOffset(1.);
    raw_corr->GetZaxis()->SetNdivisions(205);
     
      
    raw_corr->Draw("surf1");

    drawlabels_PAS_bkg(i);

    TLatex *leading_tex = new TLatex(0.74,0.72,"Leading Jet");
    leading_tex->SetName("tex04");
    leading_tex->SetNDC();
    leading_tex->SetTextFont(63);
    leading_tex->SetTextSizePixels(tspixels);
    leading_tex->Draw();


    cFigure1[i]->cd(2);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);

    temp_name.ReplaceAll("Raw_Yield","Mixed_Event");

    TH2D *mixed_bkg = (TH2D*) fbkgsummed->Get(temp_name)->Clone(temp_name);
      
    //   if(i==3){mixed_bkg->Scale(1/4.);  }

    mixed_bkg->Rebin2D(2,4);
    mixed_bkg->Scale(1/8.);


    mixed_bkg->GetXaxis()->SetRangeUser(-3.,3.);
    mixed_bkg->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.0001);
    mixed_bkg->SetLineColor(kBlack);
    mixed_bkg->GetXaxis()->SetTitle("#Delta#eta");
    mixed_bkg->GetYaxis()->SetTitle("#Delta#phi");
    mixed_bkg->GetXaxis()->CenterTitle();
    mixed_bkg->GetYaxis()->CenterTitle();
    mixed_bkg->GetZaxis()->CenterTitle();
    mixed_bkg->GetYaxis()->SetTitleSize(0.06);
    mixed_bkg->GetXaxis()->SetTitleSize(0.06);
    mixed_bkg->GetZaxis()->SetTitle("ME(#Delta#eta,#Delta#phi)");
    mixed_bkg->GetXaxis()->SetTitleOffset(1.);
    mixed_bkg->GetYaxis()->SetTitleOffset(1.);
    mixed_bkg->GetZaxis()->SetTitleOffset(1.);
    mixed_bkg->SetMinimum(mixed_min);
    mixed_bkg->SetMaximum(mixed_max);
    
    mixed_bkg->GetXaxis()->SetTitle("#Delta#eta");
    mixed_bkg->GetYaxis()->SetTitle("#Delta#phi");
    mixed_bkg->GetZaxis()->SetTitle("ME(#Delta#eta,#Delta#phi)");
    mixed_bkg->GetZaxis()->SetTitleSize(0.06);
    mixed_bkg->GetZaxis()->SetNdivisions(205);
     
    
    mixed_bkg->Draw("surf1");

    drawlabels_PAS_bkg(i);

    leading_tex->Draw();

    cFigure1[i]->cd(3);
 
    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);

    temp_name.ReplaceAll("Mixed_Event","Yield_PbPb");
    TH2D *raw_yield = (TH2D*) fbkgsummed->Get(temp_name)->Clone(temp_name);

    if(i==3){raw_yield->Scale(1/4.);}
    
    raw_yield->Rebin2D(2,4);
    raw_yield->Scale(1/dx_eta/dx_phi/8.);
      
    raw_yield->GetXaxis()->SetRangeUser(-3.,3.);
    raw_yield->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    raw_yield->GetXaxis()->SetTitle("#Delta#eta");
    raw_yield->GetYaxis()->SetTitle("#Delta#phi");
    raw_yield->GetXaxis()->CenterTitle();
    raw_yield->GetYaxis()->CenterTitle();
    raw_yield->GetZaxis()->CenterTitle();
    raw_yield->GetYaxis()->SetTitleSize(0.06);
    raw_yield->GetXaxis()->SetTitleSize(0.06);
    raw_yield->GetZaxis()->SetTitle("S(#Delta#eta,#Delta#phi) / ME(#Delta#eta,#Delta#phi)");
    raw_yield->GetZaxis()->SetTitleSize(0.06);
    raw_yield->GetXaxis()->SetTitleOffset(1.);
    raw_yield->GetYaxis()->SetTitleOffset(1.);
    raw_yield->GetZaxis()->SetTitleOffset(1.);
     
    raw_yield->SetMinimum(yield_min);
    raw_yield->SetMaximum(yield_max);
    raw_yield->GetZaxis()->SetNdivisions(205);
    raw_yield->SetLineColor(kBlack);
    raw_yield->Draw("surf1");
      
    drawlabels_PAS_bkg(i);

    leading_tex->Draw();

    cFigure1[i]->cd(4);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);


    temp_name = make_name("Raw_Yield_",0,i,3,0,centlabel,pTlabel);

    temp_name.ReplaceAll("PbPb_","");
       

    cout<<temp_name<<endl;

    TH2D *raw_corr_sub = (TH2D*)fbkgsummed_sub->Get(temp_name)->Clone(temp_name);
    
    if(i==3){raw_corr_sub->Scale(1/4.);  }
  
     
    dx_eta = raw_corr_sub->GetXaxis()->GetBinWidth(1);
    dx_phi = raw_corr_sub->GetYaxis()->GetBinWidth(1);
     
    raw_corr_sub->Rebin2D(2,4);
    raw_corr_sub->Scale(1/dx_eta/dx_phi/8.);
    raw_corr_sub->SetMinimum(raw_min);
    raw_corr_sub->SetMaximum(raw_max);
    raw_corr_sub->SetLineColor(kBlack);
    raw_corr_sub->GetXaxis()->SetRangeUser(-3.,3.);
    raw_corr_sub->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    raw_corr_sub->GetXaxis()->SetTitle("#Delta#eta");
    raw_corr_sub->GetYaxis()->SetTitle("#Delta#phi");
    raw_corr_sub->GetXaxis()->CenterTitle();
    raw_corr_sub->GetYaxis()->CenterTitle();
    raw_corr_sub->GetZaxis()->CenterTitle();
    raw_corr_sub->GetYaxis()->SetTitleSize(0.06);
    raw_corr_sub->GetXaxis()->SetTitleSize(0.06);
    raw_corr_sub->GetXaxis()->SetTitleOffset(1.);
    raw_corr_sub->GetYaxis()->SetTitleOffset(1.);
    raw_corr_sub->GetZaxis()->SetTitle("S(#Delta#eta,#Delta#phi)");
    raw_corr_sub->GetZaxis()->SetTitleSize(0.06);
    raw_corr_sub->GetZaxis()->SetTitleOffset(1.);
    raw_corr_sub->GetZaxis()->SetNdivisions(205);
     
      
    raw_corr_sub->Draw("surf1");

    drawlabels_PAS_bkg(i);

    TLatex *subleading_tex = new TLatex(0.66,0.72,"Subleading Jet");
    subleading_tex->SetName("tex04");
    subleading_tex->SetNDC();
    subleading_tex->SetTextFont(63);
    subleading_tex->SetTextSizePixels(tspixels);
    subleading_tex->Draw();

    cFigure1[i]->cd(5);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);

    temp_name.ReplaceAll("Raw_Yield","Mixed_Event");

    TH2D *mixed_bkg_sub = (TH2D*) fbkgsummed_sub->Get(temp_name)->Clone(temp_name);
      
    //   if(i==3){mixed_bkg_sub->Scale(1/4.);  }

    mixed_bkg_sub->Rebin2D(2,4);
    mixed_bkg_sub->Scale(1/8.);


    mixed_bkg_sub->GetXaxis()->SetRangeUser(-3.,3.);
    mixed_bkg_sub->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.0001);
    mixed_bkg_sub->SetLineColor(kBlack);
    mixed_bkg_sub->GetXaxis()->SetTitle("#Delta#eta");
    mixed_bkg_sub->GetYaxis()->SetTitle("#Delta#phi");
    mixed_bkg_sub->GetXaxis()->CenterTitle();
    mixed_bkg_sub->GetYaxis()->CenterTitle();
    mixed_bkg_sub->GetZaxis()->CenterTitle();
    mixed_bkg_sub->GetYaxis()->SetTitleSize(0.06);
    mixed_bkg_sub->GetXaxis()->SetTitleSize(0.06);
    mixed_bkg_sub->GetZaxis()->SetTitle("ME(#Delta#eta,#Delta#phi)");
    mixed_bkg_sub->GetXaxis()->SetTitleOffset(1.);
    mixed_bkg_sub->GetYaxis()->SetTitleOffset(1.);
    mixed_bkg_sub->GetZaxis()->SetTitleOffset(1.);
    mixed_bkg_sub->SetMinimum(mixed_min);
    mixed_bkg_sub->SetMaximum(mixed_max);
    
    mixed_bkg_sub->GetXaxis()->SetTitle("#Delta#eta");
    mixed_bkg_sub->GetYaxis()->SetTitle("#Delta#phi");
    mixed_bkg_sub->GetZaxis()->SetTitle("ME(#Delta#eta,#Delta#phi)");
    mixed_bkg_sub->GetZaxis()->SetTitleSize(0.06);
    mixed_bkg_sub->GetZaxis()->SetNdivisions(205);
     
    
    mixed_bkg_sub->Draw("surf1");

    drawlabels_PAS_bkg(i);

    subleading_tex->Draw();

    cFigure1[i]->cd(6);
 
    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);


    temp_name.ReplaceAll("Mixed_Event","Yield_PbPb");
    TH2D *raw_yield_sub = (TH2D*) fbkgsummed_sub->Get(temp_name)->Clone(temp_name);
 
    if(i==3){raw_yield_sub->Scale(1/4.);}
    
    raw_yield_sub->Rebin2D(2,4);
    raw_yield_sub->Scale(1/dx_eta/dx_phi/8.);
      
    raw_yield_sub->GetXaxis()->SetRangeUser(-3.,3.);
    raw_yield_sub->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    raw_yield_sub->GetXaxis()->SetTitle("#Delta#eta");
    raw_yield_sub->GetYaxis()->SetTitle("#Delta#phi");
    raw_yield_sub->GetXaxis()->CenterTitle();
    raw_yield_sub->GetYaxis()->CenterTitle();
    raw_yield_sub->GetZaxis()->CenterTitle();
    raw_yield_sub->GetYaxis()->SetTitleSize(0.06);
    raw_yield_sub->GetXaxis()->SetTitleSize(0.06);
    raw_yield_sub->GetZaxis()->SetTitle("S(#Delta#eta,#Delta#phi) / ME(#Delta#eta,#Delta#phi)");
    raw_yield_sub->GetZaxis()->SetTitleSize(0.06);
    raw_yield_sub->GetXaxis()->SetTitleOffset(1.);
    raw_yield_sub->GetYaxis()->SetTitleOffset(1.);
    raw_yield_sub->GetZaxis()->SetTitleOffset(1.);
     
    raw_yield_sub->SetMinimum(yield_min);
    raw_yield_sub->SetMaximum(yield_max);
    raw_yield_sub->GetZaxis()->SetNdivisions(205);
    raw_yield_sub->SetLineColor(kBlack);
    raw_yield_sub->Draw("surf1");
      
    drawlabels_PAS_bkg(i);

    subleading_tex->Draw();


    cout<<"done fig 1"<<endl;
    
    //---------------------------
    //   PAS FIGURE 2
    //--------------------------



    cFigure2[i]->cd(1);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);

    raw_yield->Draw("surf1");
    drawlabels_PAS_bkg(i);

    leading_tex->Draw();
    
  
    cFigure2[i]->cd(2);
      
    gPad->SetTopMargin(0.15);

    llimiteta = raw_yield->GetXaxis()->FindBin(-etalim+0.001);
    rlimiteta = raw_yield->GetXaxis()->FindBin(etalim-0.001);

    temp_name.ReplaceAll("Yield","Summed_bkg");
   
    TH1D *fit_bkg = (TH1D*)fbkgfit->Get(temp_name)->Clone(temp_name);
    
    TString func_name = "fitfunc4";
    func_name+=i;
    func_name+=3;
    TF1 *fitfunc = (TF1*)fbkgfit->Get(func_name)->Clone(func_name);

    cout<<func_name<<endl;

    if(i==3){fit_bkg->Scale(1/4.);}

    fit_bkg->SetMarkerColor(kBlack);
   
    fit_bkg->SetLineColor(kBlack);
    
    fit_bkg->Rebin(4);
    fit_bkg->Scale(1/dx_phi/dx_eta/4.);
    fit_bkg->SetMinimum(yield_min);
    fit_bkg->SetMaximum(yield_max);
    fit_bkg->GetXaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    fit_bkg->GetXaxis()->SetTitle("#Delta#phi");
    fit_bkg->GetXaxis()->SetLabelSize(0.045);
    fit_bkg->GetXaxis()->SetTitleSize(0.06);
    fit_bkg->GetXaxis()->SetTitleOffset(0.8);
    fit_bkg->GetYaxis()->SetLabelSize(0.045);
    fit_bkg->GetYaxis()->SetTitleSize(0.06);
    fit_bkg->GetYaxis()->SetTitleOffset(1.);
    fit_bkg->GetXaxis()->CenterTitle();
    fit_bkg->GetYaxis()->CenterTitle();
    fit_bkg->GetYaxis()->SetNdivisions(205);
    fit_bkg->GetYaxis()->SetTitle("B(#Delta#phi)");
  
    fit_bkg->Draw("hist pe1");


    float par0 = fitfunc->GetParameter(0)/dx_eta/dx_phi; 
    float par1 = fitfunc->GetParameter(1);
    float par2 = fitfunc->GetParameter(2);
    float par3 = fitfunc->GetParameter(3);
    float par4 = fitfunc->GetParameter(4)/dx_eta/dx_phi;
    float par5 = fitfunc->GetParameter(5);
    float par6 = fitfunc->GetParameter(6);
    

    if(i==3){
      par0 = fitfunc->GetParameter(0)/dx_eta/dx_phi/4.; 
      par4 = fitfunc->GetParameter(4)/dx_eta/dx_phi/4.; 
    }


    gen_gaus->FixParameter(0,par0);
    gen_gaus->FixParameter(1,par1);
    gen_gaus->FixParameter(2,par2);
    gen_gaus->FixParameter(3,par3);
    gen_gaus->FixParameter(4,par4);
    gen_gaus->FixParameter(5,par5);
    gen_gaus->FixParameter(6,par6);
  
    gen_gaus_bg->FixParameter(0,par0);
    gen_gaus_bg->FixParameter(1,par1);
    gen_gaus_bg->FixParameter(2,par2);
    gen_gaus_bg->FixParameter(3,par3);
    gen_gaus_bg->FixParameter(4,0.);
    gen_gaus_bg->FixParameter(5,0.);
    gen_gaus_bg->FixParameter(6,0.);
 
    gen_gaus_bg->SetLineColor(kRed);
    gen_gaus_bg->SetLineStyle(2);
    
    gen_gaus_bg->Draw("same");
    
    gen_gaus->Draw("same");
    
    switch(i){
    case 0:
      gen_gaus_up->SetParameter(0,par0+.18);
      gen_gaus_down->SetParameter(0,par0-.18);
      break;
    case 1:
      gen_gaus_up->SetParameter(0,par0+.12);
      gen_gaus_down->SetParameter(0,par0-.12);
      break;
    case 2:
      gen_gaus_up->SetParameter(0,par0+.05);
      gen_gaus_down->SetParameter(0,par0-.05);
      break;
    case 3:
      gen_gaus_up->SetParameter(0,par0+.02);
      gen_gaus_down->SetParameter(0,par0-.02);
      break;


    }

    gen_gaus_up->SetParameter(1,par1);
    gen_gaus_down->SetParameter(1,par1);
    gen_gaus_up->SetParameter(2,par2);
    gen_gaus_down->SetParameter(2,par2);
    gen_gaus_up->SetParameter(3,par3);
    gen_gaus_down->SetParameter(3,par3);
    gen_gaus_up->SetParameter(4,par4);
    gen_gaus_down->SetParameter(4,par4);
    gen_gaus_up->SetParameter(5,par5);
    gen_gaus_down->SetParameter(5,par5);
    gen_gaus_up->SetParameter(6,par6);
    gen_gaus_down->SetParameter(6,par6);

    gen_gaus_up->SetLineColor(kYellow);
    gen_gaus_down->SetLineColor(kYellow);
    gen_gaus_up->Draw("same");
    gen_gaus_down->Draw("same");

   
    fit_bkg->Draw("same hist pe1");
   
    TLegend *legend = new TLegend(0.18,0.5,0.6,0.65);
    legend->AddEntry(fit_bkg,"Sideband Background");
    legend->AddEntry(gen_gaus,"Background Fit");
    legend->AddEntry(gen_gaus_up,"Systematic Uncertainty");
    legend->SetLineColor(kWhite);
    legend->SetTextSize(0.045);
    
    legend->Draw();

    drawlabels_PAS_bkg2(i);
    
   

    cFigure2[i]->cd(3);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);


    temp_name.ReplaceAll("Summed_bkg","Result");
    TH2D *result = (TH2D*)fbkgfit->Get(temp_name)->Clone(temp_name);
    if(i==3){result->Scale(1/4.);}
    
    result->Rebin2D(2,4);
    result->Scale(1/8./dx_eta/dx_phi);
      
    result->GetXaxis()->SetRangeUser(-3.,3.);
    result->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    result->GetXaxis()->SetTitle("#Delta#eta");
    result->GetYaxis()->SetTitle("#Delta#phi");
    result->GetZaxis()->SetTitle("S(#Delta#eta,#Delta#phi) - B(#Delta#eta,#Delta#phi)");
    result->GetXaxis()->CenterTitle();
    result->GetYaxis()->CenterTitle();
    result->GetZaxis()->CenterTitle();
    result->GetXaxis()->SetLabelSize(0.045);
    result->GetXaxis()->SetTitleSize(0.06);
    result->GetXaxis()->SetTitleOffset(1.);
    result->GetYaxis()->SetLabelSize(0.045);
    result->GetYaxis()->SetTitleSize(0.06);
    result->GetYaxis()->SetTitleOffset(1.);
    result->GetZaxis()->SetLabelSize(0.04);
    result->GetZaxis()->SetTitleSize(0.06);
    result->GetZaxis()->SetTitleOffset(1.);
    result->GetZaxis()->SetNdivisions(205);
    result->SetMinimum(result_min);
    result->SetMaximum(result_max);
    result->SetLineColor(kBlack);
    result->Draw("surf1");
      
    drawlabels_PAS_bkg(i);

    leading_tex->Draw();
 

    cFigure2[i]->cd(4);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);

    raw_yield_sub->Draw("surf1");
    drawlabels_PAS_bkg(i);

    subleading_tex->Draw();

    cFigure2[i]->cd(5);
      
    gPad->SetTopMargin(0.15);

    llimiteta = raw_yield_sub->GetXaxis()->FindBin(-etalim+0.001);
    rlimiteta = raw_yield_sub->GetXaxis()->FindBin(etalim-0.001);

    temp_name.ReplaceAll("Result","Summed_bkg");

   
    TH1D *fit_bkg_sub = (TH1D*)fbkgfit_sub->Get(temp_name)->Clone(temp_name);
    
    func_name = "fitfunc2";
    func_name+=i;
    func_name+=3;
    fitfunc = (TF1*)fbkgfit_sub->Get(func_name)->Clone(func_name);

    fit_bkg_sub->GetFunction("gen_gaus")->SetLineColor(kWhite);


    if(i==3){fit_bkg_sub->Scale(1/4.);}

    fit_bkg_sub->SetMarkerColor(kBlack);
   
    fit_bkg_sub->SetLineColor(kBlack);
    fit_bkg_sub->Rebin(4);
    fit_bkg_sub->Scale(1/dx_phi/dx_eta/4.);
    fit_bkg_sub->SetMinimum(yield_min);
    fit_bkg_sub->SetMaximum(yield_max);
    fit_bkg_sub->GetXaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    fit_bkg_sub->GetXaxis()->SetTitle("#Delta#phi");
    fit_bkg_sub->GetXaxis()->SetLabelSize(0.045);
    fit_bkg_sub->GetXaxis()->SetTitleSize(0.06);
    fit_bkg_sub->GetXaxis()->SetTitleOffset(0.8);
    fit_bkg_sub->GetYaxis()->SetLabelSize(0.045);
    fit_bkg_sub->GetYaxis()->SetTitleSize(0.06);
    fit_bkg_sub->GetYaxis()->SetTitleOffset(1.);
    fit_bkg_sub->GetXaxis()->CenterTitle();
    fit_bkg_sub->GetYaxis()->CenterTitle();
    fit_bkg_sub->GetYaxis()->SetNdivisions(205);
    fit_bkg_sub->GetYaxis()->SetTitle("B(#Delta#phi)");
     
    fit_bkg_sub->Draw("hist pe1");
   

    par0 = fitfunc->GetParameter(0)/dx_eta/dx_phi; 
    par1 = fitfunc->GetParameter(1);
    par2 = fitfunc->GetParameter(2);
    par3 = fitfunc->GetParameter(3);
    par4 = fitfunc->GetParameter(4)/dx_eta/dx_phi;
    par5 = fitfunc->GetParameter(5);
    par6 = fitfunc->GetParameter(6);


    if(i==3){
      par0 = fitfunc->GetParameter(0)/dx_eta/dx_phi/4.; 
      par4 = fitfunc->GetParameter(4)/dx_eta/dx_phi/4.; 
    }



    sub_gen_gaus->FixParameter(0,par0);
    sub_gen_gaus->FixParameter(1,par1);
    sub_gen_gaus->FixParameter(2,par2);
    sub_gen_gaus->FixParameter(3,par3);
    sub_gen_gaus->FixParameter(4,par4);
    sub_gen_gaus->FixParameter(5,par5);
    sub_gen_gaus->FixParameter(6,par6);

    sub_gen_gaus_bg->FixParameter(0,par0);
    sub_gen_gaus_bg->FixParameter(1,par1);
    sub_gen_gaus_bg->FixParameter(2,par2);
    sub_gen_gaus_bg->FixParameter(3,par3);
    sub_gen_gaus_bg->FixParameter(4,0.);
    sub_gen_gaus_bg->FixParameter(5,0.);
    sub_gen_gaus_bg->FixParameter(6,0.);
 
    sub_gen_gaus_bg->SetLineColor(kRed);
    sub_gen_gaus_bg->SetLineStyle(2);
    

    sub_gen_gaus_bg->Draw("same");
    

    sub_gen_gaus->Draw("same");

    switch(i){
    case 0:
      sub_gen_gaus_up->SetParameter(0,par0+.18);
      sub_gen_gaus_down->SetParameter(0,par0-.18);
      break;
    case 1:
      sub_gen_gaus_up->SetParameter(0,par0+.12);
      sub_gen_gaus_down->SetParameter(0,par0-.12);
      break;
    case 2:
      sub_gen_gaus_up->SetParameter(0,par0+.05);
      sub_gen_gaus_down->SetParameter(0,par0-.05);
      break;
    case 3:
      sub_gen_gaus_up->SetParameter(0,par0+.02);
      sub_gen_gaus_down->SetParameter(0,par0-.02);
      break;


    }

    sub_gen_gaus_up->SetParameter(1,par1);
    sub_gen_gaus_down->SetParameter(1,par1);
    sub_gen_gaus_up->SetParameter(2,par2);
    sub_gen_gaus_down->SetParameter(2,par2);
    sub_gen_gaus_up->SetParameter(3,par3);
    sub_gen_gaus_down->SetParameter(3,par3);
    sub_gen_gaus_up->SetParameter(4,par4);
    sub_gen_gaus_down->SetParameter(4,par4);
    sub_gen_gaus_up->SetParameter(5,par5);
    sub_gen_gaus_down->SetParameter(5,par5);
    sub_gen_gaus_up->SetParameter(6,par6);
    sub_gen_gaus_down->SetParameter(6,par6);

    sub_gen_gaus_up->SetLineColor(kYellow);
    sub_gen_gaus_down->SetLineColor(kYellow);
    sub_gen_gaus_up->Draw("same");
    sub_gen_gaus_down->Draw("same");


   


    fit_bkg_sub->Draw("same hist pe1");
   
     
    legend->Draw();

    drawlabels_PAS_bkg2(i);
    
   

    cFigure2[i]->cd(6);

    gPad->SetTheta(60.839);
    gPad->SetPhi(38.0172);


    temp_name.ReplaceAll("Summed_bkg","Result");
    TH2D *result_sub = (TH2D*)fbkgfit_sub->Get(temp_name)->Clone(temp_name);

    if(i==3){result_sub->Scale(1/4.);}
    
    result_sub->Rebin2D(2,4);
    result_sub->Scale(1/8./dx_eta/dx_phi);
      
    result_sub->GetXaxis()->SetRangeUser(-3.,3.);
    result_sub->GetYaxis()->SetRangeUser(-1.5,3*TMath::Pi()/2-0.001);
    result_sub->GetXaxis()->SetTitle("#Delta#eta");
    result_sub->GetYaxis()->SetTitle("#Delta#phi");
    result_sub->GetZaxis()->SetTitle("S(#Delta#eta,#Delta#phi) - B(#Delta#eta,#Delta#phi)");
    result_sub->GetXaxis()->CenterTitle();
    result_sub->GetYaxis()->CenterTitle();
    result_sub->GetZaxis()->CenterTitle();
    result_sub->GetXaxis()->SetLabelSize(0.045);
    result_sub->GetXaxis()->SetTitleSize(0.06);
    result_sub->GetXaxis()->SetTitleOffset(1.);
    result_sub->GetYaxis()->SetLabelSize(0.045);
    result_sub->GetYaxis()->SetTitleSize(0.06);
    result_sub->GetYaxis()->SetTitleOffset(1.);
    result_sub->GetZaxis()->SetLabelSize(0.04);
    result_sub->GetZaxis()->SetTitleSize(0.06);
    result_sub->GetZaxis()->SetTitleOffset(1.);
    result_sub->GetZaxis()->SetNdivisions(205);
    result_sub->SetMinimum(result_min);
    result_sub->SetMaximum(result_max);
    result_sub->SetLineColor(kBlack);
    result_sub->Draw("surf1");
      
    drawlabels_PAS_bkg(i);
 
    subleading_tex->Draw();


    //Other plots go by centrality class


    int ndivisions;
    float diff_min, diff_max;

    switch(i){
    case 0: 
      result_min = -.8;
      result_max = 10.7;
      diff_min = -.8;
      diff_max = 10.7;
      ndivisions = 510;
      break;
    case 1: 
      result_min = -.8;
      result_max = 6.9;
      diff_min = -.8;
      diff_max = 6.9;
      ndivisions = 612;
      break;
    case 2: 
      result_min = -.45;
      result_max = 5.6;
      diff_min = -.45-.4;
      diff_max = 5.6-.4;
      ndivisions = 408;
      break;
    case 3: 
      result_min = -.25;
      result_max = 3.9;
      diff_min = -.25-.6;
      diff_max = 3.9-.6;
      ndivisions = 306;
      break;
    }

    for (int j=0; j<4; j++){
		   
	
      //---------------------
      //  PAS FIGURE 3
      //----------------------
	



    
	
      cFigure3[i]->cd(j+1);
      check_new_eta_rebin[0][i][j]->SetMaximum(result_max);
      check_new_eta_rebin[0][i][j]->SetMinimum(result_min);
      check_new_eta_rebin[0][i][j]->GetXaxis()->SetLabelSize(0.);
      check_new_eta_rebin[0][i][j]->GetYaxis()->SetNdivisions(ndivisions);
          
      check_new_eta_rebin[0][i][j]->Draw();
      check_new_eta_syst[0][i][j]->SetMarkerSize(1.);
      check_new_eta_syst[0][i][j]->SetMarkerStyle(20);
      check_new_eta_syst[0][i][j]->Draw("same e2");
      check_new_eta_rebin[0][i][j]->Draw("same");
      check_new_eta_rebin[1][i][j]->SetMarkerStyle(24);
      check_new_eta_rebin[1][i][j]->Draw("same");
     
      if(draw_ref) check_new_eta_ref[0][i][j]->Draw("same");
      if(draw_ref) check_new_eta_ref[1][i][j]->Draw("same");
      drawlabels(1,i,j);

	 

      if(j==3){
	tex24eta = new TLatex(textalign,texty2,phirangelabel);
	tex24eta->SetName("tex24eta");
	tex24eta->SetNDC();
	tex24eta->SetTextSizePixels(tspixels);
	tex24eta->Draw();
      }

      lineEta = new TLine(-1.5,0,1.5,0); 
      lineEta->Draw("same");
	
      if(j==0){ 
	l40 = new TLegend(textalign2-0.03,texty3-0.1,.95,texty2-0.05);
	l40->SetName("l40");
	l40->SetTextSizePixels(tspixels);
	l40->SetFillColor(kWhite);
	l40->SetLineColor(kWhite);
	l40->AddEntry(check_new_eta_syst[0][i][j],"PbPb Inclusive Jets","lpfe");
	l40->AddEntry(check_new_eta_rebin[1][i][j],"pp Inclusive Jets","lpfe");
	l40->Draw("same");
      }


      cFigure3[i]->cd(j+5);
      PbPb_pp_eta_syst[1][i][j]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV/c)^{-1}");


      PbPb_pp_eta_syst[1][i][j]->SetMinimum(diff_min);
      PbPb_pp_eta_syst[1][i][j]->SetMaximum(diff_max);
      PbPb_pp_eta_syst[1][i][j]->GetYaxis()->SetNdivisions(ndivisions);

      PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetLabelSize(ts2);
      PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetTitleSize(tstitle);
      if(j==0){ PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetLabelSize(ts2-0.01);
	PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetTitleSize(ts);
	PbPb_pp_eta_syst[1][i][j]->GetYaxis()->SetTitleSize(ts);
	PbPb_pp_eta_syst[1][i][j]->GetYaxis()->SetLabelSize(ts2);
	PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetLabelOffset(0.013);
	PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetTitleOffset(xoffset2);
      }
      //PbPb_pp_eta_syst[1][i][j]->GetXaxis()->SetLabelOffset(xoffset);
      PbPb_pp_eta_syst[1][i][j]->SetFillColor(90);
      //  PbPb_pp_eta_syst[1][i][j]->SetMarkerSize(0);
      PbPb_pp_eta_syst[1][i][j]->Draw("e2");
      PbPb_pp_eta[1][i][j]->Draw("same");

      if(draw_ref) PbPb_pp_eta_ref[1][i][j]->Draw("same");
	 
      lineEta->Draw("same");

      TPave *labelcover = new TPave(1.3,diff_min-0.8,1.6,diff_min-0.05);
      labelcover->SetLineColor(kWhite);
      labelcover->SetOption("nb");
      labelcover->SetFillColor(kWhite);
      if(j<3){ labelcover->Draw(); }
	    
	 	 
      //-----------
      // Figure 4
      //-----------


      cFigure4[i]->cd(j+1);
      check_new_phi_rebin[0][i][j]->SetMaximum(result_max);
      check_new_phi_rebin[0][i][j]->SetMinimum(result_min);
      check_new_phi_rebin[0][i][j]->GetYaxis()->SetNdivisions(ndivisions);
      check_new_phi_rebin[0][i][j]->GetXaxis()->SetLabelSize(0.);
      check_new_phi_rebin[0][i][j]->Draw();
      check_new_phi_syst[0][i][j]->Draw("same e2");
      check_new_phi_rebin[0][i][j]->Draw("same");

      if(draw_ref) check_new_phi_ref[0][i][j]->Draw("same");
      if(draw_ref) check_new_phi_ref[1][i][j]->Draw("same");

      check_new_phi_rebin[1][i][j]->Draw("same");
	 
      drawlabels(1,i,j);
	  
	  
      if(j==3){
	tex24phi = new TLatex(textalign,texty2,etarangelabel);
	tex24phi->SetName("tex24phi");
	tex24phi->SetNDC();
	tex24phi->SetTextSizePixels(tspixels);
	tex24phi->Draw();
      }


      linePhi = new TLine(-TMath::Pi()/2,0.,TMath::Pi()/2,0.);
      linePhi->Draw("same");

      if(j==0){ 
	l40->Draw("same");
      }

      cFigure4[i]->cd(j+5);
      PbPb_pp_phi_syst[1][i][j]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV/c)^{-1}");
      PbPb_pp_phi_syst[1][i][j]->SetMinimum(diff_min);
      PbPb_pp_phi_syst[1][i][j]->SetMaximum(diff_max);
      PbPb_pp_phi_syst[1][i][j]->GetYaxis()->SetNdivisions(ndivisions);
      PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetRangeUser(-1.4999,1.5);
      PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetLabelSize(ts2);
      PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetTitleSize(tstitle);
     
      if(j==0){ PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetLabelSize(ts2-0.01);
	PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetTitleSize(ts);
	PbPb_pp_phi_syst[1][i][j]->GetYaxis()->SetTitleSize(ts);
	PbPb_pp_phi_syst[1][i][j]->GetYaxis()->SetLabelSize(ts2);
	PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetTitleOffset(xoffset2);
	PbPb_pp_phi_syst[1][i][j]->GetXaxis()->SetLabelOffset(0.013);
      }
      PbPb_pp_phi_syst[1][i][j]->SetFillColor(90);
      //   PbPb_pp_phi_syst[1][i][j]->SetMarkerSize(0);
      PbPb_pp_phi_syst[1][i][j]->Draw("e2");
      PbPb_pp_phi[1][i][j]->Draw("same p");

      if(draw_ref) PbPb_pp_phi_ref[1][i][j]->Draw("same p");
      linePhi->Draw("same");


      labelcover = new TPave(1.3,diff_min-0.8,1.6,diff_min-0.05);
      labelcover->SetLineColor(kWhite);
      labelcover->SetOption("nb");
      labelcover->SetFillColor(kWhite);
      if(j<3){ labelcover->Draw(); }
      
      //---------------------
      //  PAS FIGURES 5 & 6
      //----------------------
	 
     
    

      cFigure5[i]->cd(j+1);
      check_new_eta_rebin[4][i][j]->SetMaximum(result_max);
      check_new_eta_rebin[4][i][j]->SetMinimum(result_min);
      check_new_eta_rebin[4][i][j]->GetXaxis()->SetLabelSize(0.);
      check_new_eta_rebin[4][i][j]->GetYaxis()->SetNdivisions(ndivisions);
      check_new_eta_rebin[4][i][j]->Draw();
      check_new_eta_syst[4][i][j]->SetMarkerSize(1.);
      check_new_eta_syst[4][i][j]->SetMarkerStyle(21);
      check_new_eta_syst[4][i][j]->Draw("same e2");
      check_new_eta_rebin[4][i][j]->Draw("same");
      check_new_eta_rebin[5][i][j]->SetMarkerStyle(25);
      check_new_eta_rebin[5][i][j]->Draw("same");
      if(draw_ref) check_new_eta_ref[4][i][j]->Draw("same");
      if(draw_ref) check_new_eta_ref[5][i][j]->Draw("same");
   

      drawlabels(5,i,j);

      if(j==3){
	tex24eta = new TLatex(textalign,texty2,phirangelabel);
	tex24eta->SetName("tex24eta");
	tex24eta->SetNDC();
	tex24eta->SetTextSizePixels(tspixels);
	tex24eta->Draw();
      }
	
      lineEta = new TLine(-1.5,0,1.5,0); 
      lineEta->Draw("same");
       
      if(j==0){ 
	l40 = new TLegend(textalign2-0.03,texty3-.1,.95,texty2-0.05);
	l40->SetName("l40");
	l40->SetTextSizePixels(tspixels);
	l40->SetFillColor(kWhite);
	l40->SetLineColor(kWhite);
	l40->AddEntry(check_new_eta_syst[4][i][j],"PbPb Leading Jets","lpfe");
	l40->AddEntry(check_new_eta_rebin[5][i][j],"pp Leading Jets","lpe");
	l40->Draw("same");
      }

    

      cFigure5[i]->cd(j+5);
      check_new_eta_rebin[2][i][j]->SetMaximum(result_max);
      check_new_eta_rebin[2][i][j]->SetMinimum(result_min);
      check_new_eta_rebin[2][i][j]->GetXaxis()->SetLabelSize(0.);
      check_new_eta_rebin[2][i][j]->GetYaxis()->SetNdivisions(ndivisions);
      check_new_eta_syst[2][i][j]->SetLineColor(kBlack);
      check_new_eta_rebin[2][i][j]->Draw();
    
      check_new_eta_syst[2][i][j]->SetMarkerSize(2.);
      check_new_eta_syst[2][i][j]->SetMarkerStyle(34);
      check_new_eta_syst[2][i][j]->Draw("same e2");
      check_new_eta_rebin[2][i][j]->Draw("same p");
      check_new_eta_rebin[3][i][j]->SetMarkerStyle(28);
      check_new_eta_rebin[3][i][j]->Draw("same p");
	  
      if(draw_ref) check_new_eta_ref[2][i][j]->Draw("same");
      if(draw_ref) check_new_eta_ref[3][i][j]->Draw("same");
      lineEta->Draw("same");

      
	 
      if(j==0){ 
	l41 = new TLegend(textalign2-0.03,texty3+0.05,.95,0.95);
	l41->SetName("l41");
	l41->SetTextSizePixels(tspixels);
	l41->SetFillColor(kWhite);
	l41->SetLineColor(kWhite);
	l41->AddEntry(check_new_eta_syst[2][i][j],"PbPb Subleading Jets","lpfe");
	l41->AddEntry(check_new_eta_rebin[3][i][j],"pp Subleading Jets","lpe");
	l41->Draw("same");
      }

      cFigure5[i]->cd(j+9);

      PbPb_pp_eta_syst[3][i][j]->GetYaxis()->SetTitle("Y_{PbPb}-Y_{pp}   (GeV/c)^{-1}");
      PbPb_pp_eta_syst[3][i][j]->SetMinimum(diff_min);
      PbPb_pp_eta_syst[3][i][j]->SetMaximum(diff_max);
      PbPb_pp_eta_syst[3][i][j]->GetYaxis()->SetNdivisions(ndivisions);
      PbPb_pp_eta_syst[3][i][j]->GetXaxis()->SetLabelSize(ts2);
      PbPb_pp_eta_syst[3][i][j]->GetXaxis()->SetTitleSize(tstitle);
    
      if(j==0){ PbPb_pp_eta_syst[3][i][j]->GetXaxis()->SetLabelSize(ts2-0.01);
	PbPb_pp_eta_syst[3][i][j]->GetXaxis()->SetTitleSize(ts);
	PbPb_pp_eta_syst[3][i][j]->GetYaxis()->SetTitleSize(ts);
	PbPb_pp_eta_syst[3][i][j]->GetYaxis()->SetLabelSize(ts2);
	PbPb_pp_eta_syst[3][i][j]->GetXaxis()->SetLabelOffset(0.013);
	PbPb_pp_eta_syst[3][i][j]->GetXaxis()->SetTitleOffset(xoffset2);
      }


      PbPb_pp_eta_syst[3][i][j]->SetFillColor(30);
      //  PbPb_pp_eta_syst[3][i][j]->SetMarkerSize(0);
      PbPb_pp_eta_syst[3][i][j]->Draw("e2");
      PbPb_pp_eta[3][i][j]->Draw("same");
      PbPb_pp_eta[5][i][j]->SetMarkerStyle(21);
      PbPb_pp_eta_syst[5][i][j]->SetFillColor(kOrange-2);
      PbPb_pp_eta_syst[5][i][j]->Draw("same e2");
      PbPb_pp_eta[5][i][j]->Draw("same");

     if(draw_ref) PbPb_pp_eta_ref[3][i][j]->Draw("same");
     if(draw_ref) PbPb_pp_eta_ref[5][i][j]->Draw("same");
	 
      lineEta->Draw("same");
	    
      if(j==0){ 
	l42 = new TLegend(textalign2-0.03,texty2,.95,0.95);
	l42->SetName("l42");
	l42->SetTextSizePixels(tspixels);
	l42->SetFillColor(kWhite);
	l42->SetLineColor(kWhite);
	l42->AddEntry(PbPb_pp_eta_syst[5][i][j],"Leading Jets","lpfe");
	l42->AddEntry(PbPb_pp_eta_syst[3][i][j],"Subleading Jets","lpfe");
	l42->Draw("same");
      }
      labelcover = new TPave(1.3,diff_min-0.8,1.6,diff_min-0.05);
      labelcover->SetLineColor(kWhite);
      labelcover->SetOption("nb");
      labelcover->SetFillColor(kWhite);
      labelcover->Draw();
   

      //--------dPhi-------------
      //      (Figure 6)
	  
      cFigure6[i]->cd(j+1);
      check_new_phi_rebin[4][i][j]->SetMaximum(result_max);
      check_new_phi_rebin[4][i][j]->SetMinimum(result_min);
      check_new_phi_rebin[4][i][j]->GetXaxis()->SetLabelSize(0.);
      check_new_phi_rebin[4][i][j]->GetYaxis()->SetNdivisions(ndivisions);
      check_new_phi_rebin[4][i][j]->Draw();
      check_new_phi_syst[4][i][j]->Draw("same e2");
      check_new_phi_rebin[4][i][j]->Draw("same");
      check_new_phi_rebin[5][i][j]->SetMarkerStyle(25);
      check_new_phi_rebin[5][i][j]->Draw("same");
      if(draw_ref) check_new_phi_ref[4][i][j]->Draw("same");
      if(draw_ref) check_new_phi_ref[5][i][j]->Draw("same");

      drawlabels(5,i,j);
	  
      if(j==3){
	tex24phi = new TLatex(textalign,texty2,etarangelabel);
	tex24phi->SetName("tex24phi");
	tex24phi->SetNDC();
	tex24phi->SetTextSizePixels(tspixels);
	tex24phi->Draw();
      }

      linePhi = new TLine(-1.5,0,1.5,0);
      linePhi->Draw("same");

      if(j==0){ 
	l40->Draw("same");
      }

      cFigure6[i]->cd(j+5);
      check_new_phi_rebin[2][i][j]->SetMaximum(result_max);
      check_new_phi_rebin[2][i][j]->SetMinimum(result_min);
      check_new_phi_rebin[2][i][j]->GetXaxis()->SetLabelSize(0.);
      check_new_phi_rebin[2][i][j]->GetYaxis()->SetNdivisions(ndivisions);
      check_new_phi_syst[2][i][j]->SetLineColor(kBlack);
      check_new_phi_rebin[2][i][j]->Draw();
     
      check_new_phi_syst[2][i][j]->Draw("same e2");
      check_new_phi_rebin[2][i][j]->Draw("same");
      check_new_phi_rebin[3][i][j]->SetMarkerStyle(28);
      check_new_phi_rebin[3][i][j]->Draw("same");

      if(draw_ref) check_new_phi_ref[2][i][j]->Draw("same");
      if(draw_ref) check_new_phi_ref[3][i][j]->Draw("same");
      linePhi->Draw("same");
	    
	    
      if(j==0){ 
	l41->Draw("same");
      }


      cFigure6[i]->cd(j+9);

      PbPb_pp_phi_syst[3][i][j]->GetYaxis()->SetTitle("Y_{PbPb}- Y_{pp}   (GeV/c)^{-1}");
      PbPb_pp_phi_syst[3][i][j]->SetMinimum(diff_min);
      PbPb_pp_phi_syst[3][i][j]->SetMaximum(diff_max);
      PbPb_pp_phi_syst[3][i][j]->GetYaxis()->SetNdivisions(ndivisions);
      PbPb_pp_phi_syst[3][i][j]->SetFillColor(30);
      PbPb_pp_phi_syst[3][i][j]->GetXaxis()->SetLabelSize(ts2);
      PbPb_pp_phi_syst[3][i][j]->GetXaxis()->SetTitleSize(tstitle);
    
   
      if(j==0){ PbPb_pp_phi_syst[3][i][j]->GetXaxis()->SetLabelSize(ts2-0.01);
	PbPb_pp_phi_syst[3][i][j]->GetXaxis()->SetTitleSize(ts);
	PbPb_pp_phi_syst[3][i][j]->GetYaxis()->SetTitleSize(ts);
	PbPb_pp_phi_syst[3][i][j]->GetYaxis()->SetLabelSize(ts2);
	PbPb_pp_phi_syst[3][i][j]->GetXaxis()->SetLabelOffset(0.013);
	PbPb_pp_phi_syst[3][i][j]->GetXaxis()->SetTitleOffset(xoffset2);
      }
      PbPb_pp_phi_syst[3][i][j]->SetFillColor(30);
      //  PbPb_pp_phi_syst[3][i][j]->SetMarkerSize(0);
      PbPb_pp_phi_syst[3][i][j]->Draw("e2");
      PbPb_pp_phi[3][i][j]->Draw("same");
      PbPb_pp_phi[5][i][j]->SetMarkerStyle(21);
      PbPb_pp_phi_syst[5][i][j]->SetFillColor(kOrange-2);
      PbPb_pp_phi_syst[5][i][j]->Draw("same e2");
      PbPb_pp_phi[5][i][j]->Draw("same");

      if(draw_ref) PbPb_pp_phi_ref[3][i][j]->Draw("same");
      if(draw_ref) PbPb_pp_phi_ref[5][i][j]->Draw("same");
	  
      linePhi->Draw("same");
	    
      if(j==0){ 
	l42->Draw("same");
      }
     
      labelcover = new TPave(1.3,diff_min-0.8,1.6,diff_min-0.05);
      labelcover->SetLineColor(kWhite);
      labelcover->SetOption("nb");
      labelcover->SetFillColor(kWhite);
      if(j<3){ labelcover->Draw(); }

    } //close j
  
      //--------------------------------------------
      //   Save all PAS plots except for Figure 7
      //-----------------------------------------
  

    cFigure3[i]->cd(0);
    TLatex *canvas_title = new TLatex(0.06,0.95,"CMS Preliminary");
    canvas_title->SetTextSizePixels(tspixels);
    canvas_title->SetTextFont(63);
    canvas_title->Draw();
  
    TLatex *canvas_title2 = new TLatex(0.292,0.95,"PbPb 166 #mub^{-1} (2.76 TeV)                pp 5.3 pb^{-1} (2.76 TeV)");
    canvas_title2->SetTextSizePixels(tspixels);
    canvas_title2->Draw();


    cFigure4[i]->cd(0);
    canvas_title->Draw();
    canvas_title2->Draw();

    cFigure5[i]->cd(0);
    canvas_title = new TLatex(0.06,0.965,"CMS Preliminary");
    canvas_title->SetTextSizePixels(tspixels);
    canvas_title->SetTextFont(63);
    canvas_title->Draw();

    canvas_title2 = new TLatex(0.295,0.965,"PbPb 166 #mub^{-1} (2.76 TeV)                pp 5.3 pb^{-1} (2.76 TeV)");
    canvas_title2->SetTextSizePixels(tspixels);
    canvas_title2->Draw();

    cFigure6[i]->cd(0);
    canvas_title->Draw();
    canvas_title2->Draw();

    
    cFigure1[i]->Update();
    cFigure2[i]->Update();
    cFigure3[i]->Update();
    cFigure4[i]->Update();
    cFigure5[i]->Update();
    cFigure6[i]->Update();

 
    if(draw_ref){
      figure3_name+="_WithRef";
      figure4_name+="_WithRef";
      figure5_name+="_WithRef"; 
      figure6_name+="_WithRef";
     }

    figure1_name+=".pdf";
    cFigure1[i]->SaveAs(figure1_name);
    figure1_name.ReplaceAll("pdf","png");
    cFigure1[i]->SaveAs(figure1_name);     
  
    figure2_name+=".pdf";
    cFigure2[i]->SaveAs(figure2_name);
    figure2_name.ReplaceAll("pdf","png");
    cFigure2[i]->SaveAs(figure2_name);   

    figure3_name+=".pdf";
    cFigure3[i]->SaveAs(figure3_name);
    figure3_name.ReplaceAll("pdf","png");
    cFigure3[i]->SaveAs(figure3_name);   


    figure4_name+=".pdf";
    cFigure4[i]->SaveAs(figure4_name);
    figure4_name.ReplaceAll("pdf","png");
    cFigure4[i]->SaveAs(figure4_name);   

    figure5_name+=".pdf";
    cFigure5[i]->SaveAs(figure5_name);
    figure5_name.ReplaceAll("pdf","png");
    cFigure5[i]->SaveAs(figure5_name);   
      

    figure6_name+=".pdf";
    cFigure6[i]->SaveAs(figure6_name);
    figure6_name.ReplaceAll("pdf","png");
    cFigure6[i]->SaveAs(figure6_name);   
   figure6_name.ReplaceAll("png","C");   
    cFigure6[i]->SaveAs(figure6_name);   

     
  } //closes i loop for PAS plots
    
 

  //--------------------------
  //    Figure 7
  //--------------------------
	
	
  TString figure7_name = "PAS_Figure_7";
  cFigure7 = new TCanvas(figure7_name," ",10,10,1550,500);
  cFigure7->Divide(4,1,0.,0.);

  for(int j = 0; j<4; j++){
    cFigure7->cd(j+1);

    TString leadingname = "Leading_cent";
    leadingname +=j;

    
    Integral_eta_Pt[5][j]->SetMaximum(7.6);

    
    Integral_eta_Pt[5][j]->SetMinimum(-.8);
						
    Integral_eta_Pt[5][j]->GetXaxis()->SetLabelSize(ts2);
    if(j==0){
      Integral_eta_Pt[5][j]->GetXaxis()->SetLabelSize(ts3);
      Integral_eta_Pt[5][j]->GetXaxis()->SetTitleSize(ts2);
      Integral_eta_Pt[5][j]->GetXaxis()->SetTitleOffset(xoffset+0.2);
      Integral_eta_Pt[5][j]->GetXaxis()->SetLabelOffset(0.015);
      Integral_eta_Pt[5][j]->GetYaxis()->SetLabelSize(ts2);
      Integral_eta_Pt[5][j]->GetYaxis()->SetTitleSize(ts2-0.005);
    }
      
    Integral_eta_syst[3][j]->SetMarkerSize(2.);
    Integral_eta_syst[3][j]->SetMarkerStyle(34);
    Integral_eta_syst[5][j]->SetMarkerSize(2.);
    Integral_eta_syst[5][j]->SetMarkerStyle(21);






    Integral_eta_syst[3][j]->SetMarkerColor(kCyan+3);
    Integral_eta_syst[3][j]->SetLineColor(kCyan+3);


    Integral_eta_Pt[3][j]->SetMarkerColor(kCyan+3);
    Integral_eta_Pt[3][j]->SetLineColor(kCyan+3);

    Integral_eta_Pt[3][j]->SetMarkerStyle(34);

    Integral_eta_Pt[5][j]->Draw("p");
    Integral_eta_syst[5][j]->Draw("same e2");
    Integral_eta_Pt[5][j]->Draw("same e1p");
    Integral_eta_Pt[3][j]->Draw("same p");
    Integral_eta_syst[3][j]->Draw("same e2");
     
    Integral_eta_Pt[5][j]->Draw("same e1p");
     
      
    gPad->RedrawAxis();
     
       
    if(draw_ref) Integral_eta_ref_Pt[3][j]->Draw("same p");
    if(draw_ref) Integral_eta_ref_Pt[5][j]->Draw("same p");
    
    Integral_eta_Pt[3][j]->Draw("same e1p");

    drawlabels_int_pt2(5,j);

    linePt = new TLine(1.,0,8.,0);
    linePt->SetLineStyle(2);
    linePt->SetLineWidth(2);
    linePt->Draw("same");

    if(j==3){
      tex24phi = new TLatex(textalign,texty2,TString(etarangelabel+", |#Delta#phi|<1.0"));
      tex24phi->SetName("tex24phi");
      tex24phi->SetNDC();
      tex24phi->SetTextSizePixels(tspixels);
      tex24phi->Draw();
    }


   
    if(j==0){ 
      l42 = new TLegend(textalign2-0.03,texty3-.05,.95,texty1-.05);
      l42->SetName("l42");
      l42->SetTextFont(43);
      l42->SetTextSizePixels(tspixels);
      l42->SetFillColor(kWhite);
      l42->SetLineColor(kWhite);
      l42->AddEntry(Integral_eta_syst[5][j],"Leading Jets","lpfe");
      l42->AddEntry(Integral_eta_syst[3][j],"Subleading Jets","lpfe");
      if(draw_ref) l42->AddEntry(Integral_eta_ref_Pt[5][j],"PAS","lpfe");
      l42->Draw("same");
    }

    TPave *labelcover = new TPave(7.5,-1.5,8.1,-.85);
    labelcover->SetLineColor(kWhite);
    labelcover->SetOption("nb");
    labelcover->SetFillColor(kWhite);
    if(j<3){ labelcover->Draw(); }
	    


  } //closes j loop for Figure 7

  cFigure7->Update();

  cFigure7->cd(0);
  TLatex *canvas_title = new TLatex(0.06,0.9,"CMS Preliminary");
  canvas_title->SetTextSizePixels(tspixels);
  canvas_title->SetTextFont(63);
  canvas_title->Draw();
	     
  TLatex *canvas_title2 = new TLatex(0.292,0.9,"PbPb 166 #mub^{-1} (2.76 TeV)                pp 5.3 pb^{-1} (2.76 TeV)");
  canvas_title2->SetTextSizePixels(tspixels);
  canvas_title2->Draw();
      
      
  figure7_name = "PAS_Figure_7";
 
  if(draw_ref){
    cFigure7->SaveAs((TString)(figure7_name+"_WithRef.pdf"));
    cFigure7->SaveAs((TString)(figure7_name+"_WithRef.png"));
  }else{
  

    cFigure7->SaveAs((TString)(figure7_name+".pdf"));
    cFigure7->SaveAs((TString)(figure7_name+".png"));

    for(int j = 0; j<4; j++){
      
      cFigure7->cd(j+1);

      
      Integral_eta_Pt[1][j]->SetMarkerColor(kRed);
      Integral_eta_Pt[1][j]->SetLineColor(kRed);
      Integral_eta_Pt[1][j]->Draw("same p");

      if(j==0){
	l42->AddEntry(Integral_eta_Pt[1][j],"Inclusive Jets","lpfe");
	l42->Draw();
      }
      
    }
    cFigure7->SaveAs((TString)(figure7_name+"_WithInclusive.pdf"));
    cFigure7->SaveAs((TString)(figure7_name+"_WithInclusive.png"));
  }

  return 0;
  
} //Close main loop

