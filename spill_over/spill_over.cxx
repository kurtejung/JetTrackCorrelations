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

Int_t spill_over(int gstart= 6, int gend= 12){
  // Default values run over data only.  gstart = 6, gend = 12 are values for full MC run

#include "../HIN_14_016_universals.h"

  //Set Style:

 
  gStyle->SetOptStat(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.15);
  gStyle->SetPadLeftMargin  (0.16);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadTickX       (1);
  gStyle->SetPadTickY       (1);
  gStyle->SetTextFont(43);
  gStyle->SetCanvasBorderMode(0);
 
 

  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 5;


  enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, n_data_mc_types};

  TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack"};
  int data_mc_type_code = -999;

  float CBins[nCBins+1] = {0, 10, 30, 50, 100};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%", "Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };
  TString TrkPtBin_labels[nTrkPtBins] = {"1<p_{T}^assoc.<2", "2<p_{T}^assoc.<3","3<p_{T}^assoc.<4","p_{T}^assoc.>8"};

  
  TF1 *fit0 = new TF1("fit0","[0]",-3.,3.);
 
  float x;
  TF1 *do_offset = new TF1("do_offset","-1.*[0]+x-x",-3.,3.);
  float offset;
   
  TFile *fin[12];
  TFile *fin_ref[12];
  TFile *fout[12];
  TFile *fclosures[12];

  TH2D *result[12][5][4][2];
  TH2D *result2[12][5][4][2];


  TH2D* background[12][5][4][2];
  TH1D* background_left[12][5][4][2];
  TH1D* background_right[12][5][4][2];
  TH1D* background_proj[12][5][4][2];


  TH1D *phi_proj[12][5][4][2];
  TH1D *phi_proj_rebin[12][5][4][2];
  TH1D *phi_proj_rebin2[12][5][4][2];
  TH1D *eta_proj[12][5][4][2];
  TH1D *eta_proj_rebin[12][5][4][2];
  TH1D *eta_proj_rebin2[12][5][4][2];

  TH1D *eta_proj_ref[12][5][4][2];
  TH1D *phi_proj_ref[12][5][4][2];
 
  TH1D *closure_eta[12][5][4][2];
  TH1D *closure_phi[12][5][4][2];

  TH2D *closure_2D_direct[12][5][4][2];
  TH2D *closure_2D_direct_gen[12][5][4][2];
  TH1D *closure_eta_direct[12][5][4][2];
  TH1D *closure_phi_direct[12][5][4][2];
  TH1D *closure_eta_direct_rebin[12][5][4][2];
  TH1D *closure_phi_direct_rebin[12][5][4][2];

  TH1D *closure_eta_ref[12][5][4][2];
  TH1D *closure_phi_ref[12][5][4][2];

  TH1D *HYDJET_PYTHIA_eta[12][5][4][2];
  TH1D *HYDJET_PYTHIA_phi[12][5][4][2];
 
  TCanvas *ccheckEta_wide[12][2];
  TCanvas *ccheckEta[12][2];
  TCanvas *ccheckPhi[12][2];
  
  TCanvas *cHminusP_phi[12][2];
  TCanvas *cHminusP_eta[12][2];
 
  TString in_name, plotname, outname, funcname, centlabel, datalabel, jettype,jettype2, pTlabel,checkcanvasnameEta,checkcanvasnamePhi,rawetacanvasname,HminusPphicanvasname, HminusPetacanvasname,PbPbminuscanvas_phi,PbPbminuscanvas_eta, checkcanvasnameEta_wide, EtaClosureName;

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

  vector<double> Closure_integral_eta0;
  vector<double> Closure_integral_phi0;
  vector<double> Closure_integral_eta1;
  vector<double> Closure_integral_phi1;
  vector<double> Closure_integral_eta2;
  vector<double> Closure_integral_phi2;
  vector<double> Closure_integral_eta3;
  vector<double> Closure_integral_phi3;
 
  TGraph *Closure_integral_eta_pT[12][4][2];
  TGraph *Closure_integral_phi_pT[12][4][2];
 
  Double_t check_ymax, check_ymin, dx_eta, dx_phi, bc, err, evalpt, temp1, err1;

  TF1 *gaus1d = new TF1("gaus1d","[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))");

  TF1 *gaus_phi[12][5][4];
  TF1 *gaus_eta[12][5][4];
 
  TLegend *lcheck, *leta, *lHminusP;
  
  TLine *linePhi, *lineEta;
 
  int mc_type_code;

  int llimiteta1, rlimiteta1,llimiteta2, rlimiteta2; 

  /////////////////////////

  etalim = 1.5;
  philim = 1.5;

  //-------------------------------------------------- 
  // Open data and output files
  //-------------------------------------------------
 
  // gend = 8;

  for(int g=gstart; g<gend; g++){

    //  There will only be one big "g-loop".

    switch(g){
    case 6:
      fin[g] = new TFile("../me_correct_mc/HydJet_RecoJet_GenTrack_Inclusive_Correlations.root","READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/Inclusive_Closures.root","READ");
      mc_type_code = 2;
      datalabel = "RecoJet_Inclusive";
      jettype = "";
      jettype2 = "Inclusive_";
      break;
    case 7:
      //fin[g] = new TFile("../me_correct_mc/HydJet_GenJet_GenTrack_Inclusive_Correlations.root","READ");
      fin[g] = new TFile("../me_correct_mc/Pythia_RecoJet_GenTrack_Inclusive_Correlations.root","READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/Inclusive_Closures.root","READ");
      fout[g] = new TFile("Inclusive_SpillOvers.root", "RECREATE");
      datalabel = "RecoJet_Inclusive";
      jettype = "";
      jettype2 = "Inclusive_";
      mc_type_code = 2;
      break;
    case 8:
      fin[g] = new TFile("../me_correct_mc/HydJet_RecoJet_GenTrack_SubLeading_Correlations.root","READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/SubLeading_Closures.root","READ");
      datalabel = "GenJet_SubLeading";
      jettype = "SubLeading_";
      jettype2 = "SubLeading_";
      mc_type_code = 2;
      break;
    case 9:
      //fin[g] = new TFile("../me_correct_mc/HydJet_GenJet_GenTrack_SubLeading_Correlations.root","READ");
      fin[g] = new TFile("../me_correct_mc/Pythia_RecoJet_GenTrack_SubLeading_Correlations.root","READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/SubLeading_Closures.root","READ");
      fout[g] = new TFile("SubLeading_SpillOvers.root", "RECREATE");
      datalabel = "GenJet_SubLeading";
      jettype = "SubLeading_";
      jettype2 = "SubLeading_";
     mc_type_code = 2;
      break;
    case 10:
        fin[g] = new TFile("../me_correct_mc/HydJet_RecoJet_GenTrack_Leading_Correlations.root","READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/Leading_Closures.root","READ");
      mc_type_code = 2;
      datalabel = "GenJet_Leading";
      jettype = "Leading_";
      jettype2 = "Leading_";
      break;
    case 11:
      //fin[g] = new TFile("../me_correct_mc/HydJet_GenJet_GenTrack_Leading_Correlations.root","READ");
      fin[g] = new TFile("../me_correct_mc/Pythia_RecoJet_GenTrack_Leading_Correlations.root","READ");
      fin_ref[g] = new TFile("../HIN_14_016_FROZEN_PUBLIC_PAS/analysis/Leading_Closures.root","READ");
      fout[g] = new TFile("Leading_SpillOvers.root", "RECREATE");
      datalabel = "GenJet_Leading";
      jettype = "Leading_";
      jettype2 = "Leading_";
      mc_type_code = 2;
      break;
    }

    //----------------------------------------------
    //     Start i and l loops -- set up canvases
    //------------------------------------------------

    for(int l = 0; l<1; l++){
    
      // cout << g << endl; continue;
      if(l>0&&g<6){ continue; } //no Gen/Reco for data

      //  if(l<1&&g>5){ continue; } //not interested in Gen
     
      ccheckEta[g][l] = new TCanvas(Form("EtaProjections%d%d",g,l)," ",10,10,1500,1600);
      ccheckEta[g][l]->Divide(4,4,0.,0.);
     
      ccheckEta_wide[g][l] = new TCanvas(Form("WideEtaProjections%d%d",g,l)," ",10,10,1500,1600);
      ccheckEta_wide[g][l]->Divide(4,4,0.,0.);

      ccheckPhi[g][l] = new TCanvas(Form("PhiProjections%d%d",g,l)," ",10,10,1500,1600);
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

	  TString in_name = "Yield_"; in_name+= jettype; in_name+= data_mc_type_strs[mc_type_code]; in_name+="_"; in_name+= CBin_strs[3-j]; in_name+="_"; in_name+= CBin_strs[4-j]; in_name+= "_"; in_name+=TrkPtBin_strs[i]; in_name+="_"; in_name+=TrkPtBin_strs[i+1];

	  result[g][i][j][l] = (TH2D*)fin[g]->Get(in_name)->Clone(in_name);

	  if(g%2==0){
	    llimiteta1 = result[g][i][j][l]->GetXaxis()->FindBin(-2.+0.0001);
	    rlimiteta1 = result[g][i][j][l]->GetXaxis()->FindBin(-1.0-0.0001);

	  	    
	    llimiteta2 = result[g][i][j][l]->GetXaxis()->FindBin(1.+0.0001);
	    rlimiteta2 = result[g][i][j][l]->GetXaxis()->FindBin(2.-0.0001);

	  }else{
	    llimiteta1 = result[g][i][j][l]->GetXaxis()->FindBin(-2.+0.0001);
	    rlimiteta1 = result[g][i][j][l]->GetXaxis()->FindBin(-1.0-0.0001);
	    
	    llimiteta2 = result[g][i][j][l]->GetXaxis()->FindBin(1.+0.0001);
	    rlimiteta2 = result[g][i][j][l]->GetXaxis()->FindBin(2.-0.0001);
	  }


	  background_left[g][i][j][l] = (TH1D*)result[g][i][j][l]->ProjectionY(Form("LeftSideBackground%d%d%d%d",g,i,j,l),llimiteta1,rlimiteta1);
	    
	  background_proj[g][i][j][l] = (TH1D*)result[g][i][j][l]->ProjectionY(Form("ProjectedBackground%d%d%d%d",g,i,j,l),llimiteta2,rlimiteta2);

	  background_proj[g][i][j][l]->Add(background_left[g][i][j][l]);

	  background_proj[g][i][j][l]->Scale(1./2/(rlimiteta1-llimiteta1+1));
	    
	  background[g][i][j][l] = (TH2D*)result[g][i][j][l]->Clone(Form("Background%d%d%d%d",g,i,j,l));

	  for(int k = 1;  k<result[g][i][j][l]->GetNbinsY(); k++){
	    temp1 = background_proj[g][i][j][l]->GetBinContent(k);
	    err1 = background_proj[g][i][j][l]->GetBinError(k);
	    
	    for(int m = 1;  m<result[g][i][j][l]->GetNbinsX(); m++){
	      background[g][i][j][l]->SetBinContent(m,k,temp1);
	      background[g][i][j][l]->SetBinError(m,k,err1);
		
	    }

	  }

	  result[g][i][j][l]->Add(background[g][i][j][l],-1.);


	  if(i==3){
	    result[g][i][j][l]->Scale(1./4);
	  }
	  //-------------------------------
	  //dEta projection
	  //------------------------

	  TString eta_proj_name= in_name;
	  eta_proj_name.ReplaceAll("Yield","Eta_Proj");
	    
	  llimiteta = result[g][i][j][l]->GetXaxis()->FindBin(-etalim+.001);
	  rlimiteta = result[g][i][j][l]->GetXaxis()->FindBin(etalim-.001);

	  llimitphi = result[g][i][j][l]->GetYaxis()->FindBin(-philim+.001);
	  rlimitphi = result[g][i][j][l]->GetYaxis()->FindBin(philim-.001);
	    

	  eta_proj[g][i][j][l] = result[g][i][j][l]->ProjectionX(eta_proj_name,llimitphi,rlimitphi);
	  dx_eta = eta_proj[g][i][j][l]->GetBinWidth(1);
	  eta_proj[g][i][j][l]->Scale(1/dx_eta);

	  //	  eta_proj[g][i][j][l]->Scale(0.8);
	  

	  TString eta_proj_name_rebin = eta_proj_name;
	  eta_proj_name_rebin.ReplaceAll("Eta_Proj","Eta_Proj_Rebin");
	 
	  eta_proj_rebin[g][i][j][l] = (TH1D*)Rebin_dEta(eta_proj[g][i][j][l]);
	  eta_proj_rebin[g][i][j][l]->SetName(eta_proj_name_rebin);
	 

	  //-------------------------------
	  //dPhi projection
	  //------------------------

	  TString phi_proj_name= in_name;
	  phi_proj_name.ReplaceAll("Yield","Phi_Proj");

	  phi_proj[g][i][j][l] = result[g][i][j][l]->ProjectionY(phi_proj_name,llimiteta,rlimiteta);
	  dx_phi = phi_proj[g][i][j][l]->GetBinWidth(1);
	  phi_proj[g][i][j][l]->Scale(1/dx_phi);

	  TString phi_proj_name_rebin = phi_proj_name;
	  phi_proj_name_rebin.ReplaceAll("Phi_Proj","Phi_Proj_Rebin");
	 
	  phi_proj_rebin[g][i][j][l] = (TH1D*)Rebin_dPhi(phi_proj[g][i][j][l]);
	  phi_proj_rebin[g][i][j][l]->SetName(phi_proj_name_rebin);

	  //------------------------------------------------------------------
	  //      YIELD PLOTS:  First, settings in both dEta and dPhi together
	  //-------------------------------------------------------------------
	 
	  if(g==7||g==9||g==11){

	    //  cout<<"New_check_Eta_PYTH_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"rebin"<<endl;
	    eta_proj_ref[g][i][j][l] = (TH1D*)fin_ref[g]->Get((TString)("New_check_Eta_PYTH_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"rebin"))->Clone((TString)("New_check_Eta_PYTH_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"rebin"));
	    phi_proj_ref[g][i][j][l] = (TH1D*)fin_ref[g]->Get((TString)("New_check_Phi_PYTH_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"rebin"))->Clone((TString)("New_check_Phi_PYTH_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"rebin"));
	  }else{

	    eta_proj_ref[g][i][j][l] = (TH1D*)fin_ref[g]->Get((TString)("New_check_Eta_HYD_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"rebin"))->Clone((TString)("New_check_Eta_HYD_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"rebin"));
	    phi_proj_ref[g][i][j][l] = (TH1D*)fin_ref[g]->Get((TString)("New_check_Phi_HYD_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"rebin"))->Clone((TString)("New_check_Phi_HYD_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_Pt100_Pt300_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]+"rebin"));

	  }
		  
	  eta_proj_ref[g][i][j][l]->SetLineColor(kRed);
	  eta_proj_ref[g][i][j][l]->SetMarkerColor(kRed);

	  phi_proj_ref[g][i][j][l]->SetLineColor(kRed);
	  phi_proj_ref[g][i][j][l]->SetMarkerColor(kRed);
	  
	  if(i==3){
	    eta_proj_ref[g][i][j][l]->Scale(4.);
	    phi_proj_ref[g][i][j][l]->Scale(4.);
	  }
	  
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
	    check_ymax = 10.2; 
	    check_ymin = -1.; 
	   
	    break;
	  default: 
	    check_ymax = 10.; break;
	    check_ymin = -1.; break;
	  }
	  
	  /*

	  switch(i){
	  case 0: 
	    check_ymax = 30.;
	    check_ymin = -5.; 
	    break;
	  case 1: 
	    check_ymax = 20.; 
	    check_ymin = -2.; 
	    break;
	  case 2: 
	    check_ymax = 15.; 
	    check_ymin = -1.; 
	    break;
	  case 3: 
	    check_ymax =4.; 
	    check_ymin = -0.5; 
	   
	    break;
	  default: 
	    check_ymax = 10.; break;
	    check_ymin = -1.; break;
	  }
	  */

	  switch(g){
	  case 6: 
	    eta_proj_rebin[g][i][j][l]->SetMarkerStyle(10);
	    phi_proj_rebin[g][i][j][l]->SetMarkerStyle(10);
	    break;
	  case 7:
	    eta_proj_rebin[g][i][j][l]->SetMarkerStyle(4);
	    phi_proj_rebin[g][i][j][l]->SetMarkerStyle(4);
	    break;
	  case 8: 
	    eta_proj_rebin[g][i][j][l]->SetMarkerStyle(34);
	    eta_proj_rebin[g][i][j][l]->SetMarkerSize(2);
	    phi_proj_rebin[g][i][j][l]->SetMarkerStyle(34);
	    phi_proj_rebin[g][i][j][l]->SetMarkerSize(2);
	    break;
	  case 9:
	    eta_proj_rebin[g][i][j][l]->SetMarkerStyle(28);
	    phi_proj_rebin[g][i][j][l]->SetMarkerStyle(28);
	    break;  
	  case 10: 
	    eta_proj_rebin[g][i][j][l]->SetMarkerStyle(21);
	    phi_proj_rebin[g][i][j][l]->SetMarkerStyle(21);
	    break;
	  case 11:
	    eta_proj_rebin[g][i][j][l]->SetMarkerStyle(25);
	    phi_proj_rebin[g][i][j][l]->SetMarkerStyle(25);
	    break;
	  default:
	    break;
	  }
	 
	  eta_proj_rebin[g][i][j][l]->SetMarkerColor(kBlack);
	  eta_proj_rebin[g][i][j][l]->SetLineColor(kBlack);
	  
	  phi_proj_rebin[g][i][j][l]->SetMarkerColor(kBlack);
	  phi_proj_rebin[g][i][j][l]->SetLineColor(kBlack);
	
	  //------------------
	  // dEta plotting
	  //-------------------

	  eta_proj_rebin[g][i][j][l]->SetMaximum(check_ymax);
	  eta_proj_rebin[g][i][j][l]->SetMinimum(check_ymin);

   
	  eta_proj_rebin[g][i][j][l]->SetMarkerColor(1);
	  eta_proj_rebin[g][i][j][l]->SetAxisRange(-2.5+.2,2.5-.2,"x");

	  eta_proj_rebin[g][i][j][l]->Draw();	 
	  
//	  eta_proj_ref[g][i][j][l]->Draw("same");
      
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(tstitle);
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->SetTitle("#Delta#eta");
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitle("Y#equiv1/N_{jet} d^{2}N/(d#Delta#etadp_{T}) (GeV/c)^{-1}");
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);

 
	  eta_proj_rebin[g][i][j][l]->GetXaxis()->CenterTitle();
	  eta_proj_rebin[g][i][j][l]->GetYaxis()->CenterTitle();
	  if(i<3){
	    eta_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	  }
	  if(j>0){
	    eta_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	    eta_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	  }

	  eta_proj_rebin[g][i][j][l]->Draw("same");


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

	  eta_proj_rebin2[g][i][j][l] = (TH1D*)eta_proj_rebin[g][i][j][l]->Clone(narrowetaname);
	      
	  eta_proj_rebin2[g][i][j][l]->SetAxisRange(-1.5+.1,1.5-.1,"x");

	  eta_proj_rebin2[g][i][j][l]->Draw();	 
	  eta_proj_rebin[g][i][j][l]->Draw("same");

	  drawlabels_4by4(g,i,j,datalabel);
	 
	  lineEta = new TLine(-1.5,0,1.5,0);
	  lineEta->SetLineStyle(2);
	  lineEta->Draw("same");

	  //----------------------------
	  //   Now dPhi
	  //-------------------------------


	  ccheckPhi[g][l]->cd(4*i+j+1);

	
	  phi_proj_rebin[g][i][j][l]->SetMaximum(check_ymax);
	  phi_proj_rebin[g][i][j][l]->SetMinimum(check_ymin);
	    
	  phi_proj_rebin[g][i][j][l]->Draw("same");   
      
 
	  phi_proj_rebin[g][i][j][l]->GetXaxis()->CenterTitle();
	  phi_proj_rebin[g][i][j][l]->GetYaxis()->CenterTitle();
	  phi_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(tstitle);
	  phi_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	  phi_proj_rebin[g][i][j][l]->GetXaxis()->SetTitle("#Delta#phi");
	  phi_proj_rebin[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	  phi_proj_rebin[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	  phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitle("Y#equiv1/N_{jet} d^{2}N/(d#Delta#phidp_{T}) (GeV/c)^{-1}");
	  phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	  phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);
	  if(i<3){
	    phi_proj_rebin[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	  }
	  if(j>0){
	
	    phi_proj_rebin[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	    phi_proj_rebin[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	  }
	  
	  phi_proj_rebin[g][i][j][l]->Draw("same");

	  phi_proj_ref[g][i][j][l]->Draw("same");
   
	  drawlabels_4by4(g,i,j,datalabel);

	  linePhi = new TLine(-1.51,0,1.51,0);
	  linePhi->SetLineStyle(2);
	  linePhi->Draw("same");

	 
	
	  //---------------------------------------------------------------------------
	  // WRITE ALL PLOTS (associated with correlations & correlations syst. error)
	  //---------------------------------------------------------------------------

	  if(g==7||g==9||g==11){
	   	     
	    eta_proj_rebin[g-1][i][j][l]->Write();
	    phi_proj_rebin[g][i][j][l]->Write();
	    eta_proj_rebin[g-1][i][j][l]->Write();
	    phi_proj_rebin[g][i][j][l]->Write();

	    TString HminusPnameEta = in_name;
	    HminusPnameEta.ReplaceAll("Yield","Raw_HYD_minus_PYTH_Eta");
	     
	    TString gaus_eta_name = in_name;
	    gaus_eta_name.ReplaceAll("Yield","GausFit_Eta");

	    EtaClosureName = in_name;
	    EtaClosureName.ReplaceAll("Pt100_Pt300_","");
	    EtaClosureName.ReplaceAll("Yield","SpillOvers_Eta");
	    
	    TString HminusPnamePhi = HminusPnameEta;
	    HminusPnamePhi.ReplaceAll("Eta","Phi");
	   	     
	    TString gaus_phi_name = gaus_eta_name;
	    gaus_phi_name.ReplaceAll("Eta","Phi");


	    TString PhiClosureName = EtaClosureName;
	    PhiClosureName.ReplaceAll("Eta","Phi");
	  
	    
	    closure_eta_ref[g][i][j][l] = (TH1D*)fin_ref[g]->Get((TString)("Raw_HYD_minus_PYTH_Eta_Gen_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]));

	    closure_phi_ref[g][i][j][l] = (TH1D*)fin_ref[g]->Get((TString)("Raw_HYD_minus_PYTH_Phi_Gen_"+ CBin_strs[3-j]+"_"+CBin_strs[4-j]+"_"+TrkPtBin_strs[i]+"_"+TrkPtBin_strs[i+1]));

	    closure_eta_ref[g][i][j][l]->SetLineColor(kRed);
	    closure_eta_ref[g][i][j][l]->SetMarkerColor(kRed);


	    closure_phi_ref[g][i][j][l]->SetLineColor(kRed);
	    closure_phi_ref[g][i][j][l]->SetMarkerColor(kRed);

	    if(i==3){
	      closure_eta_ref[g][i][j][l]->Scale(4.);
	      closure_phi_ref[g][i][j][l]->Scale(4.);
	    }
	    //HYDJET minus PYTHIA eta

	    cHminusP_eta[g][l]->cd(4*i+j+1);

	   

	    HYDJET_PYTHIA_eta[g][i][j][l] = (TH1D*)eta_proj_rebin[g-1][i][j][l]->Clone(HminusPnameEta);
	    HYDJET_PYTHIA_eta[g][i][j][l]->SetName(HminusPnameEta);
	    
	    HYDJET_PYTHIA_eta[g][i][j][l]->Add(eta_proj_rebin[g][i][0][l],-1.);
	      
	  
	    HYDJET_PYTHIA_eta[g][i][j][l]->SetLineColor(1);
	    HYDJET_PYTHIA_eta[g][i][j][l]->SetMarkerStyle(10);
	    HYDJET_PYTHIA_eta[g][i][j][l]->SetMarkerColor(1);
	    HYDJET_PYTHIA_eta[g][i][j][l]->SetMarkerSize(1);


	    HYDJET_PYTHIA_eta[g][i][j][l]->SetMaximum(check_ymax);
	    HYDJET_PYTHIA_eta[g][i][j][l]->SetMinimum(check_ymin);
	    HYDJET_PYTHIA_eta[g][i][j][l]->Draw();     


	    HYDJET_PYTHIA_eta[g][i][j][l]->GetXaxis()->CenterTitle();
	    HYDJET_PYTHIA_eta[g][i][j][l]->GetYaxis()->CenterTitle();
	    HYDJET_PYTHIA_eta[g][i][j][l]->GetYaxis()->SetLabelSize(tstitle);
	    HYDJET_PYTHIA_eta[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	    HYDJET_PYTHIA_eta[g][i][j][l]->GetXaxis()->SetTitle("#Delta#eta");
	    HYDJET_PYTHIA_eta[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	    HYDJET_PYTHIA_eta[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	    HYDJET_PYTHIA_eta[g][i][j][l]->GetYaxis()->SetTitle("1/N_{jet} d^{2}N/(d#Delta#phi dp_{T})");
	    HYDJET_PYTHIA_eta[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	    HYDJET_PYTHIA_eta[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);
	    if(i<3){
	      HYDJET_PYTHIA_eta[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	    }
	    if(j>0){
	      HYDJET_PYTHIA_eta[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	      HYDJET_PYTHIA_eta[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	    }
	      
	    llimiteta = HYDJET_PYTHIA_eta[g][i][j][l]->GetXaxis()->FindBin(-1.0+.0001);
	    rlimiteta = HYDJET_PYTHIA_eta[g][i][j][l]->GetXaxis()->FindBin(1.0-.0001);
	     
	    double Yield_eta = HYDJET_PYTHIA_eta[g][i][j][l]->Integral(llimiteta,rlimiteta,"width");	      
	      
	     	     
	    
	      // if(Yield_eta<0.){Yield_eta=0.;}

	      
	    gaus1d->FixParameter(0,0.);
	    gaus1d->FixParameter(1,Yield_eta);
	      
	    gaus1d->ReleaseParameter(2);
	    gaus1d->SetParameter(2,0.2);
	    gaus1d->SetParLimits(2,0.1,0.4);
	    //  if(i>1){    gaus1d->SetParLimits(2,0.05,0.2);}
	   	    	      
	    HYDJET_PYTHIA_eta[g][i][j][l]->Fit("gaus1d","","",-1.,1.);
	      
	  
	    gaus_eta[g][i][j] = new TF1(gaus_eta_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

	    for(int k= 0; k<3; k++){
	
	      double temp = gaus1d->GetParameter(k);
	      gaus_eta[g][i][j]->SetParameter(k, temp);

	      //	cout<<g<<"  "<<i<<"  "<<j<<"  "<<l<<"  "<<"  "<<k<<"  "<<temp<<endl;	
	    }

	    double err_temp= gaus1d->GetParError(2);

	
	    HYDJET_PYTHIA_eta[g][i][j][l]->SetMaximum(check_ymax);
	    HYDJET_PYTHIA_eta[g][i][j][l]->SetMinimum(check_ymin);
	      
	    //    HYDJET_PYTHIA_eta[g][i][j][l]->SetMaximum(1.5);
	    HYDJET_PYTHIA_eta[g][i][j][l]->Draw();     

	 

	    closure_eta[g][i][j][l] = (TH1D*)HYDJET_PYTHIA_eta[g][i][j][l]->Clone(EtaClosureName);
   
	    for(int k=1; k< HYDJET_PYTHIA_eta[g][i][j][l]->GetNbinsX()+1; k++){
	      evalpt = HYDJET_PYTHIA_eta[g][i][j][l]->GetBinCenter(k);
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
	    //  closure_eta[g][i][j][l]->Draw("same");

	    //	    closure_eta_ref[g][i][j][l]->Draw("same");
	    
	    //    closure_eta_direct_rebin[g][i][j][l]->Draw("same");

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

	    lHminusP->AddEntry(HYDJET_PYTHIA_eta[g][i][j][l],"HYD.-PYTH.","lpfe");
	     	
	    if(i==3){lHminusP->SetTextSize(ts2);}
	    lHminusP->Draw("same");



	    //HYDJET minus PYTHIA phi


	    cHminusP_phi[g][l]->cd(4*i+j+1);

	    HYDJET_PYTHIA_phi[g][i][j][l] = (TH1D*)phi_proj_rebin[g-1][i][j][l]->Clone(HminusPnamePhi);


	    HYDJET_PYTHIA_phi[g][i][j][l]->SetName(HminusPnamePhi);

	    HYDJET_PYTHIA_phi[g][i][j][l]->Add(phi_proj_rebin[g][i][0][l],-1.);
	  
	 
	    HYDJET_PYTHIA_phi[g][i][j][l]->SetLineColor(1);
	    HYDJET_PYTHIA_phi[g][i][j][l]->SetMarkerStyle(10);
	    HYDJET_PYTHIA_phi[g][i][j][l]->SetMarkerColor(1);
	    HYDJET_PYTHIA_phi[g][i][j][l]->SetMarkerSize(1);


	    HYDJET_PYTHIA_phi[g][i][j][l]->SetMaximum(check_ymax);
	    HYDJET_PYTHIA_phi[g][i][j][l]->SetMinimum(check_ymin);
	    HYDJET_PYTHIA_phi[g][i][j][l]->Draw();     
      
 
	    HYDJET_PYTHIA_phi[g][i][j][l]->GetXaxis()->CenterTitle();
	    HYDJET_PYTHIA_phi[g][i][j][l]->GetYaxis()->CenterTitle();
	    HYDJET_PYTHIA_phi[g][i][j][l]->GetYaxis()->SetLabelSize(tstitle);
	    HYDJET_PYTHIA_phi[g][i][j][l]->GetXaxis()->SetLabelSize(tstitle);
	    HYDJET_PYTHIA_phi[g][i][j][l]->GetXaxis()->SetTitle("#Delta#phi");
	    HYDJET_PYTHIA_phi[g][i][j][l]->GetXaxis()->SetTitleSize(tstitle);
	    HYDJET_PYTHIA_phi[g][i][j][l]->GetXaxis()->SetTitleOffset(xoffset);
	    HYDJET_PYTHIA_phi[g][i][j][l]->GetYaxis()->SetTitle("1/N_{jet} d^{2}N/(d#Delta#phi dp_{T})");
	    HYDJET_PYTHIA_phi[g][i][j][l]->GetYaxis()->SetTitleOffset(yoffset);
	    HYDJET_PYTHIA_phi[g][i][j][l]->GetYaxis()->SetTitleSize(tstitle);
	    if(i<3){
	      HYDJET_PYTHIA_phi[g][i][j][l]->GetXaxis()->SetLabelSize(0.);
	    }
	    if(j>0){
	      HYDJET_PYTHIA_phi[g][i][j][l]->GetYaxis()->SetTitleSize(0.0);
	      HYDJET_PYTHIA_phi[g][i][j][l]->GetYaxis()->SetLabelSize(0.0);
	    }

	    llimitphi = HYDJET_PYTHIA_phi[g][i][j][l]->GetXaxis()->FindBin(-1.0+.01);
	    rlimitphi = HYDJET_PYTHIA_phi[g][i][j][l]->GetXaxis()->FindBin(1.0-.001);
	   	      
	 
	    double Yield_phi = HYDJET_PYTHIA_phi[g][i][j][l]->Integral(llimitphi,rlimitphi,"width");	      
	      
	    Yield_phi = Yield_eta;
	 	      
	    //	    if(Yield_phi<0){Yield_phi=0;}
	    

	    gaus1d->SetParLimits(2,0.1,0.5);

	    if((g==6||g==4)&&i==0&&j==3){
	      double width_temp = gaus_phi[g][0][2]->GetParameter(2);
	      gaus1d->FixParameter(2,width_temp);
	    }
	  	     	      
	    HYDJET_PYTHIA_phi[g][i][j][l]->Fit("gaus1d");
	    
	    gaus_phi[g][i][j] = new TF1(gaus_phi_name,"[0]+[1]/TMath::Sqrt(2*TMath::Pi())/[2]*TMath::Exp(-0.5*TMath::Power((TMath::Abs(x)/[2]),2.))",-2.,2.);

	    for(int k= 0; k<3; k++){
	      double temp = gaus1d->GetParameter(k);
	      gaus_phi[g][i][j]->SetParameter(k, temp);
	    }

	    err_temp = gaus1d->GetParError(2);
	 
	    HYDJET_PYTHIA_phi[g][i][j][l]->SetMaximum(check_ymax);
	    HYDJET_PYTHIA_phi[g][i][j][l]->SetMinimum(check_ymin);
	    HYDJET_PYTHIA_phi[g][i][j][l]->Draw();    

	    gaus_phi[g][i][j]->SetLineColor(kBlue);
	    gaus_phi[g][i][j]->Draw("same"); 


	    closure_phi[g][i][j][l] = (TH1D*)HYDJET_PYTHIA_phi[g][i][j][l]->Clone(PhiClosureName);
      
	    for(int k=0; k< HYDJET_PYTHIA_phi[g][i][j][l]->GetNbinsX()+1; k++){
	      evalpt = HYDJET_PYTHIA_phi[g][i][j][l]->GetBinCenter(k);
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

	    //	    closure_phi_ref[g][i][j][l]->Draw("same");

	    //  closure_phi_direct_rebin[g][i][j][l]->Draw("same");
	    

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
	    
	    

	    eta_proj_rebin[g-1][i][j][l]->Write();
	    phi_proj_rebin[g-1][i][j][l]->Write();
	   
	    eta_proj_rebin[g][i][j][l]->Write();
	    phi_proj_rebin[g][i][j][l]->Write();
	    
	    
	    HYDJET_PYTHIA_eta[g][i][j][l]->Write();
	    closure_eta[g][i][j][l]->Write();
	    gaus_eta[g][i][j]->Write();

	    HYDJET_PYTHIA_phi[g][i][j][l]->Write();
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
      
      ccheckPhi[g][l]->SaveAs("Phi_Correlations_"+datalabel+".png");
      ccheckEta[g][l]->SaveAs("Eta_Correlations_"+datalabel+".png");
      ccheckEta_wide[g][l]->SaveAs("Eta_Wide_Correlations_"+datalabel+".png");
      if((g==7||g==9||g==11)){
	cHminusP_eta[g][l]->SaveAs("Eta_SpillOvers_"+jettype2+"RecoJet_GenTrack.png");
	cHminusP_phi[g][l]->SaveAs("Phi_SpillOvers_"+jettype2+"RecoJet_GenTrack.png");

	cHminusP_eta[g][l]->SaveAs("Eta_SpillOvers_"+jettype2+"RecoJet_GenTrack.pdf");
	cHminusP_phi[g][l]->SaveAs("Phi_SpillOvers_"+jettype2+"RecoJet_GenTrack.pdf");
      }
      
      /*
      
      ccheckPhi[g][l]->SaveAs("Phi_Correlations_RightJets_"+datalabel+".png");
      ccheckEta[g][l]->SaveAs("Eta_Correlations_RightJets_"+datalabel+".png");
      ccheckEta_wide[g][l]->SaveAs("Eta_Wide_Correlations_RightJets_"+datalabel+".png");
      if((g==7||g==9||g==11)){
	cHminusP_eta[g][l]->SaveAs("Eta_SpillOvers_RightJets_"+datalabel+".png");
	cHminusP_phi[g][l]->SaveAs("Phi_SpillOvers_RightJets_"+datalabel+".png");
      }
      */

    }  // Closes the [l] (gen vs. reco) loop
 
  }// closes g loop



  cout<<"starting plots for AN"<<endl;
  //**********************
  //  CLOSURE PLOTS FOR AN
  //***********************

  TCanvas *cClosuresAN_eta[12][4],*cClosuresAN_phi[12][4];
  float result_min, result_max,  val_l,val_r,err_l,err_r;
  TLegend *l40,*l41,*l42;

 for(int g=6;g<12;g++){
      
    switch(g){
    case 6: 
      continue; 	break;
    case 7: 
      datalabel = "Inclusive"; break;
    case 8: 
      continue; break;
    case 9: 
      datalabel = "Subleading"; break;
    case 10: 
      continue; break;
    case 11:
      datalabel = "Leading"; break;
    default:
      continue; 
      break;
    }

    cout<<"starting "<<datalabel<<endl;

    for(int i = 0; i<4; i++){

 
      TString pTrange;

      switch(i){
      case 0: 
	pTrange = "TrkPt1_TrkPt2"; 
	result_max = 11.;
	result_min = -1.;
	break;
      case 1: 
	pTrange = "TrkPt2_TrkPt3";
	result_max = 8.5;
	result_min = -1.; 
	break;
      case 2: 
	pTrange = "TrkPt3_TrkPt4"; 
	result_max = 6.2;
	result_min = -.45;
	break;
      case 3: 
	pTrange = "TrkPt4_TrkPt8"; 
	result_max = 4.2;
	result_min = -.25;
	break;
      }
   


      TString ClosureANeta_name = "AN_Closures_Eta_";
      ClosureANeta_name+=datalabel;
      ClosureANeta_name+=pTrange;
      cClosuresAN_eta[g][i] = new TCanvas(ClosureANeta_name," ",10,10,1500,800);
      cClosuresAN_eta[g][i]->Divide(4,2,0.,0.);
   


      TString ClosureANphi_name = "AN_Closures_Phi_";
      ClosureANphi_name+=datalabel;
      ClosureANphi_name+=pTrange;
      cClosuresAN_phi[g][i] = new TCanvas(ClosureANphi_name," ",10,10,1500,800);
      cClosuresAN_phi[g][i]->Divide(4,2,0.,0.);
	

      for(int j=0; j<4; j++){
	
	cClosuresAN_eta[g][i]->cd(j+1);

	eta_proj_rebin[g][i][j][0]->SetMaximum(result_max);
	eta_proj_rebin[g][i][j][0]->SetMinimum(result_min);
	eta_proj_rebin[g][i][j][0]->Draw();
	eta_proj_rebin[g-1][i][j][0]->Draw("same");

     
	drawlabels(g,i,j);

	if(j==3){
	  TLatex *tex24eta = new TLatex(textalign,texty2,phirangelabel);
	  tex24eta->SetName("tex24eta");
	  tex24eta->SetNDC();
	  tex24eta->SetTextSizePixels(tspixels);
	  tex24eta->Draw();
	}

	lineEta = new TLine(-1.5,0,1.5,0); 
	lineEta->Draw("same");
	
	
	if(j==0){ 

	  l40 = new TLegend(textalign2,texty2-.05,0.9,texty4-0.05);
	  l40->SetName("l40");
	  l40->SetTextSizePixels(tspixels);
	
	  l40->SetLineColor(kWhite);
	  if(g==7){ 
	    l40->AddEntry(eta_proj_rebin[g][i][j][0],"PYTHIA Inclusive","lpfe");
	    l40->AddEntry(eta_proj_rebin[g-1][i][0][0],"HYDJET Inclusive","lpfe");
	  }
	  if(g==11){ 
	    l40->AddEntry(eta_proj_rebin[g][i][j][0],"PYTHIA Leading","lpfe");
	    l40->AddEntry(eta_proj_rebin[g-1][i][0][0],"HYDJET Leading","lpfe");
	  }
	  if(g==9){ 
	    l40->AddEntry(eta_proj_rebin[g][i][j][0],"PYTHIA Subleading","lpfe");
	    l40->AddEntry(eta_proj_rebin[g-1][i][0][0],"HYDJET Subleading","lpfe");
	  }	
	  l40->SetTextSizePixels(tspixels);
	  l40->Draw("same");
	}


	cClosuresAN_eta[g][i]->cd(j+5);
	HYDJET_PYTHIA_eta[g][i][j][0]->SetMinimum(result_min);
	HYDJET_PYTHIA_eta[g][i][j][0]->SetMaximum(result_max);
	HYDJET_PYTHIA_eta[g][i][j][0]->GetXaxis()->SetLabelSize(ts2);
	HYDJET_PYTHIA_eta[g][i][j][0]->GetYaxis()->SetTitleSize(tstitle2);
	HYDJET_PYTHIA_eta[g][i][j][0]->GetYaxis()->SetTitleOffset(1.);
	HYDJET_PYTHIA_eta[g][i][j][0]->GetYaxis()->SetTitle("HYDJET-PYTHIA");
	HYDJET_PYTHIA_eta[g][i][j][0]->GetYaxis()->SetTitleSize(0);
	
	if(j==0){ HYDJET_PYTHIA_eta[g][i][j][0]->GetXaxis()->SetLabelSize(ts2-0.01);
	  HYDJET_PYTHIA_eta[g][i][j][0]->GetXaxis()->SetTitleSize(ts);
	  HYDJET_PYTHIA_eta[g][i][j][0]->GetYaxis()->SetLabelSize(ts2);
	  HYDJET_PYTHIA_eta[g][i][j][0]->GetXaxis()->SetLabelOffset(0.013);
	}
	HYDJET_PYTHIA_eta[g][i][j][0]->SetAxisRange(-1.4,1.41,"x");
	HYDJET_PYTHIA_eta[g][i][j][0]->Draw();
	gaus_eta[g][i][j]->Draw("same");


	 
	lineEta->Draw("same");

	    
	 	 
	//----------------------------------
	// Closure Plots for AN (dPhi)
	//----------------------------------

	cClosuresAN_phi[g][i]->cd(j+1);

	phi_proj_rebin[g][i][j][0]->SetMaximum(result_max);
	phi_proj_rebin[g][i][j][0]->SetMinimum(result_min);
	phi_proj_rebin[g][i][j][0]->Draw();
	phi_proj_rebin[g-1][i][j][0]->Draw("same");
     
	drawlabels(g,i,j);

	if(j==3){
	  TLatex *tex24phi = new TLatex(textalign,texty2,etarangelabel);
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


	cClosuresAN_phi[g][i]->cd(j+5);
	  
	HYDJET_PYTHIA_phi[g][i][j][0]->SetMinimum(-1.);
	HYDJET_PYTHIA_phi[g][i][j][0]->SetMaximum(result_max);
	HYDJET_PYTHIA_phi[g][i][j][0]->GetXaxis()->SetLabelSize(ts2);
	HYDJET_PYTHIA_phi[g][i][j][0]->GetYaxis()->SetTitleSize(tstitle2);
	HYDJET_PYTHIA_phi[g][i][j][0]->GetYaxis()->SetTitleOffset(1.);
	HYDJET_PYTHIA_phi[g][i][j][0]->GetYaxis()->SetTitle("HYDJET-PYTHIA");
	HYDJET_PYTHIA_phi[g][i][j][0]->GetYaxis()->SetTitleSize(0.);
	  

	if(j==0){ HYDJET_PYTHIA_phi[g][i][j][0]->GetXaxis()->SetLabelSize(ts2-0.01);
	  HYDJET_PYTHIA_phi[g][i][j][0]->GetXaxis()->SetTitleSize(ts);
	  HYDJET_PYTHIA_phi[g][i][j][0]->GetYaxis()->SetLabelSize(ts2);
	  HYDJET_PYTHIA_phi[g][i][j][0]->GetXaxis()->SetLabelOffset(0.013);
	}

	HYDJET_PYTHIA_phi[g][i][j][0]->SetAxisRange(-1.4,1.41,"x");
	HYDJET_PYTHIA_phi[g][i][j][0]->GetXaxis()->SetLabelSize(ts2);
	if(j==0){HYDJET_PYTHIA_phi[g][i][j][0]->GetXaxis()->SetLabelSize(ts2-0.01);}
	HYDJET_PYTHIA_phi[g][i][j][0]->Draw();
	gaus_phi[g][i][j]->Draw("same");


	 
	linePhi->Draw("same");
	    
	 


      }//close j


      cClosuresAN_eta[g][i]->cd(0);
 
      TLatex *canvas_title = new TLatex(0.06,0.95,"CMS Preliminary Simulation");
      canvas_title->SetTextSizePixels(tspixels);
      canvas_title->SetTextFont(63);
      canvas_title->Draw();

      TLatex *canvas_title2 = new TLatex(0.295,0.95,"PYTHIA+HYDJET");
      canvas_title2->SetTextSizePixels(tspixels);
      canvas_title2->Draw();

      ClosureANeta_name+=".pdf";
      cClosuresAN_eta[g][i]->SaveAs(ClosureANeta_name);
      ClosureANeta_name.ReplaceAll("pdf","png");
      cClosuresAN_eta[g][i]->SaveAs(ClosureANeta_name);


      cClosuresAN_phi[g][i]->cd(0);
 
      canvas_title->Draw();
      canvas_title2->Draw();

      ClosureANphi_name+=".pdf";
      cClosuresAN_phi[g][i]->SaveAs(ClosureANphi_name);
      ClosureANphi_name.ReplaceAll("pdf","png");
      cClosuresAN_phi[g][i]->SaveAs(ClosureANphi_name);
    

    }//close i 

  } //close g (MC for closure plots
 

 
  return 0;


}  //and we're done.
 
