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
#include "TMultiGraph.h"
#include "TStyle.h"
#include "TLatex.h"


#include <iostream>
#include <vector>
#include <fstream>

#include "../HIN_14_016_functions.h"


Int_t me_correct3(float etacut = 1.6, bool is_pp = kFALSE){

  //-------------Input files & set eta-----------------

  TString jetetacut, etalabel,datalabel;
  float eta_ymax;

  TLatex *centtex, *pttex;

  TF1 *pol0 = new TF1("pol0","[0]+x-x",-.3,.3);


  TFile *fin, *fout_inc, *fout_lead, *fout_sub, *fout_non, *fout_inc_pTweighted, *fout_lead_pTweighted, *fout_sub_pTweighted, *fout_non_pTweighted,  *fout_inc_pTcos, *fout_lead_pTcos, *fout_sub_pTcos, *fout_non_pTcos, *finbg;

  TString desc = "Data";
  
  if( fabs(etacut-1.6)<1e-4&& is_pp){
    fin = new TFile("../Data_Raw_Correlations/Data_pp_June24.root","READ");
    fout_inc = new TFile("pp_Inclusive_Correlations.root","RECREATE");
    fout_lead = new TFile("pp_Leading_Correlations.root","RECREATE");
    fout_sub = new TFile("pp_SubLeading_Correlations.root","RECREATE");
    //    fout_non = new TFile("pp_NonLeading_Correlations.root","RECREATE");
    jetetacut = "JetEtaCut1.6";
    etalabel = "|#eta_{jet}|<1.6";
    datalabel = "pp";
    eta_ymax = 0.5;
    
  }else if( fabs(etacut-1.6)<1e-4&& !is_pp){
    fin = new TFile("../Data_Raw_Correlations/Data2011_PbPb_June15.root","READ");
    fout_inc = new TFile("PbPb_Inclusive_Correlations.root","RECREATE");
    fout_lead = new TFile("PbPb_Leading_Correlations.root","RECREATE");
    fout_sub = new TFile("PbPb_SubLeading_Correlations.root","RECREATE");
    //   fout_non = new TFile("PbPb_NonLeading_Correlations.root","RECREATE");
    jetetacut = "JetEtaCut1.6";
    etalabel = "|#eta_{jet}|<1.6";
    datalabel = "PbPb";
    eta_ymax = 0.5;

  }else {
    cout<<etacut<<endl;
    cerr<<"No data exists for that jet cut range."<<endl;
    return -1;
  }

  //----------------------------------------------------



  TCanvas *me_proj_canvas = new TCanvas("me_proj_canvas","",0,0,1500,1600);
  me_proj_canvas->Divide(4,4,0.0000,0.0000);




  TCanvas *me_proj_leading_canvas = new TCanvas("me_proj_leading_canvas","",0,0,1500,1600);
  me_proj_leading_canvas->Divide(4,4,0.0000,0.0000);





  TCanvas *me_proj_sub_canvas = new TCanvas("me_proj_sub_canvas","",0,0,1500,1600);
  me_proj_sub_canvas->Divide(4,4,0.0000,0.0000);




  TCanvas *result_proj_canvas = new TCanvas("result_proj_canvas","",0,0,1500,1600);
  result_proj_canvas->Divide(4,4,0.0000,0.0000);



  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins;
  
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 5;

  float me00,me00_lead,me00_sub;
  int me00binl, me00binr, n_me00bins;

  float PtBins[nPtBins+1] = {100, 300};
  TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};
  

  float CBins[nCBins+1] = {0, 20, 60, 100, 200};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%","Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };
  TString TrkPtBin_labels[nTrkPtBins] = {"1<pT<2","2<pT<3","3<pT<4","4<pT<8","pT>8"};
  
  TString ref_name;
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
  
  TH1F* all_jets_corrpT[nCBins][nPtBins];
  TH1F* only_leadingjets_corrpT[nCBins][nPtBins];
  TH1F* only_subleadingjets_corrpT[nCBins][nPtBins];
  TH1F* only_nonleadingjets_corrpT[nCBins][nPtBins];
   
  TH2D* hJetTrackSignalBackground[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundLeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundSubLeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundNonLeading[nCBins][nPtBins][nTrkPtBins];
 
  TH2D* hJetTrackME[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackMELeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackMESubLeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackMENonLeading[nCBins][nPtBins][nTrkPtBins];

  TH2D *yield_inc[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_inc_pTweighted[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_inc_pTcos[nCBins][nPtBins][nTrkPtBins];
  TH1D *yield_inc_proj[nCBins][nPtBins][nTrkPtBins];
  TH1D *test_me_inc_proj[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_lead[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_lead_pTweighted[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_lead_pTcos[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_sub[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_sub_pTweighted[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_sub_pTcos[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_non[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_non_pTweighted[nCBins][nPtBins][nTrkPtBins];
  TH2D *yield_non_pTcos[nCBins][nPtBins][nTrkPtBins];
 
  TH1D *me_proj[nCBins][nPtBins][nTrkPtBins];
  TH1D *me_proj_lead[nCBins][nPtBins][nTrkPtBins];
  TH1D *me_proj_sub[nCBins][nPtBins][nTrkPtBins];
  TH1D *me_proj_non[nCBins][nPtBins][nTrkPtBins];
  


  TH2D *raw_ref_inc[nCBins][nPtBins][nTrkPtBins];
  TH2D *raw_ref_sub[nCBins][nPtBins][nTrkPtBins];
  TH2D *raw_ref_lead[nCBins][nPtBins][nTrkPtBins];


  TH1D *raw_ref_proj_inc[nCBins][nPtBins][nTrkPtBins];
  TH1D *raw_ref_proj_sub[nCBins][nPtBins][nTrkPtBins];
  TH1D *raw_ref_proj_lead[nCBins][nPtBins][nTrkPtBins];


  TH2D *new_ref_inc[nCBins][nPtBins][nTrkPtBins];
  TH2D *new_ref_sub[nCBins][nPtBins][nTrkPtBins];
  TH2D *new_ref_lead[nCBins][nPtBins][nTrkPtBins];


  TH1D *new_ref_proj_inc[nCBins][nPtBins][nTrkPtBins];
  TH1D *new_ref_proj_sub[nCBins][nPtBins][nTrkPtBins];
  TH1D *new_ref_proj_lead[nCBins][nPtBins][nTrkPtBins];



  TH1D *raw_ratio_inc[nCBins][nPtBins][nTrkPtBins];
  TH1D *raw_ratio_sub[nCBins][nPtBins][nTrkPtBins];
  TH1D *raw_ratio_lead[nCBins][nPtBins][nTrkPtBins];

 

  float norm_temp, width_temp, width_temp_x, width_temp_y, width_temp_ref, max_bin, max_cont,bc,err;
    
  //-----------------------
  // Start getting histos
  //-----------------------

  
  for (int ibin=0;ibin<nCBins;ibin++){
  
    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 

      cout<<desc + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]<<endl;

         
      all_jets_corrpT[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("all_jets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));
    
      only_leadingjets_corrpT[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_only_leadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("only_leadingjets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));
 
      only_subleadingjets_corrpT[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_only_subleadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("only_subleadingjets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));
      
      only_nonleadingjets_corrpT[ibin][ibin2] = (TH1F*)fin->Get((TString) (desc + "_only_nonleadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]))->Clone((TString) ("only_nonleadingjets_corrpT_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]));
       
      
      for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){

	cout<<desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]<<endl;
	     
	hJetTrackSignalBackground[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundSubLeading"+  CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	hJetTrackSignalBackgroundNonLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackSignalBackgroundNonLeading"+  CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Raw_Yield_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	hJetTrackME[ibin][ibin2][ibin3] = (TH2D*) fin->Get((TString)(desc + "_hJetTrackME"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));
	
	hJetTrackMELeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	hJetTrackMESubLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackMESubLeading"+  CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));

	hJetTrackMENonLeading[ibin][ibin2][ibin3] = (TH2D*)fin->Get((TString) (desc + "_hJetTrackMENonLeading"+  CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) ("Mixed_Event_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_" +TrkPtBin_strs[ibin3+1]));


      } /// ibin3

    } // ibin2

    //-------------------------
    // Normalize all!
    //------------------------


    cout<<"got hists for centbin "<<ibin<<endl;
   
    for (int ibin2=0;ibin2<nPtBins;ibin2++){

      for(int ibin3 = 0; ibin3<4; ibin3++){

	width_temp_x = 1.;
	width_temp_y = 1.;

	norm_temp = all_jets_corrpT[ibin][ibin2]->GetEntries();
	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackME[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);
   
	norm_temp = only_leadingjets_corrpT[ibin][ibin2]->GetEntries();
	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);
   
	norm_temp = only_subleadingjets_corrpT[ibin][ibin2]->GetEntries();
	hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackMESubLeading[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);

	norm_temp = only_nonleadingjets_corrpT[ibin][ibin2]->GetEntries();
	hJetTrackSignalBackgroundNonLeading[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y);
	hJetTrackMENonLeading[ibin][ibin2][ibin3]->Scale(1/norm_temp/width_temp_x/width_temp_y/50.);
	
	//     INCLUSIVE

	fout_inc->cd();
	me_proj_canvas->cd(4*(ibin3+1)-ibin);

	yield_inc[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackground[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	me_proj[ibin][ibin2][ibin3] = hJetTrackME[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",ibin,ibin2,ibin3),1,100);

	me_proj[ibin][ibin2][ibin3]->Scale (1./100);
	me_proj[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);

	me00 = 	me_proj[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);

	//	me_proj[ibin][ibin2][ibin3]->Scale(1./me00);

	cout<<me00<<endl;
	hJetTrackME[ibin][ibin2][ibin3]->Scale(1./me00);

	me_proj[ibin][ibin2][ibin3]->SetAxisRange(-.299,.299);


	
	me_proj[ibin][ibin2][ibin3]->Draw();


	if(ibin<3){
	  TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}
	if(ibin==3){
	  TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}

	if(ibin==3&&ibin3==3)	yield_inc[ibin][ibin2][ibin3]->Divide(hJetTrackME[ibin-1][ibin2][ibin3]);
	else	yield_inc[ibin][ibin2][ibin3]->Divide(hJetTrackME[ibin][ibin2][ibin3]);

	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Write();
	yield_inc[ibin][ibin2][ibin3]->Write();
	hJetTrackME[ibin][ibin2][ibin3]->Write();


  
	//  LEADING

	fout_lead->cd();
	me_proj_leading_canvas->cd(4*(ibin3+1)-ibin);


	yield_lead[ibin][ibin2][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	me_proj_lead[ibin][ibin2][ibin3] = hJetTrackMELeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_lead_temp_%d%d%d",ibin,ibin2,ibin3),1,100);

	me_proj_lead[ibin][ibin2][ibin3]->Scale (1./100);
	me_proj_lead[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);

	me00_lead = 	me_proj_lead[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);


	//   SUBLEADING

	fout_sub->cd();
	me_proj_sub_canvas->cd(4*(ibin3+1)-ibin);

	yield_sub[ibin][ibin2][ibin3] = (TH2D*) hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Clone((TString)("Yield_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	//	yield_sub_pTweighted[ibin][ibin2][ibin3] = (TH2D*) hJetTrackSignalBackgroundSubLeading_pTweighted[ibin][ibin2][ibin3]->Clone((TString)("Yield_pTweighted_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));

	//	yield_sub_pTcos[ibin][ibin2][ibin3] = (TH2D*) hJetTrackSignalBackgroundSubLeading_pTcos[ibin][ibin2][ibin3]->Clone((TString)("Yield_pTcos_"+datalabel+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+"_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));


	
	me_proj_sub[ibin][ibin2][ibin3] = hJetTrackMESubLeading[ibin][ibin2][ibin3]->ProjectionX(Form("me_proj_sub_%d%d%d",ibin,ibin2,ibin3),1,100);

	me_proj_sub[ibin][ibin2][ibin3]->Scale (1./100);
	me_proj_sub[ibin][ibin2][ibin3] -> Fit("pol0","","",-.2,.2);

	me00_sub = me_proj_sub[ibin][ibin2][ibin3]->GetFunction("pol0")->GetParameter(0);

	me00  = (me00_sub+me00_lead)/2.;


	hJetTrackMESubLeading[ibin][ibin2][ibin3]->Scale(1./me00);

	me_proj_sub[ibin][ibin2][ibin3]->SetAxisRange(-.299,.299);
	me_proj_sub[ibin][ibin2][ibin3]->Draw();


	if(ibin<3){
	  TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}
	if(ibin==3){
	  TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}


	yield_sub[ibin][ibin2][ibin3]->Divide(hJetTrackMESubLeading[ibin][ibin2][ibin3]);

	hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Write();
	yield_sub[ibin][ibin2][ibin3]->Write();
	hJetTrackMESubLeading[ibin][ibin2][ibin3]->Write();
	
	fout_lead->cd();
	me_proj_leading_canvas->cd(4*(ibin3+1)-ibin);




	hJetTrackMELeading[ibin][ibin2][ibin3]->Scale(1./me00);

	//	me_proj_lead[ibin][ibin2][ibin3]->Scale(1./me00);
	me_proj_lead[ibin][ibin2][ibin3]->SetAxisRange(-.299,.299);



	me_proj_lead[ibin][ibin2][ibin3]->Draw();


	if(ibin<3){
	  TLatex *centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}
	if(ibin==3){
	  TLatex *centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	  centtex->SetNDC();
	  centtex->Draw();
	  TLatex *pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	  pttex->SetNDC();
	  pttex->Draw();
	}

	yield_lead[ibin][ibin2][ibin3]->Divide(hJetTrackMELeading[ibin][ibin2][ibin3]);

	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Write();
	yield_lead[ibin][ibin2][ibin3] ->Write();
	hJetTrackMELeading[ibin][ibin2][ibin3]->Write();





      }//ibin3;
      
    } // ibin2

  } // ibin ( centrality ) loop

  

  me_proj_canvas->SaveAs((TString)("All_ME_Projections_Inclusive_"+datalabel+".png"));
  me_proj_leading_canvas->SaveAs((TString)("All_ME_Projections_Leading_"+datalabel+".png"));
  me_proj_sub_canvas->SaveAs((TString)("All_ME_Projections_SubLeading_"+datalabel+".png"));
   
  return 0;
} // main loop






  
     
