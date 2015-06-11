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


Int_t me_correct_mc2(bool is_pp = kFALSE, int mc_type_code = 2){

 
  const int nCBins = 4;
  const int nPtBins = 1;
  const int nTrkPtBins = 5;




  enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RightGen, RightReco, n_data_mc_types};

  TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RightGenJet_GenTrack","RightRecoJet_GenTrack"};

  int data_mc_type_code = -999;

  float CBins[nCBins+1] = {0, 10, 30, 50, 100};
  TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
  TString CBin_labels[nCBins] = {"Cent. 0-10%", "Cent. 10-30%", "Cent. 30-50%","Cent. 50-100%"};

  float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
  TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };
  TString TrkPtBin_labels[nTrkPtBins] = {"1<p_{T}^assoc.<2", "2<p_{T}^assoc.<3","3<p_{T}^assoc.<4","p_{T}^assoc.>8"};


  TCanvas *me_proj_canvas[5], *me_proj_sub_canvas[5],*me_proj_lead_canvas[5], *result_proj_canvas[5], *result_proj_sub_canvas[5],*result_proj_lead_canvas[5];

  int llimitphi,rlimitphi,llimiteta,rlimiteta,nbins;
 
  float me00;
  int me00bin;

  TString ref_name;
 
  gStyle->SetOptStat(0);  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetPadRightMargin (0.05);
  
  TH2D* hJetTrackSignalBackground[5][nCBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundLeading[5][nCBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundSubLeading[5][nCBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundNonLeading[5][nCBins][nTrkPtBins];

  TH2D* hJetTrackME[5][nCBins][nTrkPtBins];
  TH2D* hJetTrackMELeading[5][nCBins][nTrkPtBins];
  TH2D* hJetTrackMESubLeading[5][nCBins][nTrkPtBins];
  TH2D* hJetTrackMENonLeading[5][nCBins][nTrkPtBins];
TH2D *yield_inc[5][nCBins][nTrkPtBins];
  TH2D *yield_lead[5][nCBins][nTrkPtBins];
  TH2D *yield_sub[5][nCBins][nTrkPtBins];

  TH1D *yield_inc_proj[5][nCBins][nTrkPtBins];
  TH1D *yield_lead_proj[5][nCBins][nTrkPtBins];
  TH1D *yield_sub_proj[5][nCBins][nTrkPtBins];
   
  TH1D *me_proj[5][nCBins][nTrkPtBins];
  TH1D *me_proj_lead[5][nCBins][nTrkPtBins];
  TH1D *me_proj_sub[5][nCBins][nTrkPtBins];
 
  //-------------Input files & set eta-----------------

  TString jetetacut, etalabel,datalabel;
  float eta_ymax;

  TLatex *centtex, *pttex;

  TFile *fin[n_data_mc_types], *fout_inc[n_data_mc_types], *fout_lead[n_data_mc_types], *fout_sub[n_data_mc_types], *fout_non[n_data_mc_types],*fin_me[n_data_mc_types];
    
 
  if( !is_pp ){

    cout<<"../MC_Raw_Correlations/PythiaHydjet_JFFCorr_Merged_"+data_mc_type_strs[mc_type_code]+".root"<<endl;
    fin[mc_type_code] = new TFile("../MC_Raw_Correlations/PythiaHydjet_JFFCorr_Merged_"+data_mc_type_strs[mc_type_code]+".root","READ");
   
    if(mc_type_code%2!=0){
      fin_me[mc_type_code] = new TFile("../MC_Raw_Correlations/PythiaHydjet_JFFCorr_Merged_"+data_mc_type_strs[mc_type_code+1]+".root","READ");
    }

    fout_inc[mc_type_code] = new TFile("HydJet_"+data_mc_type_strs[mc_type_code]+"_Inclusive_Correlations.root","RECREATE");
    fout_lead[mc_type_code] = new TFile("HydJet_"+data_mc_type_strs[mc_type_code]+"_Leading_Correlations.root","RECREATE");
    fout_sub[mc_type_code] = new TFile("HydJet_"+data_mc_type_strs[mc_type_code]+"_SubLeading_Correlations.root","RECREATE");
  } else {
    fin[mc_type_code] = new TFile("../MC_Raw_Correlations/Pythia_JFFCorr_Merged_"+data_mc_type_strs[mc_type_code]+".root","READ");
    
    if(mc_type_code%2!=0){
      fin_me[mc_type_code] = new TFile("../MC_Raw_Correlations/Pythia_JFFCorr_Merged_"+data_mc_type_strs[mc_type_code+1]+".root","READ");
    }

    fout_inc[mc_type_code] = new TFile("Pythia_"+data_mc_type_strs[mc_type_code]+"_Inclusive_Correlations.root","RECREATE");
    fout_lead[mc_type_code] = new TFile("Pythia_"+data_mc_type_strs[mc_type_code]+"_Leading_Correlations.root","RECREATE");
    fout_sub[mc_type_code] = new TFile("Pythia_"+data_mc_type_strs[mc_type_code]+"_SubLeading_Correlations.root","RECREATE");
  }
  cout<<"opened files"<<endl;
  //-----------------------
  // Start getting histos
  //-----------------------

  
  for (int ibin=0;ibin<nCBins;ibin++){

  
    for (int ibin3=0;ibin3<nTrkPtBins-1;ibin3++){

      hJetTrackSignalBackground[mc_type_code][ibin][ibin3]=  (TH2D*)fin[mc_type_code]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackSignalBackground_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackSignalBackground_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
    
      hJetTrackSignalBackgroundLeading[mc_type_code][ibin][ibin3]= (TH2D*)fin[mc_type_code]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackSignalBackgroundLeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackSignalBackgroundLeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

      hJetTrackSignalBackgroundSubLeading[mc_type_code][ibin][ibin3]= (TH2D*)fin[mc_type_code]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackSignalBackgroundSubLeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString)(data_mc_type_strs[mc_type_code] + "_hJetTrackSignalBackgroundSubLeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));


      if(mc_type_code%2==0){

	cout<< (data_mc_type_strs[mc_type_code] + "_hJetTrackME_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1])<<endl;
	hJetTrackME[mc_type_code][ibin][ibin3]=  (TH2D*)fin[mc_type_code]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackME_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackME_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
    
	hJetTrackMELeading[mc_type_code][ibin][ibin3]= (TH2D*)fin[mc_type_code]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	hJetTrackMESubLeading[mc_type_code][ibin][ibin3]= (TH2D*)fin[mc_type_code]->Get((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString)(data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
 
	cout<<"here"<<endl;
      }else{
	hJetTrackME[mc_type_code][ibin][ibin3]=  (TH2D*)fin_me[mc_type_code]->Get((TString) (data_mc_type_strs[mc_type_code+1] + "_hJetTrackME_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackME_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
    
	hJetTrackMELeading[mc_type_code][ibin][ibin3]= (TH2D*)fin_me[mc_type_code]->Get((TString) (data_mc_type_strs[mc_type_code+1] + "_hJetTrackMELeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString) (data_mc_type_strs[mc_type_code] + "_hJetTrackMELeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));

	hJetTrackMESubLeading[mc_type_code][ibin][ibin3]= (TH2D*)fin_me[mc_type_code]->Get((TString) (data_mc_type_strs[mc_type_code+1] + "_hJetTrackMESubLeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]))->Clone((TString)(data_mc_type_strs[mc_type_code] + "_hJetTrackMESubLeading_Merged_"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]));
      }
    } /// ibin3

  } // ibin

  cout<<"here"<<endl;

  //-------------------------
  // Do ME correction
  //------------------------

  me_proj_canvas[mc_type_code] = new TCanvas((TString)("me_proj_canvas_"+data_mc_type_strs[mc_type_code]),"",0,0,1500,1600);
  me_proj_canvas[mc_type_code]->Divide(4,4,0.0000,0.0000);

  me_proj_sub_canvas[mc_type_code] = new TCanvas((TString)("me_proj_sub_canvas_"+data_mc_type_strs[mc_type_code]),"",0,0,1500,1600);
  me_proj_sub_canvas[mc_type_code]->Divide(4,4,0.0000,0.0000);

  me_proj_lead_canvas[mc_type_code] = new TCanvas((TString)("me_proj_lead_canvas_"+data_mc_type_strs[mc_type_code]),"",0,0,1500,1600);
  me_proj_lead_canvas[mc_type_code]->Divide(4,4,0.0000,0.0000);

  result_proj_canvas[mc_type_code] = new TCanvas((TString)("result_proj_canvas_"+data_mc_type_strs[mc_type_code]),"",0,0,1500,1600);
  result_proj_canvas[mc_type_code]->Divide(4,4,0.0000,0.0000);

  result_proj_lead_canvas[mc_type_code] = new TCanvas((TString)("result_proj_lead_canvas_"+data_mc_type_strs[mc_type_code]),"",0,0,1500,1600);
  result_proj_lead_canvas[mc_type_code]->Divide(4,4,0.0000,0.0000);

  result_proj_sub_canvas[mc_type_code] = new TCanvas((TString)("result_proj_sub_canvas_"+data_mc_type_strs[mc_type_code]),"",0,0,1500,1600);
  result_proj_sub_canvas[mc_type_code]->Divide(4,4,0.0000,0.0000);

  
  for (int ibin=0;ibin<nCBins;ibin++){

    for(int ibin3 = 0; ibin3<4; ibin3++){

      cout<<"I am here: "<<mc_type_code<<" "<<ibin<<" "<<ibin3<<endl;

      fout_inc[mc_type_code]->cd();

      me_proj_canvas[mc_type_code]->cd(4*(ibin3+1)-ibin);

      yield_inc[mc_type_code][ibin][ibin3] =  (TH2D*) hJetTrackSignalBackground[mc_type_code][ibin][ibin3]->Clone((TString)("Yield_"+data_mc_type_strs[mc_type_code]+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" +TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));
    
      me_proj[mc_type_code][ibin][ibin3] = (TH1D*) hJetTrackME[mc_type_code][ibin][ibin3]->ProjectionX(Form("me_proj_temp_%d%d%d",mc_type_code,ibin,ibin3),1,100);

      me_proj[mc_type_code][ibin][ibin3]->Scale (1./100);
      me_proj[mc_type_code][ibin][ibin3] -> Fit("pol0","","",-.1,.1);

      me00 = me_proj[mc_type_code][ibin][ibin3]->GetFunction("pol0")->GetParameter(0);
      me_proj[mc_type_code][ibin][ibin3]->SetAxisRange(-.49,.49);
      me_proj[mc_type_code][ibin][ibin3]->Draw();

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
      hJetTrackME[mc_type_code][ibin][ibin3]->Scale(1./me00);
      yield_inc[mc_type_code][ibin][ibin3]->Divide(hJetTrackME[mc_type_code][ibin][ibin3]);

      hJetTrackSignalBackground[mc_type_code][ibin][ibin3]->Write();
      hJetTrackME[mc_type_code][ibin][ibin3]->Write();
      yield_inc[mc_type_code][ibin][ibin3]->Write();

  
      //-------
      //Leading
      //------


      fout_lead[mc_type_code]->cd();

      me_proj_lead_canvas[mc_type_code]->cd(4*(ibin3+1)-ibin);

      yield_lead[mc_type_code][ibin][ibin3] =  (TH2D*) hJetTrackSignalBackgroundLeading[mc_type_code][ibin][ibin3]->Clone((TString)("Yield_Leading_"+data_mc_type_strs[mc_type_code]+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_"+TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));
       

      me_proj_lead[mc_type_code][ibin][ibin3] = (TH1D*)hJetTrackMELeading[mc_type_code][ibin][ibin3]->ProjectionX(Form("me_proj_lead_temp_%d%d%d",mc_type_code,ibin,ibin3),1,100);

      me_proj_lead[mc_type_code][ibin][ibin3]->Scale (1./100);
      me_proj_lead[mc_type_code][ibin][ibin3] -> Fit("pol0","","",-.1,.1);

      me00 = 	me_proj_lead[mc_type_code][ibin][ibin3]->GetFunction("pol0")->GetParameter(0);
      me_proj_lead[mc_type_code][ibin][ibin3]->SetAxisRange(-.49,.49);
      me_proj_lead[mc_type_code][ibin][ibin3]->Draw();

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

      hJetTrackMELeading[mc_type_code][ibin][ibin3]->Scale(1./me00);
      yield_lead[mc_type_code][ibin][ibin3]->Divide(hJetTrackMELeading[mc_type_code][ibin][ibin3]);

      hJetTrackSignalBackgroundLeading[mc_type_code][ibin][ibin3]->Write();
      hJetTrackMELeading[mc_type_code][ibin][ibin3]->Write();
      yield_lead[mc_type_code][ibin][ibin3]->Write();

  

      //------------
      //Subleading
      //----------
      
      fout_sub[mc_type_code]->cd();

      me_proj_sub_canvas[mc_type_code]->cd(4*(ibin3+1)-ibin);

      yield_sub[mc_type_code][ibin][ibin3] =  (TH2D*) hJetTrackSignalBackgroundSubLeading[mc_type_code][ibin][ibin3]->Clone((TString)("Yield_SubLeading_"+data_mc_type_strs[mc_type_code]+"_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+ "_" +TrkPtBin_strs[ibin3]+"_"+TrkPtBin_strs[ibin3+1]));
       
      me_proj_sub[mc_type_code][ibin][ibin3] = (TH1D*)hJetTrackMESubLeading[mc_type_code][ibin][ibin3]->ProjectionX(Form("me_proj_sub_temp_%d%d%d",mc_type_code,ibin,ibin3),1,100);


      me_proj_sub[mc_type_code][ibin][ibin3]->Scale (1./100);
      me_proj_sub[mc_type_code][ibin][ibin3] -> Fit("pol0","","",-.1,.1);

      me00 = me_proj_sub[mc_type_code][ibin][ibin3]->GetFunction("pol0")->GetParameter(0);
      me_proj_sub[mc_type_code][ibin][ibin3]->SetAxisRange(-.49,.49);
      me_proj_sub[mc_type_code][ibin][ibin3]->Draw();

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

      hJetTrackMESubLeading[mc_type_code][ibin][ibin3]->Scale(1./me00);
      yield_sub[mc_type_code][ibin][ibin3]->Divide(hJetTrackMESubLeading[mc_type_code][ibin][ibin3]);

      hJetTrackSignalBackgroundSubLeading[mc_type_code][ibin][ibin3]->Write();
      hJetTrackMESubLeading[mc_type_code][ibin][ibin3]->Write();
      yield_sub[mc_type_code][ibin][ibin3]->Write();

  

      //--------------------
      //Draw result histos
      //---------------------


      int phil = yield_inc[mc_type_code][ibin][ibin3]->GetYaxis()->FindBin(-.99999);
      int phir = yield_inc[mc_type_code][ibin][ibin3]->GetYaxis()->FindBin(.99999);

      result_proj_canvas[mc_type_code]->cd(4*(ibin3+1)-ibin);
	
      yield_inc_proj[mc_type_code][ibin][ibin3]= (TH1D*)yield_inc[mc_type_code][ibin][ibin3]->ProjectionX((TString)("yield_inc_proj_"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);
	
      yield_inc_proj[mc_type_code][ibin][ibin3]->SetMarkerStyle(10);
      yield_inc_proj[mc_type_code][ibin][ibin3]->SetMarkerSize(1);
      yield_inc_proj[mc_type_code][ibin][ibin3]->SetMarkerColor(kBlack);
      yield_inc_proj[mc_type_code][ibin][ibin3]->SetAxisRange(-2.99,2.99);
      yield_inc_proj[mc_type_code][ibin][ibin3]->GetYaxis()->SetLabelSize(0.06);
      yield_inc_proj[mc_type_code][ibin][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
      yield_inc_proj[mc_type_code][ibin][ibin3]->GetXaxis()->CenterTitle();
      yield_inc_proj[mc_type_code][ibin][ibin3]->GetXaxis()->SetTitleSize(0.06);
      yield_inc_proj[mc_type_code][ibin][ibin3]->Draw();

      if(ibin<3){
	centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->Draw();
	pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	pttex->SetNDC();
	pttex->Draw();
      }
      if(ibin==3){
	centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->Draw();
	pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	pttex->SetNDC();
	pttex->Draw();
      }


      result_proj_lead_canvas[mc_type_code]->cd(4*(ibin3+1)-ibin);
	
      yield_lead_proj[mc_type_code][ibin][ibin3]= (TH1D*)yield_lead[mc_type_code][ibin][ibin3]->ProjectionX((TString)("yield_lead_proj_"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);
	
      yield_lead_proj[mc_type_code][ibin][ibin3]->SetMarkerStyle(10);
      yield_lead_proj[mc_type_code][ibin][ibin3]->SetMarkerSize(1);
      yield_lead_proj[mc_type_code][ibin][ibin3]->SetMarkerColor(kBlack);
      yield_lead_proj[mc_type_code][ibin][ibin3]->SetAxisRange(-2.99,2.99);
      yield_lead_proj[mc_type_code][ibin][ibin3]->GetYaxis()->SetLabelSize(0.06);
      yield_lead_proj[mc_type_code][ibin][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
      yield_lead_proj[mc_type_code][ibin][ibin3]->GetXaxis()->CenterTitle();
      yield_lead_proj[mc_type_code][ibin][ibin3]->GetXaxis()->SetTitleSize(0.06);
      yield_lead_proj[mc_type_code][ibin][ibin3]->Draw();

      if(ibin<3){
	centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->Draw();
	pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	pttex->SetNDC();
	pttex->Draw();
      }
      if(ibin==3){
	centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->Draw();
	pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	pttex->SetNDC();
	pttex->Draw();
      }


      result_proj_sub_canvas[mc_type_code]->cd(4*(ibin3+1)-ibin);
	
      yield_sub_proj[mc_type_code][ibin][ibin3]= (TH1D*)yield_sub[mc_type_code][ibin][ibin3]->ProjectionX((TString)("yield_sub_proj_"+CBin_strs[ibin]+"_"+TrkPtBin_strs[ibin3]),phil,phir);
	
      yield_sub_proj[mc_type_code][ibin][ibin3]->SetMarkerStyle(10);
      yield_sub_proj[mc_type_code][ibin][ibin3]->SetMarkerSize(1);
      yield_sub_proj[mc_type_code][ibin][ibin3]->SetMarkerColor(kBlack);
      yield_sub_proj[mc_type_code][ibin][ibin3]->SetAxisRange(-2.99,2.99);
      yield_sub_proj[mc_type_code][ibin][ibin3]->GetYaxis()->SetLabelSize(0.06);
      yield_sub_proj[mc_type_code][ibin][ibin3]->GetXaxis()->SetTitle("#Delta#eta");
      yield_sub_proj[mc_type_code][ibin][ibin3]->GetXaxis()->CenterTitle();
      yield_sub_proj[mc_type_code][ibin][ibin3]->GetXaxis()->SetTitleSize(0.06);
      yield_sub_proj[mc_type_code][ibin][ibin3]->Draw();

      if(ibin<3){
	centtex = new TLatex(0.05,0.9,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->Draw();
	pttex = new TLatex(0.05,0.85,TrkPtBin_labels[ibin3]);
	pttex->SetNDC();
	pttex->Draw();
      }
      if(ibin==3){
	centtex = new TLatex(0.15,0.9,CBin_labels[ibin]);
	centtex->SetNDC();
	centtex->Draw();
	pttex = new TLatex(0.15,0.85,TrkPtBin_labels[ibin3]);
	pttex->SetNDC();
	pttex->Draw();
      }
      
	
    }//ibin3;
      
  } // ibin ( centrality ) loop
  me_proj_canvas[mc_type_code]->SaveAs((TString)("ME_Inclusive_Projections_"+data_mc_type_strs[mc_type_code]+".png"));
  me_proj_sub_canvas[mc_type_code]->SaveAs((TString)("ME_Leading_Projections_"+data_mc_type_strs[mc_type_code]+".png"));
  me_proj_lead_canvas[mc_type_code]->SaveAs((TString)("ME_SubLeading_Projections_"+data_mc_type_strs[mc_type_code]+".png"));
    
  result_proj_canvas[mc_type_code]->SaveAs((TString)("Result_Inclusive_Projections_"+data_mc_type_strs[mc_type_code]+".png"));
  result_proj_lead_canvas[mc_type_code]->SaveAs((TString)("Result_Leading_Projections_"+data_mc_type_strs[mc_type_code]+".png"));
  result_proj_sub_canvas[mc_type_code]->SaveAs((TString)("Result_SubLeading_Projections_"+data_mc_type_strs[mc_type_code]+".png"));
  
  return 0;
} // main loop
