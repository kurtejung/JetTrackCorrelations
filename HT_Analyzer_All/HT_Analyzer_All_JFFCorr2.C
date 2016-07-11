/// reco PbPb
#include <iostream>
#include "TFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TF1.h>
#include "assert.h"
#include <fstream>
#include "TMath.h"
#include "mixing_tree.h"
#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>
#include <vector>
#include "TCanvas.h"
#include "Jan_24_pp_Iterative/getTrkCorr.h"

using namespace std;

#define nCBins 4
#define nPtBins 1 
#define nTrkPtBins 5

float trkPtCut=1;

int parti = -999;
bool is_data = false;

enum enum_dataset_types {e_Data2011,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370,e_Pythia15,e_Pythia30,e_Pythia50, e_Pythia80, e_Pythia120,e_Pythia170,e_Pythia220,e_Pythia280, e_Pythia370, e_n_dataset_types};
int dataset_type_code = -999;

TString dataset_type_strs[e_n_dataset_types] = {"Data2011","Data_pp","HydJet15","HydJet30","HydJet50","HydJet80", "HydJet120", "HydJet170","HydJet220","HydJet280","HydJet370","Pythia15","Pythia30","Pythia50","Pythia80", "Pythia120", "Pythia170","Pythia220","Pythia280","Pythia370"};
//Hydjet80 = 5
//Pythia80 = 14

TString dataset_type_file_names[e_n_dataset_types] = {"ClusterData2011.txt","ClusterData_pp.txt","Hydjet15.txt","Hydjet30.txt","Hydjet50.txt","Hydjet80.txt", "Hydjet120.txt", "Hydjet170.txt","Hydjet220.txt","Hydjet280.txt","Hydjet370.txt","Pythia15.txt","Pythia30.txt","Pythia50.txt","Pythia80.txt", "Pythia120.txt", "Pythia170.txt","Pythia220.txt","Pythia280.txt","Pythia370.txt"};

int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,15,30,50,80,120,170,220,280,370,999};


enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RightGen, SpilledUnderGen, UnmatchedGen, RightReco, SpilledReco, UnmatchedReco, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,n_data_mc_types};


TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RightGenJet_GenTrack","SpilledUnderJet_GenTrack","UnmatchedGenJet_GenTrack","RightRecoJet_GenTrack","SpilledReco_GenTrack","UnmatchedReco_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};
int data_mc_type_code = -999;


float PtBins[nPtBins+1] = {100, 300};
TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};

float CBins[nCBins+1] = {0, 20, 60, 100, 200};
TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};

float TrkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};
TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt999" };


const int npt=29; 

double ptmin_pbpb[npt]={  0.5,  0.5, 0.5, 0.5, 0.5, 0.55, 0.55, 0.55, 0.55, 0.55, 0.65, 0.65, 0.65, 0.65, 0.65, 0.8, 0.8, 0.8, 0.8, 0.8,  1,  1,  1,   1,   1,  3,  3,   3,   8};
double ptmax_pbpb[npt]={ 0.55, 0.55, 0.55, 0.55, 0.55, 0.65, 0.65, 0.65, 0.65, 0.65, 0.8, 0.8, 0.8, 0.8, 0.8,  1,   1,   1,   1,   1,  3,  3,  3,   3,   3,  8,  8,   8,    300};
  
int cent_min[npt]={  0, 20, 40, 60,100,0, 20, 40, 60,100,0, 20, 40, 60,100,0, 20, 40, 60,100, 0,20,40, 60,100, 0,20, 40,     0};
int cent_max[npt]={ 20, 40, 60,100,200,20, 40, 60,100,200,20, 40, 60,100,200,20, 40, 60,100,200,20,40,60,100,200,20,40,200,200};

TFile *f_eff[npt];
TProfile *p_eff_cent[npt]; 
TProfile2D *p_eff_accept[npt]; 
TProfile *p_eff_pt[npt]; 
TProfile *p_eff_rmin[npt];
TFile *f_fake[npt];
TProfile *p_fake_cent[npt]; 
TProfile2D *p_fake_accept[npt]; 
TProfile *p_fake_pt[npt]; 
TProfile *p_fake_rmin[npt];  

const int npt_pp=4;
double ptmin_pp[]={0.5, 1, 3,  8};
double ptmax_pp[]={  1, 3, 8,300};


TFile *f_eff_pp[npt_pp];
TProfile2D *p_eff_accept_pp[npt_pp];  
TProfile *p_eff_pt_pp[npt_pp]; 
TProfile *p_eff_rmin_pp[npt_pp]; 

TFile *f_fake_pp[npt_pp];
TProfile2D *p_fake_accept_pp[npt_pp]; 
TProfile *p_fake_pt_pp[npt_pp]; 
TProfile *p_fake_rmin_pp[npt_pp]; 

TFile *f_secondary;
TH2D * hsecondary;

TFile *f_multrec;
TH2D * hmultrec;


vector <double> dijet_pt;
vector <double> dijet_eta;
vector <double> dijet_phi;
vector <double> dijet_evi;
vector <double> dijet_vz;
vector <double> dijet_cent;

#include "hist_class_def_HT.h"

//Auxillary functions defined below
void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug=false);

///***************************************************************
//     MAIN LOOP STARTS HERE!  
//*************************************************************

// arg1 = dataset type, arg2 = number of files

void HT_Analyzer_All_JFFCorr2(int datasetTypeCode = 1, int nFiles = 1){
 
  dataset_type_code = datasetTypeCode;    //// pick datasets you want to run over
  
  parti = nFiles;

  if(dataset_type_code == e_Data2011 || dataset_type_code == e_Data_pp){
    is_data = true;
    data_mc_type_code = 0;
  } else{
    is_data = false;
    data_mc_type_code = 2;
  }
    
  bool do_mixing = kFALSE;

  std::cout<<"dataset_type_code is " <<dataset_type_code<<" "<<dataset_type_strs[dataset_type_code]<<endl;
  std::cout << "Running with trkPtCut " << trkPtCut << std::endl;
    
  std::vector<TString> file_names;   file_names.clear();

  ReadFileList( file_names, dataset_type_file_names[dataset_type_code], true);
  
  cout<<"got file"<<endl;

  bool is_pp = kFALSE;
    
  if(dataset_type_code == e_Data_pp||dataset_type_code>10){is_pp = kTRUE;}

  int n_data_mc_types_used = n_data_mc_types;
  if(is_pp)n_data_mc_types_used = n_data_mc_types;  //this is because we shouldn't save sube0 scans for pythia

  hist_class *my_hists[n_data_mc_types];
 
  if(is_data){
    my_hists[data_mc_type_code] = new hist_class(data_mc_type_strs[data_mc_type_code], is_data);
  }else{
  
    for(int mc_type_i = 1; mc_type_i < n_data_mc_types_used; mc_type_i++){
      
      cout<<data_mc_type_strs[mc_type_i]<<endl;

      my_hists[mc_type_i] = new hist_class(data_mc_type_strs[mc_type_i], is_data);
    }

    data_mc_type_code = 2;
  }
   
  cout<<"made hist classes"<<endl;
  
  
  //****************************************
  //        ALL CUTS ARE HERE!
  //****************************************
 
  const double etacut = 1.6;
  const double searchetacut = 2.0;
  const double pTmaxcut = 1000.;
  const double pTmincut = 120.;
  const double leadingjetcut = 120. ;
  const double subleadingjetcut = 50. ;
  const double dphicut = 5.*(TMath::Pi())/6. ; 
  const double trketamaxcut = 2.4;

  const bool doBjets = true;
 
  
  //****************************************


  double cent, eta, pt, phi, rmin, r_reco, jeteta, jetphi, fake, eff, secondary, multrec, trkweight, trkweight_lead, trkweight_sub, vz,deta, dphi, reco_eta, gen_eta, reco_phi, gen_phi, dr, closest_dr, jet_dir_eta, jet_dir_phi;
  bool foundjet, foundjet_gen, founddijet, founddijet_gen, is_inclusive;
  int closest_j4i;
	

  double wvz = 1.;
  double wcen = 1.;
 
  TF1 *fit_cen, *fit_vz;

  int unmatched_counter = 0;

  //////////###### PTHAT SAMPLES ###########///////////////

  if(!is_data){
    TFile *f_vertex_cent = new TFile("VertexCentReweightingFits.root","READ");

    if(!is_pp)    fit_cen = (TF1*)f_vertex_cent->Get((TString)("Fit_Cent_"+dataset_type_strs[dataset_type_code]))->Clone((TString)("Fit_Cent_"+dataset_type_strs[dataset_type_code]));

    fit_vz = (TF1*)f_vertex_cent->Get((TString)("Fit_Vz_"+dataset_type_strs[dataset_type_code]))->Clone((TString)("Fit_Vz_"+dataset_type_strs[dataset_type_code]));
 
  }


  //----------------------------------------------------------------
  //    Get histograms for tracking efficiency calculation
  //-------------------------------------------------------------

  TrkCorr* trkCorr;
  if(is_pp) trkCorr = new TrkCorr("Jan_24_pp_Iterative/");
  else trkCorr = new TrkCorr("Jan18_PbPb/");

  
 
  //----------------------------------------------------------------
  //    Obtain reference PbPb Jet Spectra
  //-------------------------------------------------------------
  int pt_weight_bin;
  double pbpb_pt, pp_pt, pt_weight,  pt_weight_lead,  pt_weight_sub; 
 
  TFile *f_ref_pbpb_spectra = new TFile("./PbPb_JetSpectra_JFFCorr.root","READ");
  TFile *f_ref_pp_spectra = new TFile("./pp_JetSpectra_JFFCorr.root","READ");


  TH1D *pbpb_spectrum_inc[nCBins];
  TH1D *pbpb_spectrum_lead[nCBins];
  TH1D *pbpb_spectrum_sub[nCBins];

  TH1D *pp_spectrum_inc[nCBins];
  TH1D *pp_spectrum_lead[nCBins];
  TH1D *pp_spectrum_sub[nCBins];

 
  if(is_pp&&is_data){

    for(int ibin = 0; ibin<4; ibin++){
    
      pbpb_spectrum_inc[ibin] = (TH1D*)f_ref_pbpb_spectra->Get((TString)("PbPb_all_jets_corrpT_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"));
      pbpb_spectrum_inc[ibin]->SetName((TString)("PbPb_jet_specrum_inc_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));

      pbpb_spectrum_lead[ibin] = (TH1D*)f_ref_pbpb_spectra->Get((TString)("PbPb_only_leadingjets_corrpT_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"));
      pbpb_spectrum_lead[ibin]->SetName((TString)("PbPb_jet_specrum_lead_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));

      pbpb_spectrum_sub[ibin] = (TH1D*)f_ref_pbpb_spectra->Get((TString)("PbPb_only_subleadingjets_corrpT_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"));
      pbpb_spectrum_sub[ibin]->SetName((TString)("PbPb_jet_specrum_sub_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));


      pp_spectrum_inc[ibin] = (TH1D*)f_ref_pp_spectra->Get((TString)("pp_all_jets_corrpT_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"));
      pp_spectrum_inc[ibin]->SetName((TString)("pp_jet_specrum_inc_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));

      pp_spectrum_lead[ibin] = (TH1D*)f_ref_pp_spectra->Get((TString)("pp_only_leadingjets_corrpT_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"));
      pp_spectrum_lead[ibin]->SetName((TString)("pp_jet_specrum_lead_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));

      pp_spectrum_sub[ibin] =(TH1D*) f_ref_pp_spectra->Get((TString)("pp_only_subleadingjets_corrpT_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]+"_Pt100_Pt300"));
      pp_spectrum_sub[ibin]->SetName((TString)("pp_jet_specrum_sub_"+CBin_strs[ibin]+"_"+CBin_strs[ibin+1]));

    }
  }
 
    //--------------------------------

    //-----------------------------------------------------------------------
    //  ** START ** READING ** THE ** PRIMARY ** TREE **
    //-----------------------------------------------------------------------


    cout<<"Am I pp? "<<is_pp<<endl;

  //assert(parti <= (int) file_names.size() );
    cout << "file_names size: "<< file_names.size() << endl;

  for(int fi = 0; fi < (int) file_names.size(); fi++) {
    TFile *my_file = TFile::Open(file_names.at(fi));
    std::cout << "Current file: " << ", file_name: " << file_names.at(fi) << ", number " << fi << " of " << file_names.size() << std::endl;
    if(my_file->IsZombie()) {
      std::cout << "Is zombie" << std::endl;
    }
    
  
    TTree *inp_tree = (TTree*)my_file->Get("mixing_tree");
    mixing_tree *my_primary = new mixing_tree(inp_tree, !is_data);
    std::cout << "Successfully retrieved tree from input file!" << std::endl;
    Long64_t n_evt = my_primary->fChain->GetEntriesFast();


    TString me_file_name;
    if(is_data&&!is_pp){
      me_file_name = "/data/htrauger/PbPb_MinimumBias_7_16/Data2011_MinimumBias_Combined.root";
    }else{
      if( parti< (int) file_names.size()-1){
	me_file_name = file_names.at(parti+1);
      }else{
	me_file_name = file_names.at(0);
      }
    }
    
    cout<<me_file_name<<endl;
    
    TFile  *me_file= new TFile(me_file_name,"READ");
    TTree *inp_tree2 = (TTree*)me_file->Get("mixing_tree");
    mixing_tree *me_tree = new mixing_tree(inp_tree2);
  
    Long64_t nme = me_tree->fChain->GetEntriesFast();

   
    int meptrig = 40;
    
    gRandom->SetSeed(0);
    Long64_t  me = gRandom->Rndm()*nme;

    TH1D * centbins = new TH1D("centbins","centbins. JUST A DUMMY REALLY", 40, 0.0, 200.0);
    TH1D * vzbins = new TH1D("vzbins","vzbins. JUST A DUMMY REALLY", 30, -15., 15.);
    int jet_cent, jet_vzbin;
    vector <Int_t> me_cent;
    vector <Int_t> me_vzbin;
    vector <Int_t> me_hlt;
   

    if(do_mixing){ 
      cout<<"There are "<<nme<<" events in the MB file. First, we run through them all once."<<endl;
 
      for(int mei = 0; mei <nme; mei++){
	me_tree->fChain->GetEntry(mei);
	if(!is_pp&&!is_data){  
	  me_hlt.push_back(me_tree->HLT_HIJet80_v7);
	  me_cent.push_back((centbins->FindBin(me_tree->hiBin)));
	}else if(is_pp){
	  me_hlt.push_back(me_tree->HLT_PAJet80_NoJetID_v1);
	}else{
	  me_cent.push_back((centbins->FindBin(me_tree->hiBin)));
	  me_hlt.push_back(me_tree->HLT_PAJet80_NoJetID_v1);
	}

	me_vzbin.push_back((vzbins->FindBin(me_tree->vz->at(0))));
	if(vzbins->FindBin(me_tree->vz->at(0))>31){cout<<"THIS IS A PROBLEM!!"<<endl;}
      }
   
    }
    ///==========================   Event Loop starts ===================================
    ///==========================   Event Loop starts ===================================
  
    // n_evt = 1000;
    
    int genjet_count = 0;
    int genjet_fill_count = 0;


    cout << "total Events: "<< n_evt << endl;

    for(int evi = 0; evi < n_evt; evi++) {
     

   
      my_primary->fChain->GetEntry(evi);

      if (evi%1000==0) std::cout << " I am running on file " << fi+1 << " of " << ((int) file_names.size()) << ", evi: " << evi << " of " << n_evt << std::endl;
   

      Int_t hiBin = my_primary->hiBin;
   
      vz = my_primary->vz->at(0);

      if(!is_data){data_mc_type_code = 2;}

      my_hists[data_mc_type_code]->NEvents->Fill(hiBin/2.0);
      my_hists[data_mc_type_code]->Centrality->Fill(hiBin);
      my_hists[data_mc_type_code]->Vz->Fill(vz);
 
      wvz=1;
      wcen=1;

    
      if(!is_data){
	
     	wvz = fit_vz->Eval(vz);
	my_hists[data_mc_type_code]->Vz_new->Fill(vz,wvz);
	
	if(!is_pp){
	  wcen = fit_cen->Eval(1.*hiBin);
	  my_hists[data_mc_type_code]->Centrality_new->Fill(hiBin, wcen);
	}		       
      }

        
      int ibin2 = 0;  int ibin3=0;

      if(is_data) {
	int noise_event_selection = my_primary->pHBHENoiseFilter;
	if(noise_event_selection==0) continue;
      
	if(is_pp){
	  int event_selection = my_primary->pPAcollisionEventSelectionPA; // 2013 pp data
	  if(event_selection==0) {continue; }
	  if(my_primary->HLT_PAJet80_NoJetID_v1==0) {continue; }/// 2013 pp data
	}
	
	if(!is_pp){
	  if (my_primary->HLT_HIJet80_v1==0){ continue;}
	  int event_selection = my_primary->pcollisionEventSelection;
	  if(event_selection==0){ continue; }
	}

      } else{
	
	if(is_pp){
	  if(my_primary->HLT_PAJet80_NoJetID_v1==0) {continue; }///for Pelin's pythia
	}

	if(!is_pp){
	  if (my_primary->HLT_HIJet80_v7==0){ continue;}  
	}
	
	double evt_pthat = my_primary->pthat;
	if(evt_pthat > dataset_pthats[dataset_type_code+1]){continue; }


      }


      if(fabs(vz) > 15.) continue;      

      my_hists[data_mc_type_code]->NEvents_after_noise->Fill(hiBin/2.0);
      //      cout<<"here 0"<<endl;

      //---------------------------------------------------------------------------------------
      ///////// -------- FIND DIJETS (will fill hists along with inclusive) ------//////////
      //---------------------------------------------------------------------------------------
    
      double Aj = -99.;
      double lead_pt=0. ;
      double sublead_pt=0. ;
      int second_highest_idx=-1 ;
      int highest_idx=-1 ;

      //search for leading jet
      for(int j4i = 0; j4i < (int) my_primary->corrpt->size() ; j4i++) {
	double jet_pt= my_primary->corrpt->at(j4i);
	if(TMath::Abs(my_primary->jteta->at(j4i))>=searchetacut) continue ;
	if(jet_pt<=leadingjetcut) continue ;
  if(doBjets && my_primary->discr_csvV1->at(j4i) < 0.9) continue;
	if(jet_pt >lead_pt){
	  lead_pt=jet_pt;
	  highest_idx=j4i;
	}
      } //search for leading jet loop
    
      //search for subleading jet
      for(int ijet = 0 ; ijet < (int) my_primary->corrpt->size(); ijet++){
	if(ijet==highest_idx) continue ;
	if(TMath::Abs(my_primary->jteta->at(ijet))>= searchetacut) continue ;
	if(my_primary->corrpt->at(ijet)<=subleadingjetcut) continue ;
	if(my_primary->corrpt->at(ijet) > sublead_pt){
	  sublead_pt=my_primary->corrpt->at(ijet);
	  second_highest_idx=ijet;
	}
      }  //end of subleading jet search

      if(highest_idx< 0 || second_highest_idx< 0 ){   //only apply dijet cuts to dijets
	highest_idx = -1;
	second_highest_idx = -1;
      }
 
      if(highest_idx> -1 && second_highest_idx> -1 ){   //only apply dijet cuts to dijets
	
	//	dphi =TMath::Abs( (my_primary->jtphi->at(highest_idx))-(my_primary->jtphi->at(second_highest_idx)));
	dphi =  my_primary->jtphi->at(highest_idx) - my_primary->jtphi->at(second_highest_idx);
	if(dphi<0){dphi = -dphi;}
	if(dphi>TMath::Pi()) { dphi = 2*TMath::Pi() - dphi; }


	if((my_primary->corrpt->at(highest_idx)<= leadingjetcut )||
	   (my_primary->corrpt->at(highest_idx)>= pTmaxcut ) ||
	   (my_primary->corrpt->at(highest_idx)<= pTmincut ) ||
	   ( my_primary->corrpt->at(second_highest_idx)<= subleadingjetcut) ||
	   (TMath::Abs(my_primary->jteta->at(highest_idx)) >= etacut ) ||
	   (TMath::Abs(my_primary->jteta->at(second_highest_idx)) >= etacut )||
	   (TMath::Abs(dphi)<= dphicut)){
	  highest_idx = -1;  
	  second_highest_idx = -1; 
	}
      }else{
	founddijet = kTRUE;
      }
     
      

      //----------Then we find gen dijets if MC -----------------
      //      (This structure facilitates rightjets)
      //--------------------------------------------------------

        
      double lead_gen_pt=0. ;
      double sublead_gen_pt=0. ;
      int second_highest_idx_gen=-1 ;
      int highest_idx_gen=-1 ;

      if(!is_data){
	//search for leading jet
	for(int j4i = 0; j4i < (int) my_primary->genpt->size() ; j4i++) {
	  double jet_pt= my_primary->genpt->at(j4i);
	  if(TMath::Abs(my_primary->geneta->at(j4i))>=searchetacut) continue ;
	  if(jet_pt<=leadingjetcut) continue ;
	  if(jet_pt >lead_gen_pt){
	    lead_gen_pt=jet_pt;
	    highest_idx_gen=j4i;
	  }
	} //search for leading jet loop
    
	  //search for subleading jet
	for(int ijet = 0 ; ijet < (int) my_primary->genpt->size(); ijet++){
	  if(ijet==highest_idx_gen) continue ;
	  if(TMath::Abs(my_primary->geneta->at(ijet))>= searchetacut) continue ;
	  if(my_primary->genpt->at(ijet)<=subleadingjetcut) continue ;
	  if(my_primary->genpt->at(ijet) > sublead_gen_pt){
	    sublead_gen_pt=my_primary->genpt->at(ijet);
	    second_highest_idx_gen=ijet;
	  }
	}  //end of subleading jet search

	if(highest_idx_gen< 0 || second_highest_idx_gen< 0 ){   //only apply dijet cuts to dijets
	  highest_idx_gen = -1;
	  second_highest_idx_gen = -1;
	}
 
	if(highest_idx_gen> -1 && second_highest_idx_gen> -1 ){   //only apply dijet cuts to dijets
	
	  //	dphi =TMath::Abs( (my_primary->genphi->at(highest_idx_gen))-(my_primary->genphi->at(second_highest_idx_gen)));
	  dphi =  my_primary->genphi->at(highest_idx_gen) - my_primary->genphi->at(second_highest_idx_gen);
	  if(dphi<0){dphi = -dphi;}
	  if(dphi>TMath::Pi()) { dphi = 2*TMath::Pi() - dphi; }


	  if((my_primary->genpt->at(highest_idx_gen)<= leadingjetcut )||
	     (my_primary->genpt->at(highest_idx_gen)>= pTmaxcut ) ||
	     (my_primary->genpt->at(highest_idx_gen)<= pTmincut ) ||
	     ( my_primary->genpt->at(second_highest_idx_gen)<= subleadingjetcut) ||
	     (TMath::Abs(my_primary->geneta->at(highest_idx_gen)) >= etacut ) ||
	     (TMath::Abs(my_primary->geneta->at(second_highest_idx_gen)) >= etacut )||
	     (TMath::Abs(dphi)<= dphicut)){
	    highest_idx_gen = -1;  
	    second_highest_idx_gen = -1; 
	  }
	}else{
	  founddijet_gen = kTRUE;
	}

      } //is_data
      //----------------------------------------------------------------------------
      // Have dijet information.  Time to start filling bins.
      //----------------------------------------------------------------------------
       

      //	cout<<"here 1"<<endl;

      if(!is_data){ data_mc_type_code = 2; } //General event info we put in RecoGen, since we use this for nominal...Setting this here is actually redundant.
    
      for (int ibin=0;ibin<nCBins; ibin ++){
	if (!is_pp&&(my_primary->hiBin<CBins[ibin] || my_primary->hiBin >=CBins[ibin+1])){ continue; }
    	
	if(highest_idx > -1 && second_highest_idx > -1){ 
	  my_hists[data_mc_type_code]->NEvents_dijets->Fill(hiBin/2.0);
	  my_hists[data_mc_type_code]->dPhi_hist[ibin]->Fill(fabs(dphi));
	  Aj = (my_primary->corrpt->at(highest_idx) - my_primary->corrpt->at(second_highest_idx))/(my_primary->corrpt->at(highest_idx) + my_primary->corrpt->at(second_highest_idx));
	  my_hists[data_mc_type_code]->Aj[ibin]->Fill(Aj); 

	}
      
	/*
	  if(highest_idx > -1 && second_highest_idx > -1){ 
	  if(Aj<0.22){continue;}
	  }
	*/   
	for(int j4i = 0; j4i < (int) my_primary->corrpt->size(); j4i++) {

	  if(!is_data){ data_mc_type_code = 2; } //General event info we put in RecoGen, since we use this for nominal...Setting this here is actually redundant.
	


	  foundjet = kFALSE;
	  is_inclusive = kFALSE;
	  if( fabs(my_primary->jteta->at(j4i)) > etacut ) continue;
	  if( my_primary->corrpt->at(j4i) > pTmaxcut ) continue;
	  if(( my_primary->corrpt->at(j4i) > pTmincut )&&(my_primary->trackMax->at(j4i)/my_primary->corrpt->at(j4i) > 0.01)){
	    is_inclusive = kTRUE;  foundjet = kTRUE;
	  } 

	  
	  ibin2 = 0;  ibin3=0;
        
	  for(int pti = 0; pti < nPtBins; pti++) {
	    if (my_primary->corrpt->at(j4i) >=PtBins[pti] && my_primary->corrpt->at(j4i) < PtBins[pti+1])  ibin2 = pti ;
	  }


	  //Determine gen-jet direction once, for use later in sube0 only residual jff-jec scans.

	  jet_dir_eta = my_primary->jteta->at(j4i);
	  jet_dir_phi = my_primary->jtphi->at(j4i);
	  /*
	    closest_dr = 999.;
	    closest_j4i = -1;

	    for(int j4i_gen = 0; j4i_gen < (int) my_primary->genpt->size(); j4i_gen++) {

	    gen_phi = my_primary->genphi->at(j4i_gen);
	    gen_eta = my_primary->geneta->at(j4i_gen);
	  
	    dr = TMath::Sqrt((jet_dir_eta-gen_eta)*(jet_dir_eta-gen_eta)+(jet_dir_phi-gen_phi)*(jet_dir_phi-gen_phi));
	      
	    if(dr<closest_dr){
	    closest_j4i = j4i_gen;
	    closest_dr = dr;
	    }
	    }// j4i_gen;
	
	    if((is_inclusive||j4i==highest_idx||j4i==second_highest_idx)){
	    if(closest_dr<0.3){
	    jet_dir_eta = my_primary->geneta->at(closest_j4i);	
	    jet_dir_phi = my_primary->genphi->at(closest_j4i);
	    }else{
	    cout<<"No gen jet found: ("<<jet_dir_eta<<","<<jet_dir_phi<<") pT = "<<my_primary->corrpt->at(j4i)<<endl;
	    unmatched_counter++;

	    }
	    }

	  */
     
	  if(is_inclusive == kTRUE){
	    my_hists[data_mc_type_code]->all_jets_corrpT[ibin][ibin2]->Fill(my_primary->corrpt->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->all_jets_phi[ibin][ibin2]->Fill(my_primary->jtphi->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->all_jets_eta[ibin][ibin2]->Fill(my_primary->jteta->at(j4i), wvz*wcen); 
	  }

	  if(is_inclusive == kTRUE &&j4i!=highest_idx){
	    my_hists[data_mc_type_code]->only_nonleadingjets_corrpT[ibin][ibin2]->Fill(my_primary->corrpt->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->only_nonleadingjets_phi[ibin][ibin2]->Fill(my_primary->jtphi->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->only_nonleadingjets_eta[ibin][ibin2]->Fill(my_primary->jteta->at(j4i), wvz*wcen); 
	  }

	  if(j4i==highest_idx){
	    my_hists[data_mc_type_code]->only_leadingjets_corrpT[ibin][ibin2]->Fill(my_primary->corrpt->at(j4i), wvz*wcen);
	    my_hists[data_mc_type_code]->only_leadingjets_phi[ibin][ibin2]->Fill(my_primary->jtphi->at(j4i), wvz*wcen);
	    my_hists[data_mc_type_code]->only_leadingjets_eta[ibin][ibin2]->Fill(my_primary->jteta->at(j4i), wvz*wcen);
	    
	    if(!is_data){
	      if(closest_j4i==highest_idx_gen){
		data_mc_type_code=15;
	      }else if(closest_j4i==second_highest_idx_gen){
		data_mc_type_code=17;
	      }else{
		data_mc_type_code = 19;
	      }

	      my_hists[data_mc_type_code]->only_leadingjets_corrpT[ibin][ibin2]->Fill(my_primary->corrpt->at(j4i), wvz*wcen);
	      my_hists[data_mc_type_code]->only_leadingjets_phi[ibin][ibin2]->Fill(my_primary->jtphi->at(j4i), wvz*wcen);
	      my_hists[data_mc_type_code]->only_leadingjets_eta[ibin][ibin2]->Fill(my_primary->jteta->at(j4i), wvz*wcen);

	      data_mc_type_code = 2;
	    }
	  }

	  if(j4i==second_highest_idx){

	    my_hists[data_mc_type_code]->only_subleadingjets_corrpT[ibin][ibin2]->Fill(my_primary->corrpt->at(j4i), wvz*wcen);
	    my_hists[data_mc_type_code]->only_subleadingjets_phi[ibin][ibin2]->Fill(my_primary->jtphi->at(j4i), wvz*wcen);
	    my_hists[data_mc_type_code]->only_subleadingjets_eta[ibin][ibin2]->Fill(my_primary->jteta->at(j4i), wvz*wcen);

	    if(!is_data){
	      if(closest_j4i==second_highest_idx_gen){
		data_mc_type_code=15;
	      }else if(closest_j4i==highest_idx_gen){
		data_mc_type_code = 17;
	      }else{
		data_mc_type_code = 19;
	      }
	    
	      my_hists[data_mc_type_code]->only_subleadingjets_corrpT[ibin][ibin2]->Fill(my_primary->corrpt->at(j4i), wvz*wcen);
	      my_hists[data_mc_type_code]->only_subleadingjets_phi[ibin][ibin2]->Fill(my_primary->jtphi->at(j4i), wvz*wcen);
	      my_hists[data_mc_type_code]->only_subleadingjets_eta[ibin][ibin2]->Fill(my_primary->jteta->at(j4i), wvz*wcen);

	      data_mc_type_code = 2;
	    }

	  }

	  pt_weight = 1.;
	  pt_weight_sub = 1.;
	  pt_weight_lead = 1.;
	  
	  //------------------------------------
	  // Calculate pt_weight for all tracks
	
	  if(is_pp&&is_data){
	    pt_weight_bin = pbpb_spectrum_inc[ibin]->GetXaxis()->FindBin(my_primary->corrpt->at(j4i));
	    pbpb_pt = pbpb_spectrum_inc[ibin]->GetBinContent(pt_weight_bin);  /// PbPb data
	    pp_pt = pp_spectrum_inc[ibin]->GetBinContent(pt_weight_bin);  /// PbPb data
	    pt_weight = 1.;
	    if(pp_pt>0.0001){ pt_weight = pbpb_pt / pp_pt;
	    }

	    pt_weight_bin = pbpb_spectrum_lead[ibin]->GetXaxis()->FindBin(my_primary->corrpt->at(j4i));
	    pbpb_pt = pbpb_spectrum_lead[ibin]->GetBinContent(pt_weight_bin);  /// PbPb data
	    pp_pt = pp_spectrum_lead[ibin]->GetBinContent(pt_weight_bin);  /// PbPb data
	    pt_weight_lead = 1.;
	    if(pp_pt>0.0001){ pt_weight_lead = pbpb_pt / pp_pt;
	    }

	    pt_weight_bin = pbpb_spectrum_sub[ibin]->GetXaxis()->FindBin(my_primary->corrpt->at(j4i));
	    pbpb_pt = pbpb_spectrum_sub[ibin]->GetBinContent(pt_weight_bin);  /// PbPb data
	    pp_pt = pp_spectrum_sub[ibin]->GetBinContent(pt_weight_bin);  /// PbPb data
	    pt_weight_sub = 1.;
	    if(pp_pt>0.0001){ pt_weight_sub = pbpb_pt / pp_pt;
	    }
	  }

	  //-----------------------------------
	  if(!is_pp) {cent = my_primary->hiBin; }
	  if(!is_data){data_mc_type_code = 1; }


	  for(int tracks =0; tracks < (int) my_primary->trkPt->size(); tracks++){
	    if(fabs(my_primary->trkEta->at(tracks))>=trketamaxcut) continue;
	    if (my_primary->highPurity->at(tracks)!=1) continue;
	    if(my_primary->trkPt->at(tracks)<=trkPtCut) continue;

	  
	    for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
	      if (my_primary->trkPt->at(tracks) >=TrkPtBins[trkpti] && my_primary->trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti;
	    } /// trkpti loop
	  

	    //  Prepare for and call efficiency calculation
	 
	    eta= my_primary->trkEta->at(tracks);
	    pt= my_primary->trkPt->at(tracks);
	    phi= my_primary->trkPhi->at(tracks);
     float rmin = 999;
     for(int k = 0; k<my_primary->jtpt->size(); k++)
     {
      if(my_primary->jtpt->at(k)<50) break;
        if(/*TMath::Abs(my_primary->chargedSum->at(k)/my_primary->jtpt->at(k))<0.01 ||*/ my_primary->jteta->at(k)>2) continue;//jet quality cut
        float R = TMath::Power(my_primary->jteta->at(k)-eta,2)+TMath::Power(my_primary->jtphi->at(k)-phi,2);
        if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
      }

      double trkCorrection = trkCorr->getTrkCorr(pt, eta, phi, 0, rmin);

	    if(!is_pp){
	      secondary = 0.;
	      pt_weight = 1.;
	      pt_weight_lead = 1.;
	      pt_weight_sub = 1.;
	      multrec = 0.;
	    }  //just in case 
	  
	    trkweight = pt_weight*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);

	    trkweight_lead = pt_weight_lead*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);
	    
	    trkweight_sub = pt_weight_sub*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);
	    //---------------------------
	    // Now we are ready to fill!
	    //---------------------------
	
	   

	    my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->trkPt->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),wvz*wcen);
	    
	    my_hists[data_mc_type_code]->TrkPt_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkPt->at(tracks),trkweight*wvz*wcen);
	    my_hists[data_mc_type_code]->TrkEta_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),trkweight*wvz*wcen);
	    my_hists[data_mc_type_code]->TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),trkweight*wvz*wcen);

	   

	    if(is_inclusive == kTRUE){
	   
	    
	      deta = my_primary->jteta->at(j4i) - my_primary->trkEta->at(tracks);
	      dphi = my_primary->jtphi->at(j4i) - my_primary->trkPhi->at(tracks);
	 
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);

	    }
	  

	    if(j4i==highest_idx){
	      deta = my_primary->jteta->at(highest_idx) - my_primary->trkEta->at(tracks);
	      dphi = my_primary->jtphi->at(highest_idx) - my_primary->trkPhi->at(tracks);
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_lead*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);


	    }
	  
	    if(j4i==second_highest_idx){
	      deta = my_primary->jteta->at(second_highest_idx) - my_primary->trkEta->at(tracks);
	      dphi = my_primary->jtphi->at(second_highest_idx) - my_primary->trkPhi->at(tracks);
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_sub*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	    
	    if(is_inclusive && j4i!=highest_idx){
	      deta = my_primary->jteta->at(j4i) - my_primary->trkEta->at(tracks);
	      dphi = my_primary->jtphi->at(j4i) - my_primary->trkPhi->at(tracks);
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundNonLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	    
	  } // Track loop


	  if(!is_data){

	    data_mc_type_code =2;
	    //-------------------------------
	    //   These jets, but gen tracks
	    //-------------------------------

	    for(int tracks =0; tracks < (int) my_primary->pt->size(); tracks++){
	      if(fabs(my_primary->eta->at(tracks))>=trketamaxcut) continue;
	      if(my_primary->pt->at(tracks)<=trkPtCut) continue;
	      if(my_primary->chg->at(tracks)==0) continue;


	  
	      for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
		if (my_primary->pt->at(tracks) >=TrkPtBins[trkpti] && my_primary->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
	      } /// trkpti loop
	  
	      my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->pt->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->eta->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->phi->at(tracks),wvz*wcen);
	    
	      if(is_inclusive == kTRUE){
	   
	    
		deta = my_primary->jteta->at(j4i) - my_primary->eta->at(tracks);
		dphi = my_primary->jtphi->at(j4i) - my_primary->phi->at(tracks);
	 
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    
	      }

	      if(j4i==highest_idx){
		deta = my_primary->jteta->at(highest_idx) - my_primary->eta->at(tracks);
		dphi = my_primary->jtphi->at(highest_idx) - my_primary->phi->at(tracks);
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);


	      }
	  
	      if(j4i==second_highest_idx){
		deta = my_primary->jteta->at(second_highest_idx) - my_primary->eta->at(tracks);
		dphi = my_primary->jtphi->at(second_highest_idx) - my_primary->phi->at(tracks);
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	      }
	    
	      if(is_inclusive && j4i!=highest_idx){
		deta = my_primary->jteta->at(j4i) - my_primary->eta->at(tracks);
		dphi = my_primary->jtphi->at(j4i) - my_primary->phi->at(tracks);
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	      }
	    
	    
	      //------------------------------------------------
	      //Repeat, now subdividing into sube==0 and sube>0.   **duplicate both sube0 scan left for back-compatibility and checks**
	      //-----------------------------------------------

	      if(is_pp||my_primary->sube->at(tracks)==0) data_mc_type_code = 11;
	      else data_mc_type_code = 12;

	    
	      my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->pt->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->eta->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->phi->at(tracks),wvz*wcen);
	    
	      if(is_inclusive == kTRUE){
	   
	    
		deta = jet_dir_eta - my_primary->eta->at(tracks);
		dphi = jet_dir_phi - my_primary->phi->at(tracks);
	 
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    
	      }
	      if(j4i==highest_idx){
		deta = jet_dir_eta - my_primary->eta->at(tracks);
		dphi = jet_dir_phi- my_primary->phi->at(tracks);
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	      }
	  
	      if(j4i==second_highest_idx){
		deta = jet_dir_eta - my_primary->eta->at(tracks);
		dphi = jet_dir_phi- my_primary->phi->at(tracks);
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	      }
	    
	      if(is_inclusive && j4i!=highest_idx){
		deta = jet_dir_eta - my_primary->eta->at(tracks);
		dphi = jet_dir_phi - my_primary->phi->at(tracks);
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	      }
	      
	      if(!is_data){
		if(j4i==highest_idx){
		  if(closest_j4i==highest_idx_gen){
		    if(my_primary->sube->at(tracks)==0||is_pp) data_mc_type_code=15;
		    else data_mc_type_code = 16;
		  }else if(closest_j4i==second_highest_idx_gen){
		    if(my_primary->sube->at(tracks)==0||is_pp) data_mc_type_code=17;
		    else data_mc_type_code = 18;
		  }else{
		    if(my_primary->sube->at(tracks)==0||is_pp) data_mc_type_code=19;
		    else data_mc_type_code = 20;
		  }
  
		  deta = jet_dir_eta - my_primary->eta->at(tracks);
		  dphi = jet_dir_phi- my_primary->phi->at(tracks);
		  while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		  while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		  my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
		}
	  

		if(j4i==second_highest_idx){
		  if(closest_j4i==second_highest_idx_gen){
		    if(my_primary->sube->at(tracks)==0||is_pp) data_mc_type_code=15;
		    else data_mc_type_code = 16;
		  }else if(closest_j4i==highest_idx_gen){
		    if(my_primary->sube->at(tracks)==0||is_pp) data_mc_type_code=17;
		    else data_mc_type_code = 18;
		  }else{
		    if(my_primary->sube->at(tracks)==0||is_pp) data_mc_type_code=19;
		    else data_mc_type_code = 20;
		  }
  
		  deta = jet_dir_eta - my_primary->eta->at(tracks);
		  dphi = jet_dir_phi- my_primary->phi->at(tracks);
	
		  while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		  while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		  my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);




		  data_mc_type_code =2;

		}
	  
	      }	      

	
	    } // Gen particle loop
	    
	      //---------------
	      //Reco Rightjets
	      //---------------

	  
	      //------- Begin RightJets --------
	
	    reco_eta = my_primary->jteta->at(j4i);
	    reco_phi = my_primary->jtphi->at(j4i);

	    closest_dr = 999.;
	    closest_j4i = -1;

	    for(int j4i_gen = 0; j4i_gen < (int) my_primary->genpt->size(); j4i_gen++) {

	      gen_phi = my_primary->genphi->at(j4i_gen);
	      gen_eta = my_primary->geneta->at(j4i_gen);
	  
	      dr = TMath::Sqrt((reco_eta-gen_eta)*(reco_eta-gen_eta)+(reco_phi-gen_phi)*(reco_phi-gen_phi));
	      
	      if(dr<closest_dr){
		closest_j4i = j4i_gen;
		closest_dr = dr;
	      }
	    }// j4i_gen;
	
	
	    //------- End RightJets -------

	    if( closest_dr<0.3&&(my_primary->genpt->at(closest_j4i)>120.)){

	      data_mc_type_code = 8;
	      
	    }else if(closest_dr<0.3&&(my_primary->genpt->at(closest_j4i)<=120.)){

	      data_mc_type_code = 9;
	    }else{
	      data_mc_type_code = 10;
	    }

	    if(is_inclusive == kTRUE){
	      my_hists[data_mc_type_code]->all_jets_corrpT[ibin][ibin2]->Fill(my_primary->corrpt->at(j4i), wvz*wcen); 
	      my_hists[data_mc_type_code]->all_jets_phi[ibin][ibin2]->Fill(my_primary->jtphi->at(j4i), wvz*wcen); 
	      my_hists[data_mc_type_code]->all_jets_eta[ibin][ibin2]->Fill(my_primary->jteta->at(j4i), wvz*wcen); 


	      for(int tracks =0; tracks < (int) my_primary->pt->size(); tracks++){
		if(fabs(my_primary->eta->at(tracks))>=trketamaxcut) continue;
		if(my_primary->pt->at(tracks)<=trkPtCut) continue;
		if(my_primary->chg->at(tracks)==0) continue;
		//	if(my_primary->sube->at(tracks)!=0) continue;  //only pythi for these closures

	  
		for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
		  if (my_primary->pt->at(tracks) >=TrkPtBins[trkpti] && my_primary->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
		} /// trkpti loop
	  
		my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->pt->at(tracks),wvz*wcen);
		my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->eta->at(tracks),wvz*wcen);
		my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->phi->at(tracks),wvz*wcen);
	    
		if(is_inclusive == kTRUE){
	   
	    
		  deta = my_primary->jteta->at(j4i) - my_primary->eta->at(tracks);
		  dphi = my_primary->jtphi->at(j4i) - my_primary->phi->at(tracks);
	 
		  while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		  while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		  my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    
		}
	      
		   
	      } // Gen particle loop
	    }
	  } //!is_data 
	  

	    //----------------------------------------------------
	    //      EVENT MIXING STARTS HERE!  (For DATA ONLY)
	    //-----------------------------------------------------
       
	  if(is_data&&do_mixing){
	    jet_cent = 0;
	    if(!is_pp){ jet_cent = centbins->FindBin(my_primary->hiBin);}
	    jet_vzbin = vzbins->FindBin(my_primary->vz->at(0));

      
	    if(is_inclusive||j4i==highest_idx||j4i==second_highest_idx){  //only if we've got a trigger according to some criteria  AND we're not pp

	      //	  cout << "mixing now " << me <<" "<<nme<<endl;

	      int startovercheck = 0;
	      int mevi = 0;
	      while(mevi< meptrig &&startovercheck <2){  //
		me++;
		if(me>=nme){
		  me=0;
		  cout<<"starting over, startovercheck = "<<startovercheck<<" evi= "<<evi<<" jet_cent = "<<jet_cent<<" jet_vzbin = "<<jet_vzbin<<" mevi = "<<mevi<<" "<<me<<" "<<nme<<endl;  
		  assert(startovercheck<20);
		  startovercheck++; 
		}
	  
		if(me_hlt.at(me)==0) { continue;     }
	  
		//  Centrality matching
		if (!is_pp&&(me_cent.at(me)!=jet_cent)){ continue; }

		// Vz matching
		if(me_vzbin.at(me)==0||me_vzbin.at(me)==31){ continue; }
	   
		if(jet_vzbin!= me_vzbin.at(me)){ continue; }
	    
		me_tree->fChain->GetEntry(me);
		mevi++;

	     
		for(int tracks =0; tracks < (int) me_tree->trkPt->size(); tracks++){
		  if(fabs(me_tree->trkEta->at(tracks))>=trketamaxcut) continue;
		  if (me_tree->highPurity->at(tracks)!=1) continue;
		  if(me_tree->trkPt->at(tracks)<=trkPtCut) continue;
	  
	  
		  for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
		    if (me_tree->trkPt->at(tracks) >=TrkPtBins[trkpti] && me_tree->trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
		  } /// trkpti loop

	
		    //  Prepare efficiency inputs and calculate

		  eta= me_tree->trkEta->at(tracks);
		  pt= me_tree->trkPt->at(tracks);
		  phi= me_tree->trkPhi->at(tracks);
		  if(!is_pp){cent = me_cent.at(me);}
	  
		 float rmin = 999;
     for(int k = 0; k<my_primary->jtpt->size(); k++)
     {
      if(my_primary->jtpt->at(k)<50) break;
        if(/*TMath::Abs(my_primary->chargedSum->at(k)/my_primary->jtpt->at(k))<0.01 ||*/ my_primary->jteta->at(k)>2) continue;//jet quality cut
        float R = TMath::Power(my_primary->jteta->at(k)-eta,2)+TMath::Power(my_primary->jtphi->at(k)-phi,2);
        if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
      }

      double trkCorrection = trkCorr->getTrkCorr(pt, eta, phi, 0, rmin);
	      
		  if(!is_pp){
		    secondary = 0.;
		    pt_weight = 1.;
		    pt_weight_lead = 1.;
		    pt_weight_sub = 1.;
		    multrec = 0;
		  }  //just in case 

		  trkweight = pt_weight*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);

		  trkweight_lead = pt_weight_lead*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);
		  
		  trkweight_sub = pt_weight_sub*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);
 
		  //---------------------------
		  // Now we are ready to fill!
		  //---------------------------
		 
		  my_hists[data_mc_type_code]->ME_TrkPt[ibin][ibin2][ibin3]->Fill(me_tree->trkPt->at(tracks),wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkEta[ibin][ibin2][ibin3]->Fill(me_tree->trkEta->at(tracks),wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkPhi[ibin][ibin2][ibin3]->Fill(me_tree->trkPhi->at(tracks),wvz*wcen);
	    
		  my_hists[data_mc_type_code]->ME_TrkPt_weighted[ibin][ibin2][ibin3]->Fill(me_tree->trkPt->at(tracks),trkweight*wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkEta_weighted[ibin][ibin2][ibin3]->Fill(me_tree->trkEta->at(tracks),trkweight*wvz*wcen);
		  my_hists[data_mc_type_code]->ME_TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(me_tree->trkPhi->at(tracks),trkweight*wvz*wcen);

	    
		  if(is_inclusive){
	   	    
		    deta = my_primary->jteta->at(j4i) - me_tree->trkEta->at(tracks);
		    dphi = my_primary->jtphi->at(j4i) - me_tree->trkPhi->at(tracks);
	 
		    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		    while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		    my_hists[data_mc_type_code]->hJetTrackME[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
		    my_hists[data_mc_type_code]->hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    
		  }
	  

		  if(j4i==highest_idx){
		    deta = my_primary->jteta->at(highest_idx) - me_tree->trkEta->at(tracks);
		    dphi = my_primary->jtphi->at(highest_idx) - me_tree->trkPhi->at(tracks);
		    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		    while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		    my_hists[data_mc_type_code]->hJetTrackMELeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_lead*wvz*wcen);
		    my_hists[data_mc_type_code]->hJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
		  }
	  
		  if(j4i==second_highest_idx){
		    deta = my_primary->jteta->at(second_highest_idx) - me_tree->trkEta->at(tracks);
		    dphi = my_primary->jtphi->at(second_highest_idx) - me_tree->trkPhi->at(tracks);
		    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		    while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		    my_hists[data_mc_type_code]->hJetTrackMESubLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_sub*wvz*wcen);
		    my_hists[data_mc_type_code]->hJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
		  }

 
		  if(is_inclusive && j4i!=highest_idx){
		    deta = my_primary->jteta->at(j4i) - me_tree->trkEta->at(tracks);
		    dphi = my_primary->jtphi->at(j4i) - me_tree->trkPhi->at(tracks);
		    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		    while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		    my_hists[data_mc_type_code]->hJetTrackMENonLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
		    my_hists[data_mc_type_code]->hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
		  }

		}  //track mixed event for data

	  	  
	      } //meptrig events per trigger

	
	    }  //end only do mixed event if we've found a trigger we want to save
	
	  } //that version of mixing was just for data
	

	
	}   /// Closes jpti loop.  THIS MEANS THAT WE TAKE ALL JETS IN AN EVENT >120 GeV, not just the hardest jets.
      
	if(foundjet==kTRUE){my_hists[data_mc_type_code]->NEvents_test->Fill(hiBin/2.);}


	if(is_data){continue;}  //we're ready to move on to the next event in data

      
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////////////
  
      	  
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////   We do ALL of this again, but this time with generated jets.     //////////////
	//                (also we now handle mixing for monte carlo on by-event basis)
	///////////////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////////////////////

	data_mc_type_code = 4;

    
	//----------------------------------------------------------------------------
	// Have dijet information.  Time to start filling bins.
	//----------------------------------------------------------------------------

	//Loop over cent bins, but we pick only the right one to fill for PbPb.  We fill all cent bins (properly weighted each time) for pp.

     	
	if(highest_idx_gen > -1 && second_highest_idx_gen > -1){ 
	  my_hists[data_mc_type_code]->NEvents_dijets->Fill(hiBin/2.0);
	  my_hists[data_mc_type_code]->dPhi_hist[ibin]->Fill(fabs(dphi));
	  Aj = (my_primary->genpt->at(highest_idx_gen) - my_primary->genpt->at(second_highest_idx_gen))/(my_primary->genpt->at(highest_idx_gen) + my_primary->genpt->at(second_highest_idx_gen));
	  my_hists[data_mc_type_code]->Aj[ibin]->Fill(Aj); 
	 
	}
    
	/*
	  if(highest_idx_gen > -1 && second_highest_idx_gen > -1){ 
	  if(Aj<0.22){continue;}
	  }
	*/

	for(int j4i = 0; j4i < (int) my_primary->genpt->size(); j4i++) {

	  data_mc_type_code = 4;

	  foundjet_gen = kFALSE;
	  is_inclusive = kFALSE;
	  if( fabs(my_primary->geneta->at(j4i)) > etacut ) continue;
	  if( my_primary->genpt->at(j4i) > pTmaxcut ) continue;
	  if(my_primary->genpt->at(j4i) > pTmincut ){ is_inclusive = kTRUE; 
	    foundjet_gen = kTRUE;  	
	    genjet_count+=1;
	    
	    //	    cout<<"Here is an inclusive jet! It's the "<<genjet_count<<"th one."<<endl;
	  } 


	  ibin2 = 0;  ibin3=0;
        
	  for(int pti = 0; pti < nPtBins; pti++) {
	    if (my_primary->genpt->at(j4i) >=PtBins[pti] && my_primary->genpt->at(j4i) < PtBins[pti+1])  ibin2 = pti ;
	  }
      
	  if(is_inclusive == kTRUE){
	    
	    genjet_fill_count+=1;
	    //	    cout<<"Filling inclusive "<<data_mc_type_code<<" JetCount: "<<genjet_fill_count<<endl;

	    my_hists[data_mc_type_code]->all_jets_corrpT[ibin][ibin2]->Fill(my_primary->genpt->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->all_jets_phi[ibin][ibin2]->Fill(my_primary->genphi->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->all_jets_eta[ibin][ibin2]->Fill(my_primary->geneta->at(j4i), wvz*wcen); 
	  }

	  if(is_inclusive == kTRUE &&j4i!=highest_idx_gen){
	    my_hists[data_mc_type_code]->only_nonleadingjets_corrpT[ibin][ibin2]->Fill(my_primary->genpt->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->only_nonleadingjets_phi[ibin][ibin2]->Fill(my_primary->genphi->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->only_nonleadingjets_eta[ibin][ibin2]->Fill(my_primary->geneta->at(j4i), wvz*wcen); 
	  }

	  if(j4i==highest_idx_gen){
	    my_hists[data_mc_type_code]->only_leadingjets_corrpT[ibin][ibin2]->Fill(my_primary->genpt->at(j4i), wvz*wcen);
	    my_hists[data_mc_type_code]->only_leadingjets_phi[ibin][ibin2]->Fill(my_primary->genphi->at(j4i), wvz*wcen);
	    my_hists[data_mc_type_code]->only_leadingjets_eta[ibin][ibin2]->Fill(my_primary->geneta->at(j4i), wvz*wcen);
	
	  }

	  if(j4i==second_highest_idx_gen){
	    my_hists[data_mc_type_code]->only_subleadingjets_corrpT[ibin][ibin2]->Fill(my_primary->genpt->at(j4i),wvz*wcen);
	    my_hists[data_mc_type_code]->only_subleadingjets_phi[ibin][ibin2]->Fill(my_primary->genphi->at(j4i),wvz*wcen);
	    my_hists[data_mc_type_code]->only_subleadingjets_eta[ibin][ibin2]->Fill(my_primary->geneta->at(j4i),wvz*wcen);
	  }

	  if(!is_data){	  data_mc_type_code = 3; }
	    
	  pt_weight = 1.;

	  if(!is_pp) {cent = my_primary->hiBin; }

	  for(int tracks =0; tracks < (int) my_primary->trkPt->size(); tracks++){
	    if(fabs(my_primary->trkEta->at(tracks))>=trketamaxcut) continue;
	    if (my_primary->highPurity->at(tracks)!=1) continue;
	    if(my_primary->trkPt->at(tracks)<=trkPtCut) continue;

	  
	    for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
	      if (my_primary->trkPt->at(tracks) >=TrkPtBins[trkpti] && my_primary->trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
	    } /// trkpti loop
	  

	      //  Prepare for and call efficiency calculation

	    eta= my_primary->trkEta->at(tracks);
	    pt= my_primary->trkPt->at(tracks);
	    phi= my_primary->trkPhi->at(tracks);
	    float rmin = 999;
     for(int k = 0; k<my_primary->jtpt->size(); k++)
     {
      if(my_primary->jtpt->at(k)<50) break;
        if(/*TMath::Abs(my_primary->chargedSum->at(k)/my_primary->jtpt->at(k))<0.01 ||*/ my_primary->jteta->at(k)>2) continue;//jet quality cut
        float R = TMath::Power(my_primary->jteta->at(k)-eta,2)+TMath::Power(my_primary->jtphi->at(k)-phi,2);
        if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
      }

      double trkCorrection = trkCorr->getTrkCorr(pt, eta, phi, 0, rmin);

	    if(!is_pp){
	      secondary = 0.;
	      pt_weight = 1.;
	      pt_weight_sub = 1.;
	      pt_weight_lead = 1.;
	      multrec = 0.;
	    }  //just in case 
	  
	    trkweight = pt_weight*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);
	    trkweight_sub = pt_weight_sub*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);
	    trkweight_lead = pt_weight_lead*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);
	
 	
	    //---------------------------
	    // Now we are ready to fill!
	    //---------------------------
	
	    my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->trkPt->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),wvz*wcen);
	    
	    my_hists[data_mc_type_code]->TrkPt_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkPt->at(tracks),trkweight*wvz*wcen);
	    my_hists[data_mc_type_code]->TrkEta_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkEta->at(tracks),trkweight*wvz*wcen);
	    my_hists[data_mc_type_code]->TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(my_primary->trkPhi->at(tracks),trkweight*wvz*wcen);



	    if(is_inclusive == kTRUE){
	   
	    
	      deta = my_primary->geneta->at(j4i) - my_primary->trkEta->at(tracks);
	      dphi = my_primary->genphi->at(j4i) - my_primary->trkPhi->at(tracks);
	 
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    
	    }
	  

	    if(j4i==highest_idx_gen){
	      deta = my_primary->geneta->at(highest_idx_gen) - my_primary->trkEta->at(tracks);
	      dphi = my_primary->genphi->at(highest_idx_gen) - my_primary->trkPhi->at(tracks);
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_lead*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);


	    }
	  
	    if(j4i==second_highest_idx_gen){
	      deta = my_primary->geneta->at(second_highest_idx_gen) - my_primary->trkEta->at(tracks);
	      dphi = my_primary->genphi->at(second_highest_idx_gen) - my_primary->trkPhi->at(tracks);
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_sub*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	    
	    if(is_inclusive && j4i!=highest_idx_gen){
	      deta = my_primary->geneta->at(j4i) - my_primary->trkEta->at(tracks);
	      dphi = my_primary->genphi->at(j4i) - my_primary->trkPhi->at(tracks);
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundNonLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	    
	  } // Track loop
		

	
	    //-------------------------------
	    //   These jets, but gen tracks
	    //-------------------------------

	  data_mc_type_code = 4;
	    
	  for(int tracks =0; tracks < (int) my_primary->pt->size(); tracks++){
	    if(fabs(my_primary->eta->at(tracks))>=trketamaxcut) continue;
	    if(my_primary->pt->at(tracks)<=trkPtCut) continue;
	    if(my_primary->chg->at(tracks)==0) continue;
	    //	    if(my_primary->sube->at(tracks)!=0) continue;

	  
	    for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
	      if (my_primary->pt->at(tracks) >=TrkPtBins[trkpti] && my_primary->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
	    } /// trkpti loop
		
	    my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->pt->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->eta->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->phi->at(tracks),wvz*wcen);
	    
	
	    if(is_inclusive == kTRUE){
	   
	    
	      deta = my_primary->geneta->at(j4i) - my_primary->eta->at(tracks);
	      dphi = my_primary->genphi->at(j4i) - my_primary->phi->at(tracks);
	 
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    
	    }
	  

	    if(j4i==highest_idx_gen){
	      deta = my_primary->geneta->at(highest_idx_gen) - my_primary->eta->at(tracks);
	      dphi = my_primary->genphi->at(highest_idx_gen) - my_primary->phi->at(tracks);
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	  
	    if(j4i==second_highest_idx_gen){
	      deta = my_primary->geneta->at(second_highest_idx_gen) - my_primary->eta->at(tracks);
	      dphi = my_primary->genphi->at(second_highest_idx_gen) - my_primary->phi->at(tracks);
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	    
	    if(is_inclusive && j4i!=highest_idx_gen){
	      deta = my_primary->geneta->at(j4i) - my_primary->eta->at(tracks);
	      dphi = my_primary->genphi->at(j4i) - my_primary->phi->at(tracks);
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }


	    //-------------------------------
	    //  Same but split by sube==0
	    //-----------------------------

	    if(!is_pp){

	      if(my_primary->sube->at(tracks)==0) data_mc_type_code = 13;
	      else data_mc_type_code = 14;
	 
	  
	      my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(my_primary->pt->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(my_primary->eta->at(tracks),wvz*wcen);
	      my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(my_primary->phi->at(tracks),wvz*wcen);
	    
	      if(is_inclusive == kTRUE){
	   
	    
		deta = my_primary->geneta->at(j4i) - my_primary->eta->at(tracks);
		dphi = my_primary->genphi->at(j4i) - my_primary->phi->at(tracks);
	 
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    
	      }
	      if(j4i==highest_idx_gen){
		deta = my_primary->geneta->at(highest_idx_gen) - my_primary->eta->at(tracks);
		dphi = my_primary->genphi->at(highest_idx_gen) - my_primary->phi->at(tracks);
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);


	      }
	  
	      if(j4i==second_highest_idx_gen){
		deta = my_primary->geneta->at(second_highest_idx_gen) - my_primary->eta->at(tracks);
		dphi = my_primary->genphi->at(second_highest_idx_gen) - my_primary->phi->at(tracks);
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	      }
	    
	      if(is_inclusive && j4i!=highest_idx_gen){
		deta = my_primary->geneta->at(j4i) - my_primary->eta->at(tracks);
		dphi = my_primary->genphi->at(j4i) - my_primary->phi->at(tracks);
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
		my_hists[data_mc_type_code]->hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	      }
	    }

	    data_mc_type_code =4;
	
	    
	  } // Gen particle loop
	
	
	    //------------------------------------------------------
	    //   DEBUG.  How exactly do our genjets map onto reco?
	    //----------------------------------------------------
	 
	  gen_phi = my_primary->genphi->at(j4i);
	  gen_eta = my_primary->geneta->at(j4i);
	 
	  closest_dr = 999.;
	  closest_j4i = -1;

	  for(int j4i_reco = 0; j4i_reco < (int) my_primary->corrpt->size(); j4i_reco++) {
 
	    reco_eta = my_primary->jteta->at(j4i_reco);
	    reco_phi = my_primary->jtphi->at(j4i_reco);
	   
	    dr = TMath::Sqrt((reco_eta-gen_eta)*(reco_eta-gen_eta)+(reco_phi-gen_phi)*(reco_phi-gen_phi));
	      
	    if(dr<closest_dr){
	      closest_j4i = j4i_reco;
	      closest_dr = dr;
	    }
	  }// j4i to check all the reco jets' distance

	
	  if((closest_dr < 0.3)&&(my_primary->corrpt->at(closest_j4i)>120.)){
	    data_mc_type_code = 5;
	  }else if((closest_dr < 0.3)&&(my_primary->corrpt->at(closest_j4i)<=120.)){
	    data_mc_type_code = 6;
	  }else{
	    data_mc_type_code = 7;
	  }

	  if(is_inclusive){
	    //	    cout<<"Filling decomposed inclusive "<<data_mc_type_code<<" JetCount: "<<genjet_count<<endl;

	    my_hists[data_mc_type_code]->all_jets_corrpT[ibin][ibin2]->Fill(my_primary->genpt->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->all_jets_phi[ibin][ibin2]->Fill(my_primary->genphi->at(j4i), wvz*wcen); 
	    my_hists[data_mc_type_code]->all_jets_eta[ibin][ibin2]->Fill(my_primary->geneta->at(j4i), wvz*wcen); 

	  	    
	    for(int tracks =0; tracks < (int) my_primary->pt->size(); tracks++){
	      if(fabs(my_primary->eta->at(tracks))>=trketamaxcut) continue;
	      if(my_primary->pt->at(tracks)<=trkPtCut) continue;
	      if(my_primary->chg->at(tracks)==0) continue;
	      //      if(my_primary->sube->at(tracks)!=0) continue;
	  
	      for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
		if (my_primary->pt->at(tracks) >=TrkPtBins[trkpti] && my_primary->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
	      } /// trkpti loop
		
	      deta = my_primary->geneta->at(j4i) - my_primary->eta->at(tracks);
	      dphi = my_primary->genphi->at(j4i) - my_primary->phi->at(tracks);
	 
	      while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	      while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	      my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    
	    } //tracks

	  }
	  
	  ///--------------END DEBUG------------------------
      
	}  //end of gen jet loop.
	

	//Temporarily commented out, pending a decision to not mix MC any more (!) 
    
	//------------------------------------------------------------------
	//     MONTE CARLO EVENT MIXING STARTS HERE!  (Monte Carlo case)
	//-----------------------------------------------------------------

	//	cout<<"here 4"<<endl;

	/*

	  if((!foundjet&&!founddijet&&!foundjet_gen&&!founddijet_gen)||!do_mixing){continue; }  //Only mix if you have a jet.
    
	  int startovercheck = 0;
	  int mevi = 0;
	
  
	  jet_cent = 0;
	  if(!is_pp){ jet_cent = centbins->FindBin(my_primary->hiBin);}
	  jet_vzbin = vzbins->FindBin(my_primary->vz->at(0));

	  //	cout<<"starting matching"<<endl;
	
	  while(mevi< meptrig && startovercheck < 2){  // If we've started over once already, going through a 2nd time won't do any good.
	  me++;
	  
	  if(me>=nme){
	  me=0;
	  cout<<"MC Mixing starting over, startovercheck = "<<startovercheck<<" evi= "<<evi<<" jet_cent = "<<jet_cent<<" jet_vzbin = "<<jet_vzbin<<" mevi = "<<mevi<<" me/nme = "<<me<<"/"<<nme<<endl;  
	  //  assert(startovercheck<20);
	  startovercheck++; 
	  }
	  
	  //	  cout<<"hlt"<<endl;
	   
	  if(me_hlt.at(me)==0) {  continue;     }
	  
	  //	  cout<<"got hlt"<<endl;
	  //  Centrality matching
	  if (!is_pp&&(me_cent.at(me)!=jet_cent)){ continue; }

	  //	  cout<<"vz"<<endl;

	  // Vz matching
	  if(me_vzbin.at(me)==0||me_vzbin.at(me)==31){ continue; }
	    
	  if(jet_vzbin!= me_vzbin.at(me)){ continue; }

	  	    
	  me_tree->fChain->GetEntry(me);
	  mevi++;


	  for (int ibin=0;ibin<nCBins; ibin ++){
	  if (!is_pp&&(my_primary->hiBin<CBins[ibin] || my_primary->hiBin >=CBins[ibin+1])){ continue; }
		
	  //---------------------------------------------------------------
	  //         Matched. Now we fill recojet and genjet mixed events
	  //----------------------------------------------------------------

	  data_mc_type_code = 4;

	  for(int j4i = 0; j4i < (int) my_primary->genpt->size(); j4i++){
	      
	  is_inclusive = kFALSE;

	  if( fabs(my_primary->geneta->at(j4i)) > etacut ) continue;
	  if( my_primary->genpt->at(j4i) > pTmaxcut ) continue;
	  if( my_primary->genpt->at(j4i) > pTmincut){ is_inclusive = kTRUE; } 

		 
	  for(int tracks =0; tracks < (int) me_tree->pt->size(); tracks++){
	  if(fabs(me_tree->eta->at(tracks))>=trketamaxcut) continue;
	  if(me_tree->pt->at(tracks)<=trkPtCut) continue;
	  if(me_tree->chg->at(tracks)==0) continue;
	  if(!is_pp && me_tree->sube->at(tracks)==0) continue;
	  
	  for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
	  if (me_tree->pt->at(tracks) >=TrkPtBins[trkpti] && me_tree->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
	  } /// trkpti loop

	    //---------------------------
	    // Now we are ready to fill!  (GenGen ME)
	    //---------------------------


	    my_hists[data_mc_type_code]->ME_TrkPt[ibin][ibin2][ibin3]->Fill(me_tree->pt->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->ME_TrkEta[ibin][ibin2][ibin3]->Fill(me_tree->eta->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->ME_TrkPhi[ibin][ibin2][ibin3]->Fill(me_tree->phi->at(tracks),wvz*wcen);
	    
	    if(is_inclusive){
	   	    
	    deta = my_primary->geneta->at(j4i) - me_tree->eta->at(tracks);
	    dphi = my_primary->genphi->at(j4i) - me_tree->phi->at(tracks);
	 
	    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	    while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	    my_hists[data_mc_type_code]->hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	  

	    if(j4i==highest_idx_gen){
	    deta = my_primary->geneta->at(highest_idx_gen) - me_tree->eta->at(tracks);
	    dphi = my_primary->genphi->at(highest_idx_gen) - me_tree->phi->at(tracks);
	    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	    while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	    my_hists[data_mc_type_code]->hJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	  
	    if(j4i==second_highest_idx_gen){
	    deta = my_primary->geneta->at(second_highest_idx_gen) - me_tree->eta->at(tracks);
	    dphi = my_primary->genphi->at(second_highest_idx_gen) - me_tree->phi->at(tracks);
	    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	    while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	    my_hists[data_mc_type_code]->hJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }

 
	    if(is_inclusive && j4i!=highest_idx_gen){
	    deta = my_primary->geneta->at(j4i) - me_tree->eta->at(tracks);
	    dphi = my_primary->genphi->at(j4i) - me_tree->phi->at(tracks);
	    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	    while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	    my_hists[data_mc_type_code]->hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	      
	      
	    }/// Track loop for mixed event
		 
	    } //gen jet loop
	  
	    //---------------------------
	    // Now Reco Jets
	    //---------------------------

	    data_mc_type_code = 2;

	    for(int j4i = 0; j4i < (int) my_primary->corrpt->size(); j4i++){
	      
	    is_inclusive = kFALSE;

	    if( fabs(my_primary->jteta->at(j4i)) > etacut ) continue;
	    if( my_primary->corrpt->at(j4i) > pTmaxcut ) continue;
	    if( my_primary->corrpt->at(j4i) > pTmincut){ is_inclusive = kTRUE; } 
		 
	    for(int tracks =0; tracks < (int) me_tree->pt->size(); tracks++){
	    if(fabs(me_tree->eta->at(tracks))>=trketamaxcut) continue;
	    if(me_tree->pt->at(tracks)<=trkPtCut) continue;
	    if(me_tree->chg->at(tracks)==0) continue;
	    if(!is_pp && me_tree->sube->at(tracks)==0) continue;
	  
	    for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
	    if (me_tree->pt->at(tracks) >=TrkPtBins[trkpti] && me_tree->pt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
	    } /// trkpti loop


	    my_hists[data_mc_type_code]->ME_TrkPt[ibin][ibin2][ibin3]->Fill(me_tree->pt->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->ME_TrkEta[ibin][ibin2][ibin3]->Fill(me_tree->eta->at(tracks),wvz*wcen);
	    my_hists[data_mc_type_code]->ME_TrkPhi[ibin][ibin2][ibin3]->Fill(me_tree->phi->at(tracks),wvz*wcen);
	    
	    if(is_inclusive){
	   	    
	    deta = my_primary->jteta->at(j4i) - me_tree->eta->at(tracks);
	    dphi = my_primary->jtphi->at(j4i) - me_tree->phi->at(tracks);
	 
	    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	    while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	    my_hists[data_mc_type_code]->hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	  

	    if(j4i==highest_idx){
	    deta = my_primary->jteta->at(highest_idx) - me_tree->eta->at(tracks);
	    dphi = my_primary->jtphi->at(highest_idx) - me_tree->phi->at(tracks);
	    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	    while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	    my_hists[data_mc_type_code]->hJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	  
	    if(j4i==second_highest_idx){
	    deta = my_primary->jteta->at(second_highest_idx) - me_tree->eta->at(tracks);
	    dphi = my_primary->jtphi->at(second_highest_idx) - me_tree->phi->at(tracks);
	    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	    while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	    my_hists[data_mc_type_code]->hJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }

 
	    if(is_inclusive && j4i!=highest_idx){
	    deta = my_primary->jteta->at(j4i) - me_tree->eta->at(tracks);
	    dphi = my_primary->jtphi->at(j4i) - me_tree->phi->at(tracks);
	    while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	    while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
	    
	    my_hists[data_mc_type_code]->hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
	    }
	      
	      
	    }/// Track loop

      
		 
	    } //gen jet loop


	    //	  cout<<"here done mixing"<<endl;

	    } // cent loop for me
	  	  
	    } //meptrig events per trigger

	*/
      
      } //cent is a big loop
          
	//    cout<<"here at end"<<endl;

    } ///we do EVERYTHING one event at a time.
    
  }//FILE LOOP  (sort of a dummy, since we run on one file at a time).  
 
  cout<<"There were a total of "<<unmatched_counter<<" jets we could not match over all selections."<<endl;

  cout<<"Ready to write"<<endl;
 
  if(is_data){
    my_hists[data_mc_type_code]->Write(0);
    
  }else{
  
    for(int mc_type_i =  1; mc_type_i < n_data_mc_types_used; mc_type_i++){
      cout<<data_mc_type_strs[mc_type_i]<<endl;
      cout<<mc_type_i<<endl;
      my_hists[mc_type_i]->Write(mc_type_i);
    }
  }
  std::cout << "I am FINALLY done!!!" << std::endl;
  
} // end main


void ReadFileList(std::vector<TString> &my_file_names, TString file_of_names, bool debug)
{
  ifstream file_stream(file_of_names);
  std::string line;
  my_file_names.clear();
  if( debug ) std::cout << "Open file " << file_of_names << " to extract files to run over" << std::endl;
  if( file_stream.is_open() ) {
    if( debug ) std::cout << "Opened " << file_of_names << " for reading" << std::endl;
    int line_num = 0;
    while( file_stream >> line ) {
      if( debug ) std::cout << line_num << ": " << line << std::endl;
      TString tstring_line(line);
      if( tstring_line.CompareTo("", TString::kExact) != 0 ) my_file_names.push_back(tstring_line);
      line_num++;
    }
  } else {
    std::cout << "Error, could not open " << file_of_names << " for reading" << std::endl;
    assert(0);
  }
}
