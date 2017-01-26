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
#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>
#include <vector>
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TEnv.h"

//#include "TrkCorr_July22_Iterative_pp_eta2p4/getTrkCorr.h"
#include "TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/getTrkCorr.h"

#include "trkCorrTable/xiaoTrkCorr.h"

using namespace std;

#define nCBins 4
#define nPtBins 1 
#define nTrkPtBins 10

float trkPtCut=0.5;

int parti = -999;
bool is_data = true;

enum enum_dataset_types {e_Data2011,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370,e_Pythia15,e_Pythia30,e_Pythia50, e_Pythia80, e_Pythia120,e_Pythia170,e_Pythia220,e_Pythia280, e_Pythia370, e_n_dataset_types};
int dataset_type_code = -999;

TString dataset_type_strs[e_n_dataset_types] = {"Data2015_PbPb","Data2015_pp","HydJet15","HydJet15","HydJet15","HydJet15", "HydJet15", "HydJet170","HydJet220","HydJet280","HydJet370","Pythia15","Pythia30","Pythia50","Pythia80", "Pythia120", "Pythia170","Pythia220","Pythia280","Pythia370"};
//Hydjet80 = 5
//Pythia80 = 14

int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,15,30,50,80,120,170,220,280,370,999};


enum enum_data_mc_types {Data, RecoReco, RecoGen, GenReco, GenGen, RightGen, SpilledUnderGen, UnmatchedGen, RightReco, SpilledReco, UnmatchedReco, RecoGenSube0,RecoGenNoSube0,GenGenSube0,GenGenNoSube0,MatchedRecoGenSube0,MatchedRecoGenNoSube0,SwappedRecoGenSube0,SwappedRecoGenNoSube0, UnMatchedRecoGenSube0,UnMatchedRecoGenNoSube0,n_data_mc_types};


TString data_mc_type_strs[n_data_mc_types] = {"Data","RecoJet_RecoTrack","RecoJet_GenTrack","GenJet_RecoTrack", "GenJet_GenTrack","RightGenJet_GenTrack","SpilledUnderJet_GenTrack","UnmatchedGenJet_GenTrack","RightRecoJet_GenTrack","SpilledReco_GenTrack","UnmatchedReco_GenTrack","RecoJet_GenTrack_Sube0","RecoJet_GenTrack_NoSube0","GenJet_GenTrack_Sube0","GenJet_GenTrack_NoSube0","MatchedRecoJet_GenTrack_Sube0","MatchedRecoJet_GenTrack_NoSube0","SwappedRecoJet_GenTrack_Sube0","SwappedRecoJet_GenTrack_NoSube0","UnmatchedRecoJet_GenTrack_Sube0","UnmatchedRecoJet_GenTrack_NoSube0",};
int data_mc_type_code = -999;


float PtBins[nPtBins+1] = {100, 300};
TString PtBin_strs[nPtBins+1] = {"Pt100", "Pt300"};

float CBins[nCBins+1] = {0, 20, 60, 100, 200};
TString CBin_strs[nCBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};

float TrkPtBins[nTrkPtBins+1] = {0.5, 0.7, 1, 2, 3, 4, 8, 12, 16, 20, 999};
TString TrkPtBin_strs[nTrkPtBins+1] = {"TrkPt0p5", "TrkPt0p7", "TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt999" };


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

//****************************************
//        ALL CUTS ARE HERE!
//****************************************

const double searchetacut = 1.6;
const double pTmaxcut = 1000.;
const double pTmincut = 120.;
const double leadingjetcut = 120. ;
const double subleadingjetcut = 50. ;
const double dphicut = 5.*(TMath::Pi())/6. ; 
const double trketamaxcut = 2.4;

const bool doBjets = false;

const bool doOnlySube0 = false;
const bool doOnlySubeNot0 = false;
const bool doOnlyGluonJet = false;
const bool doOnlyQuarkJet = false;

//****************************************

//Auxillary functions defined below
void ReadFileList(std::vector<std::string> &my_file_names, TString file_of_names, bool debug=false);

bool passTrackCuts(float trkPt, float trkEta, bool highPurity, float trkChi2, int trkNdof, int trkNlayer, int trkNHit, float pfHcal, float pfEcal){
	if(abs(trkEta)>=trketamaxcut) return false;
	if(!highPurity) return false;
	if((float)trkChi2/(float)trkNdof/(float)trkNlayer > 0.15) return false;
	if(trkNHit<11 && trkPt > 0.7) return false;

	float Et = (pfHcal+pfEcal)/TMath::CosH(trkEta);
	if(!(trkPt<20 || Et > 0.5*trkPt)) return false;
	if(trkPt<=trkPtCut || trkPt > 400) return false;

	return true;

}

bool passGenTrackCuts(float trkPt, float trkEta, int chg, int sube){
    if(fabs(trkEta)>=trketamaxcut) return false;
    if(trkPt<=trkPtCut) return false;
    if(chg==0) return false;
    if(doOnlySube0 && sube!=0) return false;
    if(doOnlySubeNot0 && sube==0) return false;

    return true;

}

///***************************************************************
//     MAIN LOOP STARTS HERE!  
//*************************************************************

// arg1 = dataset type, arg2 = number of files

void HT_Analyzer_All_JFFCorr2_PbPb_Reduced(bool doCrab = 0, int jobID=0, int globalCode=2, int startfile=0, int nFiles = 1){

        bool doGenJets = false;
        bool doGenTracks = false;
        bool useOfficialTrkCorr = false;
   
        if(is_data) globalCode = 0;
        else if(doGenJets){
            if(doGenTracks) globalCode = 4;
            else globalCode = 3;
        }
        else{
            if(doGenTracks) globalCode = 2;
            else globalCode = 1;
        }

        TEnv *gEnv = new TEnv("HT_Ana_Env");
        gEnv->SetValue("TFile.AsyncPrefetching", 1);

	parti = startfile+nFiles;

	bool do_mixing = true;

	std::cout << "WARNING !!!  Running with DR-dependent tracking corrections TURNED OFF!!" << endl;

	std::cout << "Running with trkPtCut " << trkPtCut << std::endl;

	std::vector<std::string> file_names;   file_names.clear();

        bool is_pp = false;

        if(doCrab) ReadFileList( file_names, Form("job_input_file_list_%d.txt",jobID), true);
        else if(is_pp){
            if(is_data) ReadFileList( file_names, "pp_5TeV_Jan2017skim.txt",true);
            else ReadFileList( file_names, "pp5TeV_MC_Pythia6_withPFreweight.txt", true);
        }
        else{
            if(!is_data) ReadFileList( file_names, "PbPb_Hydjet_MC_Jan4Skim_newjetCorr.txt", true);
            else ReadFileList( file_names,"PbPb_5TeV_JanReskim_full.txt", true);
        }
	cout<<"got file"<<endl;

	double cent, eta, pt, phi;
	//, rmin, r_reco, jeteta, jetphi, fake, eff, vz
	double secondary, multrec, trkweight=1, trkweight_lead=1, trkweight_sub=1, deta, dphi=0, reco_eta, gen_eta, reco_phi, gen_phi, dr, closest_dr, jet_dir_eta, jet_dir_phi;
	bool foundjet=false, foundjet_gen, founddijet, founddijet_gen, is_inclusive;
	int closest_j4i=-1;


	double wvz = 1.;
	double wcen = 1.;

	TF1 *fit_cen=NULL, *fit_vz=NULL;

	int unmatched_counter = 0;
	
	hist_class *my_hists[globalCode+1];
	my_hists[globalCode] = new hist_class(data_mc_type_strs[globalCode], 1);

	//////////###### PTHAT SAMPLES ###########///////////////

	cout << "Starting reweighting " << endl;
	if(!is_data){
		//TFile *f_vertex_cent = new TFile("VertexCentReweightingFits.root","READ");
		TFile *f_vertex_cent_kurt = new TFile("VzCentReweights.root");

		//if(!is_pp)    fit_cen = (TF1*)f_vertex_cent->Get((TString)("Fit_Cent_"+dataset_type_strs[dataset_type_code]))->Clone((TString)("Fit_Cent_"+dataset_type_strs[dataset_type_code]));
		if(!is_pp) fit_cen = ((TH1D*)f_vertex_cent_kurt->Get("centReweight"))->GetFunction("pol5");
		if(!fit_cen) cout << " cant find centrality!" << endl;

		//fit_vz = (TF1*)f_vertex_cent->Get((TString)("Fit_Vz_"+dataset_type_strs[dataset_type_code]))->Clone((TString)("Fit_Vz_"+dataset_type_strs[dataset_type_code]));
		fit_vz = ((TH1D*)f_vertex_cent_kurt->Get("vzReweight_cent0"))->GetFunction("gaus");
		if(!fit_vz) cout << "cant find vz reweight" << endl;
		
	}
	else{
		fit_cen=NULL;
		fit_vz=NULL;
	}
	cout << "found reweightings" << endl;


	//----------------------------------------------------------------
	//    Get histograms for tracking efficiency calculation
	//-------------------------------------------------------------

	TrkCorr* trkCorr;
	if(is_pp) trkCorr = new TrkCorr("TrkCorr_July22_Iterative_pp_eta2p4/");
	else trkCorr = new TrkCorr("TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/");

        xiaoTrkCorr* xtc;
        if(!is_pp && !useOfficialTrkCorr){
            if(is_data) xtc = new xiaoTrkCorr("trkCorrTable/inputCorr_v13_data.root");
            else xtc = new xiaoTrkCorr("trkCorrTable/inputCorr_v13_mc.root");
        }

	//----------------------------------------------------------------
	//    Obtain reference PbPb Jet Spectra
	//-------------------------------------------------------------
	int pt_weight_bin;
	double pbpb_pt, pp_pt, pt_weight,  pt_weight_lead=1,  pt_weight_sub=1; 


	//--------------------------------

	//-----------------------------------------------------------------------
	//  ** START ** READING ** THE ** PRIMARY ** TREE **
	//-----------------------------------------------------------------------


	cout<<"Am I pp? "<<is_pp<<endl;

	//assert(parti <= (int) file_names.size() );
	cout << "file_names size: "<< file_names.size() << endl;

	//-----------------------------------------------------------------------
	//  START FILLING MIXING OBJECTS!
	//-----------------------------------------------------------------------

        int entryCalls=0;
	unsigned int meptrig = 20;

	gRandom->SetSeed(0);
	const int nCentMixBins=40;
	const int nVzMixBins=30;

	TH1D * centbins = new TH1D("centbins","centbins. JUST A DUMMY REALLY", nCentMixBins, 0.0, 200.0);
	TH1D * vzbins = new TH1D("vzbins","vzbins. JUST A DUMMY REALLY", nVzMixBins, -15., 15.);

	TChain *me_tree = new TChain("mixing_tree");
        //TTree *me_tree = 0;//(TTree*)me1->Get("mixing_tree");
        //TFile *fin;

	int randNum = gRandom->Rndm()*24 + 1;
	TStopwatch *mixingWatch = new TStopwatch();
	std::string me_file_name;
	int nmixingFiles = 1;
        if(do_mixing){
            int imix=0;
            for(int ifile=0; imix<nmixingFiles; ifile++, imix++){
                //if(ifile >=startfile && ifile < parti) continue;
                if(ifile>=(int)file_names.size()) ifile=0;
                //cout << "for mixing, adding file " << file_names.at(ifile) << endl;
                //TFile *me1 = new TFile(file_names.at(ifile).c_str());
                TFile *me1;
                if(!is_data) me1 = new TFile("/mnt/hadoop/store/user/kjung/PbPb_Pythia6Hydjet_MC_JetTrackSkim/PbPb_Pythia6Hydjet_5TeV_JetTrackCorrCrabSkim_Sept29/crab_PbPb_Pythia6Hydjet_MC_JetTrackSkim/160930_151601/0000/HydJet15_1.root");
                //else me1 = new TFile("/mnt/hadoop/store/user/kjung/PbPb_5TeV_skimJetTrack_Jan2017_looseMerge/PbPb_5TeV_JetTrackCorrCrabSkim_newJFFs/crab_PbPb_5TeV_skimJetTrack_Jan2017_looseMerge/170117_021458/0000/Data2015_311.root");
                //me_tree = (TTree*)me1->Get("mixing_tree");
                else if(is_pp) me_tree->Add(Form("/mnt/hadoop/store/user/kjung/pp_5TeV_skimJetTrack_Jan2017_looseMerge_pfCandReweight/pp_5TeV_JetTrackCorrCrabSkim_newJFFs_pfCandReweight/crab_pp_5TeV_skimJetTrack_Jan2017_looseMerge_pfCandReweight/170123_194927/0000/Data_pp_%d.root",ifile+1));
                else me_tree->Add("/mnt/hadoop/store/user/kjung/PbPb_5TeV_MinBiasSkim/Data2015_1.root");
                //me_file_name = file_names.at(ifile);
                //fin = new TFile(me_file_name.c_str());
            }
        }

	float me_vz, me_pthat;
	int me_hiBin;
	vector<float> *me_trkPt=0, *me_trkEta=0, *me_trkPhi=0, *me_jtpt=0, *me_jteta=0, *me_jtphi=0, *me_trkChi2=0, *me_pfHcal=0, *me_pfEcal=0;
	vector<bool> *me_highPurity=0;
	vector<int> *me_trkNlayer=0, *me_trkNHit=0, *me_trkNdof=0;
	vector<float> *me_pt=0, *me_eta=0, *me_phi=0, *me_genpt=0, *me_geneta=0, *me_genphi=0;
	vector<int> *me_sube=0, *me_chg=0;
	Int_t me_HBHEFilter, me_collisionEventSelection, me_phfCoincFilter;
        std::vector<std::vector<std::vector<unsigned int> > > mixing_lists;
        unsigned int jet_cent, jet_vzbin;

        if(do_mixing){
            //me_tree->SetCacheSize(5*1024*1024*1024);
            me_tree->SetBranchStatus("*",0);
            //me_tree->SetBranchStatus("calo_jteta", 1);
            //me_tree->SetBranchStatus("calo_jtphi", 1);
            //me_tree->SetBranchStatus("calo_jtpt", 1);
            if(!doGenTracks){
                me_tree->SetBranchStatus("trkPt", 1);
                me_tree->SetBranchStatus("trkEta", 1);
                me_tree->SetBranchStatus("trkPhi", 1);
                me_tree->SetBranchStatus("highPurity", 1);
                me_tree->SetBranchStatus("trkChi2", 1);
                me_tree->SetBranchStatus("trkNdof", 1);
                me_tree->SetBranchStatus("trkNlayer", 1);
                me_tree->SetBranchStatus("trkNHit", 1);
                me_tree->SetBranchStatus("pfHcal", 1);
                me_tree->SetBranchStatus("pfEcal", 1);
            }
            me_tree->SetBranchStatus("vz", 1);  
            if(!is_pp) me_tree->SetBranchStatus("hiBin", 1);
            if(doGenTracks){
                me_tree->SetBranchStatus("pt", 1);
                me_tree->SetBranchStatus("eta", 1);
                me_tree->SetBranchStatus("phi", 1);
                me_tree->SetBranchStatus("chg", 1);
                me_tree->SetBranchStatus("sube", 1);
            }
            //me_tree->SetBranchStatus("genpt", 1);
            //me_tree->SetBranchStatus("geneta", 1);
            //me_tree->SetBranchStatus("genphi", 1);
            if(!is_data) me_tree->SetBranchStatus("pthat",1);
            if(!is_pp){
                me_tree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
                me_tree->SetBranchStatus("pcollisionEventSelection",1);
                me_tree->SetBranchStatus("phfCoincFilter3",1);
            }
            else{
                me_tree->SetBranchStatus("pPAprimaryVertexFilter",1);
                me_tree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
            }

            if(!doGenTracks){
                me_tree->SetBranchAddress("trkPt", &me_trkPt);
                me_tree->SetBranchAddress("trkEta", &me_trkEta);
                me_tree->SetBranchAddress("trkPhi", &me_trkPhi);
                me_tree->SetBranchAddress("highPurity", &me_highPurity);
                me_tree->SetBranchAddress("trkChi2", &me_trkChi2);
                me_tree->SetBranchAddress("trkNdof", &me_trkNdof);
                me_tree->SetBranchAddress("trkNlayer", &me_trkNlayer);
                me_tree->SetBranchAddress("trkNHit", &me_trkNHit);
                me_tree->SetBranchAddress("pfHcal", &me_pfHcal);
                me_tree->SetBranchAddress("pfEcal", &me_pfEcal);
            }
            me_tree->SetBranchAddress("vz", &me_vz);  
            if(!is_pp) me_tree->SetBranchAddress("hiBin", &me_hiBin);
            //me_tree->SetBranchAddress("calo_jteta", &me_jteta);
            //me_tree->SetBranchAddress("calo_jtphi", &me_jtphi);
            //me_tree->SetBranchAddress("calo_jtpt", &me_jtpt);
            if(!is_pp){
                me_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&me_HBHEFilter);
                me_tree->SetBranchAddress("pcollisionEventSelection",&me_collisionEventSelection);
                me_tree->SetBranchAddress("phfCoincFilter3",&me_phfCoincFilter);
            }
            else{
                me_tree->SetBranchAddress("pPAprimaryVertexFilter",&me_collisionEventSelection);
                me_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&me_HBHEFilter);
            }
            if(!is_data){
                if(doGenTracks){
                    me_tree->SetBranchAddress("pt", &me_trkPt);
                    me_tree->SetBranchAddress("eta", &me_trkEta);
                    me_tree->SetBranchAddress("phi", &me_trkPhi);
                    me_tree->SetBranchAddress("chg", &me_chg);
                    me_tree->SetBranchAddress("sube", &me_sube);
                }
                //me_tree->SetBranchAddress("genpt", &me_genpt);
                //me_tree->SetBranchAddress("geneta", &me_geneta);
                //me_tree->SetBranchAddress("genphi", &me_genphi);

                me_tree->SetBranchAddress("pthat",&me_pthat);
            }
            for(int ivz=0; ivz<nVzMixBins; ivz++){
                vector<vector<unsigned int> > dummyVect;
                for(int icent=0; icent<nCentMixBins; icent++){
                    vector<unsigned int> dummyVect2;
                    dummyVect.push_back(dummyVect2);
                }
                mixing_lists.push_back(dummyVect);
            }

            cout << "sorting events into vz and centrality categories..." << endl;
            int totalEntries = me_tree->GetEntries();
            for(int me_evt=0; me_evt<totalEntries; me_evt++){
                if(me_evt && me_evt%5000==0) cout << "me_evt: "<< me_evt << " of " << me_tree->GetEntries() << endl;
                me_tree->GetEntry(me_evt);
                if(!is_pp && (!me_HBHEFilter || !me_collisionEventSelection || !me_phfCoincFilter)) continue;
                if(is_pp && !me_collisionEventSelection && !me_HBHEFilter) continue;
                if(abs(me_vz)>=15) continue;
                if(!is_data) if(me_pthat<80) continue;

                unsigned int me_cent = 0;
                if(!is_pp) me_cent = centbins->FindBin(me_hiBin)-1;
                unsigned int me_vzbin = vzbins->FindBin(me_vz)-1;
                //cout << "centrality: "<< me_hiBin << " vz: "<< me_vz << " cent bin: "<< me_cent << " vz bin: " << me_vzbin << endl;
                //cout << "cent bin range: "<< centbins->GetBinLowEdge(me_cent+1) << " " << centbins->GetBinLowEdge(me_cent+2) << endl;
                mixing_lists.at(me_vzbin).at(me_cent).push_back(me_evt);

            }
            cout << "sizes..." << endl;
            for(int imix=0; imix<mixing_lists.size(); imix++){
                cout << "vz bin " << imix << " " << mixing_lists.at(imix).at(0).size() << endl;
            }
            cout << "lists loaded!" << endl;
        }

	mixingWatch->Stop();
	cout << "Mixing load time: "<< mixingWatch->RealTime() << endl;

	//-----------------------------------------------------------------------
	//  Now start reading the signal files
	//-----------------------------------------------------------------------

	for(int fi = startfile; fi < (int) file_names.size() && fi<parti; fi++) {

		std::string filename = file_names.at(fi);
		TFile *my_file = TFile::Open(filename.c_str());
		std::cout << "Current file: " << ", file_name: " << filename << ", number " << fi << " of " << file_names.size() << std::endl;
		if(my_file->IsZombie()) {
			std::cout << "file " << filename << " is zombie" << std::endl;
			exit(0);
		}

		TTree *mixing_tree = (TTree*)my_file->Get("mixing_tree");
		std::cout << "Successfully retrieved tree from input file!" << std::endl;
		Long64_t n_evt = mixing_tree->GetEntriesFast();

		vector<float> *trkEta=0, *trkPhi=0, *trkPt=0, *jteta=0, *jtphi=0, *jtpt=0, *rawpt=0, *corrpt=0;
		vector<UChar_t> *trkAlgo=0;
		vector<bool> *highPurity=0;
		vector<float> *trackMax=0, *trkDxy1=0, *trkDxyError1=0, *trkDz1=0, *trkDzError1=0, *trkPtError=0, *pfEcal=0, *pfHcal=0, *trkMVALoose=0, *trkMVATight=0, *trkNdof=0, *trkChi2=0;
		vector<int> *trkNHit=0, *trkNlayer=0;
		vector<int> *sube=0, *chg=0;
		vector<float> *gpt=0, *gphi=0, *geta=0, *pPt=0, *pPhi=0, *pEta=0, *geneta=0, *genphi=0, *genpt=0;
                vector<float> *discr_csvV1=0;
                vector<int> *flavor=0;

		Int_t hiBin = -999;
		Float_t pthat = -999;
		Float_t vz = -999;
		Int_t HBHEFilter, collisionEventSelection, phfCoincFilter;

		/// higenparticles
                mixing_tree->SetBranchStatus("*",0);
		//mixing_tree->Branch("nTrk", &nTrk, "nTrk/I");
                if(!is_pp){
                    mixing_tree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
                    mixing_tree->SetBranchStatus("pcollisionEventSelection",1);
                    mixing_tree->SetBranchStatus("phfCoincFilter3",1);
                    mixing_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHEFilter);
                    mixing_tree->SetBranchAddress("pcollisionEventSelection",&collisionEventSelection);
                    mixing_tree->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
                }
		else{
                    mixing_tree->SetBranchStatus("pPAprimaryVertexFilter",1);
                    mixing_tree->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
                    mixing_tree->SetBranchAddress("pPAprimaryVertexFilter",&collisionEventSelection);
                    mixing_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHEFilter);
		}
                mixing_tree->SetBranchStatus("trkPt",1);
                mixing_tree->SetBranchStatus("trkEta",1);
                mixing_tree->SetBranchStatus("trkPhi",1);
                mixing_tree->SetBranchStatus("trkAlgo",1);
                mixing_tree->SetBranchStatus("highPurity",1);
                mixing_tree->SetBranchStatus("vz", 1);
                if(!doGenTracks){
                    mixing_tree->SetBranchAddress("trkPt", &trkPt);
                    mixing_tree->SetBranchAddress("trkEta", &trkEta);
                    mixing_tree->SetBranchAddress("trkPhi", &trkPhi);
                    mixing_tree->SetBranchAddress("trkAlgo", &trkAlgo);
                    mixing_tree->SetBranchAddress("highPurity", &highPurity);
                }
                mixing_tree->SetBranchAddress("vz", &vz);

                if(!is_pp){
                    mixing_tree->SetBranchStatus("hiBin", 1);
		    mixing_tree->SetBranchAddress("hiBin", &hiBin);
                }

                if(!doGenJets || doOnlyQuarkJet || doOnlyGluonJet){
                    mixing_tree->SetBranchStatus("calo_jteta",1);
                    mixing_tree->SetBranchStatus("calo_jtphi",1);
                    mixing_tree->SetBranchStatus("calo_jtpt",1);
                    mixing_tree->SetBranchStatus("calo_rawpt", 1);
                    mixing_tree->SetBranchStatus("calo_corrpt", 1);
                    mixing_tree->SetBranchStatus("calo_trackMax", 1);
                    if(!is_data) mixing_tree->SetBranchStatus("calo_refparton_flavor",1);

                    mixing_tree->SetBranchAddress("calo_jteta", &jteta);
                    mixing_tree->SetBranchAddress("calo_jtphi", &jtphi);
                    mixing_tree->SetBranchAddress("calo_jtpt", &jtpt);
                    mixing_tree->SetBranchAddress("calo_rawpt", &rawpt);
                    mixing_tree->SetBranchAddress("calo_corrpt", &corrpt);
                    mixing_tree->SetBranchAddress("calo_trackMax", &trackMax);
                    if(!is_data)mixing_tree->SetBranchAddress("calo_refparton_flavor",&flavor);
                }
               
                if(!is_data) mixing_tree->SetBranchStatus("pthat", 1);
                if(!doGenTracks){
                    mixing_tree->SetBranchStatus("trkDxy", 1);
                    mixing_tree->SetBranchStatus("trkDxyError",1);
                    mixing_tree->SetBranchStatus("trkDz", 1);
                    mixing_tree->SetBranchStatus("trkDzError", 1);
                    mixing_tree->SetBranchStatus("trkPtError", 1);
                    mixing_tree->SetBranchStatus("trkNHit",1);
                    mixing_tree->SetBranchStatus("trkNlayer",1);
                    mixing_tree->SetBranchStatus("trkNdof",1);
                    mixing_tree->SetBranchStatus("trkChi2",1);
                    mixing_tree->SetBranchStatus("pfEcal",1);
                    mixing_tree->SetBranchStatus("pfHcal",1);
                    mixing_tree->SetBranchStatus("trkMVALoose",1);
                    mixing_tree->SetBranchStatus("trkMVATight",1);
                }
		if(!is_data) mixing_tree->SetBranchAddress("pthat", &pthat);
                if(!doGenTracks){
                    mixing_tree->SetBranchAddress("trkDxy", &trkDxy1);
                    mixing_tree->SetBranchAddress("trkDxyError", &trkDxyError1);
                    mixing_tree->SetBranchAddress("trkDz", &trkDz1);
                    mixing_tree->SetBranchAddress("trkDzError", &trkDzError1);
                    mixing_tree->SetBranchAddress("trkPtError", &trkPtError);
                    mixing_tree->SetBranchAddress("trkNHit",&trkNHit);
                    mixing_tree->SetBranchAddress("trkNlayer",&trkNlayer);
                    mixing_tree->SetBranchAddress("trkNdof",&trkNdof);
                    mixing_tree->SetBranchAddress("trkChi2",&trkChi2);
                    mixing_tree->SetBranchAddress("pfEcal",&pfEcal);
                    mixing_tree->SetBranchAddress("pfHcal",&pfHcal);
                    mixing_tree->SetBranchAddress("trkMVALoose",&trkMVALoose);
                    mixing_tree->SetBranchAddress("trkMVATight",&trkMVATight);
                }
                mixing_tree->SetBranchStatus("calo_discr_csvV1",1);
                mixing_tree->SetBranchAddress("calo_discr_csvV1",&discr_csvV1);
                if(!is_data){
                    if(doGenTracks){
                        mixing_tree->SetBranchStatus("pt",1);
                        mixing_tree->SetBranchStatus("phi",1);
                        mixing_tree->SetBranchStatus("eta",1);
                        mixing_tree->SetBranchStatus("chg",1);
                        mixing_tree->SetBranchStatus("sube",1);
                        
                        mixing_tree->SetBranchAddress("pt", &trkPt);
                        mixing_tree->SetBranchAddress("phi", &trkPhi);
                        mixing_tree->SetBranchAddress("eta", &trkEta);
                        mixing_tree->SetBranchAddress("chg", &chg);
                        mixing_tree->SetBranchAddress("sube", &sube);
                    }
                    mixing_tree->SetBranchStatus("pPt",1);
                    mixing_tree->SetBranchStatus("pPhi",1);
                    mixing_tree->SetBranchStatus("pEta",1);
                    mixing_tree->SetBranchAddress("pPt", &pPt);
                    mixing_tree->SetBranchAddress("pPhi", &pPhi);
                    mixing_tree->SetBranchAddress("pEta", &pEta);

                    if(doGenJets && !doOnlyQuarkJet && !doOnlyGluonJet){
                        mixing_tree->SetBranchStatus("genphi", 1);
                        mixing_tree->SetBranchStatus("geneta", 1);
                        mixing_tree->SetBranchStatus("genpt", 1);

                        mixing_tree->SetBranchAddress("geneta", &jteta);
                        mixing_tree->SetBranchAddress("genphi", &jtphi);
                        mixing_tree->SetBranchAddress("genpt", &jtpt);
                        mixing_tree->SetBranchAddress("genpt",&corrpt);
                    }
                    else if(doGenJets){
                        //mixing_tree->SetBranchStatus("calo_jteta", 1);
                        //mixing_tree->SetBranchStatus("calo_jtphi", 1);
                        mixing_tree->SetBranchStatus("calo_refpt", 1);

                        cout << "genjet initialization " << endl;
                        //mixing_tree->SetBranchAddress("calo_jteta", &geneta);
                        //mixing_tree->SetBranchAddress("calo_jtphi", &genphi);
                        mixing_tree->SetBranchAddress("calo_refpt", &genpt);
                    
                    }
                }

		///==========================   Event Loop starts ===================================
		///==========================   Event Loop starts ===================================

		n_evt = mixing_tree->GetEntries();
		int genjet_count = 0;
		int genjet_fill_count = 0;

		TStopwatch evtStopwatch;
                TStopwatch fgStopwatch;
                TStopwatch mixStopwatch;
                TStopwatch mixStopwatch_p1;
                TStopwatch mixStopwatch_p2;
                TStopwatch mixStopwatch_p3;

		cout << "total Events: "<< n_evt << endl;

		evtStopwatch.Start(1);
		
		int data_mc_type_code=globalCode;

		for(int evi = 0; evi < n_evt; evi++) {

                        fgStopwatch.Start(0);
			mixing_tree->GetEntry(evi);

			if (evi%100==0) std::cout << " I am running on file " << fi+1 << " of " << ((int) file_names.size()) << ", evi: " << evi << " of " << n_evt << std::endl;

                        if(is_pp) hiBin=0;
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

			double pthatweight=1.;
			double pthatbins[10] = {15,30,50,80,120,170,220,280,370,9999};
			double xsecs[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
			double pthatEntries[9] = {0, 444180, 485654, 484515, 479933, 448829, 258987, 234835, 50644};
			double ppPthatEntries[9] = {0,0,272902,377559,467823,447683,259111,234347,50942};
                        if(!is_data){
				int ibin=0;
				while(pthat>pthatbins[ibin+1]) ibin++;
				if(!is_pp) pthatweight = (xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];
			        else pthatweight = (xsecs[ibin]-xsecs[ibin+1])/ppPthatEntries[ibin];
                        }
			if(!is_data) wvz*=pthatweight;

                        //wvz = 1.; //temporary fix

			int ibin2 = 0;  int ibin3=0;

                        if(!is_pp && (!HBHEFilter || !collisionEventSelection || !phfCoincFilter)) continue;
			if(is_pp && !collisionEventSelection && !HBHEFilter) continue;
			if(fabs(vz) > 15.) continue;      
			if(!is_data) if(pthat<80) continue;

			jet_cent = centbins->FindBin(hiBin)-1;
			if(is_pp) jet_cent=0;
                        jet_vzbin = vzbins->FindBin(vz)-1;
                        if(do_mixing){
				if(mixing_lists.at(jet_vzbin).at(jet_cent).size()<2) continue;
			}

			my_hists[data_mc_type_code]->NEvents_after_noise->Fill(hiBin/2.0);

			//---------------------------------------------------------------------------------------
			///////// -------- FIND DIJETS (will fill hists along with inclusive) ------//////////
			//---------------------------------------------------------------------------------------

			double Aj = -99.;
			double lead_pt=0. ;
			double sublead_pt=0. ;
			int second_highest_idx=-1 ;
			int highest_idx=-1 ;

			//search for leading jet
			for(int j4i = 0; j4i < (int) corrpt->size() ; j4i++) {
				double jet_pt= corrpt->at(j4i);
				if(TMath::Abs(jteta->at(j4i))>=searchetacut) continue ;
				if(jet_pt<=leadingjetcut) continue ;
				if(jet_pt >lead_pt){
					lead_pt=jet_pt;
					highest_idx=j4i;
				}
			} //search for leading jet loop
			if(highest_idx<0) continue;

			//----------------------------------------------------------------------------
			// Have dijet information.  Time to start filling bins.
			//----------------------------------------------------------------------------

			for (int ibin=0;ibin<nCBins; ibin ++){
				if(!is_pp) if((hiBin<CBins[ibin] || hiBin >=CBins[ibin+1])){ continue; }

				if(highest_idx > -1 && second_highest_idx > -1){ 
					my_hists[data_mc_type_code]->NEvents_dijets->Fill(hiBin/2.0);
					my_hists[data_mc_type_code]->dPhi_hist[ibin]->Fill(fabs(dphi));
					Aj = (corrpt->at(highest_idx) - corrpt->at(second_highest_idx))/(corrpt->at(highest_idx) + corrpt->at(second_highest_idx));
					my_hists[data_mc_type_code]->Aj[ibin]->Fill(Aj); 

				}

				/*
				   if(highest_idx > -1 && second_highest_idx > -1){ 
				   if(Aj<0.22){continue;}
				   }
				 */
				for(int j4i = 0; j4i < (int) corrpt->size(); j4i++) {

                                        //cout << " on jet " << j4i << " of " << corrpt->size() << endl;
					bool foundBjet = false;
					if(!doGenJets){
                                            if(doBjets && discr_csvV1->at(j4i) >= 0.9) foundBjet=true;
                                        }
					foundjet = kFALSE;
					is_inclusive = kFALSE;
                                        if(!is_data){
                                            if(doOnlyQuarkJet && (fabs(flavor->at(j4i))<1 || fabs(flavor->at(j4i))>6)) continue;
                                            if(doOnlyGluonJet && fabs(flavor->at(j4i))!=21) continue;
                                        }
                                        double matchedPt=0, matchedEta=0, matchedPhi=0;
                                        if(doGenJets && (doOnlyQuarkJet || doOnlyGluonJet)){
                                            matchedPt = genpt->at(j4i);
                                            matchedEta= jteta->at(j4i);
                                            matchedPhi= jtphi->at(j4i);
                                            if(matchedPt < 0) continue;
                                        }
                                        if( fabs(jteta->at(j4i)) > searchetacut ) continue;
                                        if(!doGenJets || doOnlyQuarkJet || doOnlyGluonJet){
                                            if(( corrpt->at(j4i) > pTmincut )&&(trackMax->at(j4i)/corrpt->at(j4i) > 0.01)){
                                                is_inclusive = kTRUE;  foundjet = kTRUE;
                                            } 
                                        }
                                        if(doGenJets && !doOnlyQuarkJet && !doOnlyGluonJet && corrpt->at(j4i) > pTmincut ){ is_inclusive=kTRUE; foundjet = kTRUE; }
                                        if(!foundjet) continue;

					//Determine gen-jet direction once, for use later in sube0 only residual jff-jec scans.

                                        double jet_pt, jet_eta, jet_phi;
                                        if(doGenJets && (doOnlyQuarkJet || doOnlyGluonJet)){
                                            jet_pt = matchedPt;
                                            jet_eta = matchedEta;
                                            jet_phi = matchedPhi;
                                        }
                                        else{
                                            jet_pt = corrpt->at(j4i);
                                            jet_eta = jteta->at(j4i);
                                            jet_phi = jtphi->at(j4i);
                                        }

                                        ibin2 = 0;  ibin3=0;

                                        for(int pti = 0; pti < nPtBins; pti++) {
                                            if (jet_pt >=PtBins[pti] && jet_pt < PtBins[pti+1])  ibin2 = pti ;
                                        }

					if(is_inclusive == kTRUE){
						my_hists[data_mc_type_code]->all_jets_corrpT[ibin][ibin2]->Fill(jet_pt, wvz*wcen); 
						my_hists[data_mc_type_code]->all_jets_phi[ibin][ibin2]->Fill(jet_phi, wvz*wcen); 
						my_hists[data_mc_type_code]->all_jets_eta[ibin][ibin2]->Fill(jet_eta, wvz*wcen);
						
						if(foundBjet){
							my_hists[data_mc_type_code]->all_bjets_corrpT[ibin][ibin2]->Fill(jet_pt, wvz*wcen);
							my_hists[data_mc_type_code]->all_bjets_phi[ibin][ibin2]->Fill(jet_phi, wvz*wcen);
							my_hists[data_mc_type_code]->all_bjets_eta[ibin][ibin2]->Fill(jet_eta, wvz*wcen);
						}     
					}

					if(is_inclusive == kTRUE &&j4i!=highest_idx){
						my_hists[data_mc_type_code]->only_nonleadingjets_corrpT[ibin][ibin2]->Fill(jet_pt, wvz*wcen); 
						my_hists[data_mc_type_code]->only_nonleadingjets_phi[ibin][ibin2]->Fill(jet_phi, wvz*wcen); 
						my_hists[data_mc_type_code]->only_nonleadingjets_eta[ibin][ibin2]->Fill(jet_eta, wvz*wcen);
						if(foundBjet){
							my_hists[data_mc_type_code]->only_nonleadingbjets_corrpT[ibin][ibin2]->Fill(jet_pt, wvz*wcen); 
							my_hists[data_mc_type_code]->only_nonleadingbjets_phi[ibin][ibin2]->Fill(jet_phi, wvz*wcen); 
							my_hists[data_mc_type_code]->only_nonleadingbjets_eta[ibin][ibin2]->Fill(jet_eta, wvz*wcen);
						}
					}

					if(j4i==highest_idx){
						my_hists[data_mc_type_code]->only_leadingjets_corrpT[ibin][ibin2]->Fill(jet_pt, wvz*wcen);
						my_hists[data_mc_type_code]->only_leadingjets_phi[ibin][ibin2]->Fill(jet_phi, wvz*wcen);
						my_hists[data_mc_type_code]->only_leadingjets_eta[ibin][ibin2]->Fill(jet_eta, wvz*wcen);
						if(foundBjet){
							my_hists[data_mc_type_code]->only_leadingbjets_corrpT[ibin][ibin2]->Fill(jet_pt, wvz*wcen);
							my_hists[data_mc_type_code]->only_leadingbjets_phi[ibin][ibin2]->Fill(jet_phi, wvz*wcen);
							my_hists[data_mc_type_code]->only_leadingbjets_eta[ibin][ibin2]->Fill(jet_eta, wvz*wcen);
						}

					}

					if(j4i==second_highest_idx){

						my_hists[data_mc_type_code]->only_subleadingjets_corrpT[ibin][ibin2]->Fill(jet_pt, wvz*wcen);
						my_hists[data_mc_type_code]->only_subleadingjets_phi[ibin][ibin2]->Fill(jet_phi, wvz*wcen);
						my_hists[data_mc_type_code]->only_subleadingjets_eta[ibin][ibin2]->Fill(jet_eta, wvz*wcen);
						if(foundBjet){
							my_hists[data_mc_type_code]->only_subleadingbjets_corrpT[ibin][ibin2]->Fill(jet_pt, wvz*wcen);
							my_hists[data_mc_type_code]->only_subleadingbjets_phi[ibin][ibin2]->Fill(jet_phi, wvz*wcen);
							my_hists[data_mc_type_code]->only_subleadingbjets_eta[ibin][ibin2]->Fill(jet_eta, wvz*wcen);
						}

					}

					pt_weight = 1.;
					pt_weight_sub = 1.;
					pt_weight_lead = 1.;

					//-----------------------------------
					if(!is_pp) {cent = hiBin; }

					for(int tracks =0; tracks < (int) trkPt->size(); tracks++){
                                            if(!doGenTracks){
                                                if(!passTrackCuts(trkPt->at(tracks), trkEta->at(tracks), highPurity->at(tracks), trkChi2->at(tracks), trkNdof->at(tracks), trkNlayer->at(tracks), trkNHit->at(tracks), pfHcal->at(tracks), pfEcal->at(tracks))) continue;
                                            }
                                            else{
                                                if(!passGenTrackCuts(trkPt->at(tracks), trkEta->at(tracks), chg->at(tracks), sube->at(tracks))) continue;
                                            }
						for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
							if (trkPt->at(tracks) >=TrkPtBins[trkpti] && trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti;
						} /// trkpti loop


						//  Prepare for and call efficiency calculation

						eta= trkEta->at(tracks);
						pt= trkPt->at(tracks);
						phi= trkPhi->at(tracks);
						float rmin = 999;
						/*for(unsigned int k = 0; k<jtpt->size(); k++)
						  {
						  if(jtpt->at(k)<50) break;
						  if(TMath::Abs(chargedSum->at(k)/jtpt->at(k))<0.01 || jteta->at(k)>2) continue;//jet quality cut
						  float R = TMath::Power(jteta->at(k)-eta,2)+TMath::Power(jtphi->at(k)-phi,2);
						  if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
						  }*/

						double trkCorrection;
                                                if(useOfficialTrkCorr) trkCorrection = trkCorr->getTrkCorr(pt, eta, phi, hiBin, rmin);
                                                else if(!is_pp) trkCorrection = xtc->getTrkCorr(pt, eta, phi, hiBin);
                                                else trkCorrection = trkCorr->getTrkCorr(pt, eta, phi, 1, rmin);

						if(!is_pp){
							secondary = 0.;
							pt_weight = 1.;
							pt_weight_lead = 1.;
							pt_weight_sub = 1.;
							multrec = 0.;
						}  //just in case

                                                if(doGenTracks) trkCorrection = 1.;

						trkweight = pt_weight*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);

						trkweight_lead = pt_weight_lead*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);

						trkweight_sub = pt_weight_sub*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);
						//---------------------------
						// Now we are ready to fill!
						//---------------------------



						my_hists[data_mc_type_code]->TrkPt[ibin][ibin2][ibin3]->Fill(trkPt->at(tracks),wvz*wcen);
						my_hists[data_mc_type_code]->TrkEta[ibin][ibin2][ibin3]->Fill(trkEta->at(tracks),wvz*wcen);
						my_hists[data_mc_type_code]->TrkPhi[ibin][ibin2][ibin3]->Fill(trkPhi->at(tracks),wvz*wcen);

						my_hists[data_mc_type_code]->TrkPt_weighted[ibin][ibin2][ibin3]->Fill(trkPt->at(tracks),trkweight*wvz*wcen);
						my_hists[data_mc_type_code]->TrkEta_weighted[ibin][ibin2][ibin3]->Fill(trkEta->at(tracks),trkweight*wvz*wcen);
						my_hists[data_mc_type_code]->TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(trkPhi->at(tracks),trkweight*wvz*wcen);



						if(is_inclusive == kTRUE){


							deta = jet_eta - trkEta->at(tracks);
							dphi = jet_phi - trkPhi->at(tracks);

							while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
							while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}

							my_hists[data_mc_type_code]->hJetTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
							my_hists[data_mc_type_code]->hJetTrackSignalBackground_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, pt*trkweight*wvz*wcen);
                                                        my_hists[data_mc_type_code]->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
							if(foundBjet){
								my_hists[data_mc_type_code]->hbJetTrackSignalBackground[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
								my_hists[data_mc_type_code]->hbJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
							}

						}


						if(j4i==highest_idx){
							deta = jteta->at(highest_idx) - trkEta->at(tracks);
							dphi = jtphi->at(highest_idx) - trkPhi->at(tracks);
							while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
							while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}

							my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_lead*wvz*wcen);
							my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, pt*trkweight_lead*wvz*wcen);
                                                        my_hists[data_mc_type_code]->hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
							if(foundBjet){
								my_hists[data_mc_type_code]->hbJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_lead*wvz*wcen);
								my_hists[data_mc_type_code]->hbJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
							}


						}

						if(j4i==second_highest_idx){
							deta = jteta->at(second_highest_idx) - trkEta->at(tracks);
							dphi = jtphi->at(second_highest_idx) - trkPhi->at(tracks);
							while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
							while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}

							my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_sub*wvz*wcen);
							my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, pt*trkweight_sub*wvz*wcen);
                                                        my_hists[data_mc_type_code]->hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
							if(foundBjet){
								my_hists[data_mc_type_code]->hbJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_sub*wvz*wcen);
								my_hists[data_mc_type_code]->hbJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
							}
						}

					} // Track loop
                                        fgStopwatch.Stop();
					//----------------------------------------------------
					//      EVENT MIXING STARTS HERE!  (For DATA ONLY)
					//-----------------------------------------------------

					if(do_mixing){

                                                mixStopwatch.Start(0);
						jet_cent = 0;
						if(!is_pp){ jet_cent = centbins->FindBin(hiBin)-1;}
						jet_vzbin = vzbins->FindBin(vz)-1;

                                                //cout << "mixing with cent: "<< jet_cent << " vz: " << jet_vzbin << endl;
                                                //cout << "size of mix evt: "<< mixing_lists.at(jet_vzbin).at(jet_cent).size() << endl;
							//randomize starting point of array to mix different events!
						unsigned int pct = gRandom->Rndm()*mixing_lists.at(jet_vzbin).at(jet_cent).size();

                                                //cout << "starting at " << pct << endl;
						bool circleCheck=false;
						
						unsigned int mevi = pct;
						unsigned int dummy=0;
						
                                                while(dummy < meptrig){
                                                        //cout << "mixing event " << dummy << endl;
                                                        mixStopwatch_p1.Start(0);
							if((int)mixing_lists[jet_vzbin][jet_cent].size()<2){
								//cout << "Warning! vzbin: "<< jet_vzbin-1 << ", cent bin: "<< jet_cent-1 << " is devoid of mixing events!" << endl;
								break;
							}
							if(mevi>mixing_lists.at(jet_vzbin).at(jet_cent).size()){ mevi=1; circleCheck=true; }
							
							if(mevi==0) mevi++;
							
							if(mixing_lists.at(jet_vzbin).at(jet_cent).at(mevi-1) == evi){ 
								mevi++; continue; 
							}
							mixStopwatch_p1.Stop();
                                                        mixStopwatch_p2.Start(0);
							entryCalls++;
                                                        //cout << "attempting to get entry " << mevi-1 << endl;
                                                        //cout << "total size: "<< mixing_lists.at(jet_vzbin).at(jet_cent).size() << endl;
                                                        me_tree->GetEntry(mixing_lists.at(jet_vzbin).at(jet_cent).at(mevi-1));
							//cout << " got it" << endl;
                                                        dummy++; mevi++;
							mixStopwatch_p2.Stop();

                                                        mixStopwatch_p3.Start(0);
							unsigned int me_cent=0;
							if(!is_pp){ me_cent = centbins->FindBin(me_hiBin)-1;}
							unsigned int me_vzbin = vzbins->FindBin(me_vz)-1;

								//cout << "vzbin : "<< jet_vzbin << " me vzbin: "<< me_vzbin << endl;
								//make sure i'm doing this thing correctly...
							assert(me_vzbin==jet_vzbin);
							assert(me_cent==jet_cent);
                                                        mixStopwatch_p3.Stop();
                                                        //cout << "getting tracks" << endl;
							for(int tracks =0; tracks < (int) me_trkPt->size(); tracks++){
                                                            if(!doGenTracks){    
                                                                if(!passTrackCuts(me_trkPt->at(tracks), me_trkEta->at(tracks), me_highPurity->at(tracks), me_trkChi2->at(tracks), me_trkNdof->at(tracks), me_trkNlayer->at(tracks), me_trkNHit->at(tracks), me_pfHcal->at(tracks), me_pfEcal->at(tracks))) continue;
                                                            }
                                                            else{
                                                                if(!passGenTrackCuts(me_trkPt->at(tracks), me_trkEta->at(tracks), me_chg->at(tracks), me_sube->at(tracks))) continue;
                                                            }

								for(int trkpti = 0; trkpti < nTrkPtBins; trkpti++) {
									if (me_trkPt->at(tracks) >=TrkPtBins[trkpti] && me_trkPt->at(tracks) < TrkPtBins[trkpti+1])  ibin3 = trkpti ;
								} /// trkpti loop
								
								
									//  Prepare efficiency inputs and calculate
								
								eta= me_trkEta->at(tracks);
								pt= me_trkPt->at(tracks);
								phi= me_trkPhi->at(tracks);
								if(!is_pp){cent = me_hiBin;}
								
								float rmin = 999;
									/*for(unsigned int k = 0; k<me_jtpt->size(); k++)
									  {
									  if(me_jtpt->at(k)<50) break;
									  if(TMath::Abs(chargedSum->at(k)/jtpt->at(k))<0.01 || me_jteta->at(k)>2) continue;//jet quality cut

									  float R = TMath::Power(me_jteta->at(k)-eta,2)+TMath::Power(me_jtphi->at(k)-phi,2);
									  if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
									  }*/

								double trkCorrection;
                                                                if(useOfficialTrkCorr) trkCorrection = trkCorr->getTrkCorr(pt, eta, phi, me_hiBin, rmin);
                                                                else if(!is_pp) trkCorrection = xtc->getTrkCorr(pt,eta,phi,me_hiBin);
                                                                else trkCorrection = trkCorr->getTrkCorr(pt, eta, phi, 1, rmin);

								if(!is_pp){
									secondary = 0.;
									pt_weight = 1.;
									pt_weight_lead = 1.;
									pt_weight_sub = 1.;
									multrec = 0;
								}  //just in case

                                                                if(doGenTracks) trkCorrection = 1.;
								
								
                                                                trkweight = pt_weight*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);
                                                
								trkweight_lead = pt_weight_lead*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);

								trkweight_sub = pt_weight_sub*trkCorrection; //(1-fake)*(1-secondary)/eff/(1+multrec);

								
                                                                        //---------------------------
									// Now we are ready to fill!
									//---------------------------
								
								my_hists[data_mc_type_code]->ME_TrkPt[ibin][ibin2][ibin3]->Fill(me_trkPt->at(tracks),wvz*wcen);
								my_hists[data_mc_type_code]->ME_TrkEta[ibin][ibin2][ibin3]->Fill(me_trkEta->at(tracks),wvz*wcen);
								my_hists[data_mc_type_code]->ME_TrkPhi[ibin][ibin2][ibin3]->Fill(me_trkPhi->at(tracks),wvz*wcen);
								
								my_hists[data_mc_type_code]->ME_TrkPt_weighted[ibin][ibin2][ibin3]->Fill(me_trkPt->at(tracks),trkweight*wvz*wcen);
								my_hists[data_mc_type_code]->ME_TrkEta_weighted[ibin][ibin2][ibin3]->Fill(me_trkEta->at(tracks),trkweight*wvz*wcen);
								my_hists[data_mc_type_code]->ME_TrkPhi_weighted[ibin][ibin2][ibin3]->Fill(me_trkPhi->at(tracks),trkweight*wvz*wcen);
								
								if(is_inclusive){
									
									deta = jet_eta - me_trkEta->at(tracks);
									dphi = jet_phi - me_trkPhi->at(tracks);
									
									while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
									while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
									
									my_hists[data_mc_type_code]->hJetTrackME[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
									my_hists[data_mc_type_code]->hJetTrackME_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, pt*trkweight*wvz*wcen);
                                                                        my_hists[data_mc_type_code]->hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
									if(foundBjet){
										my_hists[data_mc_type_code]->hbJetTrackME[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
										my_hists[data_mc_type_code]->hbJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
									}
									
								}
								
								
								if(j4i==highest_idx){
									deta = jteta->at(highest_idx) - me_trkEta->at(tracks);
									dphi = jtphi->at(highest_idx) - me_trkPhi->at(tracks);
									while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
									while(dphi< (-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
									
									my_hists[data_mc_type_code]->hJetTrackMELeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_lead*wvz*wcen);
									my_hists[data_mc_type_code]->hJetTrackMELeading_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, pt*trkweight_lead*wvz*wcen);
                                                                        my_hists[data_mc_type_code]->hJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
									if(foundBjet){
										my_hists[data_mc_type_code]->hbJetTrackMELeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_lead*wvz*wcen);
										my_hists[data_mc_type_code]->hbJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
									}
								}

								if(j4i==second_highest_idx){
									deta = jteta->at(second_highest_idx) - me_trkEta->at(tracks);
									dphi = jtphi->at(second_highest_idx) - me_trkPhi->at(tracks);
									while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
									while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}

									my_hists[data_mc_type_code]->hJetTrackMESubLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_sub*wvz*wcen);
									my_hists[data_mc_type_code]->hJetTrackMESubLeading_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, pt*trkweight_sub*wvz*wcen);
                                                                        my_hists[data_mc_type_code]->hJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
									if(foundBjet){
										my_hists[data_mc_type_code]->hbJetTrackMESubLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight_sub*wvz*wcen);
										my_hists[data_mc_type_code]->hbJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
									}
								}


								if(is_inclusive && j4i!=highest_idx){
									deta = jet_eta - me_trkEta->at(tracks);
									dphi = jet_phi - me_trkPhi->at(tracks);
									while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
									while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}

									my_hists[data_mc_type_code]->hJetTrackMENonLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
									my_hists[data_mc_type_code]->hJetTrackMENonLeading_pTweighted[ibin][ibin2][ibin3]->Fill(deta,dphi, pt*trkweight*wvz*wcen);
                                                                        my_hists[data_mc_type_code]->hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
									if(foundBjet){
										my_hists[data_mc_type_code]->hJetTrackMENonLeading[ibin][ibin2][ibin3]->Fill(deta,dphi, trkweight*wvz*wcen);
										my_hists[data_mc_type_code]->hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3]->Fill(deta,dphi, wvz*wcen);
									}
								}

							}  //track mixed event for data
                                                        //cout << "closed track loop" << endl;

						} //meptrig events per trigger
                                                //cout << "closed mixing loop" << endl;

                                            mixStopwatch.Stop();

					} //that version of mixing was just for data
                                        //cout << "finished mixing "<< endl;



				}   /// Closes jpti loop.  THIS MEANS THAT WE TAKE ALL JETS IN AN EVENT >120 GeV, not just the hardest jets.
                                //cout << "finished jets" << endl;
                                    
				if(foundjet==kTRUE){my_hists[data_mc_type_code]->NEvents_test->Fill(hiBin/2.);}


                  } //cent is a big loop

                        //    cout<<"here at end"<<endl;

            } ///we do EVERYTHING one event at a time.

            evtStopwatch.Stop();
            cout << "Avg event time: " << evtStopwatch.RealTime()/(double)n_evt << " sec/evt" << endl;
            cout << "Avg fg time: "<< fgStopwatch.RealTime()/(double)n_evt << " sec/evt" << endl;
            cout << "Avg mix time: "<< mixStopwatch.RealTime()/(double)n_evt << " sec/evt" << endl;
            cout << "mix time, part1: "<< mixStopwatch_p1.RealTime()/(double)n_evt << " sec/evt" << endl;
            cout << "mix time, part2: "<< mixStopwatch_p2.RealTime()/(double)n_evt << " sec/evt" << endl;
            cout << "mix time, part3: "<< mixStopwatch_p3.RealTime()/(double)n_evt << " sec/evt" << endl;
            cout << "total entry calls: " << entryCalls << endl;
            my_file->Close();
            
      }//FILE LOOP  (sort of a dummy, since we run on one file at a time).  
      
      cout<<"There were a total of "<<unmatched_counter<<" jets we could not match over all selections."<<endl;
      
      cout<<"Ready to write"<<endl;
      
      TString out_name;
      if(!is_data){
          if(!is_pp) out_name = (TString) (dataset_type_strs[globalCode+1] + "_PbPb.root");
          else out_name = (TString)"Pythia_pp.root";
      }
      else out_name = (TString) (dataset_type_strs[globalCode+1] + "_PbPb.root");
      TFile *out_file = new TFile(out_name, "RECREATE");
      out_file->cd();
      
      my_hists[globalCode]->Write(0);
      	
      out_file->Close();
      std::cout << "I am FINALLY done!!!" << std::endl;
      
      
} 


void ReadFileList(std::vector<std::string> &my_file_names, TString file_of_names, bool debug)
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
			std::string tstring_line(line);
			tstring_line.erase(std::remove(tstring_line.begin(), tstring_line.end(), '"'), tstring_line.end());
			tstring_line.erase(std::remove(tstring_line.begin(), tstring_line.end(), ','), tstring_line.end());
			tstring_line.erase(std::remove(tstring_line.begin(), tstring_line.end(), '['), tstring_line.end());
			tstring_line.erase(std::remove(tstring_line.begin(), tstring_line.end(), ']'), tstring_line.end());
			if( tstring_line != "" ) my_file_names.push_back(tstring_line);
			line_num++;
		}
	} else {
		std::cout << "Error, could not open " << file_of_names << " for reading" << std::endl;
		assert(0);
	}
}
