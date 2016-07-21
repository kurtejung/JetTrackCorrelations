#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "TF1.h"

#include "TH2F.h"
#include "TMath.h"
#include <TNtuple.h>
#include "TChain.h"
#include <TString.h>
#include <TCut.h>

#include "assert.h"
#include <fstream>
#include "TMath.h"
#include <vector>

#include "fragmentation_JEC.h"

using namespace std;


enum enum_dataset_types {e_Data2015,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370,e_Pythia15,e_Pythia30,e_Pythia50, e_Pythia80, e_Pythia120,e_Pythia170,e_Pythia220,e_Pythia280, e_Pythia370, e_n_dataset_types};
TString dataset_type_strs[e_n_dataset_types] = {"Data2015","Data_pp","HydJet15","HydJet30","HydJet50","HydJet80", "HydJet120", "HydJet170","HydJet220","HydJet280","HydJet370","Pythia15","Pythia30","Pythia50","Pythia80", "Pythia120", "Pythia170","Pythia220","Pythia280","Pythia370"};

int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,15,30,50,80,120,170,220,280,370,999};

int dataset_type_code = -999;

//arg 1 = which data set, arg 2 = output file number
void make_ntuples2(bool doCrab=0, int jobID=0, int endfile = 999, int dataset_type_code = 0 , int output_file_num = 1)
{

	bool is_data = false;

	if(dataset_type_code == 0 || dataset_type_code == 1) is_data = true;

	cout << "dataset code: " << dataset_type_code << endl;


	// assert(!is_data); //for now I'm interested in MC

	//-----------------------------
	// Set JFF-dependent corrections
	//-----------------------------

	float reco_eta, reco_phi, reco_pt, pfPt_temp, pfEta_temp, pfPhi_temp, pfId_temp, pfVsPt_temp;

	double corrected_pt, residual_corrected_pt, r;
	bool do_PbPb=1, do_pp_tracking=0;

	int radius = 4;
	bool do_residual_correction = kTRUE; 
	int nstep_residual = 3; 
	double Pf_pt_cut = 2;
	bool doFFCorrection = false;

	if(dataset_type_code== 1 || dataset_type_code > 10){do_PbPb = 0;   do_pp_tracking = 1;}

	cout<<"do_PbPb = "<<do_PbPb<<endl;
	cout<<"do_pp_tracking = "<<do_pp_tracking<<endl;


	fragmentation_JEC *FF_JEC;
	if(doFFCorrection) FF_JEC = new fragmentation_JEC(radius, do_PbPb, do_pp_tracking, do_residual_correction,nstep_residual,Pf_pt_cut); //3rd variable is only for when do_PbPb is false

	cout<<"made fragmentation JEC"<<endl;

	// For each event
	// Calculate efficiency
	if(doFFCorrection) FF_JEC->set_correction();


	vector<float> corr_pt; 

	//--------------------------------

	//////////###### centrality Bins ###########///////////////

	TTree *inp_tree;
	TTree *inp_tree2;
	TTree *inp_tree3;
	TTree *inp_tree4;
	TTree *inp_tree5;
	TTree *inp_tree6;
	TTree *inp_tree7=0;
	TTree *pftree;

	string in_file_name;

	if(doCrab){
		in_file_name = Form("job_input_file_list_%d.txt",jobID);
	}
	else if(is_data&&!do_PbPb){
		in_file_name = "HiForest_pp5TeV_eosDataset.txt";
	}else if(is_data&&do_PbPb){
		in_file_name = "PbPbForest_HiRun2015_MITForest.txt";
	}else if(dataset_type_code > 10){
		in_file_name = "/data/htrauger/Pythia_HiForest/HiForest_PYTHIA_pthat";
		in_file_name+=dataset_pthats[dataset_type_code];
		in_file_name+= ".root";
	}else if(dataset_type_code > 1&&dataset_type_code <11){
		in_file_name = "/data/mzakaria/PythiaHydjet_OfficialProduction_20150211/HiForest_PYTHIA_HYDJET_pthat";
		in_file_name+=dataset_pthats[dataset_type_code];
		in_file_name+= "_Track9_Jet30_matchEqR_merged_forest_0.root";
	}else{
		cerr<<"need to set up to run on that sample..."<<endl;
	}

	cout << "trying a filelist named "<< in_file_name << endl;

	//MC
	TString output_file_base;

	if(is_data && !do_PbPb){
		output_file_base= "./";
	}else if(is_data&&do_PbPb){
		output_file_base= "./";
	}else if(dataset_type_code > 1 &&dataset_type_code < 11){
		output_file_base= "/data/htrauger/OfficialHydjet_6_10/";
	}else if(dataset_type_code > 10){
		output_file_base= "/data/htrauger/OfficialPythia_6_24/";
	}else{
		cerr<<"nope, we can't handle that data set"<<endl;
		exit(0);
	}

	output_file_base +=dataset_type_strs[dataset_type_code];

	TString output_file_extension = "";   
	//output_file_extension += output_file_num;   
	output_file_extension += ".root";
	TFile *output_file = new TFile((TString) (output_file_base+output_file_extension), "RECREATE");
	TTree *mixing_tree = new TTree("mixing_tree", "");

	vector<float> pf_jteta, pf_jtphi, pf_jtpt, pf_rawpt, pf_corrpt, pf_trackMax;
	vector<float> calo_jteta, calo_jtphi, calo_jtpt, calo_rawpt, calo_corrpt, calo_trackMax;
	vector<int> trkAlgo;
	vector<bool> highPurity;
	vector<float> trkDxy1, trkDxyError1, trkDz1, trkDzError1, trkPtError, pfEcal, pfHcal, trkMVALoose, trkMVATight, trkChi2, trkEta, trkPhi, trkPt;
	vector<int> *pfId=0, trkNHit, trkNlayer, trkNdof;
	vector<float> *pfPt=0, *pfEta=0, *pfPhi=0, *pfVsPtInitial=0;
	vector<int> sube, chg;
	vector<float> pt, phi, eta, pPt, pPhi, pEta, geneta, genphi, genpt;
	vector<float> pf_discr_ssvHighEff, pf_discr_ssvHighPur, pf_discr_csvV1, pf_discr_prob, pf_svtxm, pf_svtxpt, pf_svtxmcorr, pf_svtxdl, pf_svtxdls;
	vector<float> calo_discr_ssvHighEff, calo_discr_ssvHighPur, calo_discr_csvV1, calo_discr_prob, calo_svtxm, calo_svtxpt, calo_svtxmcorr, calo_svtxdl, calo_svtxdls;

	Int_t HBHENoiseFilter, HBHENoiseFilterResultRun1, HBHENoiseFilterResultRun2Loose, HBHENoiseFilterResultRun2Tight, HBHEIsoNoiseFilterResult;
	Int_t phfCoincFilter1, phfCoincFilter2, phfCoincFilter3, phfCoincFilter4, phfCoincFilter5;
	Int_t eventSelection, pvFilter;
	Int_t HLT_Jet80, HLT_Jet100, HLT_Jet80_ps, HLT_Jet100_ps;
	Int_t pPAcollisionEventSelectionPA = -999;
	Int_t pClusterCompatibilityFilter, pVertexFilterCutGplus;
	Int_t hiBin = -999;
	Float_t pthat = -999;
	Float_t vz = -999;//, sumpt[15];

	/// higenparticles

	mixing_tree->Branch("HLT_ak4CaloJet80", &HLT_Jet80);
	mixing_tree->Branch("HLT_ak4CaloJet100", &HLT_Jet100);
	mixing_tree->Branch("HLT_ak4CaloJet80_Prescale", &HLT_Jet80_ps);
	mixing_tree->Branch("HLT_ak4CaloJet100_Prescale", &HLT_Jet100_ps);

	mixing_tree->Branch("HBHENoiseFilterResult",&HBHENoiseFilter);
	mixing_tree->Branch("HBHENoiseFilterResultRun1",&HBHENoiseFilterResultRun1);
	mixing_tree->Branch("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);
	mixing_tree->Branch("HBHENoiseFilterResultRun2Tight",&HBHENoiseFilterResultRun2Tight);
	mixing_tree->Branch("HBHEIsoNoiseFilterResult",&HBHEIsoNoiseFilterResult);

	if(!do_PbPb){
		mixing_tree->Branch("pPAprimaryVertexFilter",&pvFilter);
		mixing_tree->Branch("pVertexFilterCutGplus",&pVertexFilterCutGplus);
	}
	else{	
		mixing_tree->Branch("pcollisionEventSelection",&eventSelection);			
		mixing_tree->Branch("pprimaryVertexFilter",&pvFilter);
		mixing_tree->Branch("pclusterCompatibilityFilter",&pClusterCompatibilityFilter);
		mixing_tree->Branch("phfCoincFilter1",&phfCoincFilter1);
		mixing_tree->Branch("phfCoincFilter2",&phfCoincFilter2);
		mixing_tree->Branch("phfCoincFilter3",&phfCoincFilter3);
		mixing_tree->Branch("phfCoincFilter4",&phfCoincFilter4);
		mixing_tree->Branch("phfCoincFilter5",&phfCoincFilter5);
	}
	//mixing_tree->Branch("nTrk", &nTrk, "nTrk/I");
	mixing_tree->Branch("trkPt", &trkPt);
	mixing_tree->Branch("trkEta", &trkEta);
	mixing_tree->Branch("trkPhi", &trkPhi);
	mixing_tree->Branch("trkAlgo", &trkAlgo);
	mixing_tree->Branch("highPurity", &highPurity);
	mixing_tree->Branch("vz", &vz);
	//mixing_tree->Branch("pHBHENoiseFilter", &pHBHENoiseFilter, "pHBHENoiseFilter/I");
	//mixing_tree->Branch("pcollisionEventSelection", &pcollisionEventSelection, "pcollisionEventSelection/I");
	//mixing_tree->Branch("HLT_Jet80", &HLT_Jet80, "HLT_Jet80/I");
	//mixing_tree->Branch("pPAcollisionEventSelectionPA", &pPAcollisionEventSelectionPA, "pPAcollisionEventSelectionPA/I");

	mixing_tree->Branch("hiBin", &hiBin);

	mixing_tree->Branch("pf_jteta", &pf_jteta);
	mixing_tree->Branch("pf_jtphi", &pf_jtphi);
	mixing_tree->Branch("pf_jtpt", &pf_jtpt);
	mixing_tree->Branch("pf_rawpt", &pf_rawpt);
	mixing_tree->Branch("pf_corrpt", &pf_corrpt);
	mixing_tree->Branch("pf_trackMax", &pf_trackMax);

	mixing_tree->Branch("calo_jteta", &calo_jteta);
	mixing_tree->Branch("calo_jtphi", &calo_jtphi);
	mixing_tree->Branch("calo_jtpt", &calo_jtpt);
	mixing_tree->Branch("calo_rawpt", &calo_rawpt);
	mixing_tree->Branch("calo_corrpt", &calo_corrpt);
	mixing_tree->Branch("calo_trackMax", &calo_trackMax);

	if(!is_data) mixing_tree->Branch("pthat", &pthat, "pthat/F");
	mixing_tree->Branch("trkDxy", &trkDxy1);
	mixing_tree->Branch("trkDxyError", &trkDxyError1);
	mixing_tree->Branch("trkDz", &trkDz1);
	mixing_tree->Branch("trkDzError", &trkDzError1);
	mixing_tree->Branch("trkPtError", &trkPtError);
	mixing_tree->Branch("trkChi2", &trkChi2);
	mixing_tree->Branch("trkNdof", &trkNdof);
	mixing_tree->Branch("trkNHit",&trkNHit);
	mixing_tree->Branch("trkNlayer",&trkNlayer);
	mixing_tree->Branch("pfEcal",&pfEcal);
	mixing_tree->Branch("pfHcal",&pfHcal);
	mixing_tree->Branch("trkMVALoose",&trkMVALoose);
	mixing_tree->Branch("trkMVATight",&trkMVATight);

	if(!is_data){
		//mixing_tree->Branch("mult", &mult, "mult/I");
		mixing_tree->Branch("pt", &pt);
		mixing_tree->Branch("phi",  &phi);
		mixing_tree->Branch("eta",  &eta);
		mixing_tree->Branch("chg", &chg);
		mixing_tree->Branch("sube", &sube);
		//mixing_tree->Branch("nParticle", &nParticle, "nParticle/I");
		mixing_tree->Branch("pPt", &pPt);
		mixing_tree->Branch("pPhi", &pPhi);
		mixing_tree->Branch("pEta", &pEta);


		mixing_tree->Branch("geneta", &geneta);
		mixing_tree->Branch("genphi", &genphi);
		mixing_tree->Branch("genpt", &genpt);
	}

	//mixing_tree->Branch("nPFpart", &nPFpart, "nPFpart/I");
	mixing_tree->Branch("pfId", &pfId);
	mixing_tree->Branch("pfPt", &pfPt);
	mixing_tree->Branch("pfVsPtInitial", &pfVsPtInitial);
	mixing_tree->Branch("pfEta", &pfEta);
	mixing_tree->Branch("pfPhi", &pfPhi);
	//mixing_tree->Branch("sumpt", sumpt);

	//adding b-jet stuff....
	mixing_tree->Branch("pf_discr_ssvHighEff", &pf_discr_ssvHighEff);
	mixing_tree->Branch("pf_discr_ssvHighPur", &pf_discr_ssvHighPur);
	mixing_tree->Branch("pf_discr_csvV1", &pf_discr_csvV1);
	mixing_tree->Branch("pf_discr_prob", &pf_discr_prob);
	mixing_tree->Branch("pf_svtxm", &pf_svtxm);
	mixing_tree->Branch("pf_svtxpt", &pf_svtxpt);
	mixing_tree->Branch("pf_svtxmcorr", &pf_svtxmcorr);
	mixing_tree->Branch("pf_svtxdl", &pf_svtxdl);
	mixing_tree->Branch("pf_svtxdls", &pf_svtxdls);

	mixing_tree->Branch("calo_discr_ssvHighEff", &calo_discr_ssvHighEff);
	mixing_tree->Branch("calo_discr_ssvHighPur", &calo_discr_ssvHighPur);
	mixing_tree->Branch("calo_discr_csvV1", &calo_discr_csvV1);
	mixing_tree->Branch("calo_discr_prob", &calo_discr_prob);
	mixing_tree->Branch("calo_svtxm", &calo_svtxm);
	mixing_tree->Branch("calo_svtxpt", &calo_svtxpt);
	mixing_tree->Branch("calo_svtxmcorr", &calo_svtxmcorr);
	mixing_tree->Branch("calo_svtxdl", &calo_svtxdl);
	mixing_tree->Branch("calo_svtxdls", &calo_svtxdls);


	std::ifstream instr(in_file_name.c_str(), std::ifstream::in);
	if(!instr.is_open()) cout << "filelist not found!! Exiting..." << endl;
	std::string filename;
	int ifile=0;

	const int MAXJETS = 500;
	Float_t t_pf_jtpt[MAXJETS], t_pf_jteta[MAXJETS], t_pf_jtphi[MAXJETS], t_pf_discr_ssvHighEff[MAXJETS], t_pf_discr_ssvHighPur[MAXJETS], t_pf_discr_csvV1[MAXJETS], t_pf_discr_prob[MAXJETS], t_pf_svtxm[MAXJETS], t_pf_svtxpt[MAXJETS], t_pf_svtxmcorr[MAXJETS], t_pf_svtxdl[MAXJETS], t_pf_svtxdls[MAXJETS], t_pf_rawpt[MAXJETS], t_pf_trackMax[MAXJETS];
	Float_t t_calo_jtpt[MAXJETS], t_calo_jteta[MAXJETS], t_calo_jtphi[MAXJETS], t_calo_discr_ssvHighEff[MAXJETS], t_calo_discr_ssvHighPur[MAXJETS], t_calo_discr_csvV1[MAXJETS], t_calo_discr_prob[MAXJETS], t_calo_svtxm[MAXJETS], t_calo_svtxpt[MAXJETS], t_calo_svtxmcorr[MAXJETS], t_calo_svtxdl[MAXJETS], t_calo_svtxdls[MAXJETS], t_calo_rawpt[MAXJETS], t_calo_trackMax[MAXJETS];

	const int MAXPARTICLES = 60000;
	Float_t t_trkPt[MAXPARTICLES], t_trkEta[MAXPARTICLES], t_trkPhi[MAXPARTICLES], t_trkDxy1[MAXPARTICLES], t_trkDxyError1[MAXPARTICLES], t_trkDz1[MAXPARTICLES], t_trkDzError1[MAXPARTICLES], t_trkPtError[MAXPARTICLES], t_pfEcal[MAXPARTICLES], t_pfHcal[MAXPARTICLES], t_trkChi2[MAXPARTICLES];
	Bool_t t_trkMVALoose[MAXPARTICLES], t_trkMVATight[MAXPARTICLES];
	UChar_t t_trkAlgo[MAXPARTICLES], t_trkNHit[MAXPARTICLES], t_trkNlayer[MAXPARTICLES], t_trkNdof[MAXPARTICLES];
	Bool_t t_highPurity[MAXPARTICLES];

	Int_t nTrk, calo_nref, pf_nref;

	while(instr>>filename && ifile<endfile){
		filename.erase(std::remove(filename.begin(), filename.end(), '"'), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), ','), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), '['), filename.end());
		filename.erase(std::remove(filename.begin(), filename.end(), ']'), filename.end());
		cout<<"File name is "<< filename <<endl;
		ifile++;

		TFile *my_file = TFile::Open(filename.c_str());
		if(!my_file){
			int pos = filename.find_first_of('s');
			string reducedfn = filename.substr(pos+1);
			string xrdPrefix = "root://xrootd.unl.edu//";
			TFile::Open((xrdPrefix+reducedfn).c_str());
		}
		
		if(!my_file){ cout << "File cannot be found!!" << endl; exit(1); }	

		if(my_file->IsZombie()) { 
			std::cout << "Is zombie" << std::endl;
		}    

		if(do_PbPb){
			inp_tree = (TTree*)  my_file->Get(Form("akPu%dCaloJetAnalyzer/t",radius));
			pftree = (TTree*) my_file->Get(Form("akPu%dPFJetAnalyzer/t",radius));
		}else{
			inp_tree = (TTree*)  my_file->Get(Form("ak%dCaloJetAnalyzer/t",radius));
			pftree = (TTree*) my_file->Get(Form("ak%dPFJetAnalyzer/t",radius));
		}

		inp_tree2 = (TTree*)  my_file->Get("pfcandAnalyzer/pfTree");
		if(!inp_tree2){ cout << "PFCand Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree2);

		inp_tree3 = (TTree*) my_file->Get("hiEvtAnalyzer/HiTree");
		if(!inp_tree3){ cout << "Evt Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree3);

		inp_tree4 = (TTree*) my_file->Get("skimanalysis/HltTree");
		if(!inp_tree4){ cout << "Skim Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree4);

		if(do_PbPb){
			inp_tree5 = (TTree*) my_file->Get("anaTrack/trackTree");
		}else{
			inp_tree5 = (TTree*) my_file->Get("ppTrack/trackTree");
		}
		if(!inp_tree5){ cout << "track Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree5);

		inp_tree6 = (TTree*) my_file->Get("hltanalysis/HltTree");
		if(!inp_tree6){ cout << "HLT Tree not found!! Exiting..." << endl; exit(1); }
		else inp_tree->AddFriend(inp_tree6);

		if(!is_data){ 
			inp_tree7 = (TTree*) my_file->Get("HiGenParticleAna/hi");
			if(!inp_tree7){ cout << "GenPart Tree not found!! Exiting..." << endl; exit(1); }
			else inp_tree->AddFriend(inp_tree7);
		}

		inp_tree->SetBranchAddress("trkPt",t_trkPt);
		inp_tree->SetBranchAddress("trkEta",t_trkEta);
		inp_tree->SetBranchAddress("trkPhi",t_trkPhi);

		inp_tree->SetBranchAddress("trkAlgo",t_trkAlgo);
		inp_tree->SetBranchAddress("highPurity",t_highPurity);
		inp_tree->SetBranchAddress("vz",&vz);
		inp_tree->SetBranchAddress("hiBin",&hiBin);

		inp_tree->SetBranchAddress("nref",&calo_nref);
		inp_tree->SetBranchAddress("jtpt",t_calo_jtpt);
		inp_tree->SetBranchAddress("jteta",t_calo_jteta);
		inp_tree->SetBranchAddress("jtphi",t_calo_jtphi);
		inp_tree->SetBranchAddress("rawpt",t_calo_rawpt);
		inp_tree->SetBranchAddress("trackMax",t_calo_trackMax);

		pftree->SetBranchAddress("nref",&pf_nref);
		pftree->SetBranchAddress("jtpt",t_pf_jtpt);
		pftree->SetBranchAddress("jteta",t_pf_jteta);
		pftree->SetBranchAddress("jtphi",t_pf_jtphi);
		pftree->SetBranchAddress("rawpt",t_pf_rawpt);
		pftree->SetBranchAddress("trackMax",t_pf_trackMax);

		inp_tree->SetBranchAddress("nTrk",&nTrk);
		inp_tree->SetBranchAddress("trkDxy1",t_trkDxy1);
		inp_tree->SetBranchAddress("trkDxyError1",t_trkDxyError1);
		inp_tree->SetBranchAddress("trkDz1",t_trkDz1);
		inp_tree->SetBranchAddress("trkDzError1",t_trkDzError1);
		inp_tree->SetBranchAddress("trkPtError",t_trkPtError);
		inp_tree->SetBranchAddress("trkNHit",t_trkNHit);
		inp_tree->SetBranchAddress("trkNlayer",t_trkNlayer);
		inp_tree->SetBranchAddress("trkChi2",t_trkChi2);
		inp_tree->SetBranchAddress("trkNdof",t_trkNdof);
		inp_tree->SetBranchAddress("pfEcal",t_pfEcal);
		inp_tree->SetBranchAddress("pfHcal",t_pfHcal);
		if(!do_PbPb) inp_tree->SetBranchAddress("trkMVALoose",t_trkMVALoose);
		inp_tree->SetBranchAddress("trkMVATight",t_trkMVATight);

		inp_tree->SetBranchAddress("pfId",&pfId);
		inp_tree->SetBranchAddress("pfPt",&pfPt);
		inp_tree->SetBranchAddress("pfEta",&pfEta);
		inp_tree->SetBranchAddress("pfPhi",&pfPhi);
		inp_tree->SetBranchAddress("pfVsPtInitial",&pfVsPtInitial);
		//inp_tree->SetBranchAddress("sumpt",sumpt);	

		inp_tree->SetBranchAddress("discr_ssvHighEff", t_calo_discr_ssvHighEff);
		inp_tree->SetBranchAddress("discr_ssvHighPur", t_calo_discr_ssvHighPur);
		if(!do_PbPb) inp_tree->SetBranchAddress("discr_csvV1", t_calo_discr_csvV1);
		else{
			inp_tree->SetBranchAddress("discr_csvSimple", t_calo_discr_csvV1);
		}
		inp_tree->SetBranchAddress("discr_prob", t_calo_discr_prob);
		inp_tree->SetBranchAddress("svtxdl", t_calo_svtxdl);
		inp_tree->SetBranchAddress("svtxdls", t_calo_svtxdls);
		inp_tree->SetBranchAddress("svtxm", t_calo_svtxm);
		inp_tree->SetBranchAddress("svtxpt", t_calo_svtxpt);
		inp_tree->SetBranchAddress("svtxmcorr", t_calo_svtxmcorr);

		pftree->SetBranchAddress("discr_ssvHighEff", t_pf_discr_ssvHighEff);
		pftree->SetBranchAddress("discr_ssvHighPur", t_pf_discr_ssvHighPur);
		if(!do_PbPb) pftree->SetBranchAddress("discr_csvV1", t_pf_discr_csvV1);
		else{
			pftree->SetBranchAddress("discr_csvSimple", t_pf_discr_csvV1);
		}
		pftree->SetBranchAddress("discr_prob", t_pf_discr_prob);
		pftree->SetBranchAddress("svtxdl", t_pf_svtxdl);
		pftree->SetBranchAddress("svtxdls", t_pf_svtxdls);
		pftree->SetBranchAddress("svtxm", t_pf_svtxm);
		pftree->SetBranchAddress("svtxpt", t_pf_svtxpt);
		pftree->SetBranchAddress("svtxmcorr", t_pf_svtxmcorr);

		inp_tree->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilter);
		inp_tree->SetBranchAddress("HBHENoiseFilterResultRun1",&HBHENoiseFilterResultRun1);
		inp_tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);
		inp_tree->SetBranchAddress("HBHENoiseFilterResultRun2Tight",&HBHENoiseFilterResultRun2Tight);
		inp_tree->SetBranchAddress("HBHEIsoNoiseFilterResult",&HBHEIsoNoiseFilterResult);

		if(!do_PbPb){
			inp_tree->SetBranchAddress("pPAprimaryVertexFilter",&pvFilter);
			inp_tree->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutGplus);
		}
		else{	
			inp_tree->SetBranchAddress("pcollisionEventSelection",&eventSelection);			
			inp_tree->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
			inp_tree->SetBranchAddress("pclusterCompatibilityFilter",&pClusterCompatibilityFilter);
			inp_tree->SetBranchAddress("phfCoincFilter1",&phfCoincFilter1);
			inp_tree->SetBranchAddress("phfCoincFilter2",&phfCoincFilter2);
			inp_tree->SetBranchAddress("phfCoincFilter3",&phfCoincFilter3);
			inp_tree->SetBranchAddress("phfCoincFilter4",&phfCoincFilter4);
			inp_tree->SetBranchAddress("phfCoincFilter5",&phfCoincFilter5);
		}

		if(!do_PbPb){ 
			inp_tree->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1", &HLT_Jet80);
			inp_tree->SetBranchAddress("HLT_AK4PFJet100_Eta5p1_v1", &HLT_Jet100);
		}
		else{
			inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1", &HLT_Jet80);
			inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1", &HLT_Jet100);
			inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet80_Eta5p1_v1_Prescl", &HLT_Jet80_ps);
			inp_tree->SetBranchAddress("HLT_HIPuAK4CaloJet100_Eta5p1_v1_Prescl", &HLT_Jet100_ps);
		}


		int n_evt = inp_tree->GetEntriesFast();

		cout << "Entries: "<< n_evt << endl;
		for(int evi = 0; evi < n_evt; evi++) {

			inp_tree->GetEntry(evi);
			pftree->GetEntry(evi);

			if(!HLT_Jet80 && !HLT_Jet100) continue;
			if(!do_PbPb && (!pvFilter || !HBHENoiseFilter)) continue;

			// Removing event selection for PbPb skims so we can test various selection criteria
			//if(do_PbPb && (!pvFilter || !HBHENoiseFilter || !eventSelection)) continue;

			//if( evi % 1000 == 0 ) std::cout << "Filled successfully" << std::endl;

			//start calo jet loop
			for(int j4i = 0; j4i < calo_nref ; j4i++) {

				if( fabs(t_calo_jteta[j4i]) > 2. ) continue;

				//-----------------------------------------------------------------------------------------
				// Jet Energy Corrections (JFF-dependent, store final corrected values in vector corr_pt)
				//----------------------------------------------------------------------------------------


				reco_pt = t_calo_jtpt[j4i];
				reco_phi = t_calo_jtphi[j4i];
				reco_eta = t_calo_jteta[j4i];


				int npf=0;

				for(unsigned int ipf=0;ipf< pfPt->size(); ipf++){

					pfPt_temp = pfPt->at(ipf);
					pfVsPt_temp = pfVsPtInitial->at(ipf);
					pfEta_temp =  pfEta->at(ipf);
					pfPhi_temp =  pfPhi->at(ipf);
					pfId_temp = pfId->at(ipf);  //pfId == 1 for hadrons only

					//cout<<pfPt<<" "<<pfVsPt<<" "<<pfEta<<" "<<pfPhi<<" "<<pfId<<endl;
					r=sqrt(pow(reco_eta-pfEta_temp,2)+pow(acos(cos(reco_phi-pfPhi_temp)),2));

					if(do_PbPb&&r<((double)(radius)*0.1)&& pfVsPt_temp > Pf_pt_cut && pfEta_temp <2.4 && pfId_temp==1) npf++; 

					if(!do_PbPb&&r<((double)(radius)*0.1)&& pfPt_temp > Pf_pt_cut && pfEta_temp <2.4 && pfId_temp==1) npf++; 
				}

				if(do_PbPb){ 
					corrected_pt= reco_pt; //FF_JEC->get_corrected_pt(reco_pt, npf, hiBin);
				}else{ 
					corrected_pt= reco_pt; //FF_JEC->get_corrected_pt(reco_pt, npf);
				}


				if(do_residual_correction){ 
					residual_corrected_pt= reco_pt; //FF_JEC->get_residual_corrected_pt(corrected_pt,hiBin);
				}

				if( t_calo_jtpt[j4i] < 25 && residual_corrected_pt < 25) continue;

				calo_jteta.push_back(reco_eta);
				calo_jtphi.push_back(reco_phi);
				calo_jtpt.push_back(reco_pt);
				calo_corrpt.push_back(reco_pt); //residual_corrected_pt);
				calo_rawpt.push_back(t_calo_rawpt[j4i]);

				calo_discr_ssvHighEff.push_back(t_calo_discr_ssvHighEff[j4i]);
				calo_discr_ssvHighPur.push_back(t_calo_discr_ssvHighPur[j4i]);
				calo_discr_csvV1.push_back(t_calo_discr_csvV1[j4i]);
				calo_discr_prob.push_back(t_calo_discr_prob[j4i]);
				calo_svtxm.push_back(t_calo_svtxm[j4i]);
				calo_svtxpt.push_back(t_calo_svtxpt[j4i]);
				calo_svtxmcorr.push_back(t_calo_svtxmcorr[j4i]);
				calo_svtxdl.push_back(t_calo_svtxdl[j4i]);
				calo_svtxdls.push_back(t_calo_svtxdls[j4i]);

				calo_trackMax.push_back(t_calo_trackMax[j4i]);

			} /// calo jet loop


			//start pf jet loop
			for(int j4i = 0; j4i < pf_nref ; j4i++) {

				if( fabs(t_pf_jteta[j4i]) > 2. ) continue;

				//-----------------------------------------------------------------------------------------
				// Jet Energy Corrections (JFF-dependent, store final corrected values in vector corr_pt)
				//----------------------------------------------------------------------------------------


				reco_pt = t_pf_jtpt[j4i];
				reco_phi = t_pf_jtphi[j4i];
				reco_eta = t_pf_jteta[j4i];


				int npf=0;

				for(unsigned int ipf=0;ipf< pfPt->size(); ipf++){

					pfPt_temp = pfPt->at(ipf);
					pfVsPt_temp = pfVsPtInitial->at(ipf);
					pfEta_temp =  pfEta->at(ipf);
					pfPhi_temp =  pfPhi->at(ipf);
					pfId_temp = pfId->at(ipf);  //pfId == 1 for hadrons only

					//cout<<pfPt<<" "<<pfVsPt<<" "<<pfEta<<" "<<pfPhi<<" "<<pfId<<endl;
					r=sqrt(pow(reco_eta-pfEta_temp,2)+pow(acos(cos(reco_phi-pfPhi_temp)),2));

					if(do_PbPb&&r<((double)(radius)*0.1)&& pfVsPt_temp > Pf_pt_cut && pfEta_temp <2.4 && pfId_temp==1) npf++; 

					if(!do_PbPb&&r<((double)(radius)*0.1)&& pfPt_temp > Pf_pt_cut && pfEta_temp <2.4 && pfId_temp==1) npf++; 
				}

				if(do_PbPb){ 
					corrected_pt= reco_pt; //FF_JEC->get_corrected_pt(reco_pt, npf, hiBin);
				}else{ 
					corrected_pt= reco_pt; //FF_JEC->get_corrected_pt(reco_pt, npf);
				}


				if(do_residual_correction){ 
					residual_corrected_pt= reco_pt; //FF_JEC->get_residual_corrected_pt(corrected_pt,hiBin);
				}

				if( t_pf_jtpt[j4i] < 25 && residual_corrected_pt < 25) continue;

				pf_jteta.push_back(reco_eta);
				pf_jtphi.push_back(reco_phi);
				pf_jtpt.push_back(reco_pt);
				pf_corrpt.push_back(reco_pt); //residual_corrected_pt);
				pf_rawpt.push_back(t_pf_rawpt[j4i]);

				pf_discr_ssvHighEff.push_back(t_pf_discr_ssvHighEff[j4i]);
				pf_discr_ssvHighPur.push_back(t_pf_discr_ssvHighPur[j4i]);
				pf_discr_csvV1.push_back(t_pf_discr_csvV1[j4i]);
				pf_discr_prob.push_back(t_pf_discr_prob[j4i]);
				pf_svtxm.push_back(t_pf_svtxm[j4i]);
				pf_svtxpt.push_back(t_pf_svtxpt[j4i]);
				pf_svtxmcorr.push_back(t_pf_svtxmcorr[j4i]);
				pf_svtxdl.push_back(t_pf_svtxdl[j4i]);
				pf_svtxdls.push_back(t_pf_svtxdls[j4i]);

				pf_trackMax.push_back(t_pf_trackMax[j4i]);

			} /// pf jet loop


			if(!is_data){


				/* for(int j4i_gen = 0; j4i_gen < ngen ; j4i_gen++) {

				   if( fabs(geneta[j4i_gen]) > 2 ) continue;
				   if( genpt[j4i_gen] < 30 ) continue;

				   geneta.push_back(geneta[j4i_gen]);
				   genphi.push_back(genphi[j4i_gen]);
				   genpt.push_back(genpt[j4i_gen]);

				   } /// genjet loop*/

			}

			//// reco track loop
			for(int itrk=0;itrk<nTrk;itrk++){

				//very basic cuts

				if(t_trkPtError[itrk]/t_trkPt[itrk]>=0.1 || TMath::Abs(t_trkDz1[itrk]/t_trkDzError1[itrk])>=3.0 ||TMath::Abs(t_trkDxy1[itrk]/t_trkDxyError1[itrk])>=3.0) continue ;

				float eta=t_trkEta[itrk];
				if(fabs(eta)>2.4) continue; //acceptance of the tracker   

				float pt=t_trkPt[itrk];
				if(pt < 0.5) continue; //pt min

				// reco track quantities

				if(!is_data){
					/* pEta.push_back(.pEta[itrk]);
					   pPhi.push_back(.pPhi[itrk]);
					   pPt.push_back(.pPt[itrk]);*/

				}

				trkEta.push_back(t_trkEta[itrk]);
				trkPhi.push_back(t_trkPhi[itrk]);
				trkPt.push_back(t_trkPt[itrk]);
				trkAlgo.push_back(t_trkAlgo[itrk]);
				highPurity.push_back((bool)t_highPurity[itrk]);

				trkDxy1.push_back(t_trkDxy1[itrk]);
				trkDxyError1.push_back(t_trkDxyError1[itrk]);
				trkDz1.push_back(t_trkDz1[itrk]);
				trkDzError1.push_back(t_trkDzError1[itrk]);
				trkPtError.push_back(t_trkPtError[itrk]);

				trkNHit.push_back((int)t_trkNHit[itrk]);
				trkNlayer.push_back((int)t_trkNlayer[itrk]);
				trkChi2.push_back(t_trkChi2[itrk]);
				trkNdof.push_back((int)t_trkNdof[itrk]);
				pfEcal.push_back(t_pfEcal[itrk]);
				pfHcal.push_back(t_pfHcal[itrk]);
				trkMVALoose.push_back(t_trkMVALoose[itrk]);
				trkMVATight.push_back(t_trkMVATight[itrk]);

			}

			if(!is_data){

				//gen particles loop
				/* for(int ipart=0;ipart<mult;ipart++){

				   float temp_eta=t_eta[ipart];
				   if(fabs(temp_eta)>2.4) continue; //acceptance of the tracker   

				   float temp_pt= t_pt[ipart];
				   if(temp_pt < 0.5) continue; //acceptance of the tracker

				// reco track quantities

				eta.push_back(t_eta[ipart]);
				phi.push_back(t_phi[ipart]);
				pt.push_back(t_pt[ipart]);
				chg.push_back(t_chg[ipart]);
				sube.push_back(t_sube[ipart]);

				}*/
			}

			//pthat = pthat;

			///// Fill it
			mixing_tree->Fill();

			///// Reset
			trkEta.clear();
			trkPhi.clear();
			trkPt.clear();
			trkAlgo.clear();
			highPurity.clear();

			calo_jteta.clear();
			calo_jtphi.clear();
			calo_jtpt.clear();
			calo_corrpt.clear();
			calo_rawpt.clear();
			calo_trackMax.clear();

			pf_jteta.clear();
			pf_jtphi.clear();
			pf_jtpt.clear();
			pf_corrpt.clear();
			pf_rawpt.clear();
			pf_trackMax.clear();

			trkDxy1.clear();
			trkDxyError1.clear();
			trkDz1.clear();
			trkDzError1.clear();
			trkPtError.clear();
			trkChi2.clear();
			trkNdof.clear();
			pfEcal.clear();
			pfHcal.clear();
			trkMVATight.clear();
			trkMVALoose.clear();
			trkNHit.clear();
			trkNlayer.clear();

			pfId->clear();
			pfPt->clear();
			pfVsPtInitial->clear();
			pfEta->clear();
			pfPhi->clear();

			calo_discr_ssvHighEff.clear();
			calo_discr_ssvHighPur.clear();
			calo_discr_csvV1.clear();
			calo_discr_prob.clear();
			calo_svtxm.clear();
			calo_svtxpt.clear();
			calo_svtxmcorr.clear();
			calo_svtxdl.clear();
			calo_svtxdls.clear();

			pf_discr_ssvHighEff.clear();
			pf_discr_ssvHighPur.clear();
			pf_discr_csvV1.clear();
			pf_discr_prob.clear();
			pf_svtxm.clear();
			pf_svtxpt.clear();
			pf_svtxmcorr.clear();
			pf_svtxdl.clear();
			pf_svtxdls.clear();

			if(!is_data){
				pt.clear();
				phi.clear();
				eta.clear();
				chg.clear();
				sube.clear();

				pPt.clear();
				pPhi.clear();
				pEta.clear();


				geneta.clear();
				genphi.clear();
				genpt.clear();
			}

		}  ///event loop

		my_file->Close();
	}


	cout<<"writing"<<endl;

	//mixing_tree->Write();
	output_file->cd();
	mixing_tree->Write();

	output_file->Close();

	cout<<"done"<<endl;

}




