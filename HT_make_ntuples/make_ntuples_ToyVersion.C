#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"

#include "class_def/JetAna.h"
#include "class_def/Tracks.h"
#include "class_def/HLT.h"
#include "class_def/HiTree.h"
#include "class_def/Skim.h"
#include "class_def/GenParticles.h"
#include "class_def/pfcand.h"

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
using namespace std;



//arg 1 = which data set, arg 2 = output file number

enum enum_dataset_types {e_Data2011,e_Data_pp,e_HydJet15,e_HydJet30,e_HydJet50, e_HydJet80, e_HydJet120,e_HydJet170,e_HydJet220,e_HydJet280, e_HydJet370, e_n_dataset_types};
TString dataset_type_strs[e_n_dataset_types] = {"Data2011","Data_pp","HydJet15","HydJet30","HydJet50","HydJet80", "HydJet120", "HydJet170","HydJet220","HydJet280","HydJet370"};

int dataset_pthats[e_n_dataset_types+1] = {0,0,15,30,50,80,120,170,220,280,370,999};

int dataset_type_code = -999;

int main(int argc, char *argv[])
{
  assert( argc == 3);

  dataset_type_code = atoi(argv[1]);  

  bool is_data = false;
  
  if(dataset_type_code == e_Data2011 || dataset_type_code == e_Data_pp) is_data = true;
 
 
  // assert(!is_data); //for now I'm interested in MC

  int output_file_num = atoi(argv[2]);


  //////////###### centrality Bins ###########///////////////
  
  
  TTree *inp_tree;
  TTree *inp_tree2;
  TTree *inp_tree3;
  TTree *inp_tree4;
  TTree *inp_tree5;
  TTree *inp_tree6;
  TTree *inp_tree7;

  TString in_file_name;
 
  if(is_data){
    in_file_name = "/data/mzakaria/HI_53X_MB_Hydjet_STARTHI53_LV1_5_3_16_trk8_jet26_2Apr2014_1200EST_FOREST_0.root";

  }else{
    in_file_name = "/data/mzakaria/PythiaHydjet_OfficialProduction_20150211/HiForest_PYTHIA_HYDJET_pthat";
    in_file_name+=dataset_pthats[dataset_type_code];
    in_file_name+= "_Track9_Jet30_matchEqR_merged_forest_0.root";
  }
  TFile *my_file = TFile::Open(in_file_name);

  if(my_file->IsZombie()) { 
    std::cout << "Is zombie" << std::endl;
  }    

  inp_tree = (TTree*)  my_file->Get("akVs3CaloJetAnalyzer/t");
  JetAna *my_ct = new JetAna(inp_tree); 
  GenParticles *my_ct7;
 
  inp_tree7 = (TTree*) my_file->Get("HiGenParticleAna/hi"); 
  my_ct7 = new GenParticles(inp_tree7);
  
  std::cout << "Got CT" << std::endl;

  int n_evt = my_ct->fChain->GetEntriesFast();
  
    //MC
  TString output_file_base;
  output_file_base= "/data/htrauger/ToyMiniNtuples/PythiaOnly_pthat";
  output_file_base +=dataset_pthats[dataset_type_code];
  
  TString output_file_extension = "_p";   output_file_extension += output_file_num;   output_file_extension += ".root";
 
  TFile *output_file = new TFile((TString) (output_file_base + output_file_extension), "RECREATE");
  output_file->cd();
  TTree *mixing_tree = new TTree("mixing_tree", "");


  std::vector<Float_t> *pt = new std::vector<Float_t>();  pt->clear();
  std::vector<Float_t> *phi = new std::vector<Float_t>();  phi->clear();
  std::vector<Float_t> *eta = new std::vector<Float_t>();  eta->clear();
  
  std::vector<Float_t> *geneta = new std::vector<Float_t>();   geneta->clear();
  std::vector<Float_t> *genphi = new std::vector<Float_t>();   genphi->clear();
  std::vector<Float_t> *genpt = new std::vector<Float_t>();   genpt->clear();

  mixing_tree->Branch("pt", "vector<Float_t>", &pt);
  mixing_tree->Branch("phi", "vector<Float_t>", &phi);
  mixing_tree->Branch("eta", "vector<Float_t>", &eta);
   
  mixing_tree->Branch("geneta", "vector<Float_t>", &geneta);
  mixing_tree->Branch("genphi", "vector<Float_t>", &genphi);
  mixing_tree->Branch("genpt", "vector<Float_t>", &genpt);

  output_file->cd();
  
  int ev_min = 0;
  int   ev_max = 10000;

  cout << "ev_min: " << ev_min << ", Entries: " << n_evt << std::endl;

  assert( ev_min < n_evt );

  if( ev_max >= n_evt ) ev_max = n_evt;
  std::cout << "Will run from event number " << ev_min << " to " << ev_max << "\n";
  for(int evi = ev_min; evi < ev_max; evi++) {
      
    if( evi % 1000 == 0 )  std::cout << "evi: " << evi <<  " of " << n_evt << "\n";
  
    my_ct->fChain->GetEntry(evi);
    my_ct7->fChain->GetEntry(evi); 

    cout<<"got entry"<<endl;
  
    //   if( evi % 1000 == 0 ) std::cout << "Filled successfully" << std::endl;

    for(int j4i_gen = 0; j4i_gen < my_ct->ngen ; j4i_gen++) {

      if( fabs(my_ct->geneta[j4i_gen]) > 2 ) continue;
      if( my_ct->genpt[j4i_gen] < 30 ) continue;
	
      geneta->push_back(my_ct->geneta[j4i_gen]);
      genphi->push_back(my_ct->genphi[j4i_gen]);
      genpt->push_back(my_ct->genpt[j4i_gen]);
	
    } /// genjet loop

    //   cout<<"made it through jet loop"<<endl;


    for(int ipart=0;ipart<my_ct7->npart;ipart++){
 
      float temp_eta=my_ct7->eta[ipart];
      if(fabs(temp_eta)>2.4) continue; //acceptance of the tracker   

      float temp_pt= my_ct7->pt[ipart];
      if(temp_pt < 0.5) continue; //acceptance of the tracker
 
      if(my_ct7->sube[ipart]!=0) continue; //pythia only
	
      eta->push_back(my_ct7->eta[ipart]);
      phi->push_back(my_ct7->phi[ipart]);
      pt->push_back(my_ct7->pt[ipart]);
            
    }
    // cout<<"done particle loop"<<endl;  


    // cout<<eta->at(3)<<endl;
    ///// Fill it
    mixing_tree->Fill();

    cout<<"filled"<<endl;
    
    ///// Reset
  
    pt->clear();
    phi->clear();
    eta->clear();
          
    geneta->clear();
    genphi->clear();
    genpt->clear();
  
    // cout<<"done event "<<evi<<endl;

  }  ///event loop
   
  
  cout<<"writing"<<endl;
  
  //mixing_tree->Write();
  output_file->Write();
  output_file->Close();
  
  cout<<"done"<<endl;

}




