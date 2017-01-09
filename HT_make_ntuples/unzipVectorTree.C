
#include "TTree.h"
#include "TH1D.h"
#include "TFile.h"
#include <fstream>
#include <iostream>

using namespace std;

void unzipVectorTree(std::string filelist="filelist.txt", int istart, int iend){
	
	double pthatbins[10] = {15,30,50,80,120,170,220,280,370,9999};
	double xsecs[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
	double pthatEntries[9] = {0, 444180, 485654, 484515, 479933, 448829, 258987, 234835, 50644};
	
	TTree *ftree = new TTree("unzipMixTree","");
	
	double jtpt, jteta, refpt;
	int nCScand;
	ftree->Branch("jtpt", &jtpt);
	ftree->Branch("refpt", &refpt);
	ftree->Branch("jteta", &jteta);
	ftree->Branch("nCScand", &nCScand);
	ftree->Branch("weight", &weight);
	
	ifstream ifstr(filelist.c_str());
	string file;
	
	TChain *ch = new TChain("mixing_tree");
	vector<float> *t_jtpt=0, *t_refpt=0, *t_jteta=0;
	vector<int> *t_nCScand=0;
	float pthat;
	
	ch->SetBranchStatus("*",0);
	ch->SetBranchStatus("calo_jtpt",1);
	ch->SetBranchStatus("calo_refpt",1);
	ch->SetBranchStatus("calo_jteta",1);
	ch->SetBranchStatus("nCSpartGT2_id1",1);
	ch->SetBranchStatus("pthat",1);
	
	ch->SetBranchAddress("calo_jtpt",&t_jtpt);
	ch->SetBranchAddress("calo_refpt",&t_refpt);
	ch->SetBranchAddress("calo_jteta",&t_jteta);
	ch->SetBranchAddress("nCSpartGT2_id1",&t_nCScand);
	ch->SetBranchAddress("pthat",&pthat);
	
	for(int i=0; i<istart; i++){
		ifstr >> file;
	}
	for(int i=istart; i<iend; i++){
		ifstr >> file;
		ch->Add(file.c_str());
	}
	
	int entries = ch->GetEntries();
	for(int ientry=0; ientry<entries; ientry++){
		
		ch->GetEntry(ientry);
		for(int ijet=0; ijet<t_jtpt->size(); ijet++){
			jtpt = t_jtpt->at(ijet);
			jteta = t_jteta->at(ijet);
			refpt = t_refpt->at(ijet);
			nCScand = t_nCScand->at(ijet);
			
			int ibin=0;
			while(pthat>pthatbins[ibin+1]) ibin++;
			weight = (xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];

			ftree->Fill();
		}
	}
	
	TTree *fout = new TFile(Form("unzippedSkim_%d.root",iend));
	fout->cd();
	ftree->Write();
	fout->Close();
	
}