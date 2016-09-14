
#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"

#include <vector>
#include <iostream>

using namespace std;

void formatHisto(TH1D *hist, int color){

	hist->SetMarkerStyle(20);
	hist->SetLineColor(color);
	hist->SetMarkerColor(color);

}

void drawQGFraction(){

	const double PI = 3.14159;

	TChain *mix = new TChain("mixing_tree");
	mix->Add("/data/kurtjung/JetTrackCorr_skims/2p76TeV_MC_Pythia6/MergedPythia_withPartonFlavor.root");
	//mix->Add("/data/kurtjung/JetTrackCorr_skims/5TeV_MC_Pythia6/*");

	//double xsecs[11] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 3.216E-08, 0}; //pythia6 5.02 tev weights
    double xsecs[11] = {2.043e-01, 1.075E-02, 1.025E-03, 9.865E-05, 1.129E-05, 1.465E-06, 2.837E-07, 5.323E-08, 5.934e-09, 8.125e-10, 0}; //2.76 tev weights
	int recalculatedEntries[10] = {0,0,0,0,0,0,0,0,0,0};
	double pthatbins[11] = {15,30,50,80,120,170,220,280,370,460,9999};

	TFile *fout = new TFile("QGFrac_pythia6_2p76TeV.root","recreate");

	TH1D *quarkFracIncl = new TH1D("quarkFracIncl","",20,120,500); quarkFracIncl->Sumw2();
	TH1D *glueFracIncl = new TH1D("glueFracIncl","",20,120,500); glueFracIncl->Sumw2();
	TH1D *inclJets = new TH1D("inclJets","",20,120,500); inclJets->Sumw2();

	TH1D *quarkFracLead = new TH1D("quarkFracLead","",20,120,500); quarkFracLead->Sumw2();
	TH1D *glueFracLead = new TH1D("glueFracLead","",20,120,500); glueFracLead->Sumw2();
	TH1D *leadJets = new TH1D("leadJets","",20,120,500); leadJets->Sumw2();


	Int_t HBHENoiseFilterResultRun2Loose, pPAprimaryVertexFilter;
	vector<float> *calo_corrpt=0, *calo_jtphi=0;
	vector<int> *calo_refparton_flavor=0;
	float pthat;

	mix->SetBranchAddress("calo_corrpt",&calo_corrpt);
	mix->SetBranchAddress("calo_jtphi",&calo_jtphi);
	mix->SetBranchAddress("calo_refparton_flavor",&calo_refparton_flavor);
	mix->SetBranchAddress("pthat",&pthat);

	mix->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHENoiseFilterResultRun2Loose);
	mix->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);


	TH1D *pthatHisto = new TH1D("pthatHisto","",10,pthatbins);
	mix->Project("pthatHisto","pthat");
	for(int i=0; i<10; i++){
		recalculatedEntries[i] = pthatHisto->GetBinContent(i+1);
		cout << "entries between pthat " << pthatbins[i] << " and " << pthatbins[i+1] << ": " << recalculatedEntries[i] << endl;
		cout << "weight: " << (xsecs[i]-xsecs[i+1])/recalculatedEntries[i] <<endl;
	}

	int totEntries = mix->GetEntries();
	cout << "entries: "<< totEntries << endl;
	totEntries=100000;
	for(int ievt=0; ievt<totEntries; ievt++){
		mix->GetEntry(ievt);
		if(ievt && ievt%10000==0) cout << "entry: " << ievt << endl;

		//if(!HBHENoiseFilterResultRun2Loose || !pPAprimaryVertexFilter) continue;

		int ibin=0;
		double weight=0.;
		while(pthat>pthatbins[ibin]) ibin++;
		ibin--;
		weight = (xsecs[ibin]-xsecs[ibin+1])/recalculatedEntries[ibin];
		if(weight>1){ 
			cout << "xsec: "<< xsecs[ibin]-xsecs[ibin+1] << " entries: " << recalculatedEntries[ibin] << endl;
			cout << "pthat: "<< pthat << " bin " << ibin << endl;
		}
		
		for(unsigned int ijet=0; ijet<calo_corrpt->size(); ijet++){
			
			if(abs(calo_refparton_flavor->at(ijet))>0 && abs(calo_refparton_flavor->at(ijet))<6) quarkFracIncl->Fill(calo_corrpt->at(ijet), weight);
			if(abs(calo_refparton_flavor->at(ijet))==21){ glueFracIncl->Fill(calo_corrpt->at(ijet), weight); }
			inclJets->Fill(calo_corrpt->at(ijet), weight);

			if(ijet==0 && calo_corrpt->size()>1){
				double dphi = abs(calo_jtphi->at(0) - calo_jtphi->at(1));
				if(dphi>(7*PI/8.) && dphi<(9*PI/8.) && calo_corrpt->at(1)>50) {

					if(abs(calo_refparton_flavor->at(ijet))>0 && abs(calo_refparton_flavor->at(ijet))<6) quarkFracLead->Fill(calo_corrpt->at(ijet), weight);
					if(abs(calo_refparton_flavor->at(ijet))==21) glueFracLead->Fill(calo_corrpt->at(ijet), weight);
					leadJets->Fill(calo_corrpt->at(ijet), weight);

				}
			}
		}
	}

	//quarkFracLead->Divide(leadJets);
	//glueFracLead->Divide(leadJets);

	//quarkFracIncl->Divide(inclJets);
	//glueFracIncl->Divide(inclJets);

	fout->cd();

	formatHisto(quarkFracLead,1);
	formatHisto(glueFracLead,2);
	formatHisto(quarkFracIncl,4);
	formatHisto(glueFracIncl,8);

	quarkFracLead->Write();
	glueFracLead->Write();
	quarkFracIncl->Write();
	glueFracIncl->Write();
	leadJets->Write();
	inclJets->Write();

	fout->Close();

}