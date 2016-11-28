
#include "TF1.h"
#include "TH1D.h"
#include "TChain.h"
#include "TFile.h"

void deriveVzCentWeight(){

	const int nCentBins = 4;
	const int centBinScheme[nCentBins+1] = {0,10,30,50,100};
	
	TH1D *mcVz[4], *dataVz[4];
	TH1D *mcVzPthat[4];
	TH1D *mcVzPthatJetOnly[4], *dataVzPthatJetOnly[4];
	for(int i=0; i<4; i++){
		mcVz[i] = new TH1D(Form("mcVz_cent%d",i),"",60,-15,15); mcVz[i]->Sumw2();
		dataVz[i] = new TH1D(Form("dataVz_cent%d",i),"",60,-15,15); dataVz[i]->Sumw2();
		mcVzPthat[i] = new TH1D(Form("mcVzPthat_cent%d",i),"",60,-15,15); mcVzPthat[i]->Sumw2();
		mcVzPthatJetOnly[i] = new TH1D(Form("mcVzPthatJetOnly_cent%d",i),"",60,-15,15); mcVzPthatJetOnly[i]->Sumw2();
		dataVzPthatJetOnly[i] = new TH1D(Form("dataVzPthatJetOnly_cent%d",i),"",60,-15,15); dataVzPthatJetOnly[i]->Sumw2();
	}
	TH1D *mcCent = new TH1D("mcCent","",200,0,200); mcCent->Sumw2();
	TH1D *dataCent = new TH1D("dataCent","",200,0,200); dataCent->Sumw2();

	TChain *mc = new TChain("mixing_tree");
	mc->Add("/data/kurtjung/JetTrackCorr_skims/5TeV_MC_PbPb_Pythia6Hydjet/looseMerge/*.root");
	mc->SetBranchStatus("*",0);
	mc->SetBranchStatus("vz",1);
	mc->SetBranchStatus("hiBin",1);
	mc->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
	mc->SetBranchStatus("pcollisionEventSelection",1);
	mc->SetBranchStatus("phfCoincFilter3",1);
	mc->SetBranchStatus("pthat",1);
	mc->SetBranchStatus("calo_corrpt",1);

	TChain *dta = new TChain("mixing_tree");
	dta->Add("/data/kurtjung/JetTrackCorr_skims/5TeV_PbPb/looseMerge/*.root");
	dta->SetBranchStatus("*",0);
	dta->SetBranchStatus("vz",1);
	dta->SetBranchStatus("hiBin",1);
	dta->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
	dta->SetBranchStatus("pcollisionEventSelection",1);
	dta->SetBranchStatus("phfCoincFilter3",1);
	dta->SetBranchStatus("calo_corrpt",1);

	float vz, pthat;
	int hiBin;
	int HBHEFilter, collisionEventSelection, phfCoincFilter;
	vector<float> *corrpt=0;
	double pthatweight=1.;
	double pthatbins[10] = {15,30,50,80,120,170,220,280,370,9999};
	double xsecs[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
	double pthatEntries[9] = {0, 444180, 485654, 484440, 479923, 448828, 258987, 234835, 50644};
	mc->SetBranchAddress("vz",&vz);
	mc->SetBranchAddress("hiBin",&hiBin);
	mc->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHEFilter);
	mc->SetBranchAddress("pcollisionEventSelection",&collisionEventSelection);
	mc->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
	mc->SetBranchAddress("pthat",&pthat);
	mc->SetBranchAddress("calo_corrpt",&corrpt);
	cout << "starting MC..." << endl;
	for(int ientry=0; ientry<mc->GetEntries(); ientry++){
		if(ientry && ientry%10000==0) cout << " At entry : "<< ientry << " of " << mc->GetEntries() << endl;
		mc->GetEntry(ientry);

		if(!HBHEFilter || !collisionEventSelection || !phfCoincFilter) continue;

		int ibin=0;
		while(pthat>pthatbins[ibin+1]) ibin++;
		pthatweight = (xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];

		int centBin = 0;
		while(hiBin/2. > centBinScheme[centBin]) centBin++;
		mcVz[centBin-1]->Fill(vz);
		mcVzPthat[centBin-1]->Fill(vz,pthatweight);
		if(corrpt->size()>0) if(corrpt->at(0)>120) mcVzPthatJetOnly[centBin-1]->Fill(vz,pthatweight);
		mcCent->Fill(hiBin);
	}

	dta->SetBranchAddress("vz",&vz);
	dta->SetBranchAddress("hiBin",&hiBin);
	dta->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHEFilter);
	dta->SetBranchAddress("pcollisionEventSelection",&collisionEventSelection);
	dta->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
	dta->SetBranchAddress("calo_corrpt",&corrpt);
	cout << "starting data..." << endl;
	for(int ientry=0; ientry<dta->GetEntries(); ientry++){
		if(ientry && ientry%10000==0) cout << " At entry : "<< ientry << " of " << dta->GetEntries() << endl;
		dta->GetEntry(ientry);

		if(!HBHEFilter || !collisionEventSelection || !phfCoincFilter) continue;

		int centBin = 0;
		while(hiBin/2. >= centBinScheme[centBin]) centBin++;
		//cout << "centBin: "<< centBin << endl;
		dataVz[centBin-1]->Fill(vz);
		if(corrpt->size()>0) if(corrpt->at(0)>120) dataVzPthatJetOnly[centBin-1]->Fill(vz);
		dataCent->Fill(hiBin);
	}

	TFile *fout = new TFile("VzCentReweights.root","recreate");
	fout->cd();
	TH1D *vzReweight[4], *vzReweightJetOnly[4], *vzReweightPthat[4];
	for(int i=0; i<4; i++){
		mcVz[i]->Scale(1./mcVz[i]->Integral());
		mcVz[i]->Write();
		dataVz[i]->Scale(1./dataVz[i]->Integral());
		dataVz[i]->Write();
		dataVzPthatJetOnly[i]->Scale(1./dataVzPthatJetOnly[i]->Integral());
		dataVzPthatJetOnly[i]->Write();
		mcVzPthat[i]->Scale(1./mcVzPthat[i]->Integral());
		mcVzPthat[i]->Write();
		mcVzPthatJetOnly[i]->Scale(1./mcVzPthatJetOnly[i]->Integral());
		mcVzPthatJetOnly[i]->Write();

		vzReweight[i] = (TH1D*)dataVz[i]->Clone(Form("vzReweight_cent%d",i));
		vzReweight[i]->Divide(mcVz[i]);
		vzReweight[i]->Fit("gaus");
		vzReweight[i]->Write();

		vzReweightJetOnly[i] = (TH1D*)dataVzPthatJetOnly[i]->Clone(Form("vzReweightJetOnly_cent%d",i));
		vzReweightJetOnly[i]->Divide(mcVzPthatJetOnly[i]);
		vzReweightJetOnly[i]->Fit("gaus");
		vzReweightJetOnly[i]->Write();

		vzReweightPthat[i] = (TH1D*)dataVz[i]->Clone(Form("vzReweightPthat_cent%d",i));
		vzReweightPthat[i]->Divide(mcVzPthat[i]);
		vzReweightPthat[i]->Fit("gaus");
		vzReweightPthat[i]->Write();

	}
	mcCent->Scale(1./mcCent->GetEntries());
	mcCent->Write();
	dataCent->Scale(1./dataCent->GetEntries());
	dataCent->Write();

	TH1D *centReweight = (TH1D*)dataCent->Clone("centReweight");
	centReweight->Divide(mcCent);
	centReweight->Fit("pol5");
	centReweight->Write();

}