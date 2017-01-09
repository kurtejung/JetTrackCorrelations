#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include <iostream>

using namespace std;

void factorizeCorrs(){
	
	bool doSpillover = true;
	
	const double xTrkBinDouble[10] = {0.5,0.7,1.,2.,3.,4.,8.,12.,16.,20.};
	const string xCentBins[5] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
	
	TFile *f1 = new TFile("JFFcorrs_sube0_newJFFs.root");
	TFile *fsum = new TFile("JFFcorrs_fullsube_newJFFs.root");
	
	TFile *fout;
	if(doSpillover){
		fout = new TFile("JFFcorrs_spilloverFromsube0_newJFFs.root","recreate");
		fout->cd();
	}
	TH2D *hnum[10][4], *hden[10][4];
	TH1D *fp1[10][4];
	TLatex *labels[4];
	TCanvas *cc = new TCanvas("cc","",2000,1600);
	cc->Divide(4,9);
	
	TH1D *ptClosure[4];
	
	for(int j=0; j<4; j++){
		
		ptClosure[j] = new TH1D(Form("ptClosure_%d",j),"",9,xTrkBinDouble);
		
		for(int i=0; i<10; i++){
			
			hnum[i][j] = (TH2D*)fsum->Get(Form("JFFcorrs_cent%d_pt%d",j,i))->Clone(Form("JFFcorrs_cent%d_pt%d",j,i));
			hden[i][j] = (TH2D*)f1->Get(Form("JFFcorrs_cent%d_pt%d",j,i))->Clone(Form("den_cent%d_pt%d",j,i));
			if(doSpillover){
				hnum[i][j]->Add(hden[i][j],-1);
				hnum[i][j]->Write();
			}
			//else hnum[i][j] = hden[i][j];
			//hnum[i][j] = hden[i][j];

			if(i>0 && i<9){
				cc->cd(i*4+(3-j)+1-4);
				fp1[i][j] = (TH1D*)hnum[i][j]->ProjectionX();
				fp1[i][j]->GetXaxis()->SetRangeUser(-1.5,1.5);
				fp1[i][j]->GetYaxis()->SetNdivisions(505);
				fp1[i][j]->GetYaxis()->SetLabelSize(0.15);
				fp1[i][j]->GetXaxis()->SetLabelSize(0.15);
				fp1[i][j]->Rebin(2);
				fp1[i][j]->Draw();
			
				double integrErr =0.;
				double integr = fp1[i][j]->IntegralAndError(1, fp1[i][j]->GetNbinsX(), integrErr);
				
				ptClosure[j]->SetBinContent(i+1, integr/ptClosure[j]->GetBinWidth(i+1));
				ptClosure[j]->SetBinError(i+1, integrErr/ptClosure[j]->GetBinWidth(i+1));
				
				cout << "cent bin " << j << " pt bin " << i << " integral: " << fp1[i][j]->Integral() << endl;
			}
			
		}
	}
	
	TCanvas *c2 = new TCanvas("c2","",2000,400);
	c2->Divide(4,1);
	for(int j=0; j<4; j++){
		c2->cd(4-j);
		ptClosure[j]->SetMaximum(3);
		ptClosure[j]->SetMinimum(-1);
		ptClosure[j]->SetXTitle("Track p_{T} (GeV/c)");
		if(!doSpillover) ptClosure[j]->SetYTitle("JFF Correction");
		else ptClosure[j]->SetYTitle("Spillover Correction");
		ptClosure[j]->Draw();
		labels[j] = new TLatex(3,ptClosure[j]->GetMaximum()*0.85,Form("%s-%s",xCentBins[j].c_str(),xCentBins[j+1].c_str()));
		labels[j]->SetTextSize(0.06);
		labels[j]->Draw("same");
	}
	
}