

#include "Rtypes.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH1D.h"
#include "TColor.h"
#include <iostream>

using namespace std;

void drawJetShapes(){
	
	const int nCentBins=4;
	const int trkPtBins=10;
	
	const string xTrkBins[trkPtBins+1] = {"TrkPt05", "TrkPt07", "TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
	const string xCentBins[nCentBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
	
	const double mean_pts[trkPtBins] = {0.55,0.844,1.35,2.35,3.37,5.07,9.72,13.8,17.9,41.};
	const double ptWidth[trkPtBins] = {0.2,0.3,1.,1.,1.,4.,4.,4.,4.,100.};
	
	const int colors[trkPtBins] = {kBlack, kBlue-9, kYellow-9, kOrange+1, kViolet-5, kGreen+3, kRed, kRed+1, kRed+2, kRed+3};
	const double newXbins[15] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.};
	
	TH1D *PbPb_data[nCentBins][trkPtBins];
	TH1D *pp_data[trkPtBins];
	TH1D *pp_MC[trkPtBins];
	THStack *PbPbStack[nCentBins];
	THStack *ppStack;
	THStack *ppMCStack;
	
	TH1D *ppSumHisto;
	TH1D *ppMCSumHisto;
	TH1D *PbPbSumHisto[nCentBins];
	TH1D *PbPbppRatios[nCentBins];
	
	TFile *fin = new TFile("Data_bgSubtractedCorrs_dR_fullHistos.root");
	TFile *finHMC = new TFile("HallieppMC_bgSubtractedCorrs_dR_fullHistos.root");
	TFile *finMC = new TFile("ppMC_bgSubtractedCorrs_dR_fullHistos.root");
	
	for(int i=0; i<nCentBins; i++){
		for(int j=0; j<trkPtBins; j++){
			string toGet = Form("JetShape2_Yield_BkgSub_pTweightedInclusive_%s_%s_%s_%s",xCentBins[i].c_str(),xCentBins[i+1].c_str(), xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
			PbPb_data[i][j] = (TH1D*)fin->Get(toGet.c_str())->Clone(toGet.c_str());
			PbPb_data[i][j] = (TH1D*)PbPb_data[i][j]->Rebin(14,PbPb_data[i][j]->GetName(), newXbins);
			PbPb_data[i][j]->SetFillColor(colors[j]);
			PbPb_data[i][j]->GetXaxis()->SetRangeUser(0,1);
			//PbPb_data[i][j]->Scale(1./ptWidth[j]);
		}
	}
	for(int j=0; j<trkPtBins; j++){
		string toGet = Form("JetShape2_Yield_BkgSub_pTweightedInclusive_pp_%s_%s",xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
		pp_data[j] = (TH1D*)fin->Get(toGet.c_str())->Clone(toGet.c_str());
		pp_data[j] = (TH1D*)pp_data[j]->Rebin(14,pp_data[j]->GetName(), newXbins);
		pp_data[j]->SetFillColor(colors[j]);
		pp_data[j]->GetXaxis()->SetRangeUser(0,1);
		//pp_data[j]->Scale(1./ptWidth[j]);
	}
	
	for(int j=0; j<trkPtBins; j++){
		string toGet = Form("JetShape2_Yield_BkgSub_pTweightedInclusive_ppMC_%s_%s",xTrkBins[j].c_str(),xTrkBins[j+1].c_str());
		pp_MC[j] = (TH1D*)finMC->Get(toGet.c_str())->Clone(toGet.c_str());
		pp_MC[j] = (TH1D*)pp_MC[j]->Rebin(14,pp_MC[j]->GetName(), newXbins);
		pp_MC[j]->SetFillColor(colors[j]);
		pp_MC[j]->GetXaxis()->SetRangeUser(0,1);
		//pp_data[j]->Scale(1./ptWidth[j]);
	}


	//Stack and sum them...
	for(int i=0; i<nCentBins; i++){
		PbPbStack[i] = new THStack(Form("PbPbStack_Cent%d",i),"");
		for(int j=1; j<trkPtBins; j++){
			for(int ibin=1; ibin<=PbPb_data[i][j]->GetNbinsX(); ibin++){
				PbPb_data[i][j]->SetBinContent(ibin, PbPb_data[i][j]->GetBinContent(ibin)/PbPb_data[i][j]->GetBinWidth(ibin));
				PbPb_data[i][j]->SetBinError(ibin, PbPb_data[i][j]->GetBinError(ibin)/PbPb_data[i][j]->GetBinWidth(ibin));
			}
			//PbPb_data[i][j]->Scale(3./PbPb_data[i][j]->Integral(1,PbPb_data[i][j]->GetNbinsX()));
			PbPbStack[i]->Add(PbPb_data[i][j],"hist");
		}
		
		PbPbSumHisto[i] = (TH1D*)PbPb_data[i][1]->Clone(Form("PbPbSumHisto_cent%d",i));
		for(int j=2; j<trkPtBins; j++){
			PbPbSumHisto[i]->Add(PbPb_data[i][j]);
		}
		
	}
	ppStack = new THStack("ppStack","");
	for(int j=1; j<trkPtBins; j++){
		for(int ibin=1; ibin<=pp_data[j]->GetNbinsX(); ibin++){
			pp_data[j]->SetBinContent(ibin, pp_data[j]->GetBinContent(ibin)/pp_data[j]->GetBinWidth(ibin));
			pp_data[j]->SetBinError(ibin, pp_data[j]->GetBinError(ibin)/pp_data[j]->GetBinWidth(ibin));
			
			pp_MC[j]->SetBinContent(ibin, pp_MC[j]->GetBinContent(ibin)/pp_MC[j]->GetBinWidth(ibin));
			pp_MC[j]->SetBinError(ibin, pp_MC[j]->GetBinError(ibin)/pp_MC[j]->GetBinWidth(ibin));
		}
		//pp_data[j]->Scale(3./pp_data[j]->Integral(1,pp_data[j]->GetNbinsX()));
		ppStack->Add(pp_data[j],"hist");
	}
	
	ppSumHisto = (TH1D*)pp_data[1]->Clone("ppSumHisto");
	ppMCSumHisto = (TH1D*)pp_MC[1]->Clone("ppMCSumHisto");
	for(int j=2; j<trkPtBins; j++){
		ppSumHisto->Add(pp_data[j]);
		ppMCSumHisto->Add(pp_MC[j]);
	}
	
	
	for(int i=0; i<nCentBins; i++){
		PbPbppRatios[i] = (TH1D*)PbPbSumHisto[i]->Clone(Form("PbPbppRatio_cent%d",i));
		cout << " cent bin " << i << " " << PbPbppRatios[i]->GetBinContent(1) << endl;
		PbPbppRatios[i]->Divide(ppSumHisto);
	}
	
	
	TCanvas *cc = new TCanvas("cc","",1400,600);
	cc->Divide(5,2,0,0);
	for(int i=1; i<6; i++){
		cc->cd(i);
		cc->GetPad(i)->SetLogy();
		if(i==1){
			//ppStack->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
			ppSumHisto->SetMaximum(3000);
			ppStack->SetMinimum(0.2);
			//ppStack->Draw("");
			//ppStack->GetHistogram()->GetYaxis()->SetTitle("Pythia Gen-Gen P(#Deltar)");
			//ppStack->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
			ppSumHisto->SetYTitle("P(#Deltar)");
			ppSumHisto->Draw();
			ppMCSumHisto->SetMarkerColor(kBlue-3);
			ppMCSumHisto->Draw("Same");
		} 
		else{
			//PbPbStack[5-i]->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
			PbPbStack[5-i]->SetMaximum(3000);
			PbPbStack[5-i]->SetMinimum(0.2);
			//PbPbStack[5-i]->Draw("");
			//PbPbStack[5-i]->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
		}
		//cc->GetPad(i)->Update();
	}
	cc->cd(6);
	TH1D *ppDataMCRatio = (TH1D*)ppSumHisto->Clone("ppDataMCRatio");
	ppDataMCRatio->Divide(ppMCSumHisto);
	ppDataMCRatio->GetYaxis()->SetTitle("Data / Kurt MC");
	ppDataMCRatio->Draw();
	
	for(int i=7; i<=10; i++){
		cc->cd(i);
		PbPbppRatios[10-i]->SetMaximum(3.5);
		PbPbppRatios[10-i]->SetMinimum(0);
		//PbPbppRatios[10-i]->Draw();
	}
	
	//TCanvas *cc2 = new TCanvas("cc2","",1400,600);
	//cc2->Divide(1,8);
	
	
}