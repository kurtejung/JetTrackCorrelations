
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TLatex.h"


void testReflection(){

	bool doMixEvtCorr = false;
	bool dov2BGSub = false;

	TFile *fvz = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/MCVZ.root");
	TFile *fvz2 = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/mcvz_allpthat.root");
	TFile *fvz3 = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/mcvz_pthat170.root");
	
	TFile *f1 = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_fixVz_Merged.root");
	//TFile *f2 = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_fixVz_sube0_Merged.root");
	//TH1D *shift = (TH1D*)f1->Get("GenJet_GenTrack_all_jets_etaCent0_Cent10_Pt100_Pt300");

	//TH1F *shift = (TH1F*)f1->Get("GenJet_GenTrack_TrkEtaCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2");
	//TH1F *shiftMix = (TH1F*)f2->Get("GenJet_GenTrack_TrkEtaCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2");

	//TH2D *inputMix = (TH2D*)f2->Get("GenJet_GenTrack_hJetTrackME_notrkcorrCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2");

	TH2D *inputMix = (TH2D*)f1->Get("GenJet_GenTrack_hJetTrackME_notrkcorrCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2")->Clone("inputMix");
	TH2D *input = (TH2D*)f1->Get("GenJet_GenTrack_hJetTrackSignalBackground_notrkcorrCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2")->Clone("input");
	
	TH1D *mixCorr = (TH1D*)inputMix->ProjectionX("mixCorr");
	mixCorr->Scale(1./mixCorr->GetBinContent(mixCorr->FindBin(0)));
	TH2D *mixCorrProj = (TH2D*)inputMix->Clone("mixCorrProj");
	for(int ixbin=1; ixbin<=mixCorrProj->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=mixCorrProj->GetNbinsY(); iybin++){
			mixCorrProj->SetBinContent(ixbin,iybin, mixCorr->GetBinContent(ixbin));
			mixCorrProj->SetBinError(ixbin,iybin, mixCorr->GetBinError(ixbin)*sqrt(mixCorr->GetNbinsX())); //re-expand errors to match repeated binning
		}
	}
	if(doMixEvtCorr){
	 inputMix->Divide(mixCorrProj);
	 input->Divide(mixCorrProj);
	}

	int lowBin = input->GetXaxis()->FindBin(-3.5);
	int highBin = input->GetXaxis()->FindBin(-2.);
	int totBins = highBin-lowBin;
	TH1D *phiBG = (TH1D*)input->ProjectionY("phiBG",lowBin,highBin,"e");
	lowBin = input->GetXaxis()->FindBin(2.);
	highBin = input->GetXaxis()->FindBin(3.5);
	totBins+=(highBin-lowBin);
	TH1D *phiBG2 = (TH1D*)input->ProjectionY("phiBG2",lowBin,highBin,"e");
	phiBG->Add(phiBG2);
	phiBG->Scale(1./(double)(totBins+2));
	TH2D *v2Sub = (TH2D*)input->Clone("v2Sub");
	for(int ixbin=1; ixbin<=mixCorrProj->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=mixCorrProj->GetNbinsY(); iybin++){
			v2Sub->SetBinContent(ixbin,iybin, phiBG->GetBinContent(iybin));
			v2Sub->SetBinError(ixbin,iybin, phiBG->GetBinError(iybin)*sqrt(phiBG->GetNbinsX())); //re-expand errors to match repeated binning
		}
	}

	if(dov2BGSub){
		input->Add(v2Sub,-1);
	}


	TH1D *shift = (TH1D*)input->ProjectionX("shift",1,input->GetNbinsX(),"e");
	TH1D *shiftMix = (TH1D*)inputMix->ProjectionX("shiftMix",1,inputMix->GetNbinsX(),"e");
	int lowbinNear = input->GetYaxis()->FindBin(-0.785);
	int highbinNear = input->GetYaxis()->FindBin(0.785);
	TH1D *shiftNear = (TH1D*)input->ProjectionX("shiftNear",lowbinNear,highbinNear,"e");
	TH1D *shiftNearMix = (TH1D*)inputMix->ProjectionX("shiftNearMix",lowbinNear,highbinNear,"e");
	int lowbin = input->GetYaxis()->FindBin(1.378);
	int highbin = input->GetYaxis()->FindBin(1.768);
	TH1D *shiftMid = (TH1D*)input->ProjectionX("shiftMid",lowbin,highbin,"e");
	TH1D *shiftMidMix = (TH1D*)inputMix->ProjectionX("shiftMidMix",lowbin,highbin,"e");
	lowbin = input->GetYaxis()->FindBin(2.356);
	highbin = input->GetYaxis()->FindBin(3.93);
	TH1D *shiftFar = (TH1D*)input->ProjectionX("shiftFar",lowbin,highbin,"e");
	TH1D *shiftFarMix = (TH1D*)inputMix->ProjectionX("shiftFarMix",lowbin,highbin,"e");
	
	TCanvas *c2 = new TCanvas("c2","",1200,600);
	c2->Divide(2,1);
	c2->cd(1);
	TH1D *forDraw = (TH1D*)input->ProjectionX("forDraw",lowbinNear,highbinNear,"e");
	//forDraw->Rebin(4);
	forDraw->Draw();
	//input->Draw("colz");
	c2->cd(2);
	TH1D *forDrawMix = (TH1D*)inputMix->ProjectionX("forDrawMix",lowbinNear,highbinNear,"e");
	//forDrawMix->Rebin(4);
	forDrawMix->Draw("");
	//inputMix->Draw("colz");

	int rebinNo = 2;

	shift->Rebin(rebinNo);
	shiftMix->Rebin(rebinNo);
	shiftNear->Rebin(rebinNo);
	shiftNearMix->Rebin(rebinNo);
	shiftMid->Rebin(rebinNo);
	shiftMidMix->Rebin(rebinNo);
	shiftFar->Rebin(rebinNo);
	shiftFarMix->Rebin(rebinNo);

	TH1D *refl = (TH1D*)shift->Clone("refl");
	TH1D *reflMix = (TH1D*)shiftMix->Clone("reflMix");

	TH1D *reflNear = (TH1D*)shiftNear->Clone("reflNear");
	TH1D *reflNearMix = (TH1D*)shiftNearMix->Clone("reflNearMix");

	TH1D *reflMid = (TH1D*)shiftMid->Clone("reflMid");
	TH1D *reflMidMix = (TH1D*)shiftMidMix->Clone("reflMidMix");

	TH1D *reflFar = (TH1D*)shiftFar->Clone("reflFar");
	TH1D *reflFarMix = (TH1D*)shiftFarMix->Clone("reflFarMix");

	TH1F *centVz = (TH1F*)fvz2->Get("hvz80");
	centVz->Rebin(2);
	TH1F *centVz2 = (TH1F*)fvz2->Get("hvz120");
	centVz2->Rebin(2);
	TH1F *centVz3 = (TH1F*)fvz3->Get("hvz170");
	centVz3->Rebin(5);
	TH1F *centVz4 = (TH1F*)fvz2->Get("hvz220");
	centVz4->Rebin(2);
	TH1F *reflVz = (TH1F*)centVz->Clone("reflVz");
	TH1F *reflVz2 = (TH1F*)centVz2->Clone("reflVz2");
	TH1F *reflVz3 = (TH1F*)centVz3->Clone("reflVz3");
	TH1F *reflVz4 = (TH1F*)centVz4->Clone("reflVz4");

	for(int i=1; i<=centVz->GetNbinsX(); i++){
		int reflBin = centVz->FindBin(-1*(centVz->GetBinLowEdge(i)+centVz->GetBinWidth(i)/2.));
		if(centVz->GetBinContent(reflBin)>0 && centVz->GetBinContent(i)>0){
			reflVz->SetBinContent(i, centVz->GetBinContent(i)/centVz->GetBinContent(reflBin));
			reflVz->SetBinError(i, sqrt(pow(centVz->GetBinError(i)/centVz->GetBinContent(i),2)+pow(centVz->GetBinError(reflBin)/centVz->GetBinContent(reflBin),2)));

			reflVz2->SetBinContent(i, centVz2->GetBinContent(i)/centVz2->GetBinContent(reflBin));
			reflVz2->SetBinError(i, sqrt(pow(centVz2->GetBinError(i)/centVz2->GetBinContent(i),2)+pow(centVz2->GetBinError(reflBin)/centVz2->GetBinContent(reflBin),2)));

			reflVz4->SetBinContent(i, centVz4->GetBinContent(i)/centVz4->GetBinContent(reflBin));
			reflVz4->SetBinError(i, sqrt(pow(centVz4->GetBinError(i)/centVz4->GetBinContent(i),2)+pow(centVz4->GetBinError(reflBin)/centVz4->GetBinContent(reflBin),2)));
		}
	}
	for(int i=1; i<=centVz3->GetNbinsX(); i++){
		int reflBin = centVz3->FindBin(-1*(centVz3->GetBinLowEdge(i)+centVz3->GetBinWidth(i)/2.));

		if(centVz3->GetBinContent(reflBin)>0 && centVz3->GetBinContent(i)>0){
			reflVz3->SetBinContent(i, centVz3->GetBinContent(i)/centVz3->GetBinContent(reflBin));
			reflVz3->SetBinError(i, sqrt(pow(centVz3->GetBinError(i)/centVz3->GetBinContent(i),2)+pow(centVz3->GetBinError(reflBin)/centVz3->GetBinContent(reflBin),2)));
		}
	}

	for(int i=1; i<=shift->GetNbinsX(); i++){

		int reflBin = shift->FindBin(-1*(shift->GetBinLowEdge(i)+shift->GetBinWidth(i)/2.));

		if(shift->GetBinContent(reflBin)>0 && shift->GetBinContent(i)>0){
			refl->SetBinContent(i, shift->GetBinContent(i)/shift->GetBinContent(reflBin));
			refl->SetBinError(i, sqrt(pow(shift->GetBinError(i)/shift->GetBinContent(i),2)+pow(shift->GetBinError(reflBin)/shift->GetBinContent(reflBin),2)));

			reflMix->SetBinContent(i, shiftMix->GetBinContent(i)/shiftMix->GetBinContent(reflBin));
			reflMix->SetBinError(i, sqrt(pow(shiftMix->GetBinError(i)/shiftMix->GetBinContent(i),2)+pow(shiftMix->GetBinError(reflBin)/shiftMix->GetBinContent(reflBin),2)));

			reflNear->SetBinContent(i, shiftNear->GetBinContent(i)/shiftNear->GetBinContent(reflBin));
			reflNear->SetBinError(i, sqrt(pow(shiftNear->GetBinError(i)/shiftNear->GetBinContent(i),2)+pow(shiftNear->GetBinError(reflBin)/shiftNear->GetBinContent(reflBin),2)));

			reflNearMix->SetBinContent(i, shiftNearMix->GetBinContent(i)/shiftNearMix->GetBinContent(reflBin));
			reflNearMix->SetBinError(i, sqrt(pow(shiftNearMix->GetBinError(i)/shiftNearMix->GetBinContent(i),2)+pow(shiftNearMix->GetBinError(reflBin)/shiftNearMix->GetBinContent(reflBin),2)));

			reflMid->SetBinContent(i, shiftMid->GetBinContent(i)/shiftMid->GetBinContent(reflBin));
			reflMid->SetBinError(i, sqrt(pow(shiftMid->GetBinError(i)/shiftMid->GetBinContent(i),2)+pow(shiftMid->GetBinError(reflBin)/shiftMid->GetBinContent(reflBin),2)));

			reflMidMix->SetBinContent(i, shiftMidMix->GetBinContent(i)/shiftMidMix->GetBinContent(reflBin));
			reflMidMix->SetBinError(i, sqrt(pow(shiftMidMix->GetBinError(i)/shiftMidMix->GetBinContent(i),2)+pow(shiftMidMix->GetBinError(reflBin)/shiftMidMix->GetBinContent(reflBin),2)));

			reflFar->SetBinContent(i, shiftFar->GetBinContent(i)/shiftFar->GetBinContent(reflBin));
			reflFar->SetBinError(i, sqrt(pow(shiftFar->GetBinError(i)/shiftFar->GetBinContent(i),2)+pow(shiftFar->GetBinError(reflBin)/shiftFar->GetBinContent(reflBin),2)));

			reflFarMix->SetBinContent(i, shiftFarMix->GetBinContent(i)/shiftFarMix->GetBinContent(reflBin));
			reflFarMix->SetBinError(i, sqrt(pow(shiftFarMix->GetBinError(i)/shiftFarMix->GetBinContent(i),2)+pow(shiftFarMix->GetBinError(reflBin)/shiftFarMix->GetBinContent(reflBin),2)));


		}
		
	}

	TCanvas *c0 = new TCanvas("c0","",600,600);
	c0->Divide(2,2);
	c0->cd(1);
	reflVz->Draw();
	c0->cd(2);
	reflVz2->Draw();
	c0->cd(3);
	reflVz3->Draw();
	c0->cd(4);
	reflVz4->Draw();

	string titles[4] = {"Inclusive #Delta#phi","-#pi/4 < #Delta#phi < #pi/4","3#pi/8 < #Delta#phi < 5#pi/8","6#pi/8 < #Delta#phi < 10#pi/8"};
	TLatex *l1[4];

	TCanvas *c1 = new TCanvas("c1","",1200,1200);
	//c1->Divide(2,2);
	c1->cd(1);
	refl->SetYTitle("Left/Right Asymmetry");
	refl->SetXTitle("Track Eta");
	refl->SetLineColor(1);
	refl->GetYaxis()->SetRangeUser(0.85,1.15);
	refl->Draw();
	reflMix->SetLineColor(1);
	reflMix->SetMarkerStyle(21);
	reflMix->Draw("same");
	l1[0] = new TLatex(0,1.035,titles[0].c_str());
	l1[0]->Draw("same");

	/*c1->cd(2);
	reflNear->SetLineColor(2);
	reflNear->GetYaxis()->SetRangeUser(0.85,1.15);
	reflNear->Draw();
	reflNearMix->SetMarkerStyle(21);
	reflNearMix->SetMarkerColor(2);
	reflNearMix->Draw("same");
	l1[1] = new TLatex(0,1.035,titles[1].c_str());
	l1[1]->Draw("same");
	TLatex *linfo = new TLatex(-4,1.01,"50-100% Cent, 1<TrkPt<2");
	linfo->Draw("same");

	c1->cd(3);
	reflMid->SetLineColor(4);
	reflMid->GetYaxis()->SetRangeUser(0.85,1.15);
	reflMid->Draw("");
	reflMidMix->SetMarkerStyle(21);
	reflMidMix->SetMarkerColor(4);
	reflMidMix->Draw("same");
	l1[2] = new TLatex(0,1.035,titles[2].c_str());
	l1[2]->Draw("same");

	c1->cd(4);
	reflFar->SetLineColor(8);
	reflFar->GetYaxis()->SetRangeUser(0.85,1.15);
	reflFar->Draw("");
	reflFarMix->SetMarkerStyle(21);
	reflFarMix->SetMarkerColor(8);
	reflFarMix->Draw("same");
	l1[3] = new TLatex(0,1.035,titles[3].c_str());
	l1[3]->Draw("same");*/

}