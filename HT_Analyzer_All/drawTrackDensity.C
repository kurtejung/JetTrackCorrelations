#include "TFile.h"
#include "TH2D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TH1D.h"
#include <string>
#include <iostream>
#include "TCanvas.h"
#include "TLatex.h"

using namespace std;

void add2DBin(int ixbin, int iybin, int ixToAdd, int iyToAdd, TH2D *histo){
	histo->SetBinContent(ixbin,iybin,histo->GetBinContent(ixbin,iybin)+histo->GetBinContent(ixToAdd,iyToAdd));
	histo->SetBinError(ixbin,iybin,sqrt(pow(histo->GetBinError(ixbin,iybin),2) + pow(histo->GetBinError(ixToAdd,iyToAdd),2)));
}

double ReturnDivError(double num, double numErr, double den, double denErr){
	if(den==0 || num==0) return 0;
	double err2 = pow(num/den,2)*(pow(numErr/num,2)+pow(denErr/den,2));
	return sqrt(err2);
}

double ReturnSubError(double first, double firstErr, double sec, double secErr){
	double err2 = pow(first*firstErr,2)+pow(sec*secErr,2);
	return sqrt(err2);
}

TH2D* backgroundSubtract(TH2D *signal, TH2D *background, TH1F *nJets, bool doBGsub){
	
	double sidebandlow = 1.5;
	double sidebandhigh = 2.5;
	
	//Smooth background to be phi-independent
	TH1D *projEta = (TH1D*)background->ProjectionX("projEta");
	for(int ixbin=1; ixbin<=background->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=background->GetNbinsY(); iybin++){
			background->SetBinContent(ixbin,iybin, projEta->GetBinContent(ixbin));
			background->SetBinError(ixbin,iybin, projEta->GetBinError(ixbin)*sqrt(projEta->GetNbinsX()));
		}
	}
	background->Scale(1./projEta->GetBinContent(projEta->FindBin(0.)));
	
	TH1D *inclBG = (TH1D*)background->ProjectionX("inclBG",0, -1,"e");
	inclBG->Scale(1./inclBG->GetBinContent(inclBG->FindBin(0)));
	
	//Do mixed event correction
	for(int ixbin=1; ixbin<=signal->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=signal->GetNbinsY(); iybin++){
			if(inclBG->GetBinContent(ixbin)) signal->SetBinContent(ixbin,iybin, signal->GetBinContent(ixbin,iybin)/inclBG->GetBinContent(ixbin));
			else signal->SetBinContent(ixbin,iybin,0);
			if(inclBG->GetBinContent(ixbin)) signal->SetBinError(ixbin,iybin, ReturnDivError(signal->GetBinContent(ixbin,iybin), signal->GetBinError(ixbin,iybin), inclBG->GetBinContent(ixbin),inclBG->GetBinError(ixbin)));
		}
	}
	signal->Scale(1./nJets->Integral());
	//This may or may not be needed?
	signal->Scale(signal->GetNbinsX()*3./5.);
	
	//Get the dphi projection to remove v2, etc
	int lowBin = signal->GetXaxis()->FindBin(-1*sidebandhigh);
	int highBin = signal->GetXaxis()->FindBin(-1*sidebandlow);
	int totalBins = highBin-lowBin;
	TH1D *inclPhiBG = (TH1D*)signal->ProjectionY("inclPhiBG",lowBin, highBin,"e");
	
	lowBin = signal->GetXaxis()->FindBin(sidebandlow);
	highBin = signal->GetXaxis()->FindBin(sidebandhigh);
	totalBins+=(highBin-lowBin);
	TH1D *tmp1 = (TH1D*)signal->ProjectionY("tmp1",lowBin, highBin,"e");
	inclPhiBG->Add(tmp1);
	inclPhiBG->Scale(1./(double)(totalBins+2));
	
	lowBin = signal->GetYaxis()->FindBin(-5);
	highBin = signal->GetYaxis()->FindBin(5);
	totalBins = highBin-lowBin;
	TH1D *inclPhiFG = signal->ProjectionY("inclPhiFG",lowBin, highBin,"e");
	inclPhiFG->Scale(1./inclPhiFG->GetNbinsX());
	assert(signal->GetNbinsY() == inclPhiBG->GetNbinsX());
	
	if(doBGsub){
		for(int ixbin=1; ixbin<=signal->GetNbinsX(); ixbin++){
			for(int iybin=1; iybin<=signal->GetNbinsY(); iybin++){
				if(inclPhiBG->GetBinContent(iybin)) signal->SetBinContent(ixbin,iybin, signal->GetBinContent(ixbin,iybin) - inclPhiBG->GetBinContent(iybin));
				else signal->SetBinContent(ixbin,iybin,0);
				if(inclPhiBG->GetBinContent(iybin)) signal->SetBinError(ixbin,iybin, sqrt(inclPhiBG->GetNbinsX())*ReturnSubError(signal->GetBinContent(ixbin,iybin), signal->GetBinError(ixbin,iybin), inclPhiBG->GetBinContent(iybin), inclPhiBG->GetBinError(iybin)));
			}
		}
	}
	
	return signal;
}


void drawTrackDensity(string jet="Gen", string track="Reco", string jetDiv="Gen", string trackDiv="Gen", bool doBGsub=true){
	
	bool doJFFCorrs = false;
	bool applyJFFcorrs = false;
	
	string trkCorr = "_notrkcorr";
	string trkCorrDiv = "_notrkcorr";
	if(track=="Reco") trkCorr = "";
	if(trackDiv=="Reco")trkCorrDiv = "";
	if(jet=="Reco") applyJFFcorrs = true;
	//if(jet=="Reco" && jetDiv=="Reco") applyJFFcorrs = false;
		
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	
	const double xdrbins[15] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.};
	
	const int trackPtBins = 10;
	const int nCentBins = 4;
	
	const double xTrkBinDouble[trackPtBins+1] = {0.5,0.7,1.,2.,3.,4.,8.,12.,16.,20.,999.};
	const float mean_pts[trackPtBins] = {0.5,0.844,1.35,2.35,3.37,5.07,9.72,13.8,17.9,22.};
	const string xTrkBins[trackPtBins+1] = {"TrkPt0p5","TrkPt0p7", "TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt999" };
	const string xCentBins[nCentBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
	
	TFile *fJFFs = new TFile("JFFcorrs_static.root");
	TFile *fGenGen = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_%s%sReduced_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenReco = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_%s%sReduced_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	
	TH2D *GenGenSignals[trackPtBins][nCentBins];
	TH2D *GenRecoSignals[trackPtBins][nCentBins];
	TH2D *GenGenBG[trackPtBins][nCentBins];
	TH2D *GenRecoBG[trackPtBins][nCentBins];
	
	TH2D *GenGenSignalsSub[trackPtBins][nCentBins];
	TH2D *GenRecoSignalsSub[trackPtBins][nCentBins];
	
	TH1D *GenGenDRDistr[trackPtBins][nCentBins];
	TH1D *GenRecoDRDistr[trackPtBins][nCentBins];
	
	TH1D *GenGenDRDistrNorm[trackPtBins][nCentBins];
	TH1D *GenRecoDRDistrNorm[trackPtBins][nCentBins];
	
	TH1D *GenGenEtaProj[trackPtBins][nCentBins];
	TH1D *GenRecoEtaProj[trackPtBins][nCentBins];
	
	TH2D *JFFcorrs[trackPtBins];
	TH2D *JFFcorrsToApply[trackPtBins];
	
	THStack *drGenGenDensity[nCentBins];
	THStack *drGenRecoDensity[nCentBins];
	
	TH1F *nJetsGen[nCentBins];
	TH1F *nJetsReco[nCentBins];
	
	for(int i=0; i<trackPtBins; i++){
		
		string toGet = Form("JFFcorrs_pt%d",i);
		JFFcorrsToApply[i] = (TH2D*)fJFFs->Get(toGet.c_str())->Clone(toGet.c_str());
		//JFFcorrsToApply[i]->Scale(0.5);
		
		for(int j=0; j<nCentBins; j++){
			toGet = Form("%sJet_%sTrack_hJetTrackSignalBackground%s%s_%s_Pt100_Pt300_%s_%s",jet.c_str(),track.c_str(),trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			cout << "getting " << toGet << endl;
			GenGenSignals[i][j] = (TH2D*)fGenGen->Get(toGet.c_str())->Clone(toGet.c_str());
			
			toGet = Form("%sJet_%sTrack_hJetTrackSignalBackground%s%s_%s_Pt100_Pt300_%s_%s",jetDiv.c_str(),trackDiv.c_str(),trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			GenRecoSignals[i][j] = (TH2D*)fGenReco->Get(toGet.c_str())->Clone(toGet.c_str());
			
			toGet = Form("%sJet_%sTrack_hJetTrackME%s%s_%s_Pt100_Pt300_%s_%s",jet.c_str(),track.c_str(),trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			GenGenBG[i][j] = (TH2D*)fGenGen->Get(toGet.c_str())->Clone(toGet.c_str());
			
			toGet = Form("%sJet_%sTrack_hJetTrackME%s%s_%s_Pt100_Pt300_%s_%s",jetDiv.c_str(),trackDiv.c_str(),trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			GenRecoBG[i][j] = (TH2D*)fGenReco->Get(toGet.c_str())->Clone(toGet.c_str());
			
			GenGenDRDistr[i][j] = new TH1D(Form("GenGenDR_%d_%d",i,j),"",14,xdrbins);
			GenRecoDRDistr[i][j] = new TH1D(Form("GenRecoDR_%d_%d",i,j),"",14,xdrbins);
			GenGenDRDistrNorm[i][j] = new TH1D(Form("GenGenDRNorm_%d_%d",i,j),"",14,xdrbins);
			GenRecoDRDistrNorm[i][j] = new TH1D(Form("GenRecoDRNorm_%d_%d",i,j),"",14,xdrbins);
		}
	}
	
	for(int j=0; j<nCentBins; j++){
		string toGet = Form("%sJet_%sTrack_all_jets_corrpT%s_%s_Pt100_Pt300",jet.c_str(),track.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str());
		nJetsGen[j] = (TH1F*)fGenGen->Get(toGet.c_str())->Clone(toGet.c_str());
		
		toGet = Form("%sJet_%sTrack_all_jets_corrpT%s_%s_Pt100_Pt300",jetDiv.c_str(),trackDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str());
		nJetsReco[j] = (TH1F*)fGenReco->Get(toGet.c_str())->Clone(toGet.c_str());
		
		drGenGenDensity[j] = new THStack(Form("drGenGenDensity_%d",j),"");
		drGenRecoDensity[j] = new THStack(Form("drGenRecoDensity_%d",j),"");
	}
	
	for(int i=0; i<trackPtBins; i++){
		for(int j=0; j<nCentBins; j++){
			GenGenSignalsSub[i][j] = backgroundSubtract(GenGenSignals[i][j], GenGenBG[i][j], nJetsGen[j], doBGsub);
			GenRecoSignalsSub[i][j] = backgroundSubtract(GenRecoSignals[i][j], GenRecoBG[i][j], nJetsReco[j], doBGsub);
			if(applyJFFcorrs){
				if(jet=="Reco") GenGenSignalsSub[i][j]->Add(JFFcorrsToApply[i]);
				if(jetDiv=="Reco") GenRecoSignalsSub[i][j]->Add(JFFcorrsToApply[i]);				
			}
			
			//cout << "track pt" << i << " cent " << j << " gen counts: "<< GenGenSignalsSub[i][j]->GetEntries() << endl;
			//cout << "track pt" << i << " cent " << j << " counts: "<< GenRecoSignalsSub[i][j]->GetEntries() << endl;
			
			int lowBin = GenGenSignalsSub[i][j]->GetYaxis()->FindBin(-1);
			int highBin = GenGenSignalsSub[i][j]->GetYaxis()->FindBin(1);
			GenGenEtaProj[i][j] = GenGenSignalsSub[i][j]->ProjectionX(Form("GenGenEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			GenGenEtaProj[i][j]->Scale(GenGenSignalsSub[i][j]->GetYaxis()->GetBinWidth(2));
			GenGenEtaProj[i][j]->Scale(2.);
			GenRecoEtaProj[i][j] = GenRecoSignalsSub[i][j]->ProjectionX(Form("GenRecoEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			GenRecoEtaProj[i][j]->Scale(GenRecoSignalsSub[i][j]->GetYaxis()->GetBinWidth(2));
			GenRecoEtaProj[i][j]->Scale(2.);
			//GenRecoEtaProj[i][j]->Add(GenGenEtaProj[i][j],-1);
			GenRecoEtaProj[i][j]->Divide(GenGenEtaProj[i][j]);
			
			if(j==nCentBins-1 && doJFFCorrs){
				JFFcorrs[i] = (TH2D*)GenRecoSignalsSub[i][j]->Clone(Form("JFFcorrs_pt%d",i));
				JFFcorrs[i]->Add(GenGenSignalsSub[i][j], -1);
				
				//now symmetrize the corrections
				for(int ixbin=1; ixbin<=JFFcorrs[i]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<=JFFcorrs[i]->GetNbinsY(); iybin++){
						double xcent = JFFcorrs[i]->GetXaxis()->GetBinCenter(ixbin);
						double ycent = JFFcorrs[i]->GetYaxis()->GetBinCenter(iybin);
						int negX = JFFcorrs[i]->GetXaxis()->FindBin(-1*xcent);
						int negY = JFFcorrs[i]->GetYaxis()->FindBin(-1*ycent);
						
						double dr = sqrt(pow(xcent,2)+pow(ycent,2));
						if(dr>0.5) JFFcorrs[i]->SetBinContent(ixbin,iybin,0);
						else{
							if(xcent<0 || ycent<0) continue;
							add2DBin(ixbin, iybin, negX, iybin, JFFcorrs[i]);
							add2DBin(ixbin, iybin, ixbin, negY, JFFcorrs[i]);
							add2DBin(ixbin, iybin, negX, negY, JFFcorrs[i]);
							
							JFFcorrs[i]->SetBinContent(negX,iybin,JFFcorrs[i]->GetBinContent(ixbin,iybin));
							JFFcorrs[i]->SetBinError(negX,iybin,JFFcorrs[i]->GetBinError(ixbin,iybin));	
							
							JFFcorrs[i]->SetBinContent(ixbin,negY,JFFcorrs[i]->GetBinContent(ixbin,iybin));
							JFFcorrs[i]->SetBinError(ixbin,negY,JFFcorrs[i]->GetBinError(ixbin,iybin));	
							
							JFFcorrs[i]->SetBinContent(negX,negY,JFFcorrs[i]->GetBinContent(ixbin,iybin));
							JFFcorrs[i]->SetBinError(negX,negY,JFFcorrs[i]->GetBinError(ixbin,iybin));						
						}
					}
					
				}
				//take average of symmetry
				JFFcorrs[i]->Scale(1./4.);
			}
			
			for(int ixbin=1; ixbin<=GenGenSignalsSub[i][j]->GetNbinsX(); ixbin++){
				for(int iybin=1; iybin<=GenGenSignalsSub[i][j]->GetNbinsY(); iybin++){
					double xcent = GenGenSignalsSub[i][j]->GetXaxis()->GetBinCenter(ixbin);
					double ycent = GenGenSignalsSub[i][j]->GetYaxis()->GetBinCenter(iybin);
					double content = GenGenSignalsSub[i][j]->GetBinContent(GenGenSignalsSub[i][j]->GetBin(ixbin,iybin));
					//cout << "gen gen bin " << ixbin << ", " << iybin << " content " << content << endl;
					double dr = sqrt(pow(xcent,2)+pow(ycent,2));
					GenGenDRDistr[i][j]->Fill(dr,content*mean_pts[i]);
					GenGenDRDistrNorm[i][j]->Fill(dr);
					
					if(j<8) drGenGenDensity[j]->Add(GenGenDRDistr[i][j]);
										
					xcent = GenRecoSignalsSub[i][j]->GetXaxis()->GetBinCenter(ixbin);
					ycent = GenRecoSignalsSub[i][j]->GetYaxis()->GetBinCenter(iybin);
					content = GenRecoSignalsSub[i][j]->GetBinContent(GenRecoSignalsSub[i][j]->GetBin(ixbin,iybin));
					//cout << "gen reco bin " << ixbin << ", " << iybin << " content " << content << endl;
					dr = sqrt(pow(xcent,2)+pow(ycent,2));
					GenRecoDRDistr[i][j]->Fill(dr,content*mean_pts[i]);
					GenRecoDRDistrNorm[i][j]->Fill(dr);
					
					if(j<8) drGenRecoDensity[j]->Add(GenRecoDRDistr[i][j]);
					
				}
			}
			GenGenDRDistr[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			GenRecoDRDistr[i][j]->Divide(GenRecoDRDistrNorm[i][j]);
		}
	}
	
	TCanvas *ceta = new TCanvas("ceta","",1200,1600);
	TLatex *l2[trackPtBins][nCentBins];
	ceta->Divide(nCentBins,trackPtBins-2);
	for(int j=0; j<nCentBins; j++){
		for(int i=1; i<trackPtBins-1; i++){
			ceta->cd(i*nCentBins+j+1-4);
			GenRecoEtaProj[i][j]->GetYaxis()->SetNdivisions(505);
			GenRecoEtaProj[i][j]->GetYaxis()->SetLabelSize(0.15);
			GenRecoEtaProj[i][j]->GetXaxis()->SetLabelSize(0.15);
			GenRecoEtaProj[i][j]->GetXaxis()->SetRangeUser(-0.7,0.7);
			GenRecoEtaProj[i][j]->Draw();
			//GenGenEtaProj[i][j]->SetLineColor(2);
			//GenGenEtaProj[i][j]->Draw("Same");
			l2[i][j] = new TLatex(-0.5,GenRecoEtaProj[i][j]->GetMaximum()*0.85,Form("%s-%s, %g<pT<%g",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBinDouble[i],xTrkBinDouble[i+1]));
			l2[i][j]->SetTextSize(0.15);
			l2[i][j]->Draw("same");
		}
	}
	
	TCanvas *cc = new TCanvas("cc","",1200,800);
	cc->Divide(1,3);
	cc->cd(1);
	GenRecoSignalsSub[4][0]->Draw("colz");
	cc->cd(2);
	TH2D *tmp = (TH2D*)GenRecoSignalsSub[4][0]->Clone("tmp");
	tmp->Add(JFFcorrsToApply[4],-1);
	tmp->Draw("colz");
	cc->cd(3);
	JFFcorrsToApply[4]->Draw("colz");
	
	TLatex *l1[trackPtBins][nCentBins];
	TCanvas *cprint = new TCanvas("cprint","",1200,1600);
	cprint->Divide(nCentBins, trackPtBins-2);
	for(int j=0; j<nCentBins; j++){
		for(int i=1; i<trackPtBins-1; i++){
			cprint->cd(i*nCentBins+(3-j)+1-4);
			
			GenRecoDRDistr[i][j]->SetMarkerStyle(25);
			
			for(int k=1; k<=GenGenDRDistr[i][j]->GetNbinsX(); k++){
				//GenGenDRDistr[i][j]->SetBinContent(k,GenGenDRDistr[i][j]->GetBinContent(k)/GenGenDRDistr[i][j]->GetBinWidth(k));
				//GenGenDRDistr[i][j]->SetBinError(k,GenGenDRDistr[i][j]->GetBinError(k)/GenGenDRDistr[i][j]->GetBinWidth(k));
			
				//GenRecoDRDistr[i][j]->SetBinContent(k,GenRecoDRDistr[i][j]->GetBinContent(k)/GenRecoDRDistr[i][j]->GetBinWidth(k));
				//GenRecoDRDistr[i][j]->SetBinError(k,GenRecoDRDistr[i][j]->GetBinError(k)/GenRecoDRDistr[i][j]->GetBinWidth(k));
			}
			
			//GenGenDRDistr[i][j]->Draw();
			//GenRecoDRDistr[i][j]->Draw("same");
			GenGenDRDistr[i][j]->Divide(GenGenDRDistr[i][j],GenRecoDRDistr[i][j],1,1,"B");
			//GenGenDRDistr[i][j]->Divide(GenRecoDRDistr[i][j],GenGenDRDistr[i][j],1,1,"B");
			GenGenDRDistr[i][j]->SetMaximum(1.5);
			GenGenDRDistr[i][j]->SetMinimum(0.5);
			GenGenDRDistr[i][j]->GetYaxis()->SetNdivisions(505);
			GenGenDRDistr[i][j]->GetYaxis()->SetLabelSize(0.15);
			GenGenDRDistr[i][j]->GetXaxis()->SetLabelSize(0.15);
			//GenGenDRDistr[i][j]->SetYTitle("RecoReco / RecoGen #rho(#DeltaR)");
			GenGenDRDistr[i][j]->Draw("");
			l1[i][j] = new TLatex(0.2,GenGenDRDistr[i][j]->GetMaximum()*0.85,Form("%s-%s, %g<pT<%g",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBinDouble[i],xTrkBinDouble[i+1]));
			l1[i][j]->SetTextSize(0.15);
			l1[i][j]->Draw("same");
		}
	}
	string jfc = "jffCorrApplied";
	if(!applyJFFcorrs) jfc = "";
	cprint->SaveAs(Form("%s%s_div_%s%s_bgSub%d_%s_test.pdf",jet.c_str(), track.c_str(), jetDiv.c_str(), trackDiv.c_str(), doBGsub, jfc.c_str()));
	
	if(doJFFCorrs){
		TFile *fout = new TFile("JFFcorrs.root","recreate");
		fout->cd();
		for(int i=0; i<trackPtBins; i++){
			JFFcorrs[i]->Write();
		}
		fout->Close();
	}
}