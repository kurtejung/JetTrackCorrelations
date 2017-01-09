#include "TFile.h"
#include "TH2D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TH1D.h"
#include <string>
#include <iostream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraphErrors.h"

using namespace std;

void add2DBin(int ixbin, int iybin, int toAdd, TH2D *histo){
	histo->SetBinContent(ixbin,iybin,histo->GetBinContent(ixbin,iybin)+histo->GetBinContent(toAdd));
	histo->SetBinError(ixbin,iybin,sqrt(pow(histo->GetBinError(ixbin,iybin),2) + pow(histo->GetBinError(toAdd),2)));
}

double ReturnDivError(double num, double numErr, double den, double denErr){
	if(den==0 || num==0) return 0;
	double err2 = pow(num/den,2)*((pow(numErr/num,2)+pow(denErr/den,2)));
	return sqrt(err2);
}

double ReturnSubError(double firstErr, double secErr){
	double err2 = pow(firstErr,2)+pow(secErr,2);
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
	delete tmp1;
	
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
				if(inclPhiBG->GetBinContent(iybin)) signal->SetBinError(ixbin,iybin, sqrt(inclPhiBG->GetNbinsX())*ReturnSubError(signal->GetBinError(ixbin,iybin), inclPhiBG->GetBinError(iybin)));
			}
		}
	}
	
	delete inclPhiFG;
	delete inclPhiBG;
	delete inclBG;
	delete projEta;
	
	return signal;
}

void drawTrackDensity(string jet="Gen", string track="Gen", string jetDiv="Reco", string trackDiv="Gen", bool doBGsub=true){
	
	bool doJFFCorrs = true;
	bool applyJFFcorrs = false;
	bool writeOutCorrelations = false;
	bool drawUncertainties = false;
	
	string trkCorr = "_notrkcorr";
	string trkCorrDiv = "_notrkcorr";
	if(track=="Reco") trkCorr = "";
	if(trackDiv=="Reco")trkCorrDiv = "";
	if(jet=="Reco" || jetDiv=="Reco") applyJFFcorrs = true;
	if(jet=="Reco" && jetDiv=="Reco") applyJFFcorrs = false;
		
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	
	const double xdrbins[15] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.};
	
	const int trackPtBins = 10;
	const int nCentBins = 4;
	
	const double xTrkBinDouble[trackPtBins+1] = {0.5,0.7,1.,2.,3.,4.,8.,12.,16.,20.,999.};
	const float mean_pts[trackPtBins] = {0.5,0.844,1.35,2.35,3.37,5.07,9.72,13.8,17.9,22.};
	const string xTrkBins[trackPtBins+1] = {"TrkPt0p5","TrkPt0p7", "TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt999" };
	const string xCentBins[nCentBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
	
	//PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_Merged_refpt_newJetTrackCorrections_fineBin.root
	TFile *fJFFs = new TFile("JFFcorrs_sube0_newJFFs.root");
	TFile *fspillovers = new TFile("JFFcorrs_forSpillover_newJFFs.root");
	TFile *fGenGen = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_Merged_refpt_newJetTrackCorrections_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenReco = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_Merged_refpt_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenGenSube0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_sube0_Merged_refpt_newJetTrackCorrections_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenGenSubeNon0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_subeNon0_Merged_refpt_newJetTrackCorrections_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoSube0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_sube0_Merged_refpt_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenRecoSubeNon0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_subeNon0_Merged_refpt_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	/*TFile *fGenRecoQuark = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_%s%sQuarkReduced_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenRecoGluon = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_%s%sGluonReduced_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenGenQuark = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_%s%sQuarkReduced_refpt_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenGenGluon = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_%s%sGluonReduced_refpt_fineBin.root",jet.c_str(),track.c_str()));
	TFile *extra = new TFile("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_RecoGenReduced_ptCut126_fineBin.root");*/
	
	TFile *syst = new TFile(Form("relativeErrors_%s%s2%s%s.root",jet.c_str(),track.c_str(),jetDiv.c_str(),trackDiv.c_str()));
	
	TH2D *GenGenSignals[trackPtBins][nCentBins];
	TH2D *GenRecoSignals[trackPtBins][nCentBins];
	TH2D *GenGenBG[trackPtBins][nCentBins];
	TH2D *GenRecoBG[trackPtBins][nCentBins];
	
	TH2D *GenGenSignalsSub[trackPtBins][nCentBins];
	TH2D *GenRecoSignalsSub[trackPtBins][nCentBins];
	
	TH1D *GenGenDRDistr[trackPtBins][nCentBins];
	TH1D *GenRecoDRDistr[trackPtBins][nCentBins];
	TH1D *GenGenDRDistrRed[trackPtBins][nCentBins];
	TH1D *GenRecoDRDistrRed[trackPtBins][nCentBins];
	
	TH1D *GenGenDRDistrNorm[trackPtBins][nCentBins];
	TH1D *GenRecoDRDistrNorm[trackPtBins][nCentBins];
	
	TH1D *GenGenEtaProj[trackPtBins][nCentBins];
	TH1D *GenRecoEtaProj[trackPtBins][nCentBins];
	
	TH1D *GenGenPhiProj[trackPtBins][nCentBins];
	TH1D *GenRecoPhiProj[trackPtBins][nCentBins];
	
	TH1D *SpilloverEtaProj[trackPtBins][nCentBins];
	TH1D *JFFEtaProj[trackPtBins][nCentBins];
	
	TH2D *JFFcorrs[trackPtBins][nCentBins];
	TH2D *JFFcorrsToApply[trackPtBins][nCentBins];
	TH2D *spilloverCorrsToApply[trackPtBins][nCentBins];
	
	THStack *drGenGenDensity[nCentBins];
	THStack *drGenRecoDensity[nCentBins];
	
	TH1F *nJetsGen[nCentBins];
	TH1F *nJetsReco[nCentBins];
	
	TH1D *uncertainties[trackPtBins][nCentBins];
	TGraphErrors *unc[trackPtBins][nCentBins];
	
	for(int i=0; i<trackPtBins; i++){
		for(int j=0; j<nCentBins; j++){
			
			string toGet = Form("JFFcorrs_cent%d_pt%d",j,i);
			JFFcorrsToApply[i][j] = (TH2D*)fJFFs->Get(toGet.c_str())->Clone(toGet.c_str());
			spilloverCorrsToApply[i][j] = (TH2D*)fspillovers->Get(toGet.c_str())->Clone(Form("spillover_%d_%d",i,j));
		      //JFFcorrsToApply[i]->Scale(0.5);
		      
			if(drawUncertainties){
				toGet = Form("stackedRelatives_pt%d_cent%d",i,j);
				uncertainties[i][j] = (TH1D*)syst->Get(toGet.c_str())->Clone(toGet.c_str());
				unc[i][j] = new TGraphErrors(uncertainties[i][j]);
			}
			
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
			GenGenDRDistrRed[i][j] = new TH1D(Form("GenGenDR_%d_%dRed",i,j),"",14,xdrbins);
			GenRecoDRDistrRed[i][j] = new TH1D(Form("GenRecoDR_%d_%dRed",i,j),"",14,xdrbins);
			GenGenDRDistrRed[i][j]->SetMarkerColor(2);
			GenGenDRDistrRed[i][j]->SetLineColor(2);
			GenRecoDRDistrRed[i][j]->SetMarkerColor(2);
			GenRecoDRDistrRed[i][j]->SetLineColor(2);
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
			
			//This block is for debugging ONLY!!!
			/*for(int ibin=0; ibin<GenRecoSignalsSub[i][j]->GetNbinsX(); ibin++){
				for(int iybin=0; iybin<GenRecoSignalsSub[i][j]->GetNbinsY(); iybin++){
					GenRecoSignalsSub[i][j]->SetBinContent(ibin,iybin,GenGenSignalsSub[i][j]->GetBinContent(ibin,iybin));
				}
			}*/
			
			if(applyJFFcorrs){
				if(jet=="Reco"){
					//JFFcorrsToApply[i][j]->Scale(nJetsGen[j]->Integral());
					GenGenSignalsSub[i][j]->Add(JFFcorrsToApply[i][j],-1);
				}
				if(jetDiv=="Reco"){
					//JFFcorrsToApply[i][j]->Scale(nJetsReco[j]->Integral());
					GenRecoSignalsSub[i][j]->Add(JFFcorrsToApply[i][j],-1);	
				}	
			}
			
			//cout << "track pt" << i << " cent " << j << " gen counts: "<< GenGenSignalsSub[i][j]->GetEntries() << endl;
			//cout << "track pt" << i << " cent " << j << " counts: "<< GenRecoSignalsSub[i][j]->GetEntries() << endl;
			
			int lowBin = GenGenSignalsSub[i][j]->GetYaxis()->FindBin(-1);
			int highBin = GenGenSignalsSub[i][j]->GetYaxis()->FindBin(1);
			GenGenEtaProj[i][j] = GenGenSignalsSub[i][j]->ProjectionX(Form("GenGenEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			GenGenEtaProj[i][j]->Scale(1./GenGenSignalsSub[i][j]->GetXaxis()->GetBinWidth(2));
			GenGenEtaProj[i][j]->Rebin(4);
			GenGenEtaProj[i][j]->Scale(1./4.);
			GenRecoEtaProj[i][j] = GenRecoSignalsSub[i][j]->ProjectionX(Form("GenRecoEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			GenRecoEtaProj[i][j]->Scale(1./GenRecoSignalsSub[i][j]->GetXaxis()->GetBinWidth(2));
			GenRecoEtaProj[i][j]->Rebin(4);
			GenRecoEtaProj[i][j]->Scale(1./4.);
			
			lowBin = GenGenSignalsSub[i][j]->GetXaxis()->FindBin(-1);
			highBin = GenGenSignalsSub[i][j]->GetXaxis()->FindBin(1);
			GenGenPhiProj[i][j] = GenGenSignalsSub[i][j]->ProjectionY(Form("GenGenPhiProj_%d_%d",i,j),lowBin, highBin,"e");
			GenGenPhiProj[i][j]->Scale(1./GenGenSignalsSub[i][j]->GetYaxis()->GetBinWidth(2));
			GenGenPhiProj[i][j]->Rebin(4);
			GenGenPhiProj[i][j]->Scale(1./4.);
			GenRecoPhiProj[i][j] = GenRecoSignalsSub[i][j]->ProjectionY(Form("GenRecoPhiProj_%d_%d",i,j),lowBin, highBin,"e");
			GenRecoPhiProj[i][j]->Scale(1./GenRecoSignalsSub[i][j]->GetYaxis()->GetBinWidth(2));
			GenRecoPhiProj[i][j]->Rebin(4);
			GenRecoPhiProj[i][j]->Scale(1./4.);
			
			//cout << "at " << i << " " << j << " spillover int: "<< spilloverCorrsToApply[i][j]->Integral() << endl;
			//cout << " jff int: "<< JFFcorrsToApply[i][j]->Integral() << endl;
			
			lowBin = spilloverCorrsToApply[i][j]->GetYaxis()->FindBin(-1);
			highBin = spilloverCorrsToApply[i][j]->GetYaxis()->FindBin(1);
			SpilloverEtaProj[i][j] = spilloverCorrsToApply[i][j]->ProjectionX(Form("SpilloverEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			SpilloverEtaProj[i][j]->Scale(1./spilloverCorrsToApply[i][j]->GetXaxis()->GetBinWidth(2));
			SpilloverEtaProj[i][j]->Rebin(4);
			SpilloverEtaProj[i][j]->Scale(1./4.);
			JFFEtaProj[i][j] = JFFcorrsToApply[i][j]->ProjectionX(Form("JFFEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			JFFEtaProj[i][j]->Scale(1./JFFcorrsToApply[i][j]->GetXaxis()->GetBinWidth(2));
			JFFEtaProj[i][j]->Rebin(4);
			JFFEtaProj[i][j]->Scale(1./4.);
			
			//cout << "at " << i << " " << j << " spillover proj int: "<< SpilloverEtaProj[i][j]->Integral() << endl;
			//cout << " jff proj int: "<< JFFEtaProj[i][j]->Integral() << endl;
			
			if(doJFFCorrs){
				if(!applyJFFcorrs){
					JFFcorrs[i][j] = (TH2D*)GenRecoSignalsSub[i][j]->Clone(Form("JFFcorrs_cent%d_pt%d",j,i));
					JFFcorrs[i][j]->Add(GenGenSignalsSub[i][j], -1);
				}
				else{
					JFFcorrs[i][j] = JFFcorrsToApply[i][j];
				}
				//now symmetrize the corrections
				for(int ixbin=1; ixbin<=JFFcorrs[i][j]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<=JFFcorrs[i][j]->GetNbinsY(); iybin++){
						double xcent = JFFcorrs[i][j]->GetXaxis()->GetBinCenter(ixbin);
						double ycent = JFFcorrs[i][j]->GetYaxis()->GetBinCenter(iybin);
						int negX = JFFcorrs[i][j]->FindBin(-1*xcent,ycent);
						//int negY = JFFcorrs[i][j]->GetYaxis()->FindBin(-1*ycent);
						
						double dr = sqrt(pow(xcent,2)+pow(ycent,2));
						if(dr>0.8 || (dr>0.4 && i>4)){
							JFFcorrs[i][j]->SetBinContent(ixbin,iybin,0);
							JFFcorrs[i][j]->SetBinError(ixbin,iybin,0);
						}
						else{
							//add2DBin(ixbin, iybin, negX, JFFcorrs[i][j]);
							//add2DBin(ixbin, iybin, ixbin, negY, JFFcorrs[i][j]);
							//JFFcorrs[i][j]->SetBinContent(negX,iybin,JFFcorrs[i][j]->GetBinContent(ixbin,iybin));
							//JFFcorrs[i][j]->SetBinError(negX,iybin,JFFcorrs[i][j]->GetBinError(ixbin,iybin));							
						}
					}
				}
				//JFFcorrs[i][j]->Scale(0.5);
					
								
				for(int ixbin=1; ixbin<=JFFcorrs[i][j]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<=JFFcorrs[i][j]->GetNbinsY(); iybin++){
						double xcent = JFFcorrs[i][j]->GetXaxis()->GetBinCenter(ixbin);
						double ycent = JFFcorrs[i][j]->GetYaxis()->GetBinCenter(iybin);
						//int negX = JFFcorrs[i][j]->GetXaxis()->FindBin(-1*xcent);
						int negY = JFFcorrs[i][j]->FindBin(xcent,-1*ycent);
						
						double dr = sqrt(pow(xcent,2)+pow(ycent,2));
						if(dr>0.8 || (dr>0.4 && i>4)){
							JFFcorrs[i][j]->SetBinContent(ixbin,iybin,0);
							JFFcorrs[i][j]->SetBinError(ixbin,iybin,0);
						}
						else{
							//add2DBin(ixbin, iybin, negY, JFFcorrs[i][j]);
							//JFFcorrs[i][j]->SetBinContent(ixbin,negY,JFFcorrs[i][j]->GetBinContent(ixbin,iybin));
							//JFFcorrs[i][j]->SetBinError(ixbin,negY,JFFcorrs[i][j]->GetBinError(ixbin,iybin));	
												
						}
					}
					
				}
				//take average of symmetry
				//JFFcorrs[i][j]->Scale(0.5);
				
				//scale to make the corrections per jet
				//JFFcorrs[i][j]->Scale(1./nJetsReco[j]->Integral());
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
			//GenGenDRDistr[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			//GenRecoDRDistr[i][j]->Divide(GenRecoDRDistrNorm[i][j]);
			
			for(int ibin=1; ibin<=GenRecoDRDistr[i][j]->GetNbinsX(); ibin++){
				if(GenGenDRDistr[i][j]->GetBinContent(ibin) - 0.5*(GenGenDRDistr[i][j]->GetBinError(ibin)) <=0){
					GenGenDRDistr[i][j]->SetBinContent(ibin, 0);
					GenGenDRDistr[i][j]->SetBinError(ibin,0);
					GenGenDRDistrRed[i][j]->SetBinContent(ibin, 1);
					GenGenDRDistrRed[i][j]->SetBinError(ibin,0.001);
					GenRecoDRDistrRed[i][j]->SetBinContent(ibin, 1);
					GenRecoDRDistrRed[i][j]->SetBinError(ibin,0.001);
				}
				if(GenRecoDRDistr[i][j]->GetBinContent(ibin) - 0.5*(GenRecoDRDistr[i][j]->GetBinError(ibin)) <=0){
					GenRecoDRDistr[i][j]->SetBinContent(ibin, 0);
					GenRecoDRDistr[i][j]->SetBinError(ibin,0);
					GenGenDRDistrRed[i][j]->SetBinContent(ibin, 1);
					GenGenDRDistrRed[i][j]->SetBinError(ibin,0.001);
					GenRecoDRDistrRed[i][j]->SetBinContent(ibin, 1);
					GenRecoDRDistrRed[i][j]->SetBinError(ibin,0.001);
				}
			}
		}
	}
	
	TCanvas *ceta = new TCanvas("ceta","",1200,1600);
	TLatex *l2[trackPtBins][nCentBins];
	ceta->Divide(nCentBins,trackPtBins-2);
	for(int j=0; j<nCentBins; j++){
		for(int i=1; i<trackPtBins-1; i++){
			ceta->cd(i*nCentBins+(3-j)+1-4);
			GenRecoEtaProj[i][j]->GetYaxis()->SetNdivisions(505);
			GenRecoEtaProj[i][j]->GetYaxis()->SetLabelSize(0.15);
			//GenRecoEtaProj[i][j]->SetMinimum(2);//GenRecoEtaProj[i][0]->GetMinimum()*1.1);
			//GenRecoEtaProj[i][j]->SetMaximum(-2);//GenRecoEtaProj[i][0]->GetMaximum()*1.2);
			GenRecoEtaProj[i][j]->GetXaxis()->SetLabelSize(0.15);
			GenRecoEtaProj[i][j]->GetXaxis()->SetRangeUser(-2.5,2.5);
			GenRecoEtaProj[i][j]->GetYaxis()->SetRangeUser(-2,2);
			GenRecoEtaProj[i][j]->SetLineColor(2);
			GenRecoEtaProj[i][j]->Draw("");
			GenGenEtaProj[i][j]->Draw("same");
			SpilloverEtaProj[i][j]->SetLineColor(6);
			SpilloverEtaProj[i][j]->Draw("same");
			JFFEtaProj[i][j]->SetLineColor(kOrange);
			JFFEtaProj[i][j]->Draw("same");
			l2[i][j] = new TLatex(-0.5,GenRecoEtaProj[i][j]->GetMaximum()*0.85,Form("#Delta#eta, %s-%s, %g<pT<%g",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBinDouble[i],xTrkBinDouble[i+1]));
			l2[i][j]->SetTextSize(0.15);
			l2[i][j]->Draw("same");
		}
	}
	
	TCanvas *cphi = new TCanvas("cphi","",1200,1600);
	TLatex *l3[trackPtBins][nCentBins];
	cphi->Divide(nCentBins,trackPtBins-2);
	for(int j=0; j<nCentBins; j++){
		for(int i=1; i<trackPtBins-1; i++){
			cphi->cd(i*nCentBins+(3-j)+1-4);
			GenRecoPhiProj[i][j]->GetYaxis()->SetNdivisions(505);
			GenRecoPhiProj[i][j]->GetYaxis()->SetLabelSize(0.15);
			GenRecoPhiProj[i][j]->GetXaxis()->SetLabelSize(0.15);
			GenRecoPhiProj[i][j]->GetXaxis()->SetRangeUser(-1.7,3.2);
			GenRecoPhiProj[i][j]->SetLineColor(2);
			GenRecoPhiProj[i][j]->Draw("");
			GenGenPhiProj[i][j]->Draw("same");
			l3[i][j] = new TLatex(-0.5,GenRecoPhiProj[i][j]->GetMaximum()*0.85,Form("#Delta#phi, %s-%s, %g<pT<%g",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBinDouble[i],xTrkBinDouble[i+1]));
			l3[i][j]->SetTextSize(0.15);
			l3[i][j]->Draw("same");
		}
	}
	
	TCanvas *cc = new TCanvas("cc","",1200,800);
	cc->Divide(1,3);
	cc->cd(1);
	GenRecoSignalsSub[4][0]->Draw("colz");
	cc->cd(2);
	TH2D *tmp = (TH2D*)GenRecoSignalsSub[4][0]->Clone("tmp");
	tmp->Add(JFFcorrsToApply[4][0],1);
	tmp->Draw("colz");
	cc->cd(3);
	JFFcorrsToApply[4][0]->Draw("colz");
	
	TCanvas *c12 = new TCanvas("c12","",600,600);
	c12->cd();
	GenRecoEtaProj[1][3]->GetXaxis()->SetRangeUser(-3,3);
	GenRecoEtaProj[1][3]->Draw("");
	
	TFile *fout2 = new TFile("JEScorrTest.root");
	//fout2->cd();
	TH1D *RecoGenPt126[trackPtBins][nCentBins];
	TLatex *l1[trackPtBins][nCentBins];
	TCanvas *cprint = new TCanvas("cprint","",1200,1600);
	cprint->Divide(nCentBins, trackPtBins-1);
	for(int j=0; j<nCentBins; j++){
		for(int i=1; i<trackPtBins; i++){
			cprint->cd(i*nCentBins+(3-j)+1-4);
			
			GenRecoDRDistr[i][j]->SetMarkerStyle(25);
			GenGenDRDistrRed[i][j]->SetMarkerStyle(24);
			GenRecoDRDistrRed[i][j]->SetMarkerStyle(24);
			
			//RecoGenPt126[i][j] = (TH1D*)fout2->Get(Form("GenGenDR_%d_%d",i,j))->Clone(Form("GenGenDR_%d_%d",i,j));
			//RecoGenPt126[i][j]->SetMarkerStyle(25);
			
			for(int k=1; k<=GenGenDRDistr[i][j]->GetNbinsX(); k++){
				//GenGenDRDistr[i][j]->SetBinContent(k,GenGenDRDistr[i][j]->GetBinContent(k)/GenGenDRDistr[i][j]->GetBinWidth(k));
				//GenGenDRDistr[i][j]->SetBinError(k,GenGenDRDistr[i][j]->GetBinError(k)/GenGenDRDistr[i][j]->GetBinWidth(k));
			
				//GenRecoDRDistr[i][j]->SetBinContent(k,GenRecoDRDistr[i][j]->GetBinContent(k)/GenRecoDRDistr[i][j]->GetBinWidth(k));
				//GenRecoDRDistr[i][j]->SetBinError(k,GenRecoDRDistr[i][j]->GetBinError(k)/GenRecoDRDistr[i][j]->GetBinWidth(k));
			}
			
			GenGenDRDistr[i][j]->Divide(GenGenDRDistr[i][j],GenRecoDRDistr[i][j],1,1,"B");
			
			//GenGenDRDistr[i][j]->Add(GenRecoDRDistr[i][j],-1);
			
			GenGenDRDistr[i][j]->SetMaximum(1.3);
			GenGenDRDistr[i][j]->SetMinimum(0.7);
			if(i<3){
				GenGenDRDistr[i][j]->SetMaximum(1.3);
				GenGenDRDistr[i][j]->SetMinimum(0.7);
			}
			GenGenDRDistr[i][j]->GetYaxis()->SetNdivisions(505);
			GenGenDRDistr[i][j]->GetYaxis()->SetLabelSize(0.15);
			GenGenDRDistr[i][j]->GetXaxis()->SetLabelSize(0.15);
			//GenGenDRDistr[i][j]->SetYTitle("RecoReco / RecoGen #rho(#DeltaR)");
			GenGenDRDistr[i][j]->Draw("");
			//GenRecoDRDistr[i][j]->Draw("same");
			if(drawUncertainties){
				for(int ibin=0; ibin<unc[i][j]->GetN(); ibin++){
				//unc[i][j]->SetPoint(ibin, GenGenDRDistr[i][j]->GetBinCenter(ibin), GenGenDRDistr[i][j]->GetBinContent(ibin));
					unc[i][j]->SetPoint(ibin, GenGenDRDistr[i][j]->GetBinCenter(ibin+1), 1.);
					if(GenGenDRDistr[i][j]->GetBinContent(ibin+1)) unc[i][j]->SetPointError(ibin, GenGenDRDistr[i][j]->GetBinWidth(ibin+1)/2., uncertainties[i][j]->GetBinContent(ibin+1)*GenGenDRDistr[i][j]->GetBinContent(ibin+1));
					else unc[i][j]->SetPointError(ibin, GenGenDRDistr[i][j]->GetBinWidth(ibin+1)/2., uncertainties[i][j]->GetBinContent(ibin+1));
				}
				unc[i][j]->Draw("2,5,same");
			}
			GenGenDRDistr[i][j]->Draw("same");
			//RecoGenPt126[i][j]->Draw("Same");
			l1[i][j] = new TLatex(0.2,GenGenDRDistr[i][j]->GetMaximum()*0.95,Form("%s-%s, %g<pT<%g",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBinDouble[i],xTrkBinDouble[i+1]));
			l1[i][j]->SetTextSize(0.15);
			l1[i][j]->Draw("same");
			GenRecoDRDistrRed[i][j]->Draw("Same");
			
		}
	}
	fout2->Close();
	string jfc = "jffCorrApplied";
	if(!applyJFFcorrs) jfc = "";
	//cprint->SaveAs(Form("%s%s_div_%s%s_bgSub%d_%s_newCorrections.pdf",jet.c_str(), track.c_str(), jetDiv.c_str(), trackDiv.c_str(), doBGsub, jfc.c_str()));
	
	if(doJFFCorrs){
		TFile *fout = new TFile("JFFcorrs_fullsube_newJFFs.root","recreate");
		fout->cd();
		for(int i=0; i<trackPtBins; i++){
			for(int j=0; j<nCentBins; j++){
				JFFcorrs[i][j]->Write();
			}
		}
		fout->Close();
	}
	if(writeOutCorrelations){
		TFile *foutFinal = new TFile("TotalCorrelations.root","recreate");
		for(int i=0; i<trackPtBins; i++){
			for(int j=0; j<nCentBins; j++){
				GenGenSignalsSub[i][j]->Write();
				GenRecoSignalsSub[i][j]->Write();
			}
		}
		foutFinal->Close();
	}
}