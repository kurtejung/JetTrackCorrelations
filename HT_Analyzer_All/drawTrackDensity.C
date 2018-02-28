#include "TFile.h"
#include "TH2D.h"
#include "TH1F.h"
#include "THStack.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include <string>
#include <iostream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TLine.h"

using namespace std;

void refillFits(TF1 *etaFit, TF1 *phiFit, TH2D *result, int flip){
	
	TF2* ClosureFit = new TF2("ClosureFit", "[0]/2/TMath::Pi()/[1]/[2]*TMath::Exp(-1.*(x*x/[1]/[1]/2))*TMath::Exp(-1.*(y*y/[2]/[2]/2))",-5.,5.,-TMath::Pi()/2.,3*TMath::Pi()/2.);
	
	double par0;
	if(flip) par0 = etaFit->GetParameter(0);
	else par0 = etaFit->GetParameter(1);
	cout << "par0 " << par0 << endl;
	double par1= etaFit->GetParameter(2);
	cout << "par1" << par1 << endl;
	double par2= phiFit->GetParameter(2);
	cout << "par2" << par2 << endl;
	
	ClosureFit->SetParameter(0,par0);
	ClosureFit->SetParameter(1,par1);
	ClosureFit->SetParameter(2,par2);
	
	//cout << "closure" << ClosureFit->Integral(-1,1,-1,1) << endl;
	
	//cout << " refillFits 1: "<< result->ProjectionX()->Integral() << endl;
	
	for(int ixbin=1; ixbin<=result->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=result->GetNbinsY(); iybin++){
			double xcent = result->GetXaxis()->GetBinCenter(ixbin);
			double ycent = result->GetYaxis()->GetBinCenter(iybin);
			result->SetBinContent(ixbin,iybin,ClosureFit->Eval(xcent,ycent));
			result->SetBinError(ixbin, iybin, 1e-12);
		}
	}
	
	//cout << " refillFits 2: "<< result->ProjectionX()->Integral() << endl;
	
	double width_temp_x = result->GetXaxis()->GetBinWidth(1);
	double width_temp_y = result->GetYaxis()->GetBinWidth(1);
	
	result->Scale(width_temp_x*width_temp_y);
	//cout << " refillFits 3: "<< result->ProjectionX()->Integral() << endl;
}

void refillTH1Fits(TH1D *etaHist, TH1D *phiHist, TH2D *result){
	TF1 *etaFit = new TF1("etaFit","gaus",-1.5,1.5);
	TF1 *phiFit = new TF1("phiFit","gaus",-1.5,1.5);
	etaHist->Fit(etaFit,"qN0");
	phiHist->Fit(phiFit,"qN0");
	refillFits(etaFit,phiFit,result,0);
}

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

TH2D *backgroundSubtraction(TH2D *signal){
	
	double sidebandlow = 1.5;
	double sidebandhigh = 2.5;
	
	int lowBin = signal->GetXaxis()->FindBin(-1*sidebandhigh);
	int highBin = signal->GetXaxis()->FindBin(-1*sidebandlow);
	int totalBins = highBin-lowBin;
	TH1D *inclPhiBG = (TH1D*)signal->ProjectionY("inclPhiBG",lowBin, highBin,"e");
	
	lowBin = signal->GetXaxis()->FindBin(sidebandlow);
	highBin = signal->GetXaxis()->FindBin(sidebandhigh);
	totalBins+=(highBin-lowBin);
	TH1D *tmp1 = (TH1D*)signal->ProjectionY("tmp1",lowBin, highBin,"e");
	inclPhiBG->Add(tmp1);
	
	inclPhiBG->Scale(1./(double)(totalBins+2)); //SHOULD BE totalBins+2!!!!!
	//delete tmp1;
	
	lowBin = signal->GetYaxis()->FindBin(-5);
	highBin = signal->GetYaxis()->FindBin(5);
	totalBins = highBin-lowBin;
	TH1D *inclPhiFG = signal->ProjectionY("inclPhiFG",lowBin, highBin,"e");
	inclPhiFG->Scale(1./inclPhiFG->GetNbinsX());
	assert(signal->GetNbinsY() == inclPhiBG->GetNbinsX());
	
	for(int ixbin=1; ixbin<=signal->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=signal->GetNbinsY(); iybin++){
			if(inclPhiBG->GetBinContent(iybin)){
				signal->SetBinContent(ixbin,iybin, signal->GetBinContent(ixbin,iybin) - inclPhiBG->GetBinContent(iybin));
				signal->SetBinError(ixbin,iybin, ReturnSubError(signal->GetBinError(ixbin,iybin), inclPhiBG->GetBinError(iybin)));
			}
		}
	}
	
	return signal;
}


TH2D* backgroundSubtract(TH2D *signal, TH2D *background, TH1F *nJets, bool doBGsub, int doMixEvt, bool scaleJets){
	
	//Smooth background to be phi-independent
	const int nsegments = 1;
	TH1D *projEta[nsegments];
	for(int i=0; i<nsegments; i++){
		projEta[i] = (TH1D*)background->ProjectionX(Form("projEta_%d",i),(int)(1+(i*background->GetNbinsY()/(double)nsegments)),(int)(1+(i+1)*(background->GetNbinsY()/(double)nsegments)),"e");
	}
	/*for(int ixbin=1; ixbin<=background->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=background->GetNbinsY(); iybin++){
			background->SetBinContent(ixbin,iybin, projEta->GetBinContent(ixbin));
			background->SetBinError(ixbin,iybin, projEta->GetBinError(ixbin)*sqrt(background->GetNbinsY()));
		}
	}*/
	
	
	TF1 *littleFit = new TF1("littleFit","pol0",-0.3,0.3);
	if(doMixEvt==2){
		projEta[nsegments/2]->Fit(littleFit,"qN0","",-0.3,0.3);
		double normFactor = littleFit->GetParameter(0);
	      background->Scale(1./normFactor);//projEta->GetBinContent(projEta->FindBin(0.)));
	      for(int i=0; i<nsegments; i++){
	            projEta[i]->Scale(1./normFactor);//projEta->GetBinContent(projEta->FindBin(0)));
	      }
	}
	else if(doMixEvt){
		projEta[0]->Reset();
		int lowBin = signal->GetYaxis()->FindBin(1.5);
		int hiBin = signal->GetYaxis()->FindBin(2.3);
		projEta[0] = (TH1D*)signal->ProjectionX("projEta",lowBin, hiBin,"e");
		projEta[0]->Fit(littleFit,"qN0","",-0.3,0.3);
		double normFactor = littleFit->GetParameter(0);
		projEta[0]->Scale(1./normFactor);
	}
	
	//Do mixed event correction
	for(int ixbin=1; ixbin<=signal->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=signal->GetNbinsY(); iybin++){
			if(doMixEvt){
				int segment = (double)nsegments*((double)iybin/(double)signal->GetNbinsY()) - 1;
				if(doMixEvt==1) segment=0;
				if(projEta[segment]->GetBinContent(ixbin)) signal->SetBinContent(ixbin,iybin, signal->GetBinContent(ixbin,iybin)/projEta[segment]->GetBinContent(ixbin));
				else signal->SetBinContent(ixbin,iybin,0);
				if(projEta[segment]->GetBinContent(ixbin)) signal->SetBinError(ixbin,iybin, ReturnDivError(signal->GetBinContent(ixbin,iybin), signal->GetBinError(ixbin,iybin), projEta[segment]->GetBinContent(ixbin),projEta[segment]->GetBinError(ixbin)));
				
				if(projEta[segment]->GetBinContent(ixbin)) background->SetBinContent(ixbin,iybin, projEta[segment]->GetBinContent(ixbin));
				else background->SetBinContent(ixbin,iybin,0);
			}
		}
	}
	if(scaleJets) signal->Scale(1./nJets->Integral());
	
	//Get the dphi projection to remove v2, etc
	if(doBGsub) signal = backgroundSubtraction(signal);
	
	for(int i=0; i<nsegments; i++){
		//delete projEta[i];
	}
	
	return signal;
}

//Mix event 2 = full mixing correction, 1 = sideband correction, 0 = no correction
void drawTrackDensity(string jet="Gen", string track="Gen", string jetDiv="Gen", string trackDiv="Gen", bool doBGsub=false, int doMixEvt=2, int scaleJets=1){
	
	bool doJFFCorrs = false;
	bool doSpillover = false;
	bool applyJFFcorrs = false;
	bool applySpillovers = false;
	bool writeOutCorrelations = false;
	bool drawUncertainties = false;
	bool writeDRHistosForPlotting = false;
	bool doptWeight = false;
	
	bool doMixEvtShapeCorr = false;
	if(doMixEvtShapeCorr) doBGsub = false; //need to derive wing shape corrections in proportion to BG
	
	bool isData1 = false;
	bool isData2 = false;
	
	bool isBjets1 = false;
	bool isBjets2 = false;
	
	string bjetString1 = "", bjetString2="";
	if(isBjets1) bjetString1 = "b";
	if(isBjets2) bjetString2 = "b";
	string trkCorr = "_notrkcorr";
	string trkCorrDiv = "_notrkcorr";
	if(track=="Reco" || isData1) trkCorr = "";
	if(trackDiv=="Reco" || isData2) trkCorrDiv = "";
	string ptWeightString = "";
	if(doptWeight){
		ptWeightString = "_pTweighted";
		trkCorr = "";
		trkCorrDiv = "";
	}
	//if(jet=="Reco" || jetDiv=="Reco") applyJFFcorrs = true;
	//if(jet=="Reco" && jetDiv=="Reco") applyJFFcorrs = false;
		
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	
	const double xdrbins[24] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.9,1.,1.2,1.4,1.6,1.8,2.};
	
	const int trackPtBins = 10;
	const int nCentBins = 5;
	
	const double xTrkBinDouble[trackPtBins+1] = {0.5,0.7,1.,2.,3.,4.,8.,12.,16.,20.,999.};
	//const float mean_pts[trackPtBins] = {0.57,0.844,1.35,2.35,5.07,9.72,17.9};
	const float mean_pts[trackPtBins] = {0.57,0.844,1.35,2.35,3.37,5.07,9.72,13.8,17.9,22.};
	const double ptWidth[trackPtBins] = {0.2,0.3,1,1,1,4,4,4,4,200};
	
	const string xTrkBins[trackPtBins+1] = {"TrkPt0p5", "TrkPt0p7", "TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt999" };
	const string halliexTrkBins[trackPtBins+1] = {"TrkPt05", "TrkPt07", "TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
	//const string bjetxTrkBins[trackPtBins+1] = {"TrkPt0p5", "TrkPt0p7", "TrkPt1", "TrkPt2", "TrkPt4", "TrkPt8", "TrkPt16", "TrkPt999" };
	const string xCentBins[nCentBins+1] = {"Cent0","Cent10", "Cent30", "Cent50", "Cent70", "Cent100"};
	const string xOldCentBins[nCentBins+1] = {"Cent0", "Cent10", "Cent30", "Cent50", "Cent100", "Cent100"};
	//const string bjetxCentBins[nCentBins+1] = {"Cent0","Cent30","Cent100"};
	
	//PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_Merged_refpt_newJetTrackCorrections_fineBin.root
	TFile *fJFFs = new TFile("JFFcorrs_PythHydjet_sube0_finalJFFJEC_CymbalTune_withFits.root");
	TFile *fJFFsPythia = new TFile("JFFcorrs_PythiPP_newJFFJEC_withFits.root");
	TFile *fSpillovers = new TFile("JFFcorrs_PythHydjet_subeNon0_finalJFFJEC_CymbalTune_withFits.root");
	TFile *fHSpillovers = new TFile("~/Downloads/Inclusive_Hydjet_SpillOvers.root");
	TFile *fHJFFs = new TFile("~/Downloads/Inclusive_Hydjet_JFFResiduals.root");
	
	TFile *fJFFsWeighted = new TFile("JFFcorrs_PythHydjet_pTweighted_sube0_newJFFJEC_withFits.root");
	TFile *fJFFsPythiaWeighted = new TFile("JFFcorrs_PythiaOnly_pTweighted_sube0_newJFFJEC_withFits.root");
	
	//TFile *fSpillovers = new TFile("JFFcorrs_subeNon0_TestOldJFFs_withFits.root");
	
	//gengen = black, genreco = red
	TFile *fGenGen = new TFile(Form("root_output/PbPb_5TeVMC_%s%s_QuarkJet_sube0_corrHistos_merged_withMix_summed.root",jet.c_str(),track.c_str()));
	TFile *fGenReco = new TFile(Form("root_output/PbPb_5TeVMC_%s%s_GluonJet_sube0_corrHistos_merged_withMix_summed.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenGenSube0 = new TFile(Form("root_output/PbPb_5TeVMC_withMix_%s%s_sube0_CymbalTune_inclJetBinning.root",jet.c_str(),track.c_str()));
	TFile *fGenGenSubeNon0 = new TFile(Form("root_output/PbPb_5TeVMC_%s%sSube0_fineBinTrkCorrs_withTrkCorrEtaSymmV2_finalJFF_noMix.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoSube0 = new TFile(Form("root_output/PbPb_5TeVMC_%s%sSube0_fineBinTrkCorrs_withTrkCorrEtaSymmV2_finalJFF_noMix.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenRecoSubeNon0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_subeNon0_Merged_withMix_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenRecoQuarkSube0 = new TFile(Form("PbPb_5TeVMC_%s%s_QuarkJet_sube0_corrHistos_merged_withMix_summed.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenRecoGluonSube0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_GluonJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenGenQuarkSube0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_QuarkJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenGenGluonSube0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_GluonJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoQuarkSubeNon0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_QuarkJet_subeNon0_Merged_withMix_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenRecoGluonSubeNon0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_GluonJet_subeNon0_Merged_withMix_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenGenQuark = new TFile(Form("root_output/PbPb_5TeVMC_%s%sSube0_QuarkJets_fixMatch_fineBinTrkCorrs_withTrkCorrEtaSymmV2_finalJFF_noMix.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoQuark = new TFile(Form("root_output/PbPb_5TeVMC_%s%s_sube0_QuarkOnly_noMix_finalJFFs.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenGenGluon = new TFile(Form("root_output/PbPb_5TeVMC_%s%sSube0_GluonJets_fixMatch_fineBinTrkCorrs_withTrkCorrEtaSymmV2_finalJFF_noMix.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoGluon = new TFile(Form("root_output/PbPb_5TeVMC_%s%s_sube0_GluonOnly_noMix_finalJFFs.root",jetDiv.c_str(),trackDiv.c_str()));
		
	TFile *fGenGenPythia = new TFile(Form("root_output/pp_Pythia6MC_%s%s_withMix_inclJetBinning_finalJFFs.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoPythia = new TFile(Form("root_output/pp_Pythia6MC_%s%s_withMix_inclJetBinning_finalJFFs.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *extra = new TFile("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenRecoReduced_Merged_noMix_newJetTrackCorrections_pt16to20TrkCorrFix_fineBin.root");
	
	TFile *hallieGenGen = new TFile("root_output/Pythia_GenJet_GenTrack_Aug23.root");
	TFile *hallieGenReco = new TFile("root_output/HallieFiles/HydJet_GenJet_RecoTrack_Aug23.root");
	
	TFile *ppData = new TFile("root_output/pp_Data_withMix_inclJetBinning_finalJFFs_60vzMixBins.root"); //"genGen"
	TFile *PbPbData = new TFile("root_output/PbPb_5TeVData_fixMix_fineBinTrkCorrs_withTrkCorrEtaSymmV2_finalJFFs.root"); //"RecoReco"
	TFile *PbPbData2 = new TFile("root_output/PbPb_Data_withMix_CymbalTune_fullCymbalCorrs_inclJetBinning_withPtWeightME.root");
	
	TFile *hallieMix = new TFile("Hallie_PbPb_Inclusive_Correlations.root");//Mixed_Event_Cent10_Cent30_Pt100_Pt1000_TrkPt07_TrkPt1")
	
	TFile *fMix = new TFile("MixEvtFlatnessCorr.root");
	
	TFile *syst = new TFile(Form("relativeErrors_%s%s2%s%s.root",jetDiv.c_str(),trackDiv.c_str(),jet.c_str(),track.c_str()));
	
	TH2D *GenGenSignals[trackPtBins][nCentBins];
	TH2D *GenGenSignalsPtWeight[trackPtBins][nCentBins];
	TH2D *GenRecoSignals[trackPtBins][nCentBins];
	TH2D *GenGenBG[trackPtBins][nCentBins];
	TH2D *GenGenBGPtWeight[trackPtBins][nCentBins];
	TH2D *GenRecoBG[trackPtBins][nCentBins];
	
	TH2D *GenGenSignalsSub[trackPtBins][nCentBins];
	TH2D *GenGenPtWeightSignalsSub[trackPtBins][nCentBins];
	TH2D *GenRecoSignalsSub[trackPtBins][nCentBins];
	
	TH1D *GenGenDRDistr[trackPtBins][nCentBins];
	TH1D *GenRecoDRDistr[trackPtBins][nCentBins];
	TH1D *GenGenDRDistrRed[trackPtBins][nCentBins];
	TH1D *GenRecoDRDistrRed[trackPtBins][nCentBins];
	
	TH1D *GenGenDRDistrNorm[trackPtBins][nCentBins];
	TH1D *GenRecoDRDistrNorm[trackPtBins][nCentBins];
	
	TH1D *GenGenEtaProj[trackPtBins][nCentBins];
	TH1D *GenRecoEtaProj[trackPtBins][nCentBins];
	TH1D *GenGenEtaSidebandProj[trackPtBins][nCentBins];
	TH1D *GenGenPtWeightEtaSidebandProj[trackPtBins][nCentBins];
	TH1D *GenGenBGProj[trackPtBins][nCentBins];
	TH1D *GenRecoEtaSidebandProj[trackPtBins][nCentBins];
	TF1 *mixEventAfterburnerGen[trackPtBins][nCentBins];
	TF1 *mixEventAfterburnerReco[trackPtBins][nCentBins];
	TH1D *GenGenEtaSidebandXchk[trackPtBins][nCentBins];
	
	TH1D *GenGenPhiProj[trackPtBins][nCentBins];
	TH1D *GenRecoPhiProj[trackPtBins][nCentBins];
	
	TH1D *SpilloverEtaProj[trackPtBins][nCentBins];
	TH1D *HSpilloverEtaProj[trackPtBins][nCentBins];
	TH1D *JFFEtaProj[trackPtBins][nCentBins];
	TH1D *JFFDRProj[trackPtBins][nCentBins];
	
	TH1D *correctionEtaProj[trackPtBins][nCentBins];
	TH1D *correctionPhiProj[trackPtBins][nCentBins];
	TF1 *corrEtaProj[trackPtBins][nCentBins];
	TF1 *corrPhiProj[trackPtBins][nCentBins];
	
	TH2D *JFFcorrs[trackPtBins][nCentBins];
	TH2D *JFFcorrsToApply[trackPtBins][nCentBins];
	TH2D *JFFcorrsToApplyWeighted[trackPtBins][nCentBins];
	TH2D *JFFcorrsToApplyPyth[trackPtBins][nCentBins];
	TH2D *JFFcorrsToApplyPythWeighted[trackPtBins][nCentBins];
	TH2D *spilloverCorrsToApply[trackPtBins][nCentBins];
	TH2D *HspilloverCorrsToApply[trackPtBins][nCentBins];
	
	TH2D *HJFFCorrsToApply[trackPtBins][nCentBins];
	
	THStack *drGenGenDensity[nCentBins];
	THStack *drGenRecoDensity[nCentBins];
	
	TH1F *nJetsGen[nCentBins];
	TH1F *nJetsReco[nCentBins];
	
	TH1D *uncertainties[trackPtBins][nCentBins];
	TGraphErrors *unc[trackPtBins][nCentBins];
	
	for(int i=0; i<trackPtBins; i++){
		for(int j=0; j<nCentBins; j++){
			
			int jbin=j;
			if(j>3) jbin=3;
			string toGet = Form("JFFcorrs_cent%d_pt%d",jbin,i);
			if(applyJFFcorrs){
				JFFcorrsToApply[i][j] = (TH2D*)fJFFs->Get(toGet.c_str())->Clone(toGet.c_str());
				JFFcorrsToApplyWeighted[i][j] = (TH2D*)fJFFsWeighted->Get(toGet.c_str())->Clone(toGet.c_str());
				
				JFFcorrsToApplyPyth[i][j] = (TH2D*)fJFFsPythia->Get(toGet.c_str())->Clone(toGet.c_str());
				JFFcorrsToApplyPythWeighted[i][j] = (TH2D*)fJFFsPythiaWeighted->Get(toGet.c_str())->Clone(toGet.c_str());
			}
			//spilloverCorrsToApply[i][j] = new TH2D(Form("Closure_2D_%s_%s_Pt100_Pt300_%s_%s",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);
			HspilloverCorrsToApply[i][j] = new TH2D(Form("Hallie_Closure_2D_%s_%s_Pt100_Pt300_%s_%s",xOldCentBins[j].c_str(),xOldCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);
			HJFFCorrsToApply[i][j] = new TH2D(Form("Hallie_JFFClosure_2D_%s_%s_Pt100_Pt300_%s_%s",xOldCentBins[j].c_str(),xOldCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);
			if(i>0 && applyJFFcorrs && j<3){
				if(applySpillovers){
					toGet = Form("Eta_SpillOver_Fit_%s_%s_Pt100_Pt300_%s_%s",xOldCentBins[j].c_str(),xOldCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
					TF1 *etaFitIn2 = (TF1*)fHSpillovers->Get(toGet.c_str())->Clone(Form("HetaFit_%d_%d",i,j));
					toGet = Form("Phi_SpillOver_Fit_%s_%s_Pt100_Pt300_%s_%s",xOldCentBins[j].c_str(),xOldCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
					TF1 *phiFitIn2 = (TF1*)fHSpillovers->Get(toGet.c_str())->Clone(Form("HphiFit_%d_%d",i,j));
					
					refillFits(etaFitIn2, phiFitIn2, HspilloverCorrsToApply[i][j],0);
					HspilloverCorrsToApply[i][j]->Scale(ptWidth[i]);
					if(i>4 && j>2) HspilloverCorrsToApply[i][j]->Reset();
				}
				toGet = Form("Raw_JFF_Residual_Eta_%s_%s_Pt100_Pt1000_%s_%s",xOldCentBins[j].c_str(),xOldCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
				TH1D *etaFitIn3 = (TH1D*)fHJFFs->Get(toGet.c_str())->Clone(Form("JFFetaFit_%d_%d",i,j));
				etaFitIn3->Scale(1./(double)HJFFCorrsToApply[i][j]->GetNbinsX());
				toGet = Form("Raw_JFF_Residual_Phi_%s_%s_Pt100_Pt1000_%s_%s",xOldCentBins[j].c_str(),xOldCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
				TH1D *phiFitIn3 = (TH1D*)fHJFFs->Get(toGet.c_str())->Clone(Form("JFFphiFit_%d_%d",i,j));
				phiFitIn3->Scale(1./(double)HJFFCorrsToApply[i][j]->GetNbinsY());
				refillTH1Fits(etaFitIn3, phiFitIn3, HJFFCorrsToApply[i][j]);
				
				/*for(int ixbin=1; ixbin<HJFFCorrsToApply[i][j]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<HJFFCorrsToApply[i][j]->GetNbinsY(); iybin++){
						double content = (etaFitIn3->GetBinContent(etaFitIn3->FindBin(HJFFCorrsToApply[i][j]->GetXaxis()->GetBinCenter(ixbin))) + 
						phiFitIn3->GetBinContent(phiFitIn3->FindBin(HJFFCorrsToApply[i][j]->GetYaxis()->GetBinCenter(iybin))))/2.;
						HJFFCorrsToApply[i][j]->SetBinContent(ixbin,iybin,content);
						HJFFCorrsToApply[i][j]->SetBinError(ixbin,iybin,1e-12);
						
					}
				}*/
				HJFFCorrsToApply[i][j]->Scale(ptWidth[i]);
			}
			if(applySpillovers){
				int jbin=j;
				if(j>3) jbin=3;
				toGet = Form("JFFcorrs_cent%d_pt%d",jbin,i);
				spilloverCorrsToApply[i][j] = (TH2D*)fSpillovers->Get(toGet.c_str())->Clone(Form("spillover_%d_%d",i,j));
				if(i>10) spilloverCorrsToApply[i][j]->Reset();
		             //JFFcorrsToApply[i]->Scale(0.5);
		      }
		      
			if(drawUncertainties){
				toGet = Form("stackedRelatives_pt%d_cent%d",i,j);
				uncertainties[i][j] = (TH1D*)syst->Get(toGet.c_str())->Clone(toGet.c_str());
				unc[i][j] = new TGraphErrors(uncertainties[i][j]);
			}
			
			if(isData1) toGet = Form("Data_h%sJetTrackSignalBackground%s%s%s_%s_Pt100_Pt1000_%s_%s",bjetString1.c_str(),ptWeightString.c_str(), trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			else toGet = Form("%sJet_%sTrack_h%sJetTrackSignalBackground%s%s%s_%s_Pt100_Pt1000_%s_%s",jet.c_str(),track.c_str(),bjetString1.c_str(),ptWeightString.c_str(), trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			cout << "getting " << toGet << endl;
			GenGenSignals[i][j] = (TH2D*)fGenGen->Get(toGet.c_str())->Clone(toGet.c_str());
			
			if(isData1) toGet = Form("Data_h%sJetTrackSignalBackground_pTweighted%s%s_%s_Pt100_Pt1000_%s_%s",bjetString1.c_str(), trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			else toGet = Form("%sJet_%sTrack_h%sJetTrackSignalBackground%s%s%s_%s_Pt100_Pt1000_%s_%s",jet.c_str(),track.c_str(),bjetString1.c_str(),ptWeightString.c_str(), trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			cout << "getting " << toGet << endl;
			GenGenSignalsPtWeight[i][j] = (TH2D*)fGenGen->Get(toGet.c_str())->Clone(toGet.c_str());
			
			if(isData2 /*&& (i==0 || j>2)*/) toGet = Form("Data_h%sJetTrackSignalBackground%s%s%s_%s_Pt100_Pt1000_%s_%s",bjetString2.c_str(),ptWeightString.c_str(), trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			//else if(isData2) toGet = Form("Raw_Yield_%s_%s_Pt100_Pt1000_%s_%s",xCentBins[j].c_str(),xCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
			else toGet = Form("%sJet_%sTrack_h%sJetTrackSignalBackground%s%s%s_%s_Pt100_Pt1000_%s_%s",jetDiv.c_str(),trackDiv.c_str(),bjetString2.c_str(),ptWeightString.c_str(),trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			cout << "getting " << toGet << endl;
			/*if(i==0 || j>2)*/ GenRecoSignals[i][j] = (TH2D*)fGenReco->Get(toGet.c_str())->Clone(toGet.c_str());
			//else GenRecoSignals[i][j] = (TH2D*)hallieMix->Get(toGet.c_str())->Clone(toGet.c_str());
			
			if(isData1) toGet = Form("Data_h%sJetTrackME%s%s_%s_Pt100_Pt1000_%s_%s",bjetString1.c_str(),trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			else toGet = Form("%sJet_%sTrack_h%sJetTrackME%s%s_%s_Pt100_Pt1000_%s_%s",jet.c_str(),track.c_str(),bjetString1.c_str(),trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			//else toGet = Form("Data_h%sJetTrackME%s%s_%s_Pt100_Pt300_%s_%s",bjetString1.c_str(),trkCorr.c_str(),xOldCentBins[j].c_str(),xOldCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			cout << "getting " << toGet << endl;
			GenGenBG[i][j] = (TH2D*)fGenGen->Get(toGet.c_str())->Clone(toGet.c_str());
			
			if(isData2 /*&& (i==0 || j>2)*/) toGet = Form("Data_h%sJetTrackME%s%s_%s_Pt100_Pt1000_%s_%s",bjetString2.c_str(),trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			//else if(isData2) toGet = Form("Mixed_Event_%s_%s_Pt100_Pt1000_%s_%s",xCentBins[j].c_str(),xCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
			else toGet = Form("%sJet_%sTrack_h%sJetTrackME%s%s_%s_Pt100_Pt1000_%s_%s",jetDiv.c_str(),trackDiv.c_str(),bjetString2.c_str(),trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			cout << "getting " << toGet << endl;
			/*if(i==0 || j>2)*/ GenRecoBG[i][j] = (TH2D*)fGenReco->Get(toGet.c_str())->Clone(toGet.c_str());
			//else GenRecoBG[i][j] = (TH2D*)hallieMix->Get(toGet.c_str())->Clone(toGet.c_str());
			
			GenGenDRDistr[i][j] = new TH1D(Form("GenGenDR_%d_%d",i,j),"",19,xdrbins);
			GenRecoDRDistr[i][j] = new TH1D(Form("GenRecoDR_%d_%d",i,j),"",19,xdrbins);
			GenGenDRDistrRed[i][j] = new TH1D(Form("GenGenDR_%d_%dRed",i,j),"",19,xdrbins);
			GenRecoDRDistrRed[i][j] = new TH1D(Form("GenRecoDR_%d_%dRed",i,j),"",19,xdrbins);
			GenGenDRDistrRed[i][j]->SetMarkerColor(2);
			GenGenDRDistrRed[i][j]->SetLineColor(2);
			GenRecoDRDistrRed[i][j]->SetMarkerColor(2);
			GenRecoDRDistrRed[i][j]->SetLineColor(2);
			GenGenDRDistrNorm[i][j] = new TH1D(Form("GenGenDRNorm_%d_%d",i,j),"",19,xdrbins);
			GenRecoDRDistrNorm[i][j] = new TH1D(Form("GenRecoDRNorm_%d_%d",i,j),"",19,xdrbins);
			
			JFFDRProj[i][j] = new TH1D(Form("JFFDRProj_%d_%d",i,j),"",19,xdrbins); JFFDRProj[i][j]->Sumw2();
		}
	}
	
	
	for(int j=0; j<nCentBins; j++){
		string toGet;
		if(isData1) toGet = Form("Data_all_%sjets_corrpT%s_%s_Pt100_Pt1000",bjetString1.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str());
		else toGet = Form("%sJet_%sTrack_all_%sjets_corrpT%s_%s_Pt100_Pt1000",jet.c_str(),track.c_str(),bjetString1.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str());
		cout << "getting " << toGet << endl;
		nJetsGen[j] = (TH1F*)fGenGen->Get(toGet.c_str())->Clone(toGet.c_str());
		
		if(isData2) toGet = Form("Data_all_%sjets_corrpT%s_%s_Pt100_Pt1000",bjetString2.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str());
		else toGet = Form("%sJet_%sTrack_all_%sjets_corrpT%s_%s_Pt100_Pt1000",jetDiv.c_str(),trackDiv.c_str(),bjetString2.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str());
		cout << "getting " << toGet << endl;
		nJetsReco[j] = (TH1F*)fGenReco->Get(toGet.c_str())->Clone(toGet.c_str());
		
		drGenGenDensity[j] = new THStack(Form("drGenGenDensity_%d",j),"");
		drGenRecoDensity[j] = new THStack(Form("drGenRecoDensity_%d",j),"");
	}
	
	
	for(int i=0; i<trackPtBins; i++){
		for(int j=0; j<nCentBins; j++){
			cout << " bg sub GenGen for " << i << " " << j << endl;
			GenGenSignalsSub[i][j] = backgroundSubtract(GenGenSignals[i][j], GenGenBG[i][j], nJetsGen[j], doBGsub, doMixEvt, scaleJets);
			cout << " bg sub GenReco for " << i << " " << j << endl;
			GenRecoSignalsSub[i][j] = backgroundSubtract(GenRecoSignals[i][j], GenRecoBG[i][j], nJetsReco[j], doBGsub, doMixEvt, scaleJets);
			//GenRecoSignalsSub[i][j] = (TH2D*)GenRecoSignals[i][j]->Clone(Form("GenRecoSignalsSub_%d_%d",i,j));
			GenGenPtWeightSignalsSub[i][j] = backgroundSubtract(GenGenSignalsPtWeight[i][j], GenRecoBG[i][j], nJetsReco[j], doBGsub, doMixEvt, scaleJets);
			
			//This block is for debugging ONLY!!!
			/*for(int ibin=0; ibin<GenRecoSignalsSub[i][j]->GetNbinsX(); ibin++){
				for(int iybin=0; iybin<GenRecoSignalsSub[i][j]->GetNbinsY(); iybin++){
					GenRecoSignalsSub[i][j]->SetBinContent(ibin,iybin,GenGenSignalsSub[i][j]->GetBinContent(ibin,iybin));
				}
			}*/
			
			cout << " at applyJFF block" << endl;
			if(applyJFFcorrs){
				if(jet=="Reco"){
					if(!doptWeight){
						GenGenSignalsSub[i][j]->Add(JFFcorrsToApply[i][j],-1);
						if(applySpillovers && i<10) GenGenSignalsSub[i][j]->Add(HspilloverCorrsToApply[i][j],-1);
					}
					else{
						GenGenSignalsSub[i][j]->Add(JFFcorrsToApplyPythWeighted[i][j],-1);
						if(applySpillovers && i<10) GenGenSignalsSub[i][j]->Add(spilloverCorrsToApply[i][j],-1);
					}
				}
				if(jetDiv=="Reco"){
					if(!doptWeight){
						GenRecoSignalsSub[i][j]->Add(JFFcorrsToApply[i][j],-1);
						if(applySpillovers && i<10) GenRecoSignalsSub[i][j]->Add(HspilloverCorrsToApply[i][j],-1);	
					}
					else{
						GenRecoSignalsSub[i][j]->Add(JFFcorrsToApplyWeighted[i][j],-1);
						if(applySpillovers && i<10) GenRecoSignalsSub[i][j]->Add(spilloverCorrsToApply[i][j],-1);
					}
				}
			}
			
			cout << "track pt" << i << " cent " << j << " gen counts: "<< GenGenSignalsSub[i][j]->GetEntries() << endl;
			cout << "track pt" << i << " cent " << j << " counts: "<< GenRecoSignalsSub[i][j]->GetEntries() << endl;
			
			int lowBin = GenGenSignalsSub[i][j]->GetYaxis()->FindBin(-1.); //-1
			int highBin = GenGenSignalsSub[i][j]->GetYaxis()->FindBin(1.); //1
			GenGenEtaProj[i][j] = GenGenSignalsSub[i][j]->ProjectionX(Form("GenGenEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			//GenGenEtaProj[i][j]->Scale(1./GenGenSignalsSub[i][j]->GetYaxis()->GetBinWidth(2));
			GenGenEtaProj[i][j]->Rebin(5);
			GenGenEtaProj[i][j]->Scale(1./5.);
			GenRecoEtaProj[i][j] = GenRecoSignalsSub[i][j]->ProjectionX(Form("GenRecoEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			//GenRecoEtaProj[i][j]->Scale(1./GenRecoSignalsSub[i][j]->GetYaxis()->GetBinWidth(2));
			GenRecoEtaProj[i][j]->Rebin(5);
			GenRecoEtaProj[i][j]->Scale(1./5.);
			
			lowBin = GenGenSignalsSub[i][j]->GetYaxis()->FindBin(1.2); //1.4
			highBin = GenGenSignalsSub[i][j]->GetYaxis()->FindBin(2.2); //2.5
			GenGenEtaSidebandProj[i][j] = GenGenSignalsSub[i][j]->ProjectionX(Form("GenGenEtaSidebandProj_%d_%d",i,j),lowBin, highBin,"e");
			GenRecoEtaSidebandProj[i][j] = GenRecoSignalsSub[i][j]->ProjectionX(Form("GenRecoEtaSidebandProj_%d_%d",i,j),lowBin, highBin,"e");
			double sclFactor = (double)(highBin-lowBin)+1.;
			GenGenEtaSidebandProj[i][j]->Rebin(5);
			GenGenEtaSidebandProj[i][j]->Scale(1./(5.*sclFactor));
			GenRecoEtaSidebandProj[i][j]->Rebin(5);
			GenRecoEtaSidebandProj[i][j]->Scale(1./(5.*sclFactor));
			
			GenGenPtWeightEtaSidebandProj[i][j] = GenGenPtWeightSignalsSub[i][j]->ProjectionX(Form("GenGenPtWeightEtaSidebandProj_%d_%d",i,j),lowBin, highBin,"e");
			GenGenPtWeightEtaSidebandProj[i][j]->Rebin(5);
			GenGenPtWeightEtaSidebandProj[i][j]->Scale(1./(5.*sclFactor));
			
			lowBin = GenGenSignalsSub[i][j]->GetXaxis()->FindBin(-1); //-1
			highBin = GenGenSignalsSub[i][j]->GetXaxis()->FindBin(1); //1
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
			
			cout << " at spillover block" << endl;
			if(applySpillovers){
				lowBin = spilloverCorrsToApply[i][j]->GetYaxis()->FindBin(-1.5);
				highBin = spilloverCorrsToApply[i][j]->GetYaxis()->FindBin(1.5);
				SpilloverEtaProj[i][j] = spilloverCorrsToApply[i][j]->ProjectionX(Form("SpilloverEtaProj_%d_%d",i,j),lowBin, highBin,"e");
				//SpilloverEtaProj[i][j]->Scale(1./spilloverCorrsToApply[i][j]->GetYaxis()->GetBinWidth(2));
				SpilloverEtaProj[i][j]->Rebin(5);
				SpilloverEtaProj[i][j]->Scale(1./5.);
				/*HSpilloverEtaProj[i][j] = HspilloverCorrsToApply[i][j]->ProjectionX(Form("HSpilloverEtaProj_%d_%d",i,j),lowBin, highBin,"e");
				HSpilloverEtaProj[i][j]->Scale(1./HspilloverCorrsToApply[i][j]->GetYaxis()->GetBinWidth(2));
				HSpilloverEtaProj[i][j]->Rebin(5);
				HSpilloverEtaProj[i][j]->Scale(1./5.);*/
			}
			if(applyJFFcorrs){
				JFFEtaProj[i][j] = JFFcorrsToApply[i][j]->ProjectionX(Form("JFFEtaProj_%d_%d",i,j),lowBin, highBin,"e");
				//JFFEtaProj[i][j]->Scale(1./JFFcorrsToApply[i][j]->GetYaxis()->GetBinWidth(2));
				JFFEtaProj[i][j]->Rebin(5);
				JFFEtaProj[i][j]->Scale(1./5.);
			}
			//cout << "at " << i << " " << j << " spillover proj int: "<< SpilloverEtaProj[i][j]->Integral() << endl;
			//cout << " jff proj int: "<< JFFEtaProj[i][j]->Integral() << endl;
			
			if(doMixEvtShapeCorr){
				
				double shiftFactor = GenGenEtaSidebandProj[i][j]->GetBinContent(GenGenEtaSidebandProj[i][j]->FindBin(0));
				//for(int ibin=1; ibin<=GenGenEtaSidebandProj[i][j]->GetNbinsX(); ibin++){ GenGenEtaSidebandProj[i][j]->AddBinContent(ibin, shiftFactor); }
				//if(shiftFactor>0) GenGenEtaSidebandProj[i][j]->Scale(1./shiftFactor);
				mixEventAfterburnerGen[i][j] = new TF1(Form("mixEventAfterburnerGen_%d_%d",i,j),"[2]+TMath::Abs([0]*TMath::Sin([1]*x))",-2.5,2.5);
				//mixEventAfterburnerGen[i][j]->SetParLimits(0,0.,0.1);
				mixEventAfterburnerGen[i][j]->SetParLimits(1,0.7,1.3);
				mixEventAfterburnerGen[i][j]->SetParLimits(2,0.7,1.3);
				GenGenEtaSidebandProj[i][j]->Fit(mixEventAfterburnerGen[i][j],"B","",-2.5,2.5);
				
				double shiftFactor2 = GenRecoEtaSidebandProj[i][j]->GetBinContent(GenRecoEtaSidebandProj[i][j]->FindBin(0));
				//for(int ibin=1; ibin<=GenGenPtWeightEtaSidebandProj[i][j]->GetNbinsX(); ibin++){ GenGenPtWeightEtaSidebandProj[i][j]->AddBinContent(ibin, shiftFactor); }
				//if(shiftFactor2>0) GenRecoEtaSidebandProj[i][j]->Scale(1./shiftFactor2);
				mixEventAfterburnerReco[i][j] = new TF1(Form("mixEventAfterburnerReco_%d_%d",i,j),"[2]+TMath::Abs([0]*TMath::Sin([1]*x))",-2.5,2.5);
				mixEventAfterburnerReco[i][j]->SetParLimits(1,0.7,1.3);
				mixEventAfterburnerReco[i][j]->SetParLimits(2,0.7,1.3);
				GenRecoEtaSidebandProj[i][j]->Fit(mixEventAfterburnerReco[i][j],"BN","",-2.5,2.5);
				
				TH2D *tmpGenGen = (TH2D*)GenGenBG[i][j]->Clone("tempGenGen");
				TH2D *tmpGenReco = (TH2D*)GenRecoBG[i][j]->Clone("tempGenReco");
				
				for(int ixbin=1; ixbin<=tmpGenGen->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<=tmpGenGen->GetNbinsY(); iybin++){
						double binCenter = tmpGenGen->GetXaxis()->GetBinCenter(ixbin);
						tmpGenGen->SetBinContent(ixbin, iybin, mixEventAfterburnerGen[i][j]->Eval(binCenter));
						tmpGenGen->SetBinError(ixbin,iybin,1e-50);
						
						tmpGenReco->SetBinContent(ixbin, iybin, mixEventAfterburnerReco[i][j]->Eval(binCenter));
						tmpGenReco->SetBinError(ixbin,iybin,1e-50);
					}
				}
				
				if(shiftFactor && (j<3 && i<6)) GenGenBG[i][j]->Multiply(tmpGenGen);
				//GenGenBGPtWeight[i][j] = (TH2D*)GenRecoBG[i][j]->Clone(Form("Data_h%sJetTrackME%s%s_ptWeightMixEvtAfterburner_%s_Pt100_Pt1000_%s_%s",bjetString2.c_str(),trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()));
				if(shiftFactor2 && (j<3 && i<6)) GenRecoBG[i][j]->Multiply(tmpGenReco);
				
				for(int ixbin=1; ixbin<=GenGenEtaSidebandProj[i][j]->GetNbinsX(); ixbin++){
					double binCenter = GenGenEtaSidebandProj[i][j]->GetXaxis()->GetBinCenter(ixbin);
					//GenGenEtaSidebandProj[i][j]->SetBinContent(ixbin, GenGenEtaSidebandProj[i][j]->GetBinContent(ixbin)/mixEventAfterburnerGen[i][j]->Eval(binCenter));
					//GenRecoEtaSidebandProj[i][j]->SetBinContent(ixbin, GenRecoEtaSidebandProj[i][j]->GetBinContent(ixbin)/mixEventAfterburnerReco[i][j]->Eval(binCenter));
				}
			}
			
			cout << "at derive spillover block" << endl;	
			if(doJFFCorrs){

				JFFcorrs[i][j] = (TH2D*)GenRecoSignalsSub[i][j]->Clone(Form("JFFcorrs_cent%d_pt%d",j,i));
				if(!doSpillover) JFFcorrs[i][j]->Add(GenGenSignalsSub[i][j], -1);
				//JFFcorrs[i][j] = backgroundSubtraction(JFFcorrs[i][j]);

				//now symmetrize the corrections
				for(int ixbin=1; ixbin<=JFFcorrs[i][j]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<=JFFcorrs[i][j]->GetNbinsY(); iybin++){
						double xcent = JFFcorrs[i][j]->GetXaxis()->GetBinCenter(ixbin);
						double ycent = JFFcorrs[i][j]->GetYaxis()->GetBinCenter(iybin);
						int negX = JFFcorrs[i][j]->FindBin(-1*xcent,ycent);
						//int negY = JFFcorrs[i][j]->GetYaxis()->FindBin(-1*ycent);
						
						double dr = sqrt(pow(xcent,2)+pow(ycent,2));
						//if(dr>1.5 || (dr>1.5 && i<4)){
						//	JFFcorrs[i][j]->SetBinContent(ixbin,iybin,0);
						//	JFFcorrs[i][j]->SetBinError(ixbin,iybin,0);
						//}
						//else{
							//add2DBin(ixbin, iybin, negX, JFFcorrs[i][j]);
							//add2DBin(ixbin, iybin, ixbin, negY, JFFcorrs[i][j]);
							//JFFcorrs[i][j]->SetBinContent(negX,iybin,JFFcorrs[i][j]->GetBinContent(ixbin,iybin));
							//JFFcorrs[i][j]->SetBinError(negX,iybin,JFFcorrs[i][j]->GetBinError(ixbin,iybin));							
						//}
					}
				}
				//JFFcorrs[i][j]->Scale(0.5);
					
								
				/*for(int ixbin=1; ixbin<=JFFcorrs[i][j]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<=JFFcorrs[i][j]->GetNbinsY(); iybin++){
						double xcent = JFFcorrs[i][j]->GetXaxis()->GetBinCenter(ixbin);
						double ycent = JFFcorrs[i][j]->GetYaxis()->GetBinCenter(iybin);
						//int negX = JFFcorrs[i][j]->GetXaxis()->FindBin(-1*xcent);
						int negY = JFFcorrs[i][j]->FindBin(xcent,-1*ycent);
						
						double dr = sqrt(pow(xcent,2)+pow(ycent,2));
						if(dr>1.5 || (dr>1.5 && i>4)){
							JFFcorrs[i][j]->SetBinContent(ixbin,iybin,0);
							JFFcorrs[i][j]->SetBinError(ixbin,iybin,0);
						}
						else{
							//add2DBin(ixbin, iybin, negY, JFFcorrs[i][j]);
							//JFFcorrs[i][j]->SetBinContent(ixbin,negY,JFFcorrs[i][j]->GetBinContent(ixbin,iybin));
							//JFFcorrs[i][j]->SetBinError(ixbin,negY,JFFcorrs[i][j]->GetBinError(ixbin,iybin));	
												
						}
					}
					
				}*/
				cout << "starting projections" << endl;
				correctionEtaProj[i][j] = JFFcorrs[i][j]->ProjectionX(Form("correctionEtaProj%d_%d",i,j),0, -1,"e");
				correctionPhiProj[i][j] = JFFcorrs[i][j]->ProjectionY(Form("correctionPhiProj%d_%d",i,j),0, -1,"e");
				correctionEtaProj[i][j]->Rebin(5);
				correctionPhiProj[i][j]->Rebin(5);
				correctionEtaProj[i][j]->Scale(1./(5.*(double)JFFcorrs[i][j]->GetYaxis()->GetBinWidth(2)));
				correctionPhiProj[i][j]->Scale(1./(5.*(double)JFFcorrs[i][j]->GetXaxis()->GetBinWidth(2)));
				corrEtaProj[i][j] = new TF1(Form("fitcorrEtaProj_%d_%d",i,j),"gaus",-1,1);
				corrPhiProj[i][j] = new TF1(Form("fitcorrPhiProj_%d_%d",i,j),"gaus",-1,1);
				corrEtaProj[i][j]->SetParLimits(0,0,2);
				corrEtaProj[i][j]->FixParameter(1,0);
				corrEtaProj[i][j]->SetParLimits(2,9e-2,5);
				corrPhiProj[i][j]->SetParLimits(0,0,2);
				corrPhiProj[i][j]->FixParameter(1,0);
				corrPhiProj[i][j]->SetParLimits(2,9e-2,5);
				correctionEtaProj[i][j]->Fit(corrEtaProj[i][j],"B","",-1,1);
				correctionPhiProj[i][j]->Fit(corrPhiProj[i][j],"B","",-1,1);
				//take average of symmetry
				//JFFcorrs[i][j]->Scale(0.5);
				
				//scale to make the corrections per jet
				//JFFcorrs[i][j]->Scale(1./nJetsReco[j]->Integral());
			}
						
			cout << " plotting dr projections" << endl;
			for(int ixbin=1; ixbin<=GenGenSignalsSub[i][j]->GetNbinsX(); ixbin++){
				for(int iybin=1; iybin<=GenGenSignalsSub[i][j]->GetNbinsY(); iybin++){
					double xcent = GenGenSignalsSub[i][j]->GetXaxis()->GetBinCenter(ixbin);
					double ycent = GenGenSignalsSub[i][j]->GetYaxis()->GetBinCenter(iybin);
					double content = GenGenSignalsSub[i][j]->GetBinContent(GenGenSignalsSub[i][j]->GetBin(ixbin,iybin));
					//cout << "gen gen bin " << ixbin << ", " << iybin << " content " << content << endl;
					double dr = sqrt(pow(xcent,2)+pow(ycent,2));
					//if(xcent<0 || ycent<0) continue;
					GenGenDRDistr[i][j]->Fill(dr,content);
					GenGenDRDistrNorm[i][j]->Fill(dr);
					
					if(j<8) drGenGenDensity[j]->Add(GenGenDRDistr[i][j]);
										
					xcent = GenRecoSignalsSub[i][j]->GetXaxis()->GetBinCenter(ixbin);
					ycent = GenRecoSignalsSub[i][j]->GetYaxis()->GetBinCenter(iybin);
					content = GenRecoSignalsSub[i][j]->GetBinContent(GenRecoSignalsSub[i][j]->GetBin(ixbin,iybin));
					//cout << "gen reco bin " << ixbin << ", " << iybin << " content " << content << endl;
					dr = sqrt(pow(xcent,2)+pow(ycent,2));
					GenRecoDRDistr[i][j]->Fill(dr,content);
					GenRecoDRDistrNorm[i][j]->Fill(dr);
					
					if(j<8) drGenRecoDensity[j]->Add(GenRecoDRDistr[i][j]);
					
					if(applyJFFcorrs) content = JFFcorrsToApply[i][j]->GetBinContent(JFFcorrsToApply[i][j]->GetBin(ixbin, iybin));
					JFFDRProj[i][j]->Fill(dr, content);
					
				}
			}
			//GenGenDRDistr[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			//GenRecoDRDistr[i][j]->Divide(GenRecoDRDistrNorm[i][j]);
			
			cout << "at draw uncertainties" << endl;
			if(drawUncertainties){
				for(int ibin=1; ibin<=GenGenDRDistr[i][j]->GetNbinsX(); ibin++){
					if(uncertainties[i][j]->GetBinContent(ibin) > GenGenDRDistr[i][j]->GetBinContent(ibin)/(xdrbins[i+1]-xdrbins[i]) ||
					   (GenGenDRDistr[i][j]->GetBinContent(ibin-1)==0 && ibin>1)) {
						 GenGenDRDistr[i][j]->SetBinContent(ibin, 0);
						GenGenDRDistr[i][j]->SetBinError(ibin,0);
						GenGenDRDistrRed[i][j]->SetBinContent(ibin, 1);
						GenGenDRDistrRed[i][j]->SetBinError(ibin,0.001);
						GenRecoDRDistrRed[i][j]->SetBinContent(ibin, 1);
						GenRecoDRDistrRed[i][j]->SetBinError(ibin,0.001);
					}
					/*if(uncertainties[i][j]->GetBinContent(ibin) > GenRecoDRDistr[i][j]->GetBinContent(ibin)){
						GenRecoDRDistr[i][j]->SetBinContent(ibin, 0);
						GenRecoDRDistr[i][j]->SetBinError(ibin,0);
						GenGenDRDistrRed[i][j]->SetBinContent(ibin, 1);
						GenGenDRDistrRed[i][j]->SetBinError(ibin,0.001);
						GenRecoDRDistrRed[i][j]->SetBinContent(ibin, 1);
						GenRecoDRDistrRed[i][j]->SetBinError(ibin,0.001);
					}*/
				}
			}
		}
	}
	
	cout << "starting to draw..." << endl;
	TCanvas *ceta = new TCanvas("ceta","",1200,1600);
	TH1D *GenGenBGEtaProj[trackPtBins][nCentBins];
	TH1D *GenRecoBGEtaProj[trackPtBins][nCentBins];
	TLatex *l2[trackPtBins][nCentBins];
	ceta->Divide(nCentBins,(trackPtBins-1));
	TLine *lineAt1 = new TLine(-2.5,1,2.5,1);
	lineAt1->SetLineStyle(2);
	for(int j=0; j<nCentBins; j++){
		for(int i=1; i<trackPtBins; i++){
			ceta->cd(i*nCentBins+(3-j)+1-4);
			GenGenEtaProj[i][j]->GetYaxis()->SetNdivisions(505);
			GenGenEtaProj[i][j]->GetYaxis()->SetLabelSize(0.15);
			//GenGenEtaProj[i][j]->SetMinimum(0.5);//GenRecoEtaProj[i][0]->GetMinimum()*1.1);
			//GenGenEtaProj[i][j]->SetMaximum(1.5);//GenRecoEtaProj[i][0]->GetMaximum()*1.2);
			GenGenEtaProj[i][j]->GetXaxis()->SetLabelSize(0.15);
			GenGenEtaProj[i][j]->GetXaxis()->SetRangeUser(-3.5,3.5);
			//GenRecoEtaProj[i][j]->GetYaxis()->SetRangeUser(0.5,1.5);
			GenGenEtaProj[i][j]->SetLineColor(1);
			GenRecoEtaProj[i][j]->SetLineColor(2);
			GenRecoEtaProj[i][j]->SetMarkerColor(2);
			//GenGenEtaProj[i][j]->Divide(GenRecoEtaProj[i][j]);
			//GenGenEtaSidebandProj[i][j]->Draw("");
			//GenRecoEtaSidebandProj[i][j]->Draw("same");
			GenGenEtaProj[i][j]->Draw("");
			GenRecoEtaProj[i][j]->Draw("same");
			
			int loBin = GenRecoEtaProj[i][j]->FindBin(-1);
			int hiBin = GenRecoEtaProj[i][j]->FindBin(1);
			cout << "cent bin " << j << " pt bin " << i << " integral: "<< GenRecoEtaProj[i][j]->Integral(loBin,hiBin) << endl;
			
			GenGenBGEtaProj[i][j] = (TH1D*)GenGenBG[i][j]->ProjectionX()->Clone(Form("ggBGEtaProj_%d_%d",i,j));
			GenRecoBGEtaProj[i][j] = (TH1D*)GenRecoBG[i][j]->ProjectionX()->Clone(Form("grBGEtaProj_%d_%d",i,j));
			GenGenBGEtaProj[i][j]->SetLineColor(kGreen+1);
			GenGenBGEtaProj[i][j]->SetLineStyle(2);
			GenGenBGEtaProj[i][j]->Rebin(5);
			GenRecoBGEtaProj[i][j]->SetLineColor(kOrange);
			GenRecoBGEtaProj[i][j]->SetLineStyle(kOrange);
			//GenGenBGEtaProj[i][j]->Draw("same");
			int sclBin = GenGenBGEtaProj[i][j]->FindBin(-1.5);
			int sclBin2 = GenGenEtaProj[i][j]->FindBin(-1.5);
			cout << "scale factor for bg " << i << " " << j << " " << GenGenEtaProj[i][j]->GetBinContent(sclBin2) << " " << GenGenBGEtaProj[i][j]->GetBinContent(sclBin) << endl;
			GenGenBGEtaProj[i][j]->Scale(GenGenEtaProj[i][j]->GetBinContent(sclBin2)/GenGenBGEtaProj[i][j]->GetBinContent(sclBin));
			//GenGenBGEtaProj[i][j]->Draw("same");
			
			if(applySpillovers){
				SpilloverEtaProj[i][j]->SetLineColor(6);
				SpilloverEtaProj[i][j]->SetMarkerColor(6);
				SpilloverEtaProj[i][j]->SetMarkerSize(0.2);
			
			//SpilloverEtaProj[i][j]->Draw("same");
				/*HSpilloverEtaProj[i][j]->SetLineColor(kGreen+1);
				HSpilloverEtaProj[i][j]->SetMarkerColor(kGreen+1);
				HSpilloverEtaProj[i][j]->SetMarkerSize(0.2);*/
			//HSpilloverEtaProj[i][j]->Draw("same");
			}
			if(applyJFFcorrs){
				JFFEtaProj[i][j]->SetLineColor(kOrange);
				//JFFEtaProj[i][j]->Draw("same");
			}
			l2[i][j] = new TLatex(-2,GenRecoEtaProj[i][j]->GetMaximum()*0.97,Form("#Delta#eta, %s-%s, %g<pT<%g",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBinDouble[i],xTrkBinDouble[i+1]));
			l2[i][j]->SetTextSize(0.15);
			l2[i][j]->Draw("same");
			//lineAt1->Draw("Same");
			
			int lb = GenGenEtaProj[i][j]->FindBin(-0.4);
			int hb = GenGenEtaProj[i][j]->FindBin(0.4);
			double integralGenGen = GenGenEtaProj[i][j]->Integral(lb,hb);
			double integralGenReco = GenRecoEtaProj[i][j]->Integral(lb,hb);
			
			cout << "cent bin " << j << " pt: "<< i << " Scale: " << integralGenGen/integralGenReco << endl;
			
		}
	}
	
	TH1D *GenGenBGPhiProj[trackPtBins][nCentBins];
	TH1D *GenRecoBGPhiProj[trackPtBins][nCentBins];
	TCanvas *cphi = new TCanvas("cphi","",1200,1600);
	TLatex *l3[trackPtBins][nCentBins];
	cphi->Divide(nCentBins,trackPtBins-1);
	for(int j=0; j<nCentBins; j++){
		for(int i=1; i<trackPtBins; i++){
			cphi->cd(i*nCentBins+(3-j)+1-4);
			GenGenPhiProj[i][j]->GetYaxis()->SetNdivisions(505);
			GenGenPhiProj[i][j]->GetYaxis()->SetLabelSize(0.15);
			GenGenPhiProj[i][j]->GetXaxis()->SetLabelSize(0.15);
			GenGenPhiProj[i][j]->GetXaxis()->SetRangeUser(-1.7,3.2);
			GenRecoPhiProj[i][j]->SetLineColor(2);
			GenRecoPhiProj[i][j]->Draw("");
			//GenGenPhiProj[i][j]->Divide(GenRecoPhiProj[i][j]);
			GenGenPhiProj[i][j]->Draw("same");
			l3[i][j] = new TLatex(-0.5,GenRecoPhiProj[i][j]->GetMaximum()*0.97,Form("#Delta#phi, %s-%s, %g<pT<%g",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBinDouble[i],xTrkBinDouble[i+1]));
			l3[i][j]->SetTextSize(0.15);
			l3[i][j]->Draw("same");
			
			GenGenBGPhiProj[i][j] = (TH1D*)GenGenBG[i][j]->ProjectionY()->Clone(Form("ggBGPhiProj_%d_%d",i,j));
			GenRecoBGPhiProj[i][j] = (TH1D*)GenRecoBG[i][j]->ProjectionY()->Clone(Form("grBGPhiProj_%d_%d",i,j));
		}
	}
	
	TCanvas *cc = new TCanvas("cc","",1200,800);
	cc->Divide(1,3);
	cc->cd(1);
	GenRecoSignalsSub[4][0]->Draw("colz");
	cc->cd(2);
	TH2D *tmp = (TH2D*)GenRecoSignalsSub[4][0]->Clone("tmp");
	if(applyJFFcorrs){
		tmp->Add(JFFcorrsToApply[4][0],1);
		tmp->Draw("colz");
	}
	cc->cd(3);
	if(applyJFFcorrs) JFFcorrsToApply[4][0]->Draw("colz");
	
	TCanvas *c12 = new TCanvas("c12","",600,600);
	c12->cd();
	GenRecoEtaProj[1][3]->GetXaxis()->SetRangeUser(-3,3);
	GenRecoEtaProj[1][3]->Draw("");
	
	//TFile *fout2 = new TFile("JEScorrTest.root");
	//fout2->cd();
	//TH1D *RecoGenPt126[trackPtBins][nCentBins];
	TLatex *l1[trackPtBins][nCentBins];
	TCanvas *cprint = new TCanvas("cprint","",1200,1600);
	cprint->Divide(nCentBins, trackPtBins-2);
	for(int j=0; j<nCentBins; j++){
		for(int i=1; i<trackPtBins-1; i++){
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
			
			
			if(i<5) GenGenDRDistr[i][j]->Divide(GenGenDRDistr[i][j],GenRecoDRDistr[i][j],1,1,"B");
			else GenGenDRDistr[i][j]->Divide(GenGenDRDistr[i][j],GenRecoDRDistr[i][j],1,1,"B");
			
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
			GenGenDRDistr[i][j]->GetXaxis()->SetRangeUser(0,1);
			//GenGenDRDistr[i][j]->SetYTitle("RecoReco / RecoGen #rho(#DeltaR)");
			GenGenDRDistr[i][j]->Draw("");
			//GenRecoDRDistr[i][j]->Draw("same");
			//uncertainties[i][j]->Scale(mean_pts[i]*(xdrbins[i+1]-xdrbins[i]));
			//uncertainties[i][j]->Draw("same");
			if(drawUncertainties){
				for(int ibin=0; ibin<unc[i][j]->GetN(); ibin++){
				//unc[i][j]->SetPoint(ibin, GenGenDRDistr[i][j]->GetBinCenter(ibin), GenGenDRDistr[i][j]->GetBinContent(ibin));
					unc[i][j]->SetPoint(ibin, GenGenDRDistr[i][j]->GetBinCenter(ibin+1), 1.);
					if(GenGenDRDistr[i][j]->GetBinContent(ibin+1)) unc[i][j]->SetPointError(ibin, GenGenDRDistr[i][j]->GetBinWidth(ibin+1)/2., uncertainties[i][j]->GetBinContent(ibin+1));
					else unc[i][j]->SetPointError(ibin, GenGenDRDistr[i][j]->GetBinWidth(ibin+1)/2., uncertainties[i][j]->GetBinContent(ibin+1));
				}
				unc[i][j]->Draw("2,5,same");
			}
			GenGenDRDistr[i][j]->Draw("same");
			//RecoGenPt126[i][j]->Draw("Same");
			l1[i][j] = new TLatex(0.2,GenGenDRDistr[i][j]->GetMaximum()*0.85,Form("%s-%s, %g<pT<%g",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBinDouble[i],xTrkBinDouble[i+1]));
			l1[i][j]->SetTextSize(0.15);
			l1[i][j]->Draw("same");
			GenRecoDRDistrRed[i][j]->Draw("Same");
			
		}
	}
	

	if(doMixEvtShapeCorr){
		TCanvas *ccMixEvtCorr = new TCanvas("ccMixEvtCorr","",1200,1600);
		ccMixEvtCorr->Divide(nCentBins, trackPtBins-2);
		for(int j=0; j<nCentBins; j++){
			for(int i=1; i<trackPtBins-1; i++){
				ccMixEvtCorr->cd(i*nCentBins+(3-j)+1-4);
				
				GenGenEtaSidebandProj[i][j]->GetYaxis()->SetNdivisions(505);
				GenGenEtaSidebandProj[i][j]->GetYaxis()->SetLabelSize(0.15);
				GenGenEtaSidebandProj[i][j]->GetXaxis()->SetLabelSize(0.15);
				GenGenEtaSidebandProj[i][j]->GetXaxis()->SetRangeUser(-2.5,2.5);
				GenGenEtaSidebandProj[i][j]->Draw();
				//mixEventAfterburnerGen[i][j]->Draw("");
				
				GenRecoEtaSidebandProj[i][j]->SetLineColor(kMagenta+2);
				GenRecoEtaSidebandProj[i][j]->Draw("same");
				//mixEventAfterburnerReco[i][j]->Draw("same");
				
			}
		}
		
		TFile *fMixEvtCorr = new TFile("MixEvtFlatnessCorr.root","recreate");
		fMixEvtCorr->cd();
		
		for(int j=0; j<nCentBins; j++){
			nJetsGen[j]->Write();
			for(int i=0; i<trackPtBins; i++){
				GenGenSignals[i][j]->Write();
				//GenGenSignalsPtWeight[i][j]->Write();
				GenGenBG[i][j]->Write();
				//GenGenBGPtWeight[i][j]->Write();
			}
		}
		fMixEvtCorr->Close();
		
	}
	
	//fout2->Close();
	string jfc = "jffCorrApplied";
	if(!applyJFFcorrs) jfc = "";
	cprint->SaveAs(Form("%s%s_div_%s%s_bgSub%d_%s_cymbalTune_truncatePoints_finalJFFs.pdf",jet.c_str(), track.c_str(), jetDiv.c_str(), trackDiv.c_str(), doBGsub, jfc.c_str()));
	
	if(doJFFCorrs){
		string spilloverSwitch = "subeNon0";
		if(!doSpillover) spilloverSwitch = "sube0";
		TFile *fout = new TFile(Form("JFFcorrs_PythiaHydjet_%s_AllJets_finalJFFJEC_expandedRange_withFits.root",spilloverSwitch.c_str()),"recreate");
		fout->cd();
		for(int i=0; i<trackPtBins; i++){
			for(int j=0; j<nCentBins; j++){
				JFFcorrs[i][j]->Write();
				correctionEtaProj[i][j]->Write();
				correctionPhiProj[i][j]->Write();
				corrEtaProj[i][j]->Write();
				corrPhiProj[i][j]->Write();
			}
		}
		fout->Close();
	}
	if(writeOutCorrelations){
		TFile *foutFinal;
		if(isData1 || isData2) foutFinal = new TFile("TotalCorrelations_ppDataOnly_fullCymbalCorrs_withWithoutJFFs.root","recreate");
		else foutFinal = new TFile(Form("TotalCorrelations_%s%s2%s%s_GluonJets_cymbalTune_fullCymbalCorrs_finalJFFs.root",jet.c_str(),track.c_str(),jetDiv.c_str(),trackDiv.c_str()),"recreate");
		for(int i=1; i<trackPtBins-1; i++){
			for(int j=0; j<nCentBins; j++){
				GenGenSignalsSub[i][j]->Write();
				GenRecoSignalsSub[i][j]->Write();
				
				GenGenBG[i][j]->Write();
				GenRecoBG[i][j]->Write();
				GenGenEtaProj[i][j]->Write();
				GenRecoEtaProj[i][j]->Write();
				GenGenPhiProj[i][j]->Write();
				GenRecoPhiProj[i][j]->Write();
				GenRecoBGPhiProj[i][j]->Write();
			}
		}
		for(int j=0; j<nCentBins; j++){
			nJetsGen[j]->Write();
			nJetsReco[j]->Write();
		}
		foutFinal->Close();
	}
	if(writeDRHistosForPlotting){
		string drfilename = "InclJetData_bgSubtractedCorrs_dR_fullHistos_looseTrkCuts_noJFFs_xiaoTrkCorrs_CymbalTune_eta1p6.root";
		if(isBjets2) drfilename = "BJetData_bgSubtractedCorrs_dR_fullHistos_looseTrkCuts_officialJTCs_wideBinning.root";
		if(!isData1 && !isBjets1) drfilename = "PbPbMC_PbPbMCMixing_bgSubtractedCorrs_cymbalTune_fullCymbalTune_dR_QuarkVsGluonJets_fullHistos.root";
		if(!isData1 && isBjets1) drfilename = "PbPbMC_bjet_bgSubtractedCorrs_dR_fullHistos.root";
		TFile *foutdRHistos = new TFile(drfilename.c_str(),"recreate");
		foutdRHistos->cd();
		for(int i=0; i<trackPtBins; i++){
			for(int j=0; j<nCentBins; j++){
				GenRecoDRDistr[i][j]->SetName(Form("JetShape2_QuarkJets_Yield_%sjet_BkgSub%sInclusive_%s_%s_%s_%s",bjetString2.c_str(),ptWeightString.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()));
				GenRecoDRDistr[i][j]->Write();
				//HJFFCorrsToApply[i][j]->Write();
				
				GenGenDRDistr[i][j]->SetName(Form("JetShape2_GluonJets_Yield_%sjet_BkgSub%sInclusive_%s_%s_%s_%s",bjetString2.c_str(),ptWeightString.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()));
				GenGenDRDistr[i][j]->Write();
				
				GenRecoEtaProj[i][j]->SetName(Form("JetEta_QuarkJets_Yield_%sjet_BkgSub%sInclusive_%s_%s_%s_%s",bjetString2.c_str(),ptWeightString.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()));
				GenRecoEtaProj[i][j]->Write();
				
				GenGenEtaProj[i][j]->SetName(Form("JetEta_GluonJets_Yield_%sjet_BkgSub%sInclusive_%s_%s_%s_%s",bjetString2.c_str(),ptWeightString.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()));
				GenGenEtaProj[i][j]->Write();
				
				JFFDRProj[i][j]->SetName(Form("JFFCorr_%sjet_BkgSub%sInclusive_%s_%s_%s_%s",bjetString2.c_str(),ptWeightString.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()));
				JFFDRProj[i][j]->Write();
			}
			//GenGenDRDistr[i][0]->SetName(Form("JetShape2_Yield_%sjet_BkgSub%sInclusive_pp_%s_%s",bjetString1.c_str(),ptWeightString.c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()));
			//GenGenDRDistr[i][0]->Write();
			
			//GenGenEtaProj[i][0]->SetName(Form("JetEta_Yield_%sjet_BkgSub%sInclusive_pp_%s_%s",bjetString1.c_str(),ptWeightString.c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()));
			//GenGenEtaProj[i][0]->Write();
		}
		foutdRHistos->Close();
	}
}