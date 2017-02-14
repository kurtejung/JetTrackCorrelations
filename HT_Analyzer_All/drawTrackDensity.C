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
	TH1D *projEta = (TH1D*)background->ProjectionX("projEta",1,background->GetNbinsX(),"e");
	for(int ixbin=1; ixbin<=background->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=background->GetNbinsY(); iybin++){
			background->SetBinContent(ixbin,iybin, projEta->GetBinContent(ixbin));
			background->SetBinError(ixbin,iybin, projEta->GetBinError(ixbin)*sqrt(background->GetNbinsY()));
		}
	}
	TF1 *littleFit = new TF1("littleFit","pol0",-0.1,0.1);
	projEta->Fit(littleFit,"qN0","",-0.1,0.1);
	double normFactor = littleFit->GetParameter(0);
	background->Scale(1./normFactor);//projEta->GetBinContent(projEta->FindBin(0.)));
	projEta->Scale(1./normFactor);//projEta->GetBinContent(projEta->FindBin(0)));
	
	//Do mixed event correction
	for(int ixbin=1; ixbin<=signal->GetNbinsX(); ixbin++){
		for(int iybin=1; iybin<=signal->GetNbinsY(); iybin++){
			if(projEta->GetBinContent(ixbin)) signal->SetBinContent(ixbin,iybin, signal->GetBinContent(ixbin,iybin)/projEta->GetBinContent(ixbin));
			else signal->SetBinContent(ixbin,iybin,0);
			if(projEta->GetBinContent(ixbin)) signal->SetBinError(ixbin,iybin, ReturnDivError(signal->GetBinContent(ixbin,iybin), signal->GetBinError(ixbin,iybin), projEta->GetBinContent(ixbin),projEta->GetBinError(ixbin)));
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
				if(inclPhiBG->GetBinContent(iybin)) signal->SetBinError(ixbin,iybin, ReturnSubError(signal->GetBinError(ixbin,iybin), inclPhiBG->GetBinError(iybin)));
			}
		}
	}
	
	delete inclPhiFG;
	delete inclPhiBG;
	delete projEta;
	
	return signal;
}

void drawTrackDensity(string jet="Reco", string track="Reco", string jetDiv="Reco", string trackDiv="Reco", bool doBGsub=true){
	
	bool doJFFCorrs = false;
	bool applyJFFcorrs = false;
	bool writeOutCorrelations = false;
	bool drawUncertainties = false;
	bool writeDRHistosForPlotting = true;
	bool doptWeight = false;
	
	bool isData1 = true;
	bool isData2 = true;
	
	string trkCorr = "_notrkcorr";
	string trkCorrDiv = "_notrkcorr";
	if(track=="Reco") trkCorr = "";
	if(trackDiv=="Reco") trkCorrDiv = "";
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
	
	const double xdrbins[20] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.,1.2,1.4,1.6,1.8,2.};
	
	const int trackPtBins = 10;
	const int nCentBins = 4;
	
	const double xTrkBinDouble[trackPtBins+1] = {0.5,0.7,1.,2.,3.,4.,8.,12.,16.,20.,999.};
	const float mean_pts[trackPtBins] = {0.844,1.35,2.35,3.37,5.07,9.72,13.8,17.9,41.};
	const double ptWidth[trackPtBins] = {0.2,0.3,1,1,1,4,4,4,280};
	const string xTrkBins[trackPtBins+1] = {"TrkPt0p5", "TrkPt0p7", "TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt999" };
	const string halliexTrkBins[trackPtBins+1] = {"TrkPt05", "TrkPt07", "TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt300" };
	const string xCentBins[nCentBins+1] = {"Cent0", "Cent10", "Cent30","Cent50", "Cent100"};
	
	//PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_Merged_refpt_newJetTrackCorrections_fineBin.root
	TFile *fJFFs = new TFile("JFFcorrs_PythHydjet_sube0_newJFFJEC_withFits.root");
	TFile *fJFFsPythia = new TFile("JFFcorrs_PythiPP_newJFFJEC_withFits.root");
	TFile *fSpillovers = new TFile("JFFcorrs_subeNon0_newJFFs.root");
	TFile *fHSpillovers = new TFile("Inclusive_Hydjet_SpillOvers_pTweighted.root");
	TFile *fHJFFs = new TFile("Inclusive_Hydjet_Hallie_JFFResiduals.root");
	
	TFile *fJFFsWeighted = new TFile("JFFcorrs_PythHydjet_pTweighted_sube0_newJFFJEC_withFits.root");
	TFile *fJFFsPythiaWeighted = new TFile("JFFcorrs_PythiaOnly_pTweighted_sube0_newJFFJEC_withFits.root");
	
	//TFile *fSpillovers = new TFile("JFFcorrs_subeNon0_TestOldJFFs_withFits.root");
	TFile *fGenGen = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_Merged_refpt_newJetTrackCorrections_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenReco = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_Merged_refpt_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenGenSube0 = new TFile(Form("root_output/PbPb_5TeV_Pythia6Hydjet_noMix_%s%sReduced_Merged_sube0_wPtWeightHistos.root",jet.c_str(),track.c_str()));
	TFile *fGenGenSubeNon0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_subeNon0_Merged_withMix_newJetTrackCorrections_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoSube0 = new TFile(Form("root_output/PbPb_5TeV_MC_noMix_newJetTrackCorrections_fineBin_fixTrkCorrs_%s%sSube0.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenRecoSubeNon0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_subeNon0_Merged_withMix_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenRecoQuarkSube0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_QuarkJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenRecoGluonSube0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_GluonJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenGenQuarkSube0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_QuarkJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenGenGluonSube0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_GluonJet_sube0_Merged_noMix_newJetTrackCorrections_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoQuarkSubeNon0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_QuarkJet_subeNon0_Merged_withMix_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenRecoGluonSubeNon0 = new TFile(Form("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_%s%sReduced_GluonJet_subeNon0_Merged_withMix_newJetTrackCorrections_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenGenQuark = new TFile(Form("root_output/PbPb_5TeV_MC_Pythia6Hydjet_noMix_Merged_%s%sQuark_newJFFJECs.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoQuark = new TFile(Form("root_output/PbPb_5TeV_MC_Pythia6Hydjet_noMix_Merged_%s%sQuark_newJFFJECs.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *fGenGenGluon = new TFile(Form("root_output/PbPb_5TeV_MC_Pythia6Hydjet_noMix_Merged_%s%sGluon_newJFFJECs.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoGluon = new TFile(Form("root_output/PbPb_5TeV_MC_Pythia6Hydjet_noMix_Merged_%s%sGluon_newJFFJECs.root",jetDiv.c_str(),trackDiv.c_str()));
	
	TFile *fGenGenPythia = new TFile(Form("root_output/pp_5TeV_MC_Pythia_noMix_Merged_%s%sReduced_fixVz_fixTrkCorr_fineBin.root",jet.c_str(),track.c_str()));
	TFile *fGenRecoPythia = new TFile(Form("root_output/pp_5TeV_MC_Pythia_noMix_Merged_%s%sReduced_fixVz_fixTrkCorr_fineBin.root",jetDiv.c_str(),trackDiv.c_str()));
	TFile *extra = new TFile("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_GenRecoReduced_Merged_noMix_newJetTrackCorrections_pt16to20TrkCorrFix_fineBin.root");
	
	TFile *hallieGenGen = new TFile("root_output/Pythia_GenJet_GenTrack_Aug23.root");
	TFile *hallieGenReco = new TFile("root_output/HallieFiles/HydJet_GenJet_RecoTrack_Aug23.root");
	
	TFile *ppData = new TFile("root_output/ppData_5TeV_withMix_bJets_XiaoJetTrackCorrections_Merged_fineBin.root");
	TFile *PbPbData = new TFile("root_output/PbPbData_5TeV_withMix_bJets_XiaoJetTrackCorrections_Merged_fineBin.root");
	
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
	TH1D *HSpilloverEtaProj[trackPtBins][nCentBins];
	TH1D *JFFEtaProj[trackPtBins][nCentBins];
	
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
			
			string toGet = Form("JFFcorrs_cent%d_pt%d",j,i);
			JFFcorrsToApply[i][j] = (TH2D*)fJFFs->Get(toGet.c_str())->Clone(toGet.c_str());
			JFFcorrsToApplyWeighted[i][j] = (TH2D*)fJFFsWeighted->Get(toGet.c_str())->Clone(toGet.c_str());
			
			JFFcorrsToApplyPyth[i][j] = (TH2D*)fJFFsPythia->Get(toGet.c_str())->Clone(toGet.c_str());
			JFFcorrsToApplyPythWeighted[i][j] = (TH2D*)fJFFsPythiaWeighted->Get(toGet.c_str())->Clone(toGet.c_str());
			
			//spilloverCorrsToApply[i][j] = new TH2D(Form("Closure_2D_%s_%s_Pt100_Pt300_%s_%s",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);
			HspilloverCorrsToApply[i][j] = new TH2D(Form("Hallie_Closure_2D_%s_%s_Pt100_Pt300_%s_%s",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);
			HJFFCorrsToApply[i][j] = new TH2D(Form("Hallie_JFFClosure_2D_%s_%s_Pt100_Pt300_%s_%s",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str()),"",500,-5,5,200,-TMath::Pi()/2,3*TMath::Pi()/2);
			if(i>0){
				toGet = Form("Eta_SpillOver_Fit_%s_%s_Pt100_Pt300_%s_%s",xCentBins[j].c_str(),xCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
				TF1 *etaFitIn2 = (TF1*)fHSpillovers->Get(toGet.c_str())->Clone(Form("HetaFit_%d_%d",i,j));
				toGet = Form("Phi_SpillOver_Fit_%s_%s_Pt100_Pt300_%s_%s",xCentBins[j].c_str(),xCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
				TF1 *phiFitIn2 = (TF1*)fHSpillovers->Get(toGet.c_str())->Clone(Form("HphiFit_%d_%d",i,j));
				
				toGet = Form("JFF_Residual_Eta_%s_%s_Pt100_Pt300_%s_%s",xCentBins[j].c_str(),xCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
				TH1D *etaFitIn3 = (TH1D*)fHJFFs->Get(toGet.c_str())->Clone(Form("JFFetaFit_%d_%d",i,j));
				toGet = Form("JFF_Residual_Phi_%s_%s_Pt100_Pt300_%s_%s",xCentBins[j].c_str(),xCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
				TH1D *phiFitIn3 = (TH1D*)fHJFFs->Get(toGet.c_str())->Clone(Form("JFFphiFit_%d_%d",i,j));
				
				for(int ixbin=1; ixbin<HJFFCorrsToApply[i][j]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<HJFFCorrsToApply[i][j]->GetNbinsY(); iybin++){
						double content = (etaFitIn3->GetBinContent(etaFitIn3->FindBin(HJFFCorrsToApply[i][j]->GetXaxis()->GetBinCenter(ixbin))) + 
						phiFitIn3->GetBinContent(phiFitIn3->FindBin(HJFFCorrsToApply[i][j]->GetYaxis()->GetBinCenter(iybin))))/2.;
						HJFFCorrsToApply[i][j]->SetBinContent(ixbin,iybin,content);
						HJFFCorrsToApply[i][j]->SetBinError(ixbin,iybin,1e-12);
					}
				}
				//refillFits(etaFitIn, phiFitIn, spilloverCorrsToApply[i][j],1);
				refillFits(etaFitIn2, phiFitIn2, HspilloverCorrsToApply[i][j],0);
				HspilloverCorrsToApply[i][j]->Scale(ptWidth[i]);
				if(i>4) HspilloverCorrsToApply[i][j]->Reset();
			}
			toGet = Form("JFFcorrs_cent%d_pt%d",j,i);
			spilloverCorrsToApply[i][j] = (TH2D*)fSpillovers->Get(toGet.c_str())->Clone(Form("spillover_%d_%d",i,j));
			if(i>4) spilloverCorrsToApply[i][j]->Reset();
		      //JFFcorrsToApply[i]->Scale(0.5);
		      
			if(drawUncertainties){
				toGet = Form("stackedRelatives_pt%d_cent%d",i,j);
				uncertainties[i][j] = (TH1D*)syst->Get(toGet.c_str())->Clone(toGet.c_str());
				unc[i][j] = new TGraphErrors(uncertainties[i][j]);
			}
			
			if(isData1) toGet = Form("Data_hbJetTrackSignalBackground%s%s%s_%s_Pt100_Pt300_%s_%s",ptWeightString.c_str(), trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			else toGet = Form("%sJet_%sTrack_hbJetTrackSignalBackground%s%s%s_%s_Pt100_Pt300_%s_%s",jet.c_str(),track.c_str(),ptWeightString.c_str(), trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			cout << "getting " << toGet << endl;
			GenGenSignals[i][j] = (TH2D*)ppData->Get(toGet.c_str())->Clone(toGet.c_str());
			
			if(isData2) toGet = Form("Data_hbJetTrackSignalBackground%s%s%s_%s_Pt100_Pt300_%s_%s",ptWeightString.c_str(), trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			else toGet = Form("%sJet_%sTrack_hbJetTrackSignalBackground%s%s%s_%s_Pt100_Pt300_%s_%s",jetDiv.c_str(),trackDiv.c_str(),ptWeightString.c_str(),trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str());
			cout << "getting " << toGet << endl;
			GenRecoSignals[i][j] = (TH2D*)PbPbData->Get(toGet.c_str())->Clone(toGet.c_str());
			
			if(isData1) toGet = Form("Data_hbJetTrackME%s%s_%s_Pt100_Pt300_%s_%s",trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			else toGet = Form("%sJet_%sTrack_hbJetTrackME%s%s_%s_Pt100_Pt300_%s_%s",jet.c_str(),track.c_str(),trkCorr.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			cout << "getting " << toGet << endl;
			GenGenBG[i][j] = (TH2D*)ppData->Get(toGet.c_str())->Clone(toGet.c_str());
			
			if(isData2) toGet = Form("Data_hbJetTrackME%s%s_%s_Pt100_Pt300_%s_%s",trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			else toGet = Form("%sJet_%sTrack_hbJetTrackME%s%s_%s_Pt100_Pt300_%s_%s",jetDiv.c_str(),trackDiv.c_str(),trkCorrDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBins[i].c_str(),xTrkBins[i+1].c_str());
			cout << "getting " << toGet << endl;
			GenRecoBG[i][j] = (TH2D*)PbPbData->Get(toGet.c_str())->Clone(toGet.c_str());
			
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
		}
	}
	
	for(int j=0; j<nCentBins; j++){
		string toGet;
		if(isData1) toGet = Form("Data_all_bjets_corrpT%s_%s_Pt100_Pt300",xCentBins[j].c_str(),xCentBins[j+1].c_str());
		else toGet = Form("%sJet_%sTrack_all_bjets_corrpT%s_%s_Pt100_Pt300",jet.c_str(),track.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str());
		nJetsGen[j] = (TH1F*)ppData->Get(toGet.c_str())->Clone(toGet.c_str());
		
		if(isData2) toGet = Form("Data_all_bjets_corrpT%s_%s_Pt100_Pt300",xCentBins[j].c_str(),xCentBins[j+1].c_str());
		else toGet = Form("%sJet_%sTrack_all_bjets_corrpT%s_%s_Pt100_Pt300",jetDiv.c_str(),trackDiv.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str());
		nJetsReco[j] = (TH1F*)PbPbData->Get(toGet.c_str())->Clone(toGet.c_str());
		
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
					if(!doptWeight){
					GenGenSignalsSub[i][j]->Add(JFFcorrsToApplyPyth[i][j],-1);
					//GenGenSignalsSub[i][j]->Add(HspilloverCorrsToApply[i][j],-1);
					}
					else{
					GenGenSignalsSub[i][j]->Add(JFFcorrsToApplyPythWeighted[i][j],-1);
					//GenGenSignalsSub[i][j]->Add(HspilloverCorrsToApply[i][j],-1);
					}
				}
				if(jetDiv=="Reco"){
					if(!doptWeight){
						GenRecoSignalsSub[i][j]->Add(JFFcorrsToApply[i][j],-1);
						//GenGenSignalsSub[i][j]->Add(HspilloverCorrsToApply[i][j],-1);	
					}
					else{
						GenRecoSignalsSub[i][j]->Add(JFFcorrsToApplyWeighted[i][j],-1);
						//GenGenSignalsSub[i][j]->Add(HspilloverCorrsToApply[i][j],-1);
					}
				}
			}
			
			//cout << "track pt" << i << " cent " << j << " gen counts: "<< GenGenSignalsSub[i][j]->GetEntries() << endl;
			//cout << "track pt" << i << " cent " << j << " counts: "<< GenRecoSignalsSub[i][j]->GetEntries() << endl;
			
			int lowBin = GenGenSignalsSub[i][j]->GetYaxis()->FindBin(-1);
			int highBin = GenGenSignalsSub[i][j]->GetYaxis()->FindBin(1);
			GenGenEtaProj[i][j] = GenGenSignalsSub[i][j]->ProjectionX(Form("GenGenEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			//GenGenEtaProj[i][j]->Scale(1./GenGenSignalsSub[i][j]->GetYaxis()->GetBinWidth(2));
			GenGenEtaProj[i][j]->Rebin(5);
			GenGenEtaProj[i][j]->Scale(1./5.);
			GenRecoEtaProj[i][j] = GenRecoSignalsSub[i][j]->ProjectionX(Form("GenRecoEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			//GenRecoEtaProj[i][j]->Scale(1./GenRecoSignalsSub[i][j]->GetYaxis()->GetBinWidth(2));
			GenRecoEtaProj[i][j]->Rebin(5);
			GenRecoEtaProj[i][j]->Scale(1./5.);
			
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
			
			lowBin = spilloverCorrsToApply[i][j]->GetYaxis()->FindBin(-1.5);
			highBin = spilloverCorrsToApply[i][j]->GetYaxis()->FindBin(1.5);
			SpilloverEtaProj[i][j] = spilloverCorrsToApply[i][j]->ProjectionX(Form("SpilloverEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			SpilloverEtaProj[i][j]->Scale(1./spilloverCorrsToApply[i][j]->GetYaxis()->GetBinWidth(2));
			SpilloverEtaProj[i][j]->Rebin(5);
			SpilloverEtaProj[i][j]->Scale(1./5.);
			HSpilloverEtaProj[i][j] = HspilloverCorrsToApply[i][j]->ProjectionX(Form("HSpilloverEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			HSpilloverEtaProj[i][j]->Scale(1./HspilloverCorrsToApply[i][j]->GetYaxis()->GetBinWidth(2));
			HSpilloverEtaProj[i][j]->Rebin(5);
			HSpilloverEtaProj[i][j]->Scale(1./5.);
			JFFEtaProj[i][j] = JFFcorrsToApply[i][j]->ProjectionX(Form("JFFEtaProj_%d_%d",i,j),lowBin, highBin,"e");
			JFFEtaProj[i][j]->Scale(1./JFFcorrsToApply[i][j]->GetYaxis()->GetBinWidth(2));
			JFFEtaProj[i][j]->Rebin(5);
			JFFEtaProj[i][j]->Scale(1./5.);
			
			//cout << "at " << i << " " << j << " spillover proj int: "<< SpilloverEtaProj[i][j]->Integral() << endl;
			//cout << " jff proj int: "<< JFFEtaProj[i][j]->Integral() << endl;
						
			if(doJFFCorrs){

				JFFcorrs[i][j] = (TH2D*)GenRecoSignalsSub[i][j]->Clone(Form("JFFcorrs_cent%d_pt%d",j,i));
				JFFcorrs[i][j]->Add(GenGenSignalsSub[i][j], -1);

				//now symmetrize the corrections
				for(int ixbin=1; ixbin<=JFFcorrs[i][j]->GetNbinsX(); ixbin++){
					for(int iybin=1; iybin<=JFFcorrs[i][j]->GetNbinsY(); iybin++){
						double xcent = JFFcorrs[i][j]->GetXaxis()->GetBinCenter(ixbin);
						double ycent = JFFcorrs[i][j]->GetYaxis()->GetBinCenter(iybin);
						int negX = JFFcorrs[i][j]->FindBin(-1*xcent,ycent);
						//int negY = JFFcorrs[i][j]->GetYaxis()->FindBin(-1*ycent);
						
						double dr = sqrt(pow(xcent,2)+pow(ycent,2));
						if(dr>1.5 || (dr>1.5 && i>4)){
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
				correctionEtaProj[i][j] = JFFcorrs[i][j]->ProjectionX(Form("correctionEtaProj%d_%d",i,j),0, -1,"e");
				correctionPhiProj[i][j] = JFFcorrs[i][j]->ProjectionY(Form("correctionPhiProj%d_%d",i,j),0, -1,"e");
				correctionEtaProj[i][j]->Rebin(5);
				correctionPhiProj[i][j]->Rebin(5);
				correctionEtaProj[i][j]->Scale(1./(5.*(double)JFFcorrs[i][j]->GetYaxis()->GetBinWidth(2)));
				correctionPhiProj[i][j]->Scale(1./(5.*(double)JFFcorrs[i][j]->GetYaxis()->GetBinWidth(2)));
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
						
			for(int ixbin=1; ixbin<=GenGenSignalsSub[i][j]->GetNbinsX(); ixbin++){
				for(int iybin=1; iybin<=GenGenSignalsSub[i][j]->GetNbinsY(); iybin++){
					double xcent = GenGenSignalsSub[i][j]->GetXaxis()->GetBinCenter(ixbin);
					double ycent = GenGenSignalsSub[i][j]->GetYaxis()->GetBinCenter(iybin);
					double content = GenGenSignalsSub[i][j]->GetBinContent(GenGenSignalsSub[i][j]->GetBin(ixbin,iybin));
					//cout << "gen gen bin " << ixbin << ", " << iybin << " content " << content << endl;
					double dr = sqrt(pow(xcent,2)+pow(ycent,2));
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
					
				}
			}
			//GenGenDRDistr[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			//GenRecoDRDistr[i][j]->Divide(GenRecoDRDistrNorm[i][j]);
			
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
	
	TCanvas *ceta = new TCanvas("ceta","",1200,1600);
	TLatex *l2[trackPtBins][nCentBins];
	ceta->Divide(nCentBins,trackPtBins-2);
	TLine *lineAt1 = new TLine(-2.5,1,2.5,1);
	lineAt1->SetLineStyle(2);
	for(int j=0; j<nCentBins; j++){
		for(int i=1; i<trackPtBins-1; i++){
			ceta->cd(i*nCentBins+(3-j)+1-4);
			GenGenEtaProj[i][j]->GetYaxis()->SetNdivisions(505);
			GenGenEtaProj[i][j]->GetYaxis()->SetLabelSize(0.15);
			//GenGenEtaProj[i][j]->SetMinimum(0.5);//GenRecoEtaProj[i][0]->GetMinimum()*1.1);
			//GenGenEtaProj[i][j]->SetMaximum(1.5);//GenRecoEtaProj[i][0]->GetMaximum()*1.2);
			GenGenEtaProj[i][j]->GetXaxis()->SetLabelSize(0.15);
			GenGenEtaProj[i][j]->GetXaxis()->SetRangeUser(-2.5,2.5);
			//GenRecoEtaProj[i][j]->GetYaxis()->SetRangeUser(0.5,1.5);
			GenRecoEtaProj[i][j]->SetLineColor(2);
			//GenGenEtaProj[i][j]->Divide(GenRecoEtaProj[i][j]);
			GenGenEtaProj[i][j]->Draw("");
			
			GenRecoEtaProj[i][j]->Draw("same");
			SpilloverEtaProj[i][j]->SetLineColor(6);
			SpilloverEtaProj[i][j]->SetMarkerColor(6);
			SpilloverEtaProj[i][j]->SetMarkerSize(0.2);
			//SpilloverEtaProj[i][j]->Draw("same");
			HSpilloverEtaProj[i][j]->SetLineColor(kGreen+1);
			HSpilloverEtaProj[i][j]->SetMarkerColor(kGreen+1);
			HSpilloverEtaProj[i][j]->SetMarkerSize(0.2);
			//HSpilloverEtaProj[i][j]->Draw("same");
			JFFEtaProj[i][j]->SetLineColor(kOrange);
			//JFFEtaProj[i][j]->Draw("same");
			l2[i][j] = new TLatex(-2,GenRecoEtaProj[i][j]->GetMaximum()*0.97,Form("#Delta#eta, %s-%s, %g<pT<%g",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBinDouble[i],xTrkBinDouble[i+1]));
			l2[i][j]->SetTextSize(0.15);
			l2[i][j]->Draw("same");
			lineAt1->Draw("Same");
			
			int lb = GenGenEtaProj[i][j]->FindBin(-0.4);
			int hb = GenGenEtaProj[i][j]->FindBin(0.4);
			double integralGenGen = GenGenEtaProj[i][j]->Integral(lb,hb);
			double integralGenReco = GenRecoEtaProj[i][j]->Integral(lb,hb);
			
			//cout << "cent bin " << j << " pt: "<< i << " Scale: " << integralGenGen/integralGenReco << endl;
			
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
			l3[i][j] = new TLatex(-0.5,GenRecoPhiProj[i][j]->GetMaximum()*0.97,Form("#Delta#phi, %s-%s, %g<pT<%g",xCentBins[j].c_str(),xCentBins[j+1].c_str(),xTrkBinDouble[i],xTrkBinDouble[i+1]));
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
	TCanvas *cprint = new TCanvas("cprint","",1200/3.,1600);
	cprint->Divide(1, trackPtBins-1);
	for(int j=0; j<1; j++){
		for(int i=1; i<trackPtBins-1; i++){
			cprint->cd(i*1+(3-j)+1-4);
			
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
			
			
			//if(i<5) GenGenDRDistr[i][j]->Divide(GenGenDRDistr[i][j],GenRecoDRDistr[i][j],1,1);
			//else GenGenDRDistr[i][j]->Divide(GenGenDRDistr[i][j],GenRecoDRDistr[i][j],1,1);
			
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
			l1[i][j] = new TLatex(0.2,GenGenDRDistr[i][j]->GetMaximum()*0.85,Form(/*"%s-%s,*/ "%g<pT<%g",/*xCentBins[j].c_str(),xCentBins[j+1].c_str(),*/xTrkBinDouble[i],xTrkBinDouble[i+1]));
			l1[i][j]->SetTextSize(0.15);
			l1[i][j]->Draw("same");
			GenRecoDRDistrRed[i][j]->Draw("Same");
			
		}
	}
	fout2->Close();
	string jfc = "jffCorrApplied";
	if(!applyJFFcorrs) jfc = "";
	//cprint->SaveAs(Form("%s%s_div_%s%s_bgSub%d_%s_Pythiapp_truncatePoints_newCorrections_newUnc.pdf",jet.c_str(), track.c_str(), jetDiv.c_str(), trackDiv.c_str(), doBGsub, jfc.c_str()));
	
	if(doJFFCorrs){
		TFile *fout = new TFile("JFFcorrs_PythHydjet_sube0_newJFFJEC_withFits.root","recreate");
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
		if(isData1 || isData2) foutFinal = new TFile("TotalCorrelations_ppData_newJFFs.root","recreate");
		else foutFinal = new TFile(Form("TotalCorrelations_%s%s2%s%s_newJFFs.root",jet.c_str(),track.c_str(),jetDiv.c_str(),trackDiv.c_str()),"recreate");
		for(int i=0; i<trackPtBins; i++){
			for(int j=0; j<nCentBins; j++){
				GenGenSignalsSub[i][j]->Write();
				//GenRecoSignalsSub[i][j]->Write();
			}
		}
		foutFinal->Close();
	}
	if(writeDRHistosForPlotting){
		TFile *foutdRHistos = new TFile("bJetData_bgSubtractedCorrs_dR_fullHistos.root","recreate");
		foutdRHistos->cd();
		for(int i=0; i<trackPtBins; i++){
			for(int j=0; j<nCentBins; j++){
				GenRecoDRDistr[i][j]->SetName(Form("JetShape2_Yield_bjet_BkgSub%sInclusive_%s_%s_%s_%s",ptWeightString.c_str(),xCentBins[j].c_str(),xCentBins[j+1].c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str()));
				GenRecoDRDistr[i][j]->Write();
				HJFFCorrsToApply[i][j]->Write();
			}
			GenRecoDRDistr[i][0]->SetName(Form("JetShape2_Yield_bjet_BkgSub%sInclusive_HallieppMC_%s_%s",ptWeightString.c_str(),halliexTrkBins[i].c_str(),halliexTrkBins[i+1].c_str()));
			GenRecoDRDistr[i][0]->Write();
		}
		foutdRHistos->Close();
	}
}