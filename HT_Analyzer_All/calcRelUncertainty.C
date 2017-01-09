
#include "TH2D.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include <iostream>
	
using namespace std;

const double xdrbins[15] = {0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.};

const double fatdrbins[10] = {0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.};

	
const int trackPtBins = 10;
const int nCentBins = 4;

const double ptWidth[trackPtBins] = { 0.2, 0.3, 1., 1., 1., 4., 4., 4., 4., 80.};
	
const double xTrkBinDouble[trackPtBins+1] = {0.5,0.7,1.,2.,3.,4.,8.,12.,16.,20.,999.};
const string xTrkBins[trackPtBins+1] = {"TrkPt0p5","TrkPt0p7", "TrkPt1", "TrkPt2", "TrkPt3", "TrkPt4", "TrkPt8", "TrkPt12", "TrkPt16", "TrkPt20", "TrkPt999" };
const int xCentBins[nCentBins+1] = {0,10,30,50,100};
const float mean_pts[trackPtBins] = {0.5,0.844,1.35,2.35,3.37,5.07,9.72,13.8,17.9,22.};

void takeTH1Abs(TH1D* input){
	for(int i=1; i<=input->GetNbinsX(); i++){
		input->SetBinContent(i, abs(input->GetBinContent(i)));
	}
	
}
	
void calcRelUncertainty(){
		
	TFile *fin = new TFile("TotalCorrelations.root");
	TFile *fSpillovers = new TFile("JFFcorrs_spillover.root");
	TFile *fResidualJES = new TFile("JFFcorrs_sube0.root");
	
	TH2D *GenGenSignal[nCentBins][trackPtBins];
	TH2D *RecoRecoSignal[nCentBins][trackPtBins];
	TH2D *hSpillovers[nCentBins][trackPtBins];
	TH2D *hResidualJES[nCentBins][trackPtBins];
	
	TH2D *pairAccept[nCentBins][trackPtBins];
	TH2D *bgSub[nCentBins][trackPtBins];
	
	TH2D *tracking[nCentBins][trackPtBins];
	TH2D *trackingRes[nCentBins][trackPtBins];
	TH2D *otherJFF[nCentBins][trackPtBins];
	
	TH1D *etaProj[nCentBins][trackPtBins];
	TH1D *etaProjRebin[nCentBins][trackPtBins];
	TH1D *asymm[nCentBins][trackPtBins];
	
	TH1D *absoluteErrSpill[nCentBins][trackPtBins];
	TH1D *absoluteErrResidJES[nCentBins][trackPtBins];
	TH1D *absoluteErrPair[nCentBins][trackPtBins];
	TH1D *absoluteErrBG[nCentBins][trackPtBins];
	TH1D *absoluteErrTrk[nCentBins][trackPtBins];
	TH1D *absoluteErrTrkRes[nCentBins][trackPtBins];
	TH1D *absoluteErrOther[nCentBins][trackPtBins];
	
	TH1D *relErrSpill[nCentBins][trackPtBins];
	TH1D *relErrResidJES[nCentBins][trackPtBins];
	TH1D *relErrPair[nCentBins][trackPtBins];
	TH1D *relErrBG[nCentBins][trackPtBins];
	TH1D *relErrTrk[nCentBins][trackPtBins];
	TH1D *relErrTrkRes[nCentBins][trackPtBins];
	TH1D *relErrOther[nCentBins][trackPtBins];
	
	TH1D *stackedAbsolutes[nCentBins][trackPtBins];
	TH1D *stackedRelatives[nCentBins][trackPtBins];
	
	TH1D *GenGenDRDistrNorm[nCentBins][trackPtBins];
	TH1D *GenGenDRDistrNormWeight[nCentBins][trackPtBins];
	
	for(int j=0; j<trackPtBins; j++){
		for(int i=0; i<nCentBins; i++){
			absoluteErrSpill[i][j] = new TH1D(Form("absoluteErrSpill_pt%d_cent%d",j,i),"",14,xdrbins);
			absoluteErrResidJES[i][j]= new TH1D(Form("absoluteErrResidJES_pt%d_cent%d",j,i),"",14,xdrbins);
			absoluteErrPair[i][j] = new TH1D(Form("absoluteErrPair_pt%d_cent%d",j,i),"",14,xdrbins);
			absoluteErrBG[i][j] = new TH1D(Form("absoluteErrBG_pt%d_cent%d",j,i),"",14,xdrbins);
			absoluteErrTrk[i][j] = new TH1D(Form("absoluteErrTrk_pt%d_cent%d",j,i),"",14,xdrbins);
			absoluteErrTrkRes[i][j] = new TH1D(Form("absoluteErrTrkRes_pt%d_cent%d",j,i),"",14,xdrbins);
			absoluteErrOther[i][j] = new TH1D(Form("absoluteErrOther_pt%d_cent%d",j,i),"",14,xdrbins);
			
			relErrSpill[i][j] = new TH1D(Form("relErrSpill_pt%d_cent%d",j,i),"",14,xdrbins);
			relErrResidJES[i][j] = new TH1D(Form("relErrResidJES_pt%d_cent%d",j,i),"",14,xdrbins);
			relErrPair[i][j] = new TH1D(Form("relErrPair_pt%d_cent%d",j,i),"",14,xdrbins);
			relErrBG[i][j] = new TH1D(Form("relErrBG_pt%d_cent%d",j,i),"",14,xdrbins);
			relErrTrk[i][j] = new TH1D(Form("relErrTrk_pt%d_cent%d",j,i),"",14,xdrbins);
			relErrTrkRes[i][j] = new TH1D(Form("relErrTrkRes_pt%d_cent%d",j,i),"",14,xdrbins);
			relErrOther[i][j] = new TH1D(Form("relErrOther_pt%d_cent%d",j,i),"",14,xdrbins);
			
			stackedAbsolutes[i][j] = new TH1D(Form("stackedAbsolutes_pt%d_cent%d",j,i),"",14,xdrbins);
			stackedRelatives[i][j] = new TH1D(Form("stackedRelatives_pt%d_cent%d",j,i),"",14,xdrbins);
			
			GenGenDRDistrNorm[i][j] = new TH1D(Form("GenGenDRDistrNorm_pt%d_cent%d",j,i),"",14,xdrbins);
			GenGenDRDistrNormWeight[i][j] = new TH1D(Form("GenGenDRDistrNormWeight_pt%d_cent%d",j,i),"",14,xdrbins);
		}
	}
	
	for(int j=0; j<trackPtBins; j++){
		for(int i=0; i<nCentBins; i++){
						
			GenGenSignal[i][j] = (TH2D*)fin->Get(Form("GenJet_GenTrack_hJetTrackSignalBackground_notrkcorrCent%d_Cent%d_Pt100_Pt300_%s_%s",xCentBins[i],xCentBins[i+1],xTrkBins[j].c_str(),xTrkBins[j+1].c_str()))->Clone(Form("GenGenSignal_%d_%d",i,j));
			
			RecoRecoSignal[i][j] = (TH2D*)fin->Get(Form("RecoJet_GenTrack_hJetTrackSignalBackground_notrkcorrCent%d_Cent%d_Pt100_Pt300_%s_%s",xCentBins[i],xCentBins[i+1],xTrkBins[j].c_str(),xTrkBins[j+1].c_str()))->Clone(Form("GenGenSignal_%d_%d",i,j));
						
			//take only 1/2 the jff and spillovers for systematics
			hSpillovers[i][j] = (TH2D*)fSpillovers->Get(Form("JFFcorrs_cent%d_pt%d",i,j))->Clone(Form("spillovers_cent%d_pt%d",i,j));
			hSpillovers[i][j]->Scale(0.5);
			hSpillovers[i][j]->Scale(1./ptWidth[j]);
			if(j>4) hSpillovers[i][j]->Reset();
						
			//Only take 50-100% for residual JES uncertainties - they're independent of centrality
			hResidualJES[i][j] = (TH2D*)fResidualJES->Get(Form("JFFcorrs_cent3_pt%d",j))->Clone(Form("jffResid_cent%d_pt%d",i,j));
			hResidualJES[i][j]->Scale(0.5);
			hResidualJES[i][j]->Scale(1./ptWidth[j]);
						
			//Calculate maximum background deviation...
			int lowbin = RecoRecoSignal[i][j]->GetYaxis()->FindBin(-1);
			int highbin = RecoRecoSignal[i][j]->GetYaxis()->FindBin(1);
			etaProj[i][j] = RecoRecoSignal[i][j]->ProjectionX(Form("GenGenSignal_%d_%d_pfx",j,i), lowbin, highbin, "e");
			etaProj[i][j]->Scale(1./RecoRecoSignal[i][j]->GetXaxis()->GetBinWidth(2));
			etaProj[i][j]->Scale(0.5);
			
			int scaleFac = 0.2/etaProj[i][j]->GetBinWidth(2);
			etaProjRebin[i][j] = (TH1D*)etaProj[i][j]->Clone(Form("etaProjRebin_%d_%d",j,i));
			etaProjRebin[i][j]->Rebin(scaleFac);
			etaProjRebin[i][j]->Scale(1./(double)scaleFac);
			double maxDev=0;
			for(int ibin=1; ibin<etaProjRebin[i][j]->GetNbinsX(); ibin++){
				if(abs(etaProjRebin[i][j]->GetBinCenter(ibin))>1.5 && abs(etaProjRebin[i][j]->GetBinCenter(ibin))<2.5){
					if(abs(etaProjRebin[i][j]->GetBinContent(ibin)) > maxDev) maxDev = abs(etaProjRebin[i][j]->GetBinContent(ibin));
					if(i==3 && j==1) cout << " at location " << etaProjRebin[i][j]->GetBinCenter(ibin) << " content: "<< abs(etaProjRebin[i][j]->GetBinContent(ibin)) << " maxdev: "<< maxDev << endl;
				}
			}
			cout << "bin " << i << " " << j << " maxdev: "<< maxDev << endl;
			
			bgSub[i][j] = (TH2D*)RecoRecoSignal[i][j]->Clone(Form("bgSub_pt%d_cent%d",j,i));
			for(int ixbin=1; ixbin<=bgSub[i][j]->GetNbinsX(); ixbin++){
				for(int iybin=1; iybin<=bgSub[i][j]->GetNbinsY(); iybin++){
					bgSub[i][j]->SetBinContent(ixbin,iybin,maxDev);
					bgSub[i][j]->SetBinError(ixbin,iybin,0);
				}
			}
			bgSub[i][j]->Scale(1./(double)bgSub[i][j]->GetNbinsY());
			
			
			//calculate left-right asymmetry
			asymm[i][j] = (TH1D*)etaProj[i][j]->Clone(Form("pairAccept_pt%d_cent%d",j,i));
			for(int ibin=1; ibin<=etaProj[i][j]->GetNbinsX(); ibin++){
				int symmBin = etaProj[i][j]->FindBin(-1*etaProj[i][j]->GetBinCenter(ibin));
				asymm[i][j]->SetBinContent(ibin, etaProj[i][j]->GetBinContent(ibin)/etaProj[i][j]->GetBinContent(symmBin));
				double err1 = etaProj[i][j]->GetBinError(ibin)/etaProj[i][j]->GetBinContent(ibin);
				double err2 = etaProj[i][j]->GetBinError(symmBin)/etaProj[i][j]->GetBinContent(symmBin);
				asymm[i][j]->SetBinError(ibin, asymm[i][j]->GetBinContent(ibin) * sqrt(pow(err1,2)+pow(err2,2)));
			}
			TF1 *fit = new TF1("fit","pol1",-3,3);
			asymm[i][j]->Fit(fit);
			
			pairAccept[i][j] = (TH2D*)GenGenSignal[i][j]->Clone(Form("pairAccept_pt%d_cent%d",j,i));
			for(int ixbin=1; ixbin<=pairAccept[i][j]->GetNbinsX(); ixbin++){
				for(int iybin=1; iybin<=pairAccept[i][j]->GetNbinsY(); iybin++){
					pairAccept[i][j]->SetBinContent(ixbin,iybin,1.-fit->Eval(asymm[i][j]->GetBinCenter(ixbin)));
					pairAccept[i][j]->SetBinError(ixbin,iybin,0.);
				}
			}
			pairAccept[i][j]->Scale(1./(double)pairAccept[i][j]->GetNbinsY());
			
			
			//fill tracking and tracking residuals
			tracking[i][j] = (TH2D*)GenGenSignal[i][j]->Clone(Form("tracking_pt%d_cent%d",j,i));
			trackingRes[i][j] = (TH2D*)GenGenSignal[i][j]->Clone(Form("trackingRes_pt%d_cent%d",j,i));
			for(int ixbin=1; ixbin<=tracking[i][j]->GetNbinsX(); ixbin++){
				for(int iybin=1; iybin<=tracking[i][j]->GetNbinsY(); iybin++){
					tracking[i][j]->SetBinContent(ixbin, iybin, 0.05);
					tracking[i][j]->SetBinError(ixbin, iybin, 0.);
					trackingRes[i][j]->SetBinContent(ixbin, iybin, 0.05);
					trackingRes[i][j]->SetBinError(ixbin, iybin, 0.);
				}
			}
			tracking[i][j]->Scale(1./(double)tracking[i][j]->GetNbinsY());
			trackingRes[i][j]->Scale(1./(double)trackingRes[i][j]->GetNbinsY());
			
			
			//fill "other" jff uncertainties
			otherJFF[i][j] = (TH2D*)GenGenSignal[i][j]->Clone(Form("otherJFF_pt%d_cent%d",j,i));
			for(int ixbin=1; ixbin<=otherJFF[i][j]->GetNbinsX(); ixbin++){
				for(int iybin=1; iybin<=otherJFF[i][j]->GetNbinsY(); iybin++){
					otherJFF[i][j]->SetBinContent(ixbin, iybin, 0.035);
					otherJFF[i][j]->SetBinError(ixbin, iybin, 0.);
				}
			}
			otherJFF[i][j]->Scale(1./(double)otherJFF[i][j]->GetNbinsY());
			
			
			//find uncertainties!
			for(int ixbin=1; ixbin<=GenGenSignal[i][j]->GetNbinsX(); ixbin++){
				for(int iybin=1; iybin<=GenGenSignal[i][j]->GetNbinsY(); iybin++){
					double xcent = GenGenSignal[i][j]->GetXaxis()->GetBinCenter(ixbin);
					double ycent = GenGenSignal[i][j]->GetYaxis()->GetBinCenter(iybin);
					double content = GenGenSignal[i][j]->GetBinContent(ixbin,iybin);
					double dr = sqrt(pow(xcent,2)+pow(ycent,2));
					
					GenGenDRDistrNorm[i][j]->Fill(dr);
					GenGenDRDistrNormWeight[i][j]->Fill(dr,content);
			
					absoluteErrSpill[i][j]->Fill(dr, content*(hSpillovers[i][j]->GetBinContent(ixbin,iybin)));
					absoluteErrResidJES[i][j]->Fill(dr, content*(hResidualJES[i][j]->GetBinContent(ixbin,iybin)));
					absoluteErrPair[i][j]->Fill(dr, content*(pairAccept[i][j]->GetBinContent(ixbin,iybin)));
					absoluteErrBG[i][j]->Fill(dr, content*(bgSub[i][j]->GetBinContent(ixbin,iybin)));
					absoluteErrTrk[i][j]->Fill(dr, content*(tracking[i][j]->GetBinContent(ixbin,iybin)));
					absoluteErrTrkRes[i][j]->Fill(dr, content*(trackingRes[i][j]->GetBinContent(ixbin,iybin)));
					absoluteErrOther[i][j]->Fill(dr, content*(otherJFF[i][j]->GetBinContent(ixbin,iybin)));
					
					relErrSpill[i][j]->Fill(dr, (hSpillovers[i][j]->GetBinContent(ixbin,iybin)));
					relErrResidJES[i][j]->Fill(dr, (hResidualJES[i][j]->GetBinContent(ixbin,iybin)));
					relErrPair[i][j]->Fill(dr, (pairAccept[i][j]->GetBinContent(ixbin,iybin)));
					relErrBG[i][j]->Fill(dr, (bgSub[i][j]->GetBinContent(ixbin,iybin)));
					relErrTrk[i][j]->Fill(dr, (tracking[i][j]->GetBinContent(ixbin,iybin)));
					relErrTrkRes[i][j]->Fill(dr, (trackingRes[i][j]->GetBinContent(ixbin,iybin)));
					relErrOther[i][j]->Fill(dr, (otherJFF[i][j]->GetBinContent(ixbin,iybin)));
				}
			}
			
			//cout << "cent " << i << " pt " << j << " " << absoluteErrSpill[i][j]->Integral() << endl;
			
			takeTH1Abs(absoluteErrSpill[i][j]);
			takeTH1Abs(absoluteErrResidJES[i][j]);
			takeTH1Abs(absoluteErrPair[i][j]);
			takeTH1Abs(absoluteErrBG[i][j]);
			takeTH1Abs(absoluteErrTrk[i][j]);
			takeTH1Abs(absoluteErrTrkRes[i][j]);
			takeTH1Abs(absoluteErrOther[i][j]);
			
			takeTH1Abs(relErrSpill[i][j]);
			takeTH1Abs(relErrResidJES[i][j]);
			takeTH1Abs(relErrPair[i][j]);
			takeTH1Abs(relErrBG[i][j]);
			takeTH1Abs(relErrTrk[i][j]);
			takeTH1Abs(relErrTrkRes[i][j]);
			takeTH1Abs(relErrOther[i][j]);
			
			absoluteErrSpill[i][j]->Divide(GenGenDRDistrNormWeight[i][j]);
			absoluteErrResidJES[i][j]->Divide(GenGenDRDistrNormWeight[i][j]);
			absoluteErrPair[i][j]->Divide(GenGenDRDistrNormWeight[i][j]);
			absoluteErrBG[i][j]->Divide(GenGenDRDistrNormWeight[i][j]);
			absoluteErrTrk[i][j]->Divide(GenGenDRDistrNormWeight[i][j]);
			absoluteErrTrkRes[i][j]->Divide(GenGenDRDistrNormWeight[i][j]);
			absoluteErrOther[i][j]->Divide(GenGenDRDistrNormWeight[i][j]);
			
			/*absoluteErrSpill[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			absoluteErrResidJES[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			absoluteErrPair[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			absoluteErrBG[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			absoluteErrTrk[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			absoluteErrTrkRes[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			absoluteErrOther[i][j]->Divide(GenGenDRDistrNorm[i][j]);*/
			
			relErrSpill[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			relErrResidJES[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			relErrPair[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			relErrBG[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			relErrTrk[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			relErrTrkRes[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			relErrOther[i][j]->Divide(GenGenDRDistrNorm[i][j]);
			
			absoluteErrSpill[i][j]->Scale(200.);
			absoluteErrResidJES[i][j]->Scale(200.);
			absoluteErrPair[i][j]->Scale(200.);
			absoluteErrBG[i][j]->Scale(200.);
			absoluteErrTrk[i][j]->Scale(200.);
			absoluteErrTrkRes[i][j]->Scale(200.);
			absoluteErrOther[i][j]->Scale(200.);
			
			relErrSpill[i][j]->Scale(200.);
			relErrResidJES[i][j]->Scale(200.);
			relErrPair[i][j]->Scale(200.);
			relErrBG[i][j]->Scale(200.);
			relErrTrk[i][j]->Scale(200.);
			relErrTrkRes[i][j]->Scale(200.);
			relErrOther[i][j]->Scale(200.);
		}
	}
	
	//smooth out that crap at low-pt
	relErrSpill[1][1]->Add(relErrSpill[0][1]);
	relErrSpill[1][1]->Add(relErrSpill[2][1]);
	relErrSpill[1][1]->Scale(1./3.);
	
	relErrSpill[2][1]->Add(relErrSpill[1][1]);
	relErrSpill[2][1]->Add(relErrSpill[3][1]);
	relErrSpill[2][1]->Scale(1./3.);
	
	for(int j=0; j<trackPtBins; j++){
		for(int i=0; i<nCentBins; i++){
			
			//cout << "cent " << i << " pt " << j << " " << relErrResidJES[i][j]->Integral() << endl;
			
			//stack everyone!
			for(int ibin=1; ibin<=absoluteErrSpill[i][j]->GetNbinsX(); ibin++){
				stackedAbsolutes[i][j]->SetBinContent(ibin, sqrt(
				                                      pow(absoluteErrSpill[i][j]->GetBinContent(ibin),2)+
				                                      pow(absoluteErrResidJES[i][j]->GetBinContent(ibin),2)+
				                                      pow(absoluteErrPair[i][j]->GetBinContent(ibin),2)+
				                                      pow(absoluteErrBG[i][j]->GetBinContent(ibin),2)+
				                                      pow(absoluteErrTrk[i][j]->GetBinContent(ibin),2)+
				                                      pow(absoluteErrTrkRes[i][j]->GetBinContent(ibin),2)+
				                                      pow(absoluteErrOther[i][j]->GetBinContent(ibin),2)));
				stackedRelatives[i][j]->SetBinContent(ibin, sqrt(
				                                      pow(relErrSpill[i][j]->GetBinContent(ibin),2)+
				                                      pow(relErrResidJES[i][j]->GetBinContent(ibin),2)+
				                                      pow(relErrPair[i][j]->GetBinContent(ibin),2)+
				                                      pow(relErrBG[i][j]->GetBinContent(ibin),2)+
				                                      pow(relErrTrk[i][j]->GetBinContent(ibin),2)+
				                                      pow(relErrTrkRes[i][j]->GetBinContent(ibin),2)+
				                                      pow(relErrOther[i][j]->GetBinContent(ibin),2)));
			}
		}
	}
	
	TCanvas *cRelErr = new TCanvas("cRelErr","",1000,800);
	cRelErr->Divide(nCentBins,trackPtBins-1);
	TLatex *l2[trackPtBins][nCentBins];
	TLatex *l3[trackPtBins][nCentBins];
	TLegend *leg1 = new TLegend(0.1,0.1,0.9,0.9);
	leg1->AddEntry(relErrSpill[0][1], "BG Fluctuation Bias","l");
	leg1->AddEntry(relErrResidJES[0][1], "JFF Bias","l");
	leg1->AddEntry(relErrOther[0][1], "Residual JES","l");
	TLegend *leg2 = new TLegend(0.1,0.1,0.9,0.9);
	leg2->AddEntry(relErrPair[0][1], "Pair Acceptance","l");
	leg2->AddEntry(relErrBG[0][1], "Background Subtr.","l");
	TLegend *leg3 = new TLegend(0.1,0.1,0.9,0.9);
	leg3->AddEntry(relErrTrk[0][1], "Tracking Eff.","l");
	leg3->AddEntry(relErrTrkRes[0][1], "Residual Tracking Eff.","l");

	for(int i=0; i<nCentBins; i++){
		for(int j=0; j<trackPtBins-1; j++){
			cRelErr->cd(j*nCentBins+(3-i)+1);
			if(j==0){
				if(i==3) leg1->Draw();
				if(i==2) leg2->Draw();
				if(i==1) leg3->Draw();
			}
			else{
				stackedRelatives[i][j]->GetYaxis()->SetNdivisions(505);
				stackedRelatives[i][j]->GetYaxis()->SetLabelSize(0.15);
				stackedRelatives[i][j]->GetXaxis()->SetLabelSize(0.15);
				stackedRelatives[i][j]->GetYaxis()->SetTitleSize(0.12);
				stackedRelatives[i][j]->GetYaxis()->SetTitleOffset(0.6);
				stackedRelatives[i][j]->GetXaxis()->SetTitleSize(0.15);
				stackedRelatives[i][j]->SetYTitle("Relative Error");
				stackedRelatives[i][j]->SetXTitle("Jet-Track #DeltaR");
				stackedRelatives[i][j]->SetMinimum(0);
				stackedRelatives[i][j]->Draw("hist");
				relErrSpill[i][j]->SetLineColor(2);
				relErrSpill[i][j]->Draw("hist,same");
				relErrResidJES[i][j]->SetLineColor(4);
				relErrResidJES[i][j]->Draw("hist,same");
				relErrPair[i][j]->SetLineColor(8);
				relErrPair[i][j]->Draw("hist,same");
				relErrBG[i][j]->SetLineColor(kOrange+2);
				relErrBG[i][j]->Draw("hist,same");
				relErrTrk[i][j]->SetLineColor(kViolet+2);
				relErrTrk[i][j]->Draw("hist,same");
				relErrTrkRes[i][j]->SetLineColor(kCyan+1);
				relErrTrkRes[i][j]->Draw("hist,same");
				relErrOther[i][j]->SetLineColor(kGreen+2);
				relErrOther[i][j]->Draw("hist,same");
				l2[i][j] = new TLatex(0.3,stackedRelatives[i][j]->GetMaximum()*0.85,Form("Cent%d-%d, %g<pT<%g",xCentBins[i],xCentBins[i+1],xTrkBinDouble[j],xTrkBinDouble[j+1]));
				l2[i][j]->SetTextSize(0.15);
				l2[i][j]->Draw("same");
			}
		}
	}
	
	TCanvas *cAbsErr = new TCanvas("cAbsErr","",1000,800);
	cAbsErr->Divide(nCentBins,trackPtBins-1);
	for(int i=0; i<nCentBins; i++){
		for(int j=0; j<trackPtBins-1; j++){
			cAbsErr->cd(j*nCentBins+(3-i)+1);
			if(j==0){
				if(i==3) leg1->Draw();
				if(i==2) leg2->Draw();
				if(i==1) leg3->Draw();
			}
			else{
				stackedAbsolutes[i][j]->GetYaxis()->SetNdivisions(505);
				stackedAbsolutes[i][j]->GetYaxis()->SetLabelSize(0.15);
				stackedAbsolutes[i][j]->GetXaxis()->SetLabelSize(0.15);
				stackedAbsolutes[i][j]->GetYaxis()->SetTitleSize(0.12);
				stackedAbsolutes[i][j]->GetYaxis()->SetTitleOffset(0.6);
				stackedAbsolutes[i][j]->GetXaxis()->SetTitleSize(0.08);
				stackedAbsolutes[i][j]->SetYTitle("Absolute Error");
				stackedAbsolutes[i][j]->SetXTitle("Jet-Track #DeltaR");
				stackedAbsolutes[i][j]->SetMinimum(0);
				stackedAbsolutes[i][j]->Draw("hist");
				absoluteErrSpill[i][j]->SetLineColor(2);
				absoluteErrSpill[i][j]->Draw("hist,same");
				absoluteErrResidJES[i][j]->SetLineColor(4);
				absoluteErrResidJES[i][j]->Draw("hist,same");
				absoluteErrPair[i][j]->SetLineColor(8);
				absoluteErrPair[i][j]->Draw("hist,same");
				absoluteErrBG[i][j]->SetLineColor(kOrange+2);
				absoluteErrBG[i][j]->Draw("hist,same");
				absoluteErrTrk[i][j]->SetLineColor(kViolet+2);
				absoluteErrTrk[i][j]->Draw("hist,same");
				absoluteErrTrkRes[i][j]->SetLineColor(kCyan+1);
				absoluteErrTrkRes[i][j]->Draw("hist,same");
				absoluteErrOther[i][j]->SetLineColor(kGreen+2);
				absoluteErrOther[i][j]->Draw("hist,same");
				l3[i][j] = new TLatex(0.3,stackedAbsolutes[i][j]->GetMaximum()*0.85,Form("Cent%d-%d, %g<pT<%g",xCentBins[i],xCentBins[i+1],xTrkBinDouble[j],xTrkBinDouble[j+1]));
				l3[i][j]->SetTextSize(0.15);
				l3[i][j]->Draw("same");
			}
		}
	}
	
	TFile *fout = new TFile("relativeErrors_RecoReco2GenGen.root","recreate");
	fout->cd();
	for(int i=0; i<nCentBins; i++){
		for(int j=0; j<trackPtBins; j++){
			relErrSpill[i][j]->Write();
			relErrResidJES[i][j]->Write();
			relErrPair[i][j]->Write();
			relErrBG[i][j]->Write();
			relErrTrk[i][j]->Write();
			relErrTrkRes[i][j]->Write();
			relErrOther[i][j]->Write();
			stackedRelatives[i][j]->Write();
		}
	}
	
	/*TH1D *hh[nCentBins][trackPtBins];
	TCanvas *debug = new TCanvas("debug","",1000,800);
	debug->Divide(nCentBins,trackPtBins-2);
	for(int i=0; i<nCentBins; i++){
		for(int j=1; j<trackPtBins-1; j++){
			debug->cd(j*nCentBins+i+1-4);
			int lowbin = pairAccept[i][j]->GetYaxis()->FindBin(-1);
			int highbin = pairAccept[i][j]->GetYaxis()->FindBin(1);
			hh[i][j] = pairAccept[i][j]->ProjectionX(Form("hh%d_%d_pfx",j,i), lowbin, highbin, "e");
			hh[i][j]->Scale(1./2.);
			hh[i][j]->Scale(1./pairAccept[i][j]->GetXaxis()->GetBinWidth(2));
			//hh[i][j]->Draw();
			etaProjRebin[i][j]->Draw();
		}
	}*/
}