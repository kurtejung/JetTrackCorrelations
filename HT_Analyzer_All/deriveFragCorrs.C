
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TGraph.h"

using namespace std;

void deriveFragCorrs(){
	
	const int nCentBins=4;
	const int CBins[nCentBins+1] = {0, 20, 60, 100, 200};
	const int nPtBins=20;
	const double xPtBins[nPtBins+1] = {30.,40.,45.,50.,55.,60.,65.,70.,80.,85.,90.,95.,100.,110.,120.,140.,170.,220.,280.,370.,600.};
	
	const int nIters = 2;
		
	//Step1 histograms
	TH1D *csCands[nPtBins];
	TH1D *csCandsPost[nPtBins];
	TH2D *preCorrs[nCentBins][nPtBins];
	TH2D *preQCorrs[nCentBins][nPtBins];
	TH2D *preGCorrs[nCentBins][nPtBins];
	TH2D *postCorrs[nCentBins][nPtBins];
	TH2D *postQCorrs[nCentBins][nPtBins];
	TH2D *postGCorrs[nCentBins][nPtBins];
	TProfile *profs[nCentBins][nPtBins];
	TF1 *nCSfits[nCentBins][nPtBins];
	TH1D *ncsConstFitVals[nCentBins];
	TH1D *ncsSlopeFitVals[nCentBins];
	TF1 *constFit[nCentBins];
	TF1 *slopeFit[nCentBins];
	TH2D *s1recoClosure[nCentBins];
	TProfile *s1recoClosureProf[nCentBins];
	TH2D *s1QrecoClosure[nCentBins];
	TH2D *s1GrecoClosure[nCentBins];
	TH2D *s1refClosure[nCentBins];
	TH2D *s1QrefClosure[nCentBins];
	TH2D *s1GrefClosure[nCentBins];
	
	//Iterative step histograms
	TH2D *s2recoClosure[nCentBins][nIters];
	TH2D *s2QrecoClosure[nCentBins][nIters];
	TH2D *s2GrecoClosure[nCentBins][nIters];
	TH2D *s2refClosure[nCentBins][nIters];
	TH2D *s2QrefClosure[nCentBins][nIters];
	TH2D *s2GrefClosure[nCentBins][nIters];
	
	TH1D *iterClosures[nCentBins][nIters];
	TF1 *iterFits[nCentBins][nIters];
	TGraph *testGr[nCentBins][nIters];
	//TGraph *outGr[nCentBins][nIters];
	//TGraphSmooth *smoothedGr = new TGraphSmooth("normal");
	
	//Final Closure histograms
	TH2D *s3recoClosure[nCentBins];
	TH2D *s3QrecoClosure[nCentBins];
	TH2D *s3GrecoClosure[nCentBins];
	TH2D *s3refClosure[nCentBins];
	TH2D *s3QrefClosure[nCentBins];
	TH2D *s3GrefClosure[nCentBins];
	
	for(int i=0; i<nCentBins; i++){
		
		for(int k=0; k<nIters; k++){
			iterClosures[i][k] = new TH1D(Form("iterClosures_cent%d_iter%d",i,k),"",240,40,1000);
		}
		
		for(int j=0; j<nPtBins; j++){
			preCorrs[i][j] = new TH2D(Form("h_cs_preclosure_cent%d_pt%d",i,j),"",40,0.,40.,200,0.,2.);
			preCorrs[i][j]->Sumw2();
			preQCorrs[i][j] = new TH2D(Form("h_cs_q_preclosure_cent%d_pt%d",i,j),"",40,0.,40.,200,0.,2.);
			preQCorrs[i][j]->Sumw2();
			preGCorrs[i][j] = new TH2D(Form("h_cs_g_preclosure_cent%d_pt%d",i,j),"",40,0.,40.,200,0.,2.);
			preGCorrs[i][j]->Sumw2();
			
			postCorrs[i][j] = new TH2D(Form("h_cs_postclosure_cent%d_pt%d",i,j),"",40,0.,40.,200,0.,2.);
			postCorrs[i][j]->Sumw2();
			postQCorrs[i][j] = new TH2D(Form("h_cs_q_postclosure_cent%d_pt%d",i,j),"",40,0.,40.,200,0.,2.);
			postQCorrs[i][j]->Sumw2();
			postGCorrs[i][j] = new TH2D(Form("h_cs_g_postclosure_cent%d_pt%d",i,j),"",40,0.,40.,200,0.,2.);
			postGCorrs[i][j]->Sumw2();
		}
		
		ncsConstFitVals[i] = new TH1D(Form("ncsConstFitVals_cent%d",i),"",nPtBins,xPtBins); ncsConstFitVals[i]->Sumw2();
		ncsSlopeFitVals[i] = new TH1D(Form("ncsSlopeFitVals_cent%d",i),"",nPtBins,xPtBins); ncsSlopeFitVals[i]->Sumw2();
			
		s1recoClosure[i] = new TH2D(Form("s1recoClosure_cent%d",i),"",240,40,1000,200,0.,2.);
		s1recoClosure[i]->Sumw2();
		s1refClosure[i] = new TH2D(Form("s1refClosure_cent%d",i),"",240,40,1000,200,0.,2.);
		s1refClosure[i]->Sumw2();
		s1QrecoClosure[i] = new TH2D(Form("s1QrecoClosure_cent%d",i),"",240,40,1000,200,0.,2.);
		s1QrecoClosure[i]->Sumw2();
		s1QrefClosure[i] = new TH2D(Form("s1QrefClosure_cent%d",i),"",240,40,1000,200,0.,2.);
		s1QrefClosure[i]->Sumw2();
		s1GrecoClosure[i] = new TH2D(Form("s1GrecoClosure_cent%d",i),"",240,40,1000,200,0.,2.);
		s1GrecoClosure[i]->Sumw2();
		s1GrefClosure[i] = new TH2D(Form("s1GrefClosure_cent%d",i),"",240,40,1000,200,0.,2.);
		s1GrefClosure[i]->Sumw2();
		
		for(int k=0; k<nIters; k++){
			s2recoClosure[i][k] = new TH2D(Form("s2recoClosure_cent%d_iter%d",i,k),"",240,40,1000,200,0.,2.);
			s2recoClosure[i][k]->Sumw2();
			s2refClosure[i][k] = new TH2D(Form("s2refClosure_cent%d_iter%d",i,k),"",240,40,1000,200,0.,2.);
			s2refClosure[i][k]->Sumw2();
			s2QrecoClosure[i][k] = new TH2D(Form("s2QrecoClosure_cent%d_iter%d",i,k),"",240,40,1000,200,0.,2.);
			s2QrecoClosure[i][k]->Sumw2();
			s2QrefClosure[i][k] = new TH2D(Form("s2QrefClosure_cent%d_iter%d",i,k),"",240,40,1000,200,0.,2.);
			s2QrefClosure[i][k]->Sumw2();
			s2GrecoClosure[i][k] = new TH2D(Form("s2GrecoClosure_cent%d_iter%d",i,k),"",240,40,1000,200,0.,2.);
			s2GrecoClosure[i][k]->Sumw2();
			s2GrefClosure[i][k] = new TH2D(Form("s2GrefClosure_cent%d_iter%d",i,k),"",240,40,1000,200,0.,2.);
			s2GrefClosure[i][k]->Sumw2();
						
		}
		
		s3recoClosure[i] = new TH2D(Form("s3recoClosure_cent%d",i),"",115,60,500,200,0.,2.);
		s3recoClosure[i]->Sumw2();
		s3refClosure[i] = new TH2D(Form("s3refClosure_cent%d",i),"",115,60,500,200,0.,2.);
		s3refClosure[i]->Sumw2();
		s3QrecoClosure[i] = new TH2D(Form("s3QrecoClosure_cent%d",i),"",115,60,500,200,0.,2.);
		s3QrecoClosure[i]->Sumw2();
		s3QrefClosure[i] = new TH2D(Form("s3QrefClosure_cent%d",i),"",115,60,500,200,0.,2.);
		s3QrefClosure[i]->Sumw2();
		s3GrecoClosure[i] = new TH2D(Form("s3GrecoClosure_cent%d",i),"",115,60,500,200,0.,2.);
		s3GrecoClosure[i]->Sumw2();
		s3GrefClosure[i] = new TH2D(Form("s3GrefClosure_cent%d",i),"",115,60,500,200,0.,2.);
		s3GrefClosure[i]->Sumw2();
	}
	for(int j=0; j<nPtBins; j++){
		csCands[j] = new TH1D(Form("csCands_pt%d",j),"",40,0,40);
		csCands[j]->Sumw2();
		csCandsPost[j] = new TH1D(Form("csCandsPost_pt%d",j),"",40,0,40);
		csCandsPost[j]->Sumw2();
	}
	
	TFile *f1 = new TFile("/Users/kjung/Run2Validation/CScandClosure/unzippedSkimClosure.root");
	TTree *t = (TTree*)f1->Get("unzipMixTree");
	
	double jtpt, jteta, refpt, weight;
	int nCScand, bin, refparton_flavor;
	t->SetBranchAddress("corrpt", &jtpt);
	t->SetBranchAddress("refpt", &refpt);
	t->SetBranchAddress("jteta", &jteta);
	t->SetBranchAddress("nCScand", &nCScand);
	t->SetBranchAddress("hiBin",&bin);
	t->SetBranchAddress("refparton_flavor",&refparton_flavor);
	t->SetBranchAddress("weight", &weight);
	
	//Step1 - Derive the nCS-dependent corrections
	cout << "Deriving nCS corrs..." << endl;
	for(int ientry=0; ientry<t->GetEntries(); ientry++){
		if(ientry && ientry%1000000==0) cout << "at entry "<< ientry << endl;
		
		t->GetEntry(ientry);
		
		if(jtpt>600 || jtpt<30) continue;
		if(abs(jteta)>1.6) continue;
		//if(refpt<50) continue;
		//if(refparton_flavor == -999) continue;
		
		int mycbin=-1, myptbin=-1;
		for (int cbin = 0; cbin < nCentBins; cbin++){
			if (bin >= CBins[cbin] && bin < CBins[cbin+1]){
				mycbin = cbin;
			}
		}
		
		for (int ptbin = 0; ptbin < nPtBins; ptbin++){
			if (jtpt > xPtBins[ptbin] && jtpt <= xPtBins[ptbin+1]){
				myptbin = ptbin;
			}
		}
		
		csCands[myptbin]->Fill(nCScand,weight);
		preCorrs[mycbin][myptbin]->Fill(nCScand,jtpt/refpt,weight);
		s1recoClosure[mycbin]->Fill(jtpt,jtpt/refpt,weight);
		s1refClosure[mycbin]->Fill(refpt,jtpt/refpt,weight);
		if(abs(refparton_flavor)>=1 && abs(refparton_flavor)<=6){
			s1QrecoClosure[mycbin]->Fill(jtpt,jtpt/refpt,weight);
			s1QrefClosure[mycbin]->Fill(refpt,jtpt/refpt,weight);
			preQCorrs[mycbin][myptbin]->Fill(nCScand,jtpt/refpt,weight);
		}
		if(abs(refparton_flavor)==21){
			s1GrecoClosure[mycbin]->Fill(jtpt,jtpt/refpt,weight);
			s1GrefClosure[mycbin]->Fill(refpt,jtpt/refpt,weight);
			preGCorrs[mycbin][myptbin]->Fill(nCScand,jtpt/refpt,weight);
		}
	}
	
	for(int i=0; i<nCentBins; i++){
		s1recoClosureProf[i] =  s1recoClosure[i]->ProfileX();
		for(int j=0; j<nPtBins; j++){
			profs[i][j] = preCorrs[i][j]->ProfileX();
			nCSfits[i][j] = new TF1(Form("fits_%d_%d",i,j),"pol1");
			if(j<9) profs[i][j]->Fit(nCSfits[i][j],"","",2,13);
			else if(j<16) profs[i][j]->Fit(nCSfits[i][j],"","",5,20);
			else profs[i][j]->Fit(nCSfits[i][j],"","",8,30);
			
			ncsConstFitVals[i]->SetBinContent(j+1, nCSfits[i][j]->GetParameter(0));
			ncsConstFitVals[i]->SetBinError(j+1, nCSfits[i][j]->GetParError(0));
			ncsSlopeFitVals[i]->SetBinContent(j+1, nCSfits[i][j]->GetParameter(1));
			ncsSlopeFitVals[i]->SetBinError(j+1, nCSfits[i][j]->GetParError(1));
			
		}
		constFit[i] = new TF1(Form("constFit_cent%d",i),"[3]+landau");
		constFit[i]->SetParLimits(1,80,150);
		constFit[i]->SetParLimits(2,10,50);
		constFit[i]->SetParLimits(3,1,1.4);
		//constFit[i]->SetParLimits(4,1e-2,5);
		
		slopeFit[i] = new TF1(Form("slopeFit_cent%d",i),"pol3");
		ncsConstFitVals[i]->Fit(constFit[i],"","",70,600);
		ncsSlopeFitVals[i]->Fit(slopeFit[i],"","",80,600);
		
	}
	
	//Step 2 - Apply Doga's iterative method
	cout << "Deriving iterations..." << endl;
	for(int iter=0; iter<nIters+1; iter++){
		cout << "starting iteration " << iter+1 << endl;
		for(int ientry=0; ientry<t->GetEntries(); ientry++){
			if(ientry && ientry%1000000==0) cout << "at entry "<< ientry << endl;
			t->GetEntry(ientry);
			
			if(jtpt>600 || jtpt<30) continue;
			if(abs(jteta)>1.6) continue;
			//if(refpt<50) continue;
			//if(refparton_flavor == -999) continue;
			
			int mycbin=-1, myptbin=-1;
			for (int cbin = 0; cbin < nCentBins; cbin++){
				if (bin >= CBins[cbin] && bin < CBins[cbin+1]){
					mycbin = cbin;
				}
			}
			
			for (int ptbin = 0; ptbin < nPtBins; ptbin++){
				if (jtpt > xPtBins[ptbin] && jtpt <= xPtBins[ptbin+1]){
					myptbin = ptbin;
				}
			}
			
			double corrFactor;
			//if(jtpt>70) corrFactor = constFit[mycbin]->Eval(jtpt)+(slopeFit[mycbin]->Eval(jtpt)*nCScand);
			corrFactor = nCSfits[mycbin][myptbin]->Eval(nCScand);
			int ibin=s1recoClosureProf[mycbin]->FindBin(jtpt);
			corrFactor/=s1recoClosureProf[mycbin]->GetBinContent(ibin);
			double corrpt = jtpt/corrFactor;
			//double corrpt = jtpt/nCSfits[mycbin][myptbin]->Eval(nCScand);
			
			int corrptbin=0;
			while(corrpt>xPtBins[corrptbin] && corrptbin<nPtBins){
				corrptbin++;
			}
			
			if(!iter){
				csCandsPost[myptbin]->Fill(nCScand,weight);
				postCorrs[mycbin][myptbin]->Fill(nCScand,corrpt/refpt,weight);
				if(abs(refparton_flavor)>=1 && abs(refparton_flavor)<=6) postQCorrs[mycbin][myptbin]->Fill(nCScand,corrpt/refpt,weight);
				if(abs(refparton_flavor)==21) postGCorrs[mycbin][myptbin]->Fill(nCScand,corrpt/refpt,weight);
			}
			
			if(iter){
				int dummy=iter;
				while(dummy>0){

					//if(corrpt<200){
					//	corrpt = corrpt/iterFits[mycbin][dummy-1]->Eval(corrpt);
						
						//TSpline interpolation smoothes out a lot of the bumpiness from higher-order polynomial fitting
						//corrpt = corrpt/testGr[mycbin][dummy-1]->Eval(corrpt,0,"S");
					//}
					//else{
						int corrBin = iterClosures[mycbin][dummy-1]->FindBin(corrpt);
						corrpt = corrpt/iterClosures[mycbin][dummy-1]->GetBinContent(corrBin);
					//}
					dummy--;
				} 
			}
			
			//If we're one past nIters, then fill the final closures and be done
			if(iter==nIters){
				s3recoClosure[mycbin]->Fill(corrpt,corrpt/refpt,weight);
				s3refClosure[mycbin]->Fill(refpt,corrpt/refpt,weight);
				if(abs(refparton_flavor)>=1 && abs(refparton_flavor)<=6){
					s3QrecoClosure[mycbin]->Fill(corrpt,corrpt/refpt,weight);
					s3QrefClosure[mycbin]->Fill(refpt,corrpt/refpt,weight);
				}
				if(abs(refparton_flavor)==21){
					s3GrecoClosure[mycbin]->Fill(corrpt,corrpt/refpt,weight);
					s3GrefClosure[mycbin]->Fill(refpt,corrpt/refpt,weight);
				}
			}
			//Else fill the temporary corrections and reiterate...
			else{
				s2recoClosure[mycbin][iter]->Fill(corrpt,corrpt/refpt,weight);
				s2refClosure[mycbin][iter]->Fill(refpt,corrpt/refpt,weight);
				if(abs(refparton_flavor)>=1 && abs(refparton_flavor)<=6){
					s2QrecoClosure[mycbin][iter]->Fill(corrpt,corrpt/refpt,weight);
					s2QrefClosure[mycbin][iter]->Fill(refpt,corrpt/refpt,weight);
				}
				if(abs(refparton_flavor)==21){
					s2GrecoClosure[mycbin][iter]->Fill(corrpt,corrpt/refpt,weight);
					s2GrefClosure[mycbin][iter]->Fill(refpt,corrpt/refpt,weight);
				}
			}
		}
		
		//g(x) = recopt
		//f(x) = refpt
		//g(f(jtpt)*jtpt) = iterClosure
		if(iter<nIters){
			for(int i=0; i<nCentBins; i++){
				//int lowbin = s2recoClosure[i][iter]->GetYaxis()->FindBin(0.65);
				//int highbin = s2recoClosure[i][iter]->GetYaxis()->FindBin(1.35);
				TProfile *recoClos = s2recoClosure[i][iter]->ProfileX();//"_pfx",lowbin,highbin);
				TProfile *refClos = s2refClosure[i][iter]->ProfileX();//"_pfx",lowbin,highbin);
				double xbins[240];
				double ybins[240];
				for(int ibin=1; ibin<=recoClos->GetNbinsX(); ibin++){
					
					//refClos bin content = f(jtpt), recoClos center = jtpt
					double arg = refClos->GetBinContent(ibin)*recoClos->GetBinCenter(ibin);
					double gOfF = refClos->GetBinContent(refClos->FindBin(arg));
					//cout << "for bin " << ibin << " pt " << recoClos->GetBinCenter(ibin) << " gOfF = " << gOfF << endl;
					iterClosures[i][iter]->SetBinContent(ibin, gOfF);
					
					xbins[ibin-1] = iterClosures[i][iter]->GetBinCenter(ibin);
					ybins[ibin-1] = iterClosures[i][iter]->GetBinContent(ibin);
				}
				testGr[i][iter] = new TGraph(iterClosures[i][iter]->GetNbinsX(),xbins,ybins);
				//outGr[i][iter] = smoothedGr->Approx(testGr[i][iter],"linear",50,0,1,1);
				
				iterFits[i][iter] = new TF1(Form("IterFits_cent%d_iter%d",i,iter),"pol3");
				iterClosures[i][iter]->Fit(iterFits[i][iter],"","",80,200);
			}
		}
	}
	
	TFile *fout = new TFile("fullIterativeJFFClosures_closureTest.root","recreate");
	fout->cd();
	
	for(int j=0; j<nPtBins; j++){
		csCands[j]->Write();
		csCandsPost[j]->Write();
	}
	
	for(int i=0; i<nCentBins; i++){
		for(int j=0; j<nPtBins; j++){
			preCorrs[i][j]->Write();
			preQCorrs[i][j]->Write();
			preGCorrs[i][j]->Write();
			postCorrs[i][j]->Write();
			postQCorrs[i][j]->Write();
			postGCorrs[i][j]->Write();
			nCSfits[i][j]->Write();
		}
		s1recoClosure[i]->Write();
		s1QrecoClosure[i]->Write();
		s1GrecoClosure[i]->Write();
		s1refClosure[i]->Write();
		s1QrefClosure[i]->Write();
		s1GrefClosure[i]->Write();
		constFit[i]->Write();
		slopeFit[i]->Write();
		ncsConstFitVals[i]->Write();
		ncsSlopeFitVals[i]->Write();
		
		for(int k=0; k<nIters; k++){
			s2recoClosure[i][k]->Write();
			s2QrecoClosure[i][k]->Write();
			s2GrecoClosure[i][k]->Write();
			s2refClosure[i][k]->Write();
			s2QrefClosure[i][k]->Write();
			s2GrefClosure[i][k]->Write();
			iterClosures[i][k]->Write();
			iterFits[i][k]->Write();
		}
	
		s3recoClosure[i]->Write();
		s3QrecoClosure[i]->Write();
		s3GrecoClosure[i]->Write();
		s3refClosure[i]->Write();
		s3QrefClosure[i]->Write();
		s3GrefClosure[i]->Write();
	}
	
	fout->Close();
	
	TFile *fout2 = new TFile("nCScorrectionFactors.root","recreate");
	fout2->cd();
	for(int i=0; i<nCentBins; i++){
		for(int j=0; j<nPtBins; j++){
			nCSfits[i][j]->Write();

		}
		for(int k=0; k<nIters; k++){
			iterClosures[i][k]->Write();
		}
		s1recoClosureProf[i]->Write();
	}
	fout2->Close();
	
}