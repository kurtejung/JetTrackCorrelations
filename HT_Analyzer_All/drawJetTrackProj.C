

void drawJetTrackProj(){
	
	const int nTrkPtBins = 5;
	const int jetPtBins = 1;
	const int centralityBins = 1;
	const int collType = 2;

	//pp  bins are filled duplicate in every cent bin so just picking 1 is fine
	const int centBins[centralityBins+1] = {0,10};
	const int jptBins[jetPtBins+1] = {100,300};
	const int trkPtBins[nTrkPtBins+1] = {1, 2, 3, 4, 8, 999};

	double nearSideLow = -0.785; //(-pi/4)
	double nearSideHigh = 0.785; //(+pi/4)

	double awaySideLow = 2.356; //(3pi/4)
	double awaySideHigh = 3.93; //(5pi/4);
	
	TH2D *inclCorrPlots[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH2D *bjetCorrPlots[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclFG_nearPeak[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bFG_nearPeak[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclFG_awayPeak[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bFG_awayPeak[centralityBins][jetPtBins][nTrkPtBins][collType];

	TH1D *inclBG[centralityBins][jetPtBins][nTrkPtBins][collType];
	TH1D *bBG[centralityBins][jetPtBins][nTrkPtBins][collType];

	TFile *inIncl = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/pp_Data_5TeV_InclusiveJet.root");
	TFile *inB = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/pp_Data_5TeV_csv0p9Filter.root");

	TFile *inInclPbPb = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/PbPb_2015data_inclusive.root");
        TFile *inBPbPb = new TFile("/Users/kjung/JetTrackCorrelations/HT_Analyzer_All/root_output/PbPb_2015data_bJetFiltered_csv0p9.root");

	TH1D *nEvents = (TH1D*)inIncl->Get("Data_Nevents_dijets");
	int InclEvents = nEvents->GetEntries();

	TH1D *bnEvents = (TH1D*)inB->Get("Data_Nevents_dijets");
	int bEvents = bnEvents->GetEntries();

	TH1D *nEventsPbPb = (TH1D*)inInclPbPb->Get("Data_Nevents_dijets");
	int InclEvents_PbPb = nEvents->GetEntries();

	TH1D *bnEventsPbPb = (TH1D*)inBPbPb->Get("Data_Nevents_dijets");
	int bEvents_PbPb = bnEvents->GetEntries();


	for(int i=0; i<centralityBins; i++){
		for(int j=0; j<jetPtBins; j++){
			for(int k=0; k<nTrkPtBins; k++){
				inclCorrPlots[i][j][k][0] = (TH2D*)inIncl->Get(Form("Data_hJetTrackSignalBackgroundLeadingCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_pp",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				bjetCorrPlots[i][j][k][0] = (TH2D*)inB->Get(Form("Data_hJetTrackSignalBackgroundLeadingCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_pp",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));

				inclCorrPlots[i][j][k][1] = (TH2D*)inInclPbPb->Get(Form("Data_hJetTrackSignalBackgroundLeadingCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				bjetCorrPlots[i][j][k][1] = (TH2D*)inBPbPb->Get(Form("Data_hJetTrackSignalBackgroundLeadingCent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]))->Clone(Form("JetTrackCorr_LeadingJet_Cent%d_Cent%d_Pt%d_Pt%d_TrkPt%d_TrkPt%d_PbPb",centBins[i],centBins[i+1],jptBins[j],jptBins[j+1],trkPtBins[k],trkPtBins[k+1]));
				
				for(int m=0; m<collType; m++){			
	
				//Project near-side peak
					int lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(nearSideLow);
					int highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(nearSideHigh);
					inclFG_nearPeak[i][j][k][m] = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionX(Form("_pxNear_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					bFG_nearPeak[i][j][k][m] = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionX(Form("_bpxNear_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");

				//now do the middle stuff as BG
					lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(nearSideHigh);
					highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(awaySideLow);
					inclBG[i][j][k][m] = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionX(Form("_pxBG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					bBG[i][j][k][m] = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionX(Form("_bpxBG_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");


				//now do the away-side peak
					lowBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(awaySideLow);
					highBin = inclCorrPlots[i][j][k][m]->GetYaxis()->FindBin(awaySideHigh);
					inclFG_awayPeak[i][j][k][m] = (TH1D*)inclCorrPlots[i][j][k][m]->ProjectionX(Form("_pxAway_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
					bFG_awayPeak[i][j][k][m] = (TH1D*)bjetCorrPlots[i][j][k][m]->ProjectionX(Form("_bpxAway_cent%d_jtpt%d_trkpt%d_coll%d",i,j,k,m),lowBin, highBin,"e");
				}
			}
		}
	}

	int tptBin = 0;
	double bFrac = 1./0.035/0.90; //  1/bfrac/taggingPurity
	int collision = 1;

	TCanvas *cc = new TCanvas("cc","",1200,600);
	cc->Divide(3,1);
	int dummy=1;
	for(int tptBin=0; tptBin<5; tptBin+=2){
	cc->cd(dummy);
	cc->GetPad(dummy++)->SetLogy();
	//inclFG_nearPeak[0][0][tptBin]->Scale(1./(double)nEvents);
	//inclBG[0][0][tptBin]->Scale(1./(double)nEvents);
	inclFG_nearPeak[0][0][tptBin][collision]->SetXTitle("Jet #eta");
	inclFG_nearPeak[0][0][tptBin][collision]->SetYTitle("Jet Signal (peaks / yield [#pi/2,3#pi/2])");
	inclFG_nearPeak[0][0][tptBin][collision]->Divide(inclBG[0][0][tptBin][collision]);
	inclFG_nearPeak[0][0][tptBin][collision]->Draw();
	inclFG_awayPeak[0][0][tptBin][collision]->SetMarkerColor(4);
	inclFG_awayPeak[0][0][tptBin][collision]->SetLineColor(4);
	//inclFG_awayPeak[0][0][tptBin][0]->Scale(1./(double)nEvents);
	inclFG_awayPeak[0][0][tptBin][collision]->Divide(inclBG[0][0][tptBin][collision]);
	//inclFG_awayPeak[0][0][tptBin][collision]->Draw("same");

	bFG_nearPeak[0][0][tptBin][collision]->Divide(bBG[0][0][tptBin][collision]);
	bFG_nearPeak[0][0][tptBin][collision]->SetMarkerStyle(24);
	//bFG_nearPeak[0][0][tptBin][0]->Scale(bFrac);
	bFG_nearPeak[0][0][tptBin][collision]->Draw("same");

	bFG_awayPeak[0][0][tptBin][collision]->SetMarkerColor(4);
	bFG_awayPeak[0][0][tptBin][collision]->SetLineColor(4);
	bFG_awayPeak[0][0][tptBin][collision]->SetMarkerStyle(24);
	bFG_awayPeak[0][0][tptBin][collision]->Divide(bBG[0][0][tptBin][collision]);
	//bFG_awayPeak[0][0][tptBin]->Scale(bFrac);
	//bFG_awayPeak[0][0][tptBin][collision]->Draw("same");
	}


}
