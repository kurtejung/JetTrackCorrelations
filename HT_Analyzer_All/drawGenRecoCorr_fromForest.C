{
	

	string filename = "/eos/uscms/store/user/htrauger/Pythia8_Dijet80_pp502_TuneCUETP8M1/HINppWinter16DR-75X_mcRun2_asymptotic_ppAt5TeV_v3-v1_HiForest/160716_123010/0000/HiForestAOD_*.root";

	TChain *ch = new TChain("HiGenParticleAna/hi");
	ch->Add(filename.c_str());
	TChain *chReco = new TChain("ppTrack/trackTree");
	chReco->Add(filename.c_str());
	ch->AddFriend(chReco);
	TChain *chEvt = new TChain("skimanalysis/HltTree");
	chEvt->Add(filename.c_str());
	ch->AddFriend(chEvt);
	TChain *chHLT = new TChain("hltanalysis/HltTree");
	chHLT->Add(filename.c_str());
	ch->AddFriend(chHLT);

	const int MAXPARTICLES=1000;

	vector<float> *pt=0, *eta=0;
	vector<int> *chg=0;
	bool highPurity[MAXPARTICLES];
	float trackMax[MAXPARTICLES], trkDxy1[MAXPARTICLES], trkDxyError1[MAXPARTICLES], trkDz1[MAXPARTICLES], trkDzError1[MAXPARTICLES], trkPtError[MAXPARTICLES], trkMVATight[MAXPARTICLES], pfEcal[MAXPARTICLES], pfHcal[MAXPARTICLES];	
	float trkPt[MAXPARTICLES], trkEta[MAXPARTICLES];

	int pBeamScrapingFilter, pVertexFilterCutG, pPAprimaryVertexFilter, HBHENoiseFilterResult;

	int nref, nTrk;
	int HLT_Jet80;

	ch->SetBranchAddress("pt",&pt);
	ch->SetBranchAddress("eta",&eta);
	ch->SetBranchAddress("chg",&chg);

	ch->SetBranchAddress("nTrk",&nTrk);
	ch->SetBranchAddress("trkPt",trkPt);
	ch->SetBranchAddress("trkEta",trkEta);
	ch->SetBranchAddress("trkDxy1", trkDxy1);
	ch->SetBranchAddress("trkDxyError1", trkDxyError1);
	ch->SetBranchAddress("trkDz1", trkDz1);
	ch->SetBranchAddress("trkDzError1", trkDzError1);
	ch->SetBranchAddress("trkPtError", trkPtError);
	ch->SetBranchAddress("highPurity",highPurity);
	ch->SetBranchAddress("trkMVATight",trkMVATight);
	ch->SetBranchAddress("pfEcal",pfEcal);
	ch->SetBranchAddress("pfHcal",pfHcal);

	ch->SetBranchAddress("HLT_AK4CaloJet80_Eta5p1ForPPRef_v1",&HLT_Jet80);

	ch->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
	ch->SetBranchAddress("pVertexFilterCutGplus",&pVertexFilterCutG);
	ch->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
	ch->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);

	TH2D *hh = new TH2D("hh","",300,0,300,300,0,300);

	TH1D *trackPtCut = new TH1D("trackPtCut","",60,0,20); trackPtCut->Sumw2();
	TH1D *trackPtRef = new TH1D("trackPtRef","",60,0,20); trackPtRef->Sumw2();

	for(int i=0; i<ch->GetEntries(); i++){

		ch->GetEntry(i);

		//Event cuts
		if(!pBeamScrapingFilter) continue;
		if(!pPAprimaryVertexFilter) continue;
		if(!HBHENoiseFilterResult) continue;
		if(!HLT_Jet80) continue;

		int ngen=0, nreco=0;
		for(int igen=0; igen<pt->size(); igen++){
			if(abs(eta->at(igen))<2.4 && chg->at(igen)!=0 && pt->at(igen)>0.5 && pt->at(igen)<300) ngen++;
		}
		for(int ireco=0; ireco<nTrk; ireco++){
			//if(highPurity[ireco]!=1) continue;
			if(trkPtError[ireco]/trkPt[ireco]>=0.3 || TMath::Abs(trkDz1[ireco]/trkDzError1[ireco])>=3.0 ||TMath::Abs(trkDxy1[ireco]/trkDxyError1[ireco])>=3.0) continue ;
			float Et = (pfHcal[ireco]+pfEcal[ireco])/TMath::CosH(trkEta[ireco]);
			//if(!(trkPt[ireco]<20 || (Et>0.2*trkPt[ireco] && Et>trkPt[ireco]-80))) continue;

			if(abs(trkEta[ireco])>2.4) continue;
			if(trkPt[ireco]<0.5 || trkPt[ireco]>300) continue;

			nreco++;
			
		}
		for(int ireco=0; ireco<nTrk; ireco++){
			if(highPurity[ireco]!=1) continue;
			if(trkPtError[ireco]/trkPt[ireco]>=0.3 || TMath::Abs(trkDz1[ireco]/trkDzError1[ireco])>=3.0 ||TMath::Abs(trkDxy1[ireco]/trkDxyError1[ireco])>=3.0) continue ;
			float Et = (pfHcal[ireco]+pfEcal[ireco])/TMath::CosH(trkEta[ireco]);
			if(!(trkPt[ireco]<20 || (Et>0.2*trkPt[ireco] && Et>trkPt[ireco]-80))) continue;

			if(abs(trkEta[ireco])>2.4) continue;
			if(trkPt[ireco]<0.5 || trkPt[ireco]>300) continue;

			if(nreco-ngen>50){ trackPtCut->Fill(trkPt[ireco]); }
			else{ trackPtRef->Fill(trkPt[ireco]); }

		}

		hh->Fill(ngen, nreco);
		
	}

	TCanvas *c1 = new TCanvas("c1","",600,600);
	c1->cd();
	hh->Draw("colz");

	TCanvas *c2 = new TCanvas("c2","",600,600);
	c2->cd();
	trackPtCut->Scale(1./trackPtCut->Integral(30,60));
	trackPtRef->Scale(1./trackPtRef->Integral(30,60));
	trackPtRef->Draw();
	trackPtCut->SetLineColor(2);
	trackPtCut->Draw("same");


}