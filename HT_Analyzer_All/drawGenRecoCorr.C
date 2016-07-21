{
	
	TChain *ch = new TChain("mixing_tree");
	ch->Add("/data/htrauger/JetShapes2016Skims/Pythia_July19/*.root");

	vector<float> *pt=0, *eta=0, *trkPt=0, *jtpt=0;
	vector<int> *chg=0;
	vector<float> *highPurity=0;
	vector<float> *trackMax=0, *trkDxy1=0, *trkDxyError1=0, *trkDz1=0, *trkDzError1=0, *trkPtError=0, *trkMVATight=0;

	vector<float> *trkEta=0;

	ch->SetBranchAddress("pt",&pt);
	ch->SetBranchAddress("eta",&eta);
	ch->SetBranchAddress("chg",&chg);

	ch->SetBranchAddress("calo_jtpt",&jtpt);

	ch->SetBranchAddress("trkPt",&trkPt);
	ch->SetBranchAddress("trkEta",&trkEta);
	ch->SetBranchAddress("trkDxy1", &trkDxy1);
	ch->SetBranchAddress("trkDxyError1", &trkDxyError1);
	ch->SetBranchAddress("trkDz1", &trkDz1);
	ch->SetBranchAddress("trkDzError1", &trkDzError1);
	ch->SetBranchAddress("trkPtError", &trkPtError);
	ch->SetBranchAddress("highPurity",&highPurity);
	ch->SetBranchAddress("trkMVATight",&trkMVATight);
	//ch->SetBranchAddress("pfEcal",&pfEcal);
	//ch->SetBranchAddress("pfHcal",&pfHcal);

	TH2D *hh = new TH2D("hh","",300,0,300,300,0,300);

	TH1D *trackPtCut = new TH1D("trackPtCut","",60,0,20); trackPtCut->Sumw2();
	TH1D *trackPtRef = new TH1D("trackPtRef","",60,0,20); trackPtRef->Sumw2();

	TH1D *jtptCut = new TH1D("jtptCut","",100,30,500); jtptCut->Sumw2();
	TH1D *jtptRef = new TH1D("jtptRef","",100,30,500); jtptRef->Sumw2();

	for(int i=0; i<ch->GetEntries(); i++){

		ch->GetEntry(i);
		int ngen=0, nreco=0;
		for(int igen=0; igen<pt->size(); igen++){
			if(abs(eta->at(igen))<2.4 && chg->at(igen)!=0 && pt->at(igen)>5) ngen++;
		}
		for(int ireco=0; ireco<trkPt->size(); ireco++){
			if(trkPt->at(ireco)<5) continue;
			if(abs(trkEta->at(ireco))>2.4) continue;
			if(trkMVATight->at(ireco)<0.5) continue;

			nreco++;
			
		}
		for(int ireco=0; ireco<trkPt->size(); ireco++){
			if(trkPt->at(ireco)<5) continue;
			if(abs(trkEta->at(ireco))>2.4) continue;
			if(trkMVATight->at(ireco)<0.5) continue;

			if(nreco-ngen>50){ trackPtCut->Fill(trkPt->at(ireco)); }
			else{ trackPtRef->Fill(trkPt->at(ireco)); }

		}

		hh->Fill(ngen, nreco);

		for(int ijet=0; ijet<jtpt->size(); ijet++){
			if(nreco-ngen>50){ jtptCut->Fill(jtpt->at(ijet)); }
			else{ jtptRef->Fill(jtpt->at(ijet)); }
		}
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

	TCanvas *c3 = new TCanvas("c3","",600,600);
	c3->cd();
	jtptCut->Scale(1./jtptCut->Integral());
	jtptRef->Scale(1./jtptRef->Integral());
	jtptRef->Draw();
	jtptCut->SetLineColor(2);
	jtptCut->Draw("same");


}