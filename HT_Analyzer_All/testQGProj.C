{
	
	TFile *fquark = new TFile("root_output/Pythia15GenJet_GenTrack_QuarkJet_pp_5TeV.root");
	TFile *fgluon = new TFile("root_output/Pythia15GenJet_GenTrack_GluonJet_pp_5TeV.root");

	TH2D *quarkSig1 = (TH2D*)fquark->Get("RecoJet_GenTrack_Sube0_hJetTrackSignalBackground_notrkcorrCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2");
	TH2D *quarkSig4 = (TH2D*)fquark->Get("RecoJet_GenTrack_Sube0_hJetTrackSignalBackground_notrkcorrCent0_Cent10_Pt100_Pt300_TrkPt8_TrkPt999");

	TH2D *glueSig1 = (TH2D*)fgluon->Get("RecoJet_GenTrack_Sube0_hJetTrackSignalBackground_notrkcorrCent0_Cent10_Pt100_Pt300_TrkPt1_TrkPt2");
	TH2D *glueSig4 = (TH2D*)fgluon->Get("RecoJet_GenTrack_Sube0_hJetTrackSignalBackground_notrkcorrCent0_Cent10_Pt100_Pt300_TrkPt8_TrkPt999");

	double lowbin = quarkSig1->GetXaxis()->FindBin(-0.3);
	double highbin = quarkSig1->GetXaxis()->FindBin(0.3);

	cout << "scanning from " << lowbin << " to " << highbin << endl;

	TH1D *projQuarkSig1 = quarkSig1->ProjectionY("_pyPhi_quarkSig1",lowbin, highbin,"e");
	TH1D *projQuarkSig4 = quarkSig4->ProjectionY("_pyPhi_quarkSig4",lowbin, highbin,"e");

	TH1D *projGlueSig1 = glueSig1->ProjectionY("_pyPhi_glueSig1",lowbin, highbin,"e");
	TH1D *projGlueSig4 = glueSig4->ProjectionY("_pyPhi_glueSig4",lowbin, highbin,"e");

	int lbin = projQuarkSig1->FindBin(1.2);
	int hbin = projQuarkSig1->FindBin(1.6);

	projGlueSig1->Scale(projQuarkSig1->Integral(lbin,hbin)/projGlueSig1->Integral(lbin,hbin));
	projGlueSig4->Scale(projQuarkSig4->Integral(lbin,hbin)/projGlueSig4->Integral(lbin,hbin));

	TCanvas *c1 = new TCanvas("c1","",1200,600);
	c1->Divide(2,1);
	c1->cd(1);
	TLatex *eta = new TLatex(-1.2,4e-6,"|#Delta#eta|<0.3");
	TLatex *pp = new TLatex(-1.2,4e-6,"5 TeV Pythia");
	//projQuarkSig1->Divide(projGlueSig1);
	projQuarkSig1->SetMaximum(5e-6);
	projQuarkSig1->SetMinimum(0);
	projQuarkSig1->GetXaxis()->SetRangeUser(-1,1);
	projQuarkSig1->SetXTitle("Jet-Track #Delta#phi");
	projQuarkSig1->SetYTitle("Arbitrary Units");
	TLatex *l1 = new TLatex(-1.2,4e-6,"1 < Track p_{T} < 2");
	projQuarkSig1->Draw("hist");
	l1->Draw("same");
	pp->Draw("same");
	projGlueSig1->SetMarkerColor(2);
	projGlueSig1->SetLineColor(2);

	projGlueSig1->Draw("hist,same");

	c1->cd(2);
	//projQuarkSig4->Divide(projGlueSig4);
	projQuarkSig4->SetMaximum(30e-6);
	projQuarkSig4->SetMinimum(0);
	projQuarkSig4->GetXaxis()->SetRangeUser(-1,1);
	projQuarkSig4->SetXTitle("Jet-Track #Delta#phi");
	projQuarkSig4->SetYTitle("Arbitrary Units");
	TLatex *l2 = new TLatex(-1.2,14e-6,"Track p_{T} > 8");
	projQuarkSig4->Draw("hist");
	l2->Draw("same");
	eta->Draw("same");
	projGlueSig4->SetMarkerColor(2);
	projGlueSig4->SetLineColor(2);
	projGlueSig4->Draw("hist,same");


}