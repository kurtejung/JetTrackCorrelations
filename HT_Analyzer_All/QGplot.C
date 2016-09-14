{
	
	TFile *fin5 = new TFile("QGFrac_pythia6_5TeV.root");
	TFile *fin2 = new TFile("QGFrac_pythia6_2p76TeV.root");

	TH1D *quarkIncl5, *quarkIncl2, *glueIncl5, *glueIncl2;
	TH1D *quarkLead5, *quarkLead2, *glueLead5, *glueLead2;

	TH1D *inclJet5, *leadJet5;
	TH1D *inclJet2, *leadJet2;

	quarkIncl5 = (TH1D*)fin5->Get("quarkFracIncl")->Clone("quarkFracIncl5");
	glueIncl5 = (TH1D*)fin5->Get("glueFracIncl")->Clone("glueFracIncl5");
	quarkLead5 = (TH1D*)fin5->Get("quarkFracLead")->Clone("quarkFracLead5");
	glueLead5 = (TH1D*)fin5->Get("glueFracLead")->Clone("glueFracLeadl5");
	inclJet5 = (TH1D*)fin5->Get("inclJets")->Clone("inclJet5");
	leadJet5 = (TH1D*)fin5->Get("leadJets")->Clone("leadJet5");

	quarkIncl2 = (TH1D*)fin2->Get("quarkFracIncl")->Clone("quarkFracIncl2");
	glueIncl2 = (TH1D*)fin2->Get("glueFracIncl")->Clone("glueFracIncl2");
	quarkLead2 = (TH1D*)fin2->Get("quarkFracLead")->Clone("quarkFracLead2");
	glueLead2 = (TH1D*)fin2->Get("glueFracLead")->Clone("glueFracLeadl2");
	inclJet2 = (TH1D*)fin2->Get("inclJets")->Clone("inclJet2");
	leadJet2 = (TH1D*)fin2->Get("leadJets")->Clone("leadJet2");

	//quarkIncl5->Divide(quarkIncl2);
	//glueIncl5->Divide(glueIncl2);
	//quarkLead5->Divide(quarkLead2);
	//glueLead5->Divide(glueLead2);

	quarkIncl5->Divide(inclJet5);
	glueIncl5->Divide(inclJet5);
	quarkLead5->Divide(leadJet5);
	glueLead5->Divide(leadJet5);

	quarkIncl5->SetXTitle("Jet p_{T} (GeV/c)");
	//quarkIncl5->SetYTitle("5.02 TeV / 2.76 TeV");
	quarkIncl5->SetYTitle("Jet Component Fraction");

	quarkIncl5->SetMaximum(2);
	quarkIncl5->SetMinimum(0);
	quarkIncl5->Draw();
	glueIncl5->Draw("same");
	quarkLead5->Draw("same");
	glueLead5->Draw("same");

	TLegend *l1 = new TLegend(0.1,0.7,0.5,0.9);
	l1->AddEntry(quarkIncl5, "Inclusive-Jet (quark)");
	l1->AddEntry(glueIncl5, "Inclusive-Jet (gluon)");
	l1->AddEntry(quarkLead5, "Leading-Jet (quark)");
	l1->AddEntry(glueLead5, "Leading-Jet (gluon)");

	l1->Draw("Same");



}