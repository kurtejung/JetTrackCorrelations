{
	
	TFile *f1 = new TFile("/mnt/hadoop/store/user/kjung/PbPb_Pythia6Hydjet_10EvtMixed_histoFiles_Reduced_GenJetRecoTrack/PbPb_MC_Histograms/crab_PbPb_Pythia6Hydjet_10EvtMixed_histoFiles_Reduced_GenJetRecoTrack/161127_165409/0000/HydJet15_PbPb_197.root");
	TIter next2(f1->GetListOfKeys());
	TKey *key;

	while ((key = (TKey*)next2())) {
		TClass *cl = gROOT->GetClass(key->GetClassName());
		if (!cl->InheritsFrom("TH1")) continue;
		TH1 *h = (TH1*)key->ReadObj();
		cout << h->GetName() << " " << h->GetEntries() << endl;
	}

}
