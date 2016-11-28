
void fileMerger(){
	
	
	TFile *fout = new TFile("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_FullMerge_AllSteps.root","recreate");
	
	TFile *ff[4];
	TFile *ff[0] = new TFile("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_GenGenReduced.root");
	TFile *ff[1] = new TFile("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_GenRecoReduced.root");
	TFile *ff[2] = new TFile("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_RecoGenReduced.root");
	TFile *ff[3] = new TFile("root_output/PbPb_5TeV_MC_PythiaHydjet_MixHistos_Merged_RecoRecoReduced.root");
	
	fout->cd();
	for(int i=0; i<4; i++){
		
		TIter next(ff[i]->GetListOfKeys());
		TKey *key;
		
		while ((key = (TKey*)next())) {
			TClass *cl = gROOT->GetClass(key->GetClassName());
			if (!cl->InheritsFrom("TH2") && !cl->InheritsFrom("TH1")) continue;
			if(cl->InheritsFrom("TH1")){
				TH1F *h = (TH1F*)key->ReadObj();
				h->Write();
			}
			else{
				TH2D *h2 = (TH2D*)key->ReadObj();
				h2->Write();
			}
		}
	}	
}

