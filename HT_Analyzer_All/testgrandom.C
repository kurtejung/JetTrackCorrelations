{

	gRandom->SetSeed(0);

	TH1D *hh = new TH1D("hh","",100,0,1); hh->Sumw2();
	
	for(int i=0; i<1e9; i++){
		hh->Fill(gRandom->Rndm());

	}

	hh->Draw();

}