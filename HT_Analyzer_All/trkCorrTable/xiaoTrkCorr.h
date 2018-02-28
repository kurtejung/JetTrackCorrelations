

#ifndef TRKCORR_H
#define TRKCORR_H

#ifndef ROOT_TH2F
#include "TH2F.h"
#endif

#ifndef ROOT_TH1F
#include "TH1F.h"
#endif

#include <iostream>

class xiaoTrkCorr{

	public:
		int jhist;
		int khist;
		TFile* corf;
		xiaoTrkCorr(TString fname);
		float getTrkCorr(float pt, float eta, float phi, int hibin);

	private: 
		float ptbins[11];
		int hibins[47];  // we combined last 4 corr. table
		int nhibin;
		int nptbin;
		int hibinWidth;
		float etaw;
		float phiw;
		TH2F* corrTable[50][10];
};

xiaoTrkCorr::xiaoTrkCorr(TString fname){
	hibinWidth = 4;
	nhibin=46;
	nptbin =10;
	etaw = 0.1;
	phiw = 0.12;
	TFile* corf = TFile::Open(fname,"read");

	if( corf->IsOpen()){
		for(int j=0;j<nhibin;j++){
			hibins[j]=j*hibinWidth;
			for(int k=0; k<nptbin;k++){
				corrTable[j][k]= (TH2F*)corf->Get(Form("corr_%d_%d",j,k));
			}
		}
		hibins[nhibin]=200;
	}
	else cout<<"didn't open the correction file!"<<endl;

	ptbins[0] =0.5;
	ptbins[1] =0.7;
	ptbins[2] =1;
	ptbins[3] =2;
	ptbins[4] =3;
	ptbins[5] =4;
	ptbins[6] =8;
	ptbins[7] =12;
	ptbins[8] =16;
	ptbins[9] =20;
	ptbins[10] =300;
}
float xiaoTrkCorr::getTrkCorr(float pt, float eta, float phi, int hibin){
	jhist = -1;
	khist = -1;
	float corr =1;
	for(int j=0;j<nhibin;j++){
		if(hibin >= hibins[j] && hibin< hibins[j+1]) jhist=j;
	}
	if( jhist <0) {
		cout<<"hibin out of the range!"<<endl;
		return 1;
	}
	for(int k =0;k<nptbin;k++){
		if(pt >= ptbins[k] && pt< ptbins[k+1]) khist=k;
	}
	if( khist <0) {
		cout<<"pt out of the range of this correction!"<<endl;
		cout<<pt<<endl;
		return 1;
	}

	corr = corrTable[jhist][khist]->GetBinContent(corrTable[jhist][khist]->FindBin(eta,phi));

	return corr;

}
#endif
