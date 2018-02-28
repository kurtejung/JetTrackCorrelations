
#ifndef ROOT_TH2F
#include  "TH2F.h"
#endif

#ifndef ROOT_TFile
#include "TFile.h"
#endif

class xiaoTrkCorr_pp{
	public: 
		xiaoTrkCorr_pp(TString f);
		int binarySearch(float key, float* arr, int i_max, int i_min );
		float getTrkCorr(float pt, float eta, float phi, int hibin, float rmin);
	public:
		TFile* file;
		TH2F* corrTable[22];
		int nptbin = 22;
		float ptbin[23] = {0.7,0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2,
			2.5, 3, 3.5, 4,5, 6, 7, 8, 10, 12, 16, 20, 50, 999};
};

xiaoTrkCorr_pp::xiaoTrkCorr_pp(TString f){
	file = TFile::Open(f);
	for(int i=0; i<nptbin; ++i){
		corrTable[i]=(TH2F*)file->Get(Form("corr_%d",i));
	}
}

int xiaoTrkCorr_pp::binarySearch(float key, float* arr, int i_max, int i_min ){
		if(key> arr[i_max] ) return -1;
		if(key< arr[i_min] ) return -1;
		int mid = floor(float(i_max +i_min)/2);
		//	cout<<mid<<endl;
		if(mid == i_min ) return mid;
		if( arr[mid]> key) return binarySearch(key, arr, mid, i_min);
		else if( arr[mid] < key) return binarySearch(key, arr, i_max, mid);
		else return mid;
}

float xiaoTrkCorr_pp::getTrkCorr(float pt, float eta, float phi, int cent, float rmin){
	int jpt = binarySearch(pt, ptbin, nptbin,0);
	//cout<<"pt = "<<jpt<<endl;
	if(jpt <0 ) {
//		std::cout<<"error!"<<std::endl;
		std::cout<<"jpt="<<jpt<<std::endl;
		std::cout<<"pt="<<pt<<std::endl;
		return 0;
	}
	//int n1 = corrTable[jpt]->GetXaxis()->FindBin(eta);
	//int n2 = corrTable[jpt]->GetYaxis()->FindBin(phi);
	//return corrTable[jpt]->GetBinContent(n1, n2);
	return corrTable[jpt]->GetBinContent(corrTable[jpt]->FindBin(eta,phi));
}

