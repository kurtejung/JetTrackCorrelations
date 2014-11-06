#include <iostream>
#include <vector>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

struct fragmentation_JEC
{
 private:
  static const int ncent=4;
  TH2D* correction_matrix[ncent];
  TFile* correction_file;
  int radius;
  double PF_pt_cut;
  TString algo_corr; 
  int ntrkmax;
  double PF_eta_cut;
  bool doPH;
  int cent_max[ncent];
  int cent_min[ncent];
 public:
  void reset()
  { 
   for(int icent=0;icent<ncent;icent++){
    correction_matrix[icent]=NULL;
   }
   correction_file=NULL;
   ntrkmax=21;
   PF_eta_cut=2.4;
   cent_min[0]=0;
   cent_max[0]=cent_min[1]=20;
   cent_max[1]=cent_min[2]=60;
   cent_max[2]=cent_min[3]=100;
   cent_max[3]=200;
  }
  
  fragmentation_JEC(int radius, double PF_pt_cut=2,bool doPH=1)
  {
   reset();
   this->doPH=doPH;
   this->radius=radius;
   this->PF_pt_cut=PF_pt_cut;
   if(PF_pt_cut==3) ntrkmax=26;
   else if(PF_pt_cut==2) ntrkmax=31;
   else if(PF_pt_cut==1) ntrkmax=41;
  }
  
  double get_corrected_pt(double jetpt, int ntrk, int cent=0)
  {
   double correction=1;
   
   int cent_bin=0;
   if(doPH){
    for(int icent=0;icent<ncent;icent++){
     if(cent<cent_max[icent] && cent>=cent_min[icent]) cent_bin=icent;
    }
   }
   
   if(ntrk<ntrkmax) correction=correction_matrix[cent_bin]->GetBinContent(correction_matrix[cent_bin]->GetXaxis()->FindBin(ntrk),correction_matrix[cent_bin]->GetYaxis()->FindBin(jetpt));
   else correction=correction_matrix[cent_bin]->GetBinContent(correction_matrix[cent_bin]->GetXaxis()->FindBin(ntrkmax),correction_matrix[cent_bin]->GetYaxis()->FindBin(jetpt));
   
   return (1/correction)*jetpt;
  }
  
  bool passes_PF_selection(double PF_pt, double PF_eta, double PF_phi, int PF_id, double jet_eta, double jet_phi)
  { 
   if(PF_pt<PF_pt_cut || fabs(PF_eta)>PF_eta_cut || PF_id!=1) return false;
   
   double r=sqrt(pow(jet_eta-PF_eta,2)+pow(acos(cos(jet_phi-PF_phi)),2));
   
   if(r<((double)(radius)*0.1)) return true;
   else return false;
  }
  
  void set_correction()
  {
   cout<<"setting correction"<<endl;
    if(doPH){
     algo_corr=Form("akVs%dCalo",radius);
     correction_file = new TFile(Form("corrections_id1_ph/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
     for(int icent=0;icent<ncent;icent++){
	    correction_matrix[icent]=(TH2D*)correction_file->Get(Form("pNtrk_pt%d",icent));
	   } 
	  }else{
     algo_corr=Form("ak%dCalo",radius);
     correction_file = new TFile(Form("corrections_id1/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
     correction_matrix[0]=(TH2D*)correction_file->Get("pNtrk_pt");
	 }
  }
};