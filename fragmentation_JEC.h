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
  TH2D* correction_matrix;
  TFile* correction_file;
  int radius;
  double PF_pt_cut;
  TString algo_corr; 
  int ntrkmax;
  double PF_eta_cut;
 public:
  void reset()
  { 
   correction_matrix=NULL;
   correction_file=NULL;
   ntrkmax=21;
   PF_eta_cut=2.4;
  }
  
  fragmentation_JEC(int radius, double PF_pt_cut)
  {
   reset();
   
   this->radius=radius;
   this->PF_pt_cut=PF_pt_cut;
   if(PF_pt_cut==3) ntrkmax=26;
   else if(PF_pt_cut==2) ntrkmax=31;
   else if(PF_pt_cut==1) ntrkmax=41;
  }
  
  double get_corrected_pt(double jetpt, int ntrk)
  {
   double correction=1;
   
   if(ntrk<ntrkmax) correction=correction_matrix->GetBinContent(correction_matrix->GetXaxis()->FindBin(ntrk),correction_matrix->GetYaxis()->FindBin(jetpt));
   else correction=correction_matrix->GetBinContent(correction_matrix->GetXaxis()->FindBin(ntrkmax),correction_matrix->GetYaxis()->FindBin(jetpt));
   
   return (1/correction)*jetpt;
  }
  
  bool passes_PF_selection(double PF_pt, double PF_eta, double PF_phi, double jet_eta, double jet_phi)
  { 
   if(PF_pt<PF_pt_cut || fabs(PF_eta)>PF_eta_cut) return false;
   
   double r=sqrt(pow(jet_eta-PF_eta,2)+pow(acos(cos(jet_phi-PF_phi)),2));
   
   if(r<((double)(radius)*0.1)) return true;
   else return false;
  }
  
  void set_correction()
  {
    algo_corr=Form("ak%dCalo",radius);
    correction_file = new TFile(Form("corrections/FFJEC_correction_PF_%s_pt%d.root",algo_corr.Data(),(int)PF_pt_cut));
    correction_matrix=(TH2D*)correction_file->Get("pNtrk_pt");
  }
};