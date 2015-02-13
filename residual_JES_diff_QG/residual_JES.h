#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"

using namespace std;
class residual_JES
{
 private:
  static const int ncent=4;
  static const double lower_pt_cut=40;
  static const double higher_pt_cut=700;
  int radius;
  int cent_max[ncent];
  int cent_min[ncent];
  TFile *correction_file;
  TF1 * residual_correction_function[ncent];
  TString algo_corr;
  public:
  void reset()
  { 
   for(int icent=0;icent<ncent;icent++){
    residual_correction_function[icent]=NULL;
   }
   correction_file=NULL;
   cent_min[0]=0;
   cent_max[0]=cent_min[1]=20;
   cent_max[1]=cent_min[2]=60;
   cent_max[2]=cent_min[3]=100;
   cent_max[3]=200;
  }
  
  residual_JES(int radius=3)
  {
   reset();
   this->radius=radius;
   algo_corr=Form("akVs%dCalo",radius);
  }
   
  void set_correction()
  {
   correction_file = new TFile(Form("ptdependent_resJES_%s.root",algo_corr.Data()));
   for(int icent=0;icent<ncent;icent++){
	 residual_correction_function[icent]=(TF1*)correction_file->Get(Form("fit_qg_%d",icent));
   } 	
  }
  
  double get_pt_var_down(double jetpt, int cent=0)
  {
   // residual correction to correct for the effects of jet resolution in fragmentation jec with a simple centrality binned fit function
   double jetpt_down=jetpt;
   
   if(jetpt>lower_pt_cut && jetpt<higher_pt_cut){
    int cent_bin=0;
    for(int icent=0;icent<ncent;icent++){
     if(cent<cent_max[icent] && cent>=cent_min[icent]) cent_bin=icent;
    }
   
    jetpt_down=(1-fabs(residual_correction_function[cent_bin]->Eval(jetpt)))*jetpt;
   }	
   return jetpt_down;
  }
  
  double get_pt_var_up(double jetpt, int cent=0)
  {
   // residual correction to correct for the effects of jet resolution in fragmentation jec with a simple centrality binned fit function
   double jetpt_up=jetpt;
   
   if(jetpt>lower_pt_cut && jetpt<higher_pt_cut){
    int cent_bin=0;
    for(int icent=0;icent<ncent;icent++){
     if(cent<cent_max[icent] && cent>=cent_min[icent]) cent_bin=icent;
    }
    
   
    jetpt_up=(1+fabs(residual_correction_function[cent_bin]->Eval(jetpt)))*jetpt;
   }	
   return jetpt_up;
  }
  
};
