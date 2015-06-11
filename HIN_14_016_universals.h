

  //*********************************************************
  // SET INTEGRAL LIMITS AND LABELS
  //*************************************************

  double etalim = 1.0;
  double philim = 1.0;

  TString etarangelabel = "Projected |#Delta#eta|<1.0";
  TString phirangelabel = "Projected |#Delta#phi|<1.0";

  float llimitphi,rlimitphi,llimiteta,rlimiteta,nbins;
  //*********************************************************

Int_t nbounds_phi = 18;
Int_t nbounds_eta = 20;
Int_t nbounds_eta2 = 18;
 


// Double_t bin_bounds_phi[16]
//   =  {-1.50796, -1.00531,-0.879646, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319, 0.879646, 1.00531, 1.50796};

Double_t bin_bounds_phi[18]   =  {-1.50796, -1.00531,-0.879646, -.75398, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,.75398, 0.879646, 1.00531,1.50796};


//Double_t bin_bounds_phi[]   =  {-1.50796,-1.256637, -1.00531,-0.879646, -.75398, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,.75398, 0.879646, 1.00531, 1.256637,1.50796};

  Double_t bin_bounds_eta[20]
    =   {-2.5,-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.,2.5};


 Double_t bin_bounds_eta2[18]
    =   {-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.};

  Double_t pTbins[5] = {1.,2.,3.,4.,8.};

  //------------------------
  //Canvas formatting specs
  //-------------------------

  float textalign = 0.05;
  float textalign2 = 0.23;
  float ts = 0.07;
  float ts2 = 0.065;
  float ts3 = 0.057;
  float tstitle = 0.08;
float tstitle2 = 0.08;
  float xlabeloffset2 = 0.5;
  float xoffset = 0.85;
  float xoffset2 = 1.0;
  float yoffset = 0.85;
  float texty1 = 0.9;
  float texty2 = 0.8;
  float texty3 = 0.7;
  float texty4 = 0.6;
  float legendoffset = 0.65;
  float legendoffset2 = 0.7;

float tspixels = 25;

  //------------------------------
  



 

