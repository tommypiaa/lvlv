#ifndef Ana_ntuple_h
#define Ana_ntuple_h

#include <TFile.h>
#include "TH1.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TRandom.h"
#include "TEventList.h"
#include "TF1.h"
#include "TStyle.h"
#include "TChain.h"
#include <iostream>

class Ana_ntuple{

 private:
  
  float GenHiggs_Eta;
  float GenHiggs_Phi;
  float GenHiggs_Pt;
  float GenHiggs_Mass;

  float GenMET_X;
  float GenMET_Y;
  float GenMET_Pt;

  float GenFJT_Eta[2];
  float GenFJT_Phi[2];
  float GenFJT_Pt[2];
  float GenFJT_Mass[2];
  int   GenFJT_ID[2];

  float GenLepton_Eta[2];
  float GenLepton_Phi[2];
  float GenLepton_Pt[2];
  int   GenLepton_ID[2];

  float GenW_Eta[2];
  float GenW_Phi[2];
  float GenW_Pt[2];
  float GenW_Mass[2];

  float S_MET_X;
  float S_MET_Y;
  float S_MET_Pt;

  int   S_Lepton_Num;
  float S_Lepton_Pt[10];
  float S_Lepton_Eta[10];
  float S_Lepton_Phi[10];
  float S_Lepton_Energy[10];
  int   S_Lepton_ID[10];
  float S_Lepton_Ecal[10];
  float S_Lepton_Hcal[10];
  float S_Lepton_Trk[10];

  float S_FJet_Pt[2];
  float S_FJet_Eta[2];
  float S_FJet_Phi[2];
  float S_FJet_Energy[2];
  float S_FJet_RatioCone02[2];      //[S_FJet_Num]
  float S_FJet_RatioCone04[2];      //[S_FJet_Num]
  float S_FJet_Alpha[2];    //[S_FJet_Num]
  float S_FJet_Beta[2];     //[S_FJet_Num]
  float S_FJet_Bdis_BJet[2]; //[S_FJet_Num] only one combinedSV is used in the analysis ntuple
  
  int   S_CJet_Num;
  float S_CJet_Pt[20];
  float S_CJet_Eta[20];
  float S_CJet_Phi[20];
  float S_CJet_Energy[20];
  float S_CJet_RatioCone02[20];
  float S_CJet_RatioCone04[20];
  float S_CJet_Alpha[20];
  float S_CJet_Beta[20];
  float S_CJet_Bdis_BJet[20];

  TFile * outputfile;
  TTree * tree;
  
 public:
  Ana_ntuple();
  ~Ana_ntuple();
  
  void close();
  void fill();
  
  void new_event(){
    S_Lepton_Num = 0;
    S_CJet_Num   = 0;
  }
  
  void add_G_Higgs(float pt, float eta, float phi,float mass){
    GenHiggs_Eta = eta;
    GenHiggs_Phi = phi;
    GenHiggs_Pt  = pt;
    GenHiggs_Mass= mass;
  }
  

  void add_G_MET(float metx, float mety, float metpt){
    GenMET_X = metx;
    GenMET_Y = mety;
    GenMET_Pt= metpt;
  }
  void add_G_FJT(float pt0, float eta0, float phi0,float mass0, int id0,
		 float pt1, float eta1, float phi1,float mass1, int id1){
    GenFJT_Pt[0]  = pt0;
    GenFJT_Pt[1]  = pt1;
    GenFJT_Eta[0] = eta0;
    GenFJT_Eta[1] = eta1;
    GenFJT_Phi[0] = phi0;
    GenFJT_Phi[1] = phi1;
    GenFJT_Mass[0]= mass0;
    GenFJT_Mass[1]= mass1;
    GenFJT_ID[0]  = id0;
    GenFJT_ID[1]  = id1;
  }
  void add_G_Lep(float pt0, float eta0, float phi0, int id0,
		 float pt1, float eta1, float phi1, int id1){
    GenLepton_Pt[0]  = pt0;
    GenLepton_Eta[0] = eta0;
    GenLepton_Phi[0] = phi0;
    GenLepton_ID[0]  = id0;
    GenLepton_Pt[1]  = pt1;
    GenLepton_Eta[1] = eta1;
    GenLepton_Phi[1] = phi1;
    GenLepton_ID[1]  = id1;
  }
  void add_G_W(float pt0, float eta0, float phi0, float mass0,
	       float pt1, float eta1, float phi1, float mass1){
    GenW_Pt[0]  = pt0;
    GenW_Eta[0] = eta0;
    GenW_Phi[0] = phi0;
    GenW_Mass[0]= mass0;
    GenW_Pt[1]  = pt1;
    GenW_Eta[1] = eta1;
    GenW_Phi[1] = phi1;
    GenW_Mass[1]= mass1;
  }

  void add_S_MET(float metx, float mety, float metpt){
    S_MET_X = metx;
    S_MET_Y = mety;
    S_MET_Pt= metpt;
  }

  void add_S_Lep(float pt,  float eta,  float phi, float energy,  int id,
		 float trk, float ecal, float hcal){
    S_Lepton_Pt[S_Lepton_Num]     = pt;
    S_Lepton_Eta[S_Lepton_Num]    = eta;
    S_Lepton_Phi[S_Lepton_Num]    = phi;
    S_Lepton_Energy[S_Lepton_Num] = energy;
    S_Lepton_ID[S_Lepton_Num]     = id;
    S_Lepton_Ecal[S_Lepton_Num]   = ecal;
    S_Lepton_Hcal[S_Lepton_Num]   = hcal;
    S_Lepton_Trk[S_Lepton_Num]    = trk;
    S_Lepton_Num++;
  }

  void add_S_FJT(float pt0,     float eta0,    float phi0,   float energy0,
		 float cone020, float cone040, float alpha0, float beta0, float btag0,
		 float pt1,     float eta1,    float phi1,   float energy1,
		 float cone021, float cone041, float alpha1, float beta1, float btag1){
    S_FJet_Pt[0]          = pt0;
    S_FJet_Eta[0]         = eta0;
    S_FJet_Phi[0]         = phi0;
    S_FJet_Energy[0]      = energy0;
    S_FJet_RatioCone02[0] = cone020;      
    S_FJet_RatioCone04[0] = cone040;     
    S_FJet_Alpha[0]       = alpha0;   
    S_FJet_Beta[0]        = beta0;    
    S_FJet_Bdis_BJet[0]   = btag0;

    S_FJet_Pt[1]          = pt1;
    S_FJet_Eta[1]         = eta1;
    S_FJet_Phi[1]         = phi1;
    S_FJet_Energy[1]      = energy1;
    S_FJet_RatioCone02[1] = cone021;      
    S_FJet_RatioCone04[1] = cone041;     
    S_FJet_Alpha[1]       = alpha1;   
    S_FJet_Beta[1]        = beta1;    
    S_FJet_Bdis_BJet[1]   = btag1;
  }

  void add_S_CJet(float pt,     float eta,    float phi,   float energy,
		  float cone02, float cone04, float alpha, float beta, float btag){
    S_CJet_Pt[S_CJet_Num]          = pt;
    S_CJet_Eta[S_CJet_Num]         = eta;
    S_CJet_Phi[S_CJet_Num]         = phi;
    S_CJet_Energy[S_CJet_Num]      = energy;
    S_CJet_RatioCone02[S_CJet_Num] = cone02;      
    S_CJet_RatioCone04[S_CJet_Num] = cone04;     
    S_CJet_Alpha[S_CJet_Num]       = alpha;   
    S_CJet_Beta[S_CJet_Num]        = beta;    
    S_CJet_Bdis_BJet[S_CJet_Num]   = btag;
    S_CJet_Num++;
  }

  
};

#endif

#ifdef Ana_ntuple_cxx
#endif 
