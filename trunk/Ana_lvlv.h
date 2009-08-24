#ifndef Ana_lvlv_h
#define Ana_lvlv_h

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TMatrixD.h>
#include <TTree.h>
#include <TSystem.h>
#include <TString.h>
#include <TBranch.h>
#include <TMath.h>
#include "Ana_ntuple.h"
#include "Ana_ntuple.cxx"
using namespace std;


class Ana_lvlv{
  
 private:

  int     ParticleNum,ParticleID[200],ParticleST[200];
  float   ParticlePt[200],ParticleEta[200],ParticlePhi[200],ParticleMass[200];
  int     AllParticleNum,AllParticleID[5000],AllParticleST[5000];
  float   AllParticlePt[5000],AllParticleEta[5000],AllParticlePhi[5000],AllParticleMass[5000];
   
 
  float   PmetPx,PmetPy,PmetPt,PmetPhi,PmetSumEt;
  int     PJet_Number;
  float   PJet_Pt[200],PJet_Phi[200],PJet_Eta[200],PJet_Energy[200];


  float   TmetPx,TmetPy,TmetPt,TmetPhi,TmetSumEt;
  int     Jet_Number;
  float   Jet_Pt[200],Jet_Phi[200],Jet_Eta[200],Jet_Energy[200];

  int     TrkJet_Number;
  float   TrkJet_Pt[200],TrkJet_Phi[200],TrkJet_Eta[200],TrkJet_Energy[200];
  float   TrkJet_TrackInEffCor[200],TrkJet_OutConeCor[200],TrkJet_ExternalConeCor[200];

  int     TowerNum;
  float   TowerEta[8000],TowerPhi[8000],TowerHEt[8000],TowerEEt[8000],TowerEt[8000];

  int     Muon_Num;
  int     Muon_ID[100],Muon_Nhits[100];
  float   Muon_Pt[100],Muon_Eta[100],Muon_Phi[100],Muon_Energy[100];
  bool    Muon_isGlobalMuon[100],Muon_isStandAloneMuon[100],Muon_isTrackerMuon[100];
  float   Muon_vtx_x[100],Muon_vtx_y[100],Muon_vtx_z[100];
  float   Muon_qx3[100],Muon_d0[100],Muon_normalizedChi2[100];
  float   Muon_calEnergy_em[100],Muon_calEnergy_emS9[100];
  float   Muon_calEnergy_had[100],Muon_calEnergy_hadS9[100];
  float   Muon_calEnergy_tower[100],Muon_calEnergy_towerS9[100];
  int     Muon_Iso03_nTracks[100],Muon_Iso05_nTracks[100];
  float   Muon_Iso03_emEt[100],Muon_Iso03_hadEt[100],Muon_Iso03_sumPt[100];
  float   Muon_Iso05_emEt[100],Muon_Iso05_hadEt[100],Muon_Iso05_sumPt[100];

  int     Elec_Num;
  float   Elec_Pt[100],Elec_Eta[100],Elec_Phi[100],Elec_Energy[100];
  int     Elec_ID[100],Elec_id_robust[100],Elec_id_tight[100],Elec_id_loose[100];
  float   Elec_superClusterEnergy[100];
  bool    Elec_energyScaleCorrected[100];
  float   Elec_eSuperClusterOverP[100], Elec_eSeedClusterOverPout[100];
  float   Elec_caloEnergyError[100], Elec_hadOverEm[100], Elec_eSeed[100], Elec_trackMomentumOut[100], Elec_trackMomentumAtVtx[100];
  int     Elec_classification[100];
  int     Elec_nvhit[100];
  float   Elec_vtx_x[100], Elec_vtx_y[100], Elec_vtx_z[100];
  float   Elec_trk_calo_pt[100], Elec_trk_calo_eta[100], Elec_trk_calo_phi[100];
  float   Elec_trk_out_pt[100], Elec_trk_out_eta[100], Elec_trk_out_phi[100];
  float   Elec_trk_vtx_pt[100], Elec_trk_vtx_eta[100], Elec_trk_vtx_phi[100];
  float   Elec_trk_calo_x[100], Elec_trk_calo_y[100], Elec_trk_calo_z[100];
  float   Elec_trk_vtx_x[100], Elec_trk_vtx_y[100], Elec_trk_vtx_z[100];
  float   Elec_deltaEtaSeedClusterAtCalo[100], Elec_deltaPhiSeedClusterAtCalo[100];
  float   Elec_deltaEtaSuperClusterAtVtx[100], Elec_deltaPhiSuperClusterAtVtx[100];

  int     Track_Num;
  float   Track_Eta[500],Track_Phi[500],Track_Pt[500],Track_Z[500],Track_d0[500];
  float   Track_d0Er[500],Track_PhiEr[500],Track_ZEr[500],Track_EtaEr[500];
  float   Track_chi2[500],Track_ndof[500],Track_nvhit[500],Track_nlhit[500];

  int     Vtx_Num;
  float   Vtx_X[10],Vtx_Y[10],Vtx_Z[10];
  float   Vtx_Xerror[10],Vtx_Yerror[10],Vtx_Zerror[10];
  int     Vtx_Valid[10],Vtx_TrkSize[10];
  float   Vtx_ndof[10],Vtx_chi2[10];

  int     BJet_Num;
  float   BJet_Pt[100],BJet_Eta[100],BJet_Phi[100];
  float   BJet_dis_HighEffBJet[100];
  float   BJet_dis_HighPurBJet[100];
  float   BJet_dis_SV[100];
  float   BJet_dis_combinedSV[100];
  float   BJet_dis_JetBprob[100];
  float   BJet_dis_Jetprob[100];

  bool    TriggerOn;
  int     HLT_Num;
  bool    HLT_Wasrun[200],HLT_Accept[200],HLT_Error[200];


  TChain  *fChain0;//met and jet
  TChain  *fChain1;//met and jet
  TChain  *fChain2;//met and jet
  TChain  *fChain3;//lep
  TChain  *fChain4;//btag
  TChain  *fChain5;//trigger

  Ana_ntuple  *MyTree;   //mytree defination

  float    FJPt_th;

  int      Run_Trig();
  bool     Run_Gen();
  bool     Run_FJT();
  bool     Run_CJet();
  bool     Run_Lep();

  bool     lepveto(float pt, float eta, float phi);
  bool     lepveto(float pt, float eta, float phi, int type);

  // int      TagFJ_detamax(int *num1); 
  int      TagFJ_mjjmax(int *num1);
  int      TagFJ_ptmax(int *num1);
  
  float    Calc_Alpha(float pt, float eta, float phi);
  float    Calc_Beta(float pt, float eta, float phi);
  float    Calc_RatioCone(float range, float pt, float eta, float phi);
  float    Calc_Bdis(float eta, float phi);
  
  bool     Track_Quality(int num);
  bool     Track_Quality(float pt, float chi2, int ndof, float d0, float dz, int nhits);
  
  float    Calc_Elec_Ecal(float isocone, float pt, float eta, float phi, float calo_eta, float calo_phi);
  float    Calc_Elec_Hcal(float isocone, float pt, float eta, float phi);
  float    Calc_Elec_Trk(float isocone, float pt, float eta, float phi);
  
  bool     Iso_muon(float pt, int nhits, float trkiso, float ecaliso, float hcaliso);
  bool     Iso_elec(float pt, float trkiso, float ecaliso, float hcaliso, float h_e, float e_p);
  
  void     Init();
  void     Init_Trigger();

 public:
  
  Ana_lvlv();
  ~Ana_lvlv();
  int      Cut(Long64_t entry);
  void     Loop();
  void     Show(Long64_t entry = -1);
  
  void     SetTrigger(){
    TriggerOn = true; 
    Init_Trigger();
  }
  void   test(){
    std::cout<<"haha"<<endl;
    std::cout<<TmetPx<<endl;
  }
  
};
#endif
#ifdef  Ana_lvlv_cxx

Ana_lvlv::Ana_lvlv()
{
  Init();
}

Ana_lvlv::~Ana_lvlv()
{
  delete fChain0;
  delete fChain1;
  delete fChain2;
  delete fChain3;
  delete fChain4;
  delete fChain5;
  delete MyTree;
}

void
Ana_lvlv::Init(){

  TriggerOn = false;

  MyTree = new Ana_ntuple();

  fChain0 = new TChain("tree0");
  fChain1 = new TChain("tree1");
  fChain2 = new TChain("tree2");
  fChain3 = new TChain("tree3");
  fChain4 = new TChain("tree4");
  if(TriggerOn) {
    fChain5 = new TChain("tree5");
  }  

  ifstream rootlist("sig.in");
  if (!rootlist) {cout <<"No input, exit..." << endl;}
  TString rootfile;
  while(rootlist>>rootfile)
    {
      fChain0->Add(rootfile);
      fChain1->Add(rootfile);
      fChain2->Add(rootfile);
      fChain3->Add(rootfile);
      fChain4->Add(rootfile);
      if(TriggerOn) {
	fChain5->Add(rootfile);
      }       
      cout << rootfile << " added " <<endl;
    }
  
  fChain0->SetBranchAddress("ParticleNum",    &ParticleNum);
  fChain0->SetBranchAddress("ParticleID",      ParticleID);
  fChain0->SetBranchAddress("ParticleST",      ParticleST);
  fChain0->SetBranchAddress("ParticlePt",      ParticlePt);
  fChain0->SetBranchAddress("ParticleMass",    ParticleMass);
  fChain0->SetBranchAddress("ParticleEta",     ParticleEta);
  fChain0->SetBranchAddress("ParticlePhi",     ParticlePhi);
  fChain0->SetBranchAddress("AllParticleNum", &AllParticleNum);
  fChain0->SetBranchAddress("AllParticleID",   AllParticleID);
  fChain0->SetBranchAddress("AllParticleST",   AllParticleST);
  fChain0->SetBranchAddress("AllParticlePt",   AllParticlePt);
  fChain0->SetBranchAddress("AllParticleMass", AllParticleMass);
  fChain0->SetBranchAddress("AllParticleEta",  AllParticleEta);
  fChain0->SetBranchAddress("AllParticlePhi",  AllParticlePhi);
  
  fChain1->SetBranchAddress("PmetPx",       &PmetPx);
  fChain1->SetBranchAddress("PmetPy",       &PmetPy);
  fChain1->SetBranchAddress("PmetPt",       &PmetPt);
  fChain1->SetBranchAddress("PmetPhi",      &PmetPhi);
  fChain1->SetBranchAddress("PmetSumEt",    &PmetSumEt);
  fChain1->SetBranchAddress("PJet_Number",  &PJet_Number);
  fChain1->SetBranchAddress("PJet_Pt",       PJet_Pt);
  fChain1->SetBranchAddress("PJet_Phi",      PJet_Phi);
  fChain1->SetBranchAddress("PJet_Eta",      PJet_Eta);
  fChain1->SetBranchAddress("PJet_Energy",   PJet_Energy);

  fChain2->SetBranchAddress("TmetPx",       &TmetPx);
  fChain2->SetBranchAddress("TmetPy",       &TmetPy);
  fChain2->SetBranchAddress("TmetPt",       &TmetPt);
  fChain2->SetBranchAddress("TmetPhi",      &TmetPhi);
  fChain2->SetBranchAddress("TmetSumEt",    &TmetSumEt );
  fChain2->SetBranchAddress("Jet_Number",   &Jet_Number);
  fChain2->SetBranchAddress("Jet_Pt",        Jet_Pt );
  fChain2->SetBranchAddress("Jet_Phi",       Jet_Phi);
  fChain2->SetBranchAddress("Jet_Eta",       Jet_Eta);
  fChain2->SetBranchAddress("Jet_Energy",    Jet_Energy);
  fChain2->SetBranchAddress("TrkJet_Number",&TrkJet_Number);
  fChain2->SetBranchAddress("TrkJet_Pt",     TrkJet_Pt);
  fChain2->SetBranchAddress("TrkJet_Phi",    TrkJet_Phi);
  fChain2->SetBranchAddress("TrkJet_Eta",    TrkJet_Eta);
  fChain2->SetBranchAddress("TrkJet_Energy", TrkJet_Energy);
  fChain2->SetBranchAddress("TrkJet_TrackInEffCor",   TrkJet_TrackInEffCor );
  fChain2->SetBranchAddress("TrkJet_OutConeCor",      TrkJet_OutConeCor );
  fChain2->SetBranchAddress("TrkJet_ExternalConeCor", TrkJet_ExternalConeCor );
  fChain2->SetBranchAddress("TowerNum",     &TowerNum);
  fChain2->SetBranchAddress("TowerEta",      TowerEta);
  fChain2->SetBranchAddress("TowerPhi",      TowerPhi);
  fChain2->SetBranchAddress("TowerHEt",      TowerHEt);
  fChain2->SetBranchAddress("TowerEEt",      TowerEEt);
  fChain2->SetBranchAddress("TowerEt",       TowerEt);

  fChain3->SetBranchAddress("Muon_Num",             &Muon_Num);
  fChain3->SetBranchAddress("Muon_Pt",               Muon_Pt);
  fChain3->SetBranchAddress("Muon_Eta",              Muon_Eta);
  fChain3->SetBranchAddress("Muon_Phi",              Muon_Phi);
  fChain3->SetBranchAddress("Muon_Energy",           Muon_Energy);
  fChain3->SetBranchAddress("Muon_ID",               Muon_ID);
  fChain3->SetBranchAddress("Muon_isGlobalMuon",     Muon_isGlobalMuon);
  fChain3->SetBranchAddress("Muon_isStandAloneMuon", Muon_isStandAloneMuon);
  fChain3->SetBranchAddress("Muon_isTrackerMuon",    Muon_isTrackerMuon );
  fChain3->SetBranchAddress("Muon_vtx_x",            Muon_vtx_x );
  fChain3->SetBranchAddress("Muon_vtx_y",            Muon_vtx_y );
  fChain3->SetBranchAddress("Muon_vtx_z",            Muon_vtx_z );
  fChain3->SetBranchAddress("Muon_qx3",              Muon_qx3);
  fChain3->SetBranchAddress("Muon_d0",               Muon_d0);
  fChain3->SetBranchAddress("Muon_normalizedChi2",   Muon_normalizedChi2);
  fChain3->SetBranchAddress("Muon_Nhits",            Muon_Nhits); 
  fChain3->SetBranchAddress("Muon_calEnergy_em",     Muon_calEnergy_em );
  fChain3->SetBranchAddress("Muon_calEnergy_emS9",   Muon_calEnergy_emS9  );
  fChain3->SetBranchAddress("Muon_calEnergy_had",    Muon_calEnergy_had );
  fChain3->SetBranchAddress("Muon_calEnergy_hadS9",  Muon_calEnergy_hadS9 );
  fChain3->SetBranchAddress("Muon_calEnergy_tower",  Muon_calEnergy_tower );
  fChain3->SetBranchAddress("Muon_calEnergy_towerS9",Muon_calEnergy_towerS9 );
  fChain3->SetBranchAddress("Muon_Iso03_emEt",       Muon_Iso03_emEt);
  fChain3->SetBranchAddress("Muon_Iso03_hadEt",      Muon_Iso03_hadEt);
  fChain3->SetBranchAddress("Muon_Iso03_nTracks",    Muon_Iso03_nTracks);
  fChain3->SetBranchAddress("Muon_Iso03_sumPt",      Muon_Iso03_sumPt);
  fChain3->SetBranchAddress("Muon_Iso05_emEt",       Muon_Iso05_emEt);
  fChain3->SetBranchAddress("Muon_Iso05_hadEt",      Muon_Iso05_hadEt);
  fChain3->SetBranchAddress("Muon_Iso05_nTracks",    Muon_Iso05_nTracks);
  fChain3->SetBranchAddress("Muon_Iso05_sumPt",      Muon_Iso05_sumPt);

  fChain3->SetBranchAddress("Elec_Num",      &Elec_Num);
  fChain3->SetBranchAddress("Elec_Pt",        Elec_Pt);
  fChain3->SetBranchAddress("Elec_Eta",       Elec_Eta);
  fChain3->SetBranchAddress("Elec_Phi",       Elec_Phi);
  fChain3->SetBranchAddress("Elec_Energy",    Elec_Energy );
  fChain3->SetBranchAddress("Elec_ID",        Elec_ID);
  fChain3->SetBranchAddress("Elec_id_robust", Elec_id_robust);
  fChain3->SetBranchAddress("Elec_id_tight",  Elec_id_tight);
  fChain3->SetBranchAddress("Elec_id_loose",  Elec_id_loose);
  fChain3->SetBranchAddress("Elec_superClusterEnergy",   Elec_superClusterEnergy );
  fChain3->SetBranchAddress("Elec_energyScaleCorrected", Elec_energyScaleCorrected  );
  fChain3->SetBranchAddress("Elec_eSuperClusterOverP",   Elec_eSuperClusterOverP  );
  fChain3->SetBranchAddress("Elec_eSeedClusterOverPout", Elec_eSeedClusterOverPout  );
  fChain3->SetBranchAddress("Elec_caloEnergyError",      Elec_caloEnergyError  );
  fChain3->SetBranchAddress("Elec_hadOverEm",            Elec_hadOverEm );
  fChain3->SetBranchAddress("Elec_eSeed",                Elec_eSeed  );
  fChain3->SetBranchAddress("Elec_trackMomentumOut",     Elec_trackMomentumOut  );
  fChain3->SetBranchAddress("Elec_trackMomentumAtVtx",   Elec_trackMomentumAtVtx    );
  fChain3->SetBranchAddress("Elec_classification",       Elec_classification  );
  fChain3->SetBranchAddress("Elec_nvhit",                Elec_nvhit  );
  fChain3->SetBranchAddress("Elec_vtx_x",                Elec_vtx_x  );
  fChain3->SetBranchAddress("Elec_vtx_y",                Elec_vtx_y   );
  fChain3->SetBranchAddress("Elec_vtx_z",                Elec_vtx_z  );
  fChain3->SetBranchAddress("Elec_trk_calo_pt",               Elec_trk_calo_pt);
  fChain3->SetBranchAddress("Elec_trk_calo_eta",              Elec_trk_calo_eta);
  fChain3->SetBranchAddress("Elec_trk_calo_phi",              Elec_trk_calo_phi);
  fChain3->SetBranchAddress("Elec_trk_out_pt",                Elec_trk_out_pt );
  fChain3->SetBranchAddress("Elec_trk_out_eta",               Elec_trk_out_eta);
  fChain3->SetBranchAddress("Elec_trk_out_phi",               Elec_trk_out_phi);
  fChain3->SetBranchAddress("Elec_trk_vtx_pt",                Elec_trk_vtx_pt);
  fChain3->SetBranchAddress("Elec_trk_vtx_eta",               Elec_trk_vtx_eta);
  fChain3->SetBranchAddress("Elec_trk_vtx_phi",               Elec_trk_vtx_phi);
  fChain3->SetBranchAddress("Elec_trk_calo_x",                Elec_trk_calo_x );
  fChain3->SetBranchAddress("Elec_trk_calo_y",                Elec_trk_calo_y);
  fChain3->SetBranchAddress("Elec_trk_calo_z",                Elec_trk_calo_z);
  fChain3->SetBranchAddress("Elec_trk_vtx_x",                 Elec_trk_vtx_x );
  fChain3->SetBranchAddress("Elec_trk_vtx_y",                 Elec_trk_vtx_y);
  fChain3->SetBranchAddress("Elec_trk_vtx_z",                 Elec_trk_vtx_z);
  fChain3->SetBranchAddress("Elec_deltaEtaSeedClusterAtCalo", Elec_deltaEtaSeedClusterAtCalo);
  fChain3->SetBranchAddress("Elec_deltaPhiSeedClusterAtCalo", Elec_deltaPhiSeedClusterAtCalo);
  fChain3->SetBranchAddress("Elec_deltaEtaSuperClusterAtVtx", Elec_deltaEtaSuperClusterAtVtx);
  fChain3->SetBranchAddress("Elec_deltaPhiSuperClusterAtVtx", Elec_deltaPhiSuperClusterAtVtx);

  fChain3->SetBranchAddress("Track_Num",     &Track_Num);
  fChain3->SetBranchAddress("Track_Eta",      Track_Eta);
  fChain3->SetBranchAddress("Track_Phi",      Track_Phi);
  fChain3->SetBranchAddress("Track_Pt",       Track_Pt);
  fChain3->SetBranchAddress("Track_Z",        Track_Z);
  fChain3->SetBranchAddress("Track_ZEr",      Track_ZEr);
  fChain3->SetBranchAddress("Track_EtaEr",    Track_EtaEr);
  fChain3->SetBranchAddress("Track_PhiEr",    Track_PhiEr);
  fChain3->SetBranchAddress("Track_d0",       Track_d0);
  fChain3->SetBranchAddress("Track_d0Er",     Track_d0Er);
  fChain3->SetBranchAddress("Track_chi2",     Track_chi2);
  fChain3->SetBranchAddress("Track_ndof",     Track_ndof);
  fChain3->SetBranchAddress("Track_nvhit",    Track_nvhit);
  fChain3->SetBranchAddress("Track_nlhit",    Track_nlhit);

  fChain3->SetBranchAddress("Vtx_Num",       &Vtx_Num);
  fChain3->SetBranchAddress("Vtx_X",          Vtx_X);
  fChain3->SetBranchAddress("Vtx_Y",          Vtx_Y);
  fChain3->SetBranchAddress("Vtx_Z",          Vtx_Z);
  fChain3->SetBranchAddress("Vtx_Xerror",     Vtx_Xerror);
  fChain3->SetBranchAddress("Vtx_Yerror",     Vtx_Yerror);
  fChain3->SetBranchAddress("Vtx_Zerror",     Vtx_Zerror);
  fChain3->SetBranchAddress("Vtx_Valid",      Vtx_Valid);
  fChain3->SetBranchAddress("Vtx_ndof",       Vtx_ndof);
  fChain3->SetBranchAddress("Vtx_TrkSize",    Vtx_TrkSize);
  fChain3->SetBranchAddress("Vtx_chi2",       Vtx_chi2);

  fChain4->SetBranchAddress("BJet_Num",            &BJet_Num);
  fChain4->SetBranchAddress("BJet_Pt",              BJet_Pt);
  fChain4->SetBranchAddress("BJet_Eta",             BJet_Eta);
  fChain4->SetBranchAddress("BJet_Phi",             BJet_Phi);
  fChain4->SetBranchAddress("BJet_dis_HighEffBJet", BJet_dis_HighEffBJet);
  fChain4->SetBranchAddress("BJet_dis_HighPurBJet", BJet_dis_HighPurBJet);
  fChain4->SetBranchAddress("BJet_dis_SV",          BJet_dis_SV);
  fChain4->SetBranchAddress("BJet_dis_combinedSV",  BJet_dis_combinedSV);
  fChain4->SetBranchAddress("BJet_dis_JetBprob",    BJet_dis_JetBprob);
  fChain4->SetBranchAddress("BJet_dis_Jetprob",     BJet_dis_Jetprob);

}
void
Ana_lvlv::Init_Trigger(){

  fChain5 = new TChain("tree5");
  
  ifstream rootlist("sig.in");
  if (!rootlist) {cout <<"No input, exit..." << endl;}
  TString rootfile;
  while(rootlist>>rootfile)
    {
      fChain5->Add(rootfile);
      //cout << rootfile << " added " <<endl;
    }
  fChain5->SetBranchAddress("HLT_Num",     &HLT_Num);
  fChain5->SetBranchAddress("HLT_Wasrun",   HLT_Wasrun);
  fChain5->SetBranchAddress("HLT_Accept",   HLT_Accept);
  fChain5->SetBranchAddress("HLT_Error",    HLT_Error);
}


void 
Ana_lvlv::Show(Long64_t entry)
{
   fChain0->Show(entry);
   fChain1->Show(entry);
   fChain2->Show(entry);
   fChain3->Show(entry);
   fChain4->Show(entry);
   if(TriggerOn){
     fChain5->Show(entry);
   }
}

Int_t 
Ana_lvlv::Cut(Long64_t entry)
{
   return 1;
}


#endif // #ifdef makeclass_cxx
