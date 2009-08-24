#define  Ana_ntuple_cxx
#include "Ana_ntuple.h"

Ana_ntuple::Ana_ntuple(){

  outputfile = new TFile("lvlv.root","RECREATE");
  tree = new TTree("outtree13","qqH lvlv");
  
  tree->Branch("GenHiggs_Eta",      &GenHiggs_Eta,  "GenHiggs_Eta/F");
  tree->Branch("GenHiggs_Phi",      &GenHiggs_Phi,  "GenHiggs_Phi/F");
  tree->Branch("GenHiggs_Pt",       &GenHiggs_Pt,   "GenHiggs_Pt/F");
  tree->Branch("GenHiggs_Mass",     &GenHiggs_Mass, "GenHiggs_Mass/F");
  tree->Branch("GenMET_X",          &GenMET_X,      "GenMET_X/F");
  tree->Branch("GenMET_Y",          &GenMET_Y,      "GenMET_Y/F");
  tree->Branch("GenMET_Pt",         &GenMET_Pt,     "GenMET_Pt/F");
  tree->Branch("GenFJT_Eta",         GenFJT_Eta,    "GenFJT_Eta[2]/F");
  tree->Branch("GenFJT_Phi",         GenFJT_Phi,    "GenFJT_Phi[2]/F");
  tree->Branch("GenFJT_Pt",          GenFJT_Pt,     "GenFJT_Pt[2]/F");
  tree->Branch("GenFJT_Mass",        GenFJT_Mass,   "GenFJT_Mass[2]/F");
  tree->Branch("GenFJT_ID",          GenFJT_ID,     "GenFJT_ID[2]/I");
  tree->Branch("GenLepton_Eta",      GenLepton_Eta, "GenLepton_Eta[2]/F");
  tree->Branch("GenLepton_Phi",      GenLepton_Phi, "GenLepton_Phi[2]/F");
  tree->Branch("GenLepton_Pt",       GenLepton_Pt,  "GenLepton_Pt[2]/F");
  tree->Branch("GenLepton_ID",       GenLepton_ID,  "GenLepton_ID[2]/I");
  tree->Branch("GenW_Eta",           GenW_Eta,      "GenW_Eta[2]/F");
  tree->Branch("GenW_Phi",           GenW_Phi,      "GenW_Phi[2]/F");
  tree->Branch("GenW_Pt",            GenW_Pt,       "GenW_Pt[2]/F");
  tree->Branch("GenW_Mass",          GenW_Mass,     "GenW_Mass[2]/F");

  tree->Branch("S_MET_X",           &S_MET_X,       "S_MET_X/F");
  tree->Branch("S_MET_Y",           &S_MET_Y,       "S_MET_Y/F");
  tree->Branch("S_MET_Pt",          &S_MET_Pt,      "S_MET_Pt/F");

  tree->Branch("S_Lepton_Num",    &S_Lepton_Num,    "S_Lepton_Num/I");
  tree->Branch("S_Lepton_Pt",      S_Lepton_Pt,     "S_Lepton_Pt[S_Lepton_Num]/F");
  tree->Branch("S_Lepton_Eta",     S_Lepton_Eta,    "S_Lepton_Eta[S_Lepton_Num]/F");
  tree->Branch("S_Lepton_Phi",     S_Lepton_Phi,    "S_Lepton_Phi[S_Lepton_Num]/F");
  tree->Branch("S_Lepton_Energy",  S_Lepton_Energy, "S_Lepton_Energy[S_Lepton_Num]/F");
  tree->Branch("S_Lepton_ID",      S_Lepton_ID,     "S_Lepton_ID[S_Lepton_Num]/I");
  tree->Branch("S_Lepton_Trk",     S_Lepton_Trk,    "S_Lepton_Trk[S_Lepton_Num]/F");
  tree->Branch("S_Lepton_Ecal",    S_Lepton_Ecal,   "S_Lepton_Ecal[S_Lepton_Num]/F");
  tree->Branch("S_Lepton_Hcal",    S_Lepton_Hcal,   "S_Lepton_Hcal[S_Lepton_Num]/F");


  tree->Branch("S_FJet_Pt",          S_FJet_Pt,           "S_FJet_Pt[2]/F");
  tree->Branch("S_FJet_Eta",         S_FJet_Eta,          "S_FJet_Eta[2]/F");
  tree->Branch("S_FJet_Phi",         S_FJet_Phi,          "S_FJet_Phi[2]/F");
  tree->Branch("S_FJet_Energy",      S_FJet_Energy,       "S_FJet_Energy[2]/F");
  tree->Branch("S_FJet_RatioCone02", S_FJet_RatioCone02 , "S_FJet_RatioCone02[2]/F");
  tree->Branch("S_FJet_RatioCone04", S_FJet_RatioCone04 , "S_FJet_RatioCone04[2]/F");
  tree->Branch("S_FJet_Alpha",       S_FJet_Alpha,        "S_FJet_Alpha[2]/F");
  tree->Branch("S_FJet_Beta",        S_FJet_Beta,         "S_FJet_Beta[2]/F");
  tree->Branch("S_FJet_Bdis_BJet",   S_FJet_Bdis_BJet,    "S_FJet_Bdis_BJet[2]/F");

  tree->Branch("S_CJet_Num",        &S_CJet_Num,          "S_CJet_Num/I");
  tree->Branch("S_CJet_Pt",          S_CJet_Pt,           "S_CJet_Pt[S_CJet_Num]/F");
  tree->Branch("S_CJet_Eta",         S_CJet_Eta,          "S_CJet_Eta[S_CJet_Num]/F");
  tree->Branch("S_CJet_Phi",         S_CJet_Phi,          "S_CJet_Phi[S_CJet_Num]/F");
  tree->Branch("S_CJet_Energy",      S_CJet_Energy,       "S_CJet_Energy[S_CJet_Num]/F");
  tree->Branch("S_CJet_RatioCone02", S_CJet_RatioCone02 , "S_CJet_RatioCone02[S_CJet_Num]/F");
  tree->Branch("S_CJet_RatioCone04", S_CJet_RatioCone04 , "S_CJet_RatioCone04[S_CJet_Num]/F");
  tree->Branch("S_CJet_Alpha",       S_CJet_Alpha,        "S_CJet_Alpha[S_CJet_Num]/F");
  tree->Branch("S_CJet_Beta",        S_CJet_Beta,         "S_CJet_Beta[S_CJet_Num]/F");
  tree->Branch("S_CJet_Bdis_BJet",   S_CJet_Bdis_BJet,    "S_CJet_Bdis_BJet[S_CJet_Num]/F");
}

Ana_ntuple::~Ana_ntuple(){
}

void Ana_ntuple::close(){
  
  outputfile->cd();
  
  tree->Write();
  
  outputfile->Close("R");
  
}

void Ana_ntuple::fill(){
  
  tree->Fill();

}
