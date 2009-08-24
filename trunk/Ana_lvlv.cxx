#define Ana_lvlv_cxx
#include "Ana_lvlv.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

bool
Ana_lvlv::lepveto(float pt, float eta, float phi){
  return lepveto( pt, eta, phi, 1);
}
bool
Ana_lvlv::lepveto(float pt, float eta, float phi , int type){
  //vetotype = 1 jet, = 2 trk
  bool eveto = false;
  bool muveto= false;
  float vetocone = 0;
  float lepptth  = 0;
  if(type == 1)
    vetocone = 0.5;
  if(type == 2)
    vetocone = 0.015;
  if(type == 1)
    lepptth = 10;
  if(type == 2)
    lepptth = 10;
  TVector3 jet;
  jet.SetPtEtaPhi(pt,eta,phi);
  for(int i=0; i<Elec_Num; i++){
    TVector3 electron;
    electron.SetPtEtaPhi(Elec_Pt[i],Elec_Eta[i],Elec_Phi[i]);
    if( Elec_Pt[i] < lepptth ) continue;
    if( electron.DeltaR(jet) < vetocone )
      eveto = true;
  }
  for(int i=0; i<Muon_Num; i++){
    TVector3 muon;
    muon.SetPtEtaPhi(Muon_Pt[i],Muon_Eta[i],Muon_Phi[i]);
    if( Muon_Pt[i] < lepptth ) continue;
    if( muon.DeltaR(jet) < vetocone )
      muveto = true;
  }
  return (eveto||muveto);
}

int
Ana_lvlv::TagFJ_ptmax(int *num1){
  int    jetanum = -1;
  int    jetbnum = -1;
  float  ptj1 = -100;  
  float  ptj2 = -100;  
  for(int i=0 ; i<TrkJet_Number; i++){
    if( TrkJet_Pt[i] < FJPt_th ) continue;
    if( lepveto( TrkJet_Pt[i],TrkJet_Eta[i],TrkJet_Phi[i])) continue;
    if( ptj1 < TrkJet_Pt[i]){
      ptj1 = TrkJet_Pt[i];
      jetanum = i; 
    }
  }
  for(int j=0; j<TrkJet_Number; j++ ){
    if( TrkJet_Pt[j] < FJPt_th ) continue;
    if( lepveto( TrkJet_Pt[j],TrkJet_Eta[j],TrkJet_Phi[j])) continue;
    if( jetanum == j) continue;
    if( ptj2 < TrkJet_Pt[j]){
      ptj2 = TrkJet_Pt[j];
      jetbnum = j; 
    }
  }
  num1[0] = jetanum;
  num1[1] = jetbnum;
  return 1;
}  
int
Ana_lvlv::TagFJ_mjjmax(int *num1){
  int    jetpnum = -1;
  int    jetnnum = -1;
  float  dmjj = -100;  
  for(int i=0 ; i<TrkJet_Number; i++){
    if( TrkJet_Pt[i] < FJPt_th ) continue;
    if( TrkJet_Eta[i]< 0 )       continue;
    if( lepveto( TrkJet_Pt[i],TrkJet_Eta[i],TrkJet_Phi[i])) continue;
    TLorentzVector jet1;
    jet1.SetPtEtaPhiE(TrkJet_Pt[i],TrkJet_Eta[i],TrkJet_Phi[i],TrkJet_Energy[i]);
    for(int j=0; j<TrkJet_Number; j++ ){
      if( TrkJet_Pt[j] < FJPt_th ) continue;
      if( TrkJet_Eta[j]> 0  )      continue;
      if( lepveto( TrkJet_Pt[j],TrkJet_Eta[j],TrkJet_Phi[j])) continue;
      TLorentzVector jet2;
      jet2.SetPtEtaPhiE(TrkJet_Pt[j],TrkJet_Eta[j],TrkJet_Phi[j],TrkJet_Energy[j]);
      TLorentzVector JJ = jet1 + jet2;
      if(dmjj<JJ.M()){
        dmjj = JJ.M();
        jetpnum = i;
        jetnnum = j;
      }
    }
  }
  num1[0] = jetpnum;
  num1[1] = jetnnum;
  return 1;
}  

float    
Ana_lvlv::Calc_Alpha(float pt, float eta, float phi){
  float alpha = 0;
  float beta  = 0;
  for(int j=0; j<Track_Num;j++)
    {
      if(sqrt((eta-Track_Eta[j])*(eta-Track_Eta[j])+(phi-Track_Phi[j])*(phi-Track_Phi[j]))<0.5){
	beta = beta + Track_Pt[j]/pt;
	if(fabs(Track_Z[j]-Vtx_Z[0])>0.3) continue;
	alpha = alpha + Track_Pt[j]/pt;
      }
    }
  if( beta == 0 ){
    alpha = -10;
  }
  return alpha;
}
float    
Ana_lvlv::Calc_Beta(float pt, float eta, float phi){
  float alpha = 0;
  float beta  = 0;
  for(int j=0; j<Track_Num;j++)
    {
      if(sqrt((eta-Track_Eta[j])*(eta-Track_Eta[j])+(phi-Track_Phi[j])*(phi-Track_Phi[j]))<0.5){
	beta = beta + Track_Pt[j]/pt;
	if(fabs(Track_Z[j]-Vtx_Z[0])>0.3) continue;
	alpha = alpha + Track_Pt[j]/pt;
      }
    }
  if( beta == 0 ){
    beta = -1;
    alpha= -10;
  }
  beta = alpha/beta;
  return beta;
}
float    
Ana_lvlv::Calc_RatioCone(float range, float pt, float eta, float phi){
  float cone = 0;
  for(int k=0; k<TowerNum;k++)
    {
      if(sqrt((eta-TowerEta[k])*(eta-TowerEta[k])+(phi-TowerPhi[k])*(phi-TowerPhi[k]))<range)
	{
	  cone = cone + TowerEt[k]/pt;
	}
    }
  return cone;
}
float    
Ana_lvlv::Calc_Bdis(float eta, float phi){
  float bdis = 0;
  bool b_tag = false;
  for(int n=0; n<BJet_Num; n++){
    if(sqrt((BJet_Eta[n]-eta)*(BJet_Eta[n]-eta)+(BJet_Phi[n]-phi)*(BJet_Phi[n]-phi))< 0.02){
      b_tag = true;
      bdis  = BJet_dis_combinedSV[n];
    }
  }
  if(!b_tag){
    bdis = 0;
  }
  return bdis;
}

bool
Ana_lvlv::Track_Quality(int num){
  return Track_Quality(float(Track_Pt[num]),float(Track_chi2[num]),int(Track_ndof[num]),float(Track_d0[num]),float(Track_Z[num]-Vtx_Z[0]),int(Track_nvhit[num]));
}
bool   
Ana_lvlv::Track_Quality(float pt, float chi2, int ndof, float d0, float dz, int nhits)
{
  if ( pt < 1.0 )  return false;
  if ( nhits < 8 ) return false;
  if ( d0 > 1.0 )  return false;
  if ( chi2/float(ndof) > 10 ) return false;
  if ( dz > 5.00 ) return false;
  return true ;
}

float
Ana_lvlv::Calc_Elec_Ecal(float isocone, float pt, float eta, float phi, float calo_eta, float calo_phi){
  float deposit = 0;
  float vetoeta = 0.04 ;
  float vetocone= 0.07;//veto 5*5
  float hitth   = 0.1;
//   float hit;
//   int nhit;
  if(pt<10) return -1;
  
//   if(fabs(eta)<1.5){ //barrel
//     vetoeta  = 0.04;
//     vetocone = 0.045;
//     hitth    = 0.08;
//   }
//   if(fabs(eta)>1.5){ //endcap
//     vetoeta  = 0.04;
//     vetocone = 0.07;
//     hitth    = 0.3;
//   }
  
  for(int k=0; k<TowerNum; k++) {
    if(TowerEEt[k]<hitth) continue;
    float dr = sqrt((eta-TowerEta[k])*(eta-TowerEta[k])+(phi-TowerPhi[k])*(phi-TowerPhi[k]));
    float vetodr = sqrt((calo_eta-TowerEta[k])*(calo_eta-TowerEta[k])+(calo_phi-TowerPhi[k])*(calo_phi-TowerPhi[k]));
//     if(dr < isocone && fabs(calo_eta-TowerEta[k]) > vetoeta ){
//       if(hit<TowerEEt[k]) {
// 	hit = TowerEEt[k];
// 	nhit= k;
//       }
//     }
    if(dr < isocone && vetodr > vetocone && fabs(calo_eta-TowerEta[k]) > vetoeta)
      deposit+= TowerEEt[k];
  }
  return deposit;
}
float
Ana_lvlv::Calc_Elec_Hcal(float isocone, float pt, float eta, float phi){
  float deposit = 0;
  float vetocone = -0.1;
  if(pt<10) return -1;
  for(int k=0; k<TowerNum; k++) {
    float dr = sqrt((eta-TowerEta[k])*(eta-TowerEta[k])+(phi-TowerPhi[k])*(phi-TowerPhi[k]));
    if(dr < isocone && dr > vetocone)
      deposit+= TowerHEt[k];
  }
  return deposit;
}
float
Ana_lvlv::Calc_Elec_Trk(float isocone, float pt, float eta, float phi){
  float trk = 0;
  float vetocone = 0.015;
  if(pt<10) return -1;
  for(int j=0; j<Track_Num; j++) {
    if(Track_Quality(j)) {
      float dr = sqrt((Track_Eta[j]-eta)*(Track_Eta[j]-eta)+(Track_Phi[j]-phi)*(Track_Phi[j]-phi));
      if((dr < isocone)&&(dr > vetocone))
	trk += Track_Pt[j];  
    }
  }
  return trk;
}

bool     
Ana_lvlv::Iso_muon(float pt, int nhits, float trkiso, float ecaliso, float hcaliso){
  if(nhits <10 ) return false;
  if(pt>25){
    if((trkiso+ecaliso+hcaliso)<5)
      return true;
    else
      return false;
  }
  if(pt<25){
    if((trkiso+ecaliso+hcaliso)<(pt-10)/3)
      return true;
    else
      return false;
  } 
  return false;
}
bool     
Ana_lvlv::Iso_elec(float pt, float trkiso, float ecaliso, float hcaliso, float h_e, float e_p){
  if(pt>25){
    if((trkiso)<5)//+ecalis
      return true;
    else
      return false;
  }
  if(pt<25){
    if((trkiso)<(pt-10)/3)//+ecalis
      return true;
    else
      return false;
  }
  return false;
}


void
Ana_lvlv::Loop(){
  
  FJPt_th = 30;
  
  Long64_t nentries = fChain3->GetEntries();
  
  cout<<nentries<<endl;

  int numtrg_muon = 0;
  int numtrg_elec = 0;
  int numtrg_emu  = 0;
  int numlep_muon = 0;
  int numlep_elec = 0;
  int numlep_emu  = 0;
  int numfjt_muon = 0;
  int numfjt_elec = 0;
  int numfjt_emu  = 0;

  int numfjt = 0;
  int numlep = 0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) { 
    if(jentry%1000==0) 
      cout<<"jentry="<<jentry<<endl;
    
    //cout<<"=============================================="<<endl;
    //////////////////////////////////////////////////////////////////////
    //                           start                                  //
    //////////////////////////////////////////////////////////////////////
    
    fChain0->GetEntry(jentry);
    fChain1->GetEntry(jentry);
    fChain2->GetEntry(jentry);
    fChain3->GetEntry(jentry);
    fChain4->GetEntry(jentry);
    if(TriggerOn) {
      fChain5->GetEntry(jentry);
    }
    //test();
    MyTree->new_event();
    
    Run_Gen();
    Run_CJet();

    bool acc_lep = Run_Lep();
    bool acc_fjt = Run_FJT();

    if(TriggerOn) {
      bool acc_trg_muon = false;
      bool acc_trg_elec = false;
      bool acc_trg_emu  = false;
      int  resulttrg    = Run_Trig();
      if(resulttrg==1 || resulttrg==3)
	acc_trg_muon = true;
      if(resulttrg==2 || resulttrg==3)
	acc_trg_elec  = true;
      if(resulttrg==1 || resulttrg==2 || resulttrg==3)
	acc_trg_emu = true;
      if(acc_trg_muon){
	numtrg_muon++;
	if(acc_lep) {
	  numlep_muon++;
	  if(acc_fjt) {
	    numfjt_muon++;
	  }
	}
      }      
      if(acc_trg_elec){
	numtrg_elec++;
	if(acc_lep) {
	  numlep_elec++;
	  if(acc_fjt) {
	    numfjt_elec++;
	  }
	}
      }
      if(acc_trg_emu){
	numtrg_emu++;
	if(acc_lep) {
	  numlep_emu++;
	  if(acc_fjt) {
	    numfjt_emu++;
	    MyTree->fill();
	  }
	}
      }
    }

    if(!TriggerOn){
      if(acc_lep) {
	numlep++;
	if(acc_fjt) {
	  numfjt++;
	  MyTree->fill();
	}
      }
    }

  }
  
  if(TriggerOn){
    cout<<"total num      trigger num     lepiso num       jet selec num"<<endl;
    cout<<"Muon channel  "<<nentries<<"   "<<numtrg_muon<<"   "<<numlep_muon<<"   "<<numfjt_muon<<endl;
    cout<<"elec channel  "<<nentries<<"   "<<numtrg_elec<<"   "<<numlep_elec<<"   "<<numfjt_elec<<endl;
    cout<<"emu  channel  "<<nentries<<"   "<<numtrg_emu <<"   "<<numlep_emu<<"   "<<numfjt_emu<<endl;
  }
  if(!TriggerOn){
    cout<<"total num      lepiso num       jet selec num"<<endl;
    cout<<nentries<<"   "<<numlep<<"   "<<numfjt<<endl;
  }
  MyTree->close(); 
}

int     
Ana_lvlv::Run_Trig(){
  int trig = 0;
  if(HLT_Accept[82]&&!HLT_Accept[47])
    trig = 1;
  if(HLT_Accept[47]&&!HLT_Accept[82])
    trig = 2;
  if(HLT_Accept[47]&&HLT_Accept[82])
    trig = 3;
  return trig;
}


bool
Ana_lvlv::Run_Gen(){

  int G_H_num  = -1;
  int G_FJA_num= 0;
  int G_FJB_num= 0;
  int G_Wp_num = 0;
  int G_Wm_num = 0;
  int G_Lp_num = -1;
  int G_Lm_num = -1;

  for(int i = 0; i<ParticleNum; i++){
    if( ParticleID[i]==25 && ParticleST[i]==3 ) {
      G_H_num = i;
      G_FJA_num = i-2;
      G_FJB_num = i-1;
    }
    if( ParticleID[i]==24 && ParticleST[i]==2) {
      G_Wp_num = i;
    }
    if( ParticleID[i]==-24 && ParticleST[i]==2) {
      G_Wm_num = i;
    }
  }
  if(G_H_num == -1)
    return true;
  int j=-1;
  if( G_Wp_num == (G_Wm_num+1) ) j = G_Wp_num+1;
  if( G_Wm_num == (G_Wp_num+1) ) j = G_Wm_num+1;
  if(j==-1)
    return true;
  for(int k = 0,l = 0; k<4+l; k++){
    if(((ParticleID[j+k]==-11||ParticleID[j+k]==-13)&&ParticleST[j+k]==1)||ParticleID[j+k]==-15)
      G_Lp_num = j+k;
    if(((ParticleID[j+k]==11||ParticleID[j+k]==13)&&ParticleST[j+k]==1)||ParticleID[j+k]==15)
      G_Lm_num = j+k;
    if(ParticleID[j+k]==22&&ParticleST[j+k]==1)
      l++;
  }
  MyTree->add_G_Higgs(ParticlePt[G_H_num],ParticleEta[G_H_num],ParticlePhi[G_H_num],ParticleMass[G_H_num]);
  MyTree->add_G_FJT(ParticlePt[G_FJA_num],ParticleEta[G_FJA_num],ParticlePhi[G_FJA_num],ParticleMass[G_FJA_num],ParticleID[G_FJA_num],ParticlePt[G_FJB_num],ParticleEta[G_FJB_num],ParticlePhi[G_FJB_num],ParticleMass[G_FJB_num],ParticleID[G_FJB_num]);
  MyTree->add_G_W(ParticlePt[G_Wp_num],ParticleEta[G_Wp_num],ParticlePhi[G_Wp_num],ParticleMass[G_Wp_num],ParticlePt[G_Wm_num],ParticleEta[G_Wm_num],ParticlePhi[G_Wm_num],ParticleMass[G_Wm_num]);
  MyTree->add_G_Lep(ParticlePt[G_Lp_num],ParticleEta[G_Lp_num],ParticlePhi[G_Lp_num],ParticleID[G_Lp_num],ParticlePt[G_Lm_num],ParticleEta[G_Lm_num],ParticlePhi[G_Lm_num],ParticleID[G_Lm_num]);
  MyTree->add_G_MET(PmetPx,PmetPy,PmetPt);
  return true;
}

bool
Ana_lvlv::Run_FJT(){

  int num[2] = {-1,-1};

  float S_FJet_Pt[2];
  float S_FJet_Eta[2];
  float S_FJet_Phi[2];
  float S_FJet_Energy[2];
  float S_FJet_RatioCone02[2];     
  float S_FJet_RatioCone04[2];  
  float S_FJet_Alpha[2];  
  float S_FJet_Beta[2];
  float S_FJet_Bdis_BJet[2]; 

  TagFJ_ptmax( num );
  //TagFJ_detamax( num );

  if( num[0] == -1 || num[1] == -1 ) 
    return false;
  
  S_FJet_Pt[0]     = TrkJet_Pt[num[0]];
  S_FJet_Eta[0]    = TrkJet_Eta[num[0]];
  S_FJet_Phi[0]    = TrkJet_Phi[num[0]];
  S_FJet_Energy[0] = TrkJet_Energy[num[0]];
  S_FJet_Pt[1]     = TrkJet_Pt[num[1]];
  S_FJet_Eta[1]    = TrkJet_Eta[num[1]];
  S_FJet_Phi[1]    = TrkJet_Phi[num[1]];
  S_FJet_Energy[1] = TrkJet_Energy[num[1]];
  
  S_FJet_Alpha[0]     = Calc_Alpha(S_FJet_Pt[0],S_FJet_Eta[0],S_FJet_Phi[0]);
  S_FJet_Alpha[1]     = Calc_Alpha(S_FJet_Pt[1],S_FJet_Eta[1],S_FJet_Phi[1]);
  S_FJet_Beta[0]      = Calc_Beta(S_FJet_Pt[0],S_FJet_Eta[0],S_FJet_Phi[0]);
  S_FJet_Beta[1]      = Calc_Beta(S_FJet_Pt[1],S_FJet_Eta[1],S_FJet_Phi[1]);
  S_FJet_RatioCone02[0] = Calc_RatioCone(0.2, S_FJet_Pt[0],S_FJet_Eta[0],S_FJet_Phi[0]);
  S_FJet_RatioCone02[1] = Calc_RatioCone(0.2, S_FJet_Pt[1],S_FJet_Eta[1],S_FJet_Phi[1]);
  S_FJet_RatioCone04[0] = Calc_RatioCone(0.4, S_FJet_Pt[0],S_FJet_Eta[0],S_FJet_Phi[0]);
  S_FJet_RatioCone04[1] = Calc_RatioCone(0.4, S_FJet_Pt[1],S_FJet_Eta[1],S_FJet_Phi[1]);
  S_FJet_Bdis_BJet[0] = Calc_Bdis(S_FJet_Eta[0],S_FJet_Phi[0]);
  S_FJet_Bdis_BJet[1] = Calc_Bdis(S_FJet_Eta[1],S_FJet_Phi[1]);

  MyTree->add_S_FJT(S_FJet_Pt[0],S_FJet_Eta[0],S_FJet_Phi[0],S_FJet_Energy[0],
		    S_FJet_RatioCone02[0],S_FJet_RatioCone04[0],
		    S_FJet_Alpha[0],S_FJet_Beta[0],S_FJet_Bdis_BJet[0],
		    S_FJet_Pt[1],S_FJet_Eta[1],S_FJet_Phi[1],S_FJet_Energy[1],
		    S_FJet_RatioCone02[0],S_FJet_RatioCone04[0],
		    S_FJet_Alpha[1],S_FJet_Beta[1],S_FJet_Bdis_BJet[1]);
  return true;
}

bool
Ana_lvlv::Run_CJet(){
  
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
  
  int njet=0;
  
  int num[2];
  TagFJ_ptmax( num );
  
  for(int i=0; i<TrkJet_Number; i++){
    if( i == num[0] || i == num[1] ) continue;
    if( TrkJet_Pt[i]<FJPt_th ) continue;
    
    if(lepveto(TrkJet_Pt[i],TrkJet_Eta[i],TrkJet_Phi[i])) continue;
    
    S_CJet_Pt[njet]          = TrkJet_Pt[i];
    S_CJet_Eta[njet]         = TrkJet_Eta[i];
    S_CJet_Phi[njet]         = TrkJet_Phi[i];
    S_CJet_Energy[njet]      = TrkJet_Energy[i];
    S_CJet_RatioCone02[njet] = Calc_RatioCone(0.2,S_CJet_Pt[njet],S_CJet_Eta[njet],S_CJet_Phi[njet]);
    S_CJet_RatioCone04[njet] = Calc_RatioCone(0.4,S_CJet_Pt[njet],S_CJet_Eta[njet],S_CJet_Phi[njet]);
    S_CJet_Alpha[njet]       = Calc_Alpha(S_CJet_Pt[njet],S_CJet_Eta[njet],S_CJet_Phi[njet]);
    S_CJet_Beta[njet]        = Calc_Beta(S_CJet_Pt[njet],S_CJet_Eta[njet],S_CJet_Phi[njet]);
    S_CJet_Bdis_BJet[njet]   = Calc_Bdis(S_CJet_Eta[njet],S_CJet_Phi[njet]);
    
    MyTree->add_S_CJet(S_CJet_Pt[njet], S_CJet_Eta[njet],S_CJet_Phi[njet],S_CJet_Energy[njet],
		       S_CJet_RatioCone02[njet],S_CJet_RatioCone04[njet],
		       S_CJet_Alpha[njet],S_CJet_Beta[njet],S_CJet_Bdis_BJet[njet]);
    njet++;
  }
  
  S_CJet_Num = njet;
  return true;
}

bool
Ana_lvlv::Run_Lep(){

  float S_Lepton_Pt[10];
  float S_Lepton_Eta[10];
  float S_Lepton_Phi[10];
  float S_Lepton_Energy[10];
  int   S_Lepton_ID[10];
  float S_Lepton_Ecal[10];
  float S_Lepton_Hcal[10];
  float S_Lepton_Trk[10];
  
  int nlep = 0;
  float px = TmetPx;
  float py = TmetPy; 

  for(int i=0; i<Muon_Num; i++){
    TLorentzVector mu;
    mu.SetPtEtaPhiM(Muon_Pt[i],Muon_Eta[i],Muon_Phi[i],1.0566);
    px = px-mu.Px();
    py = py-mu.Py();

    if( Muon_Pt[i]<10 ) continue;      
    S_Lepton_Pt[nlep]     = Muon_Pt[i];
    S_Lepton_Eta[nlep]    = Muon_Eta[i];
    S_Lepton_Phi[nlep]    = Muon_Phi[i];
    S_Lepton_Energy[nlep] = Muon_Energy[i];
    S_Lepton_ID[nlep]     = Muon_ID[i];
    S_Lepton_Ecal[nlep]   = Muon_Iso03_emEt[i];
    S_Lepton_Hcal[nlep]   = Muon_Iso03_hadEt[i];
    S_Lepton_Trk[nlep]    = Muon_Iso03_sumPt[i];

    if(Iso_muon(S_Lepton_Pt[nlep], Muon_Nhits[i], S_Lepton_Trk[nlep], S_Lepton_Ecal[nlep], S_Lepton_Hcal[nlep])){ 
      MyTree->add_S_Lep(S_Lepton_Pt[nlep],S_Lepton_Eta[nlep],S_Lepton_Phi[nlep],S_Lepton_Energy[nlep],S_Lepton_ID[nlep],S_Lepton_Trk[nlep], S_Lepton_Ecal[nlep], S_Lepton_Hcal[nlep]);
      nlep++;
    }
  }

  MyTree->add_S_MET(px,py,sqrt(px*px+py*py));

  for(int i=0; i<Elec_Num; i++){
    if( Elec_Pt[i]<10 ) continue; 
    if( Elec_id_robust[i]!=1 ) continue;
    S_Lepton_Pt[nlep]     = Elec_Pt[i];
    S_Lepton_Eta[nlep]    = Elec_Eta[i];
    S_Lepton_Phi[nlep]    = Elec_Phi[i];
    S_Lepton_Energy[nlep] = Elec_Energy[i];
    S_Lepton_ID[nlep]     = Elec_ID[i];
    S_Lepton_Ecal[nlep]   = Calc_Elec_Ecal(0.3, S_Lepton_Pt[nlep], S_Lepton_Eta[nlep], S_Lepton_Phi[nlep], S_Lepton_Eta[nlep], S_Lepton_Phi[nlep]);
    S_Lepton_Hcal[nlep]   = Calc_Elec_Hcal(0.3, S_Lepton_Pt[nlep],S_Lepton_Eta[nlep], S_Lepton_Phi[nlep]);
    S_Lepton_Trk[nlep]    = Calc_Elec_Trk(0.3,  S_Lepton_Pt[nlep],S_Lepton_Eta[nlep], S_Lepton_Phi[nlep]);

    if(Iso_elec(S_Lepton_Pt[nlep], S_Lepton_Trk[nlep], S_Lepton_Ecal[nlep], S_Lepton_Hcal[nlep], Elec_hadOverEm[i], Elec_eSuperClusterOverP[i])){
      MyTree->add_S_Lep(S_Lepton_Pt[nlep],S_Lepton_Eta[nlep],S_Lepton_Phi[nlep],S_Lepton_Energy[nlep],S_Lepton_ID[nlep],S_Lepton_Trk[nlep], S_Lepton_Ecal[nlep], S_Lepton_Hcal[nlep]);      
      nlep++;
    }
  }

  if(nlep!=2)
    return false;
  if(S_Lepton_ID[0]*S_Lepton_ID[1]>0) 
    return false;
  if((S_Lepton_Pt[0]<20)&&(S_Lepton_Pt[1]<20))
    return false;

  return true;
}
