#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <vector>
#include <AnaAlgorithm/AnaAlgorithm.h>
#include <TH1.h>
#include <GaudiKernel/ITHistSvc.h>
#include <xAODJet/JetContainer.h>
#include <xAODTracking/TrackParticleContainer.h>
#include <xAODTruth/TruthParticleContainer.h>
#include <xAODTruth/TruthVertexContainer.h>
#include <TEfficiency.h>
#include "TCanvas.h"

class MyxAODAnalysis : public EL::AnaAlgorithm
{
public:
  // this is a standard algorithm constructor
  MyxAODAnalysis (const std::string& name, ISvcLocator* pSvcLocator);

  // these are the functions inherited from Algorithm
  virtual StatusCode initialize () override;
  virtual StatusCode execute () override;
  virtual StatusCode finalize () override;

  //Functions I defined
  Double_t decaylength(const xAOD::TruthVertex* x1, const xAOD::TruthVertex* x2);
  Float_t calcdr(const xAOD::TruthParticle* truth_p, const xAOD::TrackParticle* track_p);
  Float_t truthd0(const xAOD::TruthParticle* truth_p, const xAOD::TruthVertex* truth_v);
  Float_t truthd0(const xAOD::TruthParticle* truth_p);

  //Cuts 
  bool cut1(const xAOD::TruthParticle* truth_p, const xAOD::TrackParticle* track_p, Float_t eta_c, Float_t eta_p);
  bool cut2(Float_t mndq, Float_t cut_c);
  
private:
  // Configuration, and any other types of variables go here.
  //float m_cutValue;
  //TTree *m_myTree;
  //TH1 *m_myHist;
  int m_nonSTOP;
  std::vector<int> vector_test;
  TEfficiency* pEff;
  TCanvas* c1;


  float m_etacut;
  float m_phicut; 
  float m_drcut; 
  float m_dzcut; 

  bool m_offline_read;
  bool m_trigger_read;

  int m_cut;
  

};

#endif