#ifndef MyAnalysis_MyTrigger_H
#define MyAnalysis_MyTrigger_H

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <vector>
#include <TH1.h>
#include <GaudiKernel/ITHistSvc.h>
#include <xAODJet/JetContainer.h>
#include <xAODTracking/TrackParticleContainer.h>
#include <xAODTruth/TruthParticleContainer.h>
#include <xAODTruth/TruthVertexContainer.h>
#include <xAODMuon/MuonContainer.h>
#include <TEfficiency.h>
#include "TCanvas.h"
#include "xAODTrigMuon/L2StandAloneMuon.h"
#include "xAODTrigMuon/L2StandAloneMuonContainer.h"

class MyTrigger : public EL::AnaAlgorithm
{
public:
  // this is a standard algorithm constructor
  MyTrigger (const std::string& name, ISvcLocator* pSvcLocator);

  // these are the functions inherited from Algorithm
  virtual StatusCode initialize () override;
  virtual StatusCode execute () override;
  virtual StatusCode finalize () override;

  //functions i defined
  bool triggercond1(const xAOD::TrackParticle* LRT_track, Float_t t_d0min);
  bool triggercond2(const xAOD::TrackParticle* LRT_track1, const xAOD::TrackParticle* LRT_track2, Float_t t_d0min,  Float_t t_dphimin, Float_t t_dpt);

private:
  float m_d0min; 
  float m_dphimin; 
  float m_dpt; 
  // Configuration, and any other types of variables go here.
  //float m_cutValue;
  //TTree *m_myTree;
  //TH1 *m_myHist;
};

#endif