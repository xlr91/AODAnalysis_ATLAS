#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <vector>
#include <AnaAlgorithm/AnaAlgorithm.h>
#include <TH1.h>
#include <xAODJet/JetContainer.h>
#include <xAODTracking/TrackParticleContainer.h>
#include <xAODTruth/TruthParticleContainer.h>
#include <xAODTruth/TruthVertexContainer.h>

class MyxAODAnalysis : public EL::AnaAlgorithm
{
public:
  // this is a standard algorithm constructor
  MyxAODAnalysis (const std::string& name, ISvcLocator* pSvcLocator);

  // these are the functions inherited from Algorithm
  virtual StatusCode initialize () override;
  virtual StatusCode execute () override;
  virtual StatusCode finalize () override;
  Double_t decaylength(const xAOD::TruthVertex* x1, const xAOD::TruthVertex* x2);
  Float_t calcd0(const xAOD::TruthParticle* truth_p, const xAOD::TrackParticle* track_p);


private:
  // Configuration, and any other types of variables go here.
  //float m_cutValue;
  //TTree *m_myTree;
  //TH1 *m_myHist;
  int m_nonSTOP;
  std::vector<int> vector_test;
  

};

#endif