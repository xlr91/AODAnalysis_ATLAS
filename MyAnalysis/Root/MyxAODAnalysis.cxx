#include <AsgTools/MessageCheck.h>
#include <MyAnalysis/MyxAODAnalysis.h>
#include <xAODEventInfo/EventInfo.h>
#include <cmath>
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TMath.h"

MyxAODAnalysis :: MyxAODAnalysis (const std::string& name,
                                  ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm (name, pSvcLocator) 
    , pEff(0)
    , c1(0)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  This is also where you
  // declare all properties for your algorithm.  Note that things like
  // resetting statistics variables or booking histograms should
  // rather go into the initialize() function.

  // declareProperty( "TitleforJobOption", codetitle = 0, "Desc?");
  declareProperty( "drcut" , m_drcut  = 0.4, "Maximum value for the dr  cut");
  declareProperty( "etacut", m_etacut = 0.5, "Maximum value for the eta cut");
  declareProperty( "phicut", m_phicut = 0.3, "Maximum value for the phi cut");
  declareProperty( "dzcut" , m_dzcut  = 200, "Maximum value for the dz  cut");

  declareProperty( "TriggerRead", m_trigger_read = true, "If it reads the trigger containers");
  declareProperty( "OfflineRead", m_offline_read = true, "If it reads the offline containers");
  declareProperty( "MuonRead",    m_muon_read = true, "If it reads the muon containers");

  declareProperty( "cutnumber", m_cut = 2, "Choose what kind of cut you want");
}

Double_t MyxAODAnalysis::decaylength(const xAOD::TruthVertex* x1, const xAOD::TruthVertex* x2){
  //calculates decay length from a truth particle using their truth vertex (decay length not a method of truth particles)
  double_t result= pow((x1->x() - x2->x()), 2.0) + pow((x1->y() - x2->y()), 2.0) + pow((x1->z() - x2->z()), 2.0);
  return sqrt(result);
}

Float_t MyxAODAnalysis::calcdr(const xAOD::TruthParticle* truth_p, const xAOD::TrackParticle* track_p){
  Float_t dr2 = pow((truth_p->eta() - track_p->eta()), 2.0) + pow((truth_p->phi() - track_p->phi()), 2.0);
  return sqrt(dr2);
}

Float_t MyxAODAnalysis::calcdr(const xAOD::TruthParticle* truth_p, const xAOD::L2StandAloneMuon* standmuon_p){
  Float_t dr2 = pow((truth_p->eta() - standmuon_p->eta()), 2.0) + pow((truth_p->phi() - standmuon_p->phi()), 2.0);
  return sqrt(dr2);
}

Float_t MyxAODAnalysis::truthd0(const xAOD::TruthParticle* truth_p, const xAOD::TruthVertex* truth_v){
  //https://arxiv.org/pdf/1405.6569.pdf page 35
  Float_t ans = ((truth_p -> prodVtx()) -> x() - (truth_v -> x())) * sin(truth_p -> phi()) - ((truth_p -> prodVtx()) -> y()- (truth_v -> y())) * cos(truth_p -> phi());
  return -ans;
}

Float_t MyxAODAnalysis::truthd0(const xAOD::TruthParticle* truth_p){
  Float_t ans = ((truth_p -> prodVtx()) -> x()) * sin(truth_p -> phi()) - ((truth_p -> prodVtx()) -> y()) * cos(truth_p -> phi());
  return -ans;
}


bool MyxAODAnalysis::cut1(const xAOD::TruthParticle* truth_p, const xAOD::TrackParticle* track_p, Float_t eta_c, Float_t eta_p){
  bool ans = true;
  if( abs((truth_p->eta() - track_p->eta())) > eta_c  ||  
      abs((truth_p->phi() - track_p->phi())) > eta_p  ){
    ans = false;
  }
  return ans;
}

//more general cut function
bool MyxAODAnalysis::cut2(Float_t mndq, Float_t cut_c){
  bool ans = true;
  if(mndq > cut_c) ans = false;
  return ans;
}

StatusCode MyxAODAnalysis :: initialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  ANA_MSG_INFO ("in initialize");

  //Monitoring Histograms
  ANA_CHECK(book(TH1F("h_truthDecayLength", "STop_Decay_Length", 100, 0, 10))); //eg Stop
  ANA_CHECK(book(TH1F("h_childDecayLength", "RHadron_Decay_Length", 100, 0, 150))); //eg RHadron
  ANA_CHECK(book(TH1F("h_phiInOffline", "phi_In_Offline", 100, -3.15, -3.15))); 
  ANA_CHECK(book(TH1F("h_etaInOffline", "eta_In_Offline", 100, -5, -5)));
  ///How many times the track passes or fails the cut (1 = pass, 0 = fails)
  ANA_CHECK(book(TH1F("h_offpass", "Passed off cut", 2, 0, 2)));
  ANA_CHECK(book(TH1F("h_ftfpass", "Passed ftf cut", 2, 0, 2)));
  ANA_CHECK(book(TH1F("h_lrtpass", "Passed lrt cut", 2, 0, 2)));


  //Truth Histograms
  ANA_CHECK(book(TH1F("truth_h_dr", "dR values for truth muons", 100, 0, 0.5)));
  ANA_CHECK(book(TH1F("truth_h_d0", "d0 values for truth muons", 300, -200, 200)));
  ANA_CHECK(book(TH1F("truth_h_eta", "eta values for truth muons", 300, -5, 5)));
  ANA_CHECK(book(TH1F("truth_h_pT", "pT values for truth muons", 300, 0, 2)));
  ANA_CHECK(book(TH1F("truth_h_phi", "dphi values for truth tracks", 300, -0.001, 0.001)));
  ANA_CHECK(book(TH1F("truth_h_dphi", "dphi values for truth tracks", 300, -0.001, 0.001)));
  ANA_CHECK(book(TH1F("truth_h_deta", "deta values for truth tracks", 300, -1.2, 1.2)));
  ANA_CHECK(book(TH2F("truth_h_phivTDLength", "Phi vs Transv. DecLengths for offline", 300, 0, 20, 300, -3.2, 3.2)));
  ANA_CHECK(book(TH2F("truth_h_dphivTDLength", "dPhi vs Transv. DecLengths for truth", 300, 0, 20, 300, -0.001, 0.001)));


  //Offline Histograms
  ANA_CHECK(book(TH1F("offline_h_dr", "dR values for offline tracks", 100, 0, 0.15)));
  ANA_CHECK(book(TH1F("offline_h_d0", "d0 values for offline tracks", 300, -30, 30)));
  ANA_CHECK(book(TH1F("offline_h_eta", "eta values for offline tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("offline_h_pT", "pT values for offline tracks", 300, 0, 2)));
  ANA_CHECK(book(TH1F("offline_h_phi", "phi values for offline tracks", 300, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("offline_h_dphi", "dphi values for offline tracks", 300, -0.001, 0.001)));
  ANA_CHECK(book(TH1F("offline_h_deta", "deta values for offline tracks", 300, -1.2, 1.2)));
  ANA_CHECK(book(TH1F("offline_h_dz", "dz values for offline tracks", 300, -1.2, 1.2)));
  ANA_CHECK(book(TH2F("offline_h_phivTDLength", "Phi vs Transv. DecLengths for offline", 300, 0, 20, 300, -3.2, 3.2)));
  ANA_CHECK(book(TH2F("offline_h_dphivTDLength", "dPhi vs Transv. DecLengths for offline", 300, 0, 40, 300, -0.001, 0.001)));


  //FTF Histograms
  ANA_CHECK(book(TH1F("FTF_h_dr", "dR values for FTF tracks", 300, 0, 0.1)));
  ANA_CHECK(book(TH1F("FTF_h_dr_0-10", "dR values for FTF tracks at 0-10 mm TDL", 300, 0, 0.1)));
  ANA_CHECK(book(TH1F("FTF_h_dr_20-30", "dR values for FTF tracks at 20-30 mm TDL", 300, 0, 0.1)));
  ANA_CHECK(book(TH1F("FTF_h_d0", "d0 values for FTF tracks", 300, -100, 100))); //-30/30
  ANA_CHECK(book(TH1F("FTF_h_z0", "z0 values for LRT tracks", 300, -1000, 1000)));
  ANA_CHECK(book(TH1F("FTF_h_eta", "eta values for FTF tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("FTF_h_pT", "pT values for FTF tracks", 300, 0, 2)));
  ANA_CHECK(book(TH1F("FTF_h_phi", "phi values for FTF tracks", 300, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("FTF_h_dphi", "dphi values for FTF tracks", 300, -1, 1)));
  ANA_CHECK(book(TH1F("FTF_h_dphi_0-10", "dphi values for FTF tracks 0-10 mm TDL", 300, -1, 1)));
  ANA_CHECK(book(TH1F("FTF_h_dphi_20-30", "dphi values for FTF tracks 20-30 mm TDL", 300, -1, 1)));
  ANA_CHECK(book(TH1F("FTF_h_deta", "deta values for FTF tracks", 300, -1.2, 1.2)));
  ANA_CHECK(book(TH1F("FTF_h_dz", "dz values for FTF tracks", 300, -300, 300)));
  ANA_CHECK(book(TH1F("FTF_h_dz_0-10", "dz values for FTF tracks at 0-10 mm TDL", 300, -300, 300)));
  ANA_CHECK(book(TH1F("FTF_h_dz_20-30", "dz values for FTF tracks at 20-30 mm TDL", 300, -300, 300)));
  ANA_CHECK(book(TH2F("FTF_h_phivTDLength", "dPhi vs Transv. DecLengths for FTF", 60, 0, 20, 300, -3.2, 3.2)));
  ANA_CHECK(book(TH2F("FTF_h_dphivTDLength", "dPhi vs Transv. DecLengths for FTF", 60, 0, 40, 300, -0.001, 0.001)));
  
  ///pixel/sct hit histograms
  ANA_CHECK(book(TH1F("FTF_h_NPix", "Number of Pixel Hits", 20, 0, 20)));
  ANA_CHECK(book(TH1F("FTF_h_NSct", "Number of SCT Hits", 20, 0, 20)));
  ANA_CHECK(book(TH1F("FTF_h_Nblayer", "Number of Hits in the first layer (B-layer)", 20, 0, 20)));
  ANA_CHECK(book(TH1F("FTF_h_NcontribPix", "Number of Contributing layer of the pixel detector", 20, 0, 20)));
  ANA_CHECK(book(TH2F("FTF_h_NPixvd0", "NPix vs d0 for FTF", 20, -100, 100, 20, 0, 20)));
  ANA_CHECK(book(TH2F("FTF_h_NPixvTDLength", "NPix vs Transv. DecLengths for FTF", 40, 0, 400, 20, 0, 20)));
  ANA_CHECK(book(TH2F("FTF_h_dzvz", "dz vs z FTF", 300, -300, 300, 300, -300, 300)));
  

  //LRT Histograms
  ANA_CHECK(book(TH1F("LRT_h_dr", "dR values for LRT tracks", 300, 0, 0.1)));
  ANA_CHECK(book(TH1F("LRT_h_dr_0-10", "dR values for LRT tracks at 0-10 mm TDL", 300, 0, 0.1)));
  ANA_CHECK(book(TH1F("LRT_h_dr_20-30", "dR values for LRT tracksa at 20-30 mm TDL", 300, 0, 0.1)));
  ANA_CHECK(book(TH1F("LRT_h_d0", "d0 values for LRT tracks", 300, -100, 100)));
  ANA_CHECK(book(TH1F("LRT_h_z0", "z0 values for LRT tracks", 300, -1000, 1000)));
  ANA_CHECK(book(TH1F("LRT_h_eta", "eta values for LRT tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("LRT_h_pT", "pT values for LRT tracks", 300, 0, 2)));
  ANA_CHECK(book(TH1F("LRT_h_phi", "phi values for LRT tracks", 300, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("LRT_h_dphi", "dphi values for LRT tracks", 300, -1, 1)));
  ANA_CHECK(book(TH1F("LRT_h_dphi_0-10", "dphi values for LRT tracks 0-10 mm TDL", 300, -1, 1)));
  ANA_CHECK(book(TH1F("LRT_h_dphi_20-30", "dphi values for LRT tracks 20-30 mm TDL", 300, -1, 1)));
  ANA_CHECK(book(TH1F("LRT_h_deta", "deta values for LRT tracks", 300, -1.2, 1.2)));
  ANA_CHECK(book(TH1F("LRT_h_dz", "dz values for LRT tracks", 300, -300, 300)));
  ANA_CHECK(book(TH1F("LRT_h_dz_0-10", "dz values for LRT tracks at 0-10 mm TDL", 300, -300, 300)));
  ANA_CHECK(book(TH1F("LRT_h_dz_20-30", "dz values for LRT tracks at 20-30 mm TDL", 300, -300, 300)));
  ANA_CHECK(book(TH2F("LRT_h_phivTDLength", "dPhi vs Transv. DecLengths for LRT", 60, 0, 20, 300, -3.2, 3.2)));
  ANA_CHECK(book(TH2F("LRT_h_dphivTDLength", "dPhi vs Transv. DecLengths for LRT", 60, 0, 40, 300, -0.001, 0.001)));
  ANA_CHECK(book(TH2F("LRT_h_drvTDLength", "dR vs Transv. DecLengths for LRT", 60, 0, 400, 120, 0, 0.01))); //print out overflow 
  ANA_CHECK(book(TH2F("LRT_h_d0vTDLength", "d0 vs Transv. DecLengths for LRT", 60, 0, 400, 300, -100, 100)));
  ANA_CHECK(book(TH2F("LRT_hNoCut_drvTDLength", "dR vs Transv. DecLengths for LRT (No Cuts)", 60, 0, 400, 80, 0, 0.1))); 
  
  ///pixel/sct hit histograms
  ANA_CHECK(book(TH1F("LRT_h_NPix", "Number of Pixel Hits", 20, 0, 20)));
  ANA_CHECK(book(TH1F("LRT_h_NSct", "Number of SCT Hits", 20, 0, 20)));
  ANA_CHECK(book(TH1F("LRT_h_NCluster", "Number of SCT+Pixel Hits", 25, 0, 25)));
  ANA_CHECK(book(TH1F("LRT_h_Nblayer", "Number of Hits in the first layer (B-layer)", 20, 0, 20)));
  ANA_CHECK(book(TH1F("LRT_h_NcontribPix", "Number of Contributing layer of the pixel detector", 20, 0, 20)));
  ANA_CHECK(book(TH2F("LRT_h_NPixvd0", "NPix vs d0 for LRT", 20, -100, 100, 20, 0, 20)));
  ANA_CHECK(book(TH2F("LRT_h_NPixvTDLength", "NPix vs Transv. DecLengths for LRT", 40, 0, 400, 20, 0, 20)));
  ANA_CHECK(book(TH2F("LRT_h_NSctvd0", "NSct vs d0 for LRT", 60, -300, 300, 20, 0, 20)));
  ANA_CHECK(book(TH2F("LRT_h_NSctvTDLength", "NSct vs Transv. DecLengths for LRT", 40, 0, 400, 20, 0, 20)));
  ANA_CHECK(book(TH2F("LRT_h_NClustervd0", "NCluster vs d0 for LRT", 100, -300, 300, 25, 0, 25)));
  ANA_CHECK(book(TH2F("LRT_h_NClustervTDLength", "NCluster vs Transv. DecLengths for LRT", 40, 0, 400, 25, 0, 25)));
  ANA_CHECK(book(TH2F("LRT_h_dzvz", "dz vs z LRT", 300, -300, 300, 50, -50, 50)));

  
  //Comparison Histograms
  ANA_CHECK(book(TH1F("compare/h_d0diff_offline", "delta_d0 (Offline-truth)", 100, -2, 2)));
  ANA_CHECK(book(TH1F("compare/h_d0diff_FTF", "delta_d0 (FTF-truth)", 100, -2, 2)));
  ANA_CHECK(book(TH1F("compare/h_d0diff_LRT", "delta_d0 (LRT-truth)", 100, -2, 2)));
  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_offline", "truth_d0_vs_offline_d0", 100, -10, 10, 100, -50, 50)));
  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_FTF", "truth_d0_vs_FTF_d0", 100, -50, 50, 100, -50, 50))); // -10/ 10
  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_LRT", "truth_d0_vs_LRT_d0", 100, -50, 50, 100, -50, 50))); 
  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_All", "truth_d0_vs_All_d0", 100, -10, 10, 100, -50, 50))); 

  
  //Fakes
  ANA_CHECK(book(TH1F("FTF_h_d0fakes_n", "Fakes_function_of_d0_n_FTF", 50, -40, 40)));
  ANA_CHECK(book(TH1F("FTF_h_d0fakes_d", "Fakes_function_of_d0_d_FTF", 50, -40, 40)));
  ANA_CHECK(book(TH1F("FTF_h_d0fakes", "Fakes_function_of_d0_FTF", 50, -40, 40)));

  ANA_CHECK(book(TH1F("FTF_h_etafakes_n", "Fakes_function_of_eta_n_FTF", 50, -3, 3)));
  ANA_CHECK(book(TH1F("FTF_h_etafakes_d", "Fakes_function_of_eta_n_FTF", 50, -3, 3)));
  ANA_CHECK(book(TH1F("FTF_h_etafakes", "Fakes_function_of_eta_FTF", 50, -3, 3)));

  ANA_CHECK(book(TH1F("LRT_h_d0fakes_n", "Fakes_function_of_d0_n_LRT", 50, -300, 300)));
  ANA_CHECK(book(TH1F("LRT_h_d0fakes_d", "Fakes_function_of_d0_d_LRT", 50, -300, 300)));
  ANA_CHECK(book(TH1F("LRT_h_d0fakes", "Fakes_function_of_d0_LRT", 50, -300, 300)));

  ANA_CHECK(book(TH1F("LRT_h_etafakes_n", "Fakes_function_of_eta_n_LRT", 50, -3, 3)));
  ANA_CHECK(book(TH1F("LRT_h_etafakes_d", "Fakes_function_of_eta_d_LRT", 50, -3, 3)));
  ANA_CHECK(book(TH1F("LRT_h_etafakes", "Fakes_function_of_eta_LRT", 50, -3, 3)));


  //Efficiencies
  ///offline respect to truth
  ANA_CHECK(book(TH1F("offl_h_d0eff_n", "Offline Efficiency function of d0 n", 25, -300, 300)));
  ANA_CHECK(book(TH1F("offl_h_d0eff_d", "Offline Efficiency function of d0 d", 25, -300, 300)));
  ANA_CHECK(book(TH1F("offl_h_d0eff", "Offline Efficiency function of d0", 25, -300, 300)));

  ANA_CHECK(book(TH1F("offl_h_etaeff_n", "Offline Efficiency function of eta n", 25, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("offl_h_etaeff_d", "Offline Efficiency function of eta d", 25, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("offl_h_etaeff", "Offline Efficiency function of eta", 25, -3.2, 3.2)));

  ANA_CHECK(book(TH1F("offl_h_pTeff_n", "Offline Efficiency_function_of_pT_n", 25, 0, 150)));
  ANA_CHECK(book(TH1F("offl_h_pTeff_d", "Offline Efficiency_function_of_pT_d", 25, 0, 150)));
  ANA_CHECK(book(TH1F("offl_h_pTeff", "Offline Efficiency function of pT", 25, 0, 150)));

  ANA_CHECK(book(TH1F("offl_h_z0eff_n", "Offline Efficiency_function_of_z0_n", 100, -500, 500)));
  ANA_CHECK(book(TH1F("offl_h_z0eff_d", "Offline Efficiency_function_of_z0_d", 100, -500, 500)));
  ANA_CHECK(book(TH1F("offl_h_z0eff", "Offline Efficiency function of z0", 100, -500, 500)));
  
  ANA_CHECK(book(TH1F("offl_h_TDLeff_n", "Efficiency_function_of_TDL_n", 40, 0, 400)));
  ANA_CHECK(book(TH1F("offl_h_TDLeff_d", "Efficiency_function_of_TDL_d", 40, 0, 400)));
  ANA_CHECK(book(TH1F("offl_h_TDLeff", "Efficiency_function_of_TDL", 40, 0, 400)));

  ///FTF and LRT respect to truth
  ANA_CHECK(book(TH1F("FTF_h_d0eff_n", "Efficiency_function_of_d0_n_FTF", 25, -300, 300)));
  ANA_CHECK(book(TH1F("FTF_h_d0eff_d", "Efficiency_function_of_d0_d_FTF", 25, -300, 300)));
  ANA_CHECK(book(TH1F("FTF_h_d0eff", "Efficiency_function_of_d0_FTF", 25, -300, 300)));
  ANA_CHECK(book(TH1F("LRT_h_d0eff_n", "Efficiency_function_of_d0_n_FTF", 25, -300, 300)));
  ANA_CHECK(book(TH1F("LRT_h_d0eff_d", "Efficiency_function_of_d0_d_FTF", 25, -300, 300)));
  ANA_CHECK(book(TH1F("LRT_h_d0eff", "Efficiency_function_of_d0_FTF", 25, -300, 300)));

  ANA_CHECK(book(TH1F("FTF_h_TDLeff_n", "Efficiency_function_of_d0_n_FTF", 40, 0, 400)));
  ANA_CHECK(book(TH1F("FTF_h_TDLeff_d", "Efficiency_function_of_d0_d_FTF", 40, 0, 400)));
  ANA_CHECK(book(TH1F("FTF_h_TDLeff", "Efficiency_function_of_d0_FTF", 40, 0, 400)));
  ANA_CHECK(book(TH1F("LRT_h_TDLeff_n", "Efficiency_function_of_d0_n_FTF", 40, 0, 400)));
  ANA_CHECK(book(TH1F("LRT_h_TDLeff_d", "Efficiency_function_of_d0_d_FTF", 40, 0, 400)));
  ANA_CHECK(book(TH1F("LRT_h_TDLeff", "Efficiency_function_of_d0_FTF", 40, 0, 400)));

  ///Trigger respect to truth
  ANA_CHECK(book(TH1F("trig_h_d0eff_n", "Efficiency_function_of_d0_n", 25, -50, 50)));
  ANA_CHECK(book(TH1F("trig_h_d0eff_d", "Efficiency_function_of_d0_d", 25, -50, 50)));
  ANA_CHECK(book(TH1F("trig_h_d0eff", "Efficiency_function_of_d0", 25, -50, 50)));

  ANA_CHECK(book(TH1F("trig_h_etaeff_n", "Efficiency_function_of_eta_n", 25, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("trig_h_etaeff_d", "Efficiency_function_of_eta_d", 25, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("trig_h_etaeff", "Efficiency_function_of_eta", 25, -3.2, 3.2)));

  ANA_CHECK(book(TH1F("trig_h_pTeff_n", "Efficiency_function_of_pT_n", 25, 0, 150)));
  ANA_CHECK(book(TH1F("trig_h_pTeff_d", "Efficiency_function_of_pT_d", 25, 0, 150)));
  ANA_CHECK(book(TH1F("trig_h_pTeff", "Efficiency_function_of_pT", 25, 0, 150)));

  ////logarithmic pT values
  const Int_t nbins = 100;
  Double_t xmin = 1e-1;
  Double_t xmax = 1e2;
  Double_t logxmin = TMath::Log10(xmin);
  Double_t logxmax = TMath::Log10(xmax);
  Double_t binwidth = (logxmax-logxmin)/nbins;
  Double_t xbins[nbins+1];
  xbins[0] = xmin;
  for (Int_t i=1;i<=nbins;i++) {
    xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);
  }
  
  ANA_CHECK(book(TH1F("trig_h_pTefflog_n", "Efficiency_function_of_pT_log_n", nbins, xbins)));
  ANA_CHECK(book(TH1F("trig_h_pTefflog_d", "Efficiency_function_of_pT_log_d", nbins, xbins)));
  ANA_CHECK(book(TH1F("trig_h_pTefflog", "Efficiency_function_of_pT_log", nbins, xbins)));

  ANA_CHECK(book(TH1F("trig_h_TDLeff_n", "Efficiency_function_of_TDL_n", 40, 0, 400)));
  ANA_CHECK(book(TH1F("trig_h_TDLeff_d", "Efficiency_function_of_TDL_d", 40, 0, 400)));
  ANA_CHECK(book(TH1F("trig_h_TDLeff", "Efficiency_function_of_TDL", 40, 0, 400)));

  ANA_CHECK(book(TH1F("trig_h_z0eff_n", "Efficiency_function_of_z0_n", 100, -500, 500)));
  ANA_CHECK(book(TH1F("trig_h_z0eff_d", "Efficiency_function_of_z0_d", 100, -500, 500)));
  ANA_CHECK(book(TH1F("trig_h_z0eff", "Efficiency_function_of_z0", 100, -500, 500)));

  ///Trigger respect to offline
  ANA_CHECK(book(TH1F("tgof_h_d0eff_n", "Efficiency_function_of_d0_n", 25, -300, 300)));
  ANA_CHECK(book(TH1F("tgof_h_d0eff_d", "Efficiency_function_of_d0_d", 25, -300, 300)));
  ANA_CHECK(book(TH1F("tgof_h_d0eff", "Efficiency_function_of_d0", 25, -300, 300)));

  ANA_CHECK(book(TH1F("tgof_h_d0FTFeff_n", "Efficiency_function_of_d0_n", 25, -300, 300)));
  ANA_CHECK(book(TH1F("tgof_h_d0FTFeff_d", "Efficiency_function_of_d0_d", 25, -300, 300)));
  ANA_CHECK(book(TH1F("tgof_h_d0FTFeff", "Efficiency_function_of_d0", 25, -300, 300)));

  ANA_CHECK(book(TH1F("tgof_h_etaeff_n", "Efficiency_function_of_eta_n", 25, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("tgof_h_etaeff_d", "Efficiency_function_of_eta_d", 25, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("tgof_h_etaeff", "Efficiency_function_of_eta", 25, -3.2, 3.2))); 

  ANA_CHECK(book(TH1F("tgof_h_pTeff_n", "Efficiency_function_of_pT_n", 25, 0, 150)));
  ANA_CHECK(book(TH1F("tgof_h_pTeff_d", "Efficiency_function_of_pT_d", 25, 0, 150)));
  ANA_CHECK(book(TH1F("tgof_h_pTeff", "Efficiency_function_of_pT", 25, 0, 150)));

  ANA_CHECK(book(TH1F("tgof_h_z0eff_n", "Efficiency_function_of_z0_n", 100, -500, 500)));
  ANA_CHECK(book(TH1F("tgof_h_z0eff_d", "Efficiency_function_of_z0_d", 100, -500, 500)));
  ANA_CHECK(book(TH1F("tgof_h_z0eff", "Efficiency_function_of_z0", 100, -500, 500)));

  ANA_CHECK(book(TH1F("tgof_h_TDLeff_n", "Efficiency_function_of_TDL_n", 40, 0, 400)));
  ANA_CHECK(book(TH1F("tgof_h_TDLeff_d", "Efficiency_function_of_TDL_d", 40, 0, 400)));
  ANA_CHECK(book(TH1F("tgof_h_TDLeff", "Efficiency_function_of_TDL", 40, 0, 400)));

  ANA_CHECK(book(TH1F("tgof_h_TDLFTFeff_n", "Efficiency_function_of_TDL_n", 40, 0, 400)));
  ANA_CHECK(book(TH1F("tgof_h_TDLFTFeff_d", "Efficiency_function_of_TDL_d", 40, 0, 400)));
  ANA_CHECK(book(TH1F("tgof_h_TDLFTFeff", "Efficiency_function_of_TDL", 40, 0, 400)));



  ///Trigger-making plots (in progress):
  ///anything below this line was added rather quickly sorry for the less thought put in the names
  ANA_CHECK(book(TH1F("trigger/h_muonsig_d0max_event", "Maximum abs(d0) in each event for all tracks in truth (signal muons)", 200, 0, 400)));
  ANA_CHECK(book(TH1F("trigger/h_muonprt_d0max_event", "Maximum abs(d0) in each event for all tracks in truth (prompt muons)", 200, 0, 3))); 
  
  ANA_CHECK(book(TH1F("trigger/h_FTF_d0max_event", "Maximum abs(d0) in each event for all tracks in FTF", 200, 0, 400))); 
  ANA_CHECK(book(TH1F("trigger/h_LRT_d0max_event", "Maximum abs(d0) in each event for all tracks in LRT", 200, 0, 400))); 
  ANA_CHECK(book(TH1F("trigger/h_LRT_trigd_d0max_event", "Maximum abs(d0) in each event for LRT tracks in a dR (0.2) cone around muons", 200, 0, 400))); //repurpose this for the max event muon


  ANA_CHECK(book(TH1F("muon/h_d0eff_combined_n", "Efficiency of standalone muon as a function of d0", 25, -100, 100)));
  ANA_CHECK(book(TH1F("muon/h_d0eff_combined_d", "Efficiency of standalone muon as a function of d0", 25, -100, 100)));
  ANA_CHECK(book(TH1F("muon/h_d0eff_combined", "Efficiency of standalone muon as a function of d0", 25, -100, 100)));

  ANA_CHECK(book(TH1F("muon/h_d0eff_standalone_n", "Efficiency of standalone muon as a function of d0", 25, -100, 100)));
  ANA_CHECK(book(TH1F("muon/h_d0eff_standalone_d", "Efficiency of standalone muon as a function of d0", 25, -100, 100)));
  ANA_CHECK(book(TH1F("muon/h_d0eff_standalone", "Efficiency of standalone muon as a function of d0", 25, -100, 100)));

  //t for trigger 
  ANA_CHECK(book(TH1F("muon/t_muoncounthistogram", "number of muons in each RHadron decay", 6, 0, 6)));
  ANA_CHECK(book(TH1F("muon/t_multid0", "d0 values for muons from RHadrons, but those Rhadrons that decay to more than one muon", 50, -50, 50)));
  ANA_CHECK(book(TH1F("muon/t_multideta", "deta values for multi muons from RHadrons, RHadron as ref", 30, -2, 2)));
  ANA_CHECK(book(TH1F("muon/t_multidphi", "dphi values for multi muons from RHadrons, RHadron as ref", 30, -2, 2)));
  ANA_CHECK(book(TH2F("muon/t_multidetavdphi", "deta v dphi values for multi muons from RHadrons, RHadron as ref", 30, -2, 2, 30, -2, 2))); 
  
  ANA_CHECK(book(TH1F("muon/t_doubledeta", "deta values for double muons from RHadrons, between each other", 30, -2, 2))); //should i abs this
  ANA_CHECK(book(TH1F("muon/t_doubledphi", "dphi values for double muons from RHadrons, between each other", 30, -2, 2)));
  ANA_CHECK(book(TH2F("muon/t_doubledetavdphi", "deta v dphi values for multi muons from RHadrons, between each other", 30, -2, 2, 30, -2, 2))); 


  



  
  ANA_MSG_INFO("Offline: " << m_offline_read << " Trigger Read: " <<  m_trigger_read);
  ANA_MSG_INFO("Cut: " << m_cut);
  ANA_MSG_INFO("dr   cut: " << m_drcut);
  ANA_MSG_INFO("deta cut: " << m_etacut);
  ANA_MSG_INFO("dphi cut: " << m_phicut);
  ANA_MSG_INFO("dz   cut: " << m_dzcut);

  return StatusCode::SUCCESS;
}



StatusCode MyxAODAnalysis :: execute ()
{
  //std::vector<const xAOD::TrackParticle*> track_test;
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  ANA_MSG_INFO ("in analaysis execute");
  
  //Access InDetTrackParticles
  const xAOD::TruthParticleContainer* truthparticles = nullptr;
  const xAOD::TrackParticleContainer* offline_particles = nullptr;
  const xAOD::TrackParticleContainer* trig_FTFparticles = nullptr;
  const xAOD::TrackParticleContainer* trig_LRTparticles = nullptr;
  const xAOD::TrackParticleContainer* combinedmuon_container = nullptr;
  const xAOD::L2StandAloneMuonContainer* standalonemuons_container = nullptr;
  
  const xAOD::TrackParticle* matched_offline = nullptr;
  const xAOD::TrackParticle* matched_FTF = nullptr;
  const xAOD::TrackParticle* matched_LRT = nullptr;
  const xAOD::TrackParticle* matched_combinedmuon = nullptr;
  //const xAOD::L2StandAloneMuon* matched_standalonemuon = nullptr;
  
  Float_t mindr;
  Float_t truthd0val;
  Float_t RhTD_Length;

  Bool_t passedflag_ofl = true;
  Bool_t passedflag_ftf = true;
  Bool_t passedflag_lrt = true;
  Bool_t passedflag_muon = true;
  
  uint8_t dummy(-1);
  int NPix;
  int NSct;
  int Nblayer;
  int NcontribPix;
  
  Float_t maxedd0_muon_sig = -1;
  Float_t maxedd0_muon_prt = -1;
  Float_t maxedd0_FTF = -1;
  Float_t maxedd0_LRT = -1;
  Float_t maxedd0_LRT_trigd  = -1;

  int muoncount;


  
  
  ANA_CHECK (evtStore() -> retrieve (truthparticles, "TruthParticles"));
  ANA_MSG_INFO("Found Truth, size is " << truthparticles->size());

  if (m_offline_read) {
    ANA_CHECK (evtStore() -> retrieve (offline_particles, "InDetTrackParticles"));
    ANA_MSG_INFO("Found Offline, size is " << offline_particles->size());
  }
  if (m_trigger_read){
    ANA_CHECK (evtStore() -> retrieve (trig_FTFparticles, "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_FullScan_FTF"));
    ANA_MSG_INFO("Found FTF, size is " << trig_FTFparticles->size());
    ANA_CHECK (evtStore() -> retrieve (trig_LRTparticles, "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_FullScanlrt_FTF"));
    ANA_MSG_INFO("Found LRT, size is " << trig_LRTparticles->size());
  }

  if(m_muon_read){
    ANA_CHECK (evtStore()->retrieve (combinedmuon_container, "HLT_xAOD__TrackParticleContainer_MuonEFInfo_CombTrackParticles"));
    ANA_MSG_INFO("Found muons, size is " << combinedmuon_container -> size());
    ANA_CHECK (evtStore()->retrieve (standalonemuons_container, "HLT_xAOD__L2StandAloneMuonContainer_MuonL2SAInfo"));
    ANA_MSG_INFO("Found standalone muons, size is " << standalonemuons_container -> size());
  }

  //truth loop
  //finds only those that goes STop -> RHadron -> muon
  for (const xAOD::TruthParticle* truth : *truthparticles) {
    if (truth->absPdgId() == 1000006) {     
      ANA_MSG_INFO("TRUTH PID " << truth->pdgId());
      if (truth->nChildren() > 1) {
        for (size_t ichild=0; ichild< truth->nChildren() ; ichild++) {
          const xAOD::TruthParticle* child=truth->child(ichild);
          ANA_MSG_INFO("child PID " << child->pdgId());
          muoncount = 0;

          for (size_t igchild=0; igchild< child->nChildren() ; igchild++) {
            const xAOD::TruthParticle* gchild=child->child(igchild);
            if(gchild == nullptr){ //checks for nullptr as that would break the alg
              ANA_MSG_WARNING("Nullptr alert in gchild");
              continue;
            }

            //ANA_MSG_INFO("gchilds " << gchild -> pdgId());
            //original upgrade study 
            if (gchild->absPdgId() == 13){ 
              muoncount = muoncount + 1;
              //ANA_MSG_INFO("found muons here " << gchild -> pdgId());
              // at this point everything below are muons from RHadrons from Stops
              
              //Get the decay lengths of stop (expected to be 0) and RHadron (varies)
              const xAOD::TruthVertex* tproVtx = truth->prodVtx(); 
              const xAOD::TruthVertex* tdecVtx = truth->decayVtx(); 
              const xAOD::TruthVertex* cproVtx = child->prodVtx(); 
              const xAOD::TruthVertex* cdecVtx = child->decayVtx(); 

              if (tproVtx == nullptr || tdecVtx == nullptr || cproVtx == nullptr || cdecVtx== nullptr){ //another nullptr check grr
                 ANA_MSG_WARNING("Nullptr alert in gchild vector"); 
                 continue;
              }

              truthd0val = truthd0(gchild, cproVtx);
              RhTD_Length = sqrt(pow((cproVtx->x() - cdecVtx->x()), 2.0) + pow((cproVtx->y() - cdecVtx->y()), 2.0)); 

              
              //cuts to truth particle (reasoning is that our trigger tracks wouldnt be able to reconstruct these)
              if(gchild -> pt() < 1000 || abs(gchild -> eta()) > 2.5 || RhTD_Length > 400 || abs(truthd0val) > 300){
                  ANA_MSG_INFO("Truth particle unreconstructable");
                  continue;
              }

              //finds maximum d0 from the muons found in the signal sample (RHadrons)
              if(maxedd0_muon_sig < truthd0val){
                maxedd0_muon_sig = truthd0val; 
              }
              
              hist ("h_truthDecayLength")->Fill (decaylength(tproVtx, tdecVtx));
              hist ("h_childDecayLength")->Fill (decaylength(cproVtx, cdecVtx));
              hist ("truth_h_d0")->Fill (truthd0val);
              hist ("truth_h_eta")->Fill (gchild -> eta());                
              hist ("truth_h_pT")->Fill (gchild -> pt() / 1000);


              //Offline Tracks algorithm //note: would be nice to move this to a separate function 
              if(m_offline_read){

                mindr = 2000;
                matched_offline = nullptr;

                //Track Loop to find matched track 
                for (const xAOD::TrackParticle* offline : *offline_particles){
                  if(mindr > calcdr(gchild, offline)){
                    mindr = calcdr(gchild, offline);
                    matched_offline = offline;
                  }
                } //end of track loop
                

                switch(m_cut) {
                  case 0:
                    passedflag_ofl = true;
                    break;
                  case 1:
                    passedflag_ofl = cut1(gchild, matched_offline, m_etacut, m_phicut);
                    break;
                  case 2:
                    passedflag_ofl = cut2(mindr, m_drcut);
                    break;
                  case 3:
                    passedflag_ofl = cut2(abs(gchild->eta() - matched_offline->eta()), m_etacut);
                    break;
                  case 4:
                    passedflag_ofl = cut2(abs(gchild->phi() - matched_offline->phi()), m_phicut);
                    break;
                  case 5:
                    passedflag_ofl = cut2(abs(cdecVtx->z() - matched_offline->z0()), m_dzcut);
                    break;
                  default:
                    passedflag_ofl = true;
                }
                

                if (passedflag_ofl){
                  hist("h_offpass") -> Fill(1);
                  hist ("offline_h_dr")->Fill (mindr);
                  hist ("offline_h_d0")->Fill (matched_offline -> d0());
                  hist ("offline_h_eta")->Fill (matched_offline -> eta());
                  hist ("offline_h_phi")->Fill (matched_offline -> phi());
                  hist ("offline_h_pT")->Fill (matched_offline -> pt() / 1000);
                  hist ("compare/h_d0diff_offline")->Fill ( matched_offline -> d0() - truthd0val);
                  hist ("compare/h_d0truthvtrack_offline")->Fill (matched_offline -> d0(), truthd0val);

                  hist ("offline_h_dphi") -> Fill ( gchild->phi() - matched_offline->phi() );
                  hist ("offline_h_deta") -> Fill ( gchild->eta() - matched_offline->eta());
                  hist ("offline_h_dz") -> Fill ( cdecVtx->z() - matched_offline->z0());
            
                  hist ("offline_h_dphivTDLength")->Fill (RhTD_Length, gchild->phi() - matched_offline->phi());
                  hist ("offline_h_phivTDLength")->Fill (RhTD_Length, matched_offline->phi());

                  hist ("offl_h_d0eff_n") -> Fill(truthd0val);
                  hist ("offl_h_etaeff_n") -> Fill(gchild->eta());
                  hist ("offl_h_pTeff_n") -> Fill(gchild->pt() / 1000);
                  hist ("offl_h_TDLeff_n") -> Fill(RhTD_Length);
                  hist ("offl_h_z0eff_n") -> Fill(cdecVtx->z());
                } else {
                  hist("h_offpass") -> Fill(0);
                }

                hist ("offl_h_d0eff_d") -> Fill(truthd0val);
                hist ("offl_h_etaeff_d") -> Fill(gchild->eta());
                hist ("offl_h_pTeff_d") -> Fill(gchild->pt() / 1000);
                hist ("offl_h_TDLeff_d") -> Fill(RhTD_Length);
                hist ("offl_h_z0eff_d") -> Fill(cdecVtx->z());

              }

              //Trigger Algorithm
              if(m_trigger_read){
                //FTF First
                mindr = 2000;
                matched_FTF = nullptr;
                //finds matched track
                for (const xAOD::TrackParticle* FTF_T : *trig_FTFparticles){
                  if(mindr > calcdr(gchild, FTF_T)){
                    mindr = calcdr(gchild, FTF_T);
                    matched_FTF = FTF_T;
                  }

                  //fake plotter
                  if(calcdr(gchild, FTF_T) < 0.01){
                    hist ("FTF_h_d0fakes_d") -> Fill(truthd0val);
                    hist ("FTF_h_etafakes_d") -> Fill(gchild -> eta());
                  }

                } //end of FTF loop
    
                hist ("FTF_h_d0fakes_n") -> Fill(truthd0val);
                hist ("FTF_h_etafakes_n") -> Fill(gchild -> eta());

                //Grabbing Pixel and SCT information
                NPix        = matched_FTF->summaryValue( dummy, xAOD::numberOfPixelHits )         ? dummy :-1;
                NSct        = matched_FTF->summaryValue( dummy, xAOD::numberOfSCTHits )           ? dummy :-1;
                Nblayer     = matched_FTF->summaryValue( dummy, xAOD::numberOfBLayerHits )        ? dummy :-1;
                NcontribPix = matched_FTF->summaryValue( dummy, xAOD::numberOfContribPixelLayers )? dummy :-1;


                switch(m_cut) {
                  case 0:
                    passedflag_ftf = true;
                    break;
                  case 1:
                    passedflag_ftf = cut1(gchild, matched_FTF, m_etacut, m_phicut);
                    break;
                  case 2:
                    passedflag_ftf = cut2(mindr, m_drcut);
                    break;
                  case 3:
                    passedflag_ftf = cut2(abs(gchild->eta() - matched_FTF->eta()), m_etacut);
                    break;
                  case 4:
                    passedflag_ftf = cut2(abs(gchild->phi() - matched_FTF->phi()), m_phicut);
                    break;
                  case 5:
                    passedflag_ftf = cut2(abs(cdecVtx->z() - matched_FTF->z0()), m_dzcut);
                    break;
                  default:
                    passedflag_ftf = true;
                }

                if (RhTD_Length < 10){
                  hist("FTF_h_dphi_0-10") -> Fill ( gchild->phi() - matched_FTF->phi());
                  hist("FTF_h_dr_0-10") -> Fill (mindr);
                  hist("FTF_h_dz_0-10") -> Fill ( cdecVtx->z() - matched_FTF->z0());
                } else if (20 < RhTD_Length && RhTD_Length < 30 ){
                  hist("FTF_h_dphi_20-30") -> Fill ( gchild->phi() - matched_FTF->phi());
                  hist("FTF_h_dr_20-30") -> Fill (mindr);
                  hist("FTF_h_dz_20-30") -> Fill ( cdecVtx->z() - matched_FTF->z0());
                }


                if (passedflag_ftf){
                  hist("h_ftfpass") -> Fill(1);
                  hist ("FTF_h_dr")->Fill (mindr);
                  hist ("FTF_h_d0")->Fill (matched_FTF -> d0());
                  hist ("FTF_h_z0")->Fill (matched_FTF -> z0());
                  hist ("FTF_h_eta")->Fill (matched_FTF -> eta());
                  hist ("FTF_h_phi")->Fill (matched_FTF -> phi());
                  hist ("FTF_h_pT")->Fill (matched_FTF -> pt() / 1000);
                  hist ("compare/h_d0diff_FTF")->Fill ( matched_FTF -> d0() - truthd0val);
                  hist ("compare/h_d0truthvtrack_FTF")->Fill (matched_FTF -> d0(), truthd0val);

                  hist ("FTF_h_dphi") -> Fill ( gchild->phi() - matched_FTF->phi() );
                  hist ("FTF_h_deta") -> Fill ( gchild->eta() - matched_FTF->eta());
                  hist ("FTF_h_dz") -> Fill ( cdecVtx->z() - matched_FTF->z0());
                  hist ("FTF_h_dphivTDLength")->Fill (RhTD_Length, gchild->phi() - matched_FTF->phi());

                  hist ("FTF_h_phivTDLength")->Fill (RhTD_Length, matched_FTF->phi());    
                  hist ("FTF_h_d0eff_n") -> Fill(truthd0val);
                  hist ("FTF_h_TDLeff_n") -> Fill(RhTD_Length);

                  hist ("FTF_h_NPix") -> Fill(NPix);
                  hist ("FTF_h_NSct") -> Fill(NSct);
                  hist ("FTF_h_Nblayer") -> Fill(Nblayer);
                  hist ("FTF_h_NcontribPix") -> Fill(NcontribPix);
                  
                  hist ("FTF_h_NPixvd0") -> Fill(truthd0val, NPix);
                  hist ("FTF_h_NPixvTDLength") -> Fill(RhTD_Length, NPix);
                  hist ("FTF_h_dzvz") -> Fill(matched_FTF->z0(), cdecVtx->z() - matched_FTF->z0());
                } else {
                  hist("h_ftfpass") -> Fill(0);
                }

                hist ("FTF_h_d0eff_d") -> Fill(truthd0val);
                hist ("FTF_h_TDLeff_d") -> Fill(RhTD_Length);                

                //LRT
                mindr = 2000;
                matched_LRT = nullptr;
                for (const xAOD::TrackParticle* LRT_T : *trig_LRTparticles){
                  //finds matched track
                  if(mindr > calcdr(gchild, LRT_T)){
                    mindr = calcdr(gchild, LRT_T);
                    matched_LRT = LRT_T;
                  }

                  //fake maker
                  if(calcdr(gchild, LRT_T) < 0.01){
                    hist ("LRT_h_d0fakes_d") -> Fill(truthd0val);
                    hist ("LRT_h_etafakes_d") -> Fill(gchild -> eta());
                  }

                  //Checking for max d0 in the event 
                  //if events are within the dr cone:
                  if(calcdr(gchild, LRT_T) < 0.2){
                    if(abs(LRT_T -> d0()) > maxedd0_LRT_trigd){
                      maxedd0_LRT_trigd = abs(LRT_T -> d0());
                    }
                  }
                } //end of LRT loop

                hist ("LRT_h_d0fakes_n") -> Fill(truthd0val);
                hist ("LRT_h_etafakes_n") -> Fill(gchild -> eta());

                //grab sct/pix informatoin
                NPix        = matched_LRT->summaryValue( dummy, xAOD::numberOfPixelHits )         ? dummy :-1;
                NSct        = matched_LRT->summaryValue( dummy, xAOD::numberOfSCTHits )           ? dummy :-1;
                Nblayer     = matched_LRT->summaryValue( dummy, xAOD::numberOfBLayerHits )        ? dummy :-1;
                NcontribPix = matched_LRT->summaryValue( dummy, xAOD::numberOfContribPixelLayers )? dummy :-1;
                
                switch(m_cut) {
                  case 0:
                    passedflag_lrt = true;
                    break;
                  case 1:
                    passedflag_lrt = cut1(gchild, matched_LRT, m_etacut, m_phicut);
                    break;
                  case 2:
                    passedflag_lrt = cut2(mindr, m_drcut);
                    break;
                  case 3:
                    passedflag_lrt = cut2(abs(gchild->eta() - matched_LRT->eta()), m_etacut);
                    break;
                  case 4:
                    passedflag_lrt = cut2(abs(gchild->phi() - matched_LRT->phi()), m_phicut);
                    break;
                  case 5:
                    passedflag_lrt = cut2(abs(cdecVtx->z() - matched_LRT->z0()), m_dzcut);
                    break;
                  default:
                    passedflag_lrt = true;
                }

                if (RhTD_Length < 10){
                  hist("LRT_h_dphi_0-10") -> Fill ( gchild->phi() - matched_LRT->phi());
                  hist("LRT_h_dr_0-10") -> Fill ( mindr);
                  hist("LRT_h_dz_0-10") -> Fill ( cdecVtx->z() - matched_LRT->z0());
                } else if (20 < RhTD_Length && RhTD_Length < 30 ){
                  hist("LRT_h_dphi_20-30") -> Fill ( gchild->phi() - matched_LRT->phi());
                  hist("LRT_h_dr_20-30") -> Fill ( mindr);
                  hist("LRT_h_dz_20-30") -> Fill ( cdecVtx->z() - matched_LRT->z0());
                }


                if (passedflag_lrt){
                  hist("h_lrtpass") -> Fill(1);
                  hist ("LRT_h_dr")->Fill (mindr);
                  hist ("LRT_h_d0")->Fill (matched_LRT -> d0());
                  hist ("LRT_h_z0")->Fill (matched_LRT -> z0());
                  hist ("LRT_h_eta")->Fill (matched_LRT -> eta());
                  hist ("LRT_h_phi")->Fill (matched_LRT -> phi());
                  hist ("LRT_h_pT")->Fill (matched_LRT -> pt() / 1000);
                  
                  hist ("compare/h_d0diff_LRT")->Fill ( matched_LRT -> d0() - truthd0val);
                  hist ("compare/h_d0truthvtrack_LRT")->Fill (matched_LRT -> d0(), truthd0val);

                  hist ("LRT_h_dphi") -> Fill ( gchild->phi() - matched_LRT->phi() );
                  hist ("LRT_h_deta") -> Fill ( gchild->eta() - matched_LRT->eta());
                  hist ("LRT_h_dz") -> Fill ( cdecVtx->z() - matched_LRT->z0());
                  hist ("LRT_h_dphivTDLength")->Fill (RhTD_Length, gchild->phi() - matched_LRT->phi());
                  hist ("LRT_h_drvTDLength")->Fill (RhTD_Length, mindr);
                  hist ("LRT_h_d0vTDLength")->Fill (RhTD_Length, matched_LRT -> d0());
                  
                  hist ("LRT_h_phivTDLength")->Fill (RhTD_Length, matched_LRT->phi());
                  hist ("LRT_h_d0eff_n") -> Fill(truthd0val);
                  hist ("LRT_h_TDLeff_n") -> Fill(RhTD_Length);

                  hist ("LRT_h_NPix") -> Fill(NPix);
                  hist ("LRT_h_NSct") -> Fill(NSct);
                  hist ("LRT_h_NCluster") -> Fill(NPix + NSct);
                  hist ("LRT_h_Nblayer") -> Fill(Nblayer);
                  hist ("LRT_h_NcontribPix") -> Fill(NcontribPix);
                  
                  hist ("LRT_h_NPixvd0") -> Fill(truthd0val, NPix);
                  hist ("LRT_h_NPixvTDLength") -> Fill(RhTD_Length, NPix);
                  hist ("LRT_h_NSctvd0") -> Fill(truthd0val, NSct);
                  hist ("LRT_h_NSctvTDLength") -> Fill(RhTD_Length, NSct);
                  hist ("LRT_h_NClustervd0") -> Fill(truthd0val, NPix + NSct);
                  hist ("LRT_h_NClustervTDLength") -> Fill(RhTD_Length, NPix + NSct);                                    
                  hist ("LRT_h_dzvz") -> Fill(matched_LRT->z0(), cdecVtx->z() - matched_LRT->z0());
                } else{hist("h_lrtpass") -> Fill(0);}

                hist ("LRT_h_d0eff_d") -> Fill(truthd0val);
                hist ("LRT_h_TDLeff_d") -> Fill(RhTD_Length);
                hist ("LRT_hNoCut_drvTDLength")->Fill (RhTD_Length, mindr);




                //Efficiency plots for trigger 
                //Plots it if is passes either the FTF or LRT cuts
                if(passedflag_ftf || passedflag_lrt){
                  hist ("trig_h_d0eff_n") -> Fill(truthd0val);
                  hist ("trig_h_etaeff_n") -> Fill(gchild->eta());
                  hist ("trig_h_pTeff_n") -> Fill(gchild->pt() / 1000);
                  hist ("trig_h_pTefflog_n") -> Fill(gchild->pt() / 1000);
                  hist ("trig_h_TDLeff_n") -> Fill(RhTD_Length);
                  hist ("trig_h_z0eff_n") -> Fill(cdecVtx->z());
                }

                hist ("trig_h_d0eff_d") -> Fill(truthd0val);
                hist ("trig_h_etaeff_d") -> Fill(gchild->eta());
                hist ("trig_h_pTeff_d") -> Fill(gchild->pt() / 1000);
                hist ("trig_h_pTefflog_d") -> Fill(gchild->pt() / 1000);
                hist ("trig_h_TDLeff_d") -> Fill(RhTD_Length);
                hist ("trig_h_z0eff_d") -> Fill(cdecVtx->z());
                
              }

              //trigger with respect to offline
              if (m_offline_read && m_trigger_read){
                ///at this point it found the matched offline and matched trigger
                ///so fill in the denominator if the current offline track passed thru the cut 

                if(passedflag_ofl) {
                  hist("tgof_h_d0eff_d") -> Fill(matched_offline -> d0());
                  hist("tgof_h_d0FTFeff_d") -> Fill(matched_offline -> d0());
                  hist("tgof_h_etaeff_d") -> Fill(matched_offline->eta());
                  hist("tgof_h_pTeff_d") -> Fill(matched_offline->pt() / 1000);
                  hist("tgof_h_TDLeff_d") -> Fill(RhTD_Length);
                  hist("tgof_h_TDLFTFeff_d") -> Fill(RhTD_Length);
                  hist("tgof_h_z0eff_d") -> Fill(matched_offline->z0());

                  //if all three matches then fill it in
                  if(passedflag_ftf || passedflag_lrt) {
                    hist("tgof_h_d0eff_n") -> Fill(matched_offline -> d0());
                    hist("tgof_h_etaeff_n") -> Fill(matched_offline->eta());
                    hist("tgof_h_pTeff_n") -> Fill(matched_offline->pt() / 1000);
                    hist("tgof_h_TDLeff_n") -> Fill(RhTD_Length);
                    hist("tgof_h_z0eff_n") -> Fill(matched_offline->z0());

                    if(passedflag_ftf){
                      hist("tgof_h_d0FTFeff_n") -> Fill(matched_offline -> d0());
                      hist("tgof_h_TDLFTFeff_n") -> Fill(RhTD_Length);
                    }
                  }


                }
              }

              //muon reading
              //still in progress
              if(m_muon_read){  
                ///combined muon              
                mindr = 2000;
                matched_combinedmuon = nullptr;
                //finds matched track
                for (const xAOD::TrackParticle* muonpar : *combinedmuon_container){
                  if(mindr > calcdr(gchild, muonpar)){
                    mindr = calcdr(gchild, muonpar);
                    matched_combinedmuon = muonpar;
                  }
                } //end of track loop

                switch(m_cut) {
                  case 0:
                    passedflag_muon = true;
                    break;
                  case 1:
                    passedflag_muon = cut1(gchild, matched_combinedmuon, m_etacut, m_phicut);
                    break;
                  case 2:
                    passedflag_muon = cut2(mindr, m_drcut);
                    break;
                  case 3:
                    passedflag_muon = cut2(abs(gchild->eta() - matched_combinedmuon->eta()), m_etacut);
                    break;
                  case 4:
                    passedflag_muon = cut2(abs(gchild->phi() - matched_combinedmuon->phi()), m_phicut);
                    break;
                  case 5:
                    passedflag_muon = cut2(abs(cdecVtx->z() - matched_combinedmuon->z0()), m_dzcut);
                    break;
                  default:
                    passedflag_ofl = true;
                }

                if (passedflag_muon){
                  hist("muon/h_d0eff_combined_n") -> Fill(truthd0val);
                }
                hist("muon/h_d0eff_combined_d") -> Fill(truthd0val);

                ///standalone muon
                mindr = 2000;
                //matched_standalonemuon = nullptr;
                //find matched track
                for (const xAOD::L2StandAloneMuon* standalonemuon : *standalonemuons_container){
                  if(mindr > calcdr(gchild, standalonemuon)){
                    mindr = calcdr(gchild, standalonemuon);
                    //matched_standalonemuon = standalonemuon;
                  }
                }

                //ANA_MSG_INFO(mindr);
                switch(m_cut) {
                  case 0:
                    passedflag_muon = true;
                    break;
                  case 2:
                    passedflag_muon = cut2(mindr, m_drcut);
                    break;
                  default:
                    passedflag_ofl = true;
                }

                if (passedflag_muon){
                  hist("muon/h_d0eff_standalone_n") -> Fill(truthd0val);
                }
                hist("muon/h_d0eff_standalone_d") -> Fill(truthd0val);

              }
            }
          }// end of gchild loop

          //This part of the code is from the last week of the internship
          //studying possible trigger conditions
          ANA_MSG_INFO("muoncount " << muoncount);
          if (muoncount > 2){
            for (size_t igchild=0; igchild< child->nChildren() ; igchild++) {
              const xAOD::TruthParticle* gchild=child->child(igchild);
              if(gchild == nullptr){ //checks for nullptr as that would break the alg
                ANA_MSG_WARNING("Nullptr alert in gchild");
                continue;
              }
              if (gchild->absPdgId() == 13){
                if(gchild -> pt() < 1000) continue;
                hist("muon/t_multid0") -> Fill(truthd0(gchild, child->prodVtx()));
                hist("muon/t_multideta") -> Fill(child -> eta() - gchild -> eta());
                hist("muon/t_multidphi") -> Fill(child -> phi() - gchild ->phi());
                ANA_MSG_INFO(child -> phi() - gchild ->phi());
                hist("muon/t_multidetavdphi") -> Fill((child -> phi() - gchild ->phi()), (child -> eta() - gchild -> eta()));
              }
            }
          } 
          
          
          if (muoncount == 2){
            const xAOD::TruthParticle* gchild0 = child->child(0);
            const xAOD::TruthParticle* gchild1 = child->child(1);
            hist("muon/t_doubledeta") -> Fill(gchild0 -> eta() - gchild1 -> eta());
            hist("muon/t_doubledphi") -> Fill(gchild0 -> phi() - gchild1 -> phi());
            hist("muon/t_doubledetavdphi") -> Fill((gchild0 -> eta() - gchild1 -> eta()),(gchild0 -> phi() - gchild1 -> phi()));
          }

            
              
          if(muoncount != 0) hist("muon/t_muoncounthistogram") -> Fill(muoncount);
          //ANA_MSG_INFO("muon count " << muoncount);
        } //end of child loop
      } 
    } 

    //plots in muon maximum d0, this time coming from the prompt ones
    //signal is found somewhere up there
    if(truth -> prodVtx() != nullptr){
      if (maxedd0_muon_prt < abs(truthd0(truth))){
        maxedd0_muon_prt = abs(truthd0(truth));
      }
    }
  } //end of truth loop


  //WIP
  if(m_trigger_read){
    for (const xAOD::TrackParticle* FTF_T : *trig_FTFparticles){
      if( abs(FTF_T -> d0()) > maxedd0_FTF){
        maxedd0_FTF = abs(FTF_T -> d0());
      }
    }
    for (const xAOD::TrackParticle* LRT_T : *trig_LRTparticles){
      if( abs(LRT_T -> d0()) > maxedd0_FTF){
        maxedd0_LRT = abs(LRT_T -> d0());
      }
    }
  }
    

  
  
  hist("trigger/h_FTF_d0max_event") -> Fill(maxedd0_FTF);
  hist("trigger/h_LRT_d0max_event") -> Fill(maxedd0_LRT);

  hist("trigger/h_muonsig_d0max_event") -> Fill(maxedd0_muon_sig);
  hist("trigger/h_muonprt_d0max_event") -> Fill(maxedd0_muon_prt);

  hist("trigger/h_LRT_trigd_d0max_event") -> Fill(maxedd0_LRT_trigd );

  return StatusCode::SUCCESS;
}



StatusCode MyxAODAnalysis :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.


  //formats some necessary plots 
  hist ("FTF_h_d0eff") -> Divide(hist ("FTF_h_d0eff_n"), hist ("FTF_h_d0eff_d"));
  hist ("FTF_h_d0eff") ->SetMarkerStyle(3);
  hist ("FTF_h_d0eff") -> SetOption("P0"); 

  hist ("trig_h_d0eff") -> Divide(hist ("trig_h_d0eff_n"), hist ("trig_h_d0eff_d"));
  hist ("trig_h_d0eff") ->SetMarkerStyle(3);
  hist ("trig_h_d0eff") -> SetOption("P0"); 

  hist ("trig_h_etaeff") -> Divide(hist ("trig_h_etaeff_n"), hist ("trig_h_etaeff_d"));
  hist ("trig_h_etaeff") ->SetMarkerStyle(3);
  hist ("trig_h_etaeff") -> SetOption("P0"); 

  hist ("trig_h_pTeff") -> Divide(hist ("trig_h_pTeff_n"), hist ("trig_h_pTeff_d"));
  hist ("trig_h_pTeff") ->SetMarkerStyle(3);
  hist ("trig_h_pTeff") -> SetOption("P0");

  hist ("trig_h_pTefflog") -> Divide(hist ("trig_h_pTefflog_n"), hist ("trig_h_pTefflog_d"));
  hist ("trig_h_pTefflog") ->SetMarkerStyle(3);
  hist ("trig_h_pTefflog") -> SetOption("P0");

  hist ("FTF_h_d0fakes") -> Divide(hist ("FTF_h_d0fakes_n"), hist ("FTF_h_d0fakes_d"));
  hist ("FTF_h_d0fakes") ->SetMarkerStyle(3);

  hist ("FTF_h_etafakes") -> Divide(hist ("FTF_h_etafakes_n"), hist ("FTF_h_etafakes_d"));
  hist ("FTF_h_etafakes") ->SetMarkerStyle(3);

  hist ("LRT_h_d0fakes") -> Divide(hist ("LRT_h_d0fakes_n"), hist ("LRT_h_d0fakes_d"));
  hist ("LRT_h_d0fakes") ->SetMarkerStyle(3);

  hist ("LRT_h_etafakes") -> Divide(hist ("LRT_h_etafakes_n"), hist ("LRT_h_etafakes_d"));
  hist ("LRT_h_etafakes") ->SetMarkerStyle(3);

  hist("compare/h_d0truthvtrack_All") -> Add(hist("compare/h_d0truthvtrack_FTF"));
  hist("compare/h_d0truthvtrack_All") -> Add(hist("compare/h_d0truthvtrack_LRT"));
  hist("offline_h_dphivTDLength") -> SetOption("box");
  hist("FTF_h_dphivTDLength") -> SetOption("box");
  hist("LRT_h_dphivTDLength") -> SetOption("box");

  hist("LRT_h_drvTDLength") -> SetOption("box");
  hist("LRT_h_d0vTDLength") -> SetOption("box");
  hist("LRT_hNoCut_drvTDLength") -> SetOption("box");
  
  hist("FTF_h_NPixvd0") -> SetOption("box");
  hist("FTF_h_NPixvTDLength") -> SetOption("box");
  hist("FTF_h_dzvz") -> SetOption("box");

  hist("LRT_h_NPixvd0") -> SetOption("box");
  hist("LRT_h_NPixvTDLength") -> SetOption("box");
  hist("LRT_h_NSctvd0") -> SetOption("box");
  hist("LRT_h_NSctvTDLength") -> SetOption("box");
  hist("LRT_h_NClustervd0") -> SetOption("box");
  hist("LRT_h_NClustervTDLength") -> SetOption("box");

  hist("LRT_h_dzvz") -> SetOption("box");

  hist("muon/t_multidetavdphi") -> SetOption("box");
  hist("muon/t_doubledetavdphi") -> SetOption("box");

  ANA_MSG_INFO("Algorithm has finished running");
  return StatusCode::SUCCESS;
}