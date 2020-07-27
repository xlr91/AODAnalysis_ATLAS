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
  // declareProperty( "nonSTOP", m_nonSTOP = 0, "Desc?");
  //declareProperty( "vector_test", vector_test, "A test for the existence of vectors ykno");
  //declareProperty( "truth_vector", truth_vector, "Pointer Vector of Truth particles");
  // declareProperty( "TitleforJobOption", codetitle = 0, "Desc?");
  declareProperty( "drcut" , m_drcut  = 0.4, "Maximum value for the dr  cut");
  declareProperty( "etacut", m_etacut = 0.5, "Maximum value for the eta cut");
  declareProperty( "phicut", m_phicut = 0.3, "Maximum value for the phi cut");
  declareProperty( "dzcut" , m_dzcut  = 200, "Maximum value for the dz  cut");

  declareProperty( "TriggerRead", m_trigger_read = true, "If it reads the trigger containers");
  declareProperty( "OfflineRead", m_offline_read = true, "If it reads the offline containers");

  declareProperty( "cutnumber", m_cut = 2, "Choose what kind of cut you want");
}

Double_t MyxAODAnalysis::decaylength(const xAOD::TruthVertex* x1, const xAOD::TruthVertex* x2){

  double_t result= pow((x1->x() - x2->x()), 2.0) + pow((x1->y() - x2->y()), 2.0) + pow((x1->z() - x2->z()), 2.0);
  return sqrt(result);
}

Float_t MyxAODAnalysis::calcdr(const xAOD::TruthParticle* truth_p, const xAOD::TrackParticle* track_p){
  Float_t dr2 = pow((truth_p->eta() - track_p->eta()), 2.0) + pow((truth_p->phi() - track_p->phi()), 2.0);
  return sqrt(dr2);
}

Float_t MyxAODAnalysis::truthd0(const xAOD::TruthParticle* truth_p, const xAOD::TruthVertex* truth_v){

  //https://arxiv.org/pdf/1405.6569.pdf page 35
  Float_t ans = ((truth_p -> prodVtx()) -> x() - (truth_v -> x())) * sin(truth_p -> phi()) - ((truth_p -> prodVtx()) -> y()- (truth_v -> y())) * cos(truth_p -> phi());
  return -ans;
  //*/
  
}

Float_t MyxAODAnalysis::truthd0(const xAOD::TruthParticle* truth_p){
  //https://arxiv.org/pdf/1405.6569.pdf page 35
  Float_t ans = ((truth_p -> prodVtx()) -> x() ) * sin(truth_p -> phi()) - ((truth_p -> prodVtx()) -> y()) * cos(truth_p -> phi());
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

//eventually make it more general the switch cut

//Dr
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
  //ANA_CHECK (book (TH1F ("h_jetPt", "h_jetPt", 100, 0, 500))); // jet pt [GeV]
  //ANA_CHECK(book(TH1F("fileIdentifier", "title in graph", no. of bins, min, max)));

  //Monitoring Histograms
  ANA_CHECK(book(TH1F("h_truthDecayLength", "STop_Decay_Length", 100, 0, 10)));
  ANA_CHECK(book(TH1F("h_childDecayLength", "RHadron_Decay_Length", 100, 0, 150)));
  ANA_CHECK(book(TH1F("h_phiInOffline", "phi_In_Offline", 100, -3.15, -3.15)));
  ANA_CHECK(book(TH1F("h_etaInOffline", "eta_In_Offline", 100, -5, -5)));
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
  ANA_CHECK(book(TH1F("FTF_h_d0", "d0 values for FTF tracks", 300, -100, 100))); //-30/30
  ANA_CHECK(book(TH1F("FTF_h_eta", "eta values for FTF tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("FTF_h_pT", "pT values for FTF tracks", 300, 0, 2)));
  ANA_CHECK(book(TH1F("FTF_h_phi", "phi values for FTF tracks", 300, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("FTF_h_dphi", "dphi values for FTF tracks", 300, -1, 1)));
  ANA_CHECK(book(TH1F("FTF_h_dphi_0-10", "dphi values for FTF tracks 0-10", 300, -1, 1)));
  ANA_CHECK(book(TH1F("FTF_h_dphi_20-30", "dphi values for FTF tracks 20-30", 300, -1, 1)));
  ANA_CHECK(book(TH1F("FTF_h_deta", "deta values for FTF tracks", 300, -1.2, 1.2)));
  ANA_CHECK(book(TH1F("FTF_h_dz", "dz values for FTF tracks", 300, -300, 300)));
  ANA_CHECK(book(TH2F("FTF_h_phivTDLength", "dPhi vs Transv. DecLengths for FTF", 60, 0, 20, 300, -3.2, 3.2)));
  ANA_CHECK(book(TH2F("FTF_h_dphivTDLength", "dPhi vs Transv. DecLengths for FTF", 60, 0, 40, 300, -0.001, 0.001)));

  //LRT Histograms
  ANA_CHECK(book(TH1F("LRT_h_dr", "dR values for LRT tracks", 300, 0, 0.1)));
  ANA_CHECK(book(TH1F("LRT_h_d0", "d0 values for LRT tracks", 300, -100, 100)));
  ANA_CHECK(book(TH1F("LRT_h_eta", "eta values for LRT tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("LRT_h_pT", "pT values for LRT tracks", 300, 0, 2)));
  ANA_CHECK(book(TH1F("LRT_h_phi", "phi values for LRT tracks", 300, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("LRT_h_dphi", "dphi values for LRT tracks", 300, -1, 1)));
  ANA_CHECK(book(TH1F("LRT_h_dphi_0-10", "dphi values for LRT tracks 0-10", 300, -1, 1)));
  ANA_CHECK(book(TH1F("LRT_h_dphi_20-30", "dphi values for LRT tracks 20-30", 300, -1, 1)));
  ANA_CHECK(book(TH1F("LRT_h_deta", "deta values for LRT tracks", 300, -1.2, 1.2)));
  ANA_CHECK(book(TH1F("LRT_h_dz", "dz values for LRT tracks", 300, -300, 300)));
  ANA_CHECK(book(TH2F("LRT_h_phivTDLength", "dPhi vs Transv. DecLengths for LRT", 60, 0, 20, 300, -3.2, 3.2)));
  ANA_CHECK(book(TH2F("LRT_h_dphivTDLength", "dPhi vs Transv. DecLengths for LRT", 60, 0, 40, 300, -0.001, 0.001)));

  //Comparison Histograms
  ANA_CHECK(book(TH1F("compare/h_d0diff_offline", "delta_d0 (Offline-truth)", 100, -2, 2)));
  ANA_CHECK(book(TH1F("compare/h_d0diff_FTF", "delta_d0 (FTF-truth)", 100, -2, 2)));
  ANA_CHECK(book(TH1F("compare/h_d0diff_LRT", "delta_d0 (LRT-truth)", 100, -2, 2)));

  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_offline", "truth_d0_vs_offline_d0", 100, -10, 10, 100, -50, 50)));
  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_FTF", "truth_d0_vs_FTF_d0", 100, -50, 50, 100, -50, 50))); // -10/ 10
  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_LRT", "truth_d0_vs_LRT_d0", 100, -50, 50, 100, -50, 50))); 
  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_All", "truth_d0_vs_All_d0", 100, -10, 10, 100, -50, 50))); 

  
  //Comparison Histograms
  ANA_CHECK(book(TH1F("offl_h_d0eff_n", "Efficiency_function_of_d0_n", 50, -50, 50)));
  ANA_CHECK(book(TH1F("offl_h_d0eff_d", "Efficiency_function_of_d0_d", 50, -50, 50)));
  ANA_CHECK(book(TH1F("offl_h_d0eff", "Efficiency_function_of_d0", 50, -50, 50)));

  ANA_CHECK(book(TH1F("offl_h_etaeff_n", "Efficiency_function_of_eta_n", 50, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("offl_h_etaeff_d", "Efficiency_function_of_eta_d", 50, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("offl_h_etaeff", "Efficiency_function_of_eta", 50, -3.2, 3.2)));

  ANA_CHECK(book(TH1F("offl_h_pTeff_n", "Efficiency_function_of_pT_n", 50, 0, 2000)));
  ANA_CHECK(book(TH1F("offl_h_pTeff_d", "Efficiency_function_of_pT_d", 50, 0, 2000)));
  ANA_CHECK(book(TH1F("offl_h_pTeff", "Efficiency_function_of_pT", 50, 0, 2000)));


  ANA_CHECK(book(TH1F("FTF_h_d0eff_n", "Efficiency_function_of_d0_n_FTF", 50, -40, 40)));
  ANA_CHECK(book(TH1F("FTF_h_d0eff_d", "Efficiency_function_of_d0_d_FTF", 50, -40, 40)));
  ANA_CHECK(book(TH1F("FTF_h_d0eff", "Efficiency_function_of_d0_FTF", 50, -40, 40)));



  ANA_CHECK(book(TH1F("FTF_h_d0fakes_n", "Fakes_function_of_d0_n_FTF", 50, -40, 40)));
  ANA_CHECK(book(TH1F("FTF_h_d0fakes_d", "Fakes_function_of_d0_d_FTF", 50, -40, 40)));
  ANA_CHECK(book(TH1F("FTF_h_d0fakes", "Fakes_function_of_d0_FTF", 50, -40, 40)));

  ANA_CHECK(book(TH1F("FTF_h_etafakes_n", "Fakes_function_of_eta_n_FTF", 50, -3, 3)));
  ANA_CHECK(book(TH1F("FTF_h_etafakes_d", "Fakes_function_of_eta_n_FTF", 50, -3, 3)));
  ANA_CHECK(book(TH1F("FTF_h_etafakes", "Fakes_function_of_eta_FTF", 50, -3, 3)));

  ANA_CHECK(book(TH1F("LRT_h_d0fakes_n", "Fakes_function_of_d0_n_LRT", 50, -40, 40)));
  ANA_CHECK(book(TH1F("LRT_h_d0fakes_d", "Fakes_function_of_d0_d_LRT", 50, -40, 40)));
  ANA_CHECK(book(TH1F("LRT_h_d0fakes", "Fakes_function_of_d0_LRT", 50, -40, 40)));

  ANA_CHECK(book(TH1F("LRT_h_etafakes_n", "Fakes_function_of_eta_n_LRT", 50, -3, 3)));
  ANA_CHECK(book(TH1F("LRT_h_etafakes_d", "Fakes_function_of_eta_d_LRT", 50, -3, 3)));
  ANA_CHECK(book(TH1F("LRT_h_etafakes", "Fakes_function_of_eta_LRT", 50, -3, 3)));

  ANA_CHECK(book(TH1F("LRT_h_d0eff_n", "Efficiency_function_of_d0_n_LRT", 50, -30, 30)));
  ANA_CHECK(book(TH1F("LRT_h_d0eff_d", "Efficiency_function_of_d0_d_LRT", 50, -30, 30)));
  ANA_CHECK(book(TH1F("LRT_h_d0eff", "Efficiency_function_of_d0_LRT", 50, -30, 30)));


  ANA_CHECK(book(TH1F("trig_h_d0eff_n", "Efficiency_function_of_d0_n", 100, -50, 50)));
  ANA_CHECK(book(TH1F("trig_h_d0eff_d", "Efficiency_function_of_d0_d", 100, -50, 50)));
  ANA_CHECK(book(TH1F("trig_h_d0eff", "Efficiency_function_of_d0", 100, -50, 50)));

  ANA_CHECK(book(TH1F("trig_h_etaeff_n", "Efficiency_function_of_eta_n", 50, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("trig_h_etaeff_d", "Efficiency_function_of_eta_d", 50, -3.2, 3.2)));
  ANA_CHECK(book(TH1F("trig_h_etaeff", "Efficiency_function_of_eta", 50, -3.2, 3.2)));

  ANA_CHECK(book(TH1F("trig_h_pTeff_n", "Efficiency_function_of_pT_n", 50, 0, 150)));
  ANA_CHECK(book(TH1F("trig_h_pTeff_d", "Efficiency_function_of_pT_d", 50, 0, 150)));
  ANA_CHECK(book(TH1F("trig_h_pTeff", "Efficiency_function_of_pT", 50, 0, 150)));


    //logarithmic pT values
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

  ANA_CHECK(book(TH1F("trig_h_TDLeff_n", "Efficiency_function_of_TDL_n", 100, 0, 1000)));
  ANA_CHECK(book(TH1F("trig_h_TDLeff_d", "Efficiency_function_of_TDL_d", 100, 0, 1000)));
  ANA_CHECK(book(TH1F("trig_h_TDLeff", "Efficiency_function_of_TDL", 100, 0, 1000)));

  ANA_CHECK(book(TH1F("tgof_h_d0eff_n", "Efficiency_function_of_d0_n", 50, -40, 40)));
  ANA_CHECK(book(TH1F("tgof_h_d0eff_d", "Efficiency_function_of_d0_d", 50, -40, 40)));
  ANA_CHECK(book(TH1F("tgof_h_d0eff", "Efficiency_function_of_d0", 50, -40, 40)));

  ANA_CHECK(book(TH1F("trig_h_z0eff_n", "Efficiency_function_of_z0_n", 100, -500, 500)));
  ANA_CHECK(book(TH1F("trig_h_z0eff_d", "Efficiency_function_of_z0_d", 100, -500, 500)));
  ANA_CHECK(book(TH1F("trig_h_z0eff", "Efficiency_function_of_z0", 100, -500, 500)));



  

  

  //ANA_CHECK(regEfficiency("testthign", TEfficiency("testeff","Efficiency (Unmanaged)",300, -30, 30)));
  //ANA_CHECK(regEfficiency("test"));

  //  ANA_CHECK(pEff = new TEfficiency("Efficiency","Efficiency (Unmanaged)",300, -30, 30));

  /*
  //ANA_CHECK(book(TH1F("trig_h_d0eff", "Efficiency_function_of_d0", 300, -30, 30)));
  //ANA_CHECK(book(TH1F("h_etaeff", "Efficiency_function_of_eta", 300, 0, 5)));
  ///ANA_CHECK(book(TH1F("h_pTeff", "Efficiency_function_of_pT", 300, 0, 10)));
  */
  

  //check number of parents the muons have (just for head purposes)
  ///ANA_CHECK(book(TH1F("h_muon_parent", "MuonParents", 100, 0, 10)));

  

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
  std::vector<const xAOD::TrackParticle*> track_test;
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  ANA_MSG_INFO ("in execute");
  // retrieve the eventInfo object from the event store

  //Access InDetTrackParticles
  const xAOD::TruthParticleContainer* truthparticles;
  const xAOD::TrackParticleContainer* offline_particles;
  const xAOD::TrackParticleContainer* trig_FTFparticles;
  const xAOD::TrackParticleContainer* trig_LRTparticles;
  
  //const xAOD::TruthParticle* matched_truth;
  const xAOD::TrackParticle* matched_offline;
  const xAOD::TrackParticle* matched_FTF;
  const xAOD::TrackParticle* matched_LRT;
  Float_t mindr;
  Float_t truthd0val;
  Bool_t passedflag_ofl = true;
  Bool_t passedflag_ftf = true;
  Bool_t passedflag_lrt = true;
  Float_t RhTD_Length;
  
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

  //truth loop
  //finds only those that goes STop -> RHadron -> muon
  for (const xAOD::TruthParticle* truth : *truthparticles) {
    ANA_MSG_INFO("PDGID" << truth->absPdgId());
  
    
    if (truth->absPdgId() == 13) {     

              //Get the decay lengths of stop (expected to be 0) and RHadron (expected to be about 20 mm)
              const xAOD::TruthVertex* cdecVtx = truth->prodVtx(); 
              truthd0val = truthd0(truth);
    
              hist ("truth_h_d0")->Fill (truthd0val);
              //ANA_MSG_INFO("PV x: " << cproVtx -> x() << " y: " << cproVtx -> y()  << " z: " << cproVtx -> z() );
              hist ("truth_h_eta")->Fill (truth -> eta());                
              hist ("truth_h_pT")->Fill (truth -> pt() / 1000);


              ////Offline Tracks
              if(m_offline_read){
                mindr = 2000;
                matched_offline = nullptr;
                //Track Loop
                for (const xAOD::TrackParticle* offline : *offline_particles){
                  //finds matched track
                  if(mindr > calcdr(truth, offline)){
                    mindr = calcdr(truth, offline);
                    matched_offline = offline;
                  }
                } //end of track loop
                

                switch(m_cut) {
                  case 0:
                    passedflag_ofl = true;
                    break;
                  case 1:
                    passedflag_ofl = cut1(truth, matched_offline, m_etacut, m_phicut);
                    break;
                  case 2:
                    passedflag_ofl = cut2(mindr, m_drcut);
                    break;
                  case 3:
                    passedflag_ofl = cut2(abs(truth->eta() - matched_offline->eta()), m_etacut);
                    break;
                  case 4:
                    passedflag_ofl = cut2(abs(truth->phi() - matched_offline->phi()), m_phicut);
                    break;
                  case 5:
                    passedflag_ofl = cut2(abs(cdecVtx->z() - matched_offline->z0()), m_dzcut);
                    break;
                  default:
                    passedflag_ofl = true;
                }

                ANA_MSG_INFO("Cut Chosen is " << m_cut);
                //ANA_MSG_INFO("Passed ofl cut: " << passedflag_ofl);

                if (passedflag_ofl){
                  hist("h_offpass") -> Fill(1);
                  hist ("offline_h_dr")->Fill (mindr);
                  hist ("offline_h_d0")->Fill (matched_offline -> d0());
                  hist ("offline_h_eta")->Fill (matched_offline -> eta());
                  hist ("offline_h_phi")->Fill (matched_offline -> phi());
                  hist ("offline_h_pT")->Fill (matched_offline -> pt() / 1000);
                  hist ("compare/h_d0diff_offline")->Fill ( matched_offline -> d0() - truthd0val);
                  hist ("compare/h_d0truthvtrack_offline")->Fill (matched_offline -> d0(), truthd0val);

                  hist ("offline_h_dphi") -> Fill ( truth->phi() - matched_offline->phi() );
                  hist ("offline_h_deta") -> Fill ( truth->eta() - matched_offline->eta());
                  hist ("offline_h_dz") -> Fill ( cdecVtx->z() - matched_offline->z0());

                  hist ("offl_h_d0eff_n") -> Fill(truthd0val);
                  hist ("offl_h_etaeff_n") -> Fill(truth->eta());
                  hist ("offl_h_pTeff_n") -> Fill(truth->pt() / 1000);
                } else{hist("h_offpass") -> Fill(0); }

                hist ("offl_h_d0eff_d") -> Fill(truthd0val);
                hist ("offl_h_etaeff_d") -> Fill(truth->eta());
                hist ("offl_h_pTeff_d") -> Fill(truth->pt() / 1000);
                
              }

              //The great trigger algorithm
              ///////////FTF
              if(m_trigger_read){
                mindr = 2000;
                matched_FTF = nullptr;
                for (const xAOD::TrackParticle* FTF_T : *trig_FTFparticles){
                  //finds matched track
                  if(mindr > calcdr(truth, FTF_T)){
                    mindr = calcdr(truth, FTF_T);
                    matched_FTF = FTF_T;
                  }
                } //end of FTF loop

                //chooses what kind of cuts is used
                switch(m_cut) {
                  case 0:
                    passedflag_ftf = true;
                    break;
                  case 1:
                    passedflag_ftf = cut1(truth, matched_FTF, m_etacut, m_phicut);
                    break;
                  case 2:
                    passedflag_ftf = cut2(mindr, m_drcut);
                    break;
                  case 3:
                    passedflag_ftf = cut2(abs(truth->eta() - matched_FTF->eta()), m_etacut);
                    break;
                  case 4:
                    passedflag_ftf = cut2(abs(truth->phi() - matched_FTF->phi()), m_phicut);
                    break;
                  case 5:
                    passedflag_ftf = cut2(abs(cdecVtx->z() - matched_FTF->z0()), m_dzcut);
                    break;
                  default:
                    passedflag_ftf = true;
                }


                if (matched_FTF == nullptr) {
                  ANA_MSG_WARNING("Null ptr alert in FTF");
                  continue;
                }


                if (passedflag_ftf){
                  hist("h_ftfpass") -> Fill(1);
                  hist ("FTF_h_dr")->Fill (mindr);
                  hist ("FTF_h_d0")->Fill (matched_FTF -> d0());
                  hist ("FTF_h_eta")->Fill (matched_FTF -> eta());
                  hist ("FTF_h_phi")->Fill (matched_FTF -> phi());
                  hist ("FTF_h_pT")->Fill (matched_FTF -> pt() / 1000);
                  hist ("compare/h_d0diff_FTF")->Fill ( matched_FTF -> d0() - truthd0val);
                  hist ("compare/h_d0truthvtrack_FTF")->Fill (matched_FTF -> d0(), truthd0val);

                  hist ("FTF_h_dphi") -> Fill ( truth->phi() - matched_FTF->phi() );
                  hist ("FTF_h_deta") -> Fill ( truth->eta() - matched_FTF->eta());
                  hist ("FTF_h_dz") -> Fill ( cdecVtx->z() - matched_FTF->z0());
                  ANA_MSG_INFO("Truthz: "<< cdecVtx->z() << " Trackz: " << matched_FTF->z0());

                  hist ("FTF_h_d0eff_n") -> Fill(truthd0val);
                } else{hist("h_ftfpass") -> Fill(0);}

                hist ("FTF_h_d0eff_d") -> Fill(truthd0val);                

                ///////////LRT
                mindr = 2000;
                matched_LRT = nullptr;
                for (const xAOD::TrackParticle* LRT_T : *trig_LRTparticles){
                  if(mindr > calcdr(truth, LRT_T)){
                    mindr = calcdr(truth, LRT_T);
                    matched_LRT = LRT_T;
                  }
                } //end of LRT loop
                
          
                
                switch(m_cut) {
                  case 0:
                    passedflag_lrt = true;
                    break;
                  case 1:
                    passedflag_lrt = cut1(truth, matched_LRT, m_etacut, m_phicut);
                    break;
                  case 2:
                    passedflag_lrt = cut2(mindr, m_drcut);
                    break;
                  case 3:
                    passedflag_lrt = cut2(abs(truth->eta() - matched_LRT->eta()), m_etacut);
                    break;
                  case 4:
                    passedflag_lrt = cut2(abs(truth->phi() - matched_LRT->phi()), m_phicut);
                    break;
                  case 5:
                    passedflag_lrt = cut2(abs(cdecVtx->z() - matched_LRT->z0()), m_dzcut);
                    break;
                  default:
                    passedflag_lrt = true;
                }


                if (matched_LRT == nullptr) {
                  ANA_MSG_WARNING("Null ptr alert in lrt");
                  continue;
                }
                
                if (passedflag_lrt){
                  hist("h_lrtpass") -> Fill(1);
                  hist ("LRT_h_d0")->Fill (matched_LRT -> d0());
                  hist ("LRT_h_eta")->Fill (matched_LRT -> eta());
                  hist ("LRT_h_phi")->Fill (matched_LRT -> phi());
                  hist ("LRT_h_pT")->Fill (matched_LRT -> pt() / 1000);
                  hist ("compare/h_d0diff_LRT")->Fill ( matched_LRT -> d0() - truthd0val);
                  hist ("compare/h_d0truthvtrack_LRT")->Fill (matched_LRT -> d0(), truthd0val);

                  hist ("LRT_h_dphi") -> Fill ( truth->phi() - matched_LRT->phi() );
                  hist ("LRT_h_deta") -> Fill ( truth->eta() - matched_LRT->eta());
                  hist ("LRT_h_dz") -> Fill ( cdecVtx->z() - matched_LRT->z0());
                  hist ("LRT_h_d0eff_n") -> Fill(truthd0val);
                } else{hist("h_lrtpass") -> Fill(0);}

                hist ("LRT_h_d0eff_d") -> Fill(truthd0val);
                ///efficiency plot uses the truth values, but 'yes/no' on the booleans
                if(passedflag_ftf || passedflag_lrt){
                  hist ("trig_h_d0eff_n") -> Fill(truthd0val);
                  hist ("trig_h_etaeff_n") -> Fill(truth->eta());
                  hist ("trig_h_pTeff_n") -> Fill(truth->pt() / 1000);
                  hist ("trig_h_z0eff_n") -> Fill(cdecVtx->z());
                }

                hist ("trig_h_d0eff_d") -> Fill(truthd0val);
                hist ("trig_h_etaeff_d") -> Fill(truth->eta());
                hist ("trig_h_pTeff_d") -> Fill(truth->pt() / 1000);
                hist ("trig_h_z0eff_d") -> Fill(cdecVtx->z());
                
              }

              //trig wrt off
              if (m_offline_read && m_trigger_read){
                ///at this point it found the matched offline and matched trigger
                ///so fill in the denominator if the current offline track passed thru the cut 

                if(passedflag_ofl) {
                  hist("tgof_h_d0eff_d") -> Fill(matched_offline -> d0());
                  //if all three matches then fill it in
                  if(passedflag_ftf || passedflag_lrt) hist("tgof_h_d0eff_n") -> Fill(matched_offline -> d0());
                }
              }


    }
    
    //*/

  } //end of truth loop
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


  hist("compare/h_d0truthvtrack_All") -> Add(hist("compare/h_d0truthvtrack_FTF"));
  hist("compare/h_d0truthvtrack_All") -> Add(hist("compare/h_d0truthvtrack_LRT"));
  hist("offline_h_dphivTDLength") -> SetOption("box");
  hist("FTF_h_dphivTDLength") -> SetOption("box");
  hist("LRT_h_dphivTDLength") -> SetOption("box");



  
  //hist ("h_etaeff") ->Divide(hist ("offline_h_eta"), hist ("truth_h_eta"));
  //hist ("h_pTeff") -> Fill Divide(hist ("offline_h_pT"), hist ("truth_h_pT"));


  //  TFile* pFile = new TFile("myfile.root", "CREATE");
  //TCanvas *c2 = new TCanvas("c2"," Efficiency ",50,50,1680,1050);
  //c1 = new TCanvas("c","c",1680,1050);
  //pEff ->Draw("AP") >> hist ("h_pTeff");
  //c1 -> Print("plots/pEffThing.png");




  ANA_MSG_INFO("ENDS YEET with size of vector test " << vector_test.size());
  return StatusCode::SUCCESS;
}