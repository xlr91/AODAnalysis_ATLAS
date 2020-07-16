#include <AsgTools/MessageCheck.h>
#include <MyAnalysis/MyxAODAnalysis.h>
#include <xAODEventInfo/EventInfo.h>
#include <cmath>
#include "TEfficiency.h"
#include "TCanvas.h"

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
  declareProperty( "etacut", m_etacut = 0.1, "Maximum value for the eta cut");
  declareProperty( "phicut", m_phicut = 0.1, "Maximum value for the phi cut");
  declareProperty( "drcut" , m_drcut  = 0.002, "Maximum value for the dr  cut");

  declareProperty( "TriggerRead", m_trigger_read = true, "If it reads the trigger containers");
  declareProperty( "OfflineRead", m_offline_read = true, "If it reads the offline containers");
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
  /*
  Float_t num = (((truth_p -> prodVtx()) -> x()) * (truth_p -> px())) + (((truth_p -> prodVtx()) -> y()) * (truth_p -> py()));
  Float_t dem = pow(((truth_p -> prodVtx()) -> x()), 2.0) + pow(((truth_p -> prodVtx()) -> y()), 2.0);
  return num/(sqrt(dem));
  */

  /*
  const xAOD::TruthVertex* pv = truth_p -> prodVtx();
  Float_t num = (truth_p -> px() * pv -> y()) - (truth_p -> py() * pv -> x());
  Float_t dem = sqrt(  pow( truth_p -> px() , 2.0) + pow(truth_p -> py() , 2.0) );
  //Float_t num = (( pV-> y()) * (truth_p -> px())) - (((truth_p -> prodVtx()) -> x()) * (truth_p -> py()));
  ///Float_t dem = pow(((truth_p -> prodVtx()) -> x()), 2.0) + pow(((truth_p -> prodVtx()) -> y()), 2.0);
  //return num/(sqrt(dem));
  return num/dem;
  
  */

  ///*
  //https://arxiv.org/pdf/1405.6569.pdf page 35
  Float_t ans = ((truth_p -> prodVtx()) -> x() - (truth_v -> x())) * sin(truth_p -> phi()) - ((truth_p -> prodVtx()) -> y()- (truth_v -> y())) * cos(truth_p -> phi());
  return -ans;
  //*/
  
}


bool MyxAODAnalysis::cut1(const xAOD::TruthParticle* truth_p, const xAOD::TrackParticle* track_p, Float_t eta_c, Float_t eta_p){
  bool ans = true;
  if( abs((truth_p->eta() - track_p->eta())) > eta_c  ||  
      abs((truth_p->phi() - track_p->phi())) > eta_p  ){
    
    ans = false;
  }
  return ans;
}

bool MyxAODAnalysis::cut2(Float_t mndr, Float_t cut_c){
  bool ans = true;
  if(mndr > cut_c) ans = false;
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


  //Truth Histograms
  ANA_CHECK(book(TH1F("truth_h_dr", "dR values for truth muons", 100, 0, 0.005)));
  ANA_CHECK(book(TH1F("truth_h_d0", "d0 values for truth muons", 300, -30, 30)));
  ANA_CHECK(book(TH1F("truth_h_eta", "eta values for truth muons", 300, -5, 5)));
  ANA_CHECK(book(TH1F("truth_h_pT", "pT values for truth muons", 300, 0, 2)));
  ANA_CHECK(book(TH1F("truth_h_dphi", "phi values for truth tracks", 300, -0.001, 0.001)));
  ANA_CHECK(book(TH1F("truth_h_deta", "eta values for truth tracks", 300, -5, 5)));
  ANA_CHECK(book(TH2F("truth_h_dphivTDLength", "dPhi vs Transv. DecLengths for truth", 300, -5, 5, 300, -5, 5)));

  //Offline Histograms
  ANA_CHECK(book(TH1F("offline_h_dr", "dR values for offline tracks", 100, 0, 0.005)));
  ANA_CHECK(book(TH1F("offline_h_d0", "d0 values for offline tracks", 300, -30, 30)));
  ANA_CHECK(book(TH1F("offline_h_eta", "eta values for offline tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("offline_h_pT", "pT values for offline tracks", 300, 0, 2)));
  ANA_CHECK(book(TH1F("offline_h_dphi", "phi values for offline tracks", 300, -0.001, 0.001)));
  ANA_CHECK(book(TH1F("offline_h_deta", "eta values for offline tracks", 300, -5, 5)));
  ANA_CHECK(book(TH2F("offline_h_dphivTDLength", "dPhi vs Transv. DecLengths for offline", 300, -5, 5, 300, -5, 5)));

  //FTF Histograms
  ANA_CHECK(book(TH1F("FTF_h_dr", "dR values for FTF tracks", 100, 0, 0.005)));
  ANA_CHECK(book(TH1F("FTF_h_d0", "d0 values for FTF tracks", 300, -50, 50))); //-30/30
  ANA_CHECK(book(TH1F("FTF_h_eta", "eta values for FTF tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("FTF_h_pT", "pT values for FTF tracks", 300, 0, 2)));
  ANA_CHECK(book(TH1F("FTF_h_dphi", "phi values for FTF tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("FTF_h_deta", "eta values for FTF tracks", 300, -5, 5)));
  ANA_CHECK(book(TH2F("FTF_h_dphivTDLength", "dPhi vs Transv. DecLengths for FTF", 300, -5, 5, 300, -5, 5)));

  //LRT Histograms
  ANA_CHECK(book(TH1F("LRT_h_dr", "dR values for LRT tracks", 100, 0, 0.005)));
  ANA_CHECK(book(TH1F("LRT_h_d0", "d0 values for LRT tracks", 300, -50, 50)));
  ANA_CHECK(book(TH1F("LRT_h_eta", "eta values for LRT tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("LRT_h_pT", "pT values for LRT tracks", 300, 0, 2)));
  ANA_CHECK(book(TH1F("LRT_h_dphi", "phi values for LRT tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("LRT_h_deta", "eta values for LRT tracks", 300, -5, 5)));
  ANA_CHECK(book(TH2F("LRT_h_dphivTDLength", "dPhi vs Transv. DecLengths for LRT", 300, -5, 5, 300, -5, 5)));

  //Comparison Histograms
  ANA_CHECK(book(TH1F("compare/h_d0diff_offline", "delta_d0 (Offline-truth)", 300, -2, 2)));
  ANA_CHECK(book(TH1F("compare/h_d0diff_FTF", "delta_d0 (FTF-truth)", 300, -2, 2)));
  ANA_CHECK(book(TH1F("compare/h_d0diff_LRT", "delta_d0 (LRT-truth)", 300, -2, 2)));

  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_offline", "truth_d0_vs_offline_d0", 300, -10, 10, 300, -50, 50)));
  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_FTF", "truth_d0_vs_FTF_d0", 300, -50, 50, 300, -50, 50))); // -10/ 10
  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_LRT", "truth_d0_vs_LRT_d0", 300, -50, 50, 300, -50, 50))); 
  ANA_CHECK(book(TH2F("compare/h_d0truthvtrack_All", "truth_d0_vs_All_d0", 300, -10, 10, 300, -50, 50))); 

  
  ANA_CHECK(book(TH1F("h_d0eff_n", "Efficiency_function_of_d0_n", 300, -50, 50)));
  ANA_CHECK(book(TH1F("h_d0eff_d", "Efficiency_function_of_d0_d", 300, -50, 50)));
  ANA_CHECK(book(TH1F("h_d0eff", "Efficiency_function_of_d0", 300, -50, 50)));
  

  //ANA_CHECK(regEfficiency("testthign", TEfficiency("testeff","Efficiency (Unmanaged)",300, -30, 30)));
  //ANA_CHECK(regEfficiency("test"));

  //  ANA_CHECK(pEff = new TEfficiency("Efficiency","Efficiency (Unmanaged)",300, -30, 30));

  /*
  //ANA_CHECK(book(TH1F("h_d0eff", "Efficiency_function_of_d0", 300, -30, 30)));
  //ANA_CHECK(book(TH1F("h_etaeff", "Efficiency_function_of_eta", 300, 0, 5)));
  ///ANA_CHECK(book(TH1F("h_pTeff", "Efficiency_function_of_pT", 300, 0, 10)));
  */
  

  //check number of parents the muons have (just for head purposes)
  ///ANA_CHECK(book(TH1F("h_muon_parent", "MuonParents", 100, 0, 10)));






  
  

  ANA_MSG_INFO("Offline: " << m_offline_read << " Trigger Read: " <<  m_trigger_read);

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
    if (truth->absPdgId() == 1000006) {     
      if (truth->nChildren() > 1) {
        for (int ichild=0; ichild< truth->nChildren() ; ichild++) {
          
          const xAOD::TruthParticle* child=truth->child(ichild);
          for (int igchild=0; igchild< child->nChildren() ; igchild++) {
            
            const xAOD::TruthParticle* gchild=child->child(igchild);
            if(gchild == nullptr){
              ANA_MSG_WARNING("Nullptr alert in gchild");
              continue;
            }

            if (gchild->absPdgId() == 13){ // at this point everything below are muons from RHadrons from Stops

              //Get the decay lengths of stop (expected to be 0) and RHadron (expected to be about 20 mm)
              const xAOD::TruthVertex* tproVtx = truth->prodVtx(); 
              const xAOD::TruthVertex* tdecVtx = truth->decayVtx(); 
              const xAOD::TruthVertex* cproVtx = child->prodVtx(); 
              const xAOD::TruthVertex* cdecVtx = child->decayVtx(); 

              if (tproVtx == nullptr || tdecVtx == nullptr || cproVtx == nullptr || cdecVtx== nullptr){
                 ANA_MSG_WARNING("Nullptr alert in gchild vector"); 
                 continue;
              }


              truthd0val = truthd0(gchild, cproVtx);
              RhTD_Length = sqrt(pow((cproVtx->x() - cdecVtx->x()), 2.0) + pow((cproVtx->y() - cdecVtx->y()), 2.0)); 
              hist ("h_truthDecayLength")->Fill (decaylength(tproVtx, tdecVtx));
              hist ("h_childDecayLength")->Fill (decaylength(cproVtx, cdecVtx));


              


              ////Offline Tracks
              if(m_offline_read){
                
                mindr = 2000;
                matched_offline = nullptr;
                //Track Loop
                for (const xAOD::TrackParticle* offline : *offline_particles){
                  track_test.push_back(offline); //sanity check
                  //finds matched track
                  if(mindr > calcdr(gchild, offline)){
                    mindr = calcdr(gchild, offline);
                    matched_offline = offline;
                  }
                } //end of track loop
                
                
                //Result of matching truth muon tracks with the reco tracks
                hist ("offline_h_dr")->Fill(0);
                hist ("truth_h_d0")->Fill (truthd0val);
                //ANA_MSG_INFO("PV x: " << cproVtx -> x() << " y: " << cproVtx -> y()  << " z: " << cproVtx -> z() );
                hist ("truth_h_eta")->Fill (gchild -> eta());                
                hist ("truth_h_pT")->Fill (gchild -> pt() / 1000000);

                //Cuts
                passedflag_ofl = cut1(gchild, matched_offline, m_etacut, m_phicut);
                //passedflag_ofl = cut2(mindr, m_drcut);

                if (passedflag_ofl){
                  hist ("offline_h_dr")->Fill (mindr);
                  hist ("offline_h_d0")->Fill (matched_offline -> d0());
                  hist ("offline_h_eta")->Fill (matched_offline -> eta());
                  hist ("offline_h_pT")->Fill (matched_offline -> pt() / 1000000);
                  hist ("compare/h_d0diff_offline")->Fill ( matched_offline -> d0() - truthd0val);
                  hist ("compare/h_d0truthvtrack_offline")->Fill (matched_offline -> d0(), truthd0val);

                  hist ("offline_h_dphi") -> Fill ( gchild->phi() - matched_offline->phi() );
                  hist ("offline_h_deta") -> Fill ( gchild->eta() - matched_offline->eta());
                  hist ("offline_h_dphivTDLength")->Fill (RhTD_Length, gchild->phi() - matched_offline->phi());
                  
                }

                //pEff->Fill(true,matched_offline -> d0());
                
                ANA_MSG_INFO("Track Pointer: " << matched_offline << " Truth Pointer: " << gchild << " mindr: " << mindr);
              }

              //The great trigger algorithm
              ///////////FTF
              if(m_trigger_read){
                mindr = 2000;
                matched_FTF = nullptr;
                for (const xAOD::TrackParticle* FTF_T : *trig_FTFparticles){
                  //finds matched track
                  if(mindr > calcdr(gchild, FTF_T)){
                    mindr = calcdr(gchild, FTF_T);
                    matched_FTF = FTF_T;
                  }
                } //end of FTF loop


                //cuts
                passedflag_ftf = cut1(gchild, matched_FTF, m_etacut, m_phicut);
                //passedflag_ftf = cut2(mindr, m_drcut);
                if (passedflag_ftf){
                  hist ("FTF_h_dr")->Fill (mindr);
                  hist ("FTF_h_d0")->Fill (matched_FTF -> d0());
                  hist ("FTF_h_eta")->Fill (matched_FTF -> eta());
                  hist ("FTF_h_pT")->Fill (matched_FTF -> pt() / 1000000);
                  hist ("compare/h_d0diff_FTF")->Fill ( matched_FTF -> d0() - truthd0val);
                  hist ("compare/h_d0truthvtrack_FTF")->Fill (matched_FTF -> d0(), truthd0val);

                  hist ("FTF_h_dphi") -> Fill ( gchild->phi() - matched_FTF->phi() );
                  hist ("FTF_h_deta") -> Fill ( gchild->eta() - matched_FTF->eta());
                  hist ("FTF_h_dphivTDLength")->Fill (RhTD_Length, gchild->phi() - matched_FTF->phi());
                }
                
                hist ("h_d0eff_d") -> Fill(matched_FTF -> d0());

                ///////////LRT
                mindr = 2000;
                matched_LRT = nullptr;
                for (const xAOD::TrackParticle* LRT_T : *trig_LRTparticles){
                  //track_test.push_back(offline); //sanity check

                  //Monitoring Histograms
                  //hist ("h_phiInOffline")->Fill(offline -> phi());
                  //hist ("h_etaInOffline")->Fill(offline -> eta());
                  
                  //finds matched track
                  if(mindr > calcdr(gchild, LRT_T)){
                    mindr = calcdr(gchild, LRT_T);
                    matched_LRT = LRT_T;
                  }
                } //end of LRT loop
                
                //cuts
                passedflag_lrt = cut1(gchild, matched_LRT, m_etacut, m_phicut);
                //passedflag_lrt = cut2(mindr, m_drcut);
                if (passedflag_lrt){
                  hist ("LRT_h_dr")->Fill (mindr);
                  hist ("LRT_h_d0")->Fill (matched_LRT -> d0());
                  hist ("LRT_h_eta")->Fill (matched_LRT -> eta());
                  hist ("LRT_h_pT")->Fill (matched_LRT -> pt() / 1000000);
                  hist ("compare/h_d0diff_LRT")->Fill ( matched_LRT -> d0() - truthd0val);
                  hist ("compare/h_d0truthvtrack_LRT")->Fill (matched_LRT -> d0(), truthd0val);

                  hist ("LRT_h_dphi") -> Fill ( gchild->phi() - matched_LRT->phi() );
                  hist ("LRT_h_deta") -> Fill ( gchild->eta() - matched_LRT->eta());
                  hist ("LRT_h_dphivTDLength")->Fill (RhTD_Length, gchild->phi() - matched_LRT->phi());
                }

                ///efficiency plot uses the truth values, but 'yes/no' on the booleans
                if(passedflag_ftf || passedflag_lrt){
                  hist ("h_d0eff_n") -> Fill(truthd0val);
                }

                hist ("h_d0eff_d") -> Fill(truthd0val);
              }


            }
          }// end of gchild loop
        } //end of child loop
      } 
    } 
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

  //hist ("h_d0eff_n") -> Add(hist("LRT_h_d0"));
  //hist ("h_d0eff_n") -> Add(hist("FTF_h_d0"));


  hist ("h_d0eff") -> Divide(hist ("h_d0eff_n"), hist ("h_d0eff_d"));

  hist ("h_d0eff") ->SetMarkerStyle(3);
  //hist ("h_d0eff") -> SetOption("E1*");

  //hist ("h_d0eff") ->Divide(hist ("h_d0eff_n"), hist ("truth_h_d0"));
  hist("compare/h_d0truthvtrack_All") -> Add(hist("compare/h_d0truthvtrack_FTF"));
  hist("compare/h_d0truthvtrack_All") -> Add(hist("compare/h_d0truthvtrack_LRT"));
  //hist ("h_etaeff") ->Divide(hist ("offline_h_eta"), hist ("truth_h_eta"));
  //hist ("h_pTeff") -> Fill Divide(hist ("offline_h_pT"), hist ("truth_h_pT"));

  //TCanvas *c2 = new TCanvas("c2"," Efficiency ",50,50,1680,1050);
  //c1 = new TCanvas("c","c",1680,1050);
  //pEff ->Draw("AP") >> hist ("h_pTeff");
  //c1 -> Print("plots/pEffThing.png");

  ANA_MSG_INFO("ENDS YEET with size of vector test " << vector_test.size());
  return StatusCode::SUCCESS;
}