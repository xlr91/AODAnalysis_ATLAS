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

Float_t MyxAODAnalysis::truthd0(const xAOD::TruthParticle* truth_p){
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
  Float_t ans = ((truth_p -> prodVtx()) -> x()) * sin(truth_p -> phi()) - ((truth_p -> prodVtx()) -> y()) * cos(truth_p -> phi());
  return -ans;
  //*/
  
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

  //Triggering Histograms
  ANA_CHECK(book(TH1F("h_drvalues_offline", "dR_values", 100, 0, 0.005)));

  ANA_CHECK(book(TH1F("h_d0truth_offline", "d0_values_for_truth_muons", 300, -30, 30)));
  ANA_CHECK(book(TH1F("h_d0track_offline", "d0_values_for_matched_reco_tracks", 300, -30, 30)));
  ANA_CHECK(book(TH1F("h_d0eff", "Efficiency_function_of_d0", 300, -30, 30)));
  ANA_CHECK(book(TH1F("h_d0diff_offline", "delta_d0 (track-truth)", 300, -2, 2)));

  ANA_CHECK(book(TH1F("h_etatruth_offline", "eta_values_for_truth_muons", 300, -5, 5)));
  ANA_CHECK(book(TH1F("h_etatrack_offline", "eta_values_for_matched_reco_tracks", 300, -5, 5)));
  ANA_CHECK(book(TH1F("h_etaeff", "Efficiency_function_of_eta", 300, 0, 5)));

  ANA_CHECK(book(TH1F("h_pTtruth_offline", "pT_values_for_truth_muons", 300, 0, 2)));
  ANA_CHECK(book(TH1F("h_pTtrack_offline", "pT_values_for_matched_reco_tracks", 300, 0, 2)));
  ANA_CHECK(book(TH1F("h_pTeff", "Efficiency_function_of_pT", 300, 0, 10)));
  
  ANA_CHECK(book(TH2F("h_d0truthvtrack_offline", "truth_d0_vs_track_d0", 300, -10, 10, 300, -10, 10)));

  //check number of parents the muons have (just for head purposes)
  ANA_CHECK(book(TH1F("h_muon_parent", "MuonParents", 100, 0, 10)));

  pEff = new TEfficiency("Efficiency","Efficiency (Unmanaged)",300, -30, 30);
  

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
/*
  const xAOD::EventInfo *eventInfo = nullptr;
  ANA_CHECK (evtStore()->retrieve (eventInfo, "EventInfo"));

  // print out run and event number from retrieved object
  ANA_MSG_INFO ("in execute, runNumber = " << eventInfo->runNumber() << ", eventNumber = " << eventInfo->eventNumber());
*/

// loop over the jets in the container
  /*
  const xAOD::JetContainer* jets = nullptr;  
ANA_CHECK (evtStore()->retrieve (jets, "AntiKt4EMTopoJets"));
  for (const xAOD::Jet* jet : *jets) {
    ANA_MSG_INFO ("Jet pT " << jet->pt() );
    hist ("h_jetPt")->Fill (jet->pt() * 0.001); // GeV
  } // end for loop over jets
  */


  //// Trigger
  /*
  const xAOD::TrackParticleContainer* tracks_mu = nullptr;
  ANA_CHECK (evtStore()->retrieve (tracks_mu, "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Muon_FTF"));
    ANA_MSG_INFO ("Number of muon tracks " << tracks_mu->size() );
  for (const xAOD::TrackParticle* track : *tracks_mu) {
    ANA_MSG_INFO ("Track pT " << track->pt() );
    //    hist ("h_jetPt")->Fill (jet->pt() * 0.001); // GeV
  } // end for loop over jets

  
  const xAOD::TrackParticleContainer* tracks_lrt = nullptr;
  ANA_CHECK (evtStore()->retrieve (tracks_lrt, "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Beamspotlrt_FTF"));
    ANA_MSG_INFO ("Number of LRT tracks " << tracks_lrt->size() );
  for (const xAOD::TrackParticle* track : *tracks_lrt) {
    ANA_MSG_INFO ("Track pT " << track->pt() );
  } 

  const xAOD::TrackParticleContainer* tracks_bs = nullptr;
  ANA_CHECK (evtStore()->retrieve (tracks_bs, "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_Beamspot_FTF"));
    ANA_MSG_INFO ("Number of Beamspot tracks " << tracks_bs->size() );
  for (const xAOD::TrackParticle* track : *tracks_bs) {
    ANA_MSG_INFO ("Track pT " << track->pt() );
  } 
  */

  

  //Access InDetTrackParticles
  const xAOD::TruthParticleContainer* truthparticles;
  const xAOD::TrackParticleContainer* offlineparticles;
  const xAOD::TrackParticleContainer* trig_FTFparticles;
  const xAOD::TrackParticleContainer* trig_LRTparticles;
  const xAOD::TruthParticle* matched_truth;
  const xAOD::TrackParticle* matched_track;
  Float_t mindr;

  
  ANA_CHECK (evtStore() -> retrieve (truthparticles, "TruthParticles"));
  ANA_MSG_INFO("Found Truth, size is " << truthparticles->size());

  if (m_offline_read) {
    ANA_CHECK (evtStore() -> retrieve (trig_FTFparticles, "HLT_xAOD__TrackParticleContainer_InDetTrigTrackingxAODCnv_FullScan_FTF"));
    ANA_MSG_INFO("Found FTF, size is " << trig_FTFparticles->size());
  }
  if (m_trigger_read){
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
            if (gchild->absPdgId() == 13){ // at this point everything below are muons from RHadrons from Stops

              hist ("h_muon_parent")->Fill (gchild -> nParents()); //sanity check (unused)

              //Get the decay lengths of stop (expected to be 0) and RHadron (expected to be about 20 mm)
              const xAOD::TruthVertex* tproVtx = truth->prodVtx(); 
              const xAOD::TruthVertex* tdecVtx = truth->decayVtx(); 
              const xAOD::TruthVertex* cproVtx = child->prodVtx(); 
              const xAOD::TruthVertex* cdecVtx = child->decayVtx(); 
              hist ("h_truthDecayLength")->Fill (decaylength(tproVtx, tdecVtx));
              hist ("h_childDecayLength")->Fill (decaylength(cproVtx, cdecVtx));


              ////Offline Tracks
              if(m_offline_read){
                mindr = 2000;
                matched_track = nullptr;
                matched_truth = gchild;
  
                //Track Loop
                for (const xAOD::TrackParticle* offline : *offlineparticles){
                  track_test.push_back(offline); //sanity check

                  //Monitoring Histograms
                  hist ("h_phiInOffline")->Fill(offline -> phi());
                  hist ("h_etaInOffline")->Fill(offline -> eta());
                  
                  //finds matched track
                  if(mindr > calcdr(gchild, offline)){
                    mindr = calcdr(gchild, offline);
                    matched_track = offline;
                  }
                } //end of track loop
                
                if (mindr == 2000) matched_truth = nullptr;

                //Cuts
                Bool_t passedflag = true;
                if( abs((matched_truth->eta() - matched_track->eta())) > m_etacut  ||  
                    abs((matched_truth->phi() - matched_track->phi())) > m_phicut  ){
                  passedflag = false;
                }
                
                //Result of matching truth muon tracks with the reco tracks
                hist ("h_drvalues_offline")->Fill (mindr);

                hist ("h_d0truth_offline")->Fill (truthd0(matched_truth));
                hist ("h_d0track_offline")->Fill (matched_track -> d0());
                hist ("h_d0diff_offline")->Fill ( matched_track -> d0() - truthd0(matched_truth));
        
                hist ("h_etatruth_offline")->Fill (matched_truth -> eta());
                hist ("h_etatrack_offline")->Fill (matched_track -> eta());
                hist ("h_pTtruth_offline")->Fill (matched_truth -> pt() / 1000000);
                hist ("h_pTtrack_offline")->Fill (matched_track -> pt() / 1000000);

                hist ("h_d0truthvtrack_offline")->Fill (matched_track -> d0(), truthd0(matched_truth));

                pEff->Fill(true,matched_track -> d0());
                
                ANA_MSG_INFO("Track Pointer: " << matched_track << " Truth Pointer: " << gchild << " mindr: " << mindr);
              }


              ////FTF Tracks


              ////LRT Tracks
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



  hist ("h_d0eff") ->Divide(hist ("h_d0track_offline"), hist ("h_d0truth_offline"));
  hist ("h_etaeff") ->Divide(hist ("h_etatrack_offline"), hist ("h_etatruth_offline"));
  //hist ("h_pTeff") -> Fill Divide(hist ("h_pTtrack_offline"), hist ("h_pTtruth_offline"));

  //TCanvas *c2 = new TCanvas("c2"," Efficiency ",50,50,1680,1050);
  //c1 = new TCanvas("c","c",1680,1050);
  //pEff ->Draw("AP") >> hist ("h_pTeff");
  //c1 -> Print("plots/pEffThing.png");

  ANA_MSG_INFO("ENDS YEET with size of vector test " << vector_test.size());
  return StatusCode::SUCCESS;
}