#include <AsgTools/MessageCheck.h>
#include <MyAnalysis/MyxAODAnalysis.h>
#include <xAODEventInfo/EventInfo.h>
#include <cmath>


MyxAODAnalysis :: MyxAODAnalysis (const std::string& name,
                                  ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm (name, pSvcLocator)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  This is also where you
  // declare all properties for your algorithm.  Note that things like
  // resetting statistics variables or booking histograms should
  // rather go into the initialize() function.
  declareProperty( "nonSTOP", m_nonSTOP = 0, "Desc?");
  declareProperty( "vector_test", vector_test, "A test for the existence of vectors ykno");
  //declareProperty( "truth_vector", truth_vector, "Pointer Vector of Truth particles");
  // declareProperty( "TitleforJobOption", codetitle = 0, "Desc?");
                   

}

Double_t MyxAODAnalysis::decaylength(const xAOD::TruthVertex* x1, const xAOD::TruthVertex* x2){
  double_t result= pow((x1->x() - x2->x()), 2.0) + pow((x1->y() - x2->y()), 2.0) + pow((x1->z() - x2->z()), 2.0);
  return sqrt(result);
}

Float_t MyxAODAnalysis::calcd0(const xAOD::TruthParticle* truth_p, const xAOD::TrackParticle* track_p){
  Float_t dr2 = pow((truth_p->eta() - track_p->eta()), 2.0) + pow((truth_p->phi() - track_p->phi()), 2.0);
  return sqrt(dr2);
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
  ANA_CHECK(book(TH1F("h_truthDecayLength", "STop_Decay_Length", 100, 0, 10)));
  ANA_CHECK(book(TH1F("h_childDecayLength", "RHadron_Decay_Length", 100, 0, 150)));
  ANA_CHECK(book(TH1F("h_phiInOffline", "phi_In_Offline", 100, -3.15, -3.15)));
  ANA_CHECK(book(TH1F("h_etaInOffline", "eta_In_Offline", 100, -5, -5)));

  ANA_CHECK(book(TH1F("h_drvalues", "dR_values", 300, 0, 20)));

  //check number of parents the muons have (just for head purposes)
  ANA_CHECK(book(TH1F("h_muon_parent", "MuonParents", 100, 0, 10)));


  
  


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
  vector_test.push_back(1);
  

  //Access InDetTrackParticles
  const xAOD::TrackParticleContainer* offlineparticles;
  const xAOD::TruthParticleContainer* truthparticles;
  ANA_CHECK (evtStore() -> retrieve (truthparticles, "TruthParticles"));
  ANA_CHECK (evtStore() -> retrieve (offlineparticles, "InDetTrackParticles"));
  ANA_MSG_INFO("Found Offline, size is " << offlineparticles->size());
  ANA_MSG_INFO("Found Truth, size is " << truthparticles->size());
  m_nonSTOP = 0;

  //track loop
  for (const xAOD::TrackParticle* offline : *offlineparticles){
    track_test.push_back(offline);
    //ANA_MSG_INFO("Eta : " << offline -> eta() << "Phi : " << offline -> phi());
    Float_t mindr = 2000;
    hist ("h_phiInOffline")->Fill(offline -> phi());
    hist ("h_etaInOffline")->Fill(offline -> eta());
    
    
    ANA_CHECK (evtStore()->retrieve (truthparticles, "TruthParticles"));
    //ANA_MSG_INFO ("Number of truth particles " << truthparticles->size() );
    //truth loop
    for (const xAOD::TruthParticle* truth : *truthparticles) {
      if (truth->absPdgId() == 1000006) {
        
        ///ANA_MSG_INFO( m_nonSTOP << " number of nonstops since");
        m_nonSTOP = 0;

        ///ANA_MSG_INFO ("Truth particle pdgid pT nChildren  : " << truth->pdgId() << " " << truth->pt() << " " << truth->nChildren());    
        //ANA_MSG_INFO("Does it have prodvtx? " << truth->hasProdVtx() << " decayvtx? " << truth->hasDecayVtx());
        // check children        

        if (truth->nChildren() > 1) {
          for (int ichild=0; ichild< truth->nChildren() ; ichild++) {

            const xAOD::TruthParticle* child=truth->child(ichild);
            //ANA_MSG_INFO("child " << ichild << "  pdgid: " << child->pdgId() << " pT: " << child->pt() << " nChildren " << child->nChildren());	
          
            for (int igchild=0; igchild< child->nChildren() ; igchild++) {
              const xAOD::TruthParticle* gchild=child->child(igchild);
              if (gchild->absPdgId() == 13){
                //ANA_MSG_INFO("gchild " << igchild << "  pdgid: " << gchild->pdgId() << " nChildren " << gchild->nChildren() << " pT eta phi " << gchild->pt() << " " << gchild->eta() << " " << gchild->phi() );	
                //Decay length of parent
                const xAOD::TruthVertex* tproVtx = truth->prodVtx(); 
                const xAOD::TruthVertex* tdecVtx = truth->decayVtx(); 
                const xAOD::TruthVertex* cproVtx = child->prodVtx(); 
                const xAOD::TruthVertex* cdecVtx = child->decayVtx(); 

                

                if(mindr > calcd0(gchild, offline)){
                  mindr = calcd0(gchild, offline);
                }
                //ANA_MSG_INFO ("prod x, y, z: " << thevertex->x() << " " << thevertex->y() << " " << thevertex->z());
                //ANA_MSG_INFO ("perp: " << thevertex->perp());
                //ANA_MSG_INFO ("decay x, y, z: " << thevertex->x() << " " << thevertex->y() << " " << thevertex->z());
                //ANA_MSG_INFO ("perp: " << thevertex->perp());
                //ANA_MSG_INFO("Truth Decay length: "<< decaylength(tproVtx, tdecVtx));
                //ANA_MSG_INFO("Child Decay length: "<< decaylength(cproVtx, cdecVtx));
                hist ("h_truthDecayLength")->Fill (decaylength(tproVtx, tdecVtx));
                hist ("h_childDecayLength")->Fill (decaylength(cproVtx, cdecVtx));

                hist ("h_muon_parent")->Fill (gchild -> nParents());

                //Shove the grandchild up a vector
                //truth_vector.push_back(gchild);
        


              }
            }
          }    
        } 
      } else {
        m_nonSTOP = m_nonSTOP+1; 
      }
      /*
      if ( truth->nChildren() == 0 ) {
        ANA_MSG_INFO ("No children : Truth pdgid pT eta phi   : " << truth->pdgId() << " " << truth->pt() << " " << truth->eta() << " "  << truth->phi());
        if (truth->nParents() == 1) {

    ANA_MSG_INFO ("one parent : Truth pdgid pT eta phi   : " << truth->parent(0)->pdgId() << " " << truth->parent(0)->pt() << " " << truth->parent(0)->eta() << " "  << truth->parent(0)->phi());
        }
    } 
      */
    } // end loop over truth partciles

    //ANA_MSG_INFO("dr : " << mindr);
    hist ("h_drvalues")->Fill (mindr);
  }


  ////////////////////////////////////////////////////////////////////////////////

  ANA_MSG_INFO("Vector size of the pointer : " << track_test[0]);


  //Truth example
  


  
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

  ANA_MSG_INFO("ENDS YEET with size of vector test " << vector_test.size());
  return StatusCode::SUCCESS;
}