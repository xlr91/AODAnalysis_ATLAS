#include <AsgTools/MessageCheck.h>
#include <MyAnalysis/MyxAODAnalysis.h>
#include <xAODEventInfo/EventInfo.h>



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
  // declareProperty( "TitleforJobOption", codetitle = 0, "Desc?");
                   

}



StatusCode MyxAODAnalysis :: initialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  ANA_MSG_INFO ("in initialize");
  ANA_CHECK (book (TH1F ("h_jetPt", "h_jetPt", 100, 0, 500))); // jet pt [GeV]
  


  return StatusCode::SUCCESS;
}



StatusCode MyxAODAnalysis :: execute ()
{
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

  //Truth example
  const xAOD::TruthParticleContainer* truthparticles;
  ANA_CHECK (evtStore()->retrieve (truthparticles, "TruthParticles"));
  ANA_MSG_INFO ("Number of truth particles " << truthparticles->size() );
  
  for (const xAOD::TruthParticle* truth : *truthparticles) {
    if (truth->absPdgId() == 1000006) {
      
      ANA_MSG_INFO( m_nonSTOP << " number of nonstops since");
      m_nonSTOP = 0;

      ANA_MSG_INFO ("Truth particle pdgid pT nChildren  : " << truth->pdgId() << " " << truth->pt() << " " << truth->nChildren());    
      //ANA_MSG_INFO("Does it have prodvtx? " << truth->hasProdVtx() << " decayvtx? " << truth->hasDecayVtx());
      // check children        

      const xAOD::TruthVertex* thevertex = truth->prodVtx(); 
      //ANA_MSG_INFO ("prod x, y, z: " << thevertex->x() << " " << thevertex->y() << " " << thevertex->z());
      //ANA_MSG_INFO ("perp: " << thevertex->perp());

      thevertex = truth->decayVtx(); 
      //ANA_MSG_INFO ("decay x, y, z: " << thevertex->x() << " " << thevertex->y() << " " << thevertex->z());
      //ANA_MSG_INFO ("perp: " << thevertex->perp());
      

      if (truth->nChildren() > 1) {
        for (int ichild=0; ichild< truth->nChildren() ; ichild++) {

          const xAOD::TruthParticle* child=truth->child(ichild);
          ANA_MSG_INFO("child " << ichild << "  pdgid: " << child->pdgId() << " pT: " << child->pt() << " nChildren " << child->nChildren());	
        
          for (int igchild=0; igchild< child->nChildren() ; igchild++) {
            const xAOD::TruthParticle* gchild=child->child(igchild);
            ANA_MSG_INFO("gchild " << igchild << "  pdgid: " << gchild->pdgId() << " nChildren " << gchild->nChildren() << " pT eta phi " << gchild->pt() << " " << gchild->eta() << " " << gchild->phi() );	
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
  return StatusCode::SUCCESS;
}
