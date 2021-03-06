#include "vertexStudyLooper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include "../Tools/goodrun.h"
#include "../Tools/vtxreweight.h"
//#include "../Tools/pfjetMVAtools.h"

#include "../CORE/eventSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/susySelections.h"
#include "../CORE/muonSelections.h"
#include "../CORE/jetSelections.h"
#include "../CORE/metSelections.h"
#include "../CORE/jetSmearingTools.h"
#include "../CORE/jetcorr/JetCorrectionUncertainty.h"

bool verbose              = false;
bool doTenPercent         = false;
bool requireLeps          = true;
bool requireRecoTTbar     = true;
bool requireLepsWW        = false;
bool requireRecoLepsWW    = false;
bool doFiducialVtx        = true;

const float dzcut = 0.1;
const float drcut = 0.015;
const float dz_gen_vtx_match = 0.1;
const float fracpt_gen_vtx_match = 0.1;
const float fidvtx_cut = 15.;
const float gen_dr_match = 0.2;

using namespace std;
using namespace tas;

//--------------------------------------------------------------------

int getMotherIndex(int motherid){
  for(unsigned int i = 0; i < genps_id().size() ; i++){
    if( motherid == genps_id().at(i) ) return i;
  }

  return -1;
}

//--------------------------------------------------------------------

float vertexStudyLooper::dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 ) { 

  float dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  float deta = vec1.Eta() - vec2.Eta();

  return sqrt(dphi*dphi + deta*deta);

}

//--------------------------------------------------------------------

vertexStudyLooper::vertexStudyLooper()
{

  std::cout << " construct " << std::endl;
  g_createTree   = false;
  initialized = false;
}

//--------------------------------------------------------------------

struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

//--------------------------------------------------------------------

bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

//--------------------------------------------------------------------

std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------

// void vertexStudyLooper::InitBaby(){
// }

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void vertexStudyLooper::closeOutput()
{
  outFile->cd();
  //  outTree->Write();
  outFile->Write();
  outFile->Close();
  delete outFile;
}

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */

int vertexStudyLooper::ScanChain(TChain* chain, const TString& prefix)

{

  //  cout << "ciao " << isData << endl;

  bool isData = false;
  if( prefix.Contains("data") || prefix.Contains("2012") 
      || prefix.Contains("dimu") || prefix.Contains("diel")
      || prefix.Contains("mueg") ){
    cout << "DATA!!!" << endl;
    isData       = true;
    doTenPercent = false;
  }

  cout << "IS DATA: " << isData << endl;

  if( doTenPercent ) cout << "Processing 10% of MC" << endl;

  //------------------------------------------------------------------------------------------------------
  // set json\, vertex reweighting function and msugra cross section files
  //------------------------------------------------------------------------------------------------------
  
  if( !initialized ){

    //set json
    cout << "setting json " << g_json << endl;
    set_goodrun_file( g_json );

    //    if( prefix.Contains("ttall_massivebin") ) 
    set_vtxreweight_rootfile("vtxreweight/vtxreweight_Summer12MC_PUS10_19fb_Zselection.root",true);

    initialized = true;
  }

  //------------------------------------------------------------------------------------------------------
  // latest-and-greatest JEC
  //------------------------------------------------------------------------------------------------------

  // std::vector<std::string> jetcorr_filenames_pfL1FastJetL2L3;
  // FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3;

  // jetcorr_filenames_pfL1FastJetL2L3.clear();
  
  // string pfUncertaintyFile;
  // //string caloUncertaintyFile;

  // if ( isData ) {
  //   jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L1FastJet_AK5PF.txt");
  //   jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L2Relative_AK5PF.txt");
  //   jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L3Absolute_AK5PF.txt");
  //   jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/GR_P_V39_AN3_L2L3Residual_AK5PF.txt");

  //   pfUncertaintyFile = "jetCorrections/GR_P_V39_AN3_Uncertainty_AK5PF.txt";
  // } 
  // else {
  //   jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/START53_V21_L1FastJet_AK5PF.txt");
  //   jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/START53_V21_L2Relative_AK5PF.txt");
  //   jetcorr_filenames_pfL1FastJetL2L3.push_back  ("jetCorrections/START53_V21_L3Absolute_AK5PF.txt");
    
  //   pfUncertaintyFile = "jetCorrections/START53_V21_Uncertainty_AK5PF.txt";
  // }

  // jet_corrector_pfL1FastJetL2L3  = makeJetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

  // JetCorrectionUncertainty *pfUncertainty   = new JetCorrectionUncertainty( pfUncertaintyFile   );

  // MetCorrector *met_corrector_pfL1FastJetL2L3 = new MetCorrector(jetcorr_filenames_pfL1FastJetL2L3);

  // /*
  //  *  Jet Smearer Object to obtain the jet pt uncertainty.
  //  */

  // std::vector<std::string> list_of_file_names;
  // list_of_file_names.push_back("jetSmearData/Spring10_PtResolution_AK5PF.txt");
  // list_of_file_names.push_back("jetSmearData/Spring10_PhiResolution_AK5PF.txt");
  // list_of_file_names.push_back("jetSmearData/jet_resolutions.txt");
  // JetSmearer *jetSmearer = makeJetSmearer(list_of_file_names);
 
  // -----------------------------------------------------------

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  makeOutput(prefix);
  //  if(g_createTree) makeTree(prefix, doFakeApp, frmode);

  //  BookHistos(prefix);
  BookHistos("h");

  cout << " done with initialization "  << endl;
  
  unsigned int nEventsChain = chain->GetEntries();
  unsigned int nEventsTotal = 0;
  unsigned int nEventsPreReco = 0;
  unsigned int nEventsPass = 0;
  // map isn't needed for this purpose, vector is sufficient
  // better would be to use a struct with run, lb, event
  map<int,int> m_events;

  // loop over files
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TChainElement* currentFile = 0;

  // test chain
  if (!chain)
    {
      throw std::invalid_argument("at::ScanChain: chain is NULL!");
    }
  if (chain->GetListOfFiles()->GetEntries()<1)
    {
      throw std::invalid_argument("at::ScanChain: chain has no files!");
    }
  if (not chain->GetFile())
    {
      throw std::invalid_argument("at::ScanChain: chain has no files or file path is invalid!");
    }
  int nSkip_els_conv_dist = 0;

  while((currentFile = (TChainElement*)fileIter.Next())) {
    TFile* f = new TFile(currentFile->GetTitle());

    cout << currentFile->GetTitle() << endl;

    if (!f || f->IsZombie()) {
      throw std::runtime_error(Form("ERROR::File from TChain is invalid or corrupt: %s", currentFile->GetTitle()));
    }
    
    // get the trees in each file
    // TTree *tree = (TTree*)f->Get("Events");
    TTree *tree = dynamic_cast<TTree*>(f->Get("Events"));
    if (!tree || tree->IsZombie()) {
      throw std::runtime_error(Form("ERROR::File from TChain has an invalid TTree or is corrupt: %s", currentFile->GetTitle()));
    }

    //Matevz
    TTreeCache::SetLearnEntries(100);
    tree->SetCacheSize(128*1024*1024);

    cms2.Init(tree);
      
    unsigned int nEntries = tree->GetEntries();
    for(unsigned int z = 0; z < nEntries; ++z) {
      ++nEventsTotal;

      /////////      cout << nEventsTotal << endl;

      if( doTenPercent ){
	if( !(nEventsTotal%10==0) ) continue;
      }

      // progress feedback to user
      if (nEventsTotal % 1000 == 0){
        
        // xterm magic from L. Vacavant and A. Cerri
        if (isatty(1)){
                
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
          fflush(stdout);
        }
      }

      //Matevz
      tree->LoadTree(z);

      cms2.GetEntry(z);

      if( evt_ww_rho_vor() != evt_ww_rho_vor() ){
	cout << "Skipping event with rho = nan!!!" << endl;
	continue;
      }

      //      InitBaby();

      isdata_ = isData ? 1 : 0;

      if( verbose ){
	cout << "-------------------------------------------------------"   << endl;
	cout << "Event " << z                                               << endl;
	cout << "File  " << currentFile->GetTitle()                         << endl;
	cout << evt_dataset().at(0) << " " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << endl;
	cout << "-------------------------------------------------------"   << endl;
      }

      TString datasetname(evt_dataset().at(0));

      //---------------------------------------------
      // event cleaning and good run list
      //---------------------------------------------

      //      if( !cleaning_goodVertexApril2011() )                          continue;
      if( isData && !goodrun(evt_run(), evt_lumiBlock()) ) continue;

      //---------------------
      // skip duplicates
      //---------------------

      if( isData ) {
        DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
        if (is_duplicate(id) ){
          continue;
        }
      }

      //-------------------------------------
      // skip events with bad els_conv_dist
      //-------------------------------------

      bool skipEvent = false;
      for( unsigned int iEl = 0 ; iEl < els_conv_dist().size() ; ++iEl ){
        if( els_conv_dist().at(iEl) != els_conv_dist().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEta().at(iEl) != els_sigmaIEtaIEta().at(iEl) ){
          skipEvent = true;
        }
        if( els_sigmaIEtaIEtaSC().at(iEl) != els_sigmaIEtaIEtaSC().at(iEl) ){
          skipEvent = true;
        }
      }
             
      if( skipEvent ){
        nSkip_els_conv_dist++;
        continue;
      }
   
      int ngenleps = 0;
      int ngenleps10 = 0;
      int ngenleps_etaacc = 0;
      int npartons = 0;
      int npartons_etaacc = 0;
      // loop over genps, count e and mu with pt > 20 gev
      for (unsigned int igen = 0; igen < genps_status().size(); ++igen) {
	if (genps_status().at(igen) != 3) continue;
	if (genps_p4().at(igen).pt() < 10.) continue;
	int id = abs(genps_id().at(igen));
	if (id == 11 || id == 13) {
	  ++ngenleps10;
	}
	if (genps_p4().at(igen).pt() < 20.) continue;
	bool passeta = false;
	if (fabs(genps_p4().at(igen).eta()) < 2.4) passeta = true;
	if (id == 11 || id == 13) {
	  ++ngenleps;
	  if (passeta) ++ngenleps_etaacc;
	} // leptons
	if ((id < 6) || (id == 21)) {
	  ++npartons;
	  if (passeta) ++npartons_etaacc;
	} // partons
      }

      if (requireLeps && !ngenleps) continue;
      if (requireLepsWW && (ngenleps10 < 2) && (ngenleps < 1)) continue;

      // check gen production vtx by looking at status 3 particles
      // and make fiducial cut if requested
      float genvtx_z = genps_prod_vtx().at(2).z();
      h_genvtx_z_nocut->Fill(genvtx_z);
      if (doFiducialVtx && (fabs(genvtx_z) > fidvtx_cut)) continue;

      ++nEventsPreReco;

      //---------------------------------------------
      // reco electron selection
      //---------------------------------------------

      std::vector<int> selected_el_idx;
      int nel20 = 0;
      int nel10 = 0;
      for (unsigned int iel = 0; iel < els_p4().size(); ++iel) {
	// cut on pt and eta
	if (els_p4().at(iel).pt() < 10.) continue;
	if (fabs(els_p4().at(iel).eta()) > 2.4) continue;
	// require ID, no iso
	// !!! note that the dz, dxy cuts are commented out in CORE to avoid biases!
	if( !passElectronSelection_Stop2012_v3_NoIso( iel,true,true,false) )  continue;
	// check for gen match: match to gen particle within dR < 0.2
	bool matched = false;
	for (unsigned int igen = 0; igen < genps_id().size(); ++igen) {
	  int id = fabs(genps_id().at(igen));
	  if (id != 11) continue;
	  float dr = dRbetweenVectors(els_p4().at(iel),genps_p4().at(igen));
	  if (dr < gen_dr_match) {
	    matched = true;
	    break;
	  }
	}
	if (!matched) continue;
	++nel10;
	selected_el_idx.push_back(iel);
	if (els_p4().at(iel).pt() > 20.) ++nel20;
      }

      //---------------------------------------------
      // reco muon selection
      //---------------------------------------------

      std::vector<int> selected_mu_idx;
      int nmu20 = 0;
      int nmu10 = 0;
      for (unsigned int imu = 0; imu < mus_p4().size(); ++imu) {
	// cut on pt and eta
	if (mus_p4().at(imu).pt() < 10.) continue;
	if (fabs(mus_p4().at(imu).eta()) > 2.4) continue;
	// require ID, no iso
	// !!! note that the dz, dxy cuts are commented out in CORE to avoid biases!
	if (!muonIdNotIsolated(imu, ZMet2012_v1)) continue;
	// check for gen match: match to gen particle within dR < 0.2
	bool matched = false;
	for (unsigned int igen = 0; igen < genps_id().size(); ++igen) {
	  int id = fabs(genps_id().at(igen));
	  if (id != 13) continue;
	  float dr = dRbetweenVectors(mus_p4().at(imu),genps_p4().at(igen));
	  if (dr < gen_dr_match) {
	    matched = true;
	    break;
	  }
	}
	if (!matched) continue;
	++nmu10;
	selected_mu_idx.push_back(imu);
	if (mus_p4().at(imu).pt() > 20.) ++nmu20;
      }

      //---------------------------------------------
      // reco jet selection
      //---------------------------------------------

      int njets30 = 0;
      for (unsigned int ijet = 0; ijet < pfjets_p4().size(); ++ijet) {
	// pt
	if (pfjets_p4().at(ijet).pt() < 30.) continue;
	// eta cut for tracker
	if (fabs(pfjets_p4().at(ijet).eta()) > 2.4) continue;
	// pf jet id
	if (!passesPFJetID(ijet)) continue;
	++njets30;
      }

      //---------------------------------------------
      // reco selections
      //---------------------------------------------

      int nleps10 = nel10 + nmu10;
      int nleps20 = nel20 + nmu20;
      if (requireRecoLepsWW && ((nleps10 < 2) || (nleps20 < 1)) ) continue;
      // select ttbar: hadronic (5 jets), single lep (1l+3j), or dilep (2l+2j)
      bool pass_ttbar_had = (njets30 >= 5);
      bool pass_ttbar_1l = (nleps20 >= 1 && njets30 >= 3);
      bool pass_ttbar_2l = (nleps20 >= 2 && njets30 >= 2);
      if (requireRecoTTbar && !(pass_ttbar_had || pass_ttbar_1l || pass_ttbar_2l)) continue;

      ++nEventsPass;
      float genvtx_sumpt2 = genSumPt2();

      // make "lepton type":
      // 0 for hadronic, no leptons
      // 1 for single e
      // 2 for single mu
      // 3 for ee
      // 4 for mm
      // 5 for em
      int leptype = 0;
      if (nel20 >= 1 && nel10 >= 2) leptype = 3;
      else if (nmu20 >= 1 && nmu10 >= 2) leptype = 4;
      else if ((nmu20 >= 1 && nel10 >= 1) || (nel20 >= 1 && nmu10 >= 1)) leptype = 5;
      else if (nel20 == 1) leptype = 1;
      else if (nmu20 == 1) leptype = 2;


      //---------------------------------------------
      // loop over tracks, associate to vtx
      //---------------------------------------------

      std::vector<int> trk_bestdzvtx(trks_trk_p4().size(), -1);
      std::vector<float> vtxs_sumpt_recalc_dz(vtxs_sumpt().size(), 0.);
      std::vector<float> vtxs_sumpt_recalc_weight(vtxs_sumpt().size(), 0.);
      std::vector<float> vtxs_sumpt_hardscatter_dz(vtxs_sumpt().size(), 0.);
      std::vector<float> vtxs_sumpt_hardscatter_weight(vtxs_sumpt().size(), 0.);
      std::vector<float> vtxs_fracpt_hardscatter(vtxs_sumpt().size(), -1.);
      std::vector<float> vtxs_sumpt2_weight(vtxs_sumpt().size(), 0.);
      std::vector<float> vtxs_sumpt2_nogen_weight(vtxs_sumpt().size(), 0.);
      std::map<int,std::vector<int> > mc_idx_reco_matches;
      int mc_idx_duplicates = 0;
      float sum_hardscatter_pt = 0.;
      int ntracks = 0;

      for (unsigned int itrk = 0; itrk < trks_trk_p4().size(); ++itrk) {
	// only consider high purity tracks: bit 2 of quality mask
	if (!(trks_qualityMask().at(itrk) & (1<<2))) continue;
	++ntracks;

	// find best match vertex, from vertex algo weight and simple dz
	int bestdzvtx = associateTrackToVertex(itrk);
	//int bestdzvtx = associateTrackToVertexSimple(itrk);
	int weightvtx = trks_pvidx0().at(itrk);
	float dz = trks_dz_pv(itrk, bestdzvtx).first;
	//float dz = cms2.trks_z0().at(itrk) - vtxs_position().at(bestdzvtx).z();
	h_dz_trk_vtx->Fill(dz);
	if (fabs(dz) <= dzcut) {
	  h_trk_bestdzvtx->Fill(bestdzvtx);
	  h_trk_bestdzvtx_vs_weightvtx->Fill(trks_pvidx0().at(itrk),bestdzvtx);
	  if ((weightvtx > -9000) && (bestdzvtx != weightvtx)) {
	    h_trk_dz_vtxs->Fill(vtxs_position().at(weightvtx).z() - vtxs_position().at(bestdzvtx).z());
	  }
	  vtxs_sumpt_recalc_dz.at(bestdzvtx) += trks_trk_p4().at(itrk).pt();
	} else {
	  // closest vertex fails dzcut: set bestdzvtx to -1 for invalid
	  bestdzvtx = -1;
	}
	trk_bestdzvtx.at(itrk) = bestdzvtx;

	if (weightvtx > -9000) {
	  vtxs_sumpt_recalc_weight.at(weightvtx) += trks_trk_p4().at(itrk).pt();
	  vtxs_sumpt2_weight.at(weightvtx) += pow(trks_trk_p4().at(itrk).pt(),2);
	  if(weightvtx == 0) h_dz_trk_vtx0_weight->Fill(trks_dz_pv(itrk,0).first);
	}

	// check that track is matched to status 1 particle from hard scatter
	// apply cut on dR of reco/gen match to remove most multiply matched tracks
	if ((trk_mc_id().at(itrk) == -9999) || (trk_mcdr().at(itrk) > drcut)) {
	  // pileup track distributions
	  h_nomatch_trk_pt->Fill(trks_trk_p4().at(itrk).pt());
	  h_nomatch_trk_pt_low->Fill(trks_trk_p4().at(itrk).pt());
	  h_nomatch_trk_pt_mid->Fill(trks_trk_p4().at(itrk).pt());
	  h_nomatch_trk_pt_vs_eta->Fill(trks_trk_p4().at(itrk).eta(),trks_trk_p4().at(itrk).pt());
	  h_nomatch_trk_pt_low_vs_eta->Fill(trks_trk_p4().at(itrk).eta(),trks_trk_p4().at(itrk).pt());
	  if (weightvtx > -9000) vtxs_sumpt2_nogen_weight.at(weightvtx) += pow(trks_trk_p4().at(itrk).pt(),2);
	  if (trk_mc_id().at(itrk) != -9999) h_match_dr->Fill(trk_mcdr().at(itrk));
	  continue;
	}

	// check for tracks matched to the same gen particle..
	if (mc_idx_reco_matches.count(trk_mcidx().at(itrk)) > 0) {
	  ++mc_idx_duplicates;
	  mc_idx_reco_matches[trk_mcidx().at(itrk)].push_back(itrk);
	}
	else {
	  std::vector<int> reco_matches(1,itrk);
	  mc_idx_reco_matches[trk_mcidx().at(itrk)] = reco_matches;
	}
      } // loop over tracks

      h_ntracks->Fill(ntracks);

      //---------------------------------------------
      // loop over matched gen particles
      //---------------------------------------------

      h_mc_idx_duplicates->Fill(mc_idx_duplicates);
      // loop again to decide on duplicate matches, save/plot quantities for hard scatter
      std::map<int,std::vector<int> >::const_iterator mc_itr;
      for (mc_itr = mc_idx_reco_matches.begin(); mc_itr != mc_idx_reco_matches.end(); ++mc_itr) {
	// check for duplicate matches
	//  if unique match, just add to the hardscatter sums
	if ((*mc_itr).second.size() == 1) {
	  int itrk = (*mc_itr).second.at(0);
	  int weightvtx = trks_pvidx0().at(itrk);
	  int bestdzvtx = trk_bestdzvtx.at(itrk);

	  if (bestdzvtx > -1) vtxs_sumpt_hardscatter_dz.at(bestdzvtx) += trks_trk_p4().at(itrk).pt();
	  if (weightvtx > -9000) {
	    vtxs_sumpt_hardscatter_weight.at(weightvtx) += trks_trk_p4().at(itrk).pt();
	    if (weightvtx != 0) {
	      h_hs_trk_dz_vtxs->Fill(vtxs_position().at(weightvtx).z() - vtxs_position().at(0).z());
	      h_hs_trk_dz_vtx0->Fill(trks_dz_pv(itrk, 0).first);
	    }
	  }
	  sum_hardscatter_pt += trks_trk_p4().at(itrk).pt();
	  h_match_trk_pt->Fill(trks_trk_p4().at(itrk).pt());
	  h_match_trk_pt_low->Fill(trks_trk_p4().at(itrk).pt());
	  h_match_trk_pt_mid->Fill(trks_trk_p4().at(itrk).pt());
	} 
	// if multiple matches, loop over matched tracks and select based on most hits
	else {
	  int maxnhits = 0;
	  int maxnhits_idx = -1;
	  float mindr = 99.;
	  int mindr_idx = -1;
	  for (unsigned int imatch = 0; imatch < (*mc_itr).second.size(); ++imatch) {
	    int itrk = (*mc_itr).second.at(imatch);
	    h_mc_idx_duplicates_dr->Fill(trk_mcdr().at(itrk));
	    h_mc_idx_duplicates_pt->Fill(trks_trk_p4().at(itrk).pt());
	    h_mc_idx_duplicates_nhits->Fill(trks_validHits().at(itrk));

	    if (trks_validHits().at(itrk) > maxnhits) {
	      maxnhits = trks_validHits().at(itrk);
	      maxnhits_idx = itrk;
	    }
	    if (trk_mcdr().at(itrk) < mindr) { 
	      mindr = trk_mcdr().at(itrk);
	      mindr_idx = itrk;
	    }
	  } // loop on duplicate matches

	  int ibestmatch = maxnhits_idx;
	  int weightvtx = trks_pvidx0().at(ibestmatch);
	  int bestdzvtx = trk_bestdzvtx.at(ibestmatch);

	  if (bestdzvtx > -1) vtxs_sumpt_hardscatter_dz.at(bestdzvtx) += trks_trk_p4().at(ibestmatch).pt();
	  if (weightvtx > -9000) {
	    vtxs_sumpt_hardscatter_weight.at(weightvtx) += trks_trk_p4().at(ibestmatch).pt();
	    if (weightvtx != 0) {
	      h_hs_trk_dz_vtxs->Fill(vtxs_position().at(weightvtx).z() - vtxs_position().at(0).z());
	      h_hs_trk_dz_vtx0->Fill(trks_dz_pv(ibestmatch, 0).first);
	    }
	  }
	  sum_hardscatter_pt += trks_trk_p4().at(ibestmatch).pt();
	  h_match_trk_pt->Fill(trks_trk_p4().at(ibestmatch).pt());
	  h_match_trk_pt_low->Fill(trks_trk_p4().at(ibestmatch).pt());
	  h_match_trk_pt_mid->Fill(trks_trk_p4().at(ibestmatch).pt());
	} // if duplicate matches

      } // loop over matched gen particles

      //---------------------------------------------
      // loop over vertices
      //---------------------------------------------

      int nvtx = 0;
      int gen_match_vtx = -1;
      float mindz = 99.;
      int mindz_idx = -1;
      float mindz_eff10 = 99.;
      int mindz_eff10_idx = -1;
      for (unsigned int ivtx = 0; ivtx < vtxs_sumpt().size(); ++ivtx) {
	// count good vertices
	if (!isGoodVertex(ivtx)) continue;
	++nvtx;
	h_recovtx_z->Fill(vtxs_position().at(ivtx).z());
	h_vtx_sumpt2->Fill(vtxs_sumpt2_weight.at(ivtx));
	h_vtx_nogen_sumpt2->Fill(vtxs_sumpt2_nogen_weight.at(ivtx));

	// check to see if hard scatter vertex is among the reco vtx collection
	// if already matched to a higher sum pt^2 vertex, don't bother
	if (gen_match_vtx > -1) continue;
	// require: 
	//  |dz| < 1mm from true hard scatter vtx
	//  at least 10% of hard scatter pt
	float dz = vtxs_position().at(ivtx).z() - genvtx_z;
	if (fabs(dz) < mindz) {
	  mindz = fabs(dz);
	  mindz_idx = ivtx;
	}
	float purity_dz = vtxs_sumpt_hardscatter_dz.at(ivtx)/vtxs_sumpt_recalc_dz.at(ivtx);
	float purity_weight = vtxs_sumpt_hardscatter_weight.at(ivtx)/vtxs_sumpt_recalc_weight.at(ivtx);
	float eff_dz = vtxs_sumpt_hardscatter_dz.at(ivtx)/sum_hardscatter_pt;
	float eff_weight = vtxs_sumpt_hardscatter_weight.at(ivtx)/sum_hardscatter_pt;

	//	if (eff_dz < fracpt_gen_vtx_match) continue;

	if (fabs(dz) < mindz_eff10) {
	  mindz_eff10 = fabs(dz);
	  mindz_eff10_idx = ivtx;
	}

	if (fabs(dz) > dz_gen_vtx_match ) continue;
	//	if (vtxs_sumpt_hardscatter_dz.at(ivtx)/vtxs_sumpt_recalc_dz.at(ivtx) < fracpt_gen_vtx_match ) continue;

	h_vtx_best_purity_dz->Fill(purity_dz);
	h_vtx_best_purity_dz_vs_nvtx->Fill(nvtx,purity_dz);
	h_vtx_best_purity_weight->Fill(purity_weight);
	h_vtx_best_purity_weight_vs_nvtx->Fill(nvtx,purity_weight);

	h_vtx_best_eff_dz->Fill(eff_dz);
	h_vtx_best_eff_dz_vs_nvtx->Fill(nvtx,eff_dz);
	h_vtx_best_eff_weight->Fill(eff_weight);
	h_vtx_best_eff_weight_vs_nvtx->Fill(nvtx,eff_weight);

	h_genmatch_vtx_z->Fill(vtxs_position().at(ivtx).z());
	h_genmatch_vtx_sumpt2->Fill(vtxs_sumpt2_weight.at(ivtx));

	h_genmatch_leptype->Fill(leptype);

	// low sumpt2: check for poor electron tracks
	if (vtxs_sumpt2_weight.at(ivtx) < 500.) {
	h_genmatch_lowsumpt2_leptype->Fill(leptype);
	  for (unsigned int iel=0; iel < selected_el_idx.size(); ++iel) {
	    int trkidx = els_trkidx().at(selected_el_idx.at(iel));
	    if (trkidx >= 0) {
	      int weightvtx = trks_pvidx0().at(trkidx);
	      if (weightvtx == (int)ivtx) h_genmatch_lowsumpt2_eltrkassoc->Fill(1);
	      else {
		// check if track passes high purity selection
		if (!(trks_qualityMask().at(trkidx) & (1<<2))) h_genmatch_lowsumpt2_eltrkassoc->Fill(-1);
		else h_genmatch_lowsumpt2_eltrkassoc->Fill(0);
	      }
	    } else {
	      // invalid trkref for this electron
	      h_genmatch_lowsumpt2_eltrkassoc->Fill(-2);
	    }
	  } // selected el loop
	} // low sumpt2

	// if not matched to vtx0
	if (ivtx > 0) {
	  h_genvtx_othermatch_z->Fill(genvtx_z);
	  h_genvtx_othermatch_gensumpt2->Fill(genvtx_sumpt2);
	  h_genvtx_othermatch_sumpt2->Fill(vtxs_sumpt2_weight.at(ivtx));
	  h_genvtx_othermatch_leptype->Fill(leptype);
	}

	gen_match_vtx = ivtx;
      }
      h_nvtx->Fill(nvtx);
      h_nvtx_vs_ntrueint->Fill(puInfo_nPUvertices().at(0),nvtx);
      h_gen_match_vtx->Fill(gen_match_vtx);
      h_gen_match_vtx_vs_nvtx->Fill(nvtx,gen_match_vtx);

      h_genvtx_z->Fill(genvtx_z);
      h_genvtx_reco_dz->Fill(mindz);

      // no matched vertex
      if (gen_match_vtx == -1) {
	h_genvtx_nomatch_z->Fill(genvtx_z);
	h_genvtx_nomatch_nvtx->Fill(nvtx);
	h_genvtx_nomatch_ntracks->Fill(ntracks);
	std::cout << "-- Event with no matched vertex: " << evt_run() << " " << evt_lumiBlock() << " " << evt_event() << std::endl
		  << "     genvtx z: " << genvtx_z << ", genvtx sumpt2: " << genvtx_sumpt2 << std::endl
		  << "     nvtx: " << nvtx << ", ntracks: " << ntracks << std::endl;
	if (mindz_idx > -1) {
	  h_genvtx_nomatch_dz->Fill(genvtx_z - vtxs_position().at(mindz_idx).z());
	  std::cout << "     closest vertex dz: " << mindz << std::endl;
	  std::cout << "     ngenleps: " << ngenleps << ", in eta acc: " << ngenleps_etaacc << std::endl;
	  std::cout << "     npartons: " << npartons << ", in eta acc: " << npartons_etaacc << std::endl;

	  if (mindz < 0.2) {
	    h_genvtx_nomatch_smalldz_eff->Fill(vtxs_sumpt_hardscatter_dz.at(mindz_idx)/sum_hardscatter_pt);
	  } else {
	    // large dz: check how many leptons and partons were within pt/eta acceptance
	    h_genvtx_nomatch_largedz_nleps->Fill(ngenleps_etaacc);
	    h_genvtx_nomatch_largedz_npartons->Fill(npartons_etaacc);
	  }
	} // if mindz_idz > -1
	if (mindz_eff10_idx > -1) h_genvtx_nomatch_dz_eff10->Fill(genvtx_z - vtxs_position().at(mindz_eff10_idx).z());
	h_genvtx_nomatch_gensumpt2->Fill(genvtx_sumpt2);
      }

      // loop back over vertices to make plots for all except hard scatter vtx
      for (unsigned int ivtx = 0; ivtx < vtxs_sumpt().size(); ++ivtx) {
	// only good vertices
	if (!isGoodVertex(ivtx)) continue;
	if (int(ivtx) == gen_match_vtx) continue;
	h_vtx_nohs_sumpt2->Fill(vtxs_sumpt2_weight.at(ivtx));
	h_vtx_nogen_nohs_sumpt2->Fill(vtxs_sumpt2_nogen_weight.at(ivtx));
      }

      //---------------------------------------------
      // plots for vertex 0
      //---------------------------------------------

      const int vtx0 = 0;
      h_vtx0_genvtx_dz->Fill(genvtx_z - vtxs_position().at(vtx0).z());
      float purity_dz = vtxs_sumpt_hardscatter_dz.at(vtx0)/vtxs_sumpt_recalc_dz.at(vtx0);
      float purity_weight = vtxs_sumpt_hardscatter_weight.at(vtx0)/vtxs_sumpt_recalc_weight.at(vtx0);
      h_vtx0_purity_dz->Fill(purity_dz);
      h_vtx0_purity_dz_vs_nvtx->Fill(nvtx,purity_dz);
      h_vtx0_purity_weight->Fill(purity_weight);
      h_vtx0_purity_weight_vs_nvtx->Fill(nvtx,purity_weight);

      float eff_dz = vtxs_sumpt_hardscatter_dz.at(vtx0)/sum_hardscatter_pt;
      float eff_weight = vtxs_sumpt_hardscatter_weight.at(vtx0)/sum_hardscatter_pt;
      h_vtx0_eff_dz->Fill(eff_dz);
      h_vtx0_eff_dz_vs_nvtx->Fill(nvtx,eff_dz);
      h_vtx0_eff_weight->Fill(eff_weight);
      h_vtx0_eff_weight_vs_nvtx->Fill(nvtx,eff_weight);

      h_vtx0_hardscatter_pt_vs_sumpt->Fill(vtxs_sumpt().at(vtx0),vtxs_sumpt_hardscatter_weight.at(vtx0));
      h_vtx0_hardscatter_pt_vs_sumpt_recalc->Fill(vtxs_sumpt_recalc_weight.at(vtx0),vtxs_sumpt_hardscatter_weight.at(vtx0));
      h_vtx0_sumpt_vs_sumpt_recalc->Fill(vtxs_sumpt_recalc_weight.at(vtx0),vtxs_sumpt().at(vtx0));

      if (gen_match_vtx != vtx0) {
	h_vtx0_pu_sumpt2->Fill(vtxs_sumpt2_weight.at(vtx0));
	if (gen_match_vtx > -1) {
	  h_genvtx_othermatch_dz->Fill(vtxs_position().at(gen_match_vtx).z() - vtxs_position().at(vtx0).z());
	}
      }

      //---------------------------------------------
      // loop over selected electrons to make plots
      //---------------------------------------------

      for (unsigned int iel = 0; iel < selected_el_idx.size(); ++iel) {
	unsigned int el_idx = selected_el_idx.at(iel);
	// plot isolation with and without pileup correction, also vs vtx0 purity
	float iso = electronPFiso(el_idx);
	float iso_cor = electronPFiso(el_idx,true);
	float iso40_cor = electronPFiso40(el_idx,true);
	float trkiso = els_iso03_pf2012ext_ch().at(el_idx)/els_p4().at(el_idx).pt();
	float trkiso_abs = els_iso03_pf2012ext_ch().at(el_idx);
	float trkiso40 = els_iso04_pf2012ext_ch().at(el_idx)/els_p4().at(el_idx).pt();
	float trkiso40_abs = els_iso04_pf2012ext_ch().at(el_idx);

	h_el_iso->Fill(iso);
	h_el_iso_cor->Fill(iso_cor);
	h_el_iso_vs_vtx0_purity_dz->Fill(purity_dz,iso);
	h_el_iso_cor_vs_vtx0_purity_dz->Fill(purity_dz,iso_cor);
	h_el_iso_vs_nvtx->Fill(nvtx,iso);
	h_el_iso_cor_vs_nvtx->Fill(nvtx,iso_cor);
	h_el_iso40_cor_vs_nvtx->Fill(nvtx,iso40_cor);

	h_el_trkiso->Fill(trkiso);
	h_el_trkiso_vs_vtx0_purity_dz->Fill(purity_dz,trkiso);
	h_el_trkiso_vs_nvtx->Fill(nvtx,trkiso);
	h_el_trkiso40_vs_nvtx->Fill(nvtx,trkiso40);
	h_el_trkiso_abs->Fill(trkiso_abs);
	h_el_trkiso_abs_vs_vtx0_purity_dz->Fill(purity_dz,trkiso_abs);
	h_el_trkiso_abs_vs_nvtx->Fill(nvtx,trkiso_abs);
	h_el_trkiso40_abs_vs_nvtx->Fill(nvtx,trkiso40_abs);
      }

      //---------------------------------------------
      // loop over selected muons to make plots
      //---------------------------------------------

      for (unsigned int imu = 0; imu < selected_mu_idx.size(); ++imu) {
	unsigned int mu_idx = selected_mu_idx.at(imu);
	// plot isolation with and without pileup correction, also vs vtx0 purity
	float iso = muonPFiso(mu_idx);
	float iso_cor = muonPFiso(mu_idx,true);
	float iso40_cor = muonPFiso40(mu_idx,true);
	float trkiso = mus_isoR03_pf_ChargedHadronPt().at(mu_idx)/mus_p4().at(mu_idx).pt();
	float trkiso_abs = mus_isoR03_pf_ChargedHadronPt().at(mu_idx);
	float trkiso40 = mus_isoR04_pf_ChargedHadronPt().at(mu_idx)/mus_p4().at(mu_idx).pt();
	float trkiso40_abs = mus_isoR04_pf_ChargedHadronPt().at(mu_idx);

	h_mu_iso->Fill(iso);
	h_mu_iso_cor->Fill(iso_cor);
	h_mu_iso_vs_vtx0_purity_dz->Fill(purity_dz,iso);
	h_mu_iso_cor_vs_vtx0_purity_dz->Fill(purity_dz,iso_cor);
	h_mu_iso_vs_nvtx->Fill(nvtx,iso);
	h_mu_iso_cor_vs_nvtx->Fill(nvtx,iso_cor);
	h_mu_iso40_cor_vs_nvtx->Fill(nvtx,iso40_cor);

	h_mu_trkiso->Fill(trkiso);
	h_mu_trkiso_vs_vtx0_purity_dz->Fill(purity_dz,trkiso);
	h_mu_trkiso_vs_nvtx->Fill(nvtx,trkiso);
	h_mu_trkiso40_vs_nvtx->Fill(nvtx,trkiso40);
	h_mu_trkiso_abs->Fill(trkiso_abs);
	h_mu_trkiso_abs_vs_vtx0_purity_dz->Fill(purity_dz,trkiso_abs);
	h_mu_trkiso_abs_vs_nvtx->Fill(nvtx,trkiso_abs);
	h_mu_trkiso40_abs_vs_nvtx->Fill(nvtx,trkiso40_abs);
      }

      //---------------------------------------------
      // loop over photons
      //---------------------------------------------

      for (unsigned int iph = 0; iph < photons_p4().size(); ++iph) {
	// cut on pt and eta
	if (photons_p4().at(iph).pt() < 20.) continue;
	if (fabs(photons_p4().at(iph).eta()) > 2.4) continue;
	// require ID -- not sure what to use here
	//        if( !passElectronSelection_Stop2012_v3_NoIso( iph,true,true,false) )  continue;
	// check for gen match: match to gen particle within dR < 0.2
	bool matched = false;
	for (unsigned int igen = 0; igen < genps_id().size(); ++igen) {
	  int id = fabs(genps_id().at(igen));
	  if (id != 22) continue;
	  float dr = dRbetweenVectors(photons_p4().at(iph),genps_p4().at(igen));
	  if (dr < gen_dr_match) {
	    matched = true;
	    break;
	  }
	}
	if (!matched) continue;
	// plot (hollow) track isolation / pt
	//	float iso = photons_tkIsoHollow03().at(iph)/photons_p4().at(iph).pt();
	float iso = photonHollowTrkIso(iph)/photons_p4().at(iph).pt();

	h_ph_trkiso->Fill(iso);
	h_ph_trkiso_vs_vtx0_purity_dz->Fill(purity_dz,iso);
	h_ph_trkiso_vs_nvtx->Fill(nvtx,iso);
      }

      //---------------------------------------------
      // loop over jets
      //---------------------------------------------

      for (unsigned int ijet = 0; ijet < pfjets_p4().size(); ++ijet) {
	// pt
	if (pfjets_p4().at(ijet).pt() < 30.) continue;
	// eta cut for tracker
	if (fabs(pfjets_p4().at(ijet).eta()) > 2.4) continue;
	// check for gen match? - no requirement for now
	// require match within dR < 0.4 to both genjet and status 3 quark/glu
	bool match_genjet = false;
	bool match_parton = false;
	for (unsigned int igen = 0; igen < genps_p4().size(); ++igen) {
	  if (genps_status().at(igen) != 3) continue;
	  if ((abs(genps_id().at(igen)) > 5) && (abs(genps_id().at(igen)) != 21)) continue;
	  if (dRbetweenVectors(genps_p4().at(igen),pfjets_p4().at(ijet)) < 0.4) {
	    match_parton = true;
	    break;
	  }
	}
	for (unsigned int igen = 0; igen < genjets_p4().size(); ++igen) {
	  if (dRbetweenVectors(genjets_p4().at(igen),pfjets_p4().at(ijet)) < 0.4) {
	    match_genjet = true;
	    break;
	  }
	}
	if (!match_genjet || !match_parton) continue;

	// plot beta variable: cuts used in ewkino etc
	float beta = pfjet_beta(ijet,2,0.5);
	h_pfjet_beta->Fill(beta);
	h_pfjet_beta_vs_vtx0_purity_dz->Fill(purity_dz,beta);
	h_pfjet_beta_vs_nvtx->Fill(nvtx,beta);
      }

      //      outTree->Fill();
    
    } // entries

    delete f;
  } // currentFile

  if( nSkip_els_conv_dist > 0 )
    cout << "Skipped " << nSkip_els_conv_dist << " events due to nan in els_conv_dist" << endl;

  cout << endl;
  cout << "Sample: " << prefix << endl;
  cout << endl;
  cout << "Processed events: " << nEventsTotal << endl;
  cout << "Events before reco cuts: " << nEventsPreReco << endl;
  cout << "Passed events: " << nEventsPass << endl;
  cout << endl;

  closeOutput();
  //  if(g_createTree) closeTree();
  
  already_seen.clear();

  if (nEventsChain != nEventsTotal) 
    std::cout << "ERROR: number of events from files (" << nEventsChain 
	      << ") is not equal to total number of processed events (" << nEventsTotal << ")" << std::endl;
  
  return 0;

}


//--------------------------------------------------------------------
 
void vertexStudyLooper::BookHistos(const TString& prefix)
{
  // Prefix comes from the sample and it is passed to the scanning function
  // Suffix is "ee" "em" "em" "all" which depends on the final state
  // For example: histogram named tt_hnJet_ee would be the Njet distribution
  // for the ee final state in the ttbar sample.
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  cout << "Begin book histos..." << endl;

  // TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  // rootdir->cd();
  if (outFile) outFile->cd();

  const int max_ntracks = 5000;
  const int max_nvtx = 80;

  h_vtx0_genvtx_dz = new TH1F(Form("%s_vtx0_genvtx_dz",prefix.Data()),";dz (gen vtx, reco vtx0) [cm]",10000,-10.,10.);
  h_vtx0_hardscatter_pt_vs_sumpt = new TH2F(Form("%s_vtx0_hardscatter_pt_vs_sumpt",prefix.Data()),";vtx0 #Sigma p_{T};vtx0 #Sigma p_{T} from hard scatter",100,0,1000.,100,0,1000.);
  h_vtx0_hardscatter_pt_vs_sumpt_recalc = new TH2F(Form("%s_vtx0_hardscatter_pt_vs_sumpt_recalc",prefix.Data()),";vtx0 #Sigma p_{T};vtx0 #Sigma p_{T} from hard scatter",100,0,1000.,100,0,1000.);
  h_vtx0_sumpt_vs_sumpt_recalc = new TH2F(Form("%s_vtx0_sumpt_vs_sumpt_recalc",prefix.Data()),";vtx0 #Sigma p_{T}, recalc;vtx0 #Sigma p_{T}",100,0,1000.,100,0,1000.);
  h_vtx_best_purity_dz = new TH1F(Form("%s_vtx_best_purity_dz",prefix.Data()),";Frac of vtx0 p_{T} from hard scatter",100,0.,1.);
  h_vtx_best_purity_weight = new TH1F(Form("%s_vtx_best_purity_weight",prefix.Data()),";Frac of vtx0 p_{T} from hard scatter",100,0.,1.);
  h_vtx_best_purity_dz_vs_nvtx = new TH2F(Form("%s_vtx_best_purity_dz_vs_nvtx",prefix.Data()),";N(vtx);Frac of vtx0 p_{T} from hard scatter",max_nvtx,0,max_nvtx,100,0.,1.);
  h_vtx_best_purity_weight_vs_nvtx = new TH2F(Form("%s_vtx_best_purity_weight_vs_nvtx",prefix.Data()),";N(vtx);Frac of vtx0 p_{T} from hard scatter",max_nvtx,0,max_nvtx,100,0.,1.);
  h_vtx0_purity_dz = new TH1F(Form("%s_vtx0_purity_dz",prefix.Data()),";vtx0 hard scatter purity",100,0.,1.);
  h_vtx0_purity_weight = new TH1F(Form("%s_vtx0_purity_weight",prefix.Data()),";Frac of vtx0 p_{T} from hard scatter",100,0.,1.);
  h_vtx0_purity_dz_vs_nvtx = new TH2F(Form("%s_vtx0_purity_dz_vs_nvtx",prefix.Data()),";N(vtx);Frac of vtx0 p_{T} from hard scatter",max_nvtx,0,max_nvtx,100,0.,1.);
  h_vtx0_purity_weight_vs_nvtx = new TH2F(Form("%s_vtx0_purity_weight_vs_nvtx",prefix.Data()),";N(vtx);Frac of vtx0 p_{T} from hard scatter",max_nvtx,0,max_nvtx,100,0.,1.);
  // fraction of total hard scatter pt associated to PV0
  h_vtx_best_eff_dz = new TH1F(Form("%s_vtx_best_eff_dz",prefix.Data()),";Frac of track hard scatter p_{T} assoc to vtx0",100,0.,1.);
  h_vtx_best_eff_weight = new TH1F(Form("%s_vtx_best_eff_weight",prefix.Data()),";Frac of track hard scatter p_{T} assoc to vtx0",100,0.,1.);
  h_vtx_best_eff_dz_vs_nvtx = new TH2F(Form("%s_vtx_best_eff_dz_vs_nvtx",prefix.Data()),";N(vtx);Frac of track hard scatter p_{T} assoc to vtx0",max_nvtx,0,max_nvtx,100,0.,1.);
  h_vtx_best_eff_weight_vs_nvtx = new TH2F(Form("%s_vtx_best_eff_weight_vs_nvtx",prefix.Data()),";N(vtx);Frac of track hard scatter p_{T} assoc to vtx0",max_nvtx,0,max_nvtx,100,0.,1.);
  h_vtx0_eff_dz = new TH1F(Form("%s_vtx0_eff_dz",prefix.Data()),";Frac of track hard scatter p_{T} assoc to vtx0",100,0.,1.);
  h_vtx0_eff_weight = new TH1F(Form("%s_vtx0_eff_weight",prefix.Data()),";Frac of track hard scatter p_{T} assoc to vtx0",100,0.,1.);
  h_vtx0_eff_dz_vs_nvtx = new TH2F(Form("%s_vtx0_eff_dz_vs_nvtx",prefix.Data()),";N(vtx);Frac of track hard scatter p_{T} assoc to vtx0",max_nvtx,0,max_nvtx,100,0.,1.);
  h_vtx0_eff_weight_vs_nvtx = new TH2F(Form("%s_vtx0_eff_weight_vs_nvtx",prefix.Data()),";N(vtx);Frac of track hard scatter p_{T} assoc to vtx0",max_nvtx,0,max_nvtx,100,0.,1.);
  h_dz_trk_vtx = new TH1F(Form("%s_dz_trk_vtx",prefix.Data()),";dz(trk, best vtx) [cm]",1000,-10.,10.);
  h_dz_trk_vtx0_weight = new TH1F(Form("%s_dz_trk_vtx0_weight",prefix.Data()),";dz(trk, vtx0) [cm]",1000,-10.,10.);
  h_trk_bestdzvtx = new TH1F(Form("%s_trk_bestdzvtx",prefix.Data()),";best vtx based on dz",101,-1,100);
  h_trk_bestdzvtx_vs_weightvtx = new TH2F(Form("%s_trk_bestdzvtx_vs_weightvtx",prefix.Data()),";best weighted vtx;best vtx based on dz",100,0,100,101,-1,100);
  h_trk_dz_vtxs = new TH1F(Form("%s_trk_dz_vtxs",prefix.Data()),";dz(match vtx, closest vertex) [cm]",1000,-10.,10.);
  h_hs_trk_dz_vtxs = new TH1F(Form("%s_hs_trk_dz_vtxs",prefix.Data()),";dz(match vtx,vtx0) [cm]",1000,-10.,10.);
  h_hs_trk_dz_vtx0 = new TH1F(Form("%s_hs_trk_dz_vtx0",prefix.Data()),";dz(trk,vtx0) [cm]",1000,-10.,10.);
  h_mc_idx_duplicates = new TH1F(Form("%s_mc_idx_duplicates",prefix.Data()),";# tracks with duplicate MC matches",500,0,500);
  h_mc_idx_duplicates_dr = new TH1F(Form("%s_mc_idx_duplicates_dr",prefix.Data()),";dR(reco,gen) for tracks with duplicate matches",40,0.,0.2);
  h_mc_idx_duplicates_pt = new TH1F(Form("%s_mc_idx_duplicates_pt",prefix.Data()),";p_{T} for tracks with duplicate matches",500,0.,100.);
  h_mc_idx_duplicates_nhits = new TH1F(Form("%s_mc_idx_duplicates_nhits",prefix.Data()),";N(hits) for tracks with duplicate matches",20,0,20);
  h_match_dr = new TH1F(Form("%s_match_dr",prefix.Data()),";dR(reco,gen)",40,0.,0.2);

  h_nvtx = new TH1F(Form("%s_nvtx",prefix.Data()),";N(vtx)",max_nvtx,0,max_nvtx);
  h_nvtx_vs_ntrueint = new TH2F(Form("%s_nvtx_vs_ntrueint",prefix.Data()),";N(true int);N(vtx)",max_nvtx,0,max_nvtx,max_nvtx,0,max_nvtx);
  h_ntracks = new TH1F(Form("%s_ntracks",prefix.Data()),";N(tracks)",max_ntracks/10,0,max_ntracks);

  h_gen_match_vtx = new TH1F(Form("%s_gen_match_vtx",prefix.Data()),";reco index of true PV",51,-1,50);
  h_gen_match_vtx_vs_nvtx = new TH2F(Form("%s_gen_match_vtx_vs_nvtx",prefix.Data()),";N(vtx);reco index of true PV",max_nvtx,0,max_nvtx,51,-1,50);

  h_genvtx_z_nocut = new TH1F(Form("%s_genvtx_z_nocut",prefix.Data()),";Z (gen vtx) [cm]",50,-25.,25.);
  h_genvtx_z = new TH1F(Form("%s_genvtx_z",prefix.Data()),";Z (gen vtx) [cm]",50,-25.,25.);
  h_genmatch_vtx_sumpt2 = new TH1F(Form("%s_genmatch_vtx_sumpt2",prefix.Data()),";Vertex #Sigma gen p_{T}^{2} [GeV^{2}]",2500,0.,5000.);
  h_genmatch_vtx_z = new TH1F(Form("%s_genmatch_vtx_z",prefix.Data()),";Z (gen match vtx) [cm]",50,-25.,25.);
  h_genmatch_leptype = new TH1F(Form("%s_genmatch_leptype",prefix.Data()),";lep type",6,0.,6.);
  h_recovtx_z = new TH1F(Form("%s_recovtx_z",prefix.Data()),";Z (reco vtx) [cm]",50,-25.,25.);
  h_genvtx_reco_dz = new TH1F(Form("%s_genvtx_reco_dz",prefix.Data()),";dz (gen vtx, closest reco) [cm]",1000,-10.,10.);

  h_genmatch_lowsumpt2_leptype = new TH1F(Form("%s_genmatch_lowsumpt2_leptype",prefix.Data()),";lep type",6,0.,6.);
  h_genmatch_lowsumpt2_eltrkassoc = new TH1F(Form("%s_genmatch_lowsumpt2_eltrkassoc",prefix.Data()),"; el track assoc with vtx",4,-2.,2.);

  h_genvtx_nomatch_z = new TH1F(Form("%s_genvtx_nomatch_z",prefix.Data()),";Z (gen vtx) [cm]",50,-25.,25.);
  h_genvtx_nomatch_gensumpt2 = new TH1F(Form("%s_genvtx_nomatch_gensumpt2",prefix.Data()),";Vertex #Sigma p_{T}^{2} [GeV^{2}]",2500,0.,5000.);
  h_genvtx_nomatch_dz = new TH1F(Form("%s_genvtx_nomatch_dz",prefix.Data()),";dz (gen vtx, closest reco) [cm]",1000,-10.,10.);
  h_genvtx_nomatch_dz_eff10 = new TH1F(Form("%s_genvtx_nomatch_dz_eff10",prefix.Data()),";dz (gen vtx, closest reco) [cm]",1000,-10.,10.);
  h_genvtx_nomatch_smalldz_eff = new TH1F(Form("%s_genvtx_nomatch_smalldz_eff",prefix.Data()),";Frac of track hard scatter p_{T} assoc to vtx",100,0.,1.);
  h_genvtx_nomatch_largedz_nleps = new TH1F(Form("%s_genvtx_nomatch_largedz_nleps",prefix.Data()),";N(gen leps)",5,0.,5.);
  h_genvtx_nomatch_largedz_npartons = new TH1F(Form("%s_genvtx_nomatch_largedz_npartons",prefix.Data()),";N(gen partons)",10,0.,10.);
  h_genvtx_nomatch_nvtx = new TH1F(Form("%s_genvtx_nomatch_nvtx",prefix.Data()),";N(vtx)",max_nvtx,0,max_nvtx);
  h_genvtx_nomatch_ntracks = new TH1F(Form("%s_genvtx_nomatch_ntracks",prefix.Data()),";N(tracks)",max_ntracks/10,0,max_ntracks);

  h_genvtx_othermatch_z = new TH1F(Form("%s_genvtx_othermatch_z",prefix.Data()),";Z (gen vtx) [cm]",50,-25.,25.);
  h_genvtx_othermatch_dz = new TH1F(Form("%s_genvtx_othermatch_dz",prefix.Data()),";dz (gen vtx, highest #Sigma p_{T}^{2} vtx) [cm]",1000,-10.,10.);
  h_genvtx_othermatch_gensumpt2 = new TH1F(Form("%s_genvtx_othermatch_gensumpt2",prefix.Data()),";Vertex #Sigma gen p_{T}^{2} [GeV^{2}]",2500,0.,5000.);
  h_genvtx_othermatch_sumpt2 = new TH1F(Form("%s_genvtx_othermatch_sumpt2",prefix.Data()),";Vertex #Sigma p_{T}^{2} [GeV^{2}]",2500,0.,5000.);
  h_genvtx_othermatch_leptype = new TH1F(Form("%s_genvtx_othermatch_leptype",prefix.Data()),";lep type",6,0.,6.);

  h_vtx_sumpt2 = new TH1F(Form("%s_vtx_sumpt2",prefix.Data()),";Vertex #Sigma p_{T}^{2} [GeV^{2}]",2500,0.,5000.);
  h_vtx_nohs_sumpt2 = new TH1F(Form("%s_vtx_nohs_sumpt2",prefix.Data()),";Vertex #Sigma p_{T}^{2} [GeV^{2}]",2500,0.,5000.);
  h_vtx_nogen_sumpt2 = new TH1F(Form("%s_vtx_nogen_sumpt2",prefix.Data()),";Vertex #Sigma p_{T}^{2} [GeV^{2}]",2500,0.,5000.);
  h_vtx_nogen_nohs_sumpt2 = new TH1F(Form("%s_vtx_nogen_nohs_sumpt2",prefix.Data()),";Vertex #Sigma p_{T}^{2} [GeV^{2}]",2500,0.,5000.);
  h_vtx0_pu_sumpt2 = new TH1F(Form("%s_vtx0_pu_sumpt2",prefix.Data()),";Vertex #Sigma p_{T}^{2} [GeV^{2}]",2500,0.,5000.);

  // lepton iso, vs purity
  h_el_iso = new TH1F(Form("%s_el_iso",prefix.Data()),";electron reliso, no PU cor",200,0.,2.);
  h_el_iso_cor = new TH1F(Form("%s_el_iso_cor",prefix.Data()),";electron reliso, with PU cor",200,0.,2.);
  h_el_iso_vs_vtx0_purity_dz = new TH2F(Form("%s_el_iso_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;electron reliso, no PU cor",100,0,1.,200,0.,2.);
  h_el_iso_cor_vs_vtx0_purity_dz = new TH2F(Form("%s_el_iso_cor_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;electron reliso, with PU cor",100,0,1.,200,0.,2.);
  h_el_iso_vs_nvtx = new TH2F(Form("%s_el_iso_vs_nvtx",prefix.Data()),";nvtx;electron reliso, no PU cor",max_nvtx,0,max_nvtx,200,0.,2.);
  h_el_iso_cor_vs_nvtx = new TH2F(Form("%s_el_iso_cor_vs_nvtx",prefix.Data()),";nvtx;electron reliso, with PU cor",max_nvtx,0,max_nvtx,200,0.,2.);
  h_el_iso40_cor_vs_nvtx = new TH2F(Form("%s_el_iso40_cor_vs_nvtx",prefix.Data()),";nvtx;electron reliso, with PU cor",max_nvtx,0,max_nvtx,200,0.,2.);

  h_el_trkiso = new TH1F(Form("%s_el_trkiso",prefix.Data()),";el rel trkiso",1000,0.,2.);
  h_el_trkiso_vs_vtx0_purity_dz = new TH2F(Form("%s_el_trkiso_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;el rel trkiso",100,0,1.,1000,0.,2.);
  h_el_trkiso_vs_nvtx = new TH2F(Form("%s_el_trkiso_vs_nvtx",prefix.Data()),";nvtx;el rel trkiso",max_nvtx,0,max_nvtx,1000,0.,2.);
  h_el_trkiso40_vs_nvtx = new TH2F(Form("%s_el_trkiso40_vs_nvtx",prefix.Data()),";nvtx;el rel trkiso",max_nvtx,0,max_nvtx,1000,0.,2.);
  h_el_trkiso_abs = new TH1F(Form("%s_el_trkiso_abs",prefix.Data()),";el trkiso",1000,0.,10.);
  h_el_trkiso_abs_vs_vtx0_purity_dz = new TH2F(Form("%s_el_trkiso_abs_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;el trkiso",100,0,1.,1000,0.,10.);
  h_el_trkiso_abs_vs_nvtx = new TH2F(Form("%s_el_trkiso_abs_vs_nvtx",prefix.Data()),";nvtx;el trkiso",max_nvtx,0,max_nvtx,1000,0.,10.);
  h_el_trkiso40_abs_vs_nvtx = new TH2F(Form("%s_el_trkiso40_abs_vs_nvtx",prefix.Data()),";nvtx;el trkiso",max_nvtx,0,max_nvtx,1000,0.,10.);

  h_mu_iso = new TH1F(Form("%s_mu_iso",prefix.Data()),";mu reliso, no PU cor",200,0.,2.);
  h_mu_iso_cor = new TH1F(Form("%s_mu_iso_cor",prefix.Data()),";mu reliso, with PU cor",200,0.,2.);
  h_mu_iso_vs_vtx0_purity_dz = new TH2F(Form("%s_mu_iso_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;mu reliso, no PU cor",100,0,1.,200,0.,2.);
  h_mu_iso_cor_vs_vtx0_purity_dz = new TH2F(Form("%s_mu_iso_cor_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;mu reliso, with PU cor",100,0,1.,200,0.,2.);
  h_mu_iso_vs_nvtx = new TH2F(Form("%s_mu_iso_vs_nvtx",prefix.Data()),";nvtx;mu reliso, no PU cor",max_nvtx,0,max_nvtx,200,0.,2.);
  h_mu_iso_cor_vs_nvtx = new TH2F(Form("%s_mu_iso_cor_vs_nvtx",prefix.Data()),";nvtx;mu reliso, with PU cor",max_nvtx,0,max_nvtx,200,0.,2.);
  h_mu_iso40_cor_vs_nvtx = new TH2F(Form("%s_mu_iso40_cor_vs_nvtx",prefix.Data()),";nvtx;mu reliso, with PU cor",max_nvtx,0,max_nvtx,200,0.,2.);

  h_mu_trkiso = new TH1F(Form("%s_mu_trkiso",prefix.Data()),";mu rel trkiso",1000,0.,2.);
  h_mu_trkiso_vs_vtx0_purity_dz = new TH2F(Form("%s_mu_trkiso_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;mu rel trkiso",100,0,1.,1000,0.,2.);
  h_mu_trkiso_vs_nvtx = new TH2F(Form("%s_mu_trkiso_vs_nvtx",prefix.Data()),";nvtx;mu rel trkiso",max_nvtx,0,max_nvtx,1000,0.,2.);
  h_mu_trkiso40_vs_nvtx = new TH2F(Form("%s_mu_trkiso40_vs_nvtx",prefix.Data()),";nvtx;mu rel trkiso",max_nvtx,0,max_nvtx,1000,0.,2.);
  h_mu_trkiso_abs = new TH1F(Form("%s_mu_trkiso_abs",prefix.Data()),";mu trkiso",1000,0.,10.);
  h_mu_trkiso_abs_vs_vtx0_purity_dz = new TH2F(Form("%s_mu_trkiso_abs_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;mu trkiso",100,0,1.,1000,0.,10.);
  h_mu_trkiso_abs_vs_nvtx = new TH2F(Form("%s_mu_trkiso_abs_vs_nvtx",prefix.Data()),";nvtx;mu trkiso",max_nvtx,0,max_nvtx,1000,0.,10.);
  h_mu_trkiso40_abs_vs_nvtx = new TH2F(Form("%s_mu_trkiso40_abs_vs_nvtx",prefix.Data()),";nvtx;mu trkiso",max_nvtx,0,max_nvtx,1000,0.,10.);

  h_ph_trkiso = new TH1F(Form("%s_ph_trkiso",prefix.Data()),";ph trkiso",1000,0.,2.);
  h_ph_trkiso_vs_vtx0_purity_dz = new TH2F(Form("%s_ph_trkiso_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;ph trkiso",100,0,1.,1000,0.,2.);
  h_ph_trkiso_vs_nvtx = new TH2F(Form("%s_ph_trkiso_vs_nvtx",prefix.Data()),";nvtx;ph trkiso",max_nvtx,0,max_nvtx,1000,0.,2.);

  h_pfjet_beta = new TH1F(Form("%s_pfjet_beta",prefix.Data()),";pfjet beta",100,0.,1.);
  h_pfjet_beta_vs_vtx0_purity_dz = new TH2F(Form("%s_pfjet_beta_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;pfjet beta",100,0,1.,100,0.,1.);
  h_pfjet_beta_vs_nvtx = new TH2F(Form("%s_pfjet_beta_vs_nvtx",prefix.Data()),";nvtx;pfjet beta",max_nvtx,0,max_nvtx,100,0.,1.);

  h_match_trk_pt = new TH1F(Form("%s_match_trk_pt",prefix.Data()),";p_{T} (tracks, matched) [GeV]",500,0.,100.);
  h_match_trk_pt_low = new TH1F(Form("%s_match_trk_pt_low",prefix.Data()),";p_{T} (tracks, matched) [GeV]",100,0.,1.);
  h_match_trk_pt_mid = new TH1F(Form("%s_match_trk_pt_mid",prefix.Data()),";p_{T} (tracks, matched) [GeV]",100,1.,10.);
 
  h_nomatch_trk_pt = new TH1F(Form("%s_nomatch_trk_pt",prefix.Data()),";p_{T} (tracks, not matched) [GeV]",500,0.,100.);
  h_nomatch_trk_pt_low = new TH1F(Form("%s_nomatch_trk_pt_low",prefix.Data()),";p_{T} (tracks, not matched) [GeV]",100,0.,1.);
  h_nomatch_trk_pt_mid = new TH1F(Form("%s_nomatch_trk_pt_mid",prefix.Data()),";p_{T} (tracks, not matched) [GeV]",100,1.,10.);
 
  h_nomatch_trk_pt_vs_eta = new TH2F(Form("%s_nomatch_trk_pt_vs_eta",prefix.Data()),";#eta (tracks, not matched);p_{T} (tracks, not matched) [GeV]",50,-2.5,2.5,500,0.,100.);
  h_nomatch_trk_pt_low_vs_eta = new TH2F(Form("%s_nomatch_trk_pt_low_vs_eta",prefix.Data()),";#eta (tracks, not matched);p_{T} (tracks, not matched) [GeV]",50,-2.5,2.5,100,0.,1.);

  cout << "End book histos..." << endl;
}// CMS2::BookHistos()


void vertexStudyLooper::makeOutput(const TString& prefix){
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  TString tpsuffix = "";
  if( doTenPercent ) tpsuffix = "_tenPercent";

  outFile   = new TFile(Form("output/%s_smallTree%s.root",prefix.Data(),tpsuffix.Data()), "RECREATE");
  //  outFile   = new TFile(Form("output/%s/%s_smallTree%s%s.root",g_version,prefix.Data(),frsuffix,tpsuffix), "RECREATE");
  //  outFile   = new TFile("baby.root","RECREATE");
  outFile->cd();
  //  outTree = new TTree("t","Tree");

  //Set branch addresses
  //variables must be declared in vertexStudyLooper.h
}

//--------------------------------------------------------------------

int vertexStudyLooper::associateTrackToVertexSimple(const unsigned int itrk) {
    
  int best_vtx = -1;
  float dz = -999.;
  for (unsigned int ivtx = 0; ivtx < cms2.vtxs_position().size(); ivtx++) {
    float tmp_dz = cms2.trks_z0().at(itrk) - vtxs_position().at(ivtx).z();
    if (fabs(tmp_dz) > fabs(dz))
      continue;

    dz = tmp_dz;
    best_vtx = ivtx;
  }

  return best_vtx;
}

//--------------------------------------------------------------------

float vertexStudyLooper::electronPFiso(const unsigned int index, const bool cor) {
    
  float pt     = cms2.els_p4().at(index).pt();
  float etaAbs = fabs(cms2.els_etaSC().at(index));

  // get effective area
  float AEff = fastJetEffArea03_v1(etaAbs);
  // float AEff = 0.;
  // if (etaAbs <= 1.0) AEff = 0.10;
  // else if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.12;
  // else if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.085;
  // else if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.11;
  // else if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.12;
  // else if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.12;
  // else if (etaAbs > 2.4) AEff = 0.13;

  float pfiso_ch = cms2.els_iso03_pf2012ext_ch().at(index);
  float pfiso_em = cms2.els_iso03_pf2012ext_em().at(index);
  float pfiso_nh = cms2.els_iso03_pf2012ext_nh().at(index);
    
  // rho
  float rhoPrime = std::max(cms2.evt_kt6pf_foregiso_rho(), float(0.0));
  float pfiso_n = pfiso_em + pfiso_nh;
  if (cor) pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, float(0.0));
  float pfiso = (pfiso_ch + pfiso_n) / pt;

  return pfiso;
}

//--------------------------------------------------------------------

float vertexStudyLooper::electronPFiso40(const unsigned int index, const bool cor) {
    
  float pt     = cms2.els_p4().at(index).pt();
  float etaAbs = fabs(cms2.els_etaSC().at(index));

  // get effective area
  float AEff = fastJetEffArea04_v1(etaAbs);
  // float AEff = 0.;
  // if (etaAbs <= 1.0) AEff = 0.10;
  // else if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.12;
  // else if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.085;
  // else if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.11;
  // else if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.12;
  // else if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.12;
  // else if (etaAbs > 2.4) AEff = 0.13;

  float pfiso_ch = cms2.els_iso04_pf2012ext_ch().at(index);
  float pfiso_em = cms2.els_iso04_pf2012ext_em().at(index);
  float pfiso_nh = cms2.els_iso04_pf2012ext_nh().at(index);
    
  // rho
  float rhoPrime = std::max(cms2.evt_kt6pf_foregiso_rho(), float(0.0));
  float pfiso_n = pfiso_em + pfiso_nh;
  if (cor) pfiso_n = std::max(pfiso_em + pfiso_nh - rhoPrime * AEff, float(0.0));
  float pfiso = (pfiso_ch + pfiso_n) / pt;

  return pfiso;
}

//--------------------------------------------------------------------

float vertexStudyLooper::muonPFiso(const unsigned int imu, const bool cor) {
  float chiso = cms2.mus_isoR03_pf_ChargedHadronPt().at(imu);
  float nhiso = cms2.mus_isoR03_pf_NeutralHadronEt().at(imu);
  float emiso = cms2.mus_isoR03_pf_PhotonEt().at(imu);
  float deltaBeta = cms2.mus_isoR03_pf_PUPt().at(imu);
  float pt = cms2.mus_p4().at(imu).pt();

  float absiso = chiso + nhiso + emiso;
  if (cor) absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return (absiso / pt);

}

//--------------------------------------------------------------------

float vertexStudyLooper::muonPFiso40(const unsigned int imu, const bool cor) {
  float chiso = cms2.mus_isoR04_pf_ChargedHadronPt().at(imu);
  float nhiso = cms2.mus_isoR04_pf_NeutralHadronEt().at(imu);
  float emiso = cms2.mus_isoR04_pf_PhotonEt().at(imu);
  float deltaBeta = cms2.mus_isoR04_pf_PUPt().at(imu);
  float pt = cms2.mus_p4().at(imu).pt();

  float absiso = chiso + nhiso + emiso;
  if (cor) absiso = chiso + std::max(0.0, nhiso + emiso - 0.5 * deltaBeta);
  return (absiso / pt);

}

//--------------------------------------------------------------------

float vertexStudyLooper::photonHollowTrkIso(const unsigned int iph) {

  float dR_outer = 0.3;
  float dR_inner = 0.05;
  float iso_sum = 0.;

  for( unsigned int itrk = 0; itrk < cms2.trks_trk_p4().size(); ++itrk ){
    float dR = dRbetweenVectors( cms2.photons_p4().at(iph), cms2.trks_trk_p4().at(itrk) );
    if( dR < dR_outer && dR > dR_inner ){
      iso_sum += cms2.trks_trk_p4().at(itrk).pt();         
    }
  }

  return iso_sum;
}

//--------------------------------------------------------------------

float vertexStudyLooper::genSumPt2() {

  float sumpt2 = 0.;

  for (unsigned int igen = 0; igen < genps_p4().size(); ++igen) {
    if (genps_status().at(igen) != 3) continue;

    // skip lines up to t and tbar -- very ttbar specific?..
    if( igen < 8 ) continue;

    int id = abs(genps_id().at(igen));
    // only sum up: (udscb quarks, gluons) * 2/3?, charged leptons
    if ((id < 6) || (id == 21)) {
      // scale by 2/3 to account for energy in neutrals?
      sumpt2 += pow(genps_p4().at(igen).pt() * 2./3.,2);
    } else if ((id == 11) || (id == 13) || (id == 15)) {
      sumpt2 += pow(genps_p4().at(igen).pt(),2);
    }
  }

  return sumpt2;
}

