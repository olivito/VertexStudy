#include "vertexStudyLooper.h"
#include "TTreeCache.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include "../Tools/goodrun.h"
#include "../Tools/vtxreweight.h"
//#include "../Tools/pfjetMVAtools.h"

#include "../CORE/eventSelections.h"
#include "../CORE/trackSelections.h"
#include "../CORE/jetSelections.h"
#include "../CORE/metSelections.h"
#include "../CORE/jetSmearingTools.h"
#include "../CORE/jetcorr/JetCorrectionUncertainty.h"

bool verbose              = false;
bool doTenPercent         = false;
bool requireLeps          = false;

const float dzcut = 1.;
const float drcut = 0.015;
const float dz_gen_vtx_match = 1.;
const float fracpt_gen_vtx_match = 0.1;

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
   
      int nleps = 0;
      // loop over genps, count e and mu with pt > 15 gev
      for (unsigned int igen = 0; igen < genps_status().size(); ++igen) {
	if (genps_status().at(igen) != 3) continue;
	if (genps_p4().at(igen).pt() < 15.) continue;
	int id = abs(genps_id().at(igen));
	if (id == 11 || id == 13) ++nleps;
      }

      if (requireLeps && !nleps) continue;

      ++nEventsPass;

      //---------------------------------------------
      // loop over tracks, associate to vtx
      //---------------------------------------------

      std::vector<int> trk_bestdzvtx(trks_trk_p4().size(), -1);
      std::vector<float> vtxs_sumpt_recalc_dz(vtxs_sumpt().size(), 0.);
      std::vector<float> vtxs_sumpt_recalc_weight(vtxs_sumpt().size(), 0.);
      std::vector<float> vtxs_sumpt_hardscatter_dz(vtxs_sumpt().size(), 0.);
      std::vector<float> vtxs_sumpt_hardscatter_weight(vtxs_sumpt().size(), 0.);
      std::vector<float> vtxs_fracpt_hardscatter(vtxs_sumpt().size(), -1.);
      std::map<int,std::vector<int> > mc_idx_reco_matches;
      int mc_idx_duplicates = 0;
      float sum_hardscatter_pt = 0.;

      for (unsigned int itrk = 0; itrk < trks_trk_p4().size(); ++itrk) {
	// only consider high purity tracks: bit 2 of quality mask
	if (!(trks_qualityMask().at(itrk) & (1<<2))) continue;

	// find best match vertex, from vertex algo weight and simple dz
	int bestdzvtx = associateTrackToVertex(itrk);
	//int bestdzvtx = associateTrackToVertexSimple(itrk);
	int weightvtx = trks_pvidx0().at(itrk);
	trk_bestdzvtx.at(itrk) = bestdzvtx;
	float dz = trks_dz_pv(itrk, bestdzvtx).first;
	//float dz = cms2.trks_z0().at(itrk) - vtxs_position().at(bestdzvtx).z();
	h_dz_trk_vtx->Fill(dz);
	if (fabs(dz) <= dzcut) {
	  h_trk_bestdzvtx->Fill(bestdzvtx);
	  h_trk_bestdzvtx_vs_weightvtx->Fill(trks_pvidx0().at(itrk),bestdzvtx);
	  if ((weightvtx > -9000) && (bestdzvtx != weightvtx)) {
	    h_trk_dz_vtxs->Fill(vtxs_position().at(weightvtx).z() - vtxs_position().at(bestdzvtx).z());
	  }
	}
	vtxs_sumpt_recalc_dz.at(bestdzvtx) += trks_trk_p4().at(itrk).pt();
	if (weightvtx > -9000) {
	  vtxs_sumpt_recalc_weight.at(weightvtx) += trks_trk_p4().at(itrk).pt();
	  if(weightvtx == 0) h_dz_trk_vtx0_weight->Fill(trks_dz_pv(itrk,0).first);
	}
	// check that track is matched to status 1 particle from hard scatter
	if (trk_mc_id().at(itrk) == -9999) continue;

	h_match_dr->Fill(trk_mcdr().at(itrk));

	// apply cut on dR of reco/gen match to remove most multiply matched tracks
	if (trk_mcdr().at(itrk) > drcut) continue;

	// check for tracks matched to the same gen particle..
	if (mc_idx_reco_matches.count(trk_mcidx().at(itrk)) > 0) {
	  ++mc_idx_duplicates;
	  mc_idx_reco_matches[trk_mcidx().at(itrk)].push_back(itrk);
	}
	else {
	  std::vector<int> reco_matches(1,itrk);
	  mc_idx_reco_matches[trk_mcidx().at(itrk)] = reco_matches;
	}

	// vtxs_sumpt_hardscatter_dz.at(bestdzvtx) += trks_trk_p4().at(itrk).pt();
	// if (weightvtx > -9000) {
	//   vtxs_sumpt_hardscatter_weight.at(weightvtx) += trks_trk_p4().at(itrk).pt();
	//   if (weightvtx != 0) {
	//     h_hs_trk_dz_vtxs->Fill(vtxs_position().at(weightvtx).z() - vtxs_position().at(0).z());
	//     h_hs_trk_dz_vtx0->Fill(trks_dz_pv(itrk, 0).first);
	//   }
	// }
	// sum_hardscatter_pt += trks_trk_p4().at(itrk).pt();

      }

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

	  vtxs_sumpt_hardscatter_dz.at(bestdzvtx) += trks_trk_p4().at(itrk).pt();
	  if (weightvtx > -9000) {
	    vtxs_sumpt_hardscatter_weight.at(weightvtx) += trks_trk_p4().at(itrk).pt();
	    if (weightvtx != 0) {
	      h_hs_trk_dz_vtxs->Fill(vtxs_position().at(weightvtx).z() - vtxs_position().at(0).z());
	      h_hs_trk_dz_vtx0->Fill(trks_dz_pv(itrk, 0).first);
	    }
	  }
	  sum_hardscatter_pt += trks_trk_p4().at(itrk).pt();
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

	  vtxs_sumpt_hardscatter_dz.at(bestdzvtx) += trks_trk_p4().at(ibestmatch).pt();
	  if (weightvtx > -9000) {
	    vtxs_sumpt_hardscatter_weight.at(weightvtx) += trks_trk_p4().at(ibestmatch).pt();
	    if (weightvtx != 0) {
	      h_hs_trk_dz_vtxs->Fill(vtxs_position().at(weightvtx).z() - vtxs_position().at(0).z());
	      h_hs_trk_dz_vtx0->Fill(trks_dz_pv(ibestmatch, 0).first);
	    }
	  }
	  sum_hardscatter_pt += trks_trk_p4().at(ibestmatch).pt();
	} // if duplicate matches

      } // loop over matched gen particles

      // for (unsigned int itrk = 0; itrk < trks_trk_p4().size(); ++itrk) {
      // 	if (trk_mc_id().at(itrk) == -9999) continue;
      // 	if (mc_idxs[trk_mcidx().at(itrk)] < 2) continue;
      // 	// only consider high purity tracks: bit 2 of quality mask
      // 	if (!(trks_qualityMask().at(itrk) & (1<<2))) continue;
      // 	// apply cut on dR of reco/gen match to remove multiply matched tracks
      // 	if (trk_mcdr().at(itrk) > drcut) continue;

      // 	// make plots only for surviving tracks
      // 	h_mc_idx_duplicates_dr->Fill(trk_mcdr().at(itrk));
      // 	h_mc_idx_duplicates_pt->Fill(trks_trk_p4().at(itrk).pt());
      // 	h_mc_idx_duplicates_nhits->Fill(trks_validHits().at(itrk));
      // }

      int nvtx = 0;
      int gen_match_vtx = -1;
      for (unsigned int ivtx = 0; ivtx < vtxs_sumpt().size(); ++ivtx) {
	// count good vertices
	if (!isGoodVertex(ivtx)) continue;
	++nvtx;

	// check to see if hard scatter vertex is among the reco vtx collection
	// if already matched to a higher sum pt^2 vertex, don't bother
	if (gen_match_vtx > -1) continue;
	// require: 
	//  |dz| < 1mm from true hard scatter vtx
	//  at least 10% of hard scatter pt
	if (fabs(vtxs_position().at(ivtx).z() - genps_prod_vtx().at(2).z()) > dz_gen_vtx_match ) continue;
	//	if (vtxs_sumpt_hardscatter_dz.at(ivtx)/vtxs_sumpt_recalc_dz.at(ivtx) < fracpt_gen_vtx_match ) continue;

	float purity_dz = vtxs_sumpt_hardscatter_dz.at(ivtx)/vtxs_sumpt_recalc_dz.at(ivtx);
	float purity_weight = vtxs_sumpt_hardscatter_weight.at(ivtx)/vtxs_sumpt_recalc_weight.at(ivtx);
	float eff_dz = vtxs_sumpt_hardscatter_dz.at(ivtx)/sum_hardscatter_pt;
	float eff_weight = vtxs_sumpt_hardscatter_weight.at(ivtx)/sum_hardscatter_pt;

	if (eff_dz < fracpt_gen_vtx_match) continue;

	h_vtx_best_purity_dz->Fill(purity_dz);
	h_vtx_best_purity_dz_vs_nvtx->Fill(nvtx,purity_dz);
	h_vtx_best_purity_weight->Fill(purity_weight);
	h_vtx_best_purity_weight_vs_nvtx->Fill(nvtx,purity_weight);

	h_vtx_best_eff_dz->Fill(eff_dz);
	h_vtx_best_eff_dz_vs_nvtx->Fill(nvtx,eff_dz);
	h_vtx_best_eff_weight->Fill(eff_weight);
	h_vtx_best_eff_weight_vs_nvtx->Fill(nvtx,eff_weight);

	gen_match_vtx = ivtx;
      }
      h_nvtx->Fill(nvtx);
      h_gen_match_vtx->Fill(gen_match_vtx);

      // plots for vertex 0
      const int vtx0 = 0;
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

      // electron iso
      for (unsigned int iel = 0; iel < els_p4().size(); ++iel) {
	// cut on pt and eta
	if (els_p4().at(iel).pt() < 15.) continue;
	if (fabs(els_p4().at(iel).eta()) > 2.4) continue;
	// check for gen match
	if (abs(els_mc_id().at(iel)) != 11) continue;
	// plot isolation with and without pileup correction, also vs vtx0 purity
	float iso = electronPFiso(iel);
	float iso_cor = electronPFiso(iel,true);

	h_el_iso->Fill(iso);
	h_el_iso_cor->Fill(iso_cor);
	h_el_iso_vs_vtx0_purity_dz->Fill(purity_dz,iso);
	h_el_iso_cor_vs_vtx0_purity_dz->Fill(purity_dz,iso_cor);
	h_el_iso_vs_nvtx->Fill(nvtx,iso);
	h_el_iso_cor_vs_nvtx->Fill(nvtx,iso_cor);
      }

      // muon iso
      for (unsigned int imu = 0; imu < mus_p4().size(); ++imu) {
	// cut on pt and eta
	if (mus_p4().at(imu).pt() < 15.) continue;
	if (fabs(mus_p4().at(imu).eta()) > 2.4) continue;
	// check for gen match
	if (abs(mus_mc_id().at(imu)) != 13) continue;
	// plot isolation with and without pileup correction, also vs vtx0 purity
	float iso = muonPFiso(imu);
	float iso_cor = muonPFiso(imu,true);

	h_mu_iso->Fill(iso);
	h_mu_iso_cor->Fill(iso_cor);
	h_mu_iso_vs_vtx0_purity_dz->Fill(purity_dz,iso);
	h_mu_iso_cor_vs_vtx0_purity_dz->Fill(purity_dz,iso_cor);
	h_mu_iso_vs_nvtx->Fill(nvtx,iso);
	h_mu_iso_cor_vs_nvtx->Fill(nvtx,iso_cor);
      }

      // jet beta
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

  h_vtx0_hardscatter_pt_vs_sumpt = new TH2F(Form("%s_vtx0_hardscatter_pt_vs_sumpt",prefix.Data()),";vtx0 #Sigma p_{T};vtx0 #Sigma p_{T} from hard scatter",100,0,1000.,100,0,1000.);
  h_vtx0_hardscatter_pt_vs_sumpt_recalc = new TH2F(Form("%s_vtx0_hardscatter_pt_vs_sumpt_recalc",prefix.Data()),";vtx0 #Sigma p_{T};vtx0 #Sigma p_{T} from hard scatter",100,0,1000.,100,0,1000.);
  h_vtx0_sumpt_vs_sumpt_recalc = new TH2F(Form("%s_vtx0_sumpt_vs_sumpt_recalc",prefix.Data()),";vtx0 #Sigma p_{T}, recalc;vtx0 #Sigma p_{T}",100,0,1000.,100,0,1000.);
  h_vtx_best_purity_dz = new TH1F(Form("%s_vtx_best_purity_dz",prefix.Data()),";Frac of vtx0 p_{T} from hard scatter",100,0.,1.);
  h_vtx_best_purity_weight = new TH1F(Form("%s_vtx_best_purity_weight",prefix.Data()),";Frac of vtx0 p_{T} from hard scatter",100,0.,1.);
  h_vtx_best_purity_dz_vs_nvtx = new TH2F(Form("%s_vtx_best_purity_dz_vs_nvtx",prefix.Data()),";N(vtx);Frac of vtx0 p_{T} from hard scatter",50,0,50,100,0.,1.);
  h_vtx_best_purity_weight_vs_nvtx = new TH2F(Form("%s_vtx_best_purity_weight_vs_nvtx",prefix.Data()),";N(vtx);Frac of vtx0 p_{T} from hard scatter",50,0,50,100,0.,1.);
  h_vtx0_purity_dz = new TH1F(Form("%s_vtx0_purity_dz",prefix.Data()),";vtx0 hard scatter purity",100,0.,1.);
  h_vtx0_purity_weight = new TH1F(Form("%s_vtx0_purity_weight",prefix.Data()),";Frac of vtx0 p_{T} from hard scatter",100,0.,1.);
  h_vtx0_purity_dz_vs_nvtx = new TH2F(Form("%s_vtx0_purity_dz_vs_nvtx",prefix.Data()),";N(vtx);Frac of vtx0 p_{T} from hard scatter",50,0,50,100,0.,1.);
  h_vtx0_purity_weight_vs_nvtx = new TH2F(Form("%s_vtx0_purity_weight_vs_nvtx",prefix.Data()),";N(vtx);Frac of vtx0 p_{T} from hard scatter",50,0,50,100,0.,1.);
  // fraction of total hard scatter pt associated to PV0
  h_vtx_best_eff_dz = new TH1F(Form("%s_vtx_best_eff_dz",prefix.Data()),";Frac of track hard scatter p_{T} assoc to vtx0",100,0.,1.);
  h_vtx_best_eff_weight = new TH1F(Form("%s_vtx_best_eff_weight",prefix.Data()),";Frac of track hard scatter p_{T} assoc to vtx0",100,0.,1.);
  h_vtx_best_eff_dz_vs_nvtx = new TH2F(Form("%s_vtx_best_eff_dz_vs_nvtx",prefix.Data()),";N(vtx);Frac of track hard scatter p_{T} assoc to vtx0",50,0,50,100,0.,1.);
  h_vtx_best_eff_weight_vs_nvtx = new TH2F(Form("%s_vtx_best_eff_weight_vs_nvtx",prefix.Data()),";N(vtx);Frac of track hard scatter p_{T} assoc to vtx0",50,0,50,100,0.,1.);
  h_vtx0_eff_dz = new TH1F(Form("%s_vtx0_eff_dz",prefix.Data()),";Frac of track hard scatter p_{T} assoc to vtx0",100,0.,1.);
  h_vtx0_eff_weight = new TH1F(Form("%s_vtx0_eff_weight",prefix.Data()),";Frac of track hard scatter p_{T} assoc to vtx0",100,0.,1.);
  h_vtx0_eff_dz_vs_nvtx = new TH2F(Form("%s_vtx0_eff_dz_vs_nvtx",prefix.Data()),";N(vtx);Frac of track hard scatter p_{T} assoc to vtx0",50,0,50,100,0.,1.);
  h_vtx0_eff_weight_vs_nvtx = new TH2F(Form("%s_vtx0_eff_weight_vs_nvtx",prefix.Data()),";N(vtx);Frac of track hard scatter p_{T} assoc to vtx0",50,0,50,100,0.,1.);
  h_dz_trk_vtx = new TH1F(Form("%s_dz_trk_vtx",prefix.Data()),";dz(trk, best vtx) [mm]",1000,-10.,10.);
  h_dz_trk_vtx0_weight = new TH1F(Form("%s_dz_trk_vtx0_weight",prefix.Data()),";dz(trk, vtx0) [mm]",1000,-10.,10.);
  h_trk_bestdzvtx = new TH1F(Form("%s_trk_bestdzvtx",prefix.Data()),";best vtx based on dz",100,0,100);
  h_trk_bestdzvtx_vs_weightvtx = new TH2F(Form("%s_trk_bestdzvtx_vs_weightvtx",prefix.Data()),";best weighted vtx;best vtx based on dz",100,0,100,100,0,100);
  h_trk_dz_vtxs = new TH1F(Form("%s_trk_dz_vtxs",prefix.Data()),";dz(match vtx, closest vertex) [mm]",1000,-10.,10.);
  h_hs_trk_dz_vtxs = new TH1F(Form("%s_hs_trk_dz_vtxs",prefix.Data()),";dz(match vtx,vtx0) [mm]",1000,-10.,10.);
  h_hs_trk_dz_vtx0 = new TH1F(Form("%s_hs_trk_dz_vtx0",prefix.Data()),";dz(trk,vtx0) [mm]",1000,-10.,10.);
  h_mc_idx_duplicates = new TH1F(Form("%s_mc_idx_duplicates",prefix.Data()),";# tracks with duplicate MC matches",500,0,500);
  h_mc_idx_duplicates_dr = new TH1F(Form("%s_mc_idx_duplicates_dr",prefix.Data()),";dR(reco,gen) for tracks with duplicate matches",40,0.,0.2);
  h_mc_idx_duplicates_pt = new TH1F(Form("%s_mc_idx_duplicates_pt",prefix.Data()),";p_{T} for tracks with duplicate matches",500,0.,100.);
  h_mc_idx_duplicates_nhits = new TH1F(Form("%s_mc_idx_duplicates_nhits",prefix.Data()),";N(hits) for tracks with duplicate matches",20,0,20);
  h_match_dr = new TH1F(Form("%s_match_dr",prefix.Data()),";dR(reco,gen)",40,0.,0.2);
  h_nvtx = new TH1F(Form("%s_nvtx",prefix.Data()),";N(vtx)",50,0,50);
  h_gen_match_vtx = new TH1F(Form("%s_gen_match_vtx",prefix.Data()),";reco index of true PV",51,-1,50);

  // lepton iso, vs purity
  h_el_iso = new TH1F(Form("%s_el_iso",prefix.Data()),";electron reliso, no PU cor",100,0.,2.);
  h_el_iso_cor = new TH1F(Form("%s_el_iso_cor",prefix.Data()),";electron reliso, with PU cor",100,0.,2.);
  h_el_iso_vs_vtx0_purity_dz = new TH2F(Form("%s_el_iso_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;electron reliso, no PU cor",100,0,1.,100,0.,2.);
  h_el_iso_cor_vs_vtx0_purity_dz = new TH2F(Form("%s_el_iso_cor_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;electron reliso, with PU cor",100,0,1.,100,0.,2.);
  h_el_iso_vs_nvtx = new TH2F(Form("%s_el_iso_vs_nvtx",prefix.Data()),";nvtx;electron reliso, no PU cor",50,0,50,100,0.,2.);
  h_el_iso_cor_vs_nvtx = new TH2F(Form("%s_el_iso_cor_vs_nvtx",prefix.Data()),";nvtx;electron reliso, with PU cor",50,0,50,100,0.,2.);

  h_mu_iso = new TH1F(Form("%s_mu_iso",prefix.Data()),";mu reliso, no PU cor",100,0.,2.);
  h_mu_iso_cor = new TH1F(Form("%s_mu_iso_cor",prefix.Data()),";mu reliso, with PU cor",100,0.,2.);
  h_mu_iso_vs_vtx0_purity_dz = new TH2F(Form("%s_mu_iso_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;mu reliso, no PU cor",100,0,1.,100,0.,2.);
  h_mu_iso_cor_vs_vtx0_purity_dz = new TH2F(Form("%s_mu_iso_cor_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;mu reliso, with PU cor",100,0,1.,100,0.,2.);
  h_mu_iso_vs_nvtx = new TH2F(Form("%s_mu_iso_vs_nvtx",prefix.Data()),";nvtx;mu reliso, no PU cor",50,0,50,100,0.,2.);
  h_mu_iso_cor_vs_nvtx = new TH2F(Form("%s_mu_iso_cor_vs_nvtx",prefix.Data()),";nvtx;mu reliso, with PU cor",50,0,50,100,0.,2.);

  h_pfjet_beta = new TH1F(Form("%s_pfjet_beta",prefix.Data()),";pfjet beta",100,0.,1.);
  h_pfjet_beta_vs_vtx0_purity_dz = new TH2F(Form("%s_pfjet_beta_vs_vtx0_purity_dz",prefix.Data()),";vtx0 purity;pfjet beta",100,0,1.,100,0.,1.);
  h_pfjet_beta_vs_nvtx = new TH2F(Form("%s_pfjet_beta_vs_nvtx",prefix.Data()),";nvtx;pfjet beta",50,0,50,100,0.,1.);

 
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
  float AEff = 0.;
  if (etaAbs <= 1.0) AEff = 0.10;
  else if (etaAbs > 1.0 && etaAbs <= 1.479) AEff = 0.12;
  else if (etaAbs > 1.479 && etaAbs <= 2.0) AEff = 0.085;
  else if (etaAbs > 2.0 && etaAbs <= 2.2) AEff = 0.11;
  else if (etaAbs > 2.2 && etaAbs <= 2.3) AEff = 0.12;
  else if (etaAbs > 2.3 && etaAbs <= 2.4) AEff = 0.12;
  else if (etaAbs > 2.4) AEff = 0.13;

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

