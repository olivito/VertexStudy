#ifndef vertexStudyLooper_h
#define vertexStudyLooper_h

#include <vector>
#include <list>
#include <string>
#include <map>
#include <set>
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"

#include "TChain.h"
#include "TChainElement.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

/* ------------------------------------------------------------------------- */
/* ------------------------------------------------------------------------- */


class vertexStudyLooper
{
public: 
  vertexStudyLooper();
  ~vertexStudyLooper() {}

  int  ScanChain(TChain *chain, const TString& prefix = "" );
  void BookHistos (const TString& prefix);
  void InitBaby();
  int associateTrackToVertexSimple(const unsigned int itrk);
  float electronPFiso(const unsigned int index, const bool cor = false);
  float muonPFiso(const unsigned int imu, const bool cor = false);
  float photonHollowTrkIso(const unsigned int iph);
  float dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 );
  float genSumPt2();

  // Set globals
  void set_createTree   (bool  b)    { g_createTree   = b; }
  void set_version      (const char* v)    { g_version      = v; }
  void set_json         (const char* v)    { g_json         = v; }        

  // Baby ntuple methods
  void makeOutput (const TString& prefix);
  void closeOutput ();

private:

  // Globals
  bool  g_createTree;
  const char* g_version;
  const char* g_json;      
  bool initialized;
  bool isdata_;
  TFile* outFile;

  // histograms
  TH2F* h_vtx0_hardscatter_pt_vs_sumpt;
  TH2F* h_vtx0_hardscatter_pt_vs_sumpt_recalc;
  TH2F* h_vtx0_sumpt_vs_sumpt_recalc;
  TH1F* h_vtx_best_purity_dz;
  TH1F* h_vtx_best_purity_weight;
  TH2F* h_vtx_best_purity_dz_vs_nvtx;
  TH2F* h_vtx_best_purity_weight_vs_nvtx;
  TH1F* h_vtx0_purity_dz;
  TH1F* h_vtx0_purity_weight;
  TH2F* h_vtx0_purity_dz_vs_nvtx;
  TH2F* h_vtx0_purity_weight_vs_nvtx;
  TH1F* h_vtx_best_eff_dz;
  TH1F* h_vtx_best_eff_weight;
  TH2F* h_vtx_best_eff_dz_vs_nvtx;
  TH2F* h_vtx_best_eff_weight_vs_nvtx;
  TH1F* h_vtx0_eff_dz;
  TH1F* h_vtx0_eff_weight;
  TH2F* h_vtx0_eff_dz_vs_nvtx;
  TH2F* h_vtx0_eff_weight_vs_nvtx;
  TH1F* h_dz_trk_vtx;
  TH1F* h_dz_trk_vtx0_weight;
  TH1F* h_trk_bestdzvtx;
  TH2F* h_trk_bestdzvtx_vs_weightvtx;
  TH1F* h_trk_dz_vtxs;
  TH1F* h_hs_trk_dz_vtxs;
  TH1F* h_hs_trk_dz_vtx0;
  TH1F* h_mc_idx_duplicates;
  TH1F* h_mc_idx_duplicates_dr;
  TH1F* h_mc_idx_duplicates_pt;
  TH1F* h_mc_idx_duplicates_nhits;
  TH1F* h_match_dr;

  TH1F* h_nvtx;
  TH1F* h_ntracks;

  TH1F* h_gen_match_vtx;
  TH2F* h_gen_match_vtx_vs_nvtx;

  TH1F* h_genvtx_z_nocut;
  TH1F* h_genvtx_z;
  TH1F* h_genmatch_vtx_sumpt2;
  TH1F* h_genmatch_vtx_z;
  TH1F* h_recovtx_z;

  TH1F* h_genvtx_nomatch_z;
  TH1F* h_genvtx_nomatch_gensumpt2;
  TH1F* h_genvtx_nomatch_dz;
  TH1F* h_genvtx_nomatch_dz_eff10;
  TH1F* h_genvtx_nomatch_smalldz_eff;
  TH1F* h_genvtx_nomatch_largedz_nleps;
  TH1F* h_genvtx_nomatch_largedz_npartons;
  TH1F* h_genvtx_nomatch_nvtx;
  TH1F* h_genvtx_nomatch_ntracks;

  TH1F* h_genvtx_othermatch_z;
  TH1F* h_genvtx_othermatch_dz;
  TH1F* h_genvtx_othermatch_gensumpt2;
  TH1F* h_genvtx_othermatch_sumpt2;

  TH1F* h_vtx_sumpt2;
  TH1F* h_vtx_nohs_sumpt2;
  TH1F* h_vtx_nogen_sumpt2;
  TH1F* h_vtx_nogen_nohs_sumpt2;
  TH1F* h_vtx0_pu_sumpt2;

  TH1F* h_el_iso;
  TH1F* h_el_iso_cor;
  TH2F* h_el_iso_vs_vtx0_purity_dz;
  TH2F* h_el_iso_cor_vs_vtx0_purity_dz;
  TH2F* h_el_iso_vs_nvtx;
  TH2F* h_el_iso_cor_vs_nvtx;
  TH1F* h_el_trkiso;
  TH2F* h_el_trkiso_vs_vtx0_purity_dz;
  TH2F* h_el_trkiso_vs_nvtx;
  TH1F* h_el_trkiso_abs;
  TH2F* h_el_trkiso_abs_vs_vtx0_purity_dz;
  TH2F* h_el_trkiso_abs_vs_nvtx;

  TH1F* h_mu_iso;
  TH1F* h_mu_iso_cor;
  TH2F* h_mu_iso_vs_vtx0_purity_dz;
  TH2F* h_mu_iso_cor_vs_vtx0_purity_dz;
  TH2F* h_mu_iso_vs_nvtx;
  TH2F* h_mu_iso_cor_vs_nvtx;
  TH1F* h_mu_trkiso;
  TH2F* h_mu_trkiso_vs_vtx0_purity_dz;
  TH2F* h_mu_trkiso_vs_nvtx;
  TH1F* h_mu_trkiso_abs;
  TH2F* h_mu_trkiso_abs_vs_vtx0_purity_dz;
  TH2F* h_mu_trkiso_abs_vs_nvtx;

  TH1F* h_ph_trkiso;
  TH2F* h_ph_trkiso_vs_vtx0_purity_dz;
  TH2F* h_ph_trkiso_vs_nvtx;

  TH1F* h_pfjet_beta;
  TH2F* h_pfjet_beta_vs_vtx0_purity_dz;
  TH2F* h_pfjet_beta_vs_nvtx;

  TH1F* h_match_trk_pt;
  TH1F* h_match_trk_pt_low;
  TH1F* h_match_trk_pt_mid;

  TH1F* h_nomatch_trk_pt;
  TH1F* h_nomatch_trk_pt_low;
  TH1F* h_nomatch_trk_pt_mid;

  TH2F* h_nomatch_trk_pt_vs_eta;
  TH2F* h_nomatch_trk_pt_low_vs_eta;

};

#endif
