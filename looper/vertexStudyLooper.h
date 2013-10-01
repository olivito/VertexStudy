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
//typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > P4;
typedef std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;
//typedef map<unsigned int, unsigned int> m_uiui;

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
  float dRbetweenVectors(LorentzVector vec1, LorentzVector vec2 );

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
  TH1F* h_gen_match_vtx;

  TH1F* h_el_iso;
  TH1F* h_el_iso_cor;
  TH2F* h_el_iso_vs_vtx0_purity_dz;
  TH2F* h_el_iso_cor_vs_vtx0_purity_dz;
  TH2F* h_el_iso_vs_nvtx;
  TH2F* h_el_iso_cor_vs_nvtx;

  TH1F* h_mu_iso;
  TH1F* h_mu_iso_cor;
  TH2F* h_mu_iso_vs_vtx0_purity_dz;
  TH2F* h_mu_iso_cor_vs_vtx0_purity_dz;
  TH2F* h_mu_iso_vs_nvtx;
  TH2F* h_mu_iso_cor_vs_nvtx;

  TH1F* h_pfjet_beta;
  TH2F* h_pfjet_beta_vs_vtx0_purity_dz;
  TH2F* h_pfjet_beta_vs_nvtx;
};

#endif
