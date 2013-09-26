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

//#include "../CORE/topmass/ttdilepsolve.h" REPLACETOPMASS

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
        TH2F* h_vtx_hardscatter_pt_vs_sumpt;
        TH2F* h_vtx_hardscatter_pt_vs_sumpt_recalc;
        TH2F* h_vtx_sumpt_vs_sumpt_recalc;
        TH1F* h_vtx_fracpt_hardscatter_dz;
        TH1F* h_vtx_fracpt_hardscatter_weight;
        TH2F* h_vtx_fracpt_hardscatter_vs_nvtx;
        TH1F* h_vtx_hardscatter_fracpt_dz;
        TH1F* h_vtx_hardscatter_fracpt_weight;
        TH2F* h_vtx_hardscatter_fracpt_vs_nvtx;
        TH1F* h_dz_trk_vtx;
        TH1F* h_trk_bestdzvtx;
        TH2F* h_trk_bestdzvtx_vs_weightvtx;
        TH1F* h_trk_dz_vtxs;
        TH1F* h_hs_trk_dz_vtxs;
        TH1F* h_hs_trk_dz_vtx0;
        TH1F* h_mc_idx_duplicates;
        TH1F* h_mc_idx_duplicates_dr;

};

#endif
