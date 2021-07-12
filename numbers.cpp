//entries counter
#include "TMath.h"
#include "TH1.h"
#include "TH3.h"
#include "TFile.h"
#include "TGraph.h"
#include "TChain.h"
#include "TMath.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TGaxis.h"
#include "TImage.h"
#include "TCanvas.h"
#include "TPaveLabel.h"
#include "TSystem.h"
#include "TRandom3.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include "acceptance_c.cpp"

using namespace std;
void entries(){

  string GST = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_2261GeV/apapadop_UpdatedSchwingerRad_SuSav2_C12_2261GeV.root";

TFile *f = new TFile(GST.c_str());

TTree *tree = (TTree*) f->Get("gst");

//tree->Scan("nip:nin:nipip:nipim:nipi0:", "res && nip ==2");
Long64_t nentries = tree->GetEntries();
Long64_t nentries_qel = tree->GetEntries("qel");
Long64_t nentries_res = tree->GetEntries("res");
Long64_t nentries_mec = tree->GetEntries("mec");
Long64_t nentries_dis = tree->GetEntries("dis");
Long64_t nentries_coh = tree->GetEntries("coh");
/* EXAMPLES
Long64_t nopiTOnopi = tree->GetEntries("nipip==0 && nipi0==0 && nipim==0 && nfpip==0 && nfpi0==0 && nfpim==0");
Long64_t nopiTO1pi0 = tree->GetEntries("nipip==0 && nipi0==0 && nipim==0 && nfpip==0 && nfpi0==1 && nfpim==0");
*/

//First stage - which nucleon is hit
Long64_t nentries_hitprot = tree->GetEntries("hitnuc==2212");
Long64_t nentries_hitneut = tree->GetEntries("hitnuc==2112");
Long64_t nentries_hitcorr = tree->GetEntries("hitnuc==2000000201");

//Proton hit - before FSI 
//Indents after FSI
//pix = more than one pion
Long64_t p_p_nop  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==0 && nfn==0");

Long64_t p_p  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0");
  Long64_t p_p_p  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_p_ppi0  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_p_ppip  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_p_ppim  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_p_ppix  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t p_p_n  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfn==1 && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_p_npi0  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_p_npip = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_p_npim = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_p_npix = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t p_p_nothing = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==0 && nfn==0 && (nfpip + nfpim + nfpi0)==0");
/*
Long64_t p_ppi0  = tree->GetEntries("hitnuc==2212 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1");
  Long64_t p_ppi0_p  = tree->GetEntries("hitnuc==2212 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_ppi0_ppi0  = tree->GetEntries("hitnuc==2212 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppi0_ppip  = tree->GetEntries("hitnuc==2212 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppi0_ppim  = tree->GetEntries("hitnuc==2212 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppi0_ppix  = tree->GetEntries("hitnuc==2212 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t p_ppi0_nx  = tree->GetEntries("hitnuc==2212 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t p_ppip = tree->GetEntries("hitnuc==2212 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1");
  Long64_t p_ppip_p  = tree->GetEntries("hitnuc==2212 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_ppip_ppi0  = tree->GetEntries("hitnuc==2212 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppip_ppip  = tree->GetEntries("hitnuc==2212 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppip_ppim  = tree->GetEntries("hitnuc==2212 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppip_ppix  = tree->GetEntries("hitnuc==2212 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t p_ppip_nx  = tree->GetEntries("hitnuc==2212 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t p_ppim  = tree->GetEntries("hitnuc==2212 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1");
  Long64_t p_ppim_p  = tree->GetEntries("hitnuc==2212 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_ppim_ppi0  = tree->GetEntries("hitnuc==2212 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppim_ppip  = tree->GetEntries("hitnuc==2212 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppim_ppim  = tree->GetEntries("hitnuc==2212 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppim_ppix  = tree->GetEntries("hitnuc==2212 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t p_ppim_nx  = tree->GetEntries("hitnuc==2212 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t p_ppix = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)>1");
  Long64_t p_ppix_p  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_ppix_ppi0  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppix_ppip  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppix_ppim  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_ppix_ppix  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t p_ppix_nx  = tree->GetEntries("hitnuc==2212 && nip==1 && (nipip + nipim + nipi0)>1 && nfn==1");


Long64_t p_n  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)==0");
  Long64_t p_n_n  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_n_npi0  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_n_npip  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_n_npim  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_n_npix  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t p_n_px  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)==0 && nfp==1");

Long64_t p_npip = tree->GetEntries("hitnuc==2212 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1");
  Long64_t p_npip_n  = tree->GetEntries("hitnuc==2212 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_npip_npi0  = tree->GetEntries("hitnuc==2212 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npip_npip  = tree->GetEntries("hitnuc==2212 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npip_npim  = tree->GetEntries("hitnuc==2212 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npip_npix  = tree->GetEntries("hitnuc==2212 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t p_npip_px  = tree->GetEntries("hitnuc==2212 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1");

Long64_t p_npim  = tree->GetEntries("hitnuc==2212 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1");
  Long64_t p_npim_n  = tree->GetEntries("hitnuc==2212 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_npim_npi0  = tree->GetEntries("hitnuc==2212 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npim_npip  = tree->GetEntries("hitnuc==2212 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npim_npim  = tree->GetEntries("hitnuc==2212 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npim_npix  = tree->GetEntries("hitnuc==2212 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
 Long64_t p_npim_nx  = tree->GetEntries("hitnuc==2212 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1");


Long64_t p_npi0 = tree->GetEntries("hitnuc==2212 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1");
  Long64_t p_npi0_n  = tree->GetEntries("hitnuc==2212 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_npi0_npi0  = tree->GetEntries("hitnuc==2212 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npi0_npip  = tree->GetEntries("hitnuc==2212 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npi0_npim  = tree->GetEntries("hitnuc==2212 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npi0_npix  = tree->GetEntries("hitnuc==2212 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t p_npi0_px  = tree->GetEntries("hitnuc==2212 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1");


Long64_t p_npix = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)>1");
  Long64_t p_npix_n  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t p_npix_npi0  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npix_npip  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npix_npim  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t p_npix_npix  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t p_npix_px  = tree->GetEntries("hitnuc==2212 && nin==1 && (nipip + nipim + nipi0)>1 && nfp==1");




//Neutron hit - before FSI
Long64_t n_p  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)==0");
  Long64_t n_p_p  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t n_p_ppi0  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_p_ppip  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_p_ppim  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_p_ppix  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t n_p_nx  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)==0 && nfn==1");

Long64_t n_ppi0  = tree->GetEntries("hitnuc==2112 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1");
  Long64_t n_ppi0_p  = tree->GetEntries("hitnuc==2112 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t n_ppi0_ppi0  = tree->GetEntries("hitnuc==2112 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppi0_ppip  = tree->GetEntries("hitnuc==2112 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppi0_ppim  = tree->GetEntries("hitnuc==2112 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppi0_ppix  = tree->GetEntries("hitnuc==2112 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t n_ppi0_nx  = tree->GetEntries("hitnuc==2112 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t n_ppip = tree->GetEntries("hitnuc==2112 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1");
  Long64_t n_ppip_p  = tree->GetEntries("hitnuc==2112 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t n_ppip_ppi0  = tree->GetEntries("hitnuc==2112 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppip_ppip  = tree->GetEntries("hitnuc==2112 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppip_ppim  = tree->GetEntries("hitnuc==2112 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppip_ppix  = tree->GetEntries("hitnuc==2112 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t n_ppip_nx  = tree->GetEntries("hitnuc==2112 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t n_ppim  = tree->GetEntries("hitnuc==2112 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1");
  Long64_t n_ppim_p  = tree->GetEntries("hitnuc==2112 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t n_ppim_ppi0  = tree->GetEntries("hitnuc==2112 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppim_ppip  = tree->GetEntries("hitnuc==2112 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppim_ppim  = tree->GetEntries("hitnuc==2112 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppim_ppix  = tree->GetEntries("hitnuc==2112 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t n_ppim_nx  = tree->GetEntries("hitnuc==2112 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t n_ppix = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)>1");
  Long64_t n_ppix_p  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t n_ppix_ppi0  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppix_ppip  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppix_ppim  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_ppix_ppix  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t n_ppix_nx  = tree->GetEntries("hitnuc==2112 && nip==1 && (nipip + nipim + nipi0)>1 && nfn==1");


Long64_t n_n  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)==0");
  Long64_t n_n_n  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t n_n_npi0  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_n_npip  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_n_npim  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_n_npix  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t n_n_px  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)==0 && nfp==1");

Long64_t n_npip = tree->GetEntries("hitnuc==2112 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1");
  Long64_t n_npip_n  = tree->GetEntries("hitnuc==2112 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t n_npip_npi0  = tree->GetEntries("hitnuc==2112 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npip_npip  = tree->GetEntries("hitnuc==2112 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npip_npim  = tree->GetEntries("hitnuc==2112 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npip_npix  = tree->GetEntries("hitnuc==2112 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t n_npip_px  = tree->GetEntries("hitnuc==2112 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t n_npim  = tree->GetEntries("hitnuc==2112 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1");
  Long64_t n_npim_n  = tree->GetEntries("hitnuc==2112 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t n_npim_npi0  = tree->GetEntries("hitnuc==2112 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npim_npip  = tree->GetEntries("hitnuc==2112 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npim_npim  = tree->GetEntries("hitnuc==2112 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npim_npix  = tree->GetEntries("hitnuc==2112 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
 Long64_t n_npim_nx  = tree->GetEntries("hitnuc==2112 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1");


Long64_t n_npi0 = tree->GetEntries("hitnuc==2112 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1");
  Long64_t n_npi0_n  = tree->GetEntries("hitnuc==2112 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t n_npi0_npi0  = tree->GetEntries("hitnuc==2112 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npi0_npip  = tree->GetEntries("hitnuc==2112 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npi0_npim  = tree->GetEntries("hitnuc==2112 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npi0_npix  = tree->GetEntries("hitnuc==2112 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t n_npi0_px  = tree->GetEntries("hitnuc==2112 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1");


Long64_t n_npix = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)>1");
  Long64_t n_npix_n  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t n_npix_npi0  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npix_npip  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npix_npim  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t n_npix_npix  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t n_npix_px  = tree->GetEntries("hitnuc==2112 && nin==1 && (nipip + nipim + nipi0)>1 && nfp==1");


//Correlated hit - before FSI
Long64_t c_p  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)==0");
  Long64_t c_p_p  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t c_p_ppi0  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_p_ppip  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_p_ppim  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_p_ppix  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)==0 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t c_p_nx  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)==0 && nfn==1");

Long64_t c_ppi0  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1");
  Long64_t c_ppi0_p  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t c_ppi0_ppi0  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppi0_ppip  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppi0_ppim  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppi0_ppix  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t c_ppi0_nx  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t c_ppip = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1");
  Long64_t c_ppip_p  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t c_ppip_ppi0  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppip_ppip  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppip_ppim  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppip_ppix  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t c_ppip_nx  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t c_ppim  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1");
  Long64_t c_ppim_p  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t c_ppim_ppi0  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppim_ppip  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppim_ppim  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppim_ppix  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t c_ppim_nx  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t c_ppix = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)>1");
  Long64_t c_ppix_p  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t c_ppix_ppi0  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppix_ppip  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppix_ppim  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_ppix_ppix  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)>1 && nfp==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t c_ppix_nx  = tree->GetEntries("hitnuc==2000000201 && nip==1 && (nipip + nipim + nipi0)>1 && nfn==1");


Long64_t c_n  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)==0");
  Long64_t c_n_n  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t c_n_npi0  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_n_npip  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_n_npim  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_n_npix  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)==0 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t c_n_px  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)==0 && nfp==1");

Long64_t c_npip = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1");
  Long64_t c_npip_n  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t c_npip_npi0  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npip_npip  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npip_npim  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npip_npix  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t c_npip_px  = tree->GetEntries("hitnuc==2000000201 && nip==1 && nipip==1 && (nipip + nipim + nipi0)==1 && nfn==1");

Long64_t c_npim  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1");
  Long64_t c_npim_n  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t c_npim_npi0  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npim_npip  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npim_npim  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npim_npix  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
 Long64_t c_npim_nx  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipim==1 && (nipip + nipim + nipi0)==1 && nfp==1");


Long64_t c_npi0 = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1");
  Long64_t c_npi0_n  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t c_npi0_npi0  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npi0_npip  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npi0_npim  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npi0_npix  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t c_npi0_px  = tree->GetEntries("hitnuc==2000000201 && nin==1 && nipi0==1 && (nipip + nipim + nipi0)==1 && nfp==1");


Long64_t c_npix = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)>1");
  Long64_t c_npix_n  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1  && (nfpip + nfpim + nfpi0)==0");
  Long64_t c_npix_npi0  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npix_npip  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && nfpip==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npix_npim  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1");
  Long64_t c_npix_npix  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)>1 && nfn==1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t c_npix_px  = tree->GetEntries("hitnuc==2000000201 && nin==1 && (nipip + nipim + nipi0)>1 && nfp==1");
*/

cout << GST << endl;
cout << "number of entries" << nentries <<endl;
cout << "number of QE events" << nentries_qel <<endl;
cout << "QE prop" << int(double(nentries_qel) / double(nentries) *100.) << " \%" << endl;
cout << "number of RES events" << nentries_res <<endl;
cout << "RES prop" << int(double(nentries_res) / double(nentries) *100.) << " \%" << endl;
cout << "number of MEC events" << nentries_mec <<endl;
cout << "MEC prop" << int(double(nentries_mec) / double(nentries) *100.) << " \%" << endl;
cout << "number of DIS events" << nentries_dis <<endl;
cout << "DIS prop" << int(double(nentries_dis) / double(nentries) *100.) << " \%" << endl;
cout << "number of coherent events" << nentries_coh <<endl;
cout << "COH prop" << int(double(nentries_coh) / double(nentries) *100.) << " \%" << endl;

cout<<endl;

cout << "number of protons hit" << nentries_hitprot <<endl;
cout << "Hit protons prop" << int(double(nentries_hitprot) / double(nentries) *100.) << " \%" << endl;
cout << "number of neutrons hit" << nentries_hitneut <<endl;
cout << "Hit neutrons prop" << int(double(nentries_hitneut) / double(nentries) *100.) << " \%" << endl;
cout << "number of correlated pairs hit" << nentries_hitcorr <<endl;
cout << "Hit correlated pairs prop" << int(double(nentries_hitcorr) / double(nentries) *100.) << " \%" << endl;

cout<<endl;

cout << "1. proton-> proton" << p_p << endl; 
cout<< "2. proton -> proton -> proton" << p_p_p << endl; 
cout<< "2. proton -> proton -> proton,pi0:" << p_p_ppi0 << endl;
cout<< "2. proton -> proton -> proton,pip" << p_p_ppip << endl;
cout<< "2. proton -> proton -> proton,pim" << p_p_ppim << endl;
cout<< "2. proton -> proton -> proton,more than one pion" << p_p_ppix << endl;
cout<< "2. proton -> proton -> neutron" << p_p_n << endl;
cout<< "2. proton -> proton -> neutron,pi0:" << p_p_npi0 << endl;
cout<< "2. proton -> proton -> neutron,pip" << p_p_npip << endl;
cout<< "2. proton -> proton -> neutron,pim" << p_p_npim << endl;
cout<< "2. proton -> proton -> neutron,pix" << p_p_npix << endl;
cout<< "2. proton -> proton -> nothing" << p_p_nothing << endl;

/*
cout << "1. proton-> proton, pi0:" << p_ppi0 << endl; 
cout<< "2. proton -> proton, pi0 -> proton" << p_ppi0_p << endl; 
cout<< "2. proton -> proton, pi0 -> proton,pi0:" << p_ppi0_ppi0 << endl;
cout<< "2. proton -> proton, pi0 -> proton,pip" << p_ppi0_ppip << endl;
cout<< "2. proton -> proton, pi0 -> proton,pim" << p_ppi0_ppim << endl;
cout<< "2. proton -> proton, pi0 -> proton,more than one pion" << p_ppi0_ppix << endl;
cout<< "2. proton -> proton, pi0 -> neutron, maybe others" << p_ppi0_nx << endl;

cout << "1. proton-> proton, pip" << p_ppip << endl; 
cout<< "2. proton -> proton, pip -> proton" << p_ppip_p << endl; 
cout<< "2. proton -> proton, pip -> proton,pi0:" << p_ppip_ppi0 << endl;
cout<< "2. proton -> proton, pip -> proton,pip" << p_ppip_ppip << endl;
cout<< "2. proton -> proton, pip -> proton,pim" << p_ppip_ppim << endl;
cout<< "2. proton -> proton, pip -> proton,more than one pion" << p_ppip_ppix << endl;
cout<< "2. proton -> proton, pip -> neutron, maybe others" << p_ppip_nx << endl;

cout << "1. proton-> proton, pim" << p_ppim << endl; 
cout<< "2. proton -> proton, pim -> proton" << p_ppim_p << endl; 
cout<< "2. proton -> proton, pim -> proton,pi0:" << p_ppim_ppi0 << endl;
cout<< "2. proton -> proton, pim -> proton,pip" << p_ppim_ppip << endl;
cout<< "2. proton -> proton, pim -> proton,pim" << p_ppim_ppim << endl;
cout<< "2. proton -> proton, pim -> proton,more than one pion" << p_ppim_ppix << endl;
cout<< "2. proton -> proton, pim -> neutron, maybe others" << p_ppim_nx << endl;

cout << "1. proton-> proton, pix" << p_ppix << endl; 
cout<< "2. proton -> proton, pix -> proton" << p_ppix_p << endl; 
cout<< "2. proton -> proton, pix -> proton,pi0:" << p_ppix_ppi0 << endl;
cout<< "2. proton -> proton, pix -> proton,pip" << p_ppix_ppip << endl;
cout<< "2. proton -> proton, pix -> proton,pim" << p_ppix_ppim << endl;
cout<< "2. proton -> proton, pix -> proton,more than one pion" << p_ppix_ppix << endl;
cout<< "2. proton -> proton, pix -> neutron, maybe others" << p_ppix_nx << endl;

cout << "1. proton-> neutron" << p_n << endl; 
cout<< "2. proton -> neutron -> neutron" << p_n_n << endl; 
cout<< "2. proton -> neutron -> neutron,pi0:" << p_n_npi0 << endl;
cout<< "2. proton -> neutron -> neutron,pip" << p_n_npip << endl;
cout<< "2. proton -> neutron -> neutron,pim" << p_n_npim << endl;
cout<< "2. proton -> neutron -> neutron,more than one pion" << p_n_npix << endl;
cout<< "2. proton -> neutron -> neutron, maybe others" << p_n_px << endl;

cout << "1. proton-> neutron, pi0:" << p_npi0 << endl; 
cout<< "2. proton -> neutron, pi0 -> neutron" << p_npi0_n << endl; 
cout<< "2. proton -> neutron, pi0 -> neutron,pi0:" << p_npi0_npi0 << endl;
cout<< "2. proton -> neutron, pi0 -> neutron,pip" << p_npi0_npip << endl;
cout<< "2. proton -> neutron, pi0 -> neutron,pim" << p_npi0_npim << endl;
cout<< "2. proton -> neutron, pi0 -> neutron,more than one pion" << p_npi0_npix << endl;
cout<< "2. proton -> neutron, pi0 -> proton, maybe others" << p_npi0_px << endl;

cout << "1. proton-> neutron, pip" << p_npip << endl; 
cout<< "2. proton -> neutron, pip -> neutron" << p_npip_n << endl; 
cout<< "2. proton -> neutron, pip -> neutron,pi0:" << p_npip_npi0 << endl;
cout<< "2. proton -> neutron, pip -> neutron,pip" << p_npip_npip << endl;
cout<< "2. proton -> neutron, pip -> neutron,pim" << p_npip_npim << endl;
cout<< "2. proton -> neutron, pip -> neutron,more than one pion" << p_npip_npix << endl;
cout<< "2. proton -> neutron, pip -> proton, maybe others" << p_npip_px << endl;

cout << "1. proton-> neutron, pim" << p_npim << endl; 
cout<< "2. proton -> neutron, pim -> neutron" << p_npim_n << endl; 
cout<< "2. proton -> neutron, pim -> neutron,pi0:" << p_npim_npi0 << endl;
cout<< "2. proton -> neutron, pim -> neutron,pip" << p_npim_npip << endl;
cout<< "2. proton -> neutron, pim -> neutron,pim" << p_npim_npim << endl;
cout<< "2. proton -> neutron, pim -> neutron,more than one pion" << p_npim_npix << endl;
cout<< "2. proton -> neutron, pim -> proton, maybe others" << p_npim_nx << endl;

cout << "1. proton-> neutron, pix" << p_npix << endl; 
cout<< "2. proton -> neutron, pix -> neutron" << p_npix_n << endl; 
cout<< "2. proton -> neutron, pix -> neutron,pi0:" << p_npix_npi0 << endl;
cout<< "2. proton -> neutron, pix -> neutron,pip" << p_npix_npip << endl;
cout<< "2. proton -> neutron, pix -> neutron,pim" << p_npix_npim << endl;
cout<< "2. proton -> neutron, pix -> neutron,more than one pion" << p_npix_npix << endl;
cout<< "2. proton -> neutron, pix -> proton, maybe others" << p_npix_px << endl;

cout << "-------------------------------------------------------------------------------" << endl;

cout << "1. neutron-> proton" << n_p << endl; 
cout<< "2. neutron -> proton -> proton" << n_p_p << endl; 
cout<< "2. neutron -> proton -> proton,pi0:" << n_p_ppi0 << endl;
cout<< "2. neutron -> proton -> proton,pip" << n_p_ppip << endl;
cout<< "2. neutron -> proton -> proton,pim" << n_p_ppim << endl;
cout<< "2. neutron -> proton -> proton,more than one pion" << n_p_ppix << endl;
cout<< "2. neutron -> proton -> neutron, maybe others" << n_p_nx << endl;

cout << "1. neutron-> proton, pi0:" << n_ppi0 << endl; 
cout<< "2. neutron -> proton, pi0 -> proton" << n_ppi0_p << endl; 
cout<< "2. neutron -> proton, pi0 -> proton,pi0:" << n_ppi0_ppi0 << endl;
cout<< "2. neutron -> proton, pi0 -> proton,pip" << n_ppi0_ppip << endl;
cout<< "2. neutron -> proton, pi0 -> proton,pim" << n_ppi0_ppim << endl;
cout<< "2. neutron -> proton, pi0 -> proton,more than one pion" << n_ppi0_ppix << endl;
cout<< "2. neutron -> proton, pi0 -> neutron, maybe others" << n_ppi0_nx << endl;

cout << "1. neutron-> proton, pip" << n_ppip << endl; 
cout<< "2. neutron -> proton, pip -> proton" << n_ppip_p << endl; 
cout<< "2. neutron -> proton, pip -> proton,pi0:" << n_ppip_ppi0 << endl;
cout<< "2. neutron -> proton, pip -> proton,pip" << n_ppip_ppip << endl;
cout<< "2. neutron -> proton, pip -> proton,pim" << n_ppip_ppim << endl;
cout<< "2. neutron -> proton, pip -> proton,more than one pion" << n_ppip_ppix << endl;
cout<< "2. neutron -> proton, pip -> neutron, maybe others" << n_ppip_nx << endl;

cout << "1. neutron-> proton, pim" << n_ppim << endl; 
cout<< "2. neutron -> proton, pim -> proton" << n_ppim_p << endl; 
cout<< "2. neutron -> proton, pim -> proton,pi0:" << n_ppim_ppi0 << endl;
cout<< "2. neutron -> proton, pim -> proton,pip" << n_ppim_ppip << endl;
cout<< "2. neutron -> proton, pim -> proton,pim" << n_ppim_ppim << endl;
cout<< "2. neutron -> proton, pim -> proton,more than one pion" << n_ppim_ppix << endl;
cout<< "2. neutron -> proton, pim -> neutron, maybe others" << n_ppim_nx << endl;

cout << "1. neutron-> proton, pix" << n_ppix << endl; 
cout<< "2. neutron -> proton, pix -> proton" << n_ppix_p << endl; 
cout<< "2. neutron -> proton, pix -> proton,pi0:" << n_ppix_ppi0 << endl;
cout<< "2. neutron -> proton, pix -> proton,pip" << n_ppix_ppip << endl;
cout<< "2. neutron -> proton, pix -> proton,pim" << n_ppix_ppim << endl;
cout<< "2. neutron -> proton, pix -> proton,more than one pion" << n_ppix_ppix << endl;
cout<< "2. neutron -> proton, pix -> neutron, maybe others" << n_ppix_nx << endl;

cout << "1. neutron-> neutron" << n_n << endl; 
cout<< "2. neutron -> neutron -> neutron" << n_n_n << endl; 
cout<< "2. neutron -> neutron -> neutron,pi0:" << n_n_npi0 << endl;
cout<< "2. neutron -> neutron -> neutron,pip" << n_n_npip << endl;
cout<< "2. neutron -> neutron -> neutron,pim" << n_n_npim << endl;
cout<< "2. neutron -> neutron -> neutron,more than one pion" << n_n_npix << endl;
cout<< "2. neutron -> neutron -> neutron, maybe others" << n_n_px << endl;

cout << "1. neutron-> neutron, pi0:" << n_npi0 << endl; 
cout<< "2. neutron -> neutron, pi0 -> neutron" << n_npi0_n << endl; 
cout<< "2. neutron -> neutron, pi0 -> neutron,pi0:" << n_npi0_npi0 << endl;
cout<< "2. neutron -> neutron, pi0 -> neutron,pip" << n_npi0_npip << endl;
cout<< "2. neutron -> neutron, pi0 -> neutron,pim" << n_npi0_npim << endl;
cout<< "2. neutron -> neutron, pi0 -> neutron,more than one pion" << n_npi0_npix << endl;
cout<< "2. neutron -> neutron, pi0 -> proton, maybe others" << n_npi0_px << endl;

cout << "1. neutron-> neutron, pip" << n_npip << endl; 
cout<< "2. neutron -> neutron, pip -> neutron" << n_npip_n << endl; 
cout<< "2. neutron -> neutron, pip -> neutron,pi0:" << n_npip_npi0 << endl;
cout<< "2. neutron -> neutron, pip -> neutron,pip" << n_npip_npip << endl;
cout<< "2. neutron -> neutron, pip -> neutron,pim" << n_npip_npim << endl;
cout<< "2. neutron -> neutron, pip -> neutron,more than one pion" << n_npip_npix << endl;
cout<< "2. neutron -> neutron, pip -> proton, maybe others" << n_npip_px << endl;

cout << "1. neutron-> neutron, pim" << n_npim << endl; 
cout<< "2. neutron -> neutron, pim -> neutron" << n_npim_n << endl; 
cout<< "2. neutron -> neutron, pim -> neutron,pi0:" << n_npim_npi0 << endl;
cout<< "2. neutron -> neutron, pim -> neutron,pip" << n_npim_npip << endl;
cout<< "2. neutron -> neutron, pim -> neutron,pim" << n_npim_npim << endl;
cout<< "2. neutron -> neutron, pim -> neutron,more than one pion" << n_npim_npix << endl;
cout<< "2. neutron -> neutron, pim -> proton, maybe others" << n_npim_nx << endl;

cout << "1. neutron-> neutron, pix" << n_npix << endl; 
cout<< "2. neutron -> neutron, pix -> neutron" << n_npix_n << endl; 
cout<< "2. neutron -> neutron, pix -> neutron,pi0:" << n_npix_npi0 << endl;
cout<< "2. neutron -> neutron, pix -> neutron,pip" << n_npix_npip << endl;
cout<< "2. neutron -> neutron, pix -> neutron,pim" << n_npix_npim << endl;
cout<< "2. neutron -> neutron, pix -> neutron,more than one pion" << n_npix_npix << endl;
cout<< "2. neutron -> neutron, pix -> proton, maybe others" << n_npix_px << endl;

cout << "-------------------------------------------------------------------------------" << endl;

cout << "1. correlated-> proton" << c_p << endl; 
cout<< "2. correlated -> proton -> proton" << c_p_p << endl; 
cout<< "2. correlated -> proton -> proton,pi0:" << c_p_ppi0 << endl;
cout<< "2. correlated -> proton -> proton,pip" << c_p_ppip << endl;
cout<< "2. correlated -> proton -> proton,pim" << c_p_ppim << endl;
cout<< "2. correlated -> proton -> proton,more than one pion" << c_p_ppix << endl;
cout<< "2. correlated -> proton -> neutron, maybe others" << c_p_nx << endl;

cout << "1. correlated-> proton, pi0:" << c_ppi0 << endl; 
cout<< "2. correlated -> proton, pi0 -> proton" << c_ppi0_p << endl; 
cout<< "2. correlated -> proton, pi0 -> proton,pi0:" << c_ppi0_ppi0 << endl;
cout<< "2. correlated -> proton, pi0 -> proton,pip" << c_ppi0_ppip << endl;
cout<< "2. correlated -> proton, pi0 -> proton,pim" << c_ppi0_ppim << endl;
cout<< "2. correlated -> proton, pi0 -> proton,more than one pion" << c_ppi0_ppix << endl;
cout<< "2. correlated -> proton, pi0 -> neutron, maybe others" << c_ppi0_nx << endl;

cout << "1. correlated-> proton, pip" << c_ppip << endl; 
cout<< "2. correlated -> proton, pip -> proton" << c_ppip_p << endl; 
cout<< "2. correlated -> proton, pip -> proton,pi0:" << c_ppip_ppi0 << endl;
cout<< "2. correlated -> proton, pip -> proton,pip" << c_ppip_ppip << endl;
cout<< "2. correlated -> proton, pip -> proton,pim" << c_ppip_ppim << endl;
cout<< "2. correlated -> proton, pip -> proton,more than one pion" << c_ppip_ppix << endl;
cout<< "2. correlated -> proton, pip -> neutron, maybe others" << c_ppip_nx << endl;

cout << "1. correlated-> proton, pim" << c_ppim << endl; 
cout<< "2. correlated -> proton, pim -> proton" << c_ppim_p << endl; 
cout<< "2. correlated -> proton, pim -> proton,pi0:" << c_ppim_ppi0 << endl;
cout<< "2. correlated -> proton, pim -> proton,pip" << c_ppim_ppip << endl;
cout<< "2. correlated -> proton, pim -> proton,pim" << c_ppim_ppim << endl;
cout<< "2. correlated -> proton, pim -> proton,more than one pion" << c_ppim_ppix << endl;
cout<< "2. correlated -> proton, pim -> neutron, maybe others" << c_ppim_nx << endl;

cout << "1. correlated-> proton, pix" << c_ppix << endl; 
cout<< "2. correlated -> proton, pix -> proton" << c_ppix_p << endl; 
cout<< "2. correlated -> proton, pix -> proton,pi0:" << c_ppix_ppi0 << endl;
cout<< "2. correlated -> proton, pix -> proton,pip" << c_ppix_ppip << endl;
cout<< "2. correlated -> proton, pix -> proton,pim" << c_ppix_ppim << endl;
cout<< "2. correlated -> proton, pix -> proton,more than one pion" << c_ppix_ppix << endl;
cout<< "2. correlated -> proton, pix -> neutron, maybe others" << c_ppix_nx << endl;

cout << "1. correlated-> neutron" << c_n << endl; 
cout<< "2. correlated -> neutron -> neutron" << c_n_n << endl; 
cout<< "2. correlated -> neutron -> neutron,pi0:" << c_n_npi0 << endl;
cout<< "2. correlated -> neutron -> neutron,pip" << c_n_npip << endl;
cout<< "2. correlated -> neutron -> neutron,pim" << c_n_npim << endl;
cout<< "2. correlated -> neutron -> neutron,more than one pion" << c_n_npix << endl;
cout<< "2. correlated -> neutron -> neutron, maybe others" << c_n_px << endl;

cout << "1. correlated-> neutron, pi0:" << c_npi0 << endl; 
cout<< "2. correlated -> neutron, pi0 -> neutron" << c_npi0_n << endl; 
cout<< "2. correlated -> neutron, pi0 -> neutron,pi0:" << c_npi0_npi0 << endl;
cout<< "2. correlated -> neutron, pi0 -> neutron,pip" << c_npi0_npip << endl;
cout<< "2. correlated -> neutron, pi0 -> neutron,pim" << c_npi0_npim << endl;
cout<< "2. correlated -> neutron, pi0 -> neutron,more than one pion" << c_npi0_npix << endl;
cout<< "2. correlated -> neutron, pi0 -> proton, maybe others" << c_npi0_px << endl;

cout << "1. correlated-> neutron, pip" << c_npip << endl; 
cout<< "2. correlated -> neutron, pip -> neutron" << c_npip_n << endl; 
cout<< "2. correlated -> neutron, pip -> neutron,pi0:" << c_npip_npi0 << endl;
cout<< "2. correlated -> neutron, pip -> neutron,pip" << c_npip_npip << endl;
cout<< "2. correlated -> neutron, pip -> neutron,pim" << c_npip_npim << endl;
cout<< "2. correlated -> neutron, pip -> neutron,more than one pion" << c_npip_npix << endl;
cout<< "2. correlated -> neutron, pip -> proton, maybe others" << c_npip_px << endl;

cout << "1. correlated-> neutron, pim" << c_npim << endl; 
cout<< "2. correlated -> neutron, pim -> neutron" << c_npim_n << endl; 
cout<< "2. correlated -> neutron, pim -> neutron,pi0:" << c_npim_npi0 << endl;
cout<< "2. correlated -> neutron, pim -> neutron,pip" << c_npim_npip << endl;
cout<< "2. correlated -> neutron, pim -> neutron,pim" << c_npim_npim << endl;
cout<< "2. correlated -> neutron, pim -> neutron,more than one pion" << c_npim_npix << endl;
cout<< "2. correlated -> neutron, pim -> proton, maybe others" << c_npim_nx << endl;

cout << "1. correlated-> neutron, pix" << c_npix << endl; 
cout<< "2. correlated -> neutron, pix -> neutron" << c_npix_n << endl; 
cout<< "2. correlated -> neutron, pix -> neutron,pi0:" << c_npix_npi0 << endl;
cout<< "2. correlated -> neutron, pix -> neutron,pip" << c_npix_npip << endl;
cout<< "2. correlated -> neutron, pix -> neutron,pim" << c_npix_npim << endl;
cout<< "2. correlated -> neutron, pix -> neutron,more than one pion" << c_npix_npix << endl;
//cout<< "2. correlated -> neutron, pix -> proton, maybe others" << c_npix_px << endl;

*/
}
