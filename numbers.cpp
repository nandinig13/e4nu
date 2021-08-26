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
//--------------------------
//no pions
Long64_t i  = tree->GetEntries("(nipip + nipim + nipi0)==0");

Long64_t i_0  = tree->GetEntries("(nipip + nipim + nipi0)==0 && (nfpip + nfpim + nfpi0)==0");
  Long64_t i_pi0  = tree->GetEntries("(nipip + nipim + nipi0)==0 && (nfpip + nfpim + nfpi0)==1 && nfpi0==1");
  Long64_t i_pip  = tree->GetEntries("(nipip + nipim + nipi0)==0 && (nfpip + nfpim + nfpi0)==1 && nfpip==1");
  Long64_t i_pim  = tree->GetEntries("(nipip + nipim + nipi0)==0 && (nfpip + nfpim + nfpi0)==1 && nfpim==1");
  Long64_t i_pix  = tree->GetEntries("(nipip + nipim + nipi0)==0 && (nfpip + nfpim + nfpi0)>1");

  cout << "initially no pions : " << i << endl;
  cout << "then no pions : " << i_0 << "so prop = " << i_0 / i << endl;
  cout << "then pi0 : " << i_pi0 << "so prop = " << i_pi0 / i << endl;
  cout << "then pip : " << i_pip << "so prop = " << i_pip / i << endl;
  cout << "then pim : " << i_pim << "so prop = " << i_pim / i << endl;
  cout << "then many pions : " << i_pix << "so prop = " << i_pix / i << endl;
 
}
