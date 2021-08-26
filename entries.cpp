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

//tree->Scan("nf:pdgf");

//tree->Scan("nf:nfp:nfn:nfpip:nfpim:nfpi0");
/*
Long64_t nentries = tree->GetEntries();
Long64_t nentries_qel = tree->GetEntries("qel");
Long64_t nentries_res = tree->GetEntries("res");
Long64_t nentries_mec = tree->GetEntries("mec");
Long64_t nentries_dis = tree->GetEntries("dis");
Long64_t nentries_coh = tree->GetEntries("coh");

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


Long64_t nentries_1p1pip  = tree->GetEntries("nfp == 1 && (nipip + nfpim)==1 && nfpip==1");
Long64_t nentries_1p1pim  = tree->GetEntries("nfp == 1 && (nipip + nfpim)==1 && nfpim==1");

cout << "1p1pip: " << nentries_1p1pip << endl;
cout << "1p1pim: " << nentries_1p1pim << endl;

 Long64_t nentries_1n1pi  = tree->GetEntries("nip ==0 && && nin==1");

Long64_t nentries_1pi_1p = tree->GetEntries("nip ==0 && nin==1 && nfp == 1 && nfn == 0 && (nfpip+nfpim+nfpi0)==0");
Long64_t nentries_1pim_1p1pim = tree->GetEntries("nip ==0 && nin==1 && nfp==1 && nfn == 0 && nfpim==1 && (nfpip+nfpim+nfpi0)==1");
Long64_t nentries_1pim_1p1pip = tree->GetEntries("nip ==0 && nin==1 && nfp==1 && nfn == 0 && nfpip==1 && (nfpip+nfpim+nfpi0)==1");
Long64_t nentries_1pim_1p1pi0 = tree->GetEntries("nip ==0 && nin==1 && nfp==1 && nfn == 0 && nfpi0==1 && (nfpip+nfpim+nfpi0)==1");
Long64_t nentries_1pim_1pnpi = tree->GetEntries("nip ==0  && nin==1 && nfp==1 && nfn == 0 && (nfpip+nfpim+nfpi0)>1");
Long64_t nentries_1pim_np = tree->GetEntries("nip ==0 && nin==1 && nfp > 1 && nfn==0");
Long64_t nentries_1pim_n = tree->GetEntries("nip ==0 && nin==1 && nfp == 0 && nfn>0");
Long64_t nentries_0h = tree->GetEntries("nip ==0 && nfp == 0 && nfn==0");
Long64_t nentries_h = tree->GetEntries("nip ==0 && nfp > 0 && nfn>0");


cout << "initially 1n1pi" << nentries_1n1pi << endl;
cout<<endl;
cout << "then 1p" << nentries_1pi_1p << endl;
cout << "then 1p and pi- " << nentries_1pim_1p1pim << endl;
cout << "then 1p and pi+ " << nentries_1pim_1p1pip << endl;
cout << "then 1p and pi0: " << nentries_1pim_1p1pi0 << endl;
cout << "then 1p and many pions (pion production)" << nentries_1pim_1pnpi << endl;
cout << "then more than 1 proton (any no of pions, any =no of neutrons)" << nentries_1pim_np << endl;
cout << "then no protons, at least one neutron" << nentries_1pim_n << endl;
cout << "then no hadrons (maybe pions)" << nentries_0h << endl;
cout << "protons and neutrons and pions" << nentries_h << endl;

cout << endl;

  Long64_t nentries = tree->GetEntries("");
  Long64_t nentries_1p = tree->GetEntries("nfp == 1 && nfn == 0 && (nfpip + nfpim + nfpi0)==0");
 
  Long64_t nentries_1p1pip  = tree->GetEntries("nfp ==1 && nfpim==1 && (nfpip + nfpim + nfpi0)==1 && nfn==0");
  Long64_t nentries_1p1pi0  = tree->GetEntries("nfp ==1 && nfpi0==1 && (nfpip + nfpim + nfpi0)==1 && nfn==0");

    Long64_t nentries_np = tree->GetEntries("nfp > 1 && nfn==0 && (nfpip + nfpim + nfpi0)<2 ");

  
     Long64_t nentries_nn =  tree->GetEntries("nfn > 0 && nfp==0");

  Long64_t nentries_1pnpi = tree->GetEntries("nfn ==0 && nfp == 1 && (nfpip + nfpim + nfpi0)>1");
  Long64_t nentries_1nnpi = tree->GetEntries("nfp ==0 && nfn == 1 && (nfpip + nfpim + nfpi0)>1");


  Long64_t nentries_1p1n = tree->GetEntries("(nfpip + nfpim + nfpi0)==0 && nfn==1 && nfp==1");
  Long64_t nentries_1p1n_pi = tree->GetEntries("(nfpip + nfpim + nfpi0)>0 && nfn>=1 && nfp>=1");
 Long64_t nentries_nnnpi = tree->GetEntries("(nfpip + nfpim + nfpi0)>1 && nfn > 1");
   Long64_t nentries_npnpi = tree->GetEntries("(nfpip + nfpim + nfpi0)>1 && nfp > 1 && nfn==0");
     Long64_t nentries_nohadnopi = tree->GetEntries("(nfpip + nfpim + nfpi0)==0 && nfp ==0 && nfn == 0");


cout << "total entries" << nentries << endl;
cout<<endl;
cout << "1p" << nentries_1p << endl;
cout << "1p1pip" << nentries_1p1pip << endl;
cout << "1p1pim" << nentries_1p1pim << endl;
cout << "1p1pi0" << nentries_1p1pi0 << endl;


cout<<endl;
cout << "many protons, one or zero pions" << nentries_np << endl;
cout << "many protons, many pions(no neutrons)" << nentries_npnpi << endl;
cout << "1 proton, many pions" << nentries_1pnpi << endl;


cout<<endl;
cout << "neutrons, no protons (yes pions)" << nentries_nn << endl;

cout<<endl;
cout<<endl;

cout << "mec_1p1n (no pions)" << nentries_1p1n << endl;
cout << "neutrons and protons and pions" << nentries_1p1n_pi << endl;

cout << "no hadrons no pions" << nentries_nohadnopi << endl;



Long64_t nentries_res = tree->GetEntries("res==1");
Long64_t nentries_1p  = tree->GetEntries("res == 1 && nfp ==1 && (nfpip + nfpim + nfpi0)==0 && nfn==0");
Long64_t nentries_1n  = tree->GetEntries("res == 1 && nfn ==1 && (nfpip + nfpim + nfpi0)==0 && nfp==0");
Long64_t nentries_1p_1pi  = tree->GetEntries("res == 1 && nfp ==1 && (nfpip + nfpim + nfpi0)==1 && nfn==0");
Long64_t nentries_1n_1pi  = tree->GetEntries("res == 1 && nfp ==0 && (nfpip + nfpim + nfpi0)==1 && nfn==1");
Long64_t nentries_1p_npi  = tree->GetEntries("res == 1 && nfp ==1 && (nfpip + nfpim + nfpi0)>1 && nfn==0");
Long64_t nentries_1n_npi  = tree->GetEntries("res == 1 && nfp ==0 && (nfpip + nfpim + nfpi0)>1 && nfn==1");
Long64_t nentries_np_0n = tree->GetEntries("res == 1 && nfp > 1 && nfn==0");
Long64_t nentries_nn_0p = tree->GetEntries("res == 1 && nfp == 0 && nfn>1");
Long64_t nentries_0h = tree->GetEntries("res==1 && nfp==0 && nfn==0");
Long64_t nentries_nn_np = tree->GetEntries("res == 1 && nfp > 0 && nfn>0");

cout << "res" << nentries_res;
cout << endl;
cout << "1p" << nentries_1p << endl;
cout << "1n" << nentries_1n << endl;
cout << "1p_1pi" << nentries_1p_1pi << endl;
cout << "1n_1pi" << nentries_1n_1pi << endl;
cout << "1p_npi" << nentries_1p_npi << endl;
cout << "1n_npi" << nentries_1n_npi << endl;
cout << "np_0n" << nentries_np_0n << endl;
cout << "nn_0p" << nentries_nn_0p << endl;
cout << "nn_np" << nentries_nn_np << endl;
cout << "0h" << nentries_0h << endl;


Long64_t nentries_pip = tree->GetEntries("res == 1 && nfpip>0 && (nfpip + nfpim + nfpi0)==1");
Long64_t nentries_pim = tree->GetEntries("res == 1 && nfpim>0 && (nfpip + nfpim + nfpi0)==1");
Long64_t nentries_pi0 = tree->GetEntries("res == 1 && nfpi0>0 && (nfpip + nfpim + nfpi0)==1");

cout << "pip = " << nentries_pip << endl;
cout << "pim = " << nentries_pim << endl;
cout << "pi0 = " << nentries_pi0 << endl;

Long64_t nentries_dis = tree->GetEntries("dis == 1" );
Long64_t nentries_disn = tree->GetEntries("(dis == 1 && nfp>0) || (dis ==1 && nfn>0)");

cout << "dis = " << nentries_dis << endl;
cout << "dis npnn = " << nentries_disn << endl;



Long64_t bentries_0 = tree->GetEntries("(nipip+nipim+nipi0)==0" );
Long64_t bentries_pi0 = tree->GetEntries("(nipip+nipim+nipi0)==1 && nipi0==1" );
Long64_t bentries_pip = tree->GetEntries("(nipip+nipim+nipi0)==1 && nipip==1" );
Long64_t bentries_pim = tree->GetEntries("(nipip+nipim+nipi0)==1 && nipim==1" );
Long64_t bentries_plural = tree->GetEntries("(nipip+nipim+nipi0)>1" );

Long64_t aentries_0 = tree->GetEntries("(nfpip+nfpim+nfpi0)==0" );
Long64_t aentries_pi0 = tree->GetEntries("(nfpip+nfpim+nfpi0)==1 && nfpi0==1" );
Long64_t aentries_pip = tree->GetEntries("(nfpip+nfpim+nfpi0)==1 && nfpip==1" );
Long64_t aentries_pim = tree->GetEntries("(nfpip+nfpim+nfpi0)==1 && nfpim==1" );
Long64_t aentries_plural = tree->GetEntries("(nfpip+nfpim+nfpi0)>1" );


cout << "bentries_0" <<  bentries_0 << endl;
cout << "bentries_pi0" <<  bentries_pi0 << endl;
cout << "bentries_pip" <<  bentries_pip << endl;
cout << "bentries_pim" <<  bentries_pim << endl;
cout << "bentries_plural" <<  bentries_plural << endl;

cout << "aentries_0" <<  aentries_0 << endl;
cout << "aentries_pi0" <<  aentries_pi0 << endl;
cout << "aentries_pip" <<  aentries_pip << endl;
cout << "aentries_pim" <<  aentries_pim << endl;
cout << "aentries_plural" <<  aentries_plural << endl;
*/

Long64_t nentries_1p_1pic  = tree->GetEntries("res == 1 && nip ==1 && (nipip + nipim + nipi0)==1 && nin==0 && (nipip + nipim) == 1");
Long64_t nentries_1p_1pi0  = tree->GetEntries("res == 1 && nip ==1 && (nipip + nipim + nipi0)==1 && nin==0 && nipi0==1");

cout << "charged  " << nentries_1p_1pic << endl;
cout << "neutral  " << nentries_1p_1pi0 << endl;
}