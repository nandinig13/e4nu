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

using namespace std;
void entries(){

  string GST = "/pnfs/genie/persistent/users/apapadop/e4v_SuSav2/Exclusive/electrons/C12_1161GeV/apapadop_SuSav2_C12_1161GeV_master.root";

TFile *f = new TFile(GST.c_str());

TTree *tree = (TTree*) f->Get("gst");
//tree->Scan("nip:nin:nipip:nipim:nipi0", "res && nip ==2");
Long64_t nentries = tree->GetEntries();
Long64_t nentries_qel = tree->GetEntries("qel");
Long64_t nentries_res = tree->GetEntries("res");
Long64_t nentries_dis = tree->GetEntries("dis");

Long64_t nopiTOnopi = tree->GetEntries("nipip==0 && nipi0==0 && nipim==0 && nfpip==0 && nfpi0==0 && nfpim==0");
Long64_t nopiTO1pi0 = tree->GetEntries("nipip==0 && nipi0==0 && nipim==0 && nfpip==0 && nfpi0==1 && nfpim==0");
Long64_t nopiTO1pip = tree->GetEntries("nipip==0 && nipi0==0 && nipim==0 && nfpip==1 && nfpi0==0 && nfpim==0");
Long64_t nopiTO1pim = tree->GetEntries("nipip==0 && nipi0==0 && nipim==0 && nfpip==0 && nfpi0==0 && nfpim==1");
Long64_t nopiTO2pluspions = tree->GetEntries("nipip==0 && nipi0==0 && nipim==0 && (nfpip + nfpi0 +nfpim)>=2");

Long64_t pi0TOnopi = tree->GetEntries("nipip==0 && nipi0==1 && nipim==0 && nfpip==0 && nfpi0==0 && nfpim==0");
Long64_t pi0TO1pi0 = tree->GetEntries("nipip==0 && nipi0==1 && nipim==0 && nfpip==0 && nfpi0==1 && nfpim==0");
Long64_t pi0TO1pip = tree->GetEntries("nipip==0 && nipi0==1 && nipim==0 && nfpip==1 && nfpi0==0 && nfpim==0");
Long64_t pi0TO1pim = tree->GetEntries("nipip==0 && nipi0==1 && nipim==0 && nfpip==0 && nfpi0==0 && nfpim==1");
Long64_t pi0TO2pluspions = tree->GetEntries("nipip==0 && nipi0==1 && nipim==0 && (nfpip + nfpi0 +nfpim)>=2");

Long64_t pipTOnopi = tree->GetEntries("nipip==1 && nipi0==0 && nipim==0 && nfpip==0 && nfpi0==0 && nfpim==0");
Long64_t pipTO1pi0 = tree->GetEntries("nipip==1 && nipi0==0 && nipim==0 && nfpip==0 && nfpi0==1 && nfpim==0");
Long64_t pipTO1pip = tree->GetEntries("nipip==1 && nipi0==0 && nipim==0 && nfpip==1 && nfpi0==0 && nfpim==0");
Long64_t pipTO1pim = tree->GetEntries("nipip==1 && nipi0==0 && nipim==0 && nfpip==0 && nfpi0==0 && nfpim==1");
Long64_t pipTO2pluspions = tree->GetEntries("nipip==1 && nipi0==0 && nipim==0 && (nfpip + nfpi0 +nfpim)>=2");

 Long64_t pimTOnopi = tree->GetEntries("nipip==0 && nipi0==0 && nipim==1 && nfpip==0 && nfpi0==0 && nfpim==0");
 Long64_t pimTO1pi0 = tree->GetEntries("nipip==0 && nipi0==0 && nipim==1 && nfpip==0 && nfpi0==1 && nfpim==0");
 Long64_t pimTO1pip = tree->GetEntries("nipip==0 && nipi0==0 && nipim==1 && nfpip==1 && nfpi0==0 && nfpim==0");
 Long64_t pimTO1pim = tree->GetEntries("nipip==0 && nipi0==0 && nipim==1 && nfpip==0 && nfpi0==0 && nfpim==1");
 Long64_t pimTO2pluspions = tree->GetEntries("nipip==0 && nipi0==0 && nipim==1 && (nfpip + nfpi0 +nfpim)>=2");

 Long64_t TWOpluspionsTOnopi = tree->GetEntries("(nipip + nipi0 +nipim)>=2 && nfpip==0 && nfpi0==0 && nfpim==0");
 Long64_t TWOpluspionsTO1pi0 = tree->GetEntries("(nipip + nipi0 +nipim)>=2 && nfpip==0 && nfpi0==1 && nfpim==0");
 Long64_t TWOpluspionsTO1pip = tree->GetEntries("(nipip + nipi0 +nipim)>=2 && nfpip==1 && nfpi0==0 && nfpim==0");
 Long64_t TWOpluspionsTO1pim = tree->GetEntries("(nipip + nipi0 +nipim)>=2 && nfpip==0 && nfpi0==0 && nfpim==1");
 Long64_t TWOpluspionsTO2pluspions = tree->GetEntries("(nipip + nipi0 +nipim)>=2 && (nfpip + nfpi0 +nfpim)>=2");

cout << GST << endl;
cout << "number of entries" << nentries <<endl;
cout << "number of QE events" << nentries_qel <<endl;
cout << "number of RES events" << nentries_res <<endl;
cout << "number of DIS events" << nentries_dis <<endl;

 cout << "0pi to 0pi pi0 pi+ pi- 2+pi's " << nopiTOnopi << " & " <<nopiTO1pi0 << " & " <<nopiTO1pip<<" & "<<nopiTO1pim<<" & "<<nopiTO2pluspions<<endl;
 cout << "1pi0 to 0pi pi0 pi+ pi- 2+pi's " << pi0TOnopi << " & " <<pi0TO1pi0 << " & " <<pi0TO1pip<<" & "<<pi0TO1pim<<" & "<<pi0TO2pluspions<<endl;
 cout << "1pi+ to 0pi pi0 pi+ pi- 2+pi's " << pipTOnopi << " & " <<pipTO1pi0 << " & " <<pipTO1pip<<" & "<<pipTO1pim<<" & "<<pipTO2pluspions<<endl;
 cout << "1pi- to 0pi pi0 pi+ pi- 2+pi's " << pimTOnopi << " & " <<pimTO1pi0 << " & " <<pimTO1pip<<" & "<<pimTO1pim<<" & "<<pimTO2pluspions<<endl;
 cout << "2+pi's to 0pi pi0 pi+ pi- 2+pi's " << TWOpluspionsTOnopi << " & " <<TWOpluspionsTO1pi0 << " & " <<TWOpluspionsTO1pip<<" & "<<TWOpluspionsTO1pim<<" & "<<TWOpluspionsTO2pluspions<<endl;
}