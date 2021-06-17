#include "TH1.h"
#include "TH2.h"
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

#include "OscillationHelper.hxx"
R__LOAD_LIBRARY(WrappedProb3++.20121225.so)

using namespace std;

//PHYSICAL CONSTANTS      //can maybe one day put these in a header                                                        
const Double_t ELECTRON_MASS= 0.00051099; //Gev 
const Double_t MUON_MASS= 0.105658; //Gev     
const Double_t BEAM_ENERGY = 1.161; //GeV  //Change for different beam energies
const Double_t Eb = 0.025;//binding energies from Afro: 0.025 for 12C, 0.04 for 40Ar, 0.03 for 56Fe NEED TO CHECK VALIDITY                                                          
const Double_t PROTON_MASS =  0.938272; //mass of proton                                           
const Double_t NEUTRON_MASS =  0.939565; //mass of neutron                                        
const Double_t PION_MASS = 0.13958; //mass of charged pion 
const Double_t DELTA_MASS = 1.232; //mass of delta state                                                                 
const Double_t NEUTRAL_PION_MASS = 0.13498; //(GeV/c^2) mass of neutral pion                                                                                
const Double_t FINE_STRUCTURE_CONSTANT = 1./137.035999139; 

//GLOBAL VARIABLES FROM TREE - SEE WORD DOC
Double_t Ev; //incoming lepton energy
Double_t El; //Final state primary lepton energy
Double_t Ef; //energy of kth final state particle in hadronic system
Double_t W; // Hadronic invariant mass
Double_t x; //Bjorken x (as computed from the event record)
Double_t Q2; //Momentum transfer Q^{2} (as computed from the event record)
Bool_t res, qel, dis, mec; //is it a RES, QEL, DIS, MEC event?
//Double_t pxf, pyf, pzf; //final momenta in x, y, z of kth final state particle
Double_t pxl, pyl, pzl; //final state primary lepton momenta in x, y, z
//Final state particle variables. Final state = after intranuclear rescattering
int pdg; //pdg code of kth final state particle in hadronic system 
int nfn; //number of final state neutron and antineutron
int nfp; //number of final state proton and antiproton   
int nfpip; //number of final state pi+
int nfpim; //number of final state pi-
int nfpi0; //number of final state pi0
int nin; //number of initial state neutron 
int nip; //number of initial state proton 
int nipip; //number of inital state pi+                                                                                                   
int nipim; //number of initial state pi-                                                                                                  
int nipi0; //number of initial state pi0  
int nfkp;//number of final state K+
int nfkm;//number of final state K-     
int nfk0;//number of final state K0    
int nfother;//number of final state other heavier hadrons    
Int_t nf; //number of final state particles in hadronic system
int nfem; //number of final state e+, e-, photon
int nikp; //number of inital state K+                                                            
int nikm; //number of initial state K-                                                            
int nik0; //number of initial state K0  
int tgt; //pdg code of nuclear target
int hitnuc; //hit nucleon pdg code
bool cc; //is it a charged current event?
bool nc;
int resid;
int em;
int nres_events;
int n1in1ipi0_events;
int n1in1ipip_events;
int n1ip1ipi0_events;
int n1ip1ipip_events;
int n1ip1ipim_events;
int n1ip3ipi0_events;
int n1in1ipim1ipip1ipi0_events;
int n1in1ipim1ipip_events;
int nproton;
int nneutron;

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Double_t squared(Double_t num) {
  return TMath::Power(num, 2);
}

//CHECK WE ARE DOING ELECTRON SCATTERING, NOT NEUTRINO
bool CheckIfElectrons(TTree *tree)//
{
  bool isElectron=false;
  tree->SetBranchStatus("em", 1);
  tree->SetBranchAddress("em", &isElectron); // If this is true (electromagnetic scattering) then it isn't a neutrino
  tree->GetEntry(); // The first entry is fine, this should tell us all we need to know. Note that it would be misleading if we had a mixed neutrino and electron beam... but we don't.
  return isElectron;
}

//WHICH TARGET?
string GetTargetString(TTree *tree)
{
  int nProtons;
  string targetString="";
  tree->SetBranchStatus("Z", 1); 
  tree->SetBranchAddress("Z", &nProtons);
  tree->GetEntry();
  switch (nProtons)
  {
    case 18: targetString="argon"; break;
    case 6: targetString="carbon"; break;
    case 26: targetString="iron"; break;
    case 1: targetString="hydrogen"; break;
    case 2: targetString="helium"; break;  
    default: "unknown target";
  }
  return targetString;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//MAIN FUNCTION
int func (string inFileName, string outdir, string tester){ //not taking in any complicated parameters yet

  cout << "ROOT file being analysed:" << inFileName << endl;

//READ THE ROOT FILE AND LOAD THE GST TREE
TFile *f = new TFile(inFileName.c_str());
  f->ls();
  cout << (string(inFileName)).c_str() << endl;
  TTree *tree = (TTree*) f->Get("gst");

//NUMBER OF ENTRIES
Long64_t nentries = tree->GetEntries();
 if (tester == "test"){
     nentries = 10000;} //for debugging

  //WHICH LEPTON?
  bool isElectronMode = CheckIfElectrons(tree);
    
  //WHICH TARGET?
  string targetString = GetTargetString(tree);

  // PRINT OUT WHAT SORT OF SCATTERING IT IS
  cout<<endl;
  cout<<(isElectronMode?"Electron":"Neutrino")<< " scattering on "<<targetString<<endl;
  cout<<"Number of entries:" << nentries << endl;
  cout<<endl;



  //LOOP THROUGH ENTRIES TO FIND MAXIMUM
  tree -> SetBranchStatus("nf", 1);
  tree -> SetBranchAddress("nf", &nf);
  int maxNF = 0;
  for(Long64_t j=0;j<nentries;j++){
    tree->GetEntry(j);
    if(nf>maxNF){
      maxNF = nf;
    }
  }
  
  Int_t pdgf[maxNF];
  Double_t pxf[maxNF], pyf[maxNF], pzf[maxNF];
  Double_t Ef[maxNF];

  tree -> SetBranchStatus("*", 0); //sets all the branches off
  tree -> SetBranchStatus("Ev", 1); //turns on the one we care about
  tree -> SetBranchStatus("res", 1);
  tree -> SetBranchStatus("qel", 1);
  tree -> SetBranchStatus("dis", 1);
  tree -> SetBranchStatus("mec", 1);
  tree -> SetBranchStatus("El", 1);
  tree -> SetBranchStatus("pxl", 1);
  tree -> SetBranchStatus("pyl", 1);
  tree -> SetBranchStatus("pzl", 1);
  tree -> SetBranchStatus("Ef", 1);
  tree -> SetBranchStatus("pdgf", 1);
  tree -> SetBranchStatus("pxf", 1);
  tree -> SetBranchStatus("pyf", 1);
  tree -> SetBranchStatus("pzf", 1);
  tree -> SetBranchStatus("Q2", 1);
  tree -> SetBranchStatus("W", 1);
  tree -> SetBranchStatus("x", 1);
  tree -> SetBranchStatus("nfn", 1);
  tree -> SetBranchStatus("nfp", 1);
  tree -> SetBranchStatus("nfpip", 1);
  tree -> SetBranchStatus("nfpim", 1);
  tree -> SetBranchStatus("nfpi0", 1);
  tree -> SetBranchStatus("nfkp", 1);
  tree -> SetBranchStatus("nfkm", 1);
  tree -> SetBranchStatus("nfk0", 1); 
  tree -> SetBranchStatus("nfem", 1);
  tree -> SetBranchStatus("nip", 1);
  tree -> SetBranchStatus("nipip", 1);
  tree -> SetBranchStatus("nipim", 1);
  tree -> SetBranchStatus("nipi0", 1);
  tree -> SetBranchStatus("nin", 1);
  tree -> SetBranchStatus("nikp", 1);
  tree -> SetBranchStatus("nikm", 1);
  tree -> SetBranchStatus("nik0", 1);
  tree -> SetBranchStatus("tgt", 1);
  tree -> SetBranchStatus("hitnuc", 1);
  tree -> SetBranchStatus("cc", 1);
  tree -> SetBranchStatus("nc", 1);
  tree -> SetBranchStatus("resid",1);
  tree -> SetBranchStatus("em", 1);


  tree -> SetBranchAddress("Ev", &Ev); //puts the variable referenced into the branch called "xxx", they have the same name to make it less confusing
  tree -> SetBranchAddress("res", &res);
  tree -> SetBranchAddress("qel", &qel);
  tree -> SetBranchAddress("dis", &dis);
  tree -> SetBranchAddress("mec", &mec);
  tree -> SetBranchAddress("El", &El);  
  tree -> SetBranchAddress("pxl", &pxl);  
  tree -> SetBranchAddress("pyl", &pyl);  
  tree -> SetBranchAddress("pzl", &pzl);
  tree -> SetBranchAddress("Ef", &Ef);
  tree -> SetBranchAddress("pdgf", &pdgf);
  tree -> SetBranchAddress("pxf", &pxf);
  tree -> SetBranchAddress("pyf", &pyf);
  tree -> SetBranchAddress("pzf", &pzf);
  tree -> SetBranchAddress("Q2", &Q2);
  tree -> SetBranchAddress("W", &W);
  tree -> SetBranchAddress("x", &x);
  tree -> SetBranchAddress("nfn", &nfn);
  tree -> SetBranchAddress("nfp", &nfp);
  tree -> SetBranchAddress("nfpip", &nfpip);
  tree -> SetBranchAddress("nfpim", &nfpim);
  tree -> SetBranchAddress("nfpi0", &nfpi0);
  tree -> SetBranchAddress("nfkp", &nfkp);
  tree -> SetBranchAddress("nfkm", &nfkm);
  tree -> SetBranchAddress("nfk0", &nfk0);
  tree -> SetBranchAddress("nfem", &nfem);
  tree -> SetBranchAddress("nip", &nip);
  tree -> SetBranchAddress("nipip", &nipip);
  tree -> SetBranchAddress("nipim", &nipim);
  tree -> SetBranchAddress("nipi0", &nipi0);
  tree -> SetBranchAddress("nin", &nin);
  tree -> SetBranchAddress("nikm", &nikm);
  tree -> SetBranchAddress("nikp", &nikp);
  tree -> SetBranchAddress("nik0", &nik0);
  tree -> SetBranchAddress("tgt", &tgt);
  tree -> SetBranchAddress("hitnuc", &hitnuc);
  tree -> SetBranchAddress("cc", &cc);
  tree -> SetBranchAddress("nc", &nc);
  tree -> SetBranchAddress("resid" , &resid);
  int Bins = 280;
  double Histo_xmax = 3;
  
  //making histograms
  TH1D *hEp = new TH1D("Ep", "Proton energies", Bins, 0.0, Histo_xmax); //create a pointer to a histogram
  TH1D *hEp_primary = new TH1D("Ep_primary", "Energy of primary proton", Bins, 0.0, Histo_xmax);
  //TH1D *hEp_other = new TH1D("Ep_other", "Energy distribution of secondary protons", Bins, 0.0, Histo_xmax);
  TH1D *hEp_total = new TH1D("Ep_total", "Total proton energy", Bins, 0.0, Histo_xmax);
  TH2D *hEp_min_max_qe = new TH2D ("Ep_min_max_qe,", "comparing primary and secondary proton energies in qe events", Bins, 0.0, Histo_xmax, Bins, 0.0, Histo_xmax);
  TH2D *hEp_min_max_res = new TH2D ("Ep_min_max_res,", "comparing primary and secondary proton energies in resonant events", Bins, 0.0, Histo_xmax, Bins, 0.0, Histo_xmax);
  TH2D *hEp_min_max_mec = new TH2D ("Ep_min_max_mec,", "comparing primary and secondary proton energies in mec events", Bins, 0.0, Histo_xmax, Bins, 0.0, Histo_xmax);
  TH2D *hEp_min_max_dis = new TH2D ("Ep_min_max_dis,", "comparing primary and secondary proton energies in dis events", Bins, 0.0, Histo_xmax, Bins, 0.0, Histo_xmax);


  gRandom = new TRandom3();
  gRandom->SetSeed(10);
  
  for( Long64_t i = 0; i < nentries ; i++){
  //cout<<"Iteration: "<< i << endl;
  //cout<<"No. of final state particles: " << nf << endl;
  
      tree -> GetEntry(i);
//Initialise proton information
      double pP = 0;      //proton momentum
      double protonK = 0; //proton kinetic energy
      double pCos = 0;    //cosine of scattering angle 
      double pPhi = 0;    //phi
      int pCount = 0;
      double Ep_max = 0; //energy of the proton with max energy
      double Ep_min=999999999;
      double tot = 0; //summing up proton energies

    for(Int_t k=0; k<nf; k++) { //loop through k final particles
      int pdg = pdgf[k];
      TVector3 kVec(pxf[k], pyf[k], pzf[k]);//momentum vector for kth particle
      double kP = kVec.Mag();

//SETTING VALS IF PROTON
      if (pdg == 2212) {
	    if (kP > pP) {
	    pP = kP;
	    pCos = kVec.CosTheta(); 
	    pPhi = kVec.Phi() + TMath::Pi(); 
	    protonK = Ef[k] - PROTON_MASS;
      pCount +=1;
      tot += Ef[k];
      if (Ef[k] > Ep_max){
      Ep_max = Ef[k];
      };

      if (Ef[k] < Ep_min){
      Ep_min = Ef[k];
      }; 
 
	};
  }

};   // end loop iterating over final state particles
 
if (pCount == 2){
  if (Ep_max!=0)
    {
      cout << Ep_max << endl;
      cout << Ep_min << endl;
      cout << endl;
      hEp -> Fill(Ep_max, 1.0);
      hEp_primary -> Fill(Ep_max, 1.0);
      if (qel == true){
      hEp_min_max_qe -> Fill((Ep_max-0.938), (Ep_min-0.938));
  

      }
      if (res == true){
      hEp_min_max_res -> Fill((Ep_max-0.938), (Ep_min-0.938));

      }
      if (dis == true){
      hEp_min_max_mec -> Fill((Ep_max-0.938), (Ep_min-0.938));

      }
      if (mec == true){
      hEp_min_max_dis -> Fill((Ep_max-0.938), (Ep_min-0.938));
 

      }
      
      }
            
}
  }; //end loop iterating over entries
  
    delete gRandom;
  
  //string targetString = GetTargetString(tree);
  string  PATH = "nanPlots/" + outdir;
  string filename_nopath = inFileName.substr(inFileName.find_last_of("/")+1); //gets filename from filepath
  string outFileName = PATH+string(filename_nopath);
  TFile *output = new TFile(outFileName.c_str(),"RECREATE"); //makes the file writeable, only 1 file can be open and writeable at a time in ROOT
    hEp -> Write();
    hEp_primary -> Write();
    hEp_min_max_qe -> Write();
    hEp_min_max_res -> Write();
    hEp_min_max_mec -> Write();
    hEp_min_max_dis -> Write();
   
    
    output -> Close();
 
 cout << "Histograms filled. Output closed!" << endl;
 return 0;

}
