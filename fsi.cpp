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
#include "Constants.h"

#include "OscillationHelper.hxx"
R__LOAD_LIBRARY(WrappedProb3++.20121225.so)

using namespace std;


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
Double_t cthf;
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
enum CutID //enumerated type
{
  Q2_CUT, //Q2 > Q2_cutoff
  W_CUT, //W<WCUTOFF
  QEL_CUT, //Quasi-elastic
  BJX_CUT, 	//Bjorken x factor abs(x-1) < bjx_cutoff
  FINAL_STATE_CHARGED_PION_CUT, // final state charged pion cut
  FINAL_STATE_PROTON_CUT, //  final state proton cut
  FINAL_STATE_NEUTRON_CUT, // final state neutron cut
  FINAL_STATE_NEUTRAL_PION_CUT, // final state neutral pion cut
};

Double_t squared(Double_t num) {
  return TMath::Power(num, 2);
}


//CUT CLASS
class Cut {
private:
  string cutString_;
  string titleText_;
  CutID cutID_; //make enum for CutID where we add cuts
public:
  Cut(string _titleString, string _cutString, CutID _cutID);
  // We can add other functions here, maybe you want the efficiency of the cut for example...
  string GetTitleText();
  string GetCutString();
  bool PassesCut();//returns true or false depending on if cut id passes e.g. TRUE_RES_CUT=true then give res...
};

Cut::Cut(string _titleString, string _cutString, CutID _cutID)
{
  this->titleText_=_titleString;
  this->cutString_=_cutString;
  this->cutID_=_cutID;
}

string Cut::GetTitleText()
{
  return this->titleText_;
}

string Cut::GetCutString()
{
  return this->cutString_;
}

bool Cut::PassesCut()
{
  //bool _cutID = true; 
  switch (cutID_)
  {
    case FINAL_STATE_CHARGED_PION_CUT: return  (nfpim + nfpip > 0); 
         break; // 1 final state pion
    case FINAL_STATE_NEUTRAL_PION_CUT: return (nfpi0 ==1);
       break;
    case FINAL_STATE_PROTON_CUT: return (nfp >= 1);  //changed from ==1 because too strict?
       break; // 1 final state proton
    case FINAL_STATE_NEUTRON_CUT: return (nfn == 0);
       break;
    default: return true;
  }
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

	


int MCDecide(double prob){
  //TRandom *r3 = new TRandom3();
  //r3->SetSeed(0);
 // Double_t  num= r3->Rndm();
 auto num = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
 //cout << "The random number generated was " << num << endl;
  //double num = 0.5;
  if (num < (prob)) {
 //cout << "It was accepted" << endl;
    return 1;
  }
  else {
 cout << "It was rejected" << endl;
    return 0;
  }
  //delete r3;
  }

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


//MAIN FUNCTION
int func (string inFileName, string outdir, string tester,  string acceptance){ //not taking in any complicated parameters yet

const double_t BEAM_ENERGY = 2.261;// CHANGE AS REQUIRED

//RESOLUTIONS FOR SMEARING
    double reso_p = 0.01; // smearing for the proton
    double reso_e = 0.005; // smearing for the electron
    double reso_pi = 0.007; //smearing for pions

    if(BEAM_ENERGY == 1.161) {  //justification for this in Anjali's paper - lower energy has more error
			reso_p = 3*reso_p; reso_e = 3*reso_e; reso_pi = 3*reso_pi;
		}

  //FILEPATH = filepath ... maybe try and cout without the whole path to look nicer
  cout << "ROOT file being analysed:" << inFileName << endl;

//Could put couts here of which cuts have been specified by the command line arguments

//READ THE ROOT FILE AND LOAD THE GST TREE
TFile *f = new TFile(inFileName.c_str());
  f->ls();
  cout << (string(inFileName)).c_str() << endl;
  TTree *tree = (TTree*) f->Get("gst");

//NUMBER OF ENTRIES
Long64_t nentries = tree->GetEntries();
 if (tester == "test"){
     nentries = 7000;} //for debugging

  //WHICH LEPTON?
  bool isElectronMode = CheckIfElectrons(tree);
    
  //WHICH TARGET?
  string targetString = GetTargetString(tree);

  // PRINT OUT WHAT SORT OF SCATTERING IT IS
  cout<<endl;
  cout<<(isElectronMode?"Electron":"Neutrino")<< " scattering on "<<targetString<<endl;
  cout<<"Number of entries:" << nentries << endl;
  cout<<endl;

  //SET BINDING ENERGY ACCORDINGLY
  double_t Eb = 0;
  if (targetString == "carbon"){
    Eb = C12_bind_en;
  }
  if (targetString == "iron"){
    Eb = He4_bind_en;
  }
  if (targetString == "helium"){
    Eb = Fe_bind_en;
  }
  if (Eb == 0){
    cout << "Binding energy not set!" << endl;
  }
  
//CUTS
  vector<Cut*> cuts;
    //cuts.push_back(new Cut("Bjorken x cut","TMath::Abs(x-1) < 0.2")); 
    //cuts.push_back(new Cut("1p" , "1p", FINAL_STATE_PROTON_CUT));
    //cuts.push_back(new Cut(">= 1p" , ">= 1p", FINAL_STATE_PROTON_CUT));
    //cuts.push_back(new Cut(">1pi" , ">1pi", FINAL_STATE_CHARGED_PION_CUT)); //make sure 
   // cuts.push_back(new Cut("0n" , "0n", FINAL_STATE_NEUTRON_CUT));
   // cuts.push_back(new Cut("0pi0", "0pi0", FINAL_STATE_NEUTRAL_PION_CUT));

  string cutText="Cuts: ";
  for (int i=0; i<cuts.size();i++) // Loops all your cuts
  {
    if (i>0)cutText+=" and ";
    cutText+=cuts.at(i)->GetTitleText();
  }
  cout<<cutText<<endl;

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
  Double_t cthf[maxNF];

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
  tree -> SetBranchStatus("cthf", 1);


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
  tree -> SetBranchAddress("cthf", &cthf);

  int Bins = 300;
  double Histo_xmax = 7;
  
  TH1D *hLepP = new TH1D("LepP", "Electron momentum for accepted events", Bins, 0.0, Histo_xmax);
  TH1D *hEprot_pass = new TH1D("Eprot_pass", "Energies of protons that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEprot_fail = new TH1D("Eprot_fail", "Energies of protons that failed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEpim_pass = new TH1D("Epim_pass", "Energies of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEpim_fail = new TH1D("Epim_fail", "Energies of pion that failed acceptance", Bins, 0.0, Histo_xmax);
   TH1D *hEpip_pass = new TH1D("Epip_pass", "Energies of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEpip_fail = new TH1D("Epip_fail", "Energies of pion that failed acceptance", Bins, 0.0, Histo_xmax);
  
  TH1D *hEprot_pass_res = new TH1D("Eprot_pass_res", "Energies of protons that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEprot_fail_res = new TH1D("Eprot_fail_res", "Energies of protons that failed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEpim_pass_res = new TH1D("Epim_pass_res", "Energies of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEpim_fail_res = new TH1D("Epim_fail_res", "Energies of pion that failed acceptance", Bins, 0.0, Histo_xmax);
   TH1D *hEpip_pass_res = new TH1D("Epip_pass_res", "Energies of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEpip_fail_res = new TH1D("Epip_fail_res", "Energies of pion that failed acceptance", Bins, 0.0, Histo_xmax);
   
    TH1D *hEprot_pass_dis = new TH1D("Eprot_pass_dis", "Energies of protons that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEprot_fail_dis = new TH1D("Eprot_fail_dis", "Energies of protons that failed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEpim_pass_dis = new TH1D("Epim_pass_dis", "Energies of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEpim_fail_dis = new TH1D("Epim_fail_dis", "Energies of pion that failed acceptance", Bins, 0.0, Histo_xmax);
   TH1D *hEpip_pass_dis = new TH1D("Epip_pass_dis", "Energies of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hEpip_fail_dis = new TH1D("Epip_fail_dis", "Energies of pion that failed acceptance", Bins, 0.0, Histo_xmax);

    TH1D *hMomprot_pass = new TH1D("Momprot_pass", "Momenta of protons that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMomprot_fail = new TH1D("Momprot_fail", "Momenta of protons that failed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMompip_pass = new TH1D("Mompip_pass", "Momenta of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMompip_fail = new TH1D("Mompip_fail", "Momenta of pion that failed acceptance", Bins, 0.0, Histo_xmax);
   TH1D *hMompim_pass = new TH1D("Mompim_pass", "Momenta of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMompim_fail = new TH1D("Mompim_fail", "Momenta of pion that failed acceptance", Bins, 0.0, Histo_xmax);

   TH1D *hMomprot_pass_res = new TH1D("Momprot_pass_res", "Momenta of protons that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMomprot_fail_res = new TH1D("Momprot_fail_res", "Momenta of protons that failed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMompip_pass_res = new TH1D("Mompip_pass_res", "Momenta of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMompip_fail_res = new TH1D("Mompip_fail_res", "Momenta of pion that failed acceptance", Bins, 0.0, Histo_xmax);
   TH1D *hMompim_pass_res = new TH1D("Mompim_pass_res", "Momenta of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMompim_fail_res = new TH1D("Mompim_fail_res", "Momenta of pion that failed acceptance", Bins, 0.0, Histo_xmax);

     TH1D *hMomprot_pass_dis = new TH1D("Momprot_pass_dis", "Momenta of protons that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMomprot_fail_dis = new TH1D("Momprot_fail_dis", "Momenta of protons that failed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMompip_pass_dis = new TH1D("Mompip_pass_dis", "Momenta of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMompip_fail_dis = new TH1D("Mompip_fail_dis", "Momenta of pion that failed acceptance", Bins, 0.0, Histo_xmax);
   TH1D *hMompim_pass_dis = new TH1D("Mompim_pass_dis", "Momenta of pions that passed acceptance", Bins, 0.0, Histo_xmax);
  TH1D *hMompim_fail_dis = new TH1D("Mompim_fail_dis", "Momenta of pion that failed acceptance", Bins, 0.0, Histo_xmax);


TH1D *hctheta_prot = new TH1D("ctheta_prot", "Cos theta of proton", Bins, 0.0, Histo_xmax);
TH1D *hctheta_pip = new TH1D("ctheta_pip", "Cos theta of pip", Bins, 0.0, Histo_xmax);
TH1D *hctheta_pim = new TH1D("ctheta_pim", "Cos theta of pim", Bins, 0.0, Histo_xmax);
TH1D *hctheta_prot_pass = new TH1D("ctheta_prot_pass", "Cos theta of proton", Bins, 0.0, Histo_xmax);
TH1D *hctheta_pip_pass = new TH1D("ctheta_pip_pass", "Cos theta of pip", Bins, 0.0, Histo_xmax);
TH1D *hctheta_pim_pass = new TH1D("ctheta_pim_pass", "Cos theta of pim", Bins, 0.0, Histo_xmax);

   TH2D *hAng_pip_pass = new TH2D ("Ang_pip_pass", "Angles of pions that were rejected", Bins, 0.0, Histo_xmax, Bins, 0.0, Histo_xmax);
   TH2D *hAng_pim_pass = new TH2D ("Ang_pim_pass", "Angles of pions that were rejected", Bins, 0.0, Histo_xmax, Bins, 0.0, Histo_xmax);
  TH2D *hAng_prot_pass = new TH2D ("Ang_prot_pass", "Angles of protons that were rejected", Bins, 0.0, Histo_xmax, Bins, 0.0, Histo_xmax);
  TH2D *hAng_pip_fail = new TH2D ("Ang_pip_fail", "Angles of pions that were rejected", Bins, 0.0, Histo_xmax, Bins, 0.0, Histo_xmax);
   TH2D *hAng_pim_fail = new TH2D ("Ang_pim_fail", "Angles of pions that were rejected", Bins, 0.0, Histo_xmax, Bins, 0.0, Histo_xmax);
  TH2D *hAng_prot_fail = new TH2D ("Ang_prot_fail", "Angles of protons that were rejected", Bins, 0.0, Histo_xmax, Bins, 0.0, Histo_xmax);

TH1D *hWres_pass = new TH1D("Wres_pass", "Invariant mass for resonant events", Bins, 0.0, Histo_xmax);
TH1D *hWdis_pass = new TH1D("Wdis_pass", "Invariant mass for deep-inelastic events", Bins, 0.0, Histo_xmax);
  // --------- CLAS ACCEPTANCE MAPS -----------------
//TFile* file_acceptance_1_161 = TFile::Open("maps/e2a_maps_12C_E_1_161.root");
//TFile* file_acceptance_1_161_p = TFile::Open("maps/e2a_maps_12C_E_1_161_p.root");
//TFile* file_acceptance_1_161_pip = TFile::Open("maps/e2a_maps_12C_E_1_161_pip.root");
//TFile* file_acceptance_1_161_pim = TFile::Open("maps/e2a_maps_12C_E_1_161_pim.root");

TFile* file_acceptance_2_261 = TFile::Open("maps/e2a_maps_12C_E_2_261.root");
TFile* file_acceptance_2_261_p = TFile::Open("maps/e2a_maps_12C_E_2_261_p.root");
TFile* file_acceptance_2_261_pip = TFile::Open("maps/e2a_maps_12C_E_2_261_pip.root");
TFile* file_acceptance_2_261_pim = TFile::Open("maps/e2a_maps_12C_E_2_261_pim.root"); 
 
  //GAUSSIAN SMEARING

  gRandom = new TRandom3();
  gRandom->SetSeed(10);

  //SIGNAL EVENT COUNTER

int signal = 0;
int rescount = 0;
int true1p1pi = 0;
int qetrue1p1pi = 0;
int restrue1p1pi = 0;
int mectrue1p1pi = 0;
int distrue1p1pi = 0;


int false1p1pi = 0;
int qefalse1p1pi = 0;
int resfalse1p1pi = 0;
int disfalse1p1pi = 0;
int mecfalse1p1pi = 0;

  for( Long64_t i = 0; i < nentries ; i++)
  {
     tree -> GetEntry(i);

      vector<double> weights(maxNF, 1.0);
      vector<double> momenta(maxNF, 0.0);
      vector<double> smeared_energies(maxNF, 0.0);
      vector<double> smeared_momenta(maxNF, 0.0);
      vector<double> kinetic(maxNF, 0.0);
      vector<double> phi(maxNF, 0.0);
      //cout << "weight vector size = " << weights.size() << endl;
 
      TVector3 lP(pxl, pyl, pzl);   //get lepton information
      double leptonP = lP.Mag();
      double leptonCos = lP.CosTheta();
      double leptonPhi = lP.Phi() + TMath::Pi();
    

    //WHICH LEPTON? SETTING MASS
    double m_l; //lepton mass = muon mass for neutrino scat, lepton mass = electron mass for electron scat 
    if (isElectronMode == true){
      m_l = ELECTRON_MASS;
    }
    else {
      m_l = MUON_MASS;
      cout<<"Not in electron mode!"<<endl;
    }

    for(Int_t k=0; k<nf; k++) { //loop through k final particles
      TVector3 kVec(pxf[k], pyf[k], pzf[k]);//momentum vector for kth particle
      phi[k] = ((kVec.Phi()+TMath::Pi()));


//SMEARING ENERGY AND MOMENTA

 double smear_pipP = 0;
 double smear_pimP = 0;
 double smear_Epi = 0 ;

 double smear_lP = gRandom->Gaus(leptonP,reso_e*leptonP);
 double smear_El = TMath::Sqrt(m_l*m_l+smear_lP*smear_lP);

//ignoring neutrons and pi0
  if (pdgf[k] == 2212) { // proton
    kinetic[k] = Ef[k] - PROTON_MASS;
    momenta[k] = kVec.Mag();
    smeared_momenta[k] = gRandom->Gaus(momenta[k],reso_p*momenta[k]);
    smeared_energies[k] = TMath::Sqrt(PROTON_MASS*PROTON_MASS+smeared_momenta[k]*smeared_momenta[k]);  
    weights[k] = acceptance_c(momenta[k], cthf[k], phi[k], 2212, file_acceptance_2_261_p);
    hctheta_prot -> Fill((acos(cthf[k]) * 180/M_PI));
	}

 if (pdgf[k] == 211) { //pi+
    kinetic[k] = Ef[k] - PION_MASS;
    momenta[k] = kVec.Mag();
    smeared_momenta[k] = gRandom->Gaus(momenta[k],reso_pi*momenta[k]);
    smeared_energies[k] = TMath::Sqrt(PION_MASS*PION_MASS+smeared_momenta[k]*smeared_momenta[k]);
    weights[k] = acceptance_c(momenta[k], cthf[k], phi[k], 211, file_acceptance_2_261_pip);
    hctheta_pip -> Fill((acos(cthf[k]) * 180/M_PI));
	}

  if (pdgf[k] == -211) { //pi-
    kinetic[k] = Ef[k] - PION_MASS;
    momenta[k] = kVec.Mag();
    smeared_momenta[k] = gRandom->Gaus(momenta[k],reso_pi*momenta[k]);
    smeared_energies[k] = TMath::Sqrt(PION_MASS*PION_MASS+smeared_momenta[k]*smeared_momenta[k]); 
    weights[k] = acceptance_c(momenta[k], cthf[k], phi[k], -211, file_acceptance_2_261_pim);
    hctheta_pim -> Fill((acos(cthf[k]) * 180/M_PI));
	}

  if (pdgf[k] == 11){
    weights[k] = acceptance_c(momenta[k], cthf[k], phi[k], 11, file_acceptance_2_261);
  }
    } //end loop over final state particles


 //cout<<endl;
 cout << "----------- Entry #" << i << "----------" << endl;

 //cout <<"Number of true final state particles: "<< nf << endl;


    double smear_lP = gRandom->Gaus(leptonP,reso_e*leptonP);
    double smear_El = TMath::Sqrt(m_l*m_l+smear_lP*smear_lP);
    
    double e_acc_ratio = acceptance_c(leptonP, leptonCos, leptonPhi, 11, file_acceptance_2_261);


//COUNTING 1P1PI EVENTS 
  int protSeen = 0;
  int pipSeen = 0;
  int pimSeen = 0;

  int protAvail = 0;
  int piAvail = 0;

  //Possibly could be included in the previous loop, but don't think its a problem that its not
  for (int q = 0; q < nf ; q++){
  //cout << " - Particle " << q << " has pdg code " << pdgf[q] << ", and " << "weight: " << weights[q] << endl;
    if (pdgf[q] == 2212){
    
      protAvail+=1;
      if (MCDecide(weights[q])==1){
      protSeen += 1;
      hEprot_pass -> Fill(kinetic[q]);
      hMomprot_pass -> Fill(momenta[q]);
      
      hAng_prot_pass->Fill(((phi[q])*180/M_PI), (acos(cthf[q]) * 180/M_PI));
   
      hctheta_prot_pass ->Fill((acos(cthf[q]) * 180/M_PI));
      if (res==1){
         hEprot_pass_res -> Fill(kinetic[q]);
         hMomprot_pass_res -> Fill(momenta[q]);
      }
      if (dis==1){
        hEprot_pass_dis -> Fill(kinetic[q]);
        hMomprot_pass_dis -> Fill(momenta[q]);
      }
      }

      else{
        hEprot_fail -> Fill(kinetic[q]);
        hMomprot_fail -> Fill(momenta[q]);
        hAng_prot_fail->Fill(((phi[q])*180/M_PI), (acos(cthf[q]) * 180/M_PI));

         if (res==1){
         hEprot_fail_res -> Fill(kinetic[q]);
         hMomprot_fail_res -> Fill(momenta[q]);
         }
         if (dis==1){
        hEprot_fail_dis -> Fill(kinetic[q]);
        hMomprot_fail_dis -> Fill(momenta[q]);
         } 
      }
    }

     if (pdgf[q] == 211){
       piAvail+=1;
      if (MCDecide(weights[q])==1){
      pipSeen += 1;
      //cout << "pi +" << endl;
      hEpip_pass -> Fill(kinetic[q]);
      hMompip_pass -> Fill(momenta[q]);
     
      hAng_pip_pass->Fill(((phi[q])*180/M_PI), (acos(cthf[q]) * 180/M_PI));
      hctheta_pip_pass -> Fill((acos(cthf[q]) * 180/M_PI));
 
      if (res==1){
         hEpip_pass_res -> Fill(kinetic[q]);
         hMompip_pass_res -> Fill(momenta[q]);
      }
      if (dis==1){
        hEpip_pass_dis -> Fill(kinetic[q]);
        hMompip_pass_dis -> Fill(momenta[q]);
      }
      }

      else{
        cout << "Failed pion reason :" << (acos(cthf[q]) * 180/M_PI) << " ----" << (phi[q]) << "----" << momenta[q] << endl;
        hEpip_fail -> Fill(kinetic[q]);
        hMompip_fail -> Fill(momenta[q]);
        hAng_pip_fail->Fill((phi[q]), (acos(cthf[q]) * 180/M_PI));

         if (res==1){
         hEpip_fail_res -> Fill(kinetic[q]);
         hMompip_fail_res -> Fill(momenta[q]);
         }
         if (dis==1){
        hEpip_fail_dis -> Fill(kinetic[q]);
        hMompip_fail_dis -> Fill(momenta[q]);
         } 
      }
    }

   if (pdgf[q] == -211){
     piAvail += 1;
      if (MCDecide(weights[q])==1){
      pimSeen += 1;
      //cout << "pi -" << endl;
      hEpim_pass -> Fill(kinetic[q]);
      hMompim_pass -> Fill(momenta[q]);
      hAng_pim_pass->Fill(((phi[q])*180/M_PI), (acos(cthf[q]) * 180/M_PI));
      hctheta_pim_pass -> Fill((acos(cthf[q]) * 180/M_PI));

    

       if (res==1){
         hEpim_pass_res -> Fill(kinetic[q]);
         hMompim_pass_res -> Fill(momenta[q]);
      }
      if (dis==1){
        hEpim_pass_dis -> Fill(kinetic[q]);
        hMompim_pass_dis -> Fill(momenta[q]);
      }

      }
      else{
        hEpim_fail -> Fill(kinetic[q]);
        hMompim_fail -> Fill(momenta[q]);
        hAng_pim_fail->Fill((phi[q]), (acos(cthf[q]) * 180/M_PI));
         if (res==1){
         hEpim_fail_res -> Fill(kinetic[q]);
         hMompim_fail_res -> Fill(momenta[q]);
         }
         if (dis==1){
        hEpim_fail_dis -> Fill(kinetic[q]);
        hMompim_fail_dis -> Fill(momenta[q]);
         } 

      }
    }
    //cout << "particle:" << pdgf[q] << "// weight" << weights[q] << endl;

  } //end loop over final state particles (second)

auto piSeen = pipSeen + pimSeen;

//cout << "protseen: " << protSeen << endl;

//cout << "piseen: " << piSeen << endl;

  if (protSeen==1 && piSeen==1){
    signal+=1;
    hLepP -> Fill(smear_lP);
    //cout << i << "IS SIGNAL" << endl;
    if (res==1){
      rescount += 1;
      hWres_pass -> Fill(W);
    }

    if (dis==1){
      hWdis_pass -> Fill(W);
    }

      if (protAvail == 1 && piAvail==1){
    true1p1pi += 1;
    if (qel==1){
      qetrue1p1pi += 1;

    }
    if (res==1){
      restrue1p1pi += 1;
    }
    if (mec==1){
      mectrue1p1pi += 1;
    }
    if (dis==1){
      distrue1p1pi += 1;
    }
  }
  else if (protAvail > 1 || piAvail > 1){
    false1p1pi += 1;
     if (qel==1){
      qefalse1p1pi += 1;
    }
    if (res==1){
      resfalse1p1pi += 1;
    }
    if (mec==1){
      mecfalse1p1pi += 1;
    }
    if (dis==1){
      disfalse1p1pi += 1;
    }
  }
  }



  } // end big loop over entries
  cout << "Signal count = " << signal << endl;
  	std::cout << std::endl << "Signal proportion = " << int(double(signal) / double(nentries)*100.) << " \%" << std::endl;
	std::cout << std::endl << "Resonant proportion = " << int(double(rescount) / double(signal)*100.) << " \%" << std::endl;

cout << "1p1pi avail = " << true1p1pi <<  endl;
cout << "of which qe = " << qetrue1p1pi << endl;
cout << "of which res = " << restrue1p1pi << endl;
cout << "of which mec = " << mectrue1p1pi << endl;
cout << "of which dis = " << distrue1p1pi << endl;

cout << endl;

cout << "Accepted by losing some" << false1p1pi << endl;
cout << "of which qe = " << qefalse1p1pi << endl;
cout << "of which res = " << resfalse1p1pi << endl;
cout << "of which mec = " << mecfalse1p1pi << endl;
cout << "of which dis = " << disfalse1p1pi << endl;

    delete gRandom;

  //string targetString = GetTargetString(tree);
  string  PATH = "nanPlots/" + outdir;
  string filename_nopath = inFileName.substr(inFileName.find_last_of("/")+1); //gets filename from filepath
  string outFileName = PATH+string(filename_nopath);
  //string outFileName = PATH + "energies";
  //  TFile *output = new TFile(inFileName.c_str(),"RECREATE"); //open a file

  
  TFile *output = new TFile(outFileName.c_str(),"RECREATE"); //makes the file writeable, only 1 file can be open and writeable at a time in ROOT
  //write histograms here

 hLepP->Write();
 hEprot_pass -> Write();
 hEpip_pass -> Write();
 hEpim_pass -> Write();
 hEprot_fail -> Write();
 hEpip_fail -> Write();
 hEpim_fail -> Write();

  hEprot_pass_res -> Write();
 hEpip_pass_res -> Write();
 hEpim_pass_res -> Write();
 hEprot_fail_res -> Write();
 hEpip_fail_res -> Write();
 hEpim_fail_res -> Write();

  hEprot_pass_dis -> Write();
 hEpip_pass_dis -> Write();
 hEpim_pass_dis -> Write();
 hEprot_fail_dis -> Write();
 hEpip_fail_dis -> Write();
 hEpim_fail_dis -> Write();

 hMomprot_pass -> Write();
 hMompip_pass -> Write();
 hMompim_pass -> Write();
 hMomprot_fail -> Write();
 hMompip_fail -> Write();
 hMompim_fail -> Write();

  hMomprot_pass_res -> Write();
 hMompip_pass_res -> Write();
 hMompim_pass_res -> Write();
 hMomprot_fail_res -> Write();
 hMompip_fail_res -> Write();
 hMompim_fail_res -> Write();

  hMomprot_pass_dis -> Write();
 hMompip_pass_dis -> Write();
 hMompim_pass_dis -> Write();
 hMomprot_fail_dis -> Write();
 hMompip_fail_dis -> Write();
 hMompim_fail_dis -> Write();

 hAng_prot_pass -> Write();
 hAng_pip_pass -> Write();
 hAng_pim_pass -> Write();
 hAng_prot_fail -> Write();
 hAng_pip_fail -> Write();
 hAng_pim_fail -> Write();

hctheta_prot -> Write();
hctheta_prot_pass -> Write();
hctheta_pip -> Write();
hctheta_pip_pass -> Write();
hctheta_pim -> Write();
hctheta_pim_pass -> Write();

hWdis_pass -> Write();
hWres_pass -> Write();
  
    output -> Close();

 
 std::cout << "Histograms filled. Output closed!" << std::endl;

 
 return 0;
}
