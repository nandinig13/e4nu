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

#include "OscillationHelper.hxx"
R__LOAD_LIBRARY(WrappedProb3++.20121225.so)

using namespace std;

//Physical constants                                                              
const Double_t ELECTRON_MASS= 0.00051099; //Gev 
const Double_t MUON_MASS= 0.105658; //Gev     
const Double_t beam_energy = 2.261; //GeV  //DO I NEED TO CHANGE THIS FOR DIFFERENT BEAM ENERGY  
const Double_t Eb = 0.025;//binding energies from Afro: 0.025 for 12C, 0.04 for 40Ar, 0.03 for 56Fe NEED TO CHECK VALIDITY                                                          
const Double_t PROTON_MASS =  0.938272; //mass of proton                                           
const Double_t NEUTRON_MASS =  0.939565; //mass of neutron                                        
const Double_t PION_MASS = 0.13958; //mass of charged pion https://journals.aps.org/pr/abstract/10.1103/PhysRev.163.1451
const Double_t DELTA_MASS = 1.232; //mass of delta state                                                                 
const Double_t NEUTRAL_PION_MASS = 0.13498; //(GeV/c^2) mass of neutral pion (wikipedia, couldn't find better)                                                                                   
const Double_t FINE_STRUCTURE_CONSTANT = 1./137.035999139; //FIND SOURCE

Double_t squared(Double_t num) {
  return TMath::Power(num, 2);
}

bool CheckIfElectrons(TTree *tree);

enum CutID //enumerated type
{
  FINAL_STATE_CHARGED_PION, // final state pion cut
  FINAL_STATE_PROTON, //  final state proton cut
  FINAL_STATE_NEUTRON,
  FINAL_STATE_NEUTRAL_PION,
};

//Global variables from tree:
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
bool nc;//is it a neutral current event?    
int resid; // the ID of the baryonic resonance state if res = true
int em; // is it an electromagnetic interaction?

class Cut {
  // All this stuff is private, you can't see it from outside the class
  string cutString_; // You can do something here like with the text one
  string titleText_;                                                       //Objects
  CutID cutID_; //make enum for CutID where we add cuts
public:
  // This is a constructor - it runs when you make a new Cut object.
  // We can use it to initialize the cut
  Cut(string _titleString, string _cutString, CutID _cutID);
  // We can add other functions here, maybe you want the efficiency of the cut for example...

  string GetTitleText(); //Methods
  string GetCutString();
  bool PassesCut();//returns true or false depending on if cut id passes e.g. TRUE_RES_CUT=true then give res...
};

Cut::Cut(string _titleString, string _cutString, CutID _cutID ) //etc
{
  this->titleText_=_titleString;
  this->cutString_=_cutString;
  this->cutID_=_cutID;//
}

// Get the text and string
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
  case FINAL_STATE_CHARGED_PION: return  (nfpim + nfpip == 1); 
    break; // 1 final state pion
  case FINAL_STATE_PROTON: return (nfp >= 1); 
    break;
  case FINAL_STATE_NEUTRON: return (nfn == 0);
    break;
  case FINAL_STATE_NEUTRAL_PION: return (nfpi0 == 0);
    break;
     default: return true;
  }
}
// Returns true if it is electromagentic electron scattering, and false if not
bool CheckIfElectrons(TTree *tree)//
{
  bool isElectron=false;
  tree->SetBranchStatus("em", 1);
  tree->SetBranchAddress("em", &isElectron); // If this is true (electromagnetic scattering) then it isn't a neutrino
  tree->GetEntry(); // The first entry is fine, this should tell us all we need to know. Note that it would be misleading if we had a mixed neutrino and electron beam... but we don't.
  return isElectron;
}

// Return the text corresponding to the target material
string GetTargetString(TTree *tree)
{
  int nProtons;
  string targetString="";
  tree->SetBranchStatus("Z", 1);
  tree->SetBranchAddress("Z", &nProtons); // The number of protons in the target nucleus will tell us what element it is. Possible pitfalls - there COULD be scattering from atomic electrons but that is quite rare, and it's possible that instead of simulating from a single target material we could scatter from a complex detector made of many materials - but I don't think your simulation has that. You can check by histogramming Z - if it's the same for every entry then you are good.
  tree->GetEntry();
  switch (nProtons)
  {
  case 18: targetString="argon"; break;
  case 6: targetString="carbon"; break;
  case 26: targetString="iron"; break;
  case 1: targetString="hydrogen"; break;
  case 2: targetString="helium"; break;  
  default: "unknown target"; // You can add iron to the list...
  }
  return targetString;
}

int pions(string inFileName, string outdir) {
  //  cout<<"Cuts "<<CutID<<endl;
  cout<<"Analysing "<<inFileName<<endl;
  // Load the input root file

  TFile *f = new TFile(inFileName.c_str());
  f->ls();
  cout << (string(inFileName)).c_str() << endl;
  
  // Load the gst tree from the input file
  TTree *tree = (TTree*) f->Get("gst");
  
  Long64_t nentries = tree->GetEntries(); //eset nentries to a small number for troubleshooting
        nentries = 1000000; //for debugging only
  
    // Is it electron or neutrino mode?
  bool isElectronMode = CheckIfElectrons(tree);
    
  // Target material
  string targetString = GetTargetString(tree);
  
  // Print out what sort of scattering it is
  cout<<(isElectronMode?"Electron":"Neutrino")<<" scattering on "<<targetString<<endl;

  // CUTS
  // Make a vector of all the cuts you want to use
  // You could make groups that add (push_back) various combos of these cuts.
  // A cut group can be an object, just like a cut can be an object
  // You could then have text or other info associated with the groups
  // as well/instead of with the individual cuts.
  
  vector<Cut*> cuts;
    //cuts.push_back(new Cut("Bjorken x cut","TMath::Abs(x-1) < 0.2")); //and so on
     //cuts.push_back(new Cut("1pi" , "1pi", FINAL_STATE_CHARGED_PION));
     //cuts.push_back(new Cut("1p" , "1p", FINAL_STATE_PROTON));
     //cuts.push_back(new Cut("0pi0" , "0pi0", FINAL_STATE_NEUTRAL_PION));
     //cuts.push_back(new Cut("0n", "0n", FINAL_STATE_NEUTRON));
     //cuts.push_back(new Cut("deltaRES", "deltaRES", RESONANCE_STATE));
     cuts.push_back(new Cut(">= 1p" , ">= 1p", FINAL_STATE_PROTON));
     cuts.push_back(new Cut("1pi" , "1pi", FINAL_STATE_CHARGED_PION));
//     cout<<cuts<<endl;

  // I've put these just as text but you can use a function to define more complicated cuts if you want to - just create a Cut object and add it to your vector
  
  // Make a text string that combines all your cut texts
  
  string cutText="Cuts: ";
  for (int i=0; i<cuts.size();i++) // It loops all your cuts
  {
    if (i>0)cutText+=" and ";
    cutText+=cuts.at(i)->GetTitleText();
  }
    
  cout<<cutText<<endl;

  //loop through entries to find maximum
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

  tree -> SetBranchStatus("*", 0);
  tree -> SetBranchStatus("Ev", 1);
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

//

  tree -> SetBranchAddress("Ev", &Ev);
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
  int Bins = 250;
  double Histo_xmax = 10;
  
  TH1D *hEv = new TH1D("Ev", cutText.c_str(), Bins, 0.0, Histo_xmax); //create a pointer to a histogram
  TH1D *hcal = new TH1D("Cal", "Calorimetric energy reconstruction", Bins, 0.0, Histo_xmax); 
  TH1D *hkin = new TH1D("Kin", "Kinematic energy reconstruction", Bins, 0.0, Histo_xmax);

  TH1D *hEv_qe = new TH1D("Ev_qe", "Lepton energy in QE events (GENIE)", Bins, 0.0, Histo_xmax);
  TH1D *hEv_res = new TH1D("Ev_res", "Lepton energy in RES events (GENIE)", Bins, 0.0, Histo_xmax);
  TH1D *hEv_dis = new TH1D("Ev_dis", "Lepton energy in DIS events (GENIE)", Bins, 0.0, Histo_xmax);
  TH1D *hEv_mec = new TH1D("Ev_mec", "Lepton energy in MEC events (GENIE)", Bins, 0.0, Histo_xmax);

  TH1D *hcal_qe = new TH1D("cal_qe", "Calorimetric energy reconstruction in QE events", Bins, 0.0, Histo_xmax);
  TH1D *hcal_res = new TH1D("cal_res", "Calorimetric energy reconstruction in RES events", Bins, 0.0, Histo_xmax);
  TH1D *hcal_dis = new TH1D("cal_dis", "Calorimetric energy reconstruction in DIS events", Bins, 0.0, Histo_xmax);
  TH1D *hcal_mec = new TH1D("cal_mec", "Calorimetric energy reconstruction in MEC events", Bins, 0.0, Histo_xmax);

  TH1D *hcal_pip = new TH1D("cal_pip", "Calorimetric energy reconstruction events with 1 pi+", Bins, 0.0, Histo_xmax);
  TH1D *hcal_qe_pip = new TH1D("cal_qe_pip", "Calorimetric energy reconstruction in QE events with 1 pi+", Bins, 0.0, Histo_xmax);
  TH1D *hcal_res_pip = new TH1D("cal_res_pip", "Calorimetric energy reconstruction in RES events with 1 pi+", Bins, 0.0, Histo_xmax);
  TH1D *hcal_dis_pip = new TH1D("cal_dis_pip", "Calorimetric energy reconstruction in DIS events with 1 pi+", Bins, 0.0, Histo_xmax);
  TH1D *hcal_mec_pip = new TH1D("cal_mec_pip", "Calorimetric energy reconstruction in MEC events with 1 pi+", Bins, 0.0, Histo_xmax);
  TH1D *hcal_pim = new TH1D("cal_pim", "Calorimetric energy reconstruction events with 1 pi-", Bins, 0.0, Histo_xmax);
  TH1D *hcal_qe_pim = new TH1D("cal_qe_pim", "Calorimetric energy reconstruction in QE events with 1 pi-", Bins, 0.0, Histo_xmax);
  TH1D *hcal_res_pim = new TH1D("cal_res_pim", "Calorimetric energy reconstruction in RES events with 1 pi-", Bins, 0.0, Histo_xmax);
  TH1D *hcal_dis_pim = new TH1D("cal_dis_pim", "Calorimetric energy reconstruction in DIS events with 1 pi-", Bins, 0.0, Histo_xmax);
  TH1D *hcal_mec_pim = new TH1D("cal_mec_pim", "Calorimetric energy reconstruction in MEC events with 1 pi-", Bins, 0.0, Histo_xmax);


  TH1D *hkin_qe = new TH1D("kin_qe", "Kinematic energy reconstruction in QE events", Bins, 0.0, Histo_xmax);
  TH1D *hkin_res = new TH1D("kin_res", "Kinematic energy reconstruction in RES events", Bins, 0.0, Histo_xmax);
  TH1D *hkin_dis = new TH1D("kin_dis", "Kinematic energy reconstruction in DIS events", Bins, 0.0, Histo_xmax);
  TH1D *hkin_mec = new TH1D("kin_mec", "Kinematic energy reconstruction in MEC events", Bins, 0.0, Histo_xmax);
  TH1D *hkin_pip = new TH1D("kin_pip", "Kinematic energy reconstruction in events with 1 pi+", Bins, 0.0, Histo_xmax);  
  TH1D *hkin_qe_pip = new TH1D("kin_qe_pip", "Kinematic energy reconstruction in QE events with 1 pi+", Bins, 0.0, Histo_xmax);
  TH1D *hkin_res_pip = new TH1D("kin_res_pip", "Kinematic energy reconstruction in RES events with 1 pi+", Bins, 0.0, Histo_xmax);
  TH1D *hkin_dis_pip = new TH1D("kin_dis_pip", "Kinematic energy reconstruction in DIS events with 1 pi+", Bins, 0.0, Histo_xmax);
  TH1D *hkin_mec_pip = new TH1D("kin_mec_pip", "Kinematic energy reconstruction in MEC events with 1 pi+", Bins, 0.0, Histo_xmax);
  TH1D *hkin_pim = new TH1D("kin_pim", "Kinematic energy reconstruction in events with 1 pi-", Bins, 0.0, Histo_xmax);
  TH1D *hkin_qe_pim = new TH1D("kin_qe_pim", "Kinematic energy reconstruction in QE events with 1 pi-", Bins, 0.0, Histo_xmax);
  TH1D *hkin_res_pim = new TH1D("kin_res_pim", "Kinematic energy reconstruction in RES events with 1 pi-", Bins, 0.0, Histo_xmax);
  TH1D *hkin_dis_pim = new TH1D("kin_dis_pim", "Kinematic energy reconstruction in DIS events with 1 pi-", Bins, 0.0, Histo_xmax);
  TH1D *hkin_mec_pim = new TH1D("kin_mec_pim", "Kinematic energy reconstruction in MEC events with 1 pi-", Bins, 0.0, Histo_xmax);

  TH1D *hkin_res_Delta_pip = new TH1D("kin_res_Delta_pip", "Kinematic energy reconstruction in RES events with 1 pi+ from a (1232) RES", Bins, 0.0, Histo_xmax);
  TH1D *hkin_res_NODelta_pip = new TH1D("kin_res_NODelta_pip", "Kinematic energy reconstruction in RES events with 1 pi+ NOT from a (1232) RES", Bins, 0.0, Histo_xmax);
  TH1D *hkin_res_Delta_pim = new TH1D("kin_res_Delta_pim", "Kinematic energy reconstruction in RES events with 1 pi- from a (1232) RES", Bins, 0.0, Histo_xmax);
  TH1D *hkin_res_NODelta_pim = new TH1D("kin_res_NODelta_pim", "Kinematic energy reconstruction in RES events with 1 pi- NOT from a (1232) RES", Bins, 0.0, Histo_xmax);
    // --------- CLAS ACCEPTANCE MAPS -----------------
  TFile* file_acceptance_1_161 = TFile::Open("maps/e2a_maps_12C_E_1_161.root");
  TFile* file_acceptance_1_161_p = TFile::Open("maps/e2a_maps_12C_E_1_161_p.root");
  TFile* file_acceptance_1_161_pip = TFile::Open("maps/e2a_maps_12C_E_1_161_pip.root");
  TFile* file_acceptance_1_161_pim = TFile::Open("maps/e2a_maps_12C_E_1_161_pim.root");

  TFile* file_acceptance_2_261 = TFile::Open("maps/e2a_maps_12C_E_2_261.root");
  TFile* file_acceptance_2_261_p = TFile::Open("maps/e2a_maps_12C_E_2_261_p.root");
  TFile* file_acceptance_2_261_pip = TFile::Open("maps/e2a_maps_12C_E_2_261_pip.root");
  TFile* file_acceptance_2_261_pim = TFile::Open("maps/e2a_maps_12C_E_2_261_pim.root"); 
  
    //define random gaussian for smearing
  gRandom = new TRandom3();
  gRandom->SetSeed(10);
  
  for ( Long64_t i = 0; i < nentries ; i++) {
    
    tree -> GetEntry(i);
    TVector3 lP(pxl, pyl, pzl);   //get lepton information
    double leptonP = lP.Mag();
    double leptonCos = lP.CosTheta();
    double leptonPhi = lP.Phi() + TMath::Pi();

    double pP = 0;      //proton momentum
    double protonK = 0;   //proton kinetic energy
    double pCos = 0;    //cosine of scattering angle 
    double pPhi = 0;    //phi

    double pipP = 0;      //pion+ momentum                                                       
    double pipK = 0;   //pion+ kinetic energy                                                  
    double pipCos = 0;    //cosine of scattering angle                                             
    double pipPhi = 0;    //phi  

    double pimP = 0;      //pion- momemtum                                                                                                     
    double pimK = 0;   //pion- kinetic energy        
    double pimCos = 0;    //cosine of scattering angle                                                                                                                              
    double pimPhi = 0;    //phi                                                                                                                                                      
    double Epi = 0;

    for(Int_t k=0; k<nf; k++) { //loop through k final particles
      int pdg = pdgf[k];
      TVector3 kVec(pxf[k], pyf[k], pzf[k]);//momentum vector for kth particle
      double kP = kVec.Mag();

      if (pdg == 2212) {
	if (kP > pP) {//get primary proton
	  pP = kP;
	  pCos = kVec.CosTheta(); 
	  pPhi = kVec.Phi() + TMath::Pi(); 
	  protonK = Ef[k] - PROTON_MASS;
	}
      }
      if (pdg == 211) {
	if (kP > pipP){//get primary pion+
	  pipP = kP;
	  pipCos = kVec.CosTheta();
          pipPhi = kVec.Phi() + TMath::Pi();
          pipK = Ef[k] - PION_MASS; 
	}
      }
      if (pdg == -211) {
	if (kP > pimP){//get primary pion+                                                                              
        pimP = kP;
        pimCos = kVec.CosTheta();
        pimPhi = kVec.Phi() + TMath::Pi();
        pimK = Ef[k] - PION_MASS;
        }
      }
    }
    //Resolutions for Smearing for GENIE simulation data
    double reso_p = 0.01; // smearing for the proton
    double reso_e = 0.005; // smearing for the electron
    double reso_pi = 0.007; //smearing for pions

       //Calorimetric energy reconstruction
    double m_l; //lepton mass = muon mass for neutrino scat, lepton mass = electron mass for electron scat 
    if (isElectronMode == true){
      m_l = ELECTRON_MASS;
    }
    else {
      m_l = MUON_MASS;
    }
    double smear_pipP = 0;
    double smear_Epi = 0 ;
    double smear_pimP = 0;
    if (nfpip ==1){
      smear_pipP = gRandom->Gaus(pipP,reso_pi*pipP);
      smear_Epi = TMath::Sqrt(PION_MASS*PION_MASS+smear_pipP*smear_pipP);
      Epi = pipK + PION_MASS;
    }     
    else if (nfpim ==1){
      smear_pimP = gRandom->Gaus(pimP,reso_pi*pimP);
      smear_Epi = TMath::Sqrt(PION_MASS*PION_MASS+smear_pimP*smear_pimP);
      Epi = pimK + PION_MASS;
    }
    else {
      Epi = 0;
      smear_Epi = 0;
    }
    

    double smear_pP = gRandom->Gaus(pP,reso_p*pP);
    double smear_lP = gRandom->Gaus(leptonP,reso_e*leptonP);
    double smear_El = TMath::Sqrt(m_l*m_l+smear_lP*smear_lP);
    double smear_pK = TMath::Sqrt(PROTON_MASS*PROTON_MASS+smear_pP*smear_pP) - PROTON_MASS;

    //    double calE = smear_El + smear_pK + Eb + smear_Epi;
    double calE = smear_El + smear_pK + Eb + smear_Epi;
    double kinE = ( squared(DELTA_MASS) - squared(PROTON_MASS-Eb)-squared(m_l) + 2*(PROTON_MASS - Eb)*smear_El ) / ( 2*(PROTON_MASS - Eb - smear_El + smear_lP*leptonCos ) );
    //double kinE = ( (NEUTRON_MASS)*(NEUTRON_MASS) - (PROTON_MASS - Eb)*(PROTON_MASS-Eb) - (m_l)*(m_l) + 2*(PROTON_MASS - Eb)*El )/(2*(PROTON_MASS-Eb-El+leptonP * leptonCos ) );

    double weight = 1.0;
    double weight_pim = 1.0;
    double weight_pip = 1.0;
    int nonzero = 0;

    bool clas_acceptance = true;   
    if (clas_acceptance == true){
      double e_acc_ratio = 1.0;
      double p_acc_ratio = 1.0;
      double pip_acc_ratio = 1.0;
      double pim_acc_ratio = 1.0;

	e_acc_ratio = acceptance_c(leptonP, leptonCos, leptonPhi, 11, file_acceptance_2_261);
	p_acc_ratio = acceptance_c(pP, pCos, pPhi, 2212, file_acceptance_2_261_p);
	pip_acc_ratio = acceptance_c(pipP, pipCos, pipPhi, 211, file_acceptance_2_261_pip);
	pim_acc_ratio = acceptance_c(pimP, pimCos, pimPhi, -211, file_acceptance_2_261_pim);

  cout << nfpip << endl;
  //cout << "initial weight" << weight_pip << endl;
  //cout << "e_acc" << e_acc_ratio << endl;
  //cout << "p_acc" << p_acc_ratio << endl;
  //cout << "pip_acc" << pip_acc_ratio << endl;
  if (pip_acc_ratio != 0){
    cout << ("non zero") << endl;
  }

	weight_pip *= e_acc_ratio * p_acc_ratio * pip_acc_ratio;

	weight_pim *= e_acc_ratio * p_acc_ratio * pim_acc_ratio;

  cout << "final weight" << weight_pip << endl;
  cout << endl;

    }

  // Now make the cuts
    string cutString="";
  for (int i=0; i<cuts.size();i++) // It loops all your cuts
  {
    if (i>0)cutString+=" && ";{
    cutString+=cuts.at(i)->GetCutString();
    }
  }
    bool passescut = true;
    for (int k = 0; k < cuts.size(); k++)
      {
	if (cuts.at(k)->PassesCut() == false) {
	   passescut = false;
	 } 
      }

    if (passescut == true){
      hcal->Fill(calE, weight);
      hkin->Fill(kinE, weight);

      if (qel == true){
	hcal_qe->Fill(calE, weight);
	hkin_qe->Fill(kinE, weight);
      }
      if (res == true){
	hcal_res->Fill(calE, weight);
	hkin_res->Fill(kinE, weight);
      }
      if (dis == true){
	hcal_dis->Fill(calE, weight);
	hkin_dis->Fill(kinE, weight);
      }
      if (mec == true){
	hcal_mec->Fill(calE, weight);
	hkin_mec->Fill(kinE, weight);
      }
      if (nfpip ==1){
	hcal_pip->Fill(calE, weight_pip);
	hkin_pip->Fill(kinE, weight_pip);
             if(qel == true){
      hEv_qe -> Fill(Ev, weight_pip);
      hcal_qe_pip -> Fill(calE, weight_pip);
      hkin_qe_pip -> Fill(kinE, weight_pip);
          }if(res == true){
         hEv_res -> Fill(Ev, weight_pip);
	 hcal_res_pip -> Fill(calE, weight_pip);
	 hkin_res_pip -> Fill(kinE, weight_pip);
	 if(resid ==0){
	   hkin_res_Delta_pip->Fill(kinE, weight);
	 }
	 if(resid !=0){
	   hkin_res_NODelta_pip->Fill(kinE, weight);
	 }
    }if(mec == true){
      hEv_mec -> Fill(Ev, weight_pip);
      hcal_mec_pip -> Fill(calE, weight_pip);
      hkin_mec_pip -> Fill(kinE, weight_pip);

    } if(dis == true){
	hEv_dis -> Fill(Ev, weight_pip);
	hcal_dis_pip -> Fill(calE, weight_pip);
	hkin_dis_pip -> Fill(kinE, weight_pip);

    }
      }
      if (nfpim ==1 ){
	hcal_pim->Fill(calE, weight_pim);
	hkin_pim->Fill(kinE, weight_pim);
	if(qel == true){
	  hEv_qe -> Fill(Ev, weight_pim);
	  hcal_qe_pim -> Fill(calE, weight_pim);
	  hkin_qe_pim -> Fill(kinE, weight_pim);
	}if(res == true){
	  hEv_res -> Fill(Ev, weight_pim);
	  hcal_res_pim -> Fill(calE, weight_pim);
	  hkin_res_pim -> Fill(kinE, weight_pim);
	  if(resid ==0){
	    hkin_res_Delta_pim->Fill(kinE, weight_pim);
	  }
	  if(resid !=0){
	    hkin_res_NODelta_pim->Fill(kinE, weight);
	  }
	}if(mec == true){
	  hEv_mec -> Fill(Ev, weight_pim);
	  hcal_mec_pim -> Fill(calE, weight_pim);
	  hkin_mec_pim -> Fill(kinE, weight_pim);
	} if(dis == true){
	  hEv_dis -> Fill(Ev, weight_pim);
	  hcal_dis_pim -> Fill(calE, weight_pim);
	  hkin_dis_pim -> Fill(kinE, weight_pim);
}
    }
    }
  }
     string  PATH = "nanPlots/" + outdir;
  string filename_nopath = inFileName.substr(inFileName.find_last_of("/")+1); //gets filename from filepath
  string outFileName = PATH+string(filename_nopath);
  //  TFile *output = new TFile(inFileName.c_str(),"RECREATE"); //open a file
  TFile *output = new TFile(outFileName.c_str(),"RECREATE");
      hEv -> Write();
    hcal -> Write();
    hkin -> Write();
    
    hcal_qe -> Write();
    hcal_res -> Write();
    hcal_dis -> Write();
    hcal_mec -> Write();

    hkin_qe -> Write();
    hkin_res -> Write();
    hkin_dis -> Write();
    hkin_mec -> Write();

    hEv_qe -> Write();
    hEv_res -> Write();
    hEv_mec -> Write();
    hEv_dis -> Write();
    hcal_pip -> Write();
    hkin_pip ->Write();
    hcal_qe_pip -> Write();
    hcal_res_pip -> Write();
    hcal_dis_pip -> Write();
    hcal_mec_pip -> Write();
    hcal_pim -> Write();
    hkin_pim ->Write();
    hcal_qe_pim -> Write();
    hcal_res_pim -> Write();
    hcal_dis_pim -> Write();
    hcal_mec_pim -> Write();
    hkin_qe_pip -> Write();
    hkin_res_pip -> Write();
    hkin_dis_pip -> Write();
    hkin_mec_pip -> Write();
    hkin_qe_pim -> Write();
    hkin_res_pim -> Write();
    hkin_dis_pim -> Write();
    hkin_mec_pim -> Write();

    hkin_res_Delta_pip -> Write();
    hkin_res_Delta_pim -> Write();
    hkin_res_NODelta_pip -> Write();
    hkin_res_NODelta_pim -> Write();
     output -> Close();
 cout << "histograms filled, output closed!" << endl;
 return 0;
}