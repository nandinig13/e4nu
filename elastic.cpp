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

//PHYSICAL CONSTANTS      //can maybe one day put these in a header                                                        
const Double_t ELECTRON_MASS= 0.00051099; //Gev 
const Double_t MUON_MASS= 0.105658; //Gev     
const Double_t BEAM_ENERGY = 2.261; //GeV  //Change for different beam energies
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
  FINAL_STATE_NO_OTHER, // no other mesons in final state
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
    case FINAL_STATE_NO_OTHER: return  (nfpim + nfpip + nfpi0 + nfk0 + nfkm +nfkp == 0); 
         break; // 1 final state pion
    case FINAL_STATE_PROTON_CUT: return (nfp == 1); 
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

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

//MAIN FUNCTION
int func (string inFileName, string outdir, string tester){ //not taking in any complicated parameters yet

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

//CUTS
  vector<Cut*> cuts;
    //cuts.push_back(new Cut("Bjorken x cut","TMath::Abs(x-1) < 0.2")); 
    cuts.push_back(new Cut("no_mes" , "no_mes", FINAL_STATE_NO_OTHER));
    cuts.push_back(new Cut("1p" , "1p", FINAL_STATE_PROTON_CUT));
    cuts.push_back(new Cut("0n" , "0n", FINAL_STATE_NEUTRON_CUT));

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

  int Bins = 100;
  double Histo_xmax = 7;
  
  TH1D *hEv = new TH1D("Ev", cutText.c_str(), Bins, 0.0, Histo_xmax); //create a pointer to a histogram
  TH1D *hcal = new TH1D("Cal", "Calorimetric energy reconstruction", Bins, 0.0, Histo_xmax); 
  TH1D *hkin = new TH1D("Kin", "Kinematic energy reconstruction", Bins, 0.0, Histo_xmax);

  TH1D *hcal_qe = new TH1D("cal_qe", "Calorimetric energy reconstruction in QE events", Bins, 0.0, Histo_xmax);
  TH1D *hcal_res = new TH1D("cal_res", "Calorimetric energy reconstruction in RES events", Bins, 0.0, Histo_xmax);
  TH1D *hcal_dis = new TH1D("cal_dis", "Calorimetric energy reconstruction in DIS events", Bins, 0.0, Histo_xmax);
  TH1D *hcal_mec = new TH1D("cal_mec", "Calorimetric energy reconstruction in MEC events", Bins, 0.0, Histo_xmax);


  TH1D *hkin_qe = new TH1D("kin_qe", "Kinematic energy reconstruction in QE events", Bins, 0.0, Histo_xmax);
  TH1D *hkin_res = new TH1D("kin_res", "Kinematic energy reconstruction in RES events", Bins, 0.0, Histo_xmax);
  TH1D *hkin_dis = new TH1D("kin_dis", "Kinematic energy reconstruction in DIS events", Bins, 0.0, Histo_xmax);
  TH1D *hkin_mec = new TH1D("kin_mec", "Kinematic energy reconstruction in MEC events", Bins, 0.0, Histo_xmax);
  
/*
  // --------- CLAS ACCEPTANCE MAPS -----------------
//TFile* file_acceptance_1_161 = TFile::Open("maps/e2a_maps_12C_E_1_161.root");
//TFile* file_acceptance_1_161_p = TFile::Open("maps/e2a_maps_12C_E_1_161_p.root");
//TFile* file_acceptance_1_161_pip = TFile::Open("maps/e2a_maps_12C_E_1_161_pip.root");
//TFile* file_acceptance_1_161_pim = TFile::Open("maps/e2a_maps_12C_E_1_161_pim.root");

TFile* file_acceptance_2_261 = TFile::Open("maps/e2a_maps_12C_E_2_261.root");
TFile* file_acceptance_2_261_p = TFile::Open("maps/e2a_maps_12C_E_2_261_p.root");
TFile* file_acceptance_2_261_pip = TFile::Open("maps/e2a_maps_12C_E_2_261_pip.root");
TFile* file_acceptance_2_261_pim = TFile::Open("maps/e2a_maps_12C_E_2_261_pim.root"); 
 */
  //GAUSSIAN SMEARING

  gRandom = new TRandom3();
  gRandom->SetSeed(10);

  //SIGNAL EVENT COUNTER
  int QESignalEvents = 0;
  int MECSignalEvents = 0;
  int RESSignalEvents = 0;
  int DISSignalEvents = 0;  


  for( Long64_t i = 0; i < nentries ; i++){
  if (i%500==0){
      cout << "running" << endl;
  }
  
      tree -> GetEntry(i);
      TVector3 lP(pxl, pyl, pzl);   //get lepton information
      double leptonP = lP.Mag();
      double leptonCos = lP.CosTheta();
      double leptonPhi = lP.Phi() + TMath::Pi();
      

      double pP = 0;      //proton momentum
      double protonK = 0; //proton kinetic energy
      double pCos = 0;    //cosine of scattering angle 
      double pPhi = 0;    //phi
      int pCount = 0;


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
	}
      }
    }

//RESOLUTIONS FOR SMEARING
    double reso_p = 0.01; // smearing for the proton
    double reso_e = 0.005; // smearing for the electron
    double reso_pi = 0.007; //smearing for pions
    if(BEAM_ENERGY == 1.161) {  //justification for this in Anjali's paper - lower energy has more error
			reso_p = 3*reso_p; reso_e = 3*reso_e; reso_pi = 3*reso_pi;
		}

//STILL SETTING VALUES
    double m_l; //lepton mass = muon mass for neutrino scat, lepton mass = electron mass for electron scat 
    if (isElectronMode == true){
      m_l = ELECTRON_MASS;
    }
    else {
      m_l = MUON_MASS;
      cout<<"Not in electron mode!"<<endl;
    }

//SMEARING


    
    double smear_pP = gRandom->Gaus(pP,reso_p*pP);
    double smear_lP = gRandom->Gaus(leptonP,reso_e*leptonP);
    double smear_El = TMath::Sqrt(m_l*m_l+smear_lP*smear_lP);
    double smear_pK = TMath::Sqrt(PROTON_MASS*PROTON_MASS+smear_pP*smear_pP) - PROTON_MASS;
    
    //CALORIMETRIC ENERGY RECONSTRUCTION
    double calE = smear_El + smear_pK + Eb;

    //KINEMATIC ENERGY RECONSTRUCTION
    double kinE = ( squared(PROTON_MASS) - squared(PROTON_MASS-Eb)-squared(m_l) + 2*(PROTON_MASS - Eb)*smear_El ) / ( 2*(PROTON_MASS - Eb - smear_El + smear_lP*leptonCos ) );


/*
bool clas_acceptance = true;  //might want to pass this as an argument one day
  if (clas_acceptance == true) {
			double e_acc_ratio = 1.;
			double p_acc_ratio = 1.; 
			if (Ev <= 1.5) {
			   e_acc_ratio = acceptance_c(leptonP, leptonCos, leptonPhi, 11, file_acceptance_1_161);
			   p_acc_ratio = acceptance_c(pP, pCos, pPhi, 2212, file_acceptance_1_161_p);
			}
			else {
			   e_acc_ratio = acceptance_c(leptonP, leptonCos, leptonPhi, 11, file_acceptance_2_261);
			   p_acc_ratio = acceptance_c(pP, pCos, pPhi, 2212, file_acceptance_2_261_p);
			}
			weight *= e_acc_ratio * p_acc_ratio;
		}	
*/


// MAKING THE CUTS
    string cutString="";
   for (int i=0; i<cuts.size();i++) // It loops all your cuts
   {
    if (i>0)cutString+=" && ";{
        cutString+=cuts.at(i)->GetCutString();
    }
   }
   
    bool passescuts = true;
    for (int k = 0; k < cuts.size(); k++)
    {
	    if (cuts.at(k)->PassesCut() == false) {
	        passescuts = false;
	    } 
    }
//cout << passescuts << endl;
    if (passescuts == true){

      hcal->Fill(calE);
      hkin->Fill(kinE);
      hEv->Fill(Ev);

      if (qel == true){
    QESignalEvents += 1;
	hcal_qe->Fill(calE);
	hkin_qe->Fill(kinE);

      }
      if (res == true){
  RESSignalEvents += 1;
	hcal_res->Fill(calE);
	hkin_res->Fill(kinE);
      }

      if (dis == true){
  DISSignalEvents += 1;
	hcal_dis->Fill(calE);
	hkin_dis->Fill(kinE);

      }
      if (mec == true){
  MECSignalEvents += 1;
	hcal_mec->Fill(calE);
	hkin_mec->Fill(kinE);
      }





    } // end if passes cuts
  } // end big loop over entries

    delete gRandom;

  //PRINTING FRACTIONAL CONTRIBUTIONS
  std::cout << std::endl << "QE Fractional Contribution = " << int(double(QESignalEvents) / double(nentries) *100.) << " \%" << std::endl;
	std::cout << std::endl << "MEC Fractional Contribution = " << int(double(MECSignalEvents) / double(nentries)*100.) << " \%" << std::endl;
	std::cout << std::endl << "RES Fractional Contribution = " << int(double(RESSignalEvents) / double(nentries)*100.) << " \%" << std::endl;
	std::cout << std::endl << "DIS Fractional Contribution = " << int(double(DISSignalEvents) / double(nentries)*100.) << " \%" << std::endl;

  //string targetString = GetTargetString(tree);
  string  PATH = "nanPlots/" + outdir;
  string filename_nopath = inFileName.substr(inFileName.find_last_of("/")+1); //gets filename from filepath
  string outFileName = PATH+string(filename_nopath);
  //string outFileName = PATH + "energies";
  //  TFile *output = new TFile(inFileName.c_str(),"RECREATE"); //open a file

  
  TFile *output = new TFile(outFileName.c_str(),"RECREATE"); //makes the file writeable, only 1 file can be open and writeable at a time in ROOT
    
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

    
    output -> Close();

 
 std::cout << "Histograms filled. Output closed!" << std::endl;

 
 return 0;
}
