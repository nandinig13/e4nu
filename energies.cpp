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
#include "acceptance_c.cpp"

#include "OscillationHelper.hxx"
R__LOAD_LIBRARY(WrappedProb3++.20121225.so)

using namespace std;

Double_t squared(Double_t num) {
	return TMath::Power(num, 2);
}

// ---------- Function to find charge of a particle (Ewart) ----------
int findCharge(int pdg) {
	int apdg = abs(pdg);
	if(apdg<100) return 0; //skips bosons and other weird stuff, which are handled elsewhere
	else if(apdg>1000000000) return 0; //skip nuclei
	int q[3];
	q[0] = apdg%10000/1000;
    q[1] = apdg%1000/100;
    q[2] = apdg%100/10;
    int charge = 0;
    for(int i = 0;i<3;i++){
        if(q[i]%2==0 && q[i]!=0){
            charge+=2;
        	}
        else if(q[i]%2==1){
            charge+=1;
        	}
    	}
    switch(charge){
        case 2: return 0; //dd' quark antiquark or similar
        case 3: return pdg/apdg; //ud or ddd, charge determined by sign of pdg
        case 4: return 0; //uu' or ud'd'
        //Baryons
        case 5: return pdg/apdg; //uud'
        case 6: return 2*pdg/apdg; //uuu
        default: throw "Error! Invalid charge!";
    }
}

void textOut(string PATH, string name, TH1D* hist) {
	ofstream out;
	string savename = PATH + name;
	out.open(savename.c_str());
	for (int j=0; j < hist->GetNbinsX(); j++){
		out << hist->GetBinCenter(j+1) << "  ";
		out << hist->GetBinContent(j+1) << "\n";
	}
}

double useMigrationMatrices(TH2D* hist, double Ev) {
	int binx;
        if (Ev < 5) {		// upper limit of migration matrix is 5 GeV
		binx = hist->GetXaxis()->FindBin(Ev);
	} else {
		binx = 50;	// max bin of migration matrices
	}
	TH1D *py = new TH1D();
	py = hist->ProjectionY("py",binx,binx);
	py->ComputeIntegral();
	double recoval = py->GetRandom();
	if (Ev < 5) {
		return recoval/Ev;
	} else {
		return recoval/5;
	}	
}


int energies(string filepath, string outdir, bool clas_maps, bool do_osc, bool clas_feeddown, char* clasfile, char* hist, bool xbcut, bool qecut) { //command line inputs for the run_energies.exe 
	
	//Select file options:

        string FILEPATH = filepath; //sets a string to thefile path to the root file specified when you input filepath in the command line
	//cout << FILEPATH << endl;
	//char* NAME = filename; //source file name and save name basis
	//string SRCPATH = "/macros/energy_rec/"; //path to source file
	string save = "Fe56energies_no_smearing";
	string PATH = "energy_rec/" + save; //path to save directory, change after the + to what you want it saves as
	//const string READNAME = string(NAME)+".root/gst"; //makes filename to access GST
	const string READPATH = string(FILEPATH)+"/gst";

	Double_t m_l = 0.000510999; //0.000510999 for electron, 0.105658 for muon
	Double_t beam_energy = 2.261; //GeV

	//cout << filename << endl;
	//cout << filepath << endl; //print filepath
	//cout << (string(SRCPATH)+READNAME).c_str() << endl;
	cout << PATH << endl; //print path to save directory
	if (do_osc) { cout << "oscillating" << endl; } else { cout << "not oscillating" << endl; }
	if (qecut) { cout << "using QE cuts" << endl; } else { cout << "no QE cuts" << endl; }
	if (xbcut) { cout << "using xb cut" << endl; } else { cout << "no xb cut" << endl; }
	if (clas_maps) { cout << "using acceptance maps" << endl; } else { cout << "no maps" << endl; }
	if (clas_feeddown) {cout << "feeddown using " << clasfile << endl; } else { cout << "no feeddown" << endl; }
	
	bool clas_acceptance = clas_maps;	 // use CLAS acceptance maps
	bool osc_weights = do_osc;		 	 // use Prob3 oscillation weights
	bool clasdata = clas_feeddown;		 // use CLAS feeddown things (either Mariana's or Afro's)
	
	enum Cuts// : int
	{
		QE_fp_cut, //select for QE by requiring a single final-product particle as a proton (GENIE)
		Q2_cut,		//Q2 > Q2_cutoff
		p_cut, 		//events that have exactly one proton with momentum > 0.3 GeV/c and no other charged particles with momentum > p_cutoff
		QEL_cut, 	//events that are quasi-elastic (GENIE)
		CC_cut, 	//events that are charged current (if applicable)
		pion_cut, 	//no produced pions with momentum > pion_cutoff (GeV)
		W_cut, 		//invariant mass W < W_cutoff (GeV)
		El_cut,		// cuts on final state lepton depending on incoming energy (1.161: El > 0.4 , 2.261: El > 0.55, 4.461:  El > 1.1)
		bjx_cut, 	//Bjorken x factor abs(x-1) < bjx_cutoff
		RES_cut,	//events that are resonant (GENIE)
		Ev_cut,		//cut on a range of Ev
		CUTS_NUM
	};

	//char const* cutNames[] = {"QE_fp_cut", "Q2_cut", "p_cut", "QEL_cut", "CC_cut", "pion_cut", "W_cut", "bjx_cut"};
	bool cut_switches[CUTS_NUM] = {false};  //These flags will switch on which cuts to make on events
	bool cut_flags[CUTS_NUM] = {false}; 	//These flags will be turned true if an event passes the cuts
	char const* cut_names[CUTS_NUM] = {""};
	//Cut options: true == cut is taken
	cut_switches[CC_cut] = false;	 //charged-current events only: turn on for nu, turn off for e- events (because nc==0 && cc==0 for them)! 
	cut_switches[Q2_cut] = qecut; 	 //Q2 > Q2_cutoff (1.161: Q2 > 0.1 , 2.261: Q2 > 0.4 , 4.461: Q2 > 0.8)
	cut_switches[p_cut] = qecut; 	 //events that have exactly one proton with momentum > 0.3 GeV/c and no other charged particles with momentum > p_cutoff
	cut_switches[W_cut] = qecut;      //invariant mass W < W_cutoff (GeV)
	cut_switches[El_cut] = qecut;      // cuts on final state lepton depending on incoming energy (1.161: El > 0.4 , 2.261: El > 0.55, 4.461:  El > 1.1)
	cut_switches[pion_cut] = qecut;	 //no produced pions with momentum > pion_cutoff (GeV)
	cut_switches[bjx_cut] = xbcut;	 //Bjorken x factor abs(x-1) < bjx_cutoff
	cut_switches[QE_fp_cut] = false; //GENIE-given QEL channel by taking only single proton final particles---keep off
	cut_switches[QEL_cut] = false;   //GENIE-given QEL channel---keep off. Can turn on for CCQE run comparison (if so manually change save name to include _CCQE)
	cut_switches[RES_cut] = false;	 //only includes resonant channel
	cut_switches[Ev_cut] = true;	 //do Ev cut!
	
	cut_names[CC_cut] = "genieCCevents";	 //charged-current events only: turn on for nu, turn off for e- events (because nc==0 && cc==0 for them)! 
	cut_names[Q2_cut] = "Q2cut"; 	 //Q2 > Q2_cutoff
	cut_names[p_cut] = "1fspAboveThreshold"; 	 //events that have exactly one proton with momentum > 0.3 GeV/c and no other charged particles with momentum > p_cutoff
	cut_names[W_cut] = "Wcut";      //invariant mass W < W_cutoff (GeV)
	cut_names[pion_cut] = "noPionsAvobeThreshold";	 //no produced pions with momentum > pion_cutoff (GeV)
	cut_names[bjx_cut] = "bjxcut";	 //Bjorken x factor abs(x-1) < bjx_cutoff
	cut_names[El_cut] = "Elcut";	 //Bjorken x factor abs(x-1) < bjx_cutoff
	cut_names[QE_fp_cut] = "genie1p"; //GENIE-given QEL channel by taking only single proton final particles---keep off
	cut_names[QEL_cut] = "genieQEL";   //GENIE-given QEL channel---keep off. Can turn on for CCQE run comparison (if so manually change save name to include _CCQE)
	cut_names[RES_cut] = "genieRES";	 //only includes resonant channel
	cut_names[Ev_cut] = "EvCut";	 //do Ev cut!
	
	const double Q2_cutoff = 0.5; //GeV. 0.5 GeV for 2.2 GeV analysis, 1.0 GeV for 4.4 GeV analysis
	const double Q2_cutoff_1_161 = 0.1; //GeV. 
	const double Q2_cutoff_2_261 = 0.4; //GeV. 
	const double Q2_cutoff_4_461 = 0.8; //GeV. 
	const double El_cutoff_1_161 = 0.4; //GeV. 
	const double El_cutoff_2_261 = 0.55; //GeV. 
	const double El_cutoff_4_461 = 1.1; //GeV. 
	const double p_cutoff = 0.3;  //GeV. final-state proton momentum cutoff
	const double pion_cutoff = 0.15; //GeV. 0 to exclude pions, 70 MeV for MicroBooNE threshold, 150 MeV for JLab threshold
	const double photon_cutoff = 0.3; //GeV. 0 to exclude pions, 70 MeV for MicroBooNE threshold, 150 MeV for JLab threshold
	const double W_cutoff = 2.0;  //GeV
	const double bjx_cutoff = 0.2;
	const double Ev_lower = 0;
       	const double Ev_upper = 10;
	
	// --------- PHYSICAL CONSTANTS (from NIST unless otherwise specified) --------- 
	const Double_t Eb = 0.036;	//binding energies from Afro: 0.025 for 12C, 0.04 for 40Ar, 0.036 for 56Fe
	const Double_t m_p =  0.938272; //mass of proton
	const Double_t m_n =  0.939565; //mass of neutron
	const Double_t m_pic = 0.13958; //mass of charged pion https://journals.aps.org/pr/abstract/10.1103/PhysRev.163.1451 
	const Double_t m_pi0 = 0.13498; //(GeV/c^2) mass of neutral pion (wikipedia, couldn't find better)
	const Double_t fine_struc_const = 1./137.035999139;
	//binding energies: 34 for , 30 for 


	// --------- OSCILLATION SETUP & PARAMETERS -----------------
  	// DUNE beam dip
	double beam_dip_deg = 5.8;
	// Sin^2(Theta_12)
	double s2th12 = 0.310;
	// Sin^2(Theta_13)
	double s2th13 = 0.02241;
	// Sin^2(Theta_23)
	double s2th23 = 0.580;
	// Dm^2_21
	double dm2_21 = 7.39E-5;
	//|Dm^2_Atm|
	double dm2_atm = 2.525E-3;
	// dcp
	double dcp = 0; //TMath::Pi();
	double osc_params[] = {s2th12, s2th13, s2th23, dm2_21, dm2_atm, dcp};
	TGraph POsc; 
	OscillationHelper oh_nu;
	oh_nu.Setup_dipangle(osc_params, beam_dip_deg);
	int pdg_nu = 14;		// numu = 14, nue = 12
	oh_nu.SetOscillationChannel(pdg_nu, 12);
	
	TList *l = new TList();

	// --------- GET FEED-DOWN DATA -----------------
	TFile *fracfeed = 0;
	TH1D *clashist = 0;
	if (clasdata) {
		fracfeed = new TFile(clasfile);				//get file with electron data
		clashist = (TH1D*)fracfeed->Get(hist);		//get relevant histogram
		clashist->ComputeIntegral();				//needed for GetRandom() to work properly later
	}

	// --------- CLAS ACCEPTANCE MAPS -----------------
	TFile* file_acceptance = TFile::Open("maps/e2a_maps_12C_E_2_261.root");
	TFile* file_acceptance_p = TFile::Open("maps/e2a_maps_12C_E_2_261_p.root");
	TFile* file_acceptance_1_161 = TFile::Open("maps/e2a_maps_12C_E_1_161.root");
	TFile* file_acceptance_1_161_p = TFile::Open("maps/e2a_maps_12C_E_1_161_p.root");
	TFile* file_acceptance_2_261 = TFile::Open("maps/e2a_maps_12C_E_2_261.root");
	TFile* file_acceptance_2_261_p = TFile::Open("maps/e2a_maps_12C_E_2_261_p.root");

	TFile *matfile = new TFile("MigrationMatrices.root");
	TH2D *matrix_genie = (TH2D*)matfile->Get("MC12CMigrationMatrix");
	TH2D *matrix_data = (TH2D*)matfile->Get("DT12CMigrationMatrix");

	TChain *tree = new TChain(); //for whole-directory processing
	//tree -> Add((string(SRCPATH)+READNAME).c_str(),0); //get tree
	tree -> Add((string(READPATH)).c_str(),0); //get tree from source file
	cout << (string(READPATH)).c_str() << endl;
	Long64_t nentries = tree->GetEntries(); //For troubleshooting, can set to 10 eventsi
        //cout<< "Number of entries is" <<nentries<<endl;	
	//nentries = 10; //This is TEMPORARY for debugging!
	Double_t Ev;
	Double_t Ev_fd;
	Double_t El, pxl, pyl, pzl;		// for quasielastic reconstruction
	
	//calorimetric reconstruction variables
	tree ->	SetBranchStatus("*", 0);
	tree -> SetBranchStatus("nf", 1);
	Int_t nf;	
	int maxNF = 0;
        tree->SetBranchAddress("nf",&nf);
    	// Loop through entries to find max
	for(Long64_t j=0;j<nentries;j++){
		tree->GetEntry(j);   
	        if(nf>maxNF){
	            maxNF = nf;
	        }
      	}
	Double_t Ef[maxNF];
	Int_t pdgf[maxNF];
	Double_t pxf[maxNF], pyf[maxNF], pzf[maxNF];

	//variables for cuts and filters
	Double_t Q2, W, x;
	Bool_t res, qel, cc;
	Double_t wght;
	Bool_t dis, mec;

	tree -> SetBranchStatus("Ev", 1);
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
	tree -> SetBranchStatus("res", 1);
	tree -> SetBranchStatus("qel", 1);
	tree -> SetBranchStatus("dis", 1);
	tree -> SetBranchStatus("mec", 1);
	tree -> SetBranchStatus("cc", 1);
	tree -> SetBranchStatus("wght", 1);

	tree -> SetBranchAddress("Ev", &Ev);
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
	tree -> SetBranchAddress("res", &res);
	tree -> SetBranchAddress("qel", &qel);
	tree -> SetBranchAddress("dis", &dis);
	tree -> SetBranchAddress("mec", &mec);
	tree -> SetBranchAddress("cc", &cc);
	tree -> SetBranchAddress("wght", &wght);
	
	double HMAX = 12.0;
	double BINS = 280;
	// variable binning for 	
	int n_bins;
	double *x_values;
	double *x_qe;
	if(beam_energy>1. && beam_energy<2.){
	  //cout<< "beam energy thing" <<beam_energy <<endl;
	        n_bins=38;
		x_values = new double[n_bins+1]; 
		x_qe = new double[n_bins+1];
		for (int i=0;i<=17;i++) { 
			x_values[i]=0.4+i*0.04; 
			x_qe[i] = (x_values[i] - beam_energy) / beam_energy;
		} 
		for (int i=0;i<=20;i++) { 
			x_values[i+18]=1.08+(i+1)*0.02; 
			x_qe[i+18] = (x_values[i+18] - beam_energy) / beam_energy; 
		}
	}

	
	TH1D *hEv = new TH1D("hEv", "neutrino energy (GENIE)", BINS, 0.0, HMAX); 
	TH1D *hEQE = new TH1D("hEQE", "quasielastic energy reconstruction", BINS, 0.0, HMAX); 
	TH1D *hcal = new TH1D("hcal", "calorimetric energy reconstruction", BINS, 0.0, HMAX); 
	
	TH1D *hEv_qe = new TH1D("hEv_qe", "neutrino energy in QE events (GENIE)", BINS, 0.0, HMAX);
	TH1D *hEv_res = new TH1D("hEv_res", "neutrino energy in resonant events (GENIE)", BINS, 0.0, HMAX);
	TH1D *hEv_dis = new TH1D("hEv_dis", "neutrino energy in DIS events (GENIE)", BINS, 0.0, HMAX);
	TH1D *hEv_mec = new TH1D("hEv_mec", "neutrino energy in MEC events (GENIE)", BINS, 0.0, HMAX);
	TH1D *hEQE_qe = new TH1D("hEQE_qe", "quasielastic energy reconstruction in QE events", BINS, 0.0, HMAX);
	TH1D *hEQE_res = new TH1D("hEQE_res", "quasielastic energy reconstruction in resonant events", BINS, 0.0, HMAX);
	TH1D *hEQE_dis = new TH1D("hEQE_dis", "quasielastic energy reconstruction in DIS events", BINS, 0.0, HMAX);
	TH1D *hEQE_mec = new TH1D("hEQE_mec", "quasielastic energy reconstruction in MEC events", BINS, 0.0, HMAX);
	TH1D *hcal_qe = new TH1D("hcal_qe", "calorimetric energy reconstruction in QE events", BINS, 0.0, HMAX);
	TH1D *hcal_res = new TH1D("hcal_res", "calorimetric energy reconstruction in resonant events", BINS, 0.0, HMAX);
	TH1D *hcal_dis = new TH1D("hcal_dis", "calorimetric energy reconstruction in DIS events", BINS, 0.0, HMAX);
	TH1D *hcal_mec = new TH1D("hcal_mec", "calorimetric energy reconstruction in MEC events", BINS, 0.0, HMAX);

	TH1D *hEv_fd = new TH1D("hEv_fd", "neutrino energy with feed-down (GENIE)", BINS, 0.0, HMAX); 
	

	TH1D *hres_cal = new TH1D("hres_cal", "calorimetric reconstruction energy resolution", n_bins, x_qe); 
	TH1D *hres_fd = new TH1D("hres_fd", "feed-down reconstruction energy resolution", n_bins, x_qe); 
	TH2D *hcomp = new TH2D("hcomp", "true neutrino energy (GENIE) vs feed-down energy", BINS, 0.0, HMAX, BINS, 0.0, HMAX);
	cout << "histograms set up" << endl;
	
	gRandom = new TRandom3();
	gRandom->SetSeed(10);

	for( Long64_t i = 0; i < nentries ; i++) {
		
		tree -> GetEntry(i);

		//MAKE CUTS: determine if the ith event passes available cuts
		cut_flags[QE_fp_cut] = (nf==1 && pdgf[0]==2212);
		// cut_flags[Q2_cut] = (Q2 > Q2_cutoff);
		cut_flags[Q2_cut] = ((Ev < 1.5 && Q2 > Q2_cutoff_1_161) || (Ev >= 1.5 && Q2 > Q2_cutoff_2_261));
		cut_flags[El_cut] = ((Ev < 1.5 && El > El_cutoff_1_161) || (Ev >= 1.5 && El > El_cutoff_2_261));
		cut_flags[QEL_cut] = (qel);
		cut_flags[RES_cut] = (res);
		cut_flags[CC_cut] = (cc);
		cut_flags[W_cut] = (W < W_cutoff);
		cut_flags[bjx_cut] = (TMath::Abs(x-1) < bjx_cutoff);

		cut_flags[p_cut] = true;    //passes unless proved otherwise
		cut_flags[pion_cut] = true; //passes unless proved otherwise

		//	cout << "passes QE_fp" <<cut_flags[QE_fp_cut]<<endl;



		TVector3 lP(pxl, pyl, pzl);		//get lepton information
		double leptonP = lP.Mag();
		double leptonCos = lP.CosTheta();
		double leptonPhi = lP.Phi() + TMath::Pi();
		double E_QE = ( squared(m_n) - squared(m_p - Eb) - squared(m_l) + 2*(m_p - Eb)*El ) / ( 2*( m_p - Eb - El + leptonP * leptonCos ) );
		
		double pP = 0;			//proton momentum
		double protonK = 0;		//proton kinetic energy
		double pCos = 0;		//cosine of scattering angle 
		double pPhi = 0;		//phi
		int p_count_cut = 0;	//for single-proton cut

		for(Int_t k=0; k<nf; k++) { //loop through k final particles
			int pdg = pdgf[k];
			TVector3 kVec(pxf[k], pyf[k], pzf[k]);	//momentum vector for kth particle
			double kP = kVec.Mag();

			//charged, non-proton particle with momentum exceeding p_cutoff
			//if((abs(pdg)==0 || abs(pdg)==11 || abs(pdg)==13 || abs(pdg)==15 || abs(pdg)==24 || abs(pdg)==34 || abs(pdg)==37 || (findCharge(abs(pdg))!=0 && pdg!=2212)) && kP>p_cutoff) {
			//	cut_flags[p_cut] = false;
			//}
			//pions
			if(abs(pdg)==211) {
				if (kP > pion_cutoff) {
					cut_flags[pion_cut] = false;
				}
			}
			if(abs(pdg)==22) {
				if ( kP > photon_cutoff) {
						cut_flags[pion_cut] = false;
				}
			}

			//if proton
			if (pdg == 2212) {
				if (kP > pP) {			//get primary proton
					pP = kP;
					pCos = kVec.CosTheta(); 
					pPhi = kVec.Phi() + TMath::Pi(); 
					protonK = Ef[k] - m_p;
				}
				if (kP>p_cutoff) p_count_cut++;
			}
		}

		if (p_count_cut != 1) cut_flags[p_cut] = false;
		cut_flags[Ev_cut] = Ev_lower <= Ev && Ev <= Ev_upper;

		//Resolutions for Smearing for GENIE simulation data
		double reso_p = 0.01; // smearing for the proton
		double reso_e = 0.005; // smearing for the electrons
		//double reso_pipl = 0.007; //smearing for pions, executive decision by Larry (28.08.19)
		//double reso_pimi = 0.007; //smearing for pions, executive decision by Larry (28.08.19)
		double reso_pi = 0.007; //smearing for pions, executive decision by Larry (28.08.19)
		// Resolution defined above seems to be insufficient at 1.1 GeV -> tripled it for all particles
		if(beam_energy == 1.161) { 
			reso_p = 3*reso_p; reso_e = 3*reso_e; reso_pi = 3*reso_pi;
		}
		double smear_pP = gRandom->Gaus(pP,reso_p*pP);
		double smear_lP = gRandom->Gaus(leptonP,reso_e*leptonP);
		double smear_El = TMath::Sqrt(m_l*m_l+smear_lP*smear_lP);
		double smear_pK = TMath::Sqrt(m_p*m_p+smear_pP*smear_pP) - m_p;
		//double calE = smear_El + smear_pK + Eb; //smear calorimetric one
		double calE = El + protonK + Eb; //don't smear
		
		// ------ CALCULATE WEIGHTS ---------
		// using CLAS acceptance data
		double weight = 1;
		if (clas_acceptance == true) {
			// double e_acc_ratio = acceptance_c(leptonP, leptonCos, leptonPhi, 11, file_acceptance);
			// double p_acc_ratio = acceptance_c(pP, pCos, pPhi, 2212, file_acceptance_p);
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
		// add oscillation weights
		if (osc_weights == true) {
    			double pnu = oh_nu.GetWeight(Ev);			
				POsc.SetPoint(i, Ev, pnu);
				weight *= pnu;
		}
		if (fabs(weight) != weight) continue;
		
		// using fractional feed-down
		double scale = 1;
		if (clasdata == true) {
			//cout<<"Using clasdata data smearing "<<endl;
			scale = clashist->GetRandom() / beam_energy;  // <--- THIS IS TO USE THE FEED-DOWN STUFF
			//scale = useMigrationMatrices(matrix_genie, Ev) / Ev;
			// cout << scale << endl;
			Ev_fd = Ev * scale;
			/*double tot = 0;
			double Ev_rec;
			double wt;
			for (int b = 1; b <= clashist->GetNbinsX(); b++) {
			       Ev_rec = Ev * clashist->GetBinCenter(b) / 2.2;
			       wt = clashist->GetBinContent(b) / clashist->Integral();
			       hEv -> Fill(Ev_rec, wt);
			}*/	       
		}

		double cal_res = (calE - Ev) / Ev;
		double fd_res = (Ev_fd - Ev) / Ev;
                

		// --------- FILL HISTOGRAMS --------- 
		//Determine if event passed all implemented cuts
		bool pass = true;
		//cout <<"event number =" << i << " ";
		for(int k = 0; k < CUTS_NUM; k++) 
		  {//if(cut_switches[k]) cout << k << ":" << cut_flags[k]<<" ";  
			if(cut_switches[k]==true && cut_flags[k]==false) {
				pass = false;
			}
		}
                //cout << endl;
		//cout <<"event number =" << i <<endl;
		if (pass == true) {
                        
			//if (clasdata == false) {
				hEv -> Fill(Ev, weight);		
			//}
			if (clasdata) {
				hEv_fd -> Fill(Ev_fd, weight);
				hcomp -> Fill(Ev, Ev_fd, weight);
			}	
			hEQE -> Fill(E_QE, weight);		
			hcal -> Fill(calE, weight);
			hres_cal -> Fill(cal_res, weight);
			hres_fd -> Fill(fd_res, weight);	
			// hcomp -> Fill(Ev, calE, weight);
			
			if (qel == true) {		
				hEv_qe -> Fill(Ev, weight);
				hEQE_qe -> Fill(E_QE, weight);
				hcal_qe -> Fill(calE, weight);
			} if (res == true) {
				hEv_res-> Fill(Ev, weight);
				hEQE_res-> Fill(E_QE, weight);
				hcal_res-> Fill(calE, weight);
			} if (dis == true) {
				hEv_dis-> Fill(Ev, weight);
				hEQE_dis-> Fill(E_QE, weight);
				hcal_dis-> Fill(calE, weight);
			} if (mec == true) {
				hEv_mec-> Fill(Ev, weight);
				hEQE_mec-> Fill(E_QE, weight);
				hcal_mec-> Fill(calE, weight);
			}
		}
	}
	delete gRandom;
	//cout<<"The neutrino energy is "<<Ev<<endl;
        //cout<<x<<endl;
	// I have commented out the part where the code outputs .txt files as I will not be needing them.
      	//.txt output
	//	textOut(PATH, "EQE.txt", hEQE);
	//	textOut(PATH, "CAL.txt", hcal);
	//textOut(PATH, "NU.txt", hEv);
	//textOut(PATH, "NUFD.txt", hEv_fd);
	//textOut(PATH, "NUQE.txt", hEv_qe);
	//textOut(PATH, "NURES.txt", hEv_res);
	//textOut(PATH, "NUDIS.txt", hEv_dis);
	//textOut(PATH, "NUMEC.txt", hEv_mec);
	//textOut(PATH, "EQEQE.txt", hEQE_qe);
	//textOut(PATH, "EQERES.txt", hEQE_res);
	//textOut(PATH, "EQEDIS.txt", hEQE_dis);
	//textOut(PATH, "EQEMEC.txt", hEQE_mec);
	//textOut(PATH, "CALQE.txt", hcal_qe);
	//textOut(PATH, "CALRES.txt", hcal_res);
	//textOut(PATH, "CALDIS.txt", hcal_dis);
	//textOut(PATH, "CALMEC.txt", hcal_mec);	
	//textOut(PATH, "RESFD.txt", hres_fd);	
	//textOut(PATH, "RESCAL.txt", hres_cal);	


	//string fileName = PATH+"energy_reconstruction_"+NAME+".root";
	//rootFILEPATHNoPath=rootFILEPATH.substr(rootFILEPATH.find_last_of("/")+1);
	string FILEPATHNoPath = FILEPATH.substr(FILEPATH.find_last_of("/")+1); //gets the filename from filepath
	std::cout << "file no path: " << FILEPATHNoPath << std::endl;
        string fileName = PATH+string(FILEPATHNoPath);

	TFile *output = new TFile(fileName.c_str(),"RECREATE");
	
	hEv -> Write();
	hEQE -> Write();
	hcal -> Write();
	hcomp -> Write();
	
	if (clasdata == true) { clashist->Write(); hEv_fd->Write();} 
	
	hEv_qe -> Write();
	hEv_res -> Write();
	hEv_dis -> Write();
	hEv_mec -> Write();
	hEQE_qe -> Write();
	hEQE_res -> Write();
	hEQE_dis -> Write();
	hEQE_mec -> Write();
	hcal_qe -> Write();
	hcal_res -> Write();
	hcal_dis -> Write();
	hcal_mec -> Write();
	hres_cal -> Write();
	hres_fd -> Write();

	if (osc_weights == true) { POsc.Write("POsc");}

	output -> Close();
	cout << "histograms filled, output closed!" << endl;
	return 0;	

}
