#include "OscillationHelper.hxx"

//#undef R__LOAD_LIBRARY
//#define R__LOAD_LIBRARY(a)

#include "energies.cpp"
#include <string>
using namespace std;

int main(int argc, char *argv[]) {
	char filename[30];
	switch(argv[1][0]) {
		case('d'): strcpy(filename, "gst-nu_dune");  break;
		case('c'): strcpy(filename, "gst-nu_carbon"); break;  
		case('i'): strcpy(filename, "gst-nu_iron"); break;
		case('m'): strcpy(filename, "p14_c12_22GeV"); break;
		case('1'): strcpy(filename, "C12_nue_100k"); break;
		case('2'): strcpy(filename, "C12_nue_500k"); break;
		case('3'): strcpy(filename, "C12_nue_600k"); break;
		case('4'): strcpy(filename, "C12_nue_1M"); break;
		case('5'): strcpy(filename, "C12_nuebar_100k"); break;
		case('6'): strcpy(filename, "C12_numu_100k"); break;
		case('7'): strcpy(filename, "C12_numubar_100k"); break;
 		case('e'): strcpy(filename, "eresmaid_C12_1161_RadCorr"); break;
	        case('a'): strcpy(filename, "apapadop_SuSav2_C12_1161GeV"); break;  
	        
	}
	char clasfile[300];
	char hist[300];
	switch(argv[6][0]) {
		case('m'): 
			if (strcmp(argv[1],"carbon")) {
				strcpy(clasfile, "mariana/e2a_ep_C12_2261_neutrino6.root");
			} else {
				strcpy(clasfile, "mariana/e2a_ep_56Fe_2261_neutrino6.root");
			}
			strcpy(hist, "h_Etot_subtruct_piplpimi_prot");
			break;
		case('d'):
			if (!strcmp(argv[7],"1")){	// if bjorken x cut is happening
				strcpy(clasfile, "afroclas/xBCut/");
			} else {
				strcpy(clasfile, "afroclas/NoxBCut/");
			}				
			//if (strcmp(argv[1],"carbon") || strcmp(argv[1],"mono")) {
				strcat(clasfile, "12C_1_161_Data_Final_Plots_FSI_em.root");
			//} else {
			//	strcat(clasfile, "56Fe_4_461_Data_Final_Plots_FSI_em.root");
			//}
			strcpy(hist, "epRecoEnergy_slice_0");
			break;
		case('n'):
			strcpy(clasfile, "afrosim/old/12C_1_161_hA2018_NoRadCorr_Plots_FSI_em.root");
			strcpy(hist, "epRecoEnergy_slice_0");
			break;
		case('r'):
			//strcpy(clasfile, "afrosim/old/12C_2_261_hA2018_RadCorr_Plots_FSI_em.root");
			//if (strcmp(argv[1],"carbon")) {
				strcpy(clasfile, "afrosim/12C_1_161_hA2018_Final_RadCorr_LFGM_Plots_FSI_em.root");
			//} else {
			//	strcpy(clasfile, "afrosim/56Fe_4_461_hA2018_Final_RadCorr_LFGM_Plots_FSI_em.root");
			//}
			strcpy(hist, "epRecoEnergy_slice_0");
			break;
		case('x'):
			strcpy(clasfile, "");
			strcpy(hist, "");
			break;	
	}
	energies(filename, argv[2], !strcmp(argv[3],"1"), !strcmp(argv[4],"1"), !strcmp(argv[5],"1"), clasfile, hist, !strcmp(argv[7],"1"), !strcmp(argv[8],"1"));
}
