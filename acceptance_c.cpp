#include <TH3D.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>

using namespace std;

double acceptance_c(double p, double cost, double phi, int id,TFile* file_acceptance) {

	//Redefinition of the phi angle
	
	int redef = -30;
	//int redef = 0;

	TH3D * e_acc;
	TH3D * e_gen;

	TH3D * p_acc;
	TH3D * p_gen;

	TH3D * pip_acc;
	TH3D * pip_gen;

	TH3D * pim_acc;
	TH3D * pim_gen;

        // Electron

	if (id == 11) {

		if (redef == -30) {

			e_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
			e_gen = (TH3D*)file_acceptance->Get("Generated Particles");
		}

		if (redef == 0) {

			e_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
			e_gen = (TH3D*)file_acceptance->Get("Generated Particles");
		}		
	
		//Find number of generated events

		double e_pbin_gen = e_gen->GetXaxis()->FindBin(p);
		double e_tbin_gen = e_gen->GetYaxis()->FindBin(cost);
		double e_phibin_gen = e_gen->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
		double e_num_gen = e_gen->GetBinContent(e_pbin_gen, e_tbin_gen, e_phibin_gen);
		//Find number of accepted events

		double e_pbin_acc = e_acc->GetXaxis()->FindBin(p);
		double e_tbin_acc = e_acc->GetYaxis()->FindBin(cost);
		double e_phibin_acc = e_acc->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
		double e_num_acc = e_acc->GetBinContent(e_pbin_acc, e_tbin_acc, e_phibin_acc);

		double e_acc_ratio = (double)e_num_acc / (double)e_num_gen;
    		double e_acc_err = (double)sqrt(e_acc_ratio*(1-e_acc_ratio)) / (double)e_num_gen;
		return e_acc_ratio;
	}

        // Proton

	if (id == 2212) {

		if (redef == -30){

			p_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
			p_gen = (TH3D*)file_acceptance->Get("Generated Particles");
		}

		if (redef == 0){

			p_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
			p_gen = (TH3D*)file_acceptance->Get("Generated Particles");
		}

		//Find number of generated events

		double p_pbin_gen = p_gen->GetXaxis()->FindBin(p);
		double p_tbin_gen = p_gen->GetYaxis()->FindBin(cost);
		double p_phibin_gen = p_gen->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
		double p_num_gen = p_gen->GetBinContent(p_pbin_gen, p_tbin_gen, p_phibin_gen);
		//Find number of accepted events

		double p_pbin_acc = p_acc->GetXaxis()->FindBin(p);
		double p_tbin_acc = p_acc->GetYaxis()->FindBin(cost);
		double p_phibin_acc = p_acc->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
		double p_num_acc = p_acc->GetBinContent(p_pbin_acc, p_tbin_acc, p_phibin_acc);
		double p_acc_ratio = (double)p_num_acc / (double)p_num_gen;
		double p_acc_prr = (double)sqrt(p_acc_ratio*(1-p_acc_ratio)) / (double)p_num_gen;

		return p_acc_ratio;
	}


        // Pi+                                                                                                           

        if (id == 211) {

	  if (redef == -30){

	    pip_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
	    pip_gen = (TH3D*)file_acceptance->Get("Generated Particles");
	  }

	  if (redef == 0){

	    pip_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
	    pip_gen = (TH3D*)file_acceptance->Get("Generated Particles");
	  }

	  //Find number of generated events                                                                           

	  double pip_pbin_gen = pip_gen->GetXaxis()->FindBin(p);
	  double pip_tbin_gen = pip_gen->GetYaxis()->FindBin(cost);
	  double pip_phibin_gen = pip_gen->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
	  double pip_num_gen = pip_gen->GetBinContent(pip_pbin_gen, pip_tbin_gen, pip_phibin_gen);
	  //Find number of accepted events                                                                            

	  double pip_pbin_acc = pip_acc->GetXaxis()->FindBin(p);
	  //cout << "pip_pbin_acc" << pip_pbin_acc << endl;
	  double pip_tbin_acc = pip_acc->GetYaxis()->FindBin(cost);
	  //cout << "pip_tbin_acc" << pip_tbin_acc << endl;
	  double pip_phibin_acc = pip_acc->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
	 // cout << "pip_phibin_acc" << pip_phibin_acc << endl;
	  double pip_num_acc = pip_acc->GetBinContent(pip_pbin_acc, pip_tbin_acc, pip_phibin_acc);
	  double pip_acc_ratio = (double)pip_num_acc / (double)pip_num_gen;
	  double pip_acc_prr = (double)sqrt(pip_acc_ratio*(1-pip_acc_ratio)) / (double)pip_num_gen;
	  	//cout<<"acc"<<pip_num_acc<<endl;                                                                       
	  	//cout<<"gen"<<pip_num_gen<<endl;
	  return pip_acc_ratio;
        }
        // Pi-                                                                                                                                                                       
	if (id == -211) {

          if (redef == -30){

            pim_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
            pim_gen = (TH3D*)file_acceptance->Get("Generated Particles");
          }

          if (redef == 0){

            pim_acc = (TH3D*)file_acceptance->Get("Accepted Particles");
            pim_gen = (TH3D*)file_acceptance->Get("Generated Particles");
          }

          //Find number of generated events                                                                             

          double pim_pbin_gen = pim_gen->GetXaxis()->FindBin(p);
          double pim_tbin_gen = pim_gen->GetYaxis()->FindBin(cost);
          double pim_phibin_gen = pim_gen->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
          double pim_num_gen = pim_gen->GetBinContent(pim_pbin_gen, pim_tbin_gen, pim_phibin_gen);
          //Find number of accepted events                                                                              

          double pim_pbin_acc = pim_acc->GetXaxis()->FindBin(p);
          double pim_tbin_acc = pim_acc->GetYaxis()->FindBin(cost);
          double pim_phibin_acc = pim_acc->GetZaxis()->FindBin(phi*180/TMath::Pi()+redef);
          double pim_num_acc = pim_acc->GetBinContent(pim_pbin_acc, pim_tbin_acc, pim_phibin_acc);
          double pim_acc_ratio = (double)pim_num_acc / (double)pim_num_gen;
          double pim_acc_prr = (double)sqrt(pim_acc_ratio*(1-pim_acc_ratio)) / (double)pim_num_gen;
	  return pim_acc_ratio;
	}

return 0.;

}