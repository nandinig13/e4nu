//file to plot CLAS acceptance 
#include <string>
#include <TH3.h>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaletteAxis.h>
#include <TMath.h>
#include <TLine.h>
#include <TPad.h>
#include <TGaxis.h>
#include <iostream>
#include <vector>

using namespace std;
void acceptanceplot()
{
  gStyle->SetOptStat(false); // Turn off that ugly stats box                                                            

  // Create a TFile object and connect it to your file                                                                  
  //TFile *f = new TFile("energy_reconstruction_gst-nu_carbon.root");                                                   
  string acceptance_file  = "/genie/app/users/ngadhia/scripts/macros/maps/e2a_maps_12C_E_2_261_p.root ";
    // "f" is just a placeholder name at this point                                                                     

  TFile *f = new TFile(acceptance_file.c_str());


  f->ls(); // This will print out a list of histograms in your file     
  
   //TString Option = "xy";
    TString Option = "yz";
  //  TString Option = "xz";

   //x = momentum , y = cos theta,  z = phi
  
  TCanvas *c1 = new TCanvas("e4nu","e4nu",900,600); // We will overwrite the default titles. The numbers set the size
  
  TH3D *acc = (TH3D*) f->Get("Accepted Particles");
  TH3D *gen = (TH3D*) f->Get("Generated Particles");
  
  
  //  acc->GetYaxis()->SetRangeUser(0.6,1);
  //  gen->GetYaxis()->SetRangeUser(0.6,1);
  
  TH2D *acc2D = (TH2D*) acc->Project3D(Option);
  TH2D *gen2D = (TH2D*) gen->Project3D(Option);
    
  TH2D *accrat2D = (TH2D*) (acc2D->Clone());
  accrat2D->Divide(gen2D);
  accrat2D->SetTitle("CLAS acceptance for final state proton in electron scattering on C12 at 2.261GeV");
  //accrat2D->GetZaxis()->SetTitle("proton momentum[GeV]");
  //  accrat2D->GetYaxis()->SetTitle("proton momentum[GeV]");
  accrat2D->GetXaxis()->SetTitle("Phi(Degrees)");    
  accrat2D->GetYaxis()->SetTitle("Cos theta");
  
  gen2D->SetTitle("CLAS acceptance for electron scattering on C12 at 1.161GeV");
  gen2D->GetYaxis()->SetTitle("electron momentum[GeV]");
  gen2D->GetXaxis()->SetTitle("cos theta");
  
  gStyle->SetPalette(100); //100=kSolar colour palette
  accrat2D->Draw("colz");
  //gen2D->Draw("colz");
  //  acc2D->Draw("colz");
  c1->SaveAs("CLASacceptance.png");

}