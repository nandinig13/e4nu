#include "TColor.h"
#include <string>

using namespace std;

/* Run this from the root command line  with .x anjaliPlots.C */
void anjaliPlotsnew() // This name has to match the file name
{
  gStyle->SetOptStat(false); // Turn off that ugly stats box
    // Create a TFile object and connect it to your file
  string Filename  = "workingapapadop_UpdatedSchwingerRad_SuSav2_C12_2261GeV.root";
  TFile *f = new TFile(Filename.c_str());
  // "f" is just a placeholder name at this point
 
 
  f->ls(); // This will print out a list of histograms in your file
  
  // Get the hEv_qe histogram and store it in a TH1D object
  // (1-dimensional histogram of doubles i.e. high-precision decimals
  
   TH2D *h_EvQE = (TH2D*) f->Get("Ep_min_max_qe"); // Quasi-elastic
  // TH1D *h_EvRES = (TH1D*) f->Get("hEv_res"); // Resonant
  // TH1D *h_EvMEC = (TH1D*) f->Get("hEv_mec"); // Meson exchange currents
  // TH1D *h_EvDIS = (TH1D*) f->Get("hEv_dis"); // Deep inelastic scattering
  // Create a canvas to draw on
  TCanvas *canvas = new TCanvas("e4nu","e4nu",900,600); // We will overwrite the default titles. The numbers set the size
  THStack *hs = new THStack("hs", "Proton kinematic energies for events with 2 protons and 1 charged pion; Energy /GeV ; Event count");
  

  //Quasi-elastic
  // TH1D *h_EvQE = (TH1D*) f->Get("hEQE_qe"); // Quasi-elastic
  // h_EvQE->GetYaxis()->SetTitle("Event count");
    //h_EvQE->GetYaxis()->SetRangeUser(0,1700000);
  // h_EvQE->GetXaxis()->SetRangeUser(0,4); //Go up to 10Gev 
  //h_EvQE->GetXaxis()->SetTitle("E_{cal} (unsmeared) /GeV"); 
  //h_EvQE->SetTitle("Quasi-elastic reconstruction of electron energy for #e^- ^{12}C# at 2.261GeV (SuSav2).");
  //  h_EvQE->SetLineWidth(2);
  // h_EvQE->SetMarkerColor(kAzure-1);
  hs->Add(h_EvQE);
  //  hs->Draw("HIST");
  
  //Resonant
  TH2D *h_EvRES = (TH2D*) f->Get("Ep_min_max_res"); 
  //  h_EvRES->SetLineWidth(2);
  //   h_EvRES->SetLineColor(kPink);
   hs->Add(h_EvRES);
   // hs->Draw("HIST SAME C");
   
  //Meson-exchange
  TH2D *h_EvMEC = (TH2D*) f->Get("Ep_min_max_mec");
  // h_EvMEC->SetLineWidth(2);
  // h_EvMEC->SetMarkerColor(kViolet);  
  hs->Add(h_EvMEC);
  //hs->Draw("HIST SAME C");

  //DIS
  TH2D *h_EvDIS = (TH2D*) f->Get("Ep_min_max_dis"); 
  //h_EvDIS->SetLineWidth(2);
  // h_EvDIS->SetMarkerColor(kGreen+2);
  hs->Add(h_EvDIS);
  // hs->SetOption("HIST");
   hs->Draw("nostack");
   // hs->Draw("nostack");
   //hs->GetXaxis()->SetRangeUser(0,1);
   //hs->GetYaxis()->SetRangeUser(0,1000000);
  /*  
TLine *l = new TLine(2.2,0.0,2.2,(gPad->GetUymax())); 
  l->SetLineWidth(4);
  l->SetLineColor(kGray);
  l->SetVertical();
  gPad->Modified();
   gPad->Update();
   l->Draw();
  */
  gPad->Modified();
  gPad->Update();
 
  auto legend = new TLegend(0.75, 0.75, 0.9, 0.9); //(bottom left corner x coordinate, bottom left corner y coordinate, top right corner x coordinate, top right corner y coordinate)
  legend->SetHeader("");
  //legend->SetTextFont(60);
  legend->AddEntry(h_EvQE,"QE");
  legend->AddEntry(h_EvRES,"RES");
  legend->AddEntry(h_EvMEC,"MEC");
  legend->AddEntry(h_EvDIS,"DIS");

  legend->Draw();  
  // Save the plot you make as an image file
  canvas->SaveAs("new code energy rec(Cal)");
  

}
