
#include <string>
using namespace std;

/* Run this from the root command line  with .x anjaliPlots.C */
void anjaliPlotsnew() // This name has to match the file name
{
  gStyle->SetOptStat(false); // Turn off that ugly stats box
  
  // Create a TFile object and connect it to your file
  string Filename  = "/unix/dune/e4nu/ngadhia/xtestapapadop_SuSav2_C12_1161GeV.root.root";
  TFile *f = new TFile(Filename.c_str());
  // "f" is just a placeholder name at this point
 
 
  f->ls(); // This will print out a list of histograms in your file
  
  // Get the hEv_qe histogram and store it in a TH1D object
  // (1-dimensional histogram of doubles i.e. high-precision decimals
  
   TH1D *h_EvQE = (TH1D*) f->Get("hcal_qe"); // Quasi-elastic
  // TH1D *h_EvRES = (TH1D*) f->Get("hEv_res"); // Resonant
  // TH1D *h_EvMEC = (TH1D*) f->Get("hEv_mec"); // Meson exchange currents
  // TH1D *h_EvDIS = (TH1D*) f->Get("hEv_dis"); // Deep inelastic scattering  
  // Create a canvas to draw on
  TCanvas *canvas = new TCanvas("e4nu","e4nu",900,600); // We will overwrite the default titles. The numbers set the size
  THStack *hs = new THStack("hs", "Calorimetrically reconstructed energy for e^{-}-{}^{12}C scattering at 1.161GeV; E_{cal} /GeV ; Event count");
  

  //Quasi-elastic
  // TH1D *h_EvQE = (TH1D*) f->Get("hEQE_qe"); // Quasi-elastic
  h_EvQE->GetYaxis()->SetTitle("Event count");
  //h_EvQE->GetYaxis()->SetRangeUser(0,1700000);
  // h_EvQE->GetXaxis()->SetRangeUser(0,4); //Go up to 10Gev 
  h_EvQE->GetXaxis()->SetTitle("E_{QE} /GeV"); 
  h_EvQE->SetTitle("Quasi-elastic reconstruction of electron energy for #e^- ^{12}C# at 2.261GeV (SuSav2).");
  h_EvQE->SetLineWidth(2);
  h_EvQE->SetLineColor(kRed);
  hs->Add(h_EvQE);
  
  //Resonant
  TH1D *h_EvRES = (TH1D*) f->Get("hcal_res"); 
   h_EvRES->SetLineWidth(2);
   h_EvRES->SetLineColor(kBlue);
   hs->Add(h_EvRES);

  //Meson-exchange
  TH1D *h_EvMEC = (TH1D*) f->Get("hcal_mec");
  h_EvMEC->SetLineWidth(2);
  h_EvMEC->SetLineColor(kMagenta+2);  
  hs->Add(h_EvMEC);

  //DIS
  TH1D *h_EvDIS = (TH1D*) f->Get("hcal_dis"); 
  h_EvDIS->SetLineWidth(2);
  h_EvDIS->SetLineColor(kGreen+3);
  hs->Add(h_EvDIS);

  //True energy
  //TH1D *tru = (TH1D*)f->Get("hEv_qe");
  //tru->SetLineWidth(2);
  // tru->SetLineColor(kGray);
  //hs->Add(tru);
  // Draw the QE histogram
  //h_EvQE->Draw("HIST"); // HIST just draws a line
  // And now draw the RES one on the same axes
  //h_EvRES->Draw("HIST SAME");
  // It uses the first one you draw to set the range of the histogram
  // What will you do if a later one goes to a higher value?!
  // h_EvMEC->Draw("HIST SAME");
  
  //h_EvDIS->Draw("HIST SAME");
  hs->Draw("nostack");
  hs->GetXaxis()->SetRangeUser(0,2);
  //hs->GetYaxis()->SetRangeUser(0,1700000);
  TLine *l = new TLine(1.1,0,1.1,8000000); 
  l->SetLineWidth(4);
  l->SetLineColor(kGray);
  l->SetVertical();
  gPad->Modified();
   gPad->Update();
   l->Draw();

  // Add a legend using TLegend https://root.cern.ch/doc/master/classTLegend.html

  auto legend = new TLegend(0.75, 0.75, 0.9, 0.9); //(bottom left corner x coordinate, bottom left corner y coordinate, top right corner x coordinate, top right corner y coordinate)
  legend->AddEntry(h_EvQE,"QE");
  legend->AddEntry(h_EvRES,"RES");
  legend->AddEntry(h_EvMEC,"MEC");
  legend->AddEntry(h_EvDIS,"DIS");
  // legend->AddEntry(tru, "TRUE BEAM ENERGY");
  legend->Draw();  
  //trying to print number of entries
  ((TH1*)(hs -> GetStack() -> Last())) -> GetEntries();
  // Save the plot you make as an image file
  canvas->SaveAs("carbon 1.1 energiesCAL.png");
  

}
