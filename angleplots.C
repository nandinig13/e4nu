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

void angleplots()
{
    gStyle -> SetOptStat(false);
    string Filename = "anglesapapadop_UpdatedSchwingerRad_SuSav2_C12_2261GeV.root";
    TFile *f = new TFile(Filename.c_str());
    
    f->ls(); //print out list of histograms in file

   
    TH2D *hangPim = (TH2D*) f->Get("theta_phi_pim"); 
    TCanvas *canv = new TCanvas("e4nu", "e4nu", 900, 600);
    THStack *hs = new THStack("hs", "Theta vs phi for charged pions; cos(theta) ; Cos(theta), phi");
    gStyle->SetPalette(100);
    hs->Add(hangPim);

    TH2D *hangPip = (TH2D*) f->Get("theta_phi_pip"); 
    hs->Add(hangPip);


    hs->Draw("colz");
    hs->GetXaxis()->SetRangeUser(0,1);

    auto legend = new TLegend(0.75, 0.75, 0.9, 0.9); //(bottom left corner x coordinate, bottom left corner y coordinate, top right corner x coordinate, top right corner y coordinate)
    legend->SetHeader("");
  
    legend->AddEntry(hangPim,"Pi-");
    legend->AddEntry(hangPip,"Pi+");
    legend->Draw();  

    canv->SaveAs("angles.png");
}