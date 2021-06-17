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

void ploton()
{
    gStyle -> SetOptStat(false);
    string Filename = "workingapapadop_UpdatedSchwingerRad_SuSav2_C12_2261GeV.root";
    TFile *f = new TFile(Filename.c_str());
    
    f->ls(); //print out list of histograms in file

   
    TH2D *hEQE = (TH2D*) f->Get("Ep_min_max_qe,"); //QUASI-ELASTIC
    TCanvas *canv = new TCanvas("e4nu", "e4nu", 900, 600);
    THStack *hs = new THStack("hs", "Proton kinematic energies for events with 2 protons and 1 charged pion; Energy /GeV ; Event count");
    gStyle->SetPalette(100);
    hs->Add(hEQE);

    TH2D *hERES = (TH2D*) f->Get("Ep_min_max_res,"); //RESONANT
    hs->Add(hERES);

    TH2D *hEMEC = (TH2D*) f->Get("Ep_min_max_mec,"); //MEC
    hs->Add(hEMEC);

    TH2D *hEDIS = (TH2D*) f->Get("Ep_min_max_dis,"); //DEEP INELASTIC
    hs->Add(hEDIS);

    hs->Draw("colz");
    hs->GetXaxis()->SetRangeUser(0,1);

    auto legend = new TLegend(0.75, 0.75, 0.9, 0.9); //(bottom left corner x coordinate, bottom left corner y coordinate, top right corner x coordinate, top right corner y coordinate)
    legend->SetHeader("");
  
    legend->AddEntry(hEQE,"QE");
    legend->AddEntry(hERES,"RES");
    legend->AddEntry(hEMEC,"MEC");
    legend->AddEntry(hEDIS,"DIS");
    legend->Draw();  

    canv->SaveAs("Proton energies.png");
}