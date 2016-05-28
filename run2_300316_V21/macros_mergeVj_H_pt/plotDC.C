#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <cstdlib>
#include <iostream>
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
/*
//CMS_monoHbb_scaleF                   shape    -        -   1.0     1.0      1.0     -     - 
CMS_monoHbb_scaleR                   shape  -        -   1.0     1.0      1.0     -     - 
CMS_monoHbb_res_j                    shape    1.0   1.0    1.0     1.0     1.0    1.0  1.0  
CMS_monoHbb_scale_j                  shape    1.0   1.0    1.0     1.0     1.0    1.0  1.0  
CMS_monoHbb_metUnclusteredEn         shape    1.0   1.0    1.0     1.0     1.0    1.0  1.0  
CMS_monoHbb_ewk                      shape    -    -     1.0     1.0     -    -  -   
*/

//void drawPdf(TString syst="", TString process="")
void drawPdf()
{

  TCanvas c1;

  //  TFile *_file0 = TFile::Open("plotsmonoH_SF_shapesMET_ATLAS/MET1509999Signal1J.root");
    TFile *_file0 = TFile::Open("plotsmonoH_SF_shapesMET_CMS/MET1509999Signal1J.root");

  TH1F * nom = (TH1F*) _file0->Get("Zj");
  TH1F * pdf = (TH1F*) _file0->Get("ZjLF");

  
  TH1F * b1 = new TH1F("b1", "b1", 100, 0, 100);
  TH1F * b2 = new TH1F("b2", "b2", 100, 0, 100);
  TH1F * b3= new TH1F("b3", "b3", 100, 0, 100);



  nom->SetMarkerSize(1.5);
  nom->SetMarkerStyle(20);
  nom->SetMarkerColor(kBlack);
  nom->SetLineColor(kBlack);
  nom->SetLineStyle(1);

 nom->Draw("e1");
 nom->SetTitle("");
 nom->GetXaxis()->SetTitle("MET [GeV]");

 for(int i=9; i<100; i++){
   pdf = (TH1F*) _file0->Get(Form("Zj_CMS_monoHbb_pdf%dUp",i));
   if (pdf){
     std::cout << "pdf->GetBinContent(1) "<< pdf->GetBinContent(1) << std::endl;
     std::cout << "pdf->GetBinContent(2) "<< pdf->GetBinContent(2) << std::endl;
     std::cout << "pdf->GetBinContent(3) "<< pdf->GetBinContent(3) << std::endl;
     b1->Fill(pdf->GetBinContent(1)); 
     b2->Fill(pdf->GetBinContent(2)); 
     b3->Fill(pdf->GetBinContent(3)); 
     pdf->SetLineWidth(1.0);
     pdf->SetLineColor(kYellow);
     pdf->SetLineStyle(kDashed);
     pdf->SetFillColor(0);
     pdf->Draw("HIST SAME");
   }
 }

 pdf->SetBinContent(1, b1->GetMean());
 pdf->SetBinError(1, b1->GetRMS());

 pdf->SetBinContent(2, b2->GetMean());
 pdf->SetBinError(2, b2->GetRMS());

 pdf->SetBinContent(3, b3->GetMean());
 pdf->SetBinError(3, b3->GetRMS());


 std::cout << "pdf bin1 mean "<< b1->GetMean() << std::endl;
 std::cout << "pdf bin2 mean "<< b2->GetMean() << std::endl;
 std::cout << "pdf bin3 mean "<< b3->GetMean() << std::endl;


 pdf->SetLineWidth(2);
 pdf->SetLineStyle(1);
 pdf->SetLineColor(kBlue);
 pdf->SetMarkerSize(0);
 // pdf->SetMarkerStyle(20);
 pdf->Draw("HISTSAME E1");


 TLegend * leg = new TLegend(0.52, 0.62, 0.92, 0.92);
 leg->SetFillColor(0);
 leg->SetLineColor(0);
 leg->SetShadowColor(0);
 leg->SetTextFont(62);
 leg->SetTextSize(0.03);
 leg->AddEntry(nom, "nom Zj", "p");
 leg->AddEntry(pdf, "pdf mean and RMS", "l");

 leg->Draw("");
 c1.SaveAs("ZjPdf.pdf");



}

void drawSyst(TString syst="", TString process=""){
  TCanvas c1;
    TFile *_file0 = TFile::Open("plotsmonoH_SF_shapesMET_CMS/MET1509999Signal1J.root");
  // TFile *_file0 = TFile::Open("plotsmonoH_SF_shapesMET_ATLAS/MET1509999Signal1J.root");
  TH1F * nom = (TH1F*) _file0->Get(process);
  TH1F * systup = (TH1F*) _file0->Get(process+ "_" + syst + "Up");
  TH1F * systdown = (TH1F*) _file0->Get(process+ "_" + syst + "Down");

  if (syst == "CMS_monoHbb_ewk"){
   systup = (TH1F*) _file0->Get(process+ "_" + syst + "Down");
   systdown = (TH1F*) _file0->Get(process+ "_" + syst + "Up");
  }

  nom->SetMarkerSize(1.5);
  nom->SetMarkerStyle(20);
  // nom->SetLineStyle(1);
  nom->SetFillColor(0);
  nom->SetLineColor(kBlack);
  nom->SetMarkerColor(kBlack);
  nom->SetLineWidth(0);

  nom->Draw("e1");
  nom->SetTitle("");
  nom->GetXaxis()->SetTitle("MET [GeV]");

  systup->SetLineColor(kRed);
  systup->SetFillColor(0);
  systup->SetLineStyle(kDashed);
  systdown->SetLineColor(kBlue);
  systdown->SetFillColor(0);
   systdown->SetLineStyle(kDashed);
   systup->Draw("HIST SAME");
   systdown->Draw("HIST SAME");

 TLegend * leg = new TLegend(0.52, 0.62, 0.92, 0.92);
 leg->SetFillColor(0);
 leg->SetLineColor(0);
 leg->SetShadowColor(0);
 leg->SetTextFont(62);
 leg->SetTextSize(0.025);
 leg->AddEntry(nom,  process+ "nom", "p");
 leg->AddEntry(systup, process+syst+"Up", "l");
 leg->AddEntry(systdown, process+syst+"Down", "l");

 leg->Draw("");
 c1.SaveAs(process+syst+".pdf");



}


void plotDC(){
  gStyle->SetOptStat(0);
  gROOT->LoadMacro("tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()");
    
  TH1::SetDefaultSumw2(1);

  //  drawPdf();
 
  TString syst[6] = { "CMS_monoHbb_scaleF", "CMS_monoHbb_ewk", "CMS_monoHbb_scaleR", "CMS_monoHbb_res_j", "CMS_monoHbb_scale_j", "CMS_monoHbb_metUnclusteredEn"}; 

  for (int i=0; i<6; i++){
    drawSyst(syst[i], "Zj");
    drawSyst(syst[i], "Wj");
    drawSyst(syst[i], "TT");
    drawSyst(syst[i], "monoH_800");
  }
}



