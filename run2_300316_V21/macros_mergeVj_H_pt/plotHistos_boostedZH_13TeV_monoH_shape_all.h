#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <vector>

#include "TSystem.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "THStack.h"
#include "TCut.h"
#include "TString.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TMath.h"

//#include "XSec_8TeV19invfb.h"
//#include "XSec_13TeV20invfb.h"

#include "TLorentzVector.h"
#include "HelperBDTShape.h"




// weights correction for EWK NLO correction
double ptWeightZllH(int nGenVbosons,double GenVbosons_pt,int VtypeSim,int GenVbosons_pdgId){
  double SF = 1.;
  if (nGenVbosons ==1)
    {
      if (VtypeSim == 0 || VtypeSim == 1 || VtypeSim == 4 || VtypeSim == 5)
	{
	  if (GenVbosons_pdgId == 23)
	    {
	      //for Z options
	      if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.1808051+6.04146*(TMath::Power((GenVbosons_pt+759.098),-0.242556));
	    }
	}
      else if (VtypeSim == 2 || VtypeSim == 3)
	{
	  //for W options
	  if (GenVbosons_pdgId == 24 || GenVbosons_pdgId == -24)
	    {
	      if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.830041+7.93714*(TMath::Power((GenVbosons_pt+877.978),-0.213831));
	    }
	}
    }
  //std::cout << " SF " << " " << SF << std::endl;
  return SF>0?SF:0;
}
/*
values["WJetsHT100"      ] =     1.23 *  1.29e3 ;
values["WJetsHT200"      ] =     1.23 *  3.86e2 ;
values["WJetsHT400"      ] =      1.23 * 47.9 ;
values["WJetsHT600"      ] =       1.23 * 19.9 ;
values["WJetsIncl"       ] =    61623.000000 ;
values["ZJetsHT100"      ] =      409.860000 ;
values["ZJetsHT200"      ] =      110.880000 ;
values["ZJetsHT400"      ] =         13.189  ;
values["ZJetsHT600"      ] =        4.524300 ;
instead use the number here https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns

*/


double ptWeightQCD(double lheHT, int GenVbosons_pdgId){
  double SF = 1.;
  if (lheHT>100){

    if (GenVbosons_pdgId == 23){ // Z
    SF =   ((lheHT>100 && lheHT<200)*1.588 * (  280.35 / (409.860000) ) + (lheHT>200 && lheHT<400)*1.438 * ( 77.67 / ( 110.880000 )) + (lheHT>400 && lheHT<600)*1.494 * (10.73 / (13.189 )) + (lheHT>600)*1.139 * ( 4.116 / (4.524300) ));


    }
  

    if (abs(GenVbosons_pdgId) == 24){
      SF =   ((lheHT>100 && lheHT<200)*1.588 * (  1345 / (1.23 *  1.29e3) ) + (lheHT>200 && lheHT<400)*1.438 * ( 359.7 / ( 1.23 *  3.86e2 )) + (lheHT>400 && lheHT<600)*1.494 * (48.91 / (1.23 * 47.9 )) + (lheHT>600)*1.139 * ( 18.77 / (1.23 * 19.9) ));
      
    }
  }
  
  return SF>0?SF:0;
}


TH1F * DrawOverflow(TH1F *h)
{
  //    h->Sumw2();
  // This function paint the histogram h with an extra bin for overflows
  UInt_t nx    = h->GetNbinsX()+1;
  Double_t *xbins= new Double_t[nx+1];
  for (UInt_t i=0;i<nx;i++)
  xbins[i]=h->GetBinLowEdge(i+1);
  xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
  //  char *tempName= new char[strlen(h->GetName())+10];
  char *tempName= new char[strlen(h->GetName())];
  // sprintf(tempName,"%swtOverFlow",h->GetName());
  // Book a temporary histogram having ab extra bin for overflows
  //  TH1F *htmp = new TH1F(tempName, h->GetTitle(), nx, xbins);
  TH1F *htmp = new TH1F( h->GetName(), h->GetTitle(), nx, xbins);
  // Reset the axis labels
  
   htmp->Sumw2();
  htmp->SetXTitle(h->GetXaxis()->GetTitle());
  htmp->SetYTitle(h->GetYaxis()->GetTitle());
  // Fill the new hitogram including the extra bin for overflows
   for (UInt_t i=1; i<=nx; i++){
     //    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
     htmp->SetBinContent(i , h->GetBinContent(i));
     htmp->SetBinError( i , h->GetBinError(i));
    //   std::cout << "h->GetBinContent(i)) " << h->GetBinContent(i) << " and  h->GetBinError(i)" <<  h->GetBinError(i) << std::endl;
  }
  // Fill the underflows
   //  htmp->Fill(h->GetBinLowEdge(1)-1, h->GetBinContent(0));
   htmp->SetBinContent(0 , h->GetBinContent(0));
   htmp->SetBinError(0, h->GetBinError(0));
  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());
  // FillStyle and color
  htmp->SetFillStyle(h->GetFillStyle());
  htmp->SetFillColor(h->GetFillColor());

  
  return htmp;
}


// _____________________________________________________________________________
// This class reads the ntuples and applies basic cuts to reduce loading time.
class Events {
public:
    Events();
    ~Events();

    void read(TCut cutmc_all, TCut cutdata_all, TString processes="", TString treename = "tree");

    void set_sf(const double sf[5]) {
        sf_WjLF = sf[0];
        sf_WjHF = sf[1];
        sf_ZjLF = sf[2];
        sf_ZjHF = sf[3];
        sf_TT   = sf[4];
        return;
    }

    // Reduced TTrees
    TTree * ZH;
    TTree * ggZH;
    TTree * boostedZH;
    TTree * boostedZH_2;
    TTree * boostedZH_3;
     TTree * WH;
    TTree * WjLF;
    TTree * WjHF;
    TTree * ZjLF;
    TTree * ZjHF;

    TTree * ZjLF_HT100;
    TTree * ZjHF_HT100;
    TTree * ZjLF_HT200;
    TTree * ZjHF_HT200;
    TTree * ZjLF_HT600;
    TTree * ZjHF_HT600;


    TTree * TT;
    TTree * ST_T_s;
    TTree * ST_T_t;
    TTree * ST_T_tW;
    TTree * ST_Tbar_s;
    TTree * ST_Tbar_t;
    TTree * ST_Tbar_tW;
    TTree * s_Top;
    TTree * VV;
    TTree * VV_LF;
    TTree * VV_HF;
    //    TTree * VV_ZZ;
    TTree * QCD;
    TTree * data_obs;


    // Scale factors
    float sf_WjLF;
    float sf_WjHF;
    float sf_ZjLF;
    float sf_ZjHF;
    float sf_TT;

    int massH;  // Higgs mass

};  // Events

// _____________________________________________________________________________
// These are functions to set the styles of TH1 and TLegend, respectively
void set_style(TH1 * h, const TString& p);
void set_style_norm(TH1 * h, const TString& p);
void set_style(TLegend * leg);

// _____________________________________________________________________________
// This function makes a single plot
void MakePlot(TTree * tree, TString var, TCut cut,
              TString title, int nbinsx, double xlow, double xup,
              TString plotname="plot", TString plotdir="",
              TString options="plotLog:plotNorm") {

    std::clog << "MakePlots(): Plot var: " << var << std::endl;
    std::clog << "MakePlots(): Using cut: " << cut << std::endl;

    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));
    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    TH1F * h1         = new TH1F("h1"       , title, nbinsx, xlow, xup);
    tree->Project("h1", var, cut);

    TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
    c1->SetLogy(plotLog);
    h1->SetLineWidth  (2);
    h1->SetLineColor  (2);
    h1->SetMarkerSize (0);
    h1->GetYaxis()->SetTitle("Events");
    if (plotNorm) {
        h1->Scale(1.0 / h1->Integral());
        h1->GetYaxis()->SetTitle("arbitrary unit");
    }

    h1->Draw("hist");
    h1->Draw("same e1");

    // Setup legends
    //    TLegend * leg1 = new TLegend(0.50, 0.68, 0.72, 0.92);
    //set_style(leg1);
    //leg1->AddEntry(h1, tree->GetUserInfo()->At(0)->GetName() , "l");
    //    leg1->Draw();

    TLatex * latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(12);
    latex->SetTextFont(62);
    latex->SetTextSize(0.04);
    //latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
    latex->DrawLatex(0.19, 0.89, "CMS 2015");
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 13 TeV");


    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");

    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    return;
}

// _____________________________________________________________________________
// This function makes two plots overlaid
void MakePlot2(TTree * tree1, TTree * tree2, TString var, TCut cut,
               TString title, int nbinsx, double xlow, double xup,
               TString plotname="plot", TString plotdir="",
               TString options="plotLog:plotNorm") {

    std::clog << "MakePlots(): Plot var: " << var << std::endl;
    std::clog << "MakePlots(): Using cut: " << cut << std::endl;

    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));
    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    TH1F * h1         = new TH1F("h1"       , title, nbinsx, xlow, xup);
    TH1F * h2         = new TH1F("h2"       , title, nbinsx, xlow, xup);
    tree1->Project("h1", var, cut);
    tree2->Project("h2", var, cut);

    TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
    c1->SetLogy(plotLog);
    h1->SetLineWidth  (2);
    h1->SetLineColor  (2);
    h1->SetMarkerSize (0);
    h1->GetYaxis()->SetTitle("Events");
    h2->SetLineWidth  (2);
    h2->SetLineColor  (4);
    h2->SetMarkerSize (0);
    h2->GetYaxis()->SetTitle("Events");
    if (plotNorm) {
        h1->Scale(1.0 / h1->Integral());
        h1->GetYaxis()->SetTitle("arbitrary unit");
        h2->Scale(1.0 / h2->Integral());
        h2->GetYaxis()->SetTitle("arbitrary unit");
    }

    h1->Draw("hist");
    h1->Draw("same e1");
    h2->Draw("same hist");
    h2->Draw("same e1");


    // Setup legends
    //TLegend * leg1 = new TLegend(0.50, 0.68, 0.72, 0.92);
    //set_style(leg1);
    //leg1->AddEntry(h1, tree1->GetUserInfo()->At(0)->GetName() , "l");
    //leg1->AddEntry(h2, tree2->GetUserInfo()->At(0)->GetName() , "l");
    //    leg1->Draw();
    

    TLatex * latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(12);
    latex->SetTextFont(62);
    latex->SetTextSize(0.04);
    //latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
    latex->DrawLatex(0.19, 0.89, "CMS Simulation 2015");
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 13 TeV");


    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");



    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    //delete h2;
    return;
}

// _____________________________________________________________________________
// This function makes stacked plots
void MakePlots(const Events * ev, TString var,
               TCut cutmc, TCut cutdata,
               TString title, int nbinsx, double xlow, double xup,
               TString plotname="plot", TString plotdir="",
               TString options="printStat:plotSig:plotData:plotLog:!plotNorm:!plotBoostedZH", bool doRebin = false , int newnbinsx  = -1 , bool drawOverflow = false, TString CRname = "") {

   
  if ( newnbinsx == -1) newnbinsx = nbinsx;

    std::clog << "MakePlots(): Plot var: " << var << std::endl;
    std::clog << "MakePlots(): Using cutmc: " << cutmc << ", cutdata: " << cutdata << std::endl;

    // Parse options
    bool printStat    = options.Contains("printStat")  && (!options.Contains("!printStat"));
    //bool printCard    = options.Contains("printCard")  && (!options.Contains("!printCard"));
    //bool plotUOflow   = options.Contains("plotUOflow") && (!options.Contains("!plotUOflow"));  // FIXME: not yet implemented
    bool plotSig      = options.Contains("plotSig")    && (!options.Contains("!plotSig"));
    bool plotBoostedZH      = options.Contains("plotBoostedZH")    && (!options.Contains("!plotBoostedZH"));
    bool plotData     = options.Contains("plotData")   && (!options.Contains("!plotData"));
    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));

    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    // Book histograms
    TH1F * hZH        = new TH1F("ZH"       , title, nbinsx, xlow, xup);
    TH1F * hboostedZH        = new TH1F("boostedZH"       , title, nbinsx, xlow, xup);
    TH1F * hboostedZH_2        = new TH1F("boostedZH_2"       , title, nbinsx, xlow, xup);
    TH1F * hboostedZH_3        = new TH1F("boostedZH_3"       , title, nbinsx, xlow, xup);
    TH1F * hggZH        = new TH1F("ggZH"       , title, nbinsx, xlow, xup);
    TH1F * hWH        = new TH1F("WH"       , title, nbinsx, xlow, xup);
    TH1F * hWjLF      = new TH1F("WjLF"     , title, nbinsx, xlow, xup);
    TH1F * hWjHF      = new TH1F("WjHF"     , title, nbinsx, xlow, xup);
    TH1F * hZjLF      = new TH1F("ZjLF"     , title, nbinsx, xlow, xup);
    TH1F * hZjHF      = new TH1F("ZjHF"     , title, nbinsx, xlow, xup);
    TH1F * hTT        = new TH1F("TT"       , title, nbinsx, xlow, xup);
    TH1F * hST        = new TH1F("ST"       , title, nbinsx, xlow, xup);
    TH1F * hVV        = new TH1F("VV"       , title, nbinsx, xlow, xup);
    TH1F * hQCD        = new TH1F("QCD"       , title, nbinsx, xlow, xup);
    TH1F * hdata_obs  = new TH1F("data_obs" , title, nbinsx, xlow, xup);
    TH1F * hVH        = new TH1F("VH"       , title, nbinsx, xlow, xup);
    TH1F * hmc_exp    = new TH1F("mc_exp"   , title, nbinsx, xlow, xup);

    std::cout << "plotSig is " << plotSig <<std::endl;

    // Draw histograms, apply cut and MC event weights (= equiv_lumi * PUweight)  // no PUweight for now...
    if (plotSig) {
      //        ev->ZH->Project("ZH", var, cutmc * Form("%5f * 1.", ev->lumi_ZH));

              ev->ZH->Project("ZH", var, cutmc );
      

  std::clog << "... DONE: project ZH." << std::endl;

        ev->ggZH->Project("ggZH", var, cutmc );
        std::clog << "... DONE: project ggZH." << std::endl;

	//      ev->WH->Project("WH", var, cutmc * Form("%5f * 1.", ev->lumi_WH));
		     ev->WH->Project("WH", var, cutmc ); 
        std::clog << "... DONE: project WH." << std::endl;

    }

    if (plotBoostedZH) {
	//        ev->boostedZH->Project("boostedZH", var, cutmc * Form("%5f * 1.", ev->lumi_boostedZH));
        std::clog << "... before project boostedZH." << std::endl;
	 		//TCut effLumiboostedVH = "(3000. * 1.  * ( 0.577 * 0.1 )/ 50000.)";    //mZ'=1000, assume 1 pb and scale to 3 fb-1
        ev->boostedZH->Project("boostedZH", var, cutmc );
        std::clog << "... DONE: project boostedZH." << std::endl;

        ev->boostedZH_2->Project("boostedZH_2", var, cutmc );
        std::clog << "... DONE: project boostedZH_2." << std::endl;


        ev->boostedZH_3->Project("boostedZH_3", var, cutmc );
        std::clog << "... DONE: project boostedZH_3" << std::endl;

    }



    //    ev->WjLF->Project("WjLF", var, cutmc * Form("%5f * 1.", ev->lumi_WjLF) * Form("%f", ev->sf_WjLF));
    ev->WjLF->Project("WjLF", var, ( cutmc  * Form("%f", ev->sf_WjLF)) *  "ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId) * ptWeightQCD(lheHT, GenVbosons_pdgId)" );
    std::clog << "... DONE: project WjLF, scaled by " << ev->sf_WjLF << "." << std::endl;

    //    ev->WjLF->Project("WjHF", var, cutmc * Form("%5f * 1", ev->lumi_WjHF) * Form("%f", ev->sf_WjHF));
    ev->WjHF->Project("WjHF", var, cutmc  * Form("%f", ev->sf_WjHF) * "ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId) * ptWeightQCD(lheHT, GenVbosons_pdgId)" );
    std::clog << "... DONE: project WjHF, scaled by " << ev->sf_WjHF << "." << std::endl;

    //    ev->ZjLF->Project("ZjLF", var, cutmc * Form("%5f * 1", ev->lumi_ZjLF) * Form("%f", ev->sf_ZjLF));
    ev->ZjLF->Project("ZjLF", var,   cutmc  * Form("%f", ev->sf_ZjLF) * "ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId) * ptWeightQCD(lheHT,GenVbosons_pdgId)");
    std::clog << "... DONE: project ZjLF, scaled by " << ev->sf_ZjLF << "." << std::endl;

     //    ev->ZjHF->Project("ZjHF", var, cutmc * Form("%5f * 1", ev->lumi_ZjHF) * Form("%f", ev->sf_ZjHF));
    ev->ZjHF->Project("ZjHF", var,  cutmc  * Form("%f", ev->sf_ZjHF) * "ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId)* ptWeightQCD(lheHT,GenVbosons_pdgId)");
    std::clog << "... DONE: project ZjHF, scaled by " << ev->sf_ZjHF << "." << std::endl;    
    
    /*
    ev->ZjLF_HT100->Project("ZjLF", var, cutmc * Form("%5f * 1", ev->lumi_ZjLFHT100) * Form("%f", ev->sf_ZjLF)); 
    ev->ZjLF_HT200->Project("+ZjLF", var, cutmc * Form("%5f * 1", ev->lumi_ZjLFHT200) * Form("%f", ev->sf_ZjLF)); 
    ev->ZjLF_HT600->Project("+ZjLF", var, cutmc * Form("%5f * 1", ev->lumi_ZjLFHT600) * Form("%f", ev->sf_ZjLF)); 


    ev->ZjHF_HT100->Project("ZjHF", var, cutmc * Form("%5f * 1", ev->lumi_ZjHFHT100) * Form("%f", ev->sf_ZjHF)); 
    ev->ZjHF_HT200->Project("+ZjHF", var, cutmc * Form("%5f * 1", ev->lumi_ZjHFHT200) * Form("%f", ev->sf_ZjHF)); 
    ev->ZjHF_HT600->Project("+ZjHF", var, cutmc * Form("%5f * 1", ev->lumi_ZjHFHT600) * Form("%f", ev->sf_ZjHF)); 
    */


    //    ev->TT->Project("TT", var, cutmc * Form("%5f * 1.", ev->lumi_TT) * Form("%f", ev->sf_TT));
    ev->TT->Project("TT", var, cutmc  * Form("%f", ev->sf_TT));
    std::clog << "... DONE: project TT, scaled by " << ev->sf_TT << "." << std::endl;


    /*
    //    ev->ST_T_s    ->Project( "ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_T_s));
    ev->ST_T_s    ->Project( "ST" , var, cutmc );

    //    ev->ST_T_t    ->Project("+ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_T_t));
    ev->ST_T_t    ->Project("+ST" , var, cutmc );

    //    ev->ST_T_tW   ->Project("+ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_T_tW));
    ev->ST_T_tW   ->Project("+ST" , var, cutmc );

    //    ev->ST_Tbar_s ->Project("+ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_Tbar_s));
    ev->ST_Tbar_s ->Project("+ST" , var, cutmc );


    //    ev->ST_Tbar_t ->Project("+ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_Tbar_t));
    ev->ST_Tbar_t ->Project("+ST" , var, cutmc );

    //    ev->ST_Tbar_tW->Project("+ST" , var, cutmc * Form("%5f * 1.", ev->lumi_ST_Tbar_tW));
    ev->ST_Tbar_tW->Project("+ST" , var, cutmc );
    */
    ev->s_Top->Project("ST" , var, cutmc );                                                                                                                                    



    std::clog << "... DONE: project ST." << std::endl;


    ev->VV->Project( "VV", var, cutmc );

    /*
    ev->VV_WW->Project( "VV", var, cutmc * Form("%5f * 1.", ev->lumi_VV_WW));
    ev->VV_WZ->Project("+VV", var, cutmc * Form("%5f * 1.", ev->lumi_VV_WZ));
    ev->VV_ZZ->Project("+VV", var, cutmc * Form("%5f * 1.", ev->lumi_VV_ZZ));
    */

    std::clog << "... DONE: project VV." << std::endl;


    ev->QCD ->Project("QCD" , var, cutmc );
    std::clog << "... DONE: project QCD." << std::endl;


    if (plotData) {
        ev->data_obs->Project("data_obs", var, cutdata);
    }
    std::clog << "... DONE: project data_obs." << std::endl;

    // Sum histograms
    hVH->Sumw2();
    hVH->Add(hZH);
    hVH->Add(hWH); 
    hVH->Add(hggZH);
    std::clog << "... DONE: add ZH, WH to VH." << std::endl;
 

    // 
    hTT->Add(hST);

    // add LF and HF
    hZjHF->Add(hZjLF);
    hWjHF->Add(hWjLF);


    if (plotNorm) {
       
        TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
	c1->SetLogy(plotLog);


        hVH->Scale(1.0 / hVH->Integral());
        hVH->GetYaxis()->SetTitle("arbitrary unit");
        hggZH->Scale(1.0 / hggZH->Integral());
        hggZH->GetYaxis()->SetTitle("arbitrary unit");
        hboostedZH->Scale(1.0 / hboostedZH->Integral());
        hboostedZH->GetYaxis()->SetTitle("arbitrary unit");


        hboostedZH_2->Scale(1.0 / hboostedZH_2->Integral());
        hboostedZH_2->GetYaxis()->SetTitle("arbitrary unit");


        hboostedZH_3->Scale(1.0 / hboostedZH_3->Integral());
        hboostedZH_3->GetYaxis()->SetTitle("arbitrary unit");


        hST->Scale(1.0 / hST->Integral());
        hST->GetYaxis()->SetTitle("arbitrary unit");
        hTT->Scale(1.0 / hTT->Integral());
        hTT->GetYaxis()->SetTitle("arbitrary unit");
        hZjLF->Scale(1.0 / hZjLF->Integral());
        hZjLF->GetYaxis()->SetTitle("arbitrary unit");
        hZjHF->Scale(1.0 / hZjHF->Integral());
        hZjHF->GetYaxis()->SetTitle("arbitrary unit");
        hWjLF->Scale(1.0 / hWjLF->Integral());
        hWjLF->GetYaxis()->SetTitle("arbitrary unit");
        hWjHF->Scale(1.0 / hWjHF->Integral());
        hWjHF->GetYaxis()->SetTitle("arbitrary unit");
        hVV->Scale(1.0 / hVV->Integral());
        hVV->GetYaxis()->SetTitle("arbitrary unit");
        hQCD->Scale(1.0 / hQCD->Integral());
        hQCD->GetYaxis()->SetTitle("arbitrary unit");

    set_style_norm(hVH, "VH");
    set_style_norm(hZH, "VH");
    set_style_norm(hggZH, "ggZH");
    set_style_norm(hboostedZH, "boostedZH");

    set_style_norm(hboostedZH_2, "boostedZH_2");

    set_style_norm(hboostedZH_3, "boostedZH_3");

    set_style_norm(hWH, "VH");
    set_style_norm(hWjLF, "WjLF");
    set_style_norm(hWjHF, "WjHFb");
    set_style_norm(hZjLF, "ZjLF");
    set_style_norm(hZjHF, "ZjHFb");
    set_style_norm(hTT, "TT");
    set_style_norm(hST, "ST");
    set_style_norm(hVV, "VVHF");
    set_style_norm(hQCD, "QCD");


	double ymax = TMath::Max( hVH->GetMaximum(),1.);
	hVH->SetMaximum(ymax * 1.7 + (ymax>1 ? sqrt(ymax) : 0.));
	if (plotLog)  hVH->SetMaximum(ymax * 200 + (ymax>1 ? sqrt(ymax) : 0.));


    
	hVH->Draw("hist");
	hVH->Draw("same e1");

	//	values["ZnnH125"         ] =      (0.8696 -  0.1057) * 0.577 * 0.2 ;
	/*	values["monoH_600"           ] =         0.026 ;
	values["monoH_800"           ] =        0.0288 ;
	values["monoH_1000"           ] =       0.02337 ;
	*/

	// scale for new CMS/ATLAS numbers
	//	hboostedZH->Scale( 42.386 * 0.577 / 0.026 );
	//	hboostedZH->Scale( 372.2 * 0.577 / 0.026 );


	hboostedZH->Draw("same hist");
	hboostedZH->Draw("same e1");

	//	hboostedZH_2->Scale( 45.097 * 0.577 / 0.0288 );
	//hboostedZH_2->Scale( 230.67 * 0.577 / 0.0288 );


	hboostedZH_2->Draw("same hist");
	hboostedZH_2->Draw("same e1");

	//	hboostedZH_2->Scale( 35.444 * 0.577 / 0.02337 );
	//hboostedZH_2->Scale( 119.34 * 0.577 / 0.02337 );

	hboostedZH_3->Draw("same hist");
	hboostedZH_3->Draw("same e1");

	/*
	hggZH->Draw("same hist");
	hggZH->Draw("same e1");
	*/
	/*
	hST->Draw("same hist");
	hST->Draw("same e1");
	hTT->Draw("same hist");
	hTT->Draw("same e1");
	*/

	//	hZjLF->Draw("same hist");
	//hZjLF->Draw("same e1");
	hZjHF->Draw("same hist");
	hZjHF->Draw("same e1");
	//hWjLF->Draw("same hist");
	//hWjLF->Draw("same e1");
	hWjHF->Draw("same hist");
	hWjHF->Draw("same e1");
	hVV->Draw("same hist");
	hVV->Draw("same e1");
	hQCD->Draw("same hist");
	hQCD->Draw("same e1");

	// Setup legends
	TLegend * leg1 = new TLegend(0.50, 0.68, 0.72, 0.92);
	set_style(leg1);
	if (plotData) leg1->AddEntry(hdata_obs, "Data", "p");
	if (plotSig)  {
	  leg1->AddEntry(hVH, Form("VH(%i)", ev->massH), "f");
	  //	  leg1->AddEntry(hggZH, Form("ggZH(%i)", ev->massH), "f");
	}

	if (plotBoostedZH)  leg1->AddEntry(hboostedZH, "Z->ZH(MZ'=600), xsec * BR = 1pb", "l");

	if (plotBoostedZH)  leg1->AddEntry(hboostedZH_2, "Z->ZH(MZ'=800), xsec * BR = 1pb", "l");

	if (plotBoostedZH)  leg1->AddEntry(hboostedZH_3, "Z->ZH(MZ'=1000), xsec * BR = 1pb", "l");

	leg1->AddEntry(hTT, "top", "f");
	//	leg1->AddEntry(hST, "single top", "f");
	leg1->AddEntry(hVV, "VV", "f");

	TLegend * leg2 = new TLegend(0.72, 0.68, 0.94, 0.92);
	set_style(leg2);
	leg2->AddEntry(hWjHF, "Wj", "f");
	//	leg2->AddEntry(hWjLF, "W + LF", "f");
	leg2->AddEntry(hZjHF, "Zj", "f");
	//leg2->AddEntry(hZjLF, "Z + LF", "f");
	leg2->AddEntry(hQCD, "QCD", "f");


        leg1->Draw("same");
        leg2->Draw("same");
 
	TLatex * latex = new TLatex();
	latex->SetNDC();
	latex->SetTextAlign(12);
	latex->SetTextFont(62);
	latex->SetTextSize(0.035);
	//latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
	latex->DrawLatex(0.19, 0.89, "CMS Simulation 2015");
	latex->SetTextSize(0.035);
	latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 13 TeV");
	
	
	gPad->Print(plotdir+plotname+".png");
	gPad->Print(plotdir+plotname+".pdf");


    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    //delete h2;
	return;

   
    }



       

    hmc_exp->Sumw2();
    //    hmc_exp->Add(hWjLF);
    hmc_exp->Add(hWjHF);
    //  hmc_exp->Add(hZjLF);
    hmc_exp->Add(hZjHF);
    hmc_exp->Add(hTT);
    // hmc_exp->Add(hST);
    hmc_exp->Add(hVV);
    hmc_exp->Add(hQCD);
    hmc_exp->Add(hVH);

    std::clog << "... DONE: add all backgrounds to mc_exp." << std::endl;


    if (drawOverflow) {

      std::cout << "drawing overflow ################################" << std::endl;  
       

        hVH        = DrawOverflow(hVH ); // fixme
        hZH        = DrawOverflow(hZH             );  // fixme
        hWH        = DrawOverflow(hWH             ); //fixme
        hboostedZH = DrawOverflow(hboostedZH      );

        hboostedZH_2 = DrawOverflow(hboostedZH_2      );

        hboostedZH_3 = DrawOverflow(hboostedZH_3      );


	hggZH      = DrawOverflow(hggZH           );  // fixme
        hWjLF      = DrawOverflow(hWjLF           );
        hWjHF      = DrawOverflow(hWjHF           );
        hZjLF      = DrawOverflow(hZjLF           );
        hZjHF      = DrawOverflow(hZjHF           );
        hTT        = DrawOverflow(hTT             );
        hST        = DrawOverflow(hST             );
        hVV        = DrawOverflow(hVV             ) ;
        hQCD       = DrawOverflow(hQCD            );
      
        hmc_exp    = DrawOverflow(hmc_exp         );
        hdata_obs  = DrawOverflow(hdata_obs       );

    }




    /// Declare rebinner                                                                                                                                                                            

    if (doRebin) {
      Rebinner* rebinner = 0;
      rebinner = new Rebinner( newnbinsx , 0.25, 0.25,  xlow, xup);
      rebinner->set_signal_backgr(hboostedZH, hmc_exp);


      std::cout << "rebinning ################################" << std::endl;  
      

        hVH        = rebinner->rebin(hVH         , newnbinsx  , "VH"        ); // fixme
        hZH        = rebinner->rebin(hZH         , newnbinsx  , "ZH"        );  // fixme
        hWH        = rebinner->rebin(hWH         , newnbinsx  , "WH"        ); //fixme
        hboostedZH        = rebinner->rebin(hboostedZH         , newnbinsx  , "boostedZH"        );

        hboostedZH_2        = rebinner->rebin(hboostedZH_2         , newnbinsx  , "boostedZH_2"        );

        hboostedZH_3        = rebinner->rebin(hboostedZH_3         , newnbinsx  , "boostedZH_3"        );

	hggZH        = rebinner->rebin(hggZH     , newnbinsx  ,"ggZH"        );  // fixme
        hWjLF      = rebinner->rebin(hWjLF       , newnbinsx ,"WjLF"      );
        hWjHF      = rebinner->rebin(hWjHF       , newnbinsx ,"WjHF"      );
        hZjLF      = rebinner->rebin(hZjLF       , newnbinsx ,"ZjLF"      );
        hZjHF      = rebinner->rebin(hZjHF       , newnbinsx ,"ZjHF"      );
        hTT        = rebinner->rebin(hTT         , newnbinsx ,"TT"        );
        hST     = rebinner->rebin(hST      , newnbinsx ,"ST"     );
        hVV      = rebinner->rebin(hVV       , newnbinsx  ,"VV"      );
        hQCD       = rebinner->rebin(hQCD        , newnbinsx  ,"QCD"       );
        hmc_exp    = rebinner->rebin(hmc_exp     , newnbinsx ,"mc_exp"    );
        hdata_obs  = rebinner->rebin(hdata_obs   , newnbinsx ,"data_obs"  );
    }







    // Setting up histograms (boring stuff below)
    std::clog << "MakePlots(): Setting up histograms..." << std::endl;

    // Setup histogram stack
    THStack * hs = new THStack("hs", "");
    //    hs->Add(hST);
    hs->Add(hVV);
    hs->Add(hTT);
    //    hs->Add(hZjLF);
    hs->Add(hZjHF);
    // hs->Add(hWjLF);
    hs->Add(hWjHF);
    hs->Add(hQCD);
    if (plotSig)  hs->Add(hVH);
    //  if (plotSig)  hs->Add(hggZH);

    //    if (plotBoostedZH)  hs->Add(hboostedZH);




    if (printStat) {
        double nS = hVH->Integral();
      double nS_err = -1; 
      if (plotBoostedZH) {
	 nS = hboostedZH->IntegralAndError(1,hboostedZH->GetNbinsX(), nS_err);
      } else {
	nS = hVH->IntegralAndError(1,hVH->GetNbinsX(), nS_err);
	//nggZHS = hggZH->IntegralAndError(1,hggZH->GetNbinsX(), nggZHS_err);
      }
      double nVH_err = -1;
      double nVH = hVH->IntegralAndError(1,hVH->GetNbinsX(), nVH_err);

      double nggZH_err = -1;
      double nggZH = hggZH->IntegralAndError(1,hggZH->GetNbinsX(), nggZH_err);


      double nboostedZH_err = -1;
      double nboostedZH = hboostedZH->IntegralAndError(1,hboostedZH->GetNbinsX(), nboostedZH_err);


        //double nB = ((TH1*) hs->GetStack()->At(hs->GetHists()->GetSize()-1))->Integral();  //< integral of all added histograms in the stack
      double nB_err =-1; 
      double nB = hmc_exp->IntegralAndError(1, hmc_exp->GetNbinsX(), nB_err);
        
      double nZjHF_err =-1; 
      double nZjHF = hZjHF->IntegralAndError(1, hZjHF->GetNbinsX(), nZjHF_err);

      double nZjLF_err =-1; 
      double nZjLF = hZjLF->IntegralAndError(1, hZjLF->GetNbinsX(), nZjLF_err);
      double nWjLF_err =-1; 
      double nWjLF = hWjLF->IntegralAndError(1, hWjLF->GetNbinsX(), nWjLF_err);

      double nWjHF_err =-1;  
      double nWjHF = hWjHF->IntegralAndError( 1, hWjHF->GetNbinsX(), nWjHF_err);
      double nTT_err =-1;
      double nTT = hTT->IntegralAndError(1, hTT->GetNbinsX(), nTT_err);

      double nST_err =-1;
      double nST = hST->IntegralAndError(1, hST->GetNbinsX(), nST_err);


      double nQCD_err =-1;
      double nQCD = hQCD->IntegralAndError(1, hQCD->GetNbinsX(), nQCD_err);


      double nVV_err =-1;
      double nVV = hVV->IntegralAndError(1, hVV->GetNbinsX(), nVV_err);

        double signif_simple = nS / sqrt(nS + nB);  // only for reference
        double signif_punzi = nS / (3.0/2.0 + sqrt(nB) + 0.2*nB);
        std::cout << "MakePlots(): Printing FOM..." << std::endl;
        std::cout << Form("VH                      = %.3f +- %.3f", nVH, nVH_err) << std::endl;
        std::cout << Form("boostedZHS                      = %.3f +- %.3f", nboostedZH, nboostedZH_err) << std::endl;
        std::cout << Form("ggZHS                      = %.3f +- %.3f", nggZH, nggZH_err) << std::endl;
        std::cout << Form("B                      = %.3f +- %.3f", nB, nB_err) << std::endl;
        std::cout << Form("S/sqrt(S+B)            = %.3f ", signif_simple) << std::endl;
        std::cout << Form("S/(sqrt(B)+a/2+ 0.2*B) = %.3f", signif_punzi) << std::endl;
        std::cout << Form("ZjHF                      = %.3f +- %.3f", nZjHF, nZjHF_err) << std::endl;
        std::cout << Form("WjHF                      = %.3f +- %.3f", nWjHF, nWjHF_err) << std::endl;
        std::cout << Form("ZjLF                      = %.3f +- %.3f", nZjLF, nZjLF_err) << std::endl;
        std::cout << Form("WjLF                      = %.3f +- %.3f", nWjLF, nWjLF_err) << std::endl;
        std::cout << Form("TT                     = %.3f +- %.3f", nTT, nTT_err) << std::endl;
        std::cout << Form("ST                     = %.3f +- %.3f", nST, nST_err) << std::endl;
        std::cout << Form("VV                     = %.3f +- %.3f", nVV, nVV_err) << std::endl;
        std::cout << Form("QCD                     = %.3f +- %.3f", nQCD, nQCD_err) << std::endl;
	if (plotData) {
	  double nData = hdata_obs->IntegralAndError(1, hVV->GetNbinsX(), nVV_err);

	  std::cout << Form("data                     = %.3f ", nData) << std::endl;
	}
    }

    // Setup canvas and pads
    TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
    TPad * pad1 = new TPad("pad1", "top pad"   , 0.0, 0.3, 1.0, 1.0);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();
    TPad * pad2 = new TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.3);
    pad2->SetTopMargin(0.0);
    pad2->SetBottomMargin(0.35);
    pad2->Draw();
    pad1->cd();
    pad1->SetLogy(plotLog);

    // Setup histogram styles
    set_style(hVH, "VH");
    set_style(hZH, "VH");
    set_style(hggZH, "ggZH");
    set_style(hboostedZH, "boostedZH");
    set_style(hboostedZH_2, "boostedZH_2");
    set_style(hboostedZH_3, "boostedZH_3");
    set_style(hWH, "VH");
    set_style(hWjLF, "WjLF");
    set_style(hWjHF, "WjHFb");
    set_style(hZjLF, "ZjLF");
    set_style(hZjHF, "ZjHFb");
    set_style(hTT, "TT");
    set_style(hST, "ST");
    set_style(hVV, "VVHF");
    set_style(hQCD, "QCD");
    set_style(hdata_obs, "data_obs");

    // Setup auxiliary histograms (ratios, errors, etc)
    TH1F * staterr = (TH1F *) hmc_exp->Clone("staterr");
    staterr->Sumw2();
    //staterr->SetFillColor(kRed);
    staterr->SetFillColor(kGray+3);
    staterr->SetMarkerSize(0);
    staterr->SetFillStyle(3013);

    TH1F * ratio = (TH1F *) hdata_obs->Clone("ratio");
    ratio->Sumw2();
    ratio->SetMarkerSize(0.8);
    //ratio->SetMarkerSize(0.5);
    ratio->Divide(hdata_obs, hmc_exp, 1., 1., "");

    TH1F * ratiostaterr = (TH1F *) hmc_exp->Clone("ratiostaterr");
    ratiostaterr->Sumw2();
    ratiostaterr->SetStats(0);
    ratiostaterr->SetTitle(title);
    ratiostaterr->GetYaxis()->SetTitle("Data/MC");
    ratiostaterr->SetMaximum(2.2);
    ratiostaterr->SetMinimum(0);
    ratiostaterr->SetMarkerSize(0);
    //ratiostaterr->SetFillColor(kRed);
    ratiostaterr->SetFillColor(kGray+3);
    ratiostaterr->SetFillStyle(3013);
    ratiostaterr->GetXaxis()->SetLabelSize(0.12);
    ratiostaterr->GetXaxis()->SetTitleSize(0.14);
    ratiostaterr->GetXaxis()->SetTitleOffset(1.10);
    ratiostaterr->GetYaxis()->SetLabelSize(0.10);
    ratiostaterr->GetYaxis()->SetTitleSize(0.12);
    ratiostaterr->GetYaxis()->SetTitleOffset(0.6);
    ratiostaterr->GetYaxis()->SetNdivisions(505);

    TLine* ratiounity = new TLine(xlow,1,xup,1);
    ratiounity->SetLineStyle(2);

    for (Int_t i = 0; i < hmc_exp->GetNbinsX()+2; i++) {
      //for (Int_t i = 0; i < hmc_exp->GetNbinsX(); i++) {
        ratiostaterr->SetBinContent(i, 1.0);
        if (hmc_exp->GetBinContent(i) > 1e-6) {  //< not empty
            double binerror = hmc_exp->GetBinError(i) / hmc_exp->GetBinContent(i);
            ratiostaterr->SetBinError(i, binerror);
        } else {
            ratiostaterr->SetBinError(i, 999.);
        }
    }

    TH1F * ratiosysterr = (TH1F *) ratiostaterr->Clone("ratiosysterr");
    ratiosysterr->Sumw2();
    ratiosysterr->SetMarkerSize(0);
    ratiosysterr->SetFillColor(kYellow-4);
    //ratiosysterr->SetFillStyle(3002);
    ratiosysterr->SetFillStyle(1001);

    for (Int_t i = 0; i < hmc_exp->GetNbinsX()+2; i++) {
      //   for (Int_t i = 0; i < hmc_exp->GetNbinsX(); i++) {
        if (hmc_exp->GetBinContent(i) > 1e-6) {  //< not empty
            double binerror2 = (pow(hmc_exp->GetBinError(i), 2) +
				//                                pow(0.08 * hWjLF->GetBinContent(i), 2) +
                                pow(0.20 * hWjHF->GetBinContent(i), 2) +
				// pow(0.08 * hZjLF->GetBinContent(i), 2) +
                                pow(0.20 * hZjHF->GetBinContent(i), 2) +
                                pow(0.10 * hTT->GetBinContent(i), 2) +
                                pow(0.25 * hST->GetBinContent(i), 2) +
                                pow(0.25 * hVV->GetBinContent(i), 2)+
 	                        pow(0.50 * hQCD->GetBinContent(i), 2));

            double binerror = sqrt(binerror2);
            ratiosysterr->SetBinError(i, binerror / hmc_exp->GetBinContent(i));
        }
    }

    // Setup legends
    TLegend * leg1 = new TLegend(0.50, 0.68, 0.72, 0.92);
    set_style(leg1);
    if (plotData) leg1->AddEntry(hdata_obs, "Data", "p");
    if (plotSig)  leg1->AddEntry(hVH, Form("VH(%i)", ev->massH), "f");
    //    if (plotSig)  leg1->AddEntry(hggZH, Form("ggZH(%i)", ev->massH), "f");
    if (plotBoostedZH)  leg1->AddEntry(hboostedZH, "MZ'=600", "l");
    if (plotBoostedZH)  leg1->AddEntry(hboostedZH_2, "MZ'=800", "l");
    if (plotBoostedZH)  leg1->AddEntry(hboostedZH_3, "MZ'=1000", "l");
    leg1->AddEntry(hTT, "top", "f");
    //    leg1->AddEntry(hST, "single top", "f");


    TLegend * leg2 = new TLegend(0.72, 0.68, 0.94, 0.92);
    set_style(leg2);
    leg2->AddEntry(hVV, "VV", "f");
    leg2->AddEntry(hWjHF, "Wj", "f");
    //    leg2->AddEntry(hWjLF, "W + LF", "f");
    leg2->AddEntry(hZjHF, "Zj", "f");
    //  leg2->AddEntry(hZjLF, "Z + LF", "f");
    leg2->AddEntry(hQCD, "QCD", "f");
    leg2->AddEntry(staterr, "MC uncert. (stat)", "f");

    TLegend * ratioleg = new TLegend(0.72, 0.88, 0.94, 0.96);
    set_style(ratioleg);
    ratioleg->SetTextSize(0.07);
    ratioleg->AddEntry(ratiostaterr, "MC uncert. (stat)", "f");
    //   ratioleg->AddEntry(ratiosysterr, "MC uncert. (syst)", "f");

    TLegend * ratioleg2 = new TLegend(0.50, 0.88, 0.72, 0.96);
    set_style(ratioleg2);
    ratioleg2->SetTextSize(0.07);
    //    ratioleg->AddEntry(ratiostaterr, "MC uncert. (stat)", "f");
    ratioleg2->AddEntry(ratiosysterr, "MC uncert. (syst)", "f");


    // Draw MC signal and backgrounds
    std::clog << "MakePlots(): Drawing..." << std::endl;
    pad1->cd();
    if (plotLog) pad1->SetLogy();

    double ymax = TMath::Max(hdata_obs->GetMaximum(), hs->GetMaximum());
    hs->SetMaximum(ymax * 1.7 + (ymax>1 ? sqrt(ymax) : 0.));
    if (plotLog)  hs->SetMaximum(ymax * 200 + (ymax>1 ? sqrt(ymax) : 0.));
    hs->SetMinimum(0.01);
    hs->Draw("hist");
    hs->GetXaxis()->SetLabelSize(0);
    double binwidth = (xup - xlow) / nbinsx;
    TString ytitle = Form("Events / %.3f", binwidth);
    hs->GetYaxis()->SetTitle(ytitle);

    staterr->Draw("e2 same");
    if (plotSig) {
      /*        hVH->SetLineColor(2);
        hVH->SetLineWidth(3);
        hVH->SetFillColor(0);
        hVH->Draw("hist same");
        hggZH->SetLineColor(kOrange-2);
        hggZH->SetLineWidth(3);
        hggZH->SetFillColor(0);
        hggZH->Draw("hist same");
      */


    }
    if (plotBoostedZH){
      hboostedZH->SetLineColor(kOrange+7);
      hboostedZH->SetLineWidth(3);
      hboostedZH->SetFillColor(0);
      hboostedZH->Draw("hist same");


      hboostedZH_2->SetLineColor(kOrange-5);
      hboostedZH_2->SetLineWidth(3);
      hboostedZH_2->SetFillColor(0);
      hboostedZH_2->Draw("hist same");


      hboostedZH_3->SetLineColor(kOrange-2);
      hboostedZH_3->SetLineWidth(3);
      hboostedZH_3->SetFillColor(0);
      hboostedZH_3->Draw("hist same");


    } 
    if (plotData) {
        hdata_obs->Draw("e1 same");
    }

    // Draw legends
    leg1->Draw();
    leg2->Draw();
    TLatex * latex = new TLatex();
    latex->SetNDC();
    latex->SetTextAlign(12);
    latex->SetTextFont(62);
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.19, 0.89, "CMS Preliminary");
    //latex->DrawLatex(0.19, 0.89, "CMS Simulation 2015");
    latex->SetTextSize(0.04);
    latex->DrawLatex(0.19, 0.84, "#sqrt{s} = 13 TeV, L = 2.3 fb^{-1}");
    // NOTE: change this to your channel
    //latex->DrawLatex(0.19, 0.79, "Z(#mu#bar{#mu})H(b#bar{b})");
    //latex->DrawLatex(0.19, 0.79, "Z(e#bar{e})H(b#bar{b})");
    //latex->DrawLatex(0.19, 0.79, "W(#mu#nu)H(b#bar{b})");
    //latex->DrawLatex(0.19, 0.79, "W(e#nu)H(b#bar{b})");
    //   latex->DrawLatex(0.19, 0.79, "Z(#nu#bar{#nu})H(b#bar{b})");
      if (plotBoostedZH)  latex->DrawLatex(0.19, 0.79, "Z' #rightarrow DM + H(b#bar{b})");
      if (plotBoostedZH)  latex->DrawLatex(0.19, 0.74, CRname);

    // Draw ratio
    pad2->cd();
    pad2->SetGridy(0);
    ratiostaterr->Draw("e2");
    ratiosysterr->Draw("e2 same");
    ratiostaterr->Draw("e2 same");
    ratiounity->Draw();
    ratio->Draw("e1 same");
    ratioleg->Draw();
    ratioleg2->Draw();

    // Kolmogorov-Smirnov test and Chi2 test
    //TPaveText * pave = new TPaveText(0.18, 0.86, 0.28, 0.96, "brNDC");
    TPaveText * pave = new TPaveText(0.25, 0.86, 0.35, 0.96, "brNDC");
    if (plotData) {
        pave->SetLineColor(0);
        pave->SetFillColor(0);
        pave->SetShadowColor(0);
        pave->SetBorderSize(1);
        double nchisq = hdata_obs->Chi2Test(hmc_exp, "UWCHI2/NDF");  // MC uncert. (stat)
        double kolprob = hdata_obs->KolmogorovTest(hmc_exp);  // MC uncert. (stat)
        //TText * text = pave->AddText(Form("#chi_{#nu}^{2} = %.3f", nchisq));
               TText * text = pave->AddText(Form("#chi_{#nu}^{2} = %.3f, K_{s} = %.3f", nchisq, kolprob));
        text->SetTextFont(62);
        text->SetTextSize(0.07);
        //text->SetTextSize(0.06);
	// pave->Draw();
    }

    // Print
    std::clog << "MakePlots(): Printing..." << std::endl;
    pad1->cd();
    gPad->RedrawAxis();
    gPad->Modified();
    gPad->Update();
    pad2->cd();
    gPad->RedrawAxis();
    gPad->Modified();
    gPad->Update();
    c1->cd();

    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");

    TFile* outrootfile = TFile::Open(plotdir+plotname+".root", "RECREATE");
    hZH->Write();
    hVH->Write();
    hboostedZH->Write("monoH");
    hboostedZH_2->Write("monoH_2");
    hboostedZH_3->Write("monoH_3");
    hggZH->Write();
    hWH->Write();
    hWjLF->Write();
    hWjHF->Write("Wj");
    hZjLF->Write();
    hZjHF->Write("Zj");
    hTT->Write();
    hST->Write();
    hVV->Write();
    hQCD->Write();
    hmc_exp->Write();
    hdata_obs->Write();
    outrootfile->Close();

    // Clean up
    delete staterr;
    delete ratio;
    delete ratiostaterr;
    delete ratiosysterr;
    delete leg1;
    delete leg2;
    delete ratioleg;
    delete ratioleg2;
    delete latex;
    delete pave;
    delete hs;
    delete pad1;
    delete pad2;
    delete c1;

    delete hZH;
    delete hggZH;
    delete hWH;
    delete hWjLF;
    delete hWjHF;
    delete hZjLF;
    delete hZjHF;
    delete hTT;
    delete hST;
    delete hVV;
    delete hQCD;
    delete hdata_obs;
    delete hVH;
    delete hboostedZH;
    delete hboostedZH_2;
    delete hboostedZH_3;
    delete hmc_exp;

    std::clog << "MakePlots(): DONE!" << std::endl;

    return;
}

//______________________________________________________________________________
void MakeDatacard(TString channel, TString dcname, TString wsname,
                  bool useshapes,
                  TString options="unblind:SplusB") {


    std::clog << "MakeDatacard(): starting " << std::endl;

    // Parse options
    bool unblind    = options.Contains("unblind")    && (!options.Contains("!unblind"));
    bool SplusB     = options.Contains("SplusB")     && (!options.Contains("!SplusB"));
    channel += "_13TeV";
    
    if (useshapes && !unblind) {
        std::clog << "MakeDatacard(): I only support 'unblind' mode when using shapes." << std::endl;
        unblind = true;
    }
    

    std::ofstream dc;
    dc.setf(ios::fixed,ios::floatfield);
    dc.precision(3);
    dc.open(dcname.Data());



    TString wsname_new = dcname(0,dcname.Sizeof()-5)+".root";
    gSystem->Exec(Form("cp %s %s", wsname.Data(), wsname_new.Data()));
    TFile * wsfile = TFile::Open(wsname);


    std::clog << "MakeDatacard(): getting the histos " << std::endl;

      double nVH_err = -1;
      double nVH =  ((TH1 *) wsfile->Get("VH"))->IntegralAndError(1,  999 , nVH_err);




      double nZH_err = -1;
      double nZH =  ((TH1 *) wsfile->Get("ZH"))->IntegralAndError(1,  999 , nZH_err);


    std::clog << "MakeDatacard(): got ZH " << std::endl;

      double nWH_err = 0;
      double nWH =  ((TH1 *) wsfile->Get("WH"))->IntegralAndError(1,  999 , nWH_err);

    std::clog << "MakeDatacard(): got WH " << std::endl;

      double nggZH_err = -1;
      double nggZH =  ((TH1 *) wsfile->Get("ggZH"))->IntegralAndError(1,  999 , nggZH_err);

    std::clog << "MakeDatacard(): got ggZH " << std::endl;

      double nboostedZH_err = -1;
      double nboostedZH =  ((TH1 *) wsfile->Get("monoH"))->IntegralAndError(1,  999 , nboostedZH_err);

    std::clog << "MakeDatacard(): got boostedZH " << std::endl;

      double nZjHF_err =-1; 
      double nZjHF =  ((TH1 *) wsfile->Get("Zj"))->IntegralAndError(1,  999 , nZjHF_err);

      //      double nZjLF_err =-1; 
      // double nZjLF =  ((TH1 *) wsfile->Get("ZjLF"))->IntegralAndError(1,  999 , nZjLF_err);


      double nWjHF_err =-1; 
      double nWjHF =  ((TH1 *) wsfile->Get("Wj"))->IntegralAndError(1,  999 , nWjHF_err);

      // double nWjLF_err =-1; 
      //double nWjLF =  ((TH1 *) wsfile->Get("WjLF"))->IntegralAndError(1,  999 , nWjLF_err);

      double nTT_err =-1; 
      double nTT =  ((TH1 *) wsfile->Get("TT"))->IntegralAndError(1,  999 , nTT_err);



      double nST_err =-1; 
      double nST =  ((TH1 *) wsfile->Get("ST"))->IntegralAndError(1,  999 , nST_err);

      double nVV_err =-1; 
      double nVV =  ((TH1 *) wsfile->Get("VV"))->IntegralAndError(1,  999 , nVV_err);


      double nQCD_err =-1; 
      double nQCD=  ((TH1 *) wsfile->Get("QCD"))->IntegralAndError(1,  999 , nQCD_err);



    int jmax = 6;
    dc << "imax 1 number of channels" << std::endl;
    dc << "jmax " << jmax << " number of processes minus 1 ('*' = automatic)" << std::endl;
    dc << "kmax * number of nuisance parameters (sources of systematical uncertainties)" << std::endl;
    dc << "-----------------------------------" << std::endl;
    if (useshapes) {
        //dc << "shapes * * " << wsname << " $CHANNEL:$PROCESS $CHANNEL:$PROCESS_$SYSTEMATIC" << std::endl;
        dc << "shapes * * " << wsname_new << " $PROCESS $PROCESS_$SYSTEMATIC" << std::endl;
        dc << "-----------------------------------" << std::endl;
    }
    dc << "bin         " << channel << std::endl;
    if (unblind) {
        dc << "observation " << ((TH1 *) wsfile->Get("data_obs"))->Integral() << std::endl;
    } else if (SplusB) {
      dc << "observation " << ((TH1 *) wsfile->Get("mc_exp"))->Integral() + ((TH1 *) wsfile->Get("VH"))->Integral()  + ((TH1 *) wsfile->Get("monoH"))->Integral()  << std::endl;
    } else {
        dc << "observation " << ((TH1 *) wsfile->Get("mc_exp"))->Integral() << std::endl;
    }
    dc << "#prediction " << ((TH1 *) wsfile->Get("mc_exp"))->Integral() << std::endl;
    dc << "-----------------------------------" << std::endl;

    dc << "### all sognal rates        " << setw(10) << right << ((TH1 *) wsfile->Get("monoH"))->Integral() << " " 
                                          << setw(10) << right << ((TH1 *) wsfile->Get("monoH_2"))->Integral() << " " 
       << setw(10) << right << ((TH1 *) wsfile->Get("monoH_3"))->Integral() << std::endl;

    dc << "bin         "; for (int j=0; j!=jmax+1; j++)  dc << channel << "   "; dc << std::endl;
    dc << "process     monoH       VH       Wj         Zj          TT         VV       QCD  " << std::endl;
    dc << "process     -1          1        2          3            4          5        6   " << std::endl;
    dc << "rate        " << setw(10) << right << ((TH1 *) wsfile->Get("monoH"))->Integral() << " "
      //                         << setw(10) << right << ((TH1 *) wsfile->Get("ggZH"))->Integral() << " "
      //                   << setw(10) << right << ((TH1 *) wsfile->Get("ZH"))->Integral() << " "
                   << setw(10) << right << ((TH1 *) wsfile->Get("VH"))->Integral() << " "
      //                   << setw(10) << right << ((TH1 *) wsfile->Get("WH"))->Integral() << " "
      //                         << setw(10) << right << ((TH1 *) wsfile->Get("WjLF"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("Wj"))->Integral() << " "
      //                 << setw(10) << right << ((TH1 *) wsfile->Get("ZjLF"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("Zj"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("TT"))->Integral() << " "
      // << setw(10) << right << ((TH1 *) wsfile->Get("ST"))->Integral() << " "
                         << setw(10) << right << ((TH1 *) wsfile->Get("VV"))->Integral() << "  " 
                         << setw(10) << right << ((TH1 *) wsfile->Get("QCD"))->Integral() <<std::endl;
    dc.precision(2);
    dc << "-----------------------------------" << std::endl;
    dc << "" << std::endl;
    dc << "# This datacard is for demonstration purposes only!" << std::endl;
    dc << "#################################### #####  monoH   VH     Wj    Zj     TT    VV    QCD " << std::endl;
    dc << "lumi_13TeV                           lnN    1.026   1.026  -     -      -     1.026 1.026" << std::endl;
    dc << "pdf_qqbar                            lnN    1.01    1.01   -     -      -     -     -" << std::endl;
    dc << "pdf_gg                               lnN    -        -     -     -      -     -     1.50  "  << std::endl;
    dc << "QCDscale_VH                          lnN    -        1.04  -     -      -     -     - " << std::endl;
    //    dc << "QCDscale_ggZH                        lnN    -      1.18   -  -     -     -        -     -     -     -     " << std::endl;
    //    dc << "QCDscale_ttbar                       lnN    -      -      -     -     -         -     1.06  -     -     -     " << std::endl;
    dc << "QCDscale_VV                          lnN     -      -      -     -      -    1.04  -     " << std::endl;
    //    dc << "CMS_monoHbb_boost_EWK_13TeV             lnN    1.05      -      1.05  1.10  -         -     -     -     -     -     " << std::endl;
    //  dc << "CMS_monoHbb_boost_QCD_13TeV             lnN    1.10      1.10   1.10  1.10  -      -     -     -     -     -     " << std::endl;
    //    dc << "CMS_monoHbb_ST                          lnN    -      -      -      -     -        -     -    1.25  -     -     " << std::endl;
    dc << "CMS_monoHbb_VV                       lnN     -       -     -    -      -      1.25  -     " << std::endl;
    dc << "CMS_monoHbb_eff_b                    lnN    1.07   1.07    -     -     -   1.07  1.07  " << std::endl;
    dc << "CMS_monoHbb_fake_b_13TeV             lnN     -       -     -     -        -       1.03     1.03" << std::endl;
    dc << "CMS_monoHbb_res_j                    lnN    1.03   1.03    -     -     -    1.03  1.03  " << std::endl;
    dc << "CMS_monoHbb_scale_j                  lnN    1.05   1.05    -     -     -    1.05  1.05  1.05  " << std::endl;
    dc << "CMS_monoHbb_WjHFUnc_SF               lnN    -      -      1.10   -     -     -      -     " << std::endl;
    //   dc << "CMS_monoHbb_Zj_SF                     lnN    -      -     -     -     -      1.10  -     -     -     -      -     " << std::endl;
    //    dc << "CMS_monoHbb_ZjHF_SF                     lnN    -      -     -     -     -     -     -     1.50  -     -     -      -     " << std::endl;
    // dc << "CMS_monoHbb_TT_SF                       lnN    -      -     -         -     -     -     1.10  -     -      -     " << std::endl;

    dc << "CMS_monoHbb_monoH_stat                lnN    "<< 1.+nboostedZH_err/nboostedZH<<"  -        -     -     -     -     -     " << std::endl;
    //    dc << "CMS_monoHbb_ggZH_stat                   lnN    -     "<< 1.+nggZH_err/nggZH<<"  -     -        -     -     -     -     -     -     " << std::endl;
    //  dc << "CMS_monoHbb_ZH_stat                     lnN    -     -     "<< 1.+nZH_err/nZH<<"  -     -         -     -     -     -     -     " << std::endl;
    // dc << "CMS_monoHbb_WH_stat                     lnN    -     -     -     "<< 1.+nWH_err/nWH<<"  -        -     -     -     -     -     " << std::endl;
    dc << "CMS_monoHbb_VH_stat                   lnN    -      "<< 1.+nVH_err/nWH<<"  -         -     -     -     -     " << std::endl;
    dc << "CMS_monoHbb_Wj_stat                   lnN    -     -     "<< 1.+nWjHF_err/nWjHF<<"  -     -     -     -     -     -         " << std::endl;
    //    dc << "CMS_monoHbb_WjHF_stat                   lnN    -     -     -     -     -     "<<1.+nWjHF_err/nWjHF<<"  -     -     -     -     -     -              " << std::endl;
    // dc << "CMS_monoHbb_ZjLF_stat                   lnN    -     -     -     -     -     -     "<<1.+nZjLF_err/nZjLF<<"  -     -     -     -     -                   " << std::endl;
    dc << "CMS_monoHbb_Zj_stat                   lnN    -     -       -     "<<1.+nZjHF_err/nZjHF<<"  -       -     -                         " << std::endl;
    dc << "CMS_monoHbb_TT_stat                     lnN    -     -     -     -     "<<1.+nTT_err/nTT<<"      -     -                               " << std::endl;
    //    dc << "CMS_monoHbb_ST_stat                     lnN    -     -     -     -     "<<1.+nST_err/nST<<"  -     -                               " << std::endl;
    dc << "CMS_monoHbb_VV_stat                     lnN     -     -     -     -     -    "<<1.+nVV_err/nVV<<"  -                               " << std::endl;
    dc << "CMS_monoHbb_QCD_stat                    lnN     -     -     -     -     -    -    "<<1.+nQCD_err/nQCD<<"                                " << std::endl;


    dc << "CMS_monoHbb_trigger_MET                 lnN    1.03   1.03   -        -     -     1.03  1.03 " << std::endl;
    dc << "#################################### #####  monoH    VH    Wj  Zj  TT    VV    QCD  " << std::endl;
    dc.close();

    wsfile->Close();

    std::clog << "MakeDatacard(): DONE!" << std::endl;

    return;
}




//______________________________________________________________________________
Events::Events()
  : ZH(0),
    ggZH(0),
    WH(0),
    WjLF(0),
    WjHF(0),
    ZjLF(0),
    ZjHF(0),
    TT(0),
    ST_T_s(0),
    ST_T_t(0),
    ST_T_tW(0),
    ST_Tbar_s(0),
    ST_Tbar_t(0),
    ST_Tbar_tW(0),
    VV(0),
    //    VV_WZ(0),
    // VV_ZZ(0),
    data_obs(0),
  sf_WjLF(1.32), sf_WjHF(1.32), sf_ZjLF(1.), sf_ZjHF(1.), sf_TT(0.87),  
    massH(125) {}

Events::~Events() {
    delete ZH;
    delete ggZH;
    delete WH;
    delete WjLF;
    delete WjHF;
    delete ZjLF;
    delete ZjHF;
    delete TT;
    delete ST_T_s;
    delete ST_T_t;
    delete ST_T_tW;
    delete ST_Tbar_s;
    delete ST_Tbar_t;
    delete ST_Tbar_tW;
    delete VV;
    //    delete VV_WW;
    //  delete VV_WZ;
    //delete VV_ZZ;
    delete data_obs;
}

void Events::read( TCut cutmc_all, TCut cutdata_all, TString processes , TString treename ) {
    //TString indir = "/uscms_data/d2/jiafu/CMSDAS2014/VHbbAnalysis/skim/";
    //TString indir = "/eos/uscms/store/user/cmsdas/2014/Hbb/Step2/";
  //    TString indir = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Phys14_PU20bx25/skimV11/step3/";
  // TString indir = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Spring15_PU20bx25/skimV12_v2/step3/skim_ZnnH_classification/step4/"; //Step4_ZnnH125.root
TString indir = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/V21_FromSean/step3/"; //Step4_ZnnH125.root


    TString prefix = "Step3_";
    TString suffix = ".root";
    //     TString treename = "tree_ZnunuHighPt_test";
     //TString treename = "tree";

    TCut Vbb  = "Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numBHadrons>0)>=2";
    TCut Vb  = "Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numBHadrons>0)==1";
    TCut Vcc  = TCut("Sum$(GenJet_pt>20 && abs(GenJet_eta)<2.4 && GenJet_numCHadrons>0)>=1") && !Vb && !Vbb;
    TCut Vudsg = !(Vbb || Vb || Vcc);


    TCut cutHF = Vb || Vbb ;//"(abs(Jet_hadFlavour[hJCidx[0]])==5 || abs(Jet_mcFlavour[hJCidx[1]])==5)";  //< for b quarks, pdgId = 5 
    TCut cutLF = !(Vbb || Vb) ;//"(abs(Jet_mcFlavour[hJCidx[0]])!=5 && abs(Jet_mcFlavour[hJCidx[1]])!=5)";  //< for non-b quarks

    //    std::map<std::string, float> lumis = GetLumis();

    // parse processes
    bool loadAll  = processes == "" || processes.Contains("ALL") || processes.Contains("All") || processes.Contains("all");

    bool loadVH   = processes.Contains("VH");
    bool loadggZH   = processes.Contains("VH") || loadAll;
    bool loadboostedZH   = processes.Contains("boostedZH") || loadAll ;
    bool loadZH   = processes.Contains("ZH") || loadVH || loadAll;
    bool loadWH   = processes.Contains("WH") || loadVH || loadAll;
    bool loadWJ   = processes.Contains("WJ") || processes.Contains("WjLF") || processes.Contains("WjHF") || loadAll;
    bool loadZJ   = processes.Contains("ZJ") ||   processes.Contains("ZjLF") ||  processes.Contains("ZjHF") ||   loadAll;
    bool loadTT   = processes.Contains("TT") || loadAll;
    bool loadST   = processes.Contains("ST") || loadAll;
    bool loadVV   = processes.Contains("VV") || loadAll;
    bool loadQCD   = processes.Contains("QCD") || loadAll;
    bool loadData = processes.Contains("DATA") || processes.Contains("Data") || processes.Contains("data") || loadAll;


    // Monte Carlo______________________________________________________________
    // NOTE: for Zll, use the default "ZllH"
    // NOTE: for Znn, change "ZllH" to "ZnnH"
    if (loadZH) {
        TChain ZH_(treename);
        //ZH_.Add(indir + prefix + Form("ZllH%i", massH) + suffix);
        //ZH = (TTree *) ZH_.CopyTree(cutmc_all);
        //lumi_ZH = lumis[Form("ZllH%i", massH)];
        ZH_.Add(indir + prefix + Form("ZnnH%i", massH) + suffix);
        ZH = (TTree *) ZH_.CopyTree(cutmc_all);
	//        lumi_ZH = lumis[Form("ZnnH%i", massH)];
        std::clog << "... DONE: ZH copy tree. N=" << ZH->GetEntries() << std::endl;
    }

    if (loadboostedZH) {
        TChain boostedZH_(treename);
        boostedZH_.Add(indir + prefix + "monoH_600" + suffix);
        boostedZH = (TTree *) boostedZH_.CopyTree(cutmc_all);
        std::clog << "... DONE: boostedZH copy tree. N=" << boostedZH->GetEntries() << std::endl;

        TChain boostedZH_2_(treename);
        boostedZH_2_.Add(indir + prefix + "monoH_800" + suffix);
        boostedZH_2 = (TTree *) boostedZH_2_.CopyTree(cutmc_all);
        std::clog << "... DONE: boostedZH_2 copy tree. N=" << boostedZH_2->GetEntries() << std::endl;


        TChain boostedZH_3_(treename);
        boostedZH_3_.Add(indir + prefix + "monoH_1000" + suffix);
        boostedZH_3 = (TTree *) boostedZH_3_.CopyTree(cutmc_all);
        std::clog << "... DONE: boostedZH_3 copy tree. N=" << boostedZH_3->GetEntries() << std::endl;



    }


    if (loadggZH) {
        TChain ggZH_(treename);
        ggZH_.Add(indir + prefix + Form("ggZH%i", massH) + suffix);
        ggZH = (TTree *) ggZH_.CopyTree(cutmc_all);
	//        lumi_ggZH = lumis[Form("ggZH%i", massH)];
        std::clog << "... DONE: ggZH copy tree. N=" << ggZH->GetEntries() << std::endl;
    }


    if (loadWH) {
       TChain WH_(treename);
        WH_.Add(indir + prefix + Form("WlnH%i", massH) + suffix);
       WH = (TTree *) WH_.CopyTree(cutmc_all);
	//        lumi_WH = lumis[Form("WlnH%i", massH)];
              std::clog << "... DONE: WH copy tree. N=" << WH->GetEntries() << std::endl;
    }

    if (loadWJ) {
        TChain WjLF_(treename);
	//        WjLF_.Add(indir + prefix + "WJetsPtW100" + suffix);
        WjLF_.Add(indir + prefix + "WJets" + suffix);
        WjLF = (TTree *) WjLF_.CopyTree( cutmc_all + cutLF  );
	//        lumi_WjLF = lumis["WJetsIncl"];
        std::clog << "... DONE: WjLF copy tree. N=" << WjLF->GetEntries() << std::endl;

        TChain WjHF_(treename);
        WjHF_.Add(indir + prefix + "WJets" + suffix);
        WjHF = (TTree *) WjHF_.CopyTree( cutmc_all + cutHF);
	//        lumi_WjHF = lumis["WJetsIncl"];
        std::clog << "... DONE: WjHF copy tree. N=" << WjHF->GetEntries() << std::endl;
    }

    if (loadZJ) {
 
        TChain ZjLF_(treename);
	//        ZjLF_.Add(indir + prefix + "WJetsPtW100" + suffix);
        ZjLF_.Add(indir + prefix + "ZJets" + suffix);
        ZjLF = (TTree *) ZjLF_.CopyTree( cutmc_all + cutLF  );
	//        lumi_ZjLF = lumis["ZJetsHT800"];
        std::clog << "... DONE: ZjLF copy tree. N=" << ZjLF->GetEntries() << std::endl;

        TChain ZjHF_(treename);
        ZjHF_.Add(indir + prefix + "ZJets" + suffix);
        ZjHF = (TTree *) ZjHF_.CopyTree( cutmc_all + cutHF);
	//        lumi_ZjHF = lumis["ZJetsHT800"];
        std::clog << "... DONE: ZjHF copy tree. N=" << ZjHF->GetEntries() << std::endl;

      /*

        TChain ZjLF_HT100_(treename);
        ZjLF_HT100_.Add(indir + prefix + "ZJetsHT100" + suffix);
        ZjLF_HT100 = (TTree *) ZjLF_HT100_.CopyTree(cutmc_all + cutLF);
        lumi_ZjLFHT100 = lumis["ZJetsHT100"];
        std::clog << "... DONE: ZjLF HT100 copy tree. N=" << ZjLF_HT100->GetEntries() << std::endl;


        TChain ZjLF_HT200_(treename);
        ZjLF_HT200_.Add(indir + prefix + "ZJetsHT200" + suffix);
        ZjLF_HT200 = (TTree *) ZjLF_HT200_.CopyTree(cutmc_all + cutLF);
        lumi_ZjLFHT200 = lumis["ZJetsHT200"];
        std::clog << "... DONE: ZjLF HT200 copy tree. N=" << ZjLF_HT200->GetEntries() << std::endl;

	TChain ZjLF_HT800_(treename);
        ZjLF_HT800_.Add(indir + prefix + "ZJetsHT800" + suffix);
        ZjLF_HT800 = (TTree *) ZjLF_HT800_.CopyTree(cutmc_all + cutLF);
        lumi_ZjLFHT800 = lumis["ZJetsHT800"];
        std::clog << "... DONE: ZjLF HT800 copy tree. N=" << ZjLF_HT800->GetEntries() << std::endl;



        TChain ZjHF_HT100_(treename);
        ZjHF_HT100_.Add(indir + prefix + "ZJetsHT100" + suffix);
        ZjHF_HT100 = (TTree *) ZjHF_HT100_.CopyTree(cutmc_all + cutHF);
        lumi_ZjHFHT100 = lumis["ZJetsHT100"];
        std::clog << "... DONE: ZjHF HT100 copy tree. N=" << ZjHF_HT100->GetEntries() << std::endl;


        TChain ZjHF_HT200_(treename);
        ZjHF_HT200_.Add(indir + prefix + "ZJetsHT200" + suffix);
        ZjHF_HT200 = (TTree *) ZjHF_HT200_.CopyTree(cutmc_all + cutHF);
        lumi_ZjHFHT200 = lumis["ZJetsHT200"];
        std::clog << "... DONE: ZjHF HT200 copy tree. N=" << ZjHF_HT200->GetEntries() << std::endl;

	TChain ZjHF_HT800_(treename);
        ZjHF_HT800_.Add(indir + prefix + "ZJetsHT800" + suffix);
        ZjHF_HT800 = (TTree *) ZjHF_HT800_.CopyTree(cutmc_all + cutHF);
        lumi_ZjHFHT800 = lumis["ZJetsHT800"];
        std::clog << "... DONE: ZjHF HT800 copy tree. N=" << ZjHF_HT800->GetEntries() << std::endl;

      */



    }

    /* if (loadZjLF) { */
    /*     // NOTE: for Zll, use the default "DYJetsPtZ100" */
    /*     // NOTE: for Znn, change "DYJetsPtZ100" to "ZJetsPtZ100" */
    /*     TChain ZjLF_(treename); */
    /*     //ZjLF_.Add(indir + prefix + "DYJetsPtZ100" + suffix); */
    /*     //ZjLF = (TTree *) ZjLF_.CopyTree(cutmc_all + cutLF); */
    /*     //lumi_ZjLF = lumis["DYJetsPtZ100"]; */
    /*     ZjLF_.Add(indir + prefix + "ZJets" + suffix); */
    /*     ZjLF = (TTree *) ZjLF_.CopyTree(cutmc_all + cutLF); */
    /*     lumi_ZjLF = lumis["ZJetsHT800"]; */
    /*     std::clog << "... DONE: ZjLF copy tree. N=" << ZjLF->GetEntries() << std::endl; */

    /* } */


    /* if (loadZjHF) { */
    /*     // NOTE: for Zll, use the default "DYJetsPtZ100" */
    /*     // NOTE: for Znn, change "DYJetsPtZ100" to "ZJetsPtZ100" */
    /*     TChain ZjHF_(treename); */
    /*     //ZjLF_.Add(indir + prefix + "DYJetsPtZ100" + suffix); */
    /*     //ZjLF = (TTree *) ZjLF_.CopyTree(cutmc_all + cutLF); */
    /*     //lumi_ZjLF = lumis["DYJetsPtZ100"]; */
    /*     ZjHF_.Add(indir + prefix + "ZJets" + suffix); */
    /*     ZjHF = (TTree *) ZjHF_.CopyTree(cutmc_all + cutHF); */
    /*     lumi_ZjHF = lumis["ZJetsHT800"]; */
    /*     std::clog << "... DONE: ZjHF copy tree. N=" << ZjHF->GetEntries() << std::endl; */

    /* } */




    if (loadTT) {
        TChain TT_(treename);
	TT_.Add(indir + prefix + "TTMad" + suffix);
	//	    	TT_.Add(indir + prefix + "TTMadv2" + suffix);
        TT = (TTree *) TT_.CopyTree(cutmc_all);
	//       	lumi_TT = lumis["TTPythia8"];
	// lumi_TT = lumis["TTMadv2"];
        std::clog << "... DONE: TT copy tree. N=" << TT->GetEntries() << std::endl;

    }

    if (loadST) {
      /*
        TChain ST_T_s_(treename);
        ST_T_s_.Add(indir + prefix + "T_s" + suffix);
        ST_T_s = (TTree *) ST_T_s_.CopyTree(cutmc_all);
        lumi_ST_T_s = lumis["T_s"];
        std::clog << "... DONE: ST_T_s copy tree. N=" << ST_T_s->GetEntries() << std::endl;

        TChain ST_T_t_(treename);
        ST_T_t_.Add(indir + prefix + "T_t" + suffix);
        ST_T_t = (TTree *) ST_T_t_.CopyTree(cutmc_all);
        lumi_ST_T_t = lumis["T_t"];
        std::clog << "... DONE: ST_T_t copy tree. N="  << ST_T_t->GetEntries() << std::endl;

        TChain ST_T_tW_(treename);
        ST_T_tW_.Add(indir + prefix + "T_tW" + suffix);
        ST_T_tW = (TTree *) ST_T_tW_.CopyTree(cutmc_all);
        lumi_ST_T_tW = lumis["T_tW"];
        std::clog << "... DONE: ST_T_tW copy tree. N=" << ST_T_tW->GetEntries() << std::endl;

        TChain ST_Tbar_s_(treename);
        ST_Tbar_s_.Add(indir + prefix + "Tbar_s" + suffix);
        ST_Tbar_s = (TTree *) ST_Tbar_s_.CopyTree(cutmc_all);
        lumi_ST_Tbar_s = lumis["Tbar_s"];
        std::clog << "... DONE: ST_Tbar_s copy tree. N=" << ST_Tbar_s->GetEntries() << std::endl;

        TChain ST_Tbar_t_(treename);
        ST_Tbar_t_.Add(indir + prefix + "Tbar_t" + suffix);
        ST_Tbar_t = (TTree *) ST_Tbar_t_.CopyTree(cutmc_all);
        lumi_ST_Tbar_t = lumis["Tbar_t"];
        std::clog << "... DONE: ST_Tbar_t copy tree. N=" << ST_Tbar_t->GetEntries() << std::endl;

        TChain ST_Tbar_tW_(treename);
        ST_Tbar_tW_.Add(indir + prefix + "Tbar_tW" + suffix);
        ST_Tbar_tW = (TTree *) ST_Tbar_tW_.CopyTree(cutmc_all);
        lumi_ST_Tbar_tW = lumis["Tbar_tW"];
        std::clog << "... DONE: ST_Tbar_tW copy tree. N=" << ST_Tbar_tW->GetEntries() << std::endl; 
      */
                                                                                                                                                                                               
      TChain s_Top_(treename);                                                                                                                                                                
      s_Top_.Add(indir + prefix + "s_Top" + suffix);                                                                                                                                        
      s_Top = (TTree *) s_Top_.CopyTree(cutmc_all);                                                                                                                                      
      //      lumi_s_Top = lumis["Tbar_tW"];                                                                                                                                                          
      std::clog << "... DONE: s_Top copy tree. N=" << s_Top->GetEntries() << std::endl; 

    }

    if (loadVV) {


        TChain VV_(treename);
	VV_.Add(indir + prefix + "VV" + suffix);
        VV = (TTree *) VV_.CopyTree(  cutmc_all );
	//  lumi_VV_ = lumis["WW"];
        std::clog << "... DONE: VV copy tree. N=" << VV ->GetEntries() << std::endl;



      /*
        TChain VV_WW_(treename);
        VV_WW_.Add(indir + prefix + "WW" + suffix);
        VV_WW = (TTree *) VV_WW_.CopyTree(cutmc_all);
        lumi_VV_WW = lumis["WW"];
        std::clog << "... DONE: VV_WW copy tree. N=" << VV_WW->GetEntries() << std::endl;

        TChain VV_WZ_(treename);
        VV_WZ_.Add(indir + prefix + "WZ" + suffix);
        VV_WZ = (TTree *) VV_WZ_.CopyTree(cutmc_all);
        lumi_VV_WZ = lumis["WZ"];
        std::clog << "... DONE: VV_WZ copy tree. N=" << VV_WZ->GetEntries() << std::endl;

        TChain VV_ZZ_(treename);
        VV_ZZ_.Add(indir + prefix + "ZZ" + suffix);
        VV_ZZ = (TTree *) VV_ZZ_.CopyTree(cutmc_all);
        lumi_VV_ZZ = lumis["ZZ"];
        std::clog << "... DONE: VV_ZZ copy tree. N=" << VV_ZZ->GetEntries() << std::endl;
      */

    }



    if (loadQCD) {
        TChain QCD_(treename);
	QCD_.Add(indir + prefix + "QCD" + suffix);
        QCD = (TTree *) QCD_.CopyTree(cutmc_all);
	//        lumi_QCD = lumis["QCD_HT800"];
        std::clog << "... DONE: QCD copy tree. N=" << QCD->GetEntries() << std::endl;
    }





    // Data_____________________________________________________________________
    // NOTE: for Zmm and Wmn, use the default "SingleMu"
    // NOTE: for Wen, change both "SingleMu" to "SingleEl"
    // NOTE: for Zee, change both "SingleMu" to "DoubleEl"
    // NOTE: for Znn, change both "SingleMu" to "MET"
    if (loadData) {
        TChain data_obs_(treename);
	//       	data_obs_.Add(indir + prefix + "MET" + suffix);
	data_obs_.Add(indir + "monoHbbHighPt_V21_MET.root");
	//	data_obs_.Add("/afs/cern.ch/user/d/degrutto/eos/cms/store/group/phys_higgs/hbb/ntuples/V21/user/arizzi/VHBBHeppyV21/MET/VHBB_HEPPY_V21_MET__Run2015D-16Dec2015-v1/160317_131113/0000/tree_*.root");
	//        data_obs_.Add(indir +  "skim_Step3_MET__Run2015C-05Oct2015" + suffix );
        //data_obs_.Add(indir +  "skim_Step3_MET__Run2015D-05Oct2015" + suffix );
        //data_obs_.Add(indir +  "skim_Step3_MET__Run2015D-PromptReco" + suffix );



        //data_obs_.Add(indir + prefix + "SingleMu_Prompt" + suffix);
        //data_obs_.Add(indir + prefix + "SingleEl_ReReco" + suffix);
        //data_obs_.Add(indir + prefix + "SingleEl_Prompt" + suffix);
        //data_obs_.Add(indir + prefix + "DoubleEl_ReReco" + suffix);
        //data_obs_.Add(indir + prefix + "DoubleEl_Prompt" + suffix);

	/*        data_obs_.Add(indir + prefix + "MET_ReReco" + suffix);
        data_obs_.Add(indir + prefix + "MET_Prompt" + suffix);
      	*/

	data_obs = (TTree *) data_obs_.CopyTree(cutdata_all);
        std::clog << "... DONE: data_obs copy tree. N=" << data_obs->GetEntries() << std::endl;

    }

    return;
}


void set_style(TH1 * h, const TString& p) {
    if (p == "VH") {
        h->SetFillColor  (2);
        h->SetMarkerColor(2);
    } else if (p == "data_obs" || p == "Data") {
        h->SetMarkerSize(0.8);
        h->SetMarkerStyle(20);
    } else {
        h->SetLineColor(kBlack);
        if (p == "WjLF") {
            h->SetFillColor  (814);
            h->SetMarkerColor(814);
        } else if (p == "WjHFc") {
            h->SetFillColor  (816);
            h->SetMarkerColor(816);
        } else if (p == "WjHFb") {
            h->SetFillColor  (820);
            h->SetMarkerColor(820);
        } else if (p == "ZjLF") {
            h->SetFillColor  (401);
            h->SetMarkerColor(401);
        } else if (p == "ZjHFc") {
            h->SetFillColor  (41);
            h->SetMarkerColor(41);
        } else if (p == "ZjHFb") {
            h->SetFillColor  (5);
            h->SetMarkerColor(5);
        } else if (p == "TT") {
            h->SetFillColor  (596);
            h->SetMarkerColor(596);
        } else if (p == "ST") {
            h->SetFillColor  (840);
            h->SetMarkerColor(840);
        } else if (p == "VV") {
            h->SetFillColor  (922);
            h->SetMarkerColor(922);
        } else if (p == "VVHF") {
            h->SetFillColor  (920);
            h->SetMarkerColor(920);
        } else if (p == "QCD") {
            h->SetFillColor(616);
            h->SetMarkerColor(616);
        } else if (p == "boostedZH") {
            h->SetFillColor  (kOrange+7);
            h->SetMarkerColor(kOrange+7);
        } else if (p == "boostedZH_2") {
            h->SetFillColor  (kOrange-5);
            h->SetMarkerColor(kOrange-5);
        } else if (p == "boostedZH_3") {
            h->SetFillColor  (kOrange-2);
            h->SetMarkerColor(kOrange-2);
        }
	else if (p == "ggZH") {
            h->SetFillColor  (kOrange-2);
            h->SetMarkerColor(kOrange-2);
        }
    }
    return;
}


void set_style_norm(TH1 * h, const TString& p) {
    if (p == "VH") {
        h->SetMaximum(10.0); 
        h->SetLineColor(2);
        h->SetLineWidth(3);
        h->SetMarkerColor(2);
    } else if (p == "data_obs" || p == "Data") {
        h->SetMaximum(10.0); 
        h->SetMarkerSize(0.8);
        h->SetLineWidth(3);
        h->SetMarkerStyle(20);
    } else {
        h->SetLineColor(kBlack);
        if (p == "WjLF") {
        h->SetMaximum(10.0); 
	  h->SetLineColor  (814);
        h->SetLineWidth(3);
            h->SetMarkerColor(814);
        } else if (p == "WjHFc") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (816);
            h->SetMarkerColor(816);
        } else if (p == "WjHFb") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (820);
            h->SetMarkerColor(820);
        } else if (p == "ZjLF") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (401);
            h->SetMarkerColor(401);
        } else if (p == "ZjHFc") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (41);
            h->SetMarkerColor(41);
        } else if (p == "ZjHFb") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (5);
            h->SetMarkerColor(5);
        } else if (p == "TT") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (596);
            h->SetMarkerColor(596);
        } else if (p == "ST") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (840);
            h->SetMarkerColor(840);
        } else if (p == "VV") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (922);
            h->SetMarkerColor(922);
        } else if (p == "VVHF") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (920);
            h->SetMarkerColor(920);
        } else if (p == "QCD") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
          h->SetLineColor  (616);
          h->SetMarkerColor(616);
        } else if (p == "boostedZH") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor(kOrange+7);
            h->SetMarkerColor(kOrange+7);
        } else if (p == "boostedZH_2") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor(kOrange-5);
            h->SetMarkerColor(kOrange-5);
        } else if (p == "boostedZH_3") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor(kOrange-2);
            h->SetMarkerColor(kOrange-2);
        } else if (p == "ggZH") {
        h->SetMaximum(10.0);
        h->SetLineWidth(3); 
            h->SetLineColor  (kOrange-2);
            h->SetMarkerColor(kOrange-2);
        }

    }
    return;
}


void set_style(TLegend* leg) {
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->SetShadowColor(0);
    leg->SetTextFont(62);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(1);
    return;
}





