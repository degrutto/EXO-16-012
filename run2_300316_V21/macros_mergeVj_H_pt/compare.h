#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <vector>

#include "TSystem.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TTree.h"
#include "TH1.h"
#include "THStack.h"
#include "TCut.h"
#include "TString.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif




// _____________________________________________________________________________
// These are style functions
void set_legstyle(TLegend* leg) {
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextFont(62);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(1);
  return;
}




void set_style(TCanvas * c, TH1 * h, int style, int fillcolor =3001, int styleChange=0) {

   
   /*   h->Fit("gaus", "", "",  h->GetMean()-12, h->GetMean()+18);*/
  //  h->SetLineStyle( style) ;
  if (fillcolor!=3001)  h->SetLineColor( fillcolor) ;
  //  h->SetMarkerStyle(20 + style) ;
  // h->SetFillColor(style) ;
  //h->SetFillStyle(3004) ;
  if (styleChange==0)   h->SetFillStyle(0) ;
  if (styleChange==0)   h->SetLineStyle(3) ;


   c->Modified(); 
   c->Update();
 
   /* TPaveStats * stats = (TPaveStats*) c->GetPrimitive("stats");
  if (stats)
    { stats->SetName(h->GetTitle());
      stats->SetY1NDC(.4 + 10 * style);
      stats->SetY2NDC(.6 + 10 * style);

      stats->SetX1NDC(.5 + 10 * style);
      stats->SetX2NDC(.8 + 10 * style);
      //stats->Draw(); 
       c->Modified();
       c->Update();

    }
   */
 
 
  return;
}



void set_style(TEfficiency * h, int style, int fillcolor =3001) {
  h->SetLineStyle( style) ;
   if(fillcolor != 3001) {
    h->SetLineColor(0) ;
    h->SetFillColor(1) ;
    h->SetFillStyle(3004) ;
  }
  return;
}



void set_style(TProfile * h, int style, int fillcolor =3001) {
  h->SetLineStyle( style) ;
   if(fillcolor != 3001) {
    h->SetLineColor(0) ;
    h->SetFillColor(1) ;
    h->SetFillStyle(3004) ;
  }
  return;
}



void set_style(TGraphAsymmErrors * h, int style, int fillcolor =3001) {
  //  h->SetLineStyle( style) ;
  h->SetMarkerStyle(20 + style) ;
   if(fillcolor != 3001) {
    h->SetLineColor(0) ;
    h->SetFillColor(1) ;
    h->SetFillStyle(3004) ;
  
  }
  return;
}




// _____________________________________________________________________________
// This function makes 1 plot
void  MakePlot( TCanvas *c1, TLegend * leg1, TTree * tree1, TString var1,  TCut cut,
               TString title, int nbinsx, double xlow, double xup,
               TString plotname="plot", TString plotdir="",
		 TString options="plotLog:plotNorm", TString lab1 = "" , int style =1, int fillcolor = 3001 , TString labY = "Events") {


    std::clog << "MakePlots(): Plot vars: " << var1 << " " <<   std::endl;
    std::clog << "MakePlots(): Using cut: " << cut << std::endl;

    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));
    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    TH1F * h1         = new TH1F( lab1       , title, nbinsx, xlow, xup);
    tree1->Project( lab1, var1, cut);


    c1->SetLogy(plotLog);

    h1->SetLineWidth  (2);
    h1->SetLineColor  (2);
    //    h1->SetMarkerSize (0);
    h1->GetYaxis()->SetTitle(labY);
    h1->GetXaxis()->SetTitle(title);
    if (plotNorm) {
        h1->Scale(1.0 / h1->Integral());
        h1->GetYaxis()->SetTitle("arbitrary unit");

    }

    
    set_style( c1, h1, style , fillcolor);



   
    h1->Draw("sames hist");
    h1->Draw("sames e1");


    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");

    // Setup legends
    set_legstyle(leg1);
    if (fillcolor == 3001) {
      leg1->AddEntry(h1,  lab1, "l p");
    } else {
      leg1->AddEntry(h1,  lab1, "f");
    }

    leg1->Draw("same");       
    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    //delete h2;
    return ;
}





// _____________________________________________________________________________
// This function makes two plots overlaid
void  MakePlot2( TCanvas *c1, TLegend * leg1, TTree * tree1, TTree * tree2, TString var1, TString var2,  TCut cut, TCut cut2,
               TString title, int nbinsx, double xlow, double xup,
               TString plotname="plot", TString plotdir="",
		 TString options="plotLog:plotNorm", TString lab1 = "", TString lab2 = "" , int style =1, int styleChange =0 , int fillcolor = 3001 , TString labY = "Events") {


    std::clog << "MakePlots(): Plot vars: " << var1 << " " <<  var2 <<  std::endl;
    std::clog << "MakePlots(): Using cuts: " << cut << "and "<< cut2 << std::endl;

    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));
    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    TH1F * h1         = new TH1F( lab1       , title, nbinsx, xlow, xup);
    TH1F * h2         = new TH1F( lab2       , title, nbinsx, xlow, xup);
    tree1->Project( lab1, var1, cut);
    tree2->Project( lab2, var2, cut2);


    c1->SetLogy(plotLog);

    h1->SetLineWidth  (2);
    h1->SetLineColor  (2);
    h1->SetMarkerSize (0);
    h1->GetYaxis()->SetTitle(labY);
    h1->GetXaxis()->SetTitle(title);
    h2->SetLineWidth  (2);
    h2->SetLineColor  (4);
    h2->SetMarkerSize (0);
    h2->GetYaxis()->SetTitle(labY);
    h2->GetXaxis()->SetTitle(title);
    if (plotNorm) {
        h1->Scale(1.0 / h1->Integral());
        h1->GetYaxis()->SetTitle("arbitrary unit");
        h2->Scale(1.0 / h2->Integral());
        h2->GetYaxis()->SetTitle("arbitrary unit");
    }




   

    h2->Draw("sames hist");
    h2->Draw("sames e1");

    set_style(c1, h2, style  , fillcolor , !styleChange );

    h1->Draw("sames hist");
    h1->Draw("sames e1");
    set_style(c1, h1, style , fillcolor, !styleChange);


    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");

    // Setup legends
   set_legstyle(leg1);
      leg1->AddEntry(h1,  lab1, "l");
      leg1->AddEntry(h2,  lab2, "l");

    leg1->Draw("same");       
    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    //delete h2;
    return ;
}



// _____________________________________________________________________________
// This function makes 1 profile plot
void  MakeProfile( TCanvas *c1, TLegend * leg1, TTree * tree1,  TString var1,   TCut cut,
               TString title, int nbinsx, double xlow, double xup,
               TString plotname="plot", TString plotdir="",
		 TString options="plotLog:plotNorm", TString lab1 = "",  int style =1, TString labY = "Events") {


    std::clog << "MakePlots(): Plot vars: " << var1 << " "<<  std::endl;
    std::clog << "MakePlots(): Using cuts: " << cut  <<   std::endl;

    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));
    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    TProfile * h1         = new TProfile( lab1       , title, nbinsx, xlow, xup, 0 , 1.);
    tree1->Project( lab1, var1, cut , "profs");


    c1->SetLogy(plotLog);

    h1->SetLineWidth  (2);
    h1->SetLineColor  (2);
    h1->SetMarkerSize (0);
    h1->GetYaxis()->SetTitle(labY);
    h1->GetXaxis()->SetTitle(title);
    if (plotNorm) {
        h1->Scale(1.0 / h1->Integral());
        h1->GetYaxis()->SetTitle("arbitrary unit");
    }





   
    h1->Draw("sames hist");
    h1->Draw("sames e1");


    set_style( h1, style);

    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");

    // Setup legends
    set_legstyle(leg1);
    leg1->AddEntry(h1,  lab1, "l");

    leg1->Draw("sames");       
    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    //delete h2;
    return ;
}




// _____________________________________________________________________________
// This function makes two plots overlaid
void  MakeProfile2( TCanvas *c1, TLegend * leg1, TTree * tree1, TTree * tree2, TString var1, TString var2,  TCut cut,
               TString title, int nbinsx, double xlow, double xup,
               TString plotname="plot", TString plotdir="",
		 TString options="plotLog:plotNorm", TString lab1 = "", TString lab2 = "" , int style =1, TString labY = "Events") {


    std::clog << "MakePlots(): Plot vars: " << var1 << " " <<  var2 <<  std::endl;
    std::clog << "MakePlots(): Using cut: " << cut << std::endl;

    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));
    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    TProfile * h1         = new TProfile( lab1       , title, nbinsx, xlow, xup, 0, 1);
    TProfile * h2         = new TProfile( lab2       , title, nbinsx, xlow, xup, 0, 1);
    tree1->Project( lab1,   var1,  cut , "profs");
    tree2->Project(  lab2, var2,  cut, "profs");


    c1->SetLogy(plotLog);

    h1->SetLineWidth  (2);
    h1->SetLineColor  (2);
    h1->SetMarkerSize (0);
    h1->GetYaxis()->SetTitle(labY);
    h1->GetXaxis()->SetTitle(title);
    h2->SetLineWidth  (2);
    h2->SetLineColor  (4);
    h2->SetMarkerSize (0);
    h2->GetYaxis()->SetTitle(labY);
    h2->GetXaxis()->SetTitle(title);
    if (plotNorm) {
        h1->Scale(1.0 / h1->Integral());
        h1->GetYaxis()->SetTitle("arbitrary unit");
        h2->Scale(1.0 / h2->Integral());
        h2->GetYaxis()->SetTitle("arbitrary unit");
    }





   
    h1->Draw("sames hist");
    h1->Draw("sames e1");
    set_style(h1, style);

    h2->Draw("sames hist");
    h2->Draw("sames e1");



    set_style(h2, style);

    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");

    // Setup legends
    set_legstyle(leg1);
    leg1->AddEntry(h1,  lab1, "l");
    leg1->AddEntry(h2,  lab2, "l");

    leg1->Draw("same");       
    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    //delete h2;
    return ;
}





// _____________________________________________________________________________
// This function makes two plots overlaid
void  MakeEfficiencyPlot( TCanvas *c1, TLegend * leg1, TTree * tree1, TString var1,  TCut cutnum, TCut cutden,
               TString title, int nbinsx, double xlow, double xup,
               TString plotname="plot", TString plotdir="",
		 TString options="plotLog:plotNorm", TString lab1 = "" , int style =1, TString labY = "Events") {


    std::clog << "MakePlots(): Plot var: " << var1 << " "  <<  std::endl;
    std::clog << "MakePlots(): Using cuts: " << cutnum << " "<< cutden <<  std::endl;

    bool plotLog      = options.Contains("plotLog")    && (!options.Contains("!plotLog"));
    bool plotNorm     = options.Contains("plotNorm")   && (!options.Contains("!plotNorm"));

    TH1F * h1den         = new TH1F( "h1den"      , title, nbinsx, xlow, xup);
    TH1F * h2num         = new TH1F( "h2num"       , title, nbinsx, xlow, xup);
    tree1->Project( "h1den", var1, cutden );
    tree1->Project( "h2num",  var1, cutnum);
    //    TEfficiency * h1 = new TEfficiency(*h2num,*h1den, "CP");
    h1den->Sumw2();
    h2num->Sumw2();

    TGraphAsymmErrors * h1= new TGraphAsymmErrors (h2num , h1den, "cp") ;
    //    h1 = Divide(h_after_selection,h_before_selection,"cl=0.683 b(1,1) mode")

    double eff = h2num->GetEntries() / h1den ->GetEntries();
    double eff_err = sqrt(eff * (1. - eff) / h1den ->GetEntries());
    c1->SetLogy(plotLog);

    h1->SetLineWidth  (1);
    //    h1->SetLineColor  (2);
    h1->SetMarkerSize (2);
    h1->GetYaxis()->SetTitle(labY);
    h1->SetMaximum(1.2);
    h1->GetXaxis()->SetTitle(title);




    

    if (style ==1) {     h1->Draw("sames APE");}
    else{
    h1->Draw("sames PE");
    }

    set_style(h1, style);

    gPad->Print(plotdir+plotname+".png");
    gPad->Print(plotdir+plotname+".pdf");

    // Setup legends
    set_legstyle(leg1);
    leg1->AddEntry(h1,  lab1 + Form("; tot eff= %.3f #pm %.3f", eff, eff_err), "p");


    leg1->Draw("same");       
    // If not deleted, they will remain in the memory
    //delete c1;
    //delete h1;
    //delete h2;
    return ;
}

