#include "TROOT.h"
#include "TStyle.h"


#include "TLorentzVector.h"
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


#include <iostream>
#include <fstream>




#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline double deltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  while (result > TMath::Pi()) result -= 2* TMath::Pi();
  while (result <= -TMath::Pi()) result += 2* TMath::Pi();
  return result;
}

inline double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return sqrt(deta*deta + dphi*dphi);
}


inline float triggercorrMET(float met) {
  if (met < 110.)  return 0.940;
  if (met < 120.)  return 0.977;
  if (met < 130.)  return 0.984;
  if (met < 140.)  return 0.9820;
  if (met < 150.)  return 1.000;
  return 1.;
}


inline double evalHMETMassiveMt(double m0, double pt0, double phi0, double m1, double pt1, double phi1)
{

  return TMath::Sqrt(m0*m0 + m1*m1 + 2.0 * (TMath::Sqrt((m0*m0+pt0*pt0)*(m1*m1+pt1*pt1)) - pt0*pt1*TMath::Cos(deltaPhi(phi0, phi1)) ));

}



inline double evalHMETPt(double m0, double pt0, double phi0, double eta0, double m1, double pt1, double phi1, double eta1)

{

  TLorentzVector v1,v2;
  v1.SetPtEtaPhiM(pt0,eta0, phi0, m0);
  v2.SetPtEtaPhiM(pt1,eta1, phi1, m1);
  return (v1+v2).Pt();

}



inline double evalHMETPhi(double m0, double pt0, double phi0, double eta0, double m1, double pt1, double phi1, double eta1)

{

  TLorentzVector v1,v2;
  v1.SetPtEtaPhiM(pt0,eta0, phi0, m0);
  v2.SetPtEtaPhiM(pt1,eta1, phi1, m1);
  return (v1+v2).Phi();

}


inline float deltaPhiMETjets(float metphi, float jetphi, float jetpt, float jeteta, float minpt=25, float maxeta=4.5) {
  if(jetpt <= minpt || TMath::Abs(jeteta) >= maxeta)
    return 999.;
  else
    return fabs(deltaPhi(metphi, jetphi));
}


//______________________________________________________________________________                                                                                                                    
                                                                                                                                                                                                   



// To run, do "root -l plotHistos_Znn.C+"

void SignalAcceptance() {
    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    TH1::SetDefaultSumw2(1);
    //gROOT->SetBatch(1);

    TCanvas c1;

    TString plotdir = "plotsmonoH/";
    if (gSystem->AccessPathName(plotdir))
        gSystem->mkdir(plotdir);


    ////////////////////////////////////////////////////////////////////////////
    // Task 1 (a)                                                             //
    // - please comment out other tasks                                       //
    ////////////////////////////////////////////////////////////////////////////

    // Read from ntuples





     // TCut cutmc_all = "MinIf$(abs(deltaPhi(V_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>1.5";

     TCut cutMinimal = "(Vtype==2 || Vtype==3 || Vtype==4) &&  HCSV_pt>0 && met_pt>150  && Jet_btagCSV[hJidx[1]]>0.3 && min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30";

         TCut cutQCD = "met_pt>150 &&  MinIf$( abs( deltaPhi(met_phi,Jet_phi)), Jet_pt>30 && Jet_puId &&  abs(Jet_eta)<4.5 )> 0.5 && abs( deltaPhi(met_phi,tkMet_phi))<0.7 &&   tkMet_pt>0  && min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30";


	 // TCut cutQCD = "";



     TCut addCenJet30m0 = "(Sum$(Jet_pt>30 && Jet_puId &&  abs(Jet_eta)<4.5)-2)>0";
     TCut addCenJet30e0 = "(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==0";
     TCut addCenJet30le1 = "(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)<=1";
     TCut addCenJet30e1 = "(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==1";
     TCut naddGoodLeptons10m0= "(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 ||  aLeptons_jetDR>0.3 ))+Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3 )))>0";
     TCut naddGoodLeptons10e0= "(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 ||  aLeptons_jetDR>0.3 ))+Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3 )))==0";
     TCut naddGoodTaus20e0= "Sum$(TauGood_idDecayMode >= 1 && TauGood_idCI3hit >= 1 && TauGood_pt > 20. && abs(TauGood_eta) < 2.3)==0";




     TCut cutSignal0J = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  Jet_btagCSV[hJidx[1]]>0.605 && (HCSV_mass>90 && HCSV_mass<150)" + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1;


    TCut weightmc = "";
    // TCut weightmc = "(efflumi * (1.26/3.)*  sign(genWeight))";
    TCut weightdata = "";

    TCut met = "met_pt>150";
    TCut jets = "min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30";
    TCut btag = "min(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]] )>0.46";
    TCut mass = "(H_mass>90 && H_mass<150)";
    //    TCut mt = "evalHMETMassiveMt(HCSV_mass, HCSV_pt, HCSV_phi, 91.18,  met_pt, met_phi)>500";
 
    TH1D * Acceptance = new TH1D("Acceptance", "Acceptance", 7, 0, 7);
    
    Acceptance->GetXaxis()->SetBinLabel(1, "met>150+trigger");
    Acceptance->GetXaxis()->SetBinLabel(2, "pt_{1,2}>30,|#eta_{1,2}|<2.4");
    Acceptance->GetXaxis()->SetBinLabel(3, "QCD cuts");
    Acceptance->GetXaxis()->SetBinLabel(4, "0 additional leptons");
    Acceptance->GetXaxis()->SetBinLabel(5, "<=1 additional jets");
    Acceptance->GetXaxis()->SetBinLabel(6, "btag");
    Acceptance->GetXaxis()->SetBinLabel(7, "mass cut");
    //    Acceptance->GetXaxis()->SetBinLabel(8, "mt cut");


    TFile * file = TFile::Open("/afs/cern.ch/user/d/degrutto/scratch3/VHbbRun2/CMSSW_7_6_3_patch2/src/run2_300316_V21/macros_mergeVj_H_pt/eos/cms/store/group/cmst3/user/degrutto/VHBBHeppyV21_add/ZprimeToA0hToA0chichihbb_2HDM_MZp-800_MA0-300_13TeV-madgraph.root");
    TFile * file_600 = TFile::Open("/afs/cern.ch/user/d/degrutto/scratch3/VHbbRun2/CMSSW_7_6_3_patch2/src/run2_300316_V21/macros_mergeVj_H_pt/eos/cms/store/group/cmst3/user/degrutto/VHBBHeppyV21_add/ZprimeToA0hToA0chichihbb_2HDM_MZp-600_MA0-300_13TeV-madgraph.root");
    TFile * file_1000 = TFile::Open("/afs/cern.ch/user/d/degrutto/scratch3/VHbbRun2/CMSSW_7_6_3_patch2/src/run2_300316_V21/macros_mergeVj_H_pt/eos/cms/store/group/cmst3/user/degrutto/VHBBHeppyV21_add/ZprimeToA0hToA0chichihbb_2HDM_MZp-1000_MA0-300_13TeV-madgraph.root");


    TTree * tree = (TTree *) file->Get("tree"); 
    TH1D * count = (TH1D *) file->Get("Count"); 
    std::cout << "total count is " << count->GetBinContent(1) << std::endl;

    int nStart = count->GetBinContent(1) ;
    TH1D * htemp = new TH1D("htemp", "htemp", 10000, 0, 10000) ;
    tree->Project("htemp", "met_pt" , met ) ;
   
    std::cout << "met " << htemp->Integral() << std::endl;

    Acceptance->SetBinContent(1, double (htemp->Integral() / nStart ) );

    tree->Project("htemp", "met_pt" ,  met  + jets  ) ;

    std::cout << "after j1j2 |#eta|<2.4, pt>30 " << htemp->Integral() << std::endl;

    Acceptance->SetBinContent(2, double (htemp->Integral() / nStart ) );




    tree->Project("htemp", "met_pt" ,  met  + jets + cutQCD  ) ;

    std::cout << "QCD cuts " << htemp->Integral() << std::endl;

    Acceptance->SetBinContent(3, double (htemp->Integral() / nStart ) );

    tree->Project("htemp", "met_pt" ,   met  + jets + naddGoodLeptons10e0 + naddGoodTaus20e0  + cutQCD ) ;

    std::cout << "0 additionla lepton " << htemp->Integral() << std::endl;


    Acceptance->SetBinContent(4, double (htemp->Integral() / nStart ) );

    tree->Project("htemp", "met_pt" ,  met  + jets  + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

    std::cout << "at max 1  additionljets " << htemp->Integral() << std::endl;

    Acceptance->SetBinContent(5, double (htemp->Integral() / nStart ) );

    tree->Project("htemp", "met_pt" , met + jets + btag  + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

    std::cout << "btag " << htemp->Integral() << std::endl;

    Acceptance->SetBinContent(6, double (htemp->Integral() / nStart ) );
    tree->Project("htemp", "met_pt" ,  met + jets + btag +  mass +   naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

    std::cout << "mass " << htemp->Integral() << std::endl;
    Acceptance->SetBinContent(7, double (htemp->Integral() / nStart ) );


    tree->Project("htemp", "met_pt" ,   met + jets + btag +  mass   +  naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;


    Acceptance->SetBinContent(8, double (htemp->Integral() / nStart ) );







    TH1D * Acceptance_600 = new TH1D("Acceptance_600", "Acceptance_600", 7, 0, 7);
    
    Acceptance_600->GetXaxis()->SetBinLabel(1, "met>150+trigger");
    Acceptance_600->GetXaxis()->SetBinLabel(2, "pt_{1,2}>30,|#eta_{1,2}|<2.4,hpt>150");
    Acceptance_600->GetXaxis()->SetBinLabel(3, "QCD cuts");
    Acceptance_600->GetXaxis()->SetBinLabel(4, "0 additional leptons");
    Acceptance_600->GetXaxis()->SetBinLabel(5, "<=1 additional jets");
    Acceptance_600->GetXaxis()->SetBinLabel(6, "btag");
    Acceptance_600->GetXaxis()->SetBinLabel(7, "mass cut");
    //    Acceptance_600->GetXaxis()->SetBinLabel(8, "mt cut");



    TTree * tree_600 = (TTree *) file_600->Get("tree"); 
    TH1D * count_600 = (TH1D *) file_600->Get("Count"); 
    std::cout << "total count is " << count_600->GetBinContent(1) << std::endl;

    nStart = count_600->GetBinContent(1) ;
    tree_600->Project("htemp", "met_pt" , met ) ;
   
    std::cout << "met " << htemp->Integral() << std::endl;

    Acceptance_600->SetBinContent(1, double (htemp->Integral() / nStart ) );

    tree_600->Project("htemp", "met_pt" ,  met  + jets  ) ;

    std::cout << "after j1j2 |#eta|<2.4, pt>30 " << htemp->Integral() << std::endl;

    Acceptance_600->SetBinContent(2, double (htemp->Integral() / nStart ) );




    tree_600->Project("htemp", "met_pt" ,  met  + jets + cutQCD  ) ;

    std::cout << "QCD cuts " << htemp->Integral() << std::endl;

    Acceptance_600->SetBinContent(3, double (htemp->Integral() / nStart ) );

    tree_600->Project("htemp", "met_pt" ,   met  + jets + naddGoodLeptons10e0 + naddGoodTaus20e0  + cutQCD ) ;

    std::cout << "0 additionla lepton " << htemp->Integral() << std::endl;


    Acceptance_600->SetBinContent(4, double (htemp->Integral() / nStart ) );

    tree_600->Project("htemp", "met_pt" ,  met  + jets  + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

    std::cout << "at max 1  additionljets " << htemp->Integral() << std::endl;

    Acceptance_600->SetBinContent(5, double (htemp->Integral() / nStart ) );

    tree_600->Project("htemp", "met_pt" , met + jets + btag  + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

    std::cout << "btag " << htemp->Integral() << std::endl;

    Acceptance_600->SetBinContent(6, double (htemp->Integral() / nStart ) );
    tree_600->Project("htemp", "met_pt" ,  met + jets + btag +  mass +   naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

    std::cout << "mass " << htemp->Integral() << std::endl;
    Acceptance_600->SetBinContent(7, double (htemp->Integral() / nStart ) );


    tree_600->Project("htemp", "met_pt" ,   met + jets + btag +  mass   +  naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

       std::cout << "mt " <<     htemp->Integral()  << std::endl;
    Acceptance_600->SetBinContent(8, double (htemp->Integral() / nStart ) );



    TH1D * Acceptance_1000 = new TH1D("Acceptance_1000", "Acceptance_1000", 7, 0, 7);
    
    Acceptance_1000->GetXaxis()->SetBinLabel(1, "met>150+trigger");
    Acceptance_1000->GetXaxis()->SetBinLabel(2, "pt_{1,2}>30,|#eta_{1,2}|<2.4,hpt>150");
    Acceptance_1000->GetXaxis()->SetBinLabel(3, "QCD cuts");
    Acceptance_1000->GetXaxis()->SetBinLabel(4, "0 additional leptons");
    Acceptance_1000->GetXaxis()->SetBinLabel(5, "<=1 additional jets");
    Acceptance_1000->GetXaxis()->SetBinLabel(6, "btag");
    Acceptance_1000->GetXaxis()->SetBinLabel(7, "mass cut");
    //    Acceptance_1000->GetXaxis()->SetBinLabel(8, "mt cut");



    TTree * tree_1000 = (TTree *) file_1000->Get("tree"); 
    TH1D * count_1000 = (TH1D *) file_1000->Get("Count"); 
    std::cout << "total count is " << count_1000->GetBinContent(1) << std::endl;

    nStart = count_1000->GetBinContent(1) ;
    tree_1000->Project("htemp", "met_pt" , met ) ;
   
    std::cout << "met " << htemp->Integral() << std::endl;

    Acceptance_1000->SetBinContent(1, double (htemp->Integral() / nStart ) );

    tree_1000->Project("htemp", "met_pt" ,  met  + jets  ) ;

    std::cout << "after j1j2 |#eta|<2.4, pt>30 " << htemp->Integral() << std::endl;

    Acceptance_1000->SetBinContent(2, double (htemp->Integral() / nStart ) );




    tree_1000->Project("htemp", "met_pt" ,  met  + jets + cutQCD  ) ;

    std::cout << "QCD cuts " << htemp->Integral() << std::endl;

    Acceptance_1000->SetBinContent(3, double (htemp->Integral() / nStart ) );

    tree_1000->Project("htemp", "met_pt" ,   met  + jets + naddGoodLeptons10e0 + naddGoodTaus20e0  + cutQCD ) ;

    std::cout << "0 additionla lepton " << htemp->Integral() << std::endl;


    Acceptance_1000->SetBinContent(4, double (htemp->Integral() / nStart ) );

    tree_1000->Project("htemp", "met_pt" ,  met  + jets  + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

    std::cout << "at max 1  additionljets " << htemp->Integral() << std::endl;

    Acceptance_1000->SetBinContent(5, double (htemp->Integral() / nStart ) );

    tree_1000->Project("htemp", "met_pt" , met + jets + btag  + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

    std::cout << "btag " << htemp->Integral() << std::endl;

    Acceptance_1000->SetBinContent(6, double (htemp->Integral() / nStart ) );
    tree_1000->Project("htemp", "met_pt" ,  met + jets + btag +  mass +   naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

    std::cout << "mass " << htemp->Integral() << std::endl;
    Acceptance_1000->SetBinContent(7, double (htemp->Integral() / nStart ) );


    tree_1000->Project("htemp", "met_pt" ,   met + jets + btag +  mass   +  naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1 + cutQCD) ;

    std::cout << "mt " <<     htemp->Integral()  << std::endl;
    Acceptance_1000->SetBinContent(8, double (htemp->Integral() / nStart ) );


    //    Acceptance->Draw("TEXT00");
    Acceptance_1000->SetLineColor(kBlue);
    Acceptance_1000->Draw("TEXT00");
    Acceptance_1000->SetMinimum(0);
    Acceptance_1000->SetMaximum(1);

   Acceptance->Draw("TEXT00SAME");
    Acceptance_600->SetLineColor(kRed);
    Acceptance_600->Draw("TEXT00SAME");


    // Setup legends                                                                                                                    
    TLegend * leg1 = new TLegend(0.70, 0.68, 0.82, 0.92);

 
    leg1->AddEntry(Acceptance_600, "MZp600", "l");
    leg1->AddEntry(Acceptance, "MZp800", "l");
    leg1->AddEntry(Acceptance_1000, "MZp1000", "l");

    leg1->Draw("same");


  

    c1.SaveAs("CutFlowEfficiency.pdf");


}


