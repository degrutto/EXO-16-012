#include "compare.h"
#include "TMath.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TFileCollection.h"

// To run, do "root -l plotHistos_Znn.C+"

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



inline double evalm1m2Pt(double m0, double pt0, double phi0, double eta0, double m1, double pt1, double phi1, double eta1)

{

  TLorentzVector v1,v2;
  v1.SetPtEtaPhiM(pt0,eta0, phi0, m0);
  v2.SetPtEtaPhiM(pt1,eta1, phi1, m1);
  return (v1+v2).Pt();

}


inline double evalm1m2Phi(double m0, double pt0, double phi0, double eta0, double m1, double pt1, double phi1, double eta1)

{

  TLorentzVector v1,v2;
  v1.SetPtEtaPhiM(pt0,eta0, phi0, m0);
  v2.SetPtEtaPhiM(pt1,eta1, phi1, m1);
  return (v1+v2).Phi();

}




void compare() {
     gROOT->LoadMacro("tdrstyle.C");
     gROOT->ProcessLine("setTDRStyle()");

    TH1::SetDefaultSumw2(1);
    //gROOT->SetBatch(1);

     TString plotdir = "plots/";
    


    if (gSystem->AccessPathName(plotdir))
        gSystem->mkdir(plotdir);


  TCanvas * c1 = new TCanvas("c1", "c1", 700, 700);
  TLegend * leg = new TLegend(0.70, 0.68, 0.9, 0.92);
 

  //      TFile * file1 = TFile::Open("/afs/cern.ch/user/d/degrutto/eos/cms/store/group/cmst3/user/degrutto/MonoHiggsWP_M800/MonoHiggsWP_M800.root");
      TFile * file1 = TFile::Open("/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_v14/step3/Step3_ZJets.root");

//  TFile * file2 = TFile::Open("/afs/cern.ch/user/d/degrutto/eos/cms/store/group/cmst3/user/degrutto/MonoHiggsWP_M1500/MonoHiggsWP_M1500.root");
      ///File * file2 = TFile::Open("/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_v14/step3/Step3_ZJets.root");
	 //TFile * file2 = TFile::Open("M1500.root");
      //	 TFile * file2 = TFile::Open("/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_v14/SingleMuon__Run2015*/skim_SingleMuon__Run2015*tree_*");


      TChain * tree2 = new TChain("tree");

      //      tree2->AddFile("/afs/cern.ch/user/d/degrutto/eos/cms/store/group/phys_higgs/hbb/ntuples/V14/SingleMuon/VHBB_HEPPY_V14_SingleMuon__Run2015*/*/*/tree_*.root");
      //      tree2->AddFileInfoList("/afs/cern.ch/user/d/degrutto/eos/cms/store/group/phys_higgs/hbb/ntuples/V14/SingleMuon/VHBB_HEPPY_V14_SingleMuon__Run2015*/*/*/tree_*.root");


      TFileCollection * fc =new TFileCollection("dum","","singleMuonFiles.txt");
      tree2->AddFileInfoList(  ( TCollection * ) fc->GetList());


      TTree * tree1 = (TTree *) file1->Get("tree"); 


  

      TCut cutQCD = "HCSV_pt>150 && met_pt>150 &&  MinIf$( abs( deltaPhi(met_phi,Jet_phi)), Jet_pt>30 && Jet_puId &&  abs(Jet_eta)<4.5 )> 0.7 && abs( deltaPhi(met_phi,tkMet_phi))<0.7 &&   tkMet_pt>0  && Jet_btagCSV[hJCidx[1]]>0.3 && min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 && tkMet_pt>30";


      TCut cutQCDMM = "HCSV_pt>150 && fakeMET_pt>150";


      TCut cutAntiQCD = "HCSV_pt>150 && met_pt>150 &&   MinIf$( abs( deltaPhi(met_phi,Jet_phi)), Jet_pt>30 && Jet_puId &&  abs(Jet_eta)<4.5 )< 0.7  && abs( deltaPhi(met_phi,tkMet_phi))>0.7 && Jet_btagCSV[hJCidx[1]]>0.3 && min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 && tkMet_pt<30";


      TCut addCenJet30m0 = "(Sum$(Jet_pt>30 && Jet_puId &&  abs(Jet_eta)<4.5)-2)>0";
      TCut addCenJet30e0 = "(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==0";
      TCut addCenJet30le1 = "(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)<=1";
      TCut addCenJet30e1 = "(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==1";
      TCut naddGoodLeptons10m0= "(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 ||  aLeptons_jetDR>0.3 ))+Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3 )))>0";
      TCut naddGoodLeptons10e0= "(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 ||  aLeptons_jetDR>0.3 ))+Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3 )))==0";
      TCut naddGoodTaus20e0= "Sum$(TauGood_idDecayMode >= 1 && TauGood_idCI3hit >= 1 && TauGood_pt > 20. && abs(TauGood_eta) < 2.3)==0";


     TCut cutSignalLoose = "Vtype==4 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605 && (HCSV_mass<90 || HCSV_mass>150)" + naddGoodLeptons10e0 + naddGoodTaus20e0;

     TCut cutSignal0J = "Vtype==4 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605" + naddGoodLeptons10e0 + naddGoodTaus20e0;

     TCut cut1b0J = "Vtype==4 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605 && Jet_btagCSV[hJCidx[0]]<0.97 &&  (HCSV_mass>90 && HCSV_mass<150)" + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30e0;


     TCut cutSignal1J = "Vtype==4 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605 && (HCSV_mass>90 && HCSV_mass<150)" + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30e1;
     TCut cutSignalSB = "Vtype==4 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605 && (HCSV_mass<90 || HCSV_mass>150)" + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30e0;


     TCut cutSignalMM0J = "json==1 && Vtype==0 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605" ;


     TCut cutSignalMMSB = "json==1 && Vtype==0 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605 && !(HCSV_mass>90 && HCSV_mass<150)" ;
     
     
     
     TString options  = "!plotLog:plotNorm";    // use "plotLog" to plot on log-y scale,                                                                                                          


     TCut cutdata = "json==1 && ( HLT_BIT_HLT_PFMET170_NoiseCleaned_v ||  HLT_BIT_HLT_PFMET120_PFMHT120_IDTight_v)";
     //     TCut trigMC = "(HLT_BIT_HLT_PFMET170_NoiseCleaned_v || HLT_BIT_HLT_PFMET120_PFMHT120_IDLoose_v)";
     TCut trigMC = "";

     TString channel  = "monoH_Zp800";

   

     TCut weightmc = "(efflumi *  sign(genWeight)) * puWeight * weightTrig";
     // TCut weightmc = "(efflumi * (1.26/3.)*  sign(genWeight))";
     TCut weightdata = "";




   
      // use "!plotLog" to plot on linear scale;                                                                                                           
      // use "plotNorm" to plot normalized plots,                                                                                                          
      // use "!plotNorm" otherwise.                              


    // Using "ev->ZH" for ZH and "ev->ZjHF" for Z+HF
      //      MakePlot(c1, tree2, var, cut, title, nbinsx, xlow, xup, plotname+"1", plotdir, options );
      


      
      TString var = "evalHMETMassiveMt(HCSV_mass, HCSV_pt, HCSV_phi, 91.18,  met_pt, V_phi)";
      TString title    = "MT [GeV]";          // the title of the histogram
      TString plotname = "mt";      // the name of the image file
      int nbinsx       = 15;                      // number of bins
      double xlow      = 0;                    // the low edge of x-axis
      double xup       = 1500;                   // the upper edge of x-axis
      //      MakePlot2(c1, leg, tree1, tree2, var,  var, (cutQCD + trigMC + cutSignal0J) * weightmc ,  (cutQCD + trigMC + cutSignalSB) * weightmc  , title, nbinsx, xlow, xup, plotname + "compareMTshape", plotdir, options, "ZJets,SR0J", "ZJets,SB0J");

  TLegend * leg2 = new TLegend(0.70, 0.48, 0.9, 0.72);

  //  MakePlot2(c1, leg, tree1, tree1, var,  var, (cutQCD + trigMC + cutSignal0J) * weightmc ,  (cutQCD + trigMC + cut1b0J) * weightmc  , title, nbinsx, xlow, xup, plotname + "compareMTshape", plotdir, options, "ZnunuJets_SR0J", "ZnunuJets_1b0J");



 TString  varMM = "evalHMETMassiveMt(HCSV_mass, HCSV_pt, HCSV_phi, 0., fakeMET_pt, fakeMET_phi )";

 //MakePlot2(c1, leg, tree1, tree2, var,  varMM,   (  trigMC + cutSignal0J) * weightmc   ,  cutQCDMM + cutSignalMM0J   , title, nbinsx, xlow, xup, plotname + "compareMTshape", plotdir, options, "ZnunuJets_SR", "doubleMuon_fakeMET_ZmmJets");

     title    = "MET [GeV]";          // the title of the histogram
       plotname = "met";      // the name of the image file
       nbinsx       = 10;                      // number of bins
       xlow      = 150;                    // the low edge of x-axis
       xup       = 450;                   // the upper edge of x-axis


       MakePlot2(c1, leg, tree1, tree2, "met_pt",  "fakeMET_pt",   (  trigMC + cutSignal0J) * weightmc   ,  cutQCDMM + cutSignalMM0J   , title, nbinsx, xlow, xup, plotname + "compareMETshape", plotdir, options, "ZnunuJets_SR", "doubleMuon_fakeMET_ZmmJets");






}
