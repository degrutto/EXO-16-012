#include "plotHistos_boostedZH_13TeV_monoH.h"
#include "TROOT.h"
#include "TStyle.h"







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

void plotHistos_test() {
    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    TH1::SetDefaultSumw2(1);
    //gROOT->SetBatch(1);

    TString plotdir = "plotsmonoH/";

    if (gSystem->AccessPathName(plotdir))
        gSystem->mkdir(plotdir);


    ////////////////////////////////////////////////////////////////////////////
    // Task 1 (a)                                                             //
    // - please comment out other tasks                                       //
    ////////////////////////////////////////////////////////////////////////////

    // Read from ntuples





     // TCut cutmc_all = "MinIf$(abs(deltaPhi(V_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>1.5";

     TCut cutMinimal = "(Vtype==2 || Vtype==3 || Vtype==4) &&  HCSV_pt>150 && met_pt>150  && Jet_btagCSV[hJCidx[1]]>0.3 && min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 && HCSV_mass<300";

     TCut cutQCD = "HCSV_pt>150 && met_pt>150 &&  MinIf$( abs( deltaPhi(met_phi,Jet_phi)), Jet_pt>30 && Jet_puId &&  abs(Jet_eta)<4.5 )> 0.7 && abs( deltaPhi(met_phi,tkMet_phi))<0.7 &&   tkMet_pt>0 && Jet_btagCSV[hJCidx[1]]>0.3 && min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 && tkMet_pt>30";


     TCut cutAntiQCD = "HCSV_pt>150 && met_pt>150 &&   MinIf$( abs( deltaPhi(met_phi,Jet_phi)), Jet_pt>30 && Jet_puId &&  abs(Jet_eta)<4.5 )< 0.7  && abs( deltaPhi(met_phi,tkMet_phi))>0.7 && Jet_btagCSV[hJCidx[1]]>0.3 && min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 && tkMet_pt<30";




     TCut addCenJet30m0 = "(Sum$(Jet_pt>30 && Jet_puId &&  abs(Jet_eta)<4.5)-2)>0";
     TCut addCenJet30e0 = "(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==0";
     TCut addCenJet30le1 = "(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)<=1";
     TCut addCenJet30e1 = "(Sum$(Jet_pt>30 && Jet_puId && abs(Jet_eta)<4.5)-2)==1";
     TCut naddGoodLeptons10m0= "(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 ||  aLeptons_jetDR>0.3 ))+Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3 )))>0";
     TCut naddGoodLeptons10e0= "(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 ||  aLeptons_jetDR>0.3 ))+Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3 )))==0";
     TCut naddGoodTaus20e0= "Sum$(TauGood_idDecayMode >= 1 && TauGood_idCI3hit >= 1 && TauGood_pt > 20. && abs(TauGood_eta) < 2.3)==0";



     TCut cutSignalLoose = "Vtype==4 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605 && (HCSV_mass<90 || HCSV_mass>150)" + naddGoodLeptons10e0 + naddGoodTaus20e0;

     TCut cutSignal0J = "Vtype==4 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605 && (HCSV_mass>90 && HCSV_mass<150)" + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30le1;
     TCut cutSignal1J = "Vtype==4 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605 && (HCSV_mass>90 && HCSV_mass<150)" + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30e1;
     TCut cutSignalSB = "Vtype==4 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605 && (HCSV_mass<90 || HCSV_mass>150)" + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30e0;

     TCut cutSignal1b = "Vtype==4 &&  min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30 &&  Jet_btagCSV[hJCidx[1]]>0.605 && Jet_btagCSV[hJCidx[0]]<0.97 && ((HCSV_mass>90 && HCSV_mass<150))" + naddGoodLeptons10e0 + naddGoodTaus20e0 + addCenJet30e0;


     TCut cutTT = "(Vtype==2 || Vtype==3) && vLeptons_pt>30" + addCenJet30m0 + /*naddGoodLeptons5m0 + */"Jet_btagCSV[hJCidx[0]]>0.97 && Jet_btagCSV[hJCidx[1]]<0.97";

     
     TCut cutZlight = "Vtype==4" + addCenJet30e0 + naddGoodLeptons10e0 +  "Jet_btagCSV[hJCidx[0]]<0.97";



     TCut cutZbb = "Vtype==4 && (HCSV_mass<100 || HCSV_mass>140)" + addCenJet30e0 +  naddGoodLeptons10e0 +  "Jet_btagCSV[hJCidx[1]]>0.89";

     TCut cutWlight = "(Vtype==2 || Vtype==3) && vLeptons_pt>30" + addCenJet30e0 + /*naddGoodLeptons5m0 +*/ "Jet_btagCSV[hJCidx[0]]<0.97";
     TCut cutWbb = "(Vtype==2 || Vtype==3) && vLeptons_pt>30" + addCenJet30e0 + /*naddGoodLeptons5m0 +*/ "Jet_btagCSV[hJCidx[1]]>0.89";
  
     //TCut cutmc_all = "Vtype==4  &&  t[hJCidx[0]]>30 ))  &&  Sum$(aLeptons_pt>5  && aLeptons_relIso03<1. && ( deltaR(Jet_eta[hJCidx[0]], aLeptons_eta, Jet_phi[hJCidx[0]], aLeptons_phi)>0.4 &&  deltaR(Jet_eta[hJCidx[1]], aLeptons_eta, Jet_phi[hJCidx[1]], aLeptons_phi)>0.4))>=0 && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId==1)<=99  && abs(deltaPhi(HCSV_phi, V_phi))>2.5";
     /*
     Higgs candidate:             Vtype >=0
       MET>150:                       met_pt>150
       Higgs pT>150:                HCSV_reg_pt>150
       CSV2>0.2:                      Jet_btagCSV[hJCidx[1]]>0.2
       DeltaPhi(H,MET)>1.5      abs(TVector2::Phi_mpi_pi(HCSV_reg_phi-met_phi))>1.5
       Higgs jet with pt>30       min(Jet_pt[hJCidx[0]],Jet_pt[hJCidx[1]])>30
       H mass <300                 HCSV_reg_mass<300 && HCSV_mass>0
     */

       /*


Nothing = 1 && <!Cuts|FlagsMET!>
NoQCD1 = <!Cuts|Nothing!> && HCSV_reg_pt>150 && minDeltaPhiJet2Met>0.7
NoQCD2 = <!Cuts|NoQCD1!> && minDeltaPhiCenJetNoPU30Met >0.5
NoQCD3 = <!Cuts|NoQCD2!> && tkMetPVchs_pt>25 && tkMet_pt>25
       ;metType1p2_pt/sqrt(met_sumEt)>5 && DeltaPhiJet1Jet2<2.5



CR_Nothing = <!Cuts|NoQCD3!>
    CR_TTbar = <!Cuts|NoQCD3!> && naddGoodLeptons5>=1 && (addCenJet30)>=1 && Jet_btagCSV[hJCidx[0]]>0.94
    CR_WLight = <!Cuts|NoQCD3!> && naddGoodLeptons5>=1 && (addCenJet30)<=1 && Jet_btagCSV[hJCidx[0]]<0.94
CR_Wbb = <!Cuts|NoQCD3!> && naddGoodLeptons5>=1 && (addCenJet30)<=1 && Jet_btagCSV[hJCidx[0]]>0.94									 
CR_ZLight = <!Cuts|NoQCD3!> && naddGoodLeptons5==0 && (addCenJet30)<=1 && Jet_btagCSV[hJCidx[0]]<0.94	
CR_Zbb = <!Cuts|NoQCD3!> && naddGoodLeptons5==0 && (addCenJet30)<=1 && Jet_btagCSV[hJCidx[0]]>0.94 && (H.mass<110 || H.mass>140)
CR_QCD = <!Cuts|Nothing!> && minDeltaPhiJet2Met<0.7

       */
     TCut cutdata = "(json==1 && ( HLT_BIT_HLT_PFMET170_NoiseCleaned_v ||  HLT_BIT_HLT_PFMET120_PFMHT120_IDTight_v))";
     //     TCut trigMC = "(HLT_BIT_HLT_PFMET170_NoiseCleaned_v || HLT_BIT_HLT_PFMET120_PFMHT120_IDLoose_v)";
     TCut trigMC = "";

    TString channel  = "monoH_Zp800";

    //    ev->read(   cutMinimal   ,  cutMinimal, "" , "tree"); 


    TCut weightmc = "(efflumi * (1.28/3.) * sign(genWeight)) * puWeight * weightTrig";
    // TCut weightmc = "(efflumi * (1.26/3.)*  sign(genWeight))";
    TCut weightdata = "";



    // Read from ntuples
    Events * ev = new Events();
    TString process  = "data_obs";         // can try other processes: VH, WJ, ZJ, TT
    ev->read(cutQCD, cutQCD, process);

    TString plotname = process + "_met";      // the name of the image file
    int nbinsx       = 15;                      // number of bins
    double xlow      = 150.0;                    // the low edge of x-axis
    double xup       = 300.0;                   // the upper edge of x-axis
    TString options  = "!plotLog:!plotNorm";    // use "plotLog" to plot on log-y scale,
                                                // use "!plotLog" to plot on linear scale;
                                                // use "plotNorm" to plot normalized plots,
                            


      //   MakePlots(ev, varN[i], (cutQCD + trigMC + cutZlight ) * weightmc , (cutQCD + cutZlight + cutdata) , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "ZlightLin", plotdir, options);

    
      //MakePlot(ev->ZH, "H.pt", cut, "; p_{T}(jj) [GeV]", 16, 0, 240., process + "_Hpt", plotdir, options);

         MakePlot(ev->data_obs, "met_pt", cutdata + cutZlight + cutQCD, ";met [GeV]", 20, 150, 350., "met", plotdir , options);




}

//  LocalWords:  deltaPhi
