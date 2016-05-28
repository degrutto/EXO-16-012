#include "plotHistos_boostedZH_13TeV_monoH_shape.h"
#include "TROOT.h"
#include "TStyle.h"
#include "puWeight.h"

#include "BTagCalibrationStandalone.h"
// setup calibration readers
BTagCalibration calib("CSVv2", "CSVv2.csv");
BTagCalibrationReader readerHF(&calib,               // calibration instance
                             BTagEntry::OP_LOOSE,  // operating point
                             "mujets",               // measurement type
                             "central");           // systematics type

BTagCalibrationReader readerLF(&calib,               // calibration instance
                             BTagEntry::OP_LOOSE,  // operating point
                             "incl",               // measurement type
                             "central");           // systematics type



inline double weightBtag(double pt, double eta, unsigned int flav) { // assuming CSVM
  
  double SF=1.0; 
  float MaxBJetPt = 670., MaxLJetPt = 1000.;
  if (pt>MaxBJetPt && (flav==5 || flav==4) )  { // use MaxLJetPt for  light jets
    pt = MaxBJetPt -0.01; 
  }
  if (pt>MaxLJetPt && !(flav==5 || flav==4))  { // use MaxLJetPt for  light jets
    pt = MaxLJetPt - 0.01; 
  }
  /*  enum JetFlavor {
    FLAV_B=0,
    FLAV_C=1,
    FLAV_UDSG=2,
  };
  */
  switch(flav)
    {
    case 5:
      SF = readerHF.eval( BTagEntry::FLAV_B , eta, pt); 
      break;
    case 4:
      SF = readerHF.eval( BTagEntry::FLAV_C , eta, pt); 
      break; 
    
    default:
      SF = readerLF.eval( BTagEntry::FLAV_UDSG , eta, pt); 
    }
  //  std::cout << "SF for pt, " << pt <<" eta, " << eta << " flav " << flav << " is " << SF << std::endl; 
  return SF ;


}
  /*
  if (btag>0.89) { // batg eff sf for CSVM
    if (pt<50)
      return 0.944;
    else if (pt<70) 
      return 0.939;
    else if  (pt<100)
      return 0.878;
    else 
      return 0.899;
  } 
  else {
    return 1.1;
  }
 
  return 1.0;
  */







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


inline double evalHMETMassiveMt(double m0, double pt0, double phi0, double eta0,double m1, double pt1, double eta1, double phi1)
{

  return TMath::Sqrt(m0*m0 + m1*m1 + 2.0 * (TMath::Sqrt((m0*m0+pt0*pt0)*(m1*m1+pt1*pt1)) - pt0*pt1*TMath::Cos(deltaPhi(phi0, phi1)) ));
  //TLorentzVector v1,v2;
  // v1.SetPtEtaPhiM(pt0,eta0, phi0, m0);
  //v2.SetPtEtaPhiM(pt1,eta1, phi1, m1);
  // return (v1+v2).Mt();


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

void plotHistos_boostedZH_13TeV_monoH_shape() {
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
    Events * ev = new Events();




     // TCut cutmc_all = "MinIf$(abs(deltaPhi(V_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>1.5";

     TCut cutMinimal = "mhtJet30>120 && Jet_btagCSV[hJCidx[1]]>0.4 && HCSV_mass<500 && abs(deltaPhi(HCSV_phi,met_phi))>0.7 && json && (Vtype==2 || Vtype==3 || Vtype==4) &&  H_reg_pt>0 && met_pt>150  && min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 && Jet_btagCSV[hJidx[1]]>0.46 && Jet_btagCSV[hJidx[0]]>0.46 &&  Jet_id[hJidx[0]]>3 && Jet_id[hJidx[1]]>3";

     TCut cutQCD = "H_reg_pt>0 && met_pt>150 &&   MinIf$( abs(deltaPhi(met_phi,Jet_phi)), Jet_pt>30 && Jet_puId>3 &&  abs(Jet_eta)<4.5 )> 0.7 && Jet_btagCSV[hJidx[0]]>0.3 && Jet_btagCSV[hJidx[1]]>0.3 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30";


     TCut cutAntiQCD = "H_reg_pt>0 && met_pt>150   && Jet_btagCSV[hJidx[0]]>0.3 &&  Jet_btagCSV[hJidx[1]]>0.3 && min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30";




     TCut addCenJet30m0 = "(Sum$(Jet_pt>30 && Jet_puId>3 &&  abs(Jet_eta)<4.5)-2)>0";
     TCut addCenJet30e0 = "(Sum$(Jet_pt>30 && Jet_puId>3 && abs(Jet_eta)<4.5)-2)==0";
     TCut addCenJet30le1 = "(Sum$(Jet_pt>30 && Jet_puId>3 && abs(Jet_eta)<4.5)-2)<=1";
     TCut addBCenJet30e0 = "(Sum$(Jet_pt>30 && Jet_puId>3 && abs(Jet_eta)<2.5 && Jet_btagCSV>0.46)-2)==0";

     TCut addCenJet30e1 = "(Sum$(Jet_pt>30 && Jet_puId>3 && abs(Jet_eta)<4.5)-2)==1";
     TCut naddGoodLeptons10m0= "(Sum$(aLeptons_pt>10 && ( aLeptons_relIso03<0.4 && aLeptons_looseIdSusy!=0 ))+Sum$(vLeptons_pt>10 && (vLeptons_relIso03<0.4 && vLeptons_looseIdSusy!=0 )))>0";
     //     TCut naddGoodLeptons10e0= "(Sum$(aLeptons_pt>10 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 ||  aLeptons_jetDR>0.3 ))+Sum$(vLeptons_pt>10 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3 )))==0";
     TCut naddGoodLeptons10e0= "(Sum$(aLeptons_pt>10 && (aLeptons_relIso03<0.4 && aLeptons_looseIdSusy!=0 ))+Sum$(vLeptons_pt>10 && (vLeptons_relIso03<0.4 && vLeptons_looseIdSusy!=0  )))==0";
     TCut naddGoodTaus20e0= "Sum$(TauGood_idDecayMode >= 1 && TauGood_idCI3hit >= 1 && TauGood_pt > 20. && abs(TauGood_eta) < 2.3)==0";



     TCut cutSignalLoose = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  Jet_btagCSV[hJidx[0]]>0.46 &&  Jet_btagCSV[hJidx[1]]>0.46 && (H_reg_mass<100 || H_reg_mass>140)";

     TCut cutSignal0J = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  min(Jet_btagCSV[hJidx[0]], Jet_btagCSV[hJidx[1]])>0.46 && (H_reg_mass>100 && H_reg_mass<140)" + addCenJet30e0 + naddGoodLeptons10e0 + naddGoodTaus20e0 + addBCenJet30e0;
     TCut cutSignal1J = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  max(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]] )>0.80 && min(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]] )>0.46  && (H_reg_mass>100 && H_reg_mass<140)" + addCenJet30le1 + naddGoodLeptons10e0 + naddGoodTaus20e0 + addBCenJet30e0;
     TCut cutSignalSB = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  min(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])>0.46 && (H_reg_mass<100 || H_reg_mass>140)" + addCenJet30le1 + naddGoodLeptons10e0 + naddGoodTaus20e0 + addBCenJet30e0; // le1 -> e0

     TCut cutSignal1b = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  max(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])>0.97 && ((H_reg_mass>100 && H_reg_mass<140))" + addCenJet30e0 ;


     TCut cutTT = "(Vtype==2 || Vtype==3) && vLeptons_pt>30" + addCenJet30m0 + /*naddGoodLeptons5m0 + */" max(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])>0.46&& min(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])>0.46";
     TCut cutZlight = "Vtype==4" + addCenJet30e0;



     TCut cutZbb = "Vtype==4 && (H_reg_mass<100 || H_reg_mass>140)" + addCenJet30e0 +    "min(Jet_btagCSV[hJidx[1]],Jet_btagCSV[hJidx[0]])>0.46";

     //     TCut cutWlight = "(Vtype==2 || Vtype==3) && vLeptons_pt>30" + addCenJet30e0 + /*naddGoodLeptons5m0 +*/ "Jet_btagCSV[hJidx[0]]<0.46";
     TCut cutWlight =   addCenJet30e0 + "(Vtype==2 || Vtype==3) && vLeptons_pt>30" ; /*naddGoodLeptons5m0 +*/;

     TCut cutWbb =  addCenJet30e0 + "(Vtype==2 || Vtype==3) && vLeptons_pt>30" + /*naddGoodLeptons5m0 +*/ "min(Jet_btagCSV[hJidx[1]],Jet_btagCSV[hJidx[0]])>0.46"; 

//&& (nhttCandidates==0 || ( nhttCandidates==1 && (httCandidates_mass[0]<120 || httCandidates_mass[0]>190 )))";
  
     //TCut cutmc_all = "Vtype==4  &&  t[hJidx[0]]>30 ))  &&  Sum$(aLeptons_pt>5  && aLeptons_relIso03<1. && ( deltaR(Jet_eta[hJidx[0]], aLeptons_eta, Jet_phi[hJidx[0]], aLeptons_phi)>0.4 &&  deltaR(Jet_eta[hJidx[1]], aLeptons_eta, Jet_phi[hJidx[1]], aLeptons_phi)>0.4))>=0 && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId==1)<=99  && abs(deltaPhi(H_phi, V_phi))>2.5";
     /*
     Higgs candidate:             Vtype >=0
       MET>150:                       met_pt>150
       Higgs pT>150:                H_reg_pt>150
       CSV2>0.2:                      Jet_btagCSV[hJidx[1]]>0.2
       DeltaPhi(H,MET)>1.5      abs(TVector2::Phi_mpi_pi(H_reg_phi-met_phi))>1.5
       Higgs jet with pt>30       min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30
       H mass <300                 H_reg_mass<300 && H_reg_mass>0
     */

       /*


Nothing = 1 && <!Cuts|FlagsMET!>
NoQCD1 = <!Cuts|Nothing!> && H_reg_pt>150 && minDeltaPhiJet2Met>0.7
NoQCD2 = <!Cuts|NoQCD1!> && minDeltaPhiCenJetNoPU30Met >0.5
NoQCD3 = <!Cuts|NoQCD2!> && tkMetPVchs_pt>25 && tkMet_pt>25
       ;metType1p2_pt/sqrt(met_sumEt)>5 && DeltaPhiJet1Jet2<2.5



CR_Nothing = <!Cuts|NoQCD3!>
    CR_TTbar = <!Cuts|NoQCD3!> && naddGoodLeptons5>=1 && (addCenJet30)>=1 && Jet_btagCSV[hJidx[0]]>0.94
    CR_WLight = <!Cuts|NoQCD3!> && naddGoodLeptons5>=1 && (addCenJet30)<=1 && Jet_btagCSV[hJidx[0]]<0.94
CR_Wbb = <!Cuts|NoQCD3!> && naddGoodLeptons5>=1 && (addCenJet30)<=1 && Jet_btagCSV[hJidx[0]]>0.94									 
CR_ZLight = <!Cuts|NoQCD3!> && naddGoodLeptons5==0 && (addCenJet30)<=1 && Jet_btagCSV[hJidx[0]]<0.94	
CR_Zbb = <!Cuts|NoQCD3!> && naddGoodLeptons5==0 && (addCenJet30)<=1 && Jet_btagCSV[hJidx[0]]>0.94 && (H.mass<110 || H.mass>140)
CR_QCD = <!Cuts|Nothing!> && minDeltaPhiJet2Met<0.7

       */
     //     TCut cutdata = "json==1 && ( HLT_BIT_HLT_PFMET170_NoiseCleaned_v ||  HLT_BIT_HLT_PFMET90_PFMHT90_IDTight_v) && Flag_hbheFilterNew &&Flag_hbheIsoFilter && Flag_goodVertices &&Flag_eeBadScFilter &&Flag_CSCTightHaloFilter";
     TCut cutdata = "json==1 && ( HLT_BIT_HLT_PFMET170_NoiseCleaned_v ||  HLT_BIT_HLT_PFMET90_PFMHT90_IDTight_v) && Flag_hbheFilterNew && Flag_goodVertices &&Flag_eeBadScFilter";
     //     TCut trigMC = "(HLT_BIT_HLT_PFMET170_NoiseCleaned_v || HLT_BIT_HLT_PFMET120_PFMHT120_IDLoose_v)";
     //   TCut trigMC = "(HLT_BIT_HLT_PFMET170_NoiseCleaned_v ||  HLT_BIT_HLT_PFMET90_PFMHT90_IDLoose_v)";

      TCut trigMC = "(HLT_BIT_HLT_PFMET170_NoiseCleaned_v ||  HLT_BIT_HLT_PFMET90_PFMHT90_IDTight_v)";

    TString channel  = "monoH_Zp600";

    ev->read(   cutMinimal   ,  cutMinimal, "" , "tree"); 


    TCut weightmc = "efflumi * (2.32/3.) * sign(genWeight) * puWeight  * weightBtag(Jet_pt[hJidx[0]], Jet_eta[hJidx[0]] , abs(Jet_hadronFlavour[hJidx[0]])) * weightBtag(Jet_pt[hJidx[1]], Jet_eta[hJidx[1]], abs(Jet_hadronFlavour[hJidx[1]])) * weightTrig"  ;
    // TCut weightmc = "(efflumi * (1.26/3.)*  sign(genWeight))";
    TCut weightdata = "";

    /*
    TString var      = "H_reg_mass";                // the variable to plot
    TString title    = ";Hmass [GeV]";          // the title of the histogram
    TString plotname = "massBefore";      // the name of the image file
    int nbinsx       = 10;                      // number of bins
    double xlow      = 0;                    // the low edge of x-axis
    double xup       = 150;                   // the upper edge of x-axis
    */
    TString options  = "printStat:plotSig:plotData:!plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,

   
    TString varN[20] =     { "run==1 ? met_shifted_JetResUp_pt: met_pt" , "run==1? met_shifted_JetResDown_pt : met_pt" , "run==1?met_shifted_JetEnUp_pt:met_pt" , "run==1? met_shifted_JetEnDown_pt:met_pt",   "V_pt", "H_reg_mass",    "H_reg_pt"     ,"Sum$(Jet_pt>30 && Jet_puId>3 &&  abs(Jet_eta)<2.5)" , "min( abs( deltaPhi(met_phi,Jet_phi[hJidx[0]] )), abs(deltaPhi(met_phi,Jet_phi[hJidx[1]])) )", "htJet30" , "max(Jet_pt[hJidx[0]], Jet_pt[hJidx[1]])", "min(Jet_pt[hJidx[0]], Jet_pt[hJidx[1]])", "Jet_btagCSV[hJidx[0]]", "Jet_btagCSV[hJidx[1]]", "met_pt", "nPVs" , "MinIf$( abs( deltaPhi(met_phi,Jet_phi)), Jet_pt>30 && Jet_puId>3 &&  abs(Jet_eta)<2.5)" , "tkMet_pt" , "abs( deltaPhi(met_phi, tkMet_phi))", "abs(deltaPhi(met_phi, H_phi))"};
    TString titleN[20] =   {  ";metJerUp [GeV]" , ";metJerDown [GeV]" , ";metJecUp [GeV]", ";metJecDown [GeV]",  ";Vpt [GeV]" , ";Hmass [GeV]", ";H_reg_pt[GeV]",";NaddJets",  ";min#Delta#Phi(MET,hJets)" , ";HT [GeV]", ";H_jetPt1 [GeV]",  ";H_jetPt2",  ";btagCSV1" , ";btagCSV2", ";MET[GEV]", ";nPVs", ";min#Delta#Phi(MET,Jets30)" , ";tkMet [GeV]", ";#Delta#Phi(MET,tkMET)" , ";#Delta#Phi(MET,H)"};
    TString plotnameN[20] ={ "metJerUp", "metJerDown" , "metJecUp", "metJecDown" ,  "vpt" , "hmass",        "hpt"         ,"NaddJets",  "minDPhiMJ2" , "ht30", "hj1", "hj2", "csv1" , "csv2", "met" , "nPVs", "minDPhiMJ" , "tkMet", "DPhiMETtkMET", "DPhiMETH" };
    int nbinsxN[20] =      {60, 60 , 60 , 60,  50, 50,             50           , 8 ,    32,   50,  50, 50, 20, 20, 60, 30 , 32 , 30, 32, 32};
    double xlowN[20] =     { 0, 0 ,0 ,0 ,0,  0,              0            , 0,     0   , 0 ,  0 , 0 , 0 , 0 , 150 , 0 , 0, 0 , 0 , 0 };
    double xupN[20] =      { 900, 900, 900, 900, 800, 500,            800          , 8     ,3.2 , 800, 800, 800, 1 , 1, 900 , 30 , 3.2, 300, 3.2, 3.2};

    /*
    for (int i=0;  i<20; i++ ) {
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutSignalLoose ) * weightmc ,  cutdata + cutQCD + cutSignalLoose   , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "SignalLooseLin", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutSignal0J ) * weightmc ,  cutdata + cutQCD + cutSignal0J , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "Signal0JLin", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutSignal1J ) * weightmc ,  cutdata +  cutQCD + cutSignal1J  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "Signal1JLin", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutSignalSB ) * weightmc ,  cutdata + cutQCD + cutSignalSB  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "SignalSBLin", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutSignal1b ) * weightmc ,  cutdata  +  cutQCD + cutSignal1b  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "Signal1bLin", plotdir, options, 0, -1, 1 );
  

    MakePlots(ev, varN[i], (cutQCD + trigMC + cutTT ) * weightmc ,  cutdata + cutQCD + cutTT  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "TTbarLin", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutZlight ) * weightmc ,  cutdata + cutQCD + cutZlight  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "ZlightLin", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutZbb ) * weightmc ,  cutdata + cutQCD + cutZbb  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "ZbbLin", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutWlight ) * weightmc ,   cutdata + cutQCD + cutWlight  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "WlightLin", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutWbb ) * weightmc ,  cutdata + cutQCD + cutWbb  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "WbbLin", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutAntiQCD + trigMC ) * weightmc ,  cutdata + cutAntiQCD   , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "QCDLin", plotdir, options, 0, -1, 1 );

    }
    
 
    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,  

    for (int i=0;  i<20; i++ ) {
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutSignalLoose ) * weightmc , cutdata + cutQCD + cutSignalLoose  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "SignalLooseLog", plotdir, options, 0, -1, 1 );

    MakePlots(ev, varN[i], (cutQCD + trigMC + cutSignal0J ) * weightmc , cutdata + cutQCD + cutSignal0J   , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "Signal0JLog", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutSignal1J ) * weightmc , cutdata + cutQCD + cutSignal1J  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "Signal1JLog", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutSignalSB ) * weightmc , cutdata + cutQCD + cutSignalSB  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "SignalSBLog", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutSignal1b ) * weightmc , cutdata + cutQCD + cutSignal1b   , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "Signal1bLog", plotdir, options, 0, -1, 1 );

    MakePlots(ev, varN[i], (cutQCD + trigMC + cutTT ) * weightmc ,  cutdata + cutQCD + cutTT  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "TTbarLog", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutZlight ) * weightmc ,  cutdata + cutQCD + cutZlight  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "ZlightLog", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutZbb ) * weightmc ,  cutdata + cutQCD + cutZbb  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "ZbbLog", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutWlight ) * weightmc ,  cutdata + cutQCD + cutWlight  , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "WlightLog", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutQCD + trigMC + cutWbb ) * weightmc ,   cutdata + cutQCD + cutWbb , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "WbbLog", plotdir, options, 0, -1, 1 );
    MakePlots(ev, varN[i], (cutAntiQCD + trigMC ) * weightmc ,  cutdata + cutAntiQCD   , titleN[i], nbinsxN[i], xlowN[i], xupN[i], plotnameN[i] + "QCDLog", plotdir, options, 0, -1, 1 );

    }
    

    */

    //    TString var = "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi,   met_pt, 0 , met_phi)";
    TString var = "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)";


    TString title    = ";MT [GeV]";          // the title of the histogram
    TString plotname = "mt";      // the name of the image file
    int nbinsx       = 11;                      // number of bins
    double xlow      = 0;                    // the low edge of x-axis
    double xup       = 1100;                   // the upper edge of x-axis
    bool doRebin = false;
    bool doOverFlow = false;


    int metedge[4] =      { 200, 300, 400, 9999 };
    for (int i =0; i<3; i++) {
    TCut cutmet = Form("met_pt>%d && met_pt<%d",  metedge[i],metedge[i+1]);


    options  = "printStat:plotSig:plotData:!plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
    //  MakePlots(ev, var, (cutQCD + trigMC + cutSignalLoose + "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)>0") * weightmc ,  cutdata + cutQCD + cutSignalLoose + "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)>0" , title, nbinsx, xlow, xup, plotname + "SignalLoose", plotdir, options,     doRebin, 9);
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet + "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)>0") * weightmc ,  cutdata + cutQCD + cutSignal0J + cutmet + "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)>0" , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow);

    TString dcname    = Form("zp600met%d%d_mt_0J.txt", metedge[i],metedge[i+1]) ;   // the datacard name
    TString wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J.root";                // the workspace name
    bool useshapes = false;
    TString options1  = "unblind:SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);


    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet + "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)>0") * weightmc ,  cutdata + cutQCD + cutSignal1J + cutmet+ "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)>0"  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow);

    dcname    = Form("zp600met%d%d_mt_1J.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) +"Signal1J.root";                // the workspace name
    useshapes = false;
    options1  = "unblind:SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);


    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet + "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)>0" ) * weightmc ,  cutdata + cutQCD + cutSignalSB + cutmet + "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)>0"  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow);

    dcname    =  Form("zp600met%d%d_mt_SB.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB.root";                // the workspace name
    useshapes = false;
    options1  = "unblind:SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);


    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet + "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)>0" ) * weightmc ,  cutdata + cutQCD + cutSignal1b + cutmet + "evalHMETMassiveMt(H_reg_mass, H_reg_pt, H_reg_eta, H_phi, 91.2,  met_pt, 0 , met_phi)>0" , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow);

    dcname    =  Form("zp600met%d%d_mt_1b.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1b.root";                // the workspace name
    useshapes = false;
    options1  = "unblind:SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);

    //    nbinsx       = 10;                      // number of bins
    // xlow      = 0;                    // the low edge of x-axis
    // xup      = 600;                    // the low edge of x-axis
    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc ,  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow);




    dcname    = Form("zp600met%d%d_mt_TT.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar.root";                // the workspace name
    useshapes = false;
    options1  = "unblind:SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);






    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc ,  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options);

    dcname    = Form("zp600met%d%d_mt_Wlight.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight.root";                // the workspace name
    useshapes = false;


    options1  = "unblind:SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc ,   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow);
    //    MakePlots(ev, var, (cutAntiQCD + trigMC + cutWbb ) * weightmc ,  cutdata + cutAntiQCD + cutWbb  , title, nbinsx, xlow, xup, plotname + "QCD", plotdir, options);
    dcname    = Form("zp600met%d%d_mt_Wbb.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) +"Wbb.root";                // the workspace name
    useshapes = false;


    options1  = "unblind:SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);




    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc ,   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow);

    //    MakePlots(ev, var, (cutQCD + trigMC + cutTT ) * weightmc ,  cutdata + cutQCD + cutTT  , title, nbinsx, xlow, xup, plotname + "TTbar", plotdir, options);





    dcname    = Form("zp600met%d%d_mt_Zlight.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Zlight.root";                // the workspace name
    useshapes = false;
    options1  = "unblind:SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);



    }
    


}

//  LocalWords:  deltaPhi
