#include "plotHistos_boostedZH_13TeV_monoH_shape_all_syst.h"
#include "TROOT.h"
#include "TStyle.h"
#include "puWeight.h"

#include "BTagCalibrationStandalone.h"
// setup calibration readers
BTagCalibration calib("CSVv2", "CSVv2.csv");
BTagCalibrationReader readerHF(&calib,               // calibration instance
                             BTagEntry::OP_MEDIUM,  // operating point
                             "mujets",               // measurement type
                             "central");           // systematics type

BTagCalibrationReader readerLF(&calib,               // calibration instance
                             BTagEntry::OP_MEDIUM,  // operating point
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

void plotHistos_boostedZH_13TeV_monoH_shape_all_syst() {
    gROOT->LoadMacro("tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    TH1::SetDefaultSumw2(1);
    //gROOT->SetBatch(1);

    TString plotdir = "plotsmonoH_SF_shapesMET_CMS/";

    if (gSystem->AccessPathName(plotdir))
        gSystem->mkdir(plotdir);


    ////////////////////////////////////////////////////////////////////////////
    // Task 1 (a)                                                             //
    // - please comment out other tasks                                       //
    ////////////////////////////////////////////////////////////////////////////

    // Read from ntuples
    Events * ev = new Events();




     // TCut cutmc_all = "MinIf$(abs(deltaPhi(V_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>1.5";

     TCut cutMinimal = "mhtJet30>150 && HCSV_pt>100 && Jet_btagCSV[hJCidx[1]]>0.4 && HCSV_mass<500 && abs(deltaPhi(HCSV_phi,met_phi))>0.7 && json && (Vtype==2 || Vtype==3 || Vtype==4) &&  H_reg_pt>0 && met_pt>150  && min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 && Jet_btagCSV[hJidx[1]]>0.46 && Jet_btagCSV[hJidx[0]]>0.46 &&  Jet_id[hJidx[0]]>3 && Jet_id[hJidx[1]]>3";

     TCut cutQCD = "H_reg_pt>0 && met_pt>170 &&   MinIf$( abs(deltaPhi(met_phi,Jet_phi)), Jet_pt>30 && Jet_puId>3 &&  abs(Jet_eta)<4.5 )> 0.7 && abs( deltaPhi(met_phi, tkMet_phi))<0.5 && Jet_btagCSV[hJidx[0]]>0.3 && Jet_btagCSV[hJidx[1]]>0.3 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30";


     TCut cutAntiQCD = "H_reg_pt>0 && met_pt>150   && Jet_btagCSV[hJidx[0]]>0.3 &&  Jet_btagCSV[hJidx[1]]>0.3 && min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30";




     TCut addCenJet30m0 = "(Sum$(Jet_pt>30 && Jet_puId>3 &&  abs(Jet_eta)<4.5)-2)>0";
     TCut addCenJet30e0 =  "(Sum$(Jet_pt>30 && Jet_puId>3 && abs(Jet_eta)<4.5)-2)==0";
     TCut addCenJet30le1 = "(Sum$(Jet_pt>30 && Jet_puId>3 && abs(Jet_eta)<4.5)-2)<=1";
     TCut addBCenJet30e0 = "(Sum$(Jet_pt>30 && Jet_puId>3 && abs(Jet_eta)<2.5 && Jet_btagCSV>0.46)-2)==0";

     TCut addCenJet30e1 = "(Sum$(Jet_pt>30 && Jet_puId>3 && abs(Jet_eta)<4.5)-2)==1";
     //     TCut naddGoodLeptons10m0= "(Sum$(aLeptons_pt>10 && ( aLeptons_relIso03<0.4 && aLeptons_looseIdSusy!=0 ))+Sum$(vLeptons_pt>10 && (vLeptons_relIso03<0.4 && vLeptons_looseIdSusy!=0 )))>0";
     TCut naddGoodLeptons10m0= "(Sum$(aLeptons_pt>7 && ( aLeptons_relIso03<0.4 && aLeptons_looseIdSusy!=0 ))+Sum$(vLeptons_pt>7 && (vLeptons_relIso03<0.4 && vLeptons_looseIdSusy!=0 )))>0";
       TCut naddGoodLeptons10e0= "(Sum$(aLeptons_pt>7 && (aLeptons_jetBTagCSV<0.25 || aLeptons_relIso03<0.4 || aLeptons_looseIdSusy!=0 ||  aLeptons_jetDR>0.3 ))+Sum$(vLeptons_pt>7 && (vLeptons_jetBTagCSV<0.25 || vLeptons_relIso03<0.4 || vLeptons_looseIdSusy!=0 || vLeptons_jetDR>0.3 )))==0";
	       //     TCut naddGoodLeptons10e0= "(Sum$(aLeptons_pt>10 && (aLeptons_relIso03<0.4 && aLeptons_looseIdSusy!=0 ))+Sum$(vLeptons_pt>10 && (vLeptons_relIso03<0.4 && vLeptons_looseIdSusy!=0  )))==0";
     //     TCut naddGoodLeptons10e0= "(Sum$(aLeptons_pt>7 && (aLeptons_relIso03<0.4 ))+Sum$(vLeptons_pt>7 && (vLeptons_relIso03<0.4  )))==0";
     TCut naddGoodTaus20e0= "Sum$(TauGood_idDecayMode >= 1 && TauGood_idCI3hit >= 1 && TauGood_pt > 20. && abs(TauGood_eta) < 2.3)==0";



     TCut cutSignalLoose = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  Jet_btagCSV[hJidx[0]]>0.46 &&  Jet_btagCSV[hJidx[1]]>0.46 && (HaddJetsdR08_mass<100 || HaddJetsdR08_mass>140)";

     TCut cutSignal0J = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  min(Jet_btagCSV[hJidx[0]], Jet_btagCSV[hJidx[1]])>0.46 && deltaR_jj<2.0 && (HaddJetsdR08_mass>100 && HaddJetsdR08_mass<140)" + addCenJet30e0 + naddGoodLeptons10e0 + naddGoodTaus20e0 + addBCenJet30e0;
     //     TCut cutSignal1J = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  max(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]] )>0.935 && min(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]] )>0.46" + addCenJet30le1 + naddGoodLeptons10e0 + naddGoodTaus20e0 + addBCenJet30e0;
     TCut cutSignal1J = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  max(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]] )>0.8 && min(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]] )>0.8 && deltaR_jj<2.0 &&  (HaddJetsdR08_mass>100 && HaddJetsdR08_mass<140)" + addCenJet30le1 + naddGoodLeptons10e0 + naddGoodTaus20e0 + addBCenJet30e0;


     TCut cutSignalSB = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  min(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])>0.8 && (HaddJetsdR08_mass<100 || HaddJetsdR08_mass>140)" + addCenJet30le1 + naddGoodLeptons10e0 + naddGoodTaus20e0 + addBCenJet30e0; // le1 -> e0

     //     TCut cutSignal1b = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  max(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])<0.935 && min(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])>0.46" + addCenJet30e0 ;
     TCut cutSignal1b = "Vtype==4 &&  min(Jet_pt[hJidx[0]],Jet_pt[hJidx[1]])>30 &&  max(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])<0.935 && min(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])>0.46" + addCenJet30e0 ;


     TCut cutTT = "(Vtype==2 || Vtype==3) && vLeptons_pt>30" + addCenJet30m0 + /*naddGoodLeptons5m0 + */" max(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])>0.46&& min(Jet_btagCSV[hJidx[0]],Jet_btagCSV[hJidx[1]])>0.46";
     TCut cutZlight = "Vtype==4" + addCenJet30e0;



     TCut cutZbb = "Vtype==4 && (HaddJetsdR08_mass<100 || HaddJetsdR08_mass>140)" + addCenJet30e0 +    "min(Jet_btagCSV[hJidx[1]],Jet_btagCSV[hJidx[0]])>0.46";

     //     TCut cutWlight = "(Vtype==2 || Vtype==3) && vLeptons_pt>30" + addCenJet30e0 + /*naddGoodLeptons5m0 +*/ "Jet_btagCSV[hJidx[0]]<0.46";
     TCut cutWlight =   addCenJet30e0 + "(Vtype==2 || Vtype==3) && vLeptons_pt>30" ; /*naddGoodLeptons5m0 +*/;

     TCut cutWbb =  addCenJet30e0 + "(Vtype==2 || Vtype==3) && vLeptons_pt>30" + /*naddGoodLeptons5m0 +*/ "min(Jet_btagCSV[hJidx[1]],Jet_btagCSV[hJidx[0]])>0.46"; 


     TCut cutdata = "json==1 && ( HLT_BIT_HLT_PFMET170_NoiseCleaned_v ||  HLT_BIT_HLT_PFMET90_PFMHT90_IDTight_v) && Flag_hbheFilterNew && Flag_goodVertices &&Flag_eeBadScFilter";

      TCut trigMC = "(HLT_BIT_HLT_PFMET170_NoiseCleaned_v ||  HLT_BIT_HLT_PFMET90_PFMHT90_IDTight_v)";

    TString channel  = "monoH_Zp600";

    ev->read(   cutMinimal   ,  cutMinimal, "" , "tree"); 


    TCut weightmc = "efflumi * (2.32/3.) * sign(genWeight) * puWeight  * weightBtag(Jet_pt[hJidx[0]], Jet_eta[hJidx[0]] , abs(Jet_hadronFlavour[hJidx[0]])) * weightBtag(Jet_pt[hJidx[1]], Jet_eta[hJidx[1]], abs(Jet_hadronFlavour[hJidx[1]])) * weightTrig"  ;
    // TCut weightmc = "(efflumi * (1.26/3.)*  sign(genWeight))";
    TCut weightdata = "";

    TString options  = "printStat:plotSig:plotData:!plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,


    TString var = "met_pt";


    TString title    = ";MET [GeV]";          // the title of the histogram
    TString plotname = "MET";      // the name of the image file
    int nbinsx       = 2;                      // number of bins
    double xlow      = 150;                    // the low edge of x-axis
    double xup       = 350;                   // the upper edge of x-axis
    bool doRebin = false;
    bool doOverFlow = true;


    int metedge[2] =      { 150,  9999 };
    for (int i =0; i<1; i++) {



      // nominal  

  TCut cutmet = Form("met_pt>%d && met_pt<%d",  metedge[i],metedge[i+1]);


    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   


    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc ,  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)");

    TString dcname    = Form("zp600metshape%d%d_0J.txt", metedge[i],metedge[i+1]) ;   // the datacard name
    TString wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J.root";                // the workspace name
    bool useshapes = true;
    TString options1  = "unblind:!SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);


    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc ,  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)");

    dcname    = Form("zp600metshape%d%d_1J.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) +"Signal1J.root";                // the workspace name
    useshapes = true;
    options1  = "unblind:!SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);


    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)");

    dcname    =  Form("zp600metshape%d%d_SB.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB.root";                // the workspace name
    useshapes = true;
    options1  = "unblind:!SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);


    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)");

    dcname    =  Form("zp600metshape%d%d_1b.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1b.root";                // the workspace name
    useshapes = true;
    options1  = "unblind:!SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);

    //    nbinsx       = 10;                      // number of bins
    // xlow      = 0;                    // the low edge of x-axis
    // xup      = 1000;                    // the low edge of x-axis
    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc ,  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)");




    dcname    = Form("zp600metshape%d%d_TT.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar.root";                // the workspace name
    useshapes = true;
    options1  = "unblind:!SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);






    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc ,  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, "W+light CR (resolved)");

    dcname    = Form("zp600metshape%d%d_Wlight.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight.root";                // the workspace name
    useshapes = true;


    options1  = "unblind:!SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc ,   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)");
    //    MakePlots(ev, var, (cutAntiQCD + trigMC + cutWbb ) * weightmc ,  cutdata + cutAntiQCD + cutWbb  , title, nbinsx, xlow, xup, plotname + "QCD", plotdir, options);
    dcname    = Form("zp600metshape%d%d_Wbb.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) +"Wbb.root";                // the workspace name
    useshapes = true;


    options1  = "unblind:!SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);




    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc ,   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)");

    //    MakePlots(ev, var, (cutQCD + trigMC + cutTT ) * weightmc ,  cutdata + cutQCD + cutTT  , title, nbinsx, xlow, xup, plotname + "TTbar", plotdir, options);





    dcname    = Form("zp600metshape%d%d_Zlight.txt", metedge[i],metedge[i+1]);   // the datacard name
    wsname    = plotdir + plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Zlight.root";                // the workspace name
    useshapes = true;
    options1  = "unblind:!SplusB";
    MakeDatacard(channel, dcname, wsname, useshapes, options1);




    // JEC up
      // nominal  
    //    "run==1 ? met_shifted_JetResUp_pt: met_pt" , "run==1? met_shifted_JetResDown_pt : met_pt" , "run==1?met_shifted_JetEnUp_pt:met_pt" , "run==1? met_shifted_JetEnDown_pt:met_pt"

     var = "run==1?met_shifted_JetEnUp_pt:met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   
    //CMS_monoHbb_res_j                    lnN    1.03   1.03    -     -     -    1.03  1.03  
    //CMS_monoHbb_scale_j   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc ,  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_scale_jUp");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc ,  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_scale_jUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_scale_jUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_scale_jUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc ,  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_scale_jUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc ,  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_scale_jUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc ,   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_scale_jUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc ,   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_scale_jUp");



    // JEC down
      // nominal  
    //    "run==1 ? met_shifted_JetResUp_pt: met_pt" , "run==1? met_shifted_JetResDown_pt : met_pt" , "run==1?met_shifted_JetEnUp_pt:met_pt" , "run==1? met_shifted_JetEnDown_pt:met_pt"

     var = "run==1?met_shifted_JetEnDown_pt:met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc ,  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_scale_jDown");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc ,  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_scale_jDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_scale_jDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_scale_jDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc ,  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_scale_jDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc ,  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_scale_jDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc ,   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_scale_jDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc ,   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_scale_jDown");




    // JEC down
      // nominal  
    //    "run==1 ? met_shifted_JetResUp_pt: met_pt" , "run==1? met_shifted_JetResDown_pt : met_pt" , "run==1?met_shifted_JetEnUp_pt:met_pt" , "run==1? met_shifted_JetEnDown_pt:met_pt"

     var = "run==1?met_shifted_JetResUp_pt:met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc ,  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_res_jUp");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc ,  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_res_jUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_res_jUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_res_jUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc ,  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_res_jUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc ,  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_res_jUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc ,   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_res_jUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc ,   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_res_jUp");



    // JEC down
      // nominal  
    //    "run==1 ? met_shifted_JetResUp_pt: met_pt" , "run==1? met_shifted_JetResDown_pt : met_pt" , "run==1?met_shifted_JetEnUp_pt:met_pt" , "run==1? met_shifted_JetEnDown_pt:met_pt"

     var = "run==1?met_shifted_JetResDown_pt:met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc ,  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_res_jDown");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc ,  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_res_jDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_res_jDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_res_jDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc ,  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_res_jDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc ,  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_res_jDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc ,   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_res_jDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc ,   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_res_jDown");


    //CMS_monoHbb_metUnclusteredEn
     var = "run==1?met_shifted_UnclusteredEnUp_pt:met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc ,  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_metUnclusteredEnUp");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc ,  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_metUnclusteredEnUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_metUnclusteredEnUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_metUnclusteredEnUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc ,  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_metUnclusteredEnUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc ,  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_metUnclusteredEnUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc ,   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_metUnclusteredEnUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc ,   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_metUnclusteredEnUp");




     var = "run==1?met_shifted_UnclusteredEnDown_pt:met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc ,  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_metUnclusteredEnDown");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc ,  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_metUnclusteredEnDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_metUnclusteredEnDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_metUnclusteredEnDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc ,  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_metUnclusteredEnDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc ,  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_metUnclusteredEnDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc ,   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_metUnclusteredEnDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc ,   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_metUnclusteredEnDown");


    
    //     var = "run==1?met_pt*LHE_weights_scale_wgt[0]:met_pt";
     var = "met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc * "LHE_weights_scale_wgt[0]" ,  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_scaleFUp");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc * "LHE_weights_scale_wgt[0]",  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_scaleFUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc  * "LHE_weights_scale_wgt[0]",  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_scaleFUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc * "LHE_weights_scale_wgt[0]",  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_scaleFUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc * "LHE_weights_scale_wgt[0]",  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_scaleFUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc * "LHE_weights_scale_wgt[0]",  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_scaleFUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc * "LHE_weights_scale_wgt[0]",   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_scaleFUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc * "LHE_weights_scale_wgt[0]" ,   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_scaleFUp");



    

     var = "met_pt";

    

     options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
     
    
     MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc * "LHE_weights_scale_wgt[1]",  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_scaleFDown");
    
   
     MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc * "LHE_weights_scale_wgt[1]",  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_scaleFDown");

    
    
    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc * "LHE_weights_scale_wgt[1]",  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_scaleFDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc * "LHE_weights_scale_wgt[1]",  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_scaleFDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc * "LHE_weights_scale_wgt[1]",  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_scaleFDown");

    

    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc * "LHE_weights_scale_wgt[1]",  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_scaleFDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc * "LHE_weights_scale_wgt[1]",   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_scaleFDown");

    
    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc * "LHE_weights_scale_wgt[1]",   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_scaleFDown");





    



     var = "met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc * "LHE_weights_scale_wgt[2]",  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_scaleRUp");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc  * "LHE_weights_scale_wgt[2]",  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_scaleRUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc  * "LHE_weights_scale_wgt[2]",  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_scaleRUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc  * "LHE_weights_scale_wgt[2]",  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_scaleRUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc  * "LHE_weights_scale_wgt[2]",  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_scaleRUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc  * "LHE_weights_scale_wgt[2]",  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_scaleRUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc  * "LHE_weights_scale_wgt[2]",   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_scaleRUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc  * "LHE_weights_scale_wgt[2]",   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_scaleRUp");





     var = "met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc  * "LHE_weights_scale_wgt[3]",  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_scaleRDown");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc  * "LHE_weights_scale_wgt[3]",  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_scaleRDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc  * "LHE_weights_scale_wgt[3]",  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_scaleRDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc  * "LHE_weights_scale_wgt[3]",  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_scaleRDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc  * "LHE_weights_scale_wgt[3]",  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_scaleRDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc  * "LHE_weights_scale_wgt[3]",  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_scaleRDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc  * "LHE_weights_scale_wgt[3]",   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_scaleRDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc  * "LHE_weights_scale_wgt[3]",   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_scaleRDown");




    //    ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId) 

    var = "met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc * "1/(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))" ,  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_ewkDown");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc * "1/(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_ewkDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc * "1/(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_ewkDown");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc * "1/(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_ewkDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc * "1/(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_ewkDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc * "1/(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_ewkDown");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc * "1/(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_ewkDown");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc * "1/(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_ewkDown");



    //    ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId) 

    var = "met_pt";

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc * "(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", "_CMS_monoHbb_ewkUp");

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc * "(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", "_CMS_monoHbb_ewkUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc * "(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", "_CMS_monoHbb_ewkUp");

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc * "(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", "_CMS_monoHbb_ewkUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc * "(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", "_CMS_monoHbb_ewkUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc * "(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", "_CMS_monoHbb_ewkUp");



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc * "(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", "_CMS_monoHbb_ewkUp");


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc * "(ptWeightZllH(nGenVbosons,GenVbosons_pt,VtypeSim,GenVbosons_pdgId))",   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", "_CMS_monoHbb_ewkUp");

    
    // pdf now
    /// 
    /*
    for (int j=9; j<110; j++){

      var = Form("run==1?met_pt*LHE_weights_pdf_wgt[%d]:met_pt", j);

  

    options  = "printStat:plotSig:plotData:plotLog:!plotNorm:plotBoostedZH";    // use "!plotLog" to plot on log-y scale,
   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal0J + cutmet )  *   weightmc ,  cutdata + cutQCD + cutSignal0J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal0J", plotdir, options,     doRebin, 6, doOverFlow, "SR 0J (resolved)", Form("_CMS_monoHbb_pdf%dUp",j));

   
    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1J + cutmet ) * weightmc ,  cutdata + cutQCD + cutSignal1J + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Signal1J", plotdir, options,  doRebin, 6, doOverFlow, "SR (resolved)", Form("_CMS_monoHbb_pdf%dUp",j));

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignalSB + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignalSB + cutmet   , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "SignalSB", plotdir, options,  doRebin, 6, doOverFlow, "mass sideband CR (resolved)", Form("_CMS_monoHbb_pdf%dUp",j));

   

    MakePlots(ev, var, (cutQCD + trigMC + cutSignal1b + cutmet  ) * weightmc ,  cutdata + cutQCD + cutSignal1b + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Signal1b", plotdir, options,  doRebin, 6, doOverFlow, "SR 1b (resolved)", Form("_CMS_monoHbb_pdf%dUp",j));


    MakePlots(ev, var, (cutQCD + trigMC + cutTT + cutmet ) * weightmc ,  cutdata + cutQCD + cutTT  + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "TTbar", plotdir, options,  doRebin, 6, doOverFlow, "TT CR (resolved)", Form("_CMS_monoHbb_pdf%dUp",j));



    MakePlots(ev, var, (cutQCD + trigMC + cutWlight + cutmet ) * weightmc ,  cutdata + cutQCD + cutWlight + cutmet  , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) + "Wlight", plotdir, options, doRebin, 6, doOverFlow, "W+light CR (resolved)", Form("_CMS_monoHbb_pdf%dUp",j));



    MakePlots(ev, var, (cutQCD + trigMC + cutWbb + cutmet) * weightmc ,   cutdata + cutQCD + cutWbb + cutmet , title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Wbb", plotdir, options,  doRebin, 6, doOverFlow, "j CR (resolved)", Form("_CMS_monoHbb_pdf%dUp",j));


    MakePlots(ev, var, (cutQCD + trigMC + cutZlight + cutmet ) * weightmc ,   cutdata + cutQCD + cutZlight + cutmet, title, nbinsx, xlow, xup, plotname + Form("%d%d", metedge[i],metedge[i+1]) +  "Zlight", plotdir, options,  doRebin, 6, doOverFlow, "Z+light CR (resolved)", Form("_CMS_monoHbb_pdf%dUp",j));


    }

    */
    


    }


    }

//  LocalWords:  deltaPhi
