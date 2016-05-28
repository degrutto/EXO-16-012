#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TString.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TMath.h"
//#include "TPRegexp.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeFormula.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

#include "HelperTMVA.h"
//#include "HelperFunctions.h"
#include "HelperNtuples.h"
#include "HelperVHbbDataFormats.h"

#include "EquationSolver.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/QuantFuncMathCore.h"

namespace math {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > XYZTLorentzVector;
}

#define MAXJ 100
#define MAXL 110

//#define CSVSYST
#ifdef CSVSYST
#include "BTagReshaping.h"  // FIXME: local version is out of date with V6 Step 2
#endif

//#define JECFWLITE
#ifdef JECFWLITE
#include "JECFWLiteStandalone.h"
#endif

//#define STITCH

// FIXME: also patch V.pt

//#define PATCHMETTYPE1CORR


// initial
float weightTrig(float met)
{
  if (met<150.) return 0.;
		  //		  "exponential", "[2]*(1e0-exp(-[0]*(x-[1])))",minx,500);
  /* 1  p0           4.37294e-02   2.73746e-02  -2.31457e-03  -8.39905e-04
   2  p1           1.09515e+02   2.90232e+01   4.17730e-04   2.71771e-04
   3  p2           9.39967e-01   2.44644e-02   2.44644e-02   9.50151e-05
  */
  if (met>=150 && met<500.) return  (9.99e-01 * (1e0-exp(-4.37294e-02*(met-1.09515e+02))));
  if (met>=500.) return 9.99e-01;
  return 0.;
}



// FIXME: change it to output the delta, not met+delta
double shift_met_by_Nvtx(double met, double metphi, int Nvtx, int EVENT_run)
{
    // from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py?revision=1.6&view=markup
    double mex = met * cos(metphi);
    double mey = met * sin(metphi);
    double px = 0.0, py = 0.0;
    if (EVENT_run != 1) {
        //pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data
        px = +0.2661 + 0.3217*Nvtx;
        py = -0.2251 - 0.1747*Nvtx;
    } else {
        //pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc
        px = +0.1166 + 0.0200*Nvtx;
        py = +0.2764 - 0.1280*Nvtx;
    }
    mex -= px;
    mey -= py;
    return TMath::Sqrt(mex * mex + mey * mey);
}

// FIXME: change it to output the delta, not met+delta
double shift_metphi_by_Nvtx(double met, double metphi, int Nvtx, int EVENT_run)
{
    // from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py?revision=1.6&view=markup
    double mex = met * cos(metphi);
    double mey = met * sin(metphi);
    double px = 0.0, py = 0.0;
    if (EVENT_run != 1) {
        //pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data
        px = +0.2661 + 0.3217 * Nvtx;
        py = -0.2251 - 0.1747 * Nvtx;
    } else {
        //pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc
        px = +0.1166 + 0.0200 * Nvtx;
        py = +0.2764 - 0.1280 * Nvtx;
    }
    mex -= px;
    mey -= py;
    if (mex == 0.0 && mey == 0.0)
        return 0.0;
    double phi1 = TMath::ATan2(mey, mex);
    double phi2 = TMath::ATan2(mey, mex) - 2.0 * M_PI;
    if (TMath::Abs(phi1 - metphi) < TMath::Abs(phi2 - metphi) + 0.5 * M_PI)
        return phi1;
    else
        return phi2;
}

TLorentzVector makePtEtaPhiE(double ptCor, double pt, double eta, double phi, double e, 
                             bool rescaleEnergy=true)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(ptCor, eta, phi, (rescaleEnergy ? (e * ptCor / pt) : e));
    return j;
}


TLorentzVector makePtEtaPhiM(double ptCor, double pt, double eta, double phi, double m, 
                             bool print=false)
{
    TLorentzVector j;
    if (print ) { 
      std::cout << "evaluating vector with ptCor, pt, eta, phi, m " << ptCor << " , " <<  pt << " , " << eta << " , " << phi << " , " << m << std::endl;
    }
    j.SetPtEtaPhiM(ptCor, eta, phi,  m );
    return j;
}


// Copied from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/TopQuarkAnalysis/SingleTop/src/TopProducer.cc?revision=1.9&view=markup
TLorentzVector getNu4Momentum(const TLorentzVector& TLepton, const TLorentzVector& TMET)
{
  const math::XYZTLorentzVector Lepton(TLepton.Px(), TLepton.Py(), TLepton.Pz(), TLepton.E());
  const math::XYZTLorentzVector MET(TMET.Px(), TMET.Py(), 0., TMET.E());

  double  mW = 80.38;

  //std::vector<math::XYZTLorentzVector> result;
  std::vector<TLorentzVector> result;

  //  double Wmt = sqrt(pow(Lepton.et()+MET.pt(),2) - pow(Lepton.px()+MET.px(),2) - pow(Lepton.py()+MET.py(),2) );
    
  double MisET2 = (MET.px()*MET.px() + MET.py()*MET.py());
  double mu = (mW*mW)/2 + MET.px()*Lepton.px() + MET.py()*Lepton.py();
  double a  = (mu*Lepton.pz())/(Lepton.energy()*Lepton.energy() - Lepton.pz()*Lepton.pz());
  double a2 = TMath::Power(a,2);
  double b  = (TMath::Power(Lepton.energy(),2.)*(MisET2) - TMath::Power(mu,2.))/(TMath::Power(Lepton.energy(),2) - TMath::Power(Lepton.pz(),2));
  double pz1(0),pz2(0),pznu(0);
  int nNuSol(0);

  //math::XYZTLorentzVector p4nu_rec;
  TLorentzVector p4nu_rec;
  math::XYZTLorentzVector p4W_rec;
  math::XYZTLorentzVector p4b_rec;
  math::XYZTLorentzVector p4Top_rec;
  math::XYZTLorentzVector p4lep_rec;

  p4lep_rec.SetPxPyPzE(Lepton.px(),Lepton.py(),Lepton.pz(),Lepton.energy());
  
  //math::XYZTLorentzVector p40_rec(0,0,0,0);

  if(a2-b > 0 ){
    //if(!usePositiveDeltaSolutions_)
    //  {
    //    result.push_back(p40_rec);
    //    return result;
    //  }
    double root = sqrt(a2-b);
    pz1 = a + root;
    pz2 = a - root;
    nNuSol = 2;

    //if(usePzPlusSolutions_)pznu = pz1;    
    //if(usePzMinusSolutions_)pznu = pz2;
    //if(usePzAbsValMinimumSolutions_){
      pznu = pz1;
      if(fabs(pz1)>fabs(pz2)) pznu = pz2;
    //}

    double Enu = sqrt(MisET2 + pznu*pznu);

    p4nu_rec.SetPxPyPzE(MET.px(), MET.py(), pznu, Enu);

    result.push_back(p4nu_rec);

  }
  else{
    //if(!useNegativeDeltaSolutions_){
    //  result.push_back(p40_rec);
    //  return result;
    //}
    //    double xprime = sqrt(mW;

    double ptlep = Lepton.pt(),pxlep=Lepton.px(),pylep=Lepton.py(),metpx=MET.px(),metpy=MET.py();

    double EquationA = 1;
    double EquationB = -3*pylep*mW/(ptlep);
    double EquationC = mW*mW*(2*pylep*pylep)/(ptlep*ptlep)+mW*mW-4*pxlep*pxlep*pxlep*metpx/(ptlep*ptlep)-4*pxlep*pxlep*pylep*metpy/(ptlep*ptlep);
    double EquationD = 4*pxlep*pxlep*mW*metpy/(ptlep)-pylep*mW*mW*mW/ptlep;

    std::vector<long double> solutions = EquationSolve<long double>((long double)EquationA,(long double)EquationB,(long double)EquationC,(long double)EquationD);

    std::vector<long double> solutions2 = EquationSolve<long double>((long double)EquationA,-(long double)EquationB,(long double)EquationC,-(long double)EquationD);

    double deltaMin = 14000*14000;
    double zeroValue = -mW*mW/(4*pxlep); 
    double minPx=0;
    double minPy=0;

    //    std::cout<<"a "<<EquationA << " b " << EquationB  <<" c "<< EquationC <<" d "<< EquationD << std::endl; 
      
    //if(usePxMinusSolutions_){
      for( int i =0; i< (int)solutions.size();++i){
      if(solutions[i]<0 ) continue;
      double p_x = (solutions[i]*solutions[i]-mW*mW)/(4*pxlep); 
      double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x -mW*ptlep*solutions[i])/(2*pxlep*pxlep);
      double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy); 

      //      std::cout<<"intermediate solution1 met x "<<metpx << " min px " << p_x  <<" met y "<<metpy <<" min py "<< p_y << std::endl; 

      if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
      minPx=p_x;
      minPy=p_y;}
      //     std::cout<<"solution1 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
      }
    //}

    //if(usePxPlusSolutions_){
      for( int i =0; i< (int)solutions2.size();++i){
        if(solutions2[i]<0 ) continue;
        double p_x = (solutions2[i]*solutions2[i]-mW*mW)/(4*pxlep); 
        double p_y = ( mW*mW*pylep + 2*pxlep*pylep*p_x +mW*ptlep*solutions2[i])/(2*pxlep*pxlep);
        double Delta2 = (p_x-metpx)*(p_x-metpx)+(p_y-metpy)*(p_y-metpy); 
        //  std::cout<<"intermediate solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
        if(Delta2< deltaMin && Delta2 > 0){deltaMin = Delta2;
          minPx=p_x;
          minPy=p_y;
        }
        //        std::cout<<"solution2 met x "<<metpx << " min px " << minPx  <<" met y "<<metpy <<" min py "<< minPy << std::endl; 
      }
    //}

    double pyZeroValue= ( mW*mW*pxlep + 2*pxlep*pylep*zeroValue);
    double delta2ZeroValue= (zeroValue-metpx)*(zeroValue-metpx) + (pyZeroValue-metpy)*(pyZeroValue-metpy);

    if(deltaMin==14000*14000) return TLorentzVector(0,0,0,0);
    //if(deltaMin==14000*14000) return result.front();
    //    else std::cout << " test " << std::endl;

    if(delta2ZeroValue < deltaMin){
      deltaMin = delta2ZeroValue;
      minPx=zeroValue;
      minPy=pyZeroValue;}

    //    std::cout<<" MtW2 from min py and min px "<< sqrt((minPy*minPy+minPx*minPx))*ptlep*2 -2*(pxlep*minPx + pylep*minPy)  <<std::endl;
    ///    ////Y part   

    double mu_Minimum = (mW*mW)/2 + minPx*pxlep + minPy*pylep;
    double a_Minimum  = (mu_Minimum*Lepton.pz())/(Lepton.energy()*Lepton.energy() - Lepton.pz()*Lepton.pz());
    pznu = a_Minimum;

    //if(!useMetForNegativeSolutions_){
      double Enu = sqrt(minPx*minPx+minPy*minPy + pznu*pznu);
      p4nu_rec.SetPxPyPzE(minPx, minPy, pznu , Enu);
    //}
    //else{
    //  pznu = a;
    //  double Enu = sqrt(metpx*metpx+metpy*metpy + pznu*pznu);
    //  p4nu_rec.SetPxPyPzE(metpx, metpy, pznu , Enu);
    //}
    result.push_back(p4nu_rec);
  }
  return result.front();
}

void ResetDeleteBranches(TTree* tree)
{
    TObjArray* branches = tree->GetListOfBranches();
    Int_t nb = branches->GetEntriesFast();
    for (Int_t i = 0; i < nb; ++i) {
        TBranch* br = (TBranch*) branches->UncheckedAt(i);
        if (br->InheritsFrom(TBranchElement::Class())) {
            ((TBranchElement*) br)->ResetDeleteObject();
        }
    }
}


////////////////////////////////////////////////////////////////////////////////
/// Main                                                                     ///
////////////////////////////////////////////////////////////////////////////////
void GrowTree_noReg(TString process, TString outname,  Long64_t beginEntry=0, Long64_t endEntry=-1)
{

  //    const float lumi = 20000.0;                                                                                                                                                                        

    gROOT->SetBatch(1);
    TH1::SetDefaultSumw2(1);
    gROOT->LoadMacro("HelperFunctions.h");  //< make functions visible to TTreeFormula

    TStopwatch sw;
    sw.Start();

 
    //  const TString indir   = "/afs/cern.ch/work/s/swang373/public/V21/";
    //    const TString indir   = "dcache:/pnfs/cms/WAX/resilient/jiafu/ZnunuHbb/skim_ZnnH_baseline/";
       const TString indir   = "/afs/cern.ch/user/d/degrutto/eos/cms/store/group/cmst3/user/degrutto/VHBBHeppyV21_add/";
    const TString outdir  = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/V21_FromSean/step3/";
    //const TString prefix  = "ZvvHighPt_V21_";
       const TString prefix  = "";
    const TString suffix  = ".root";

    TFile *input = TFile::Open(indir + prefix + process + suffix);
    if (!input) {
        std::cout << "ERROR: Could not open input file." << std::endl;
        exit(1);
    }
    /// Make output directory if it doesn't exist
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);

    std::cout << "--- GrowTree                 : Using input file: " << input->GetName() << std::endl;
    TTree *inTree = (TTree *) input->Get("tree");
    TH1F  *hcount = (TH1F *) input->Get("Count");
    TFile *output(0);
    if (beginEntry == 0 && endEntry == -1)
        output = TFile::Open(outdir + "Step3_" + outname + suffix, "RECREATE");
    else
        output = TFile::Open(outdir + "Step3_" + outname + TString::Format("_%Li_%Li", beginEntry, endEntry) + suffix, "RECREATE");
    TTree *outTree = inTree->CloneTree(0); // Do no copy the data yet
    /// The clone should not delete any shared i/o buffers.
    ResetDeleteBranches(outTree);
   
    float met_pt;

    inTree->SetBranchStatus("*", 1);
    inTree->SetBranchAddress("met_pt", &met_pt);                                                                                                                                                               

    ///-- Make new branches ----------------------------------------------------
    float  lumiNom;
    float efflumi,  xsec, weightTrigger; 

    outTree->Branch("lumiNom", &lumiNom, "lumiNom/F");
    outTree->Branch("efflumi", &efflumi, "efflumi/F");
    outTree->Branch("xsec", &xsec ,"xsec/F");
    outTree->Branch("weightTrig", &weightTrigger ,"weightTrigger/F");

    /// Get effective lumis
    std::map < std::string, float > xsecs = GetXsec();
    xsec = xsecs[outname.Data()];
    assert(xsec > 0);

    lumiNom       = lumiToScale * xsec / (hcount->GetBinContent(1));
    efflumi    =  lumiToScale * xsec / (hcount->GetBinContent(2)  -  hcount->GetBinContent(3));  // pos - def

    if (hcount->GetBinContent(2) == 0) {
      efflumi    =  lumiToScale * xsec / (hcount->GetBinContent(1));  // pos - def  
    }

    std::cout << "sample, efflumi, lumiToScale, xsec,  hcount->GetBinContent(1) , hcount->GetBinContent(2), hcount->GetBinContent(3) " << std::endl;
    std::cout << process.Data() << ", " << efflumi << " ," <<  lumiToScale<<" ," <<  xsec <<", " <<  hcount->GetBinContent(1) <<" , "  <<  hcount->GetBinContent(2) <<" , " <<  hcount->GetBinContent(3)  << std::endl;


    const Long64_t nentries = inTree->GetEntries();                                                                                                                                                 
    if (endEntry < 0)  endEntry = nentries;                                                                                                                                                          


    Long64_t ievt = 0;
    for (ievt=TMath::Max(ievt, beginEntry); ievt<TMath::Min(nentries, endEntry); ievt++) {
        if (ievt % 10000 == 0)
       
            std::cout << "--- ... Processing event: " << ievt << std::endl;
    
        const Long64_t local_entry = inTree->LoadTree(ievt);  // faster, but only for TTreeFormula
        if (local_entry < 0)  break;
        inTree->GetEntry(ievt);  // same event as received by LoadTree()
	if (met_pt < 150.) continue;
        weightTrigger = weightTrig(met_pt); 
	/*
        nalep_pt10_Znn = 0;
          
        /// Loop on additional leptons
        for (Int_t ial = 0; ial < naLeptons; ial++) {
	  if (aLeptons_pt[ial] > 10.  && aLeptons_looseIdPOG[ial] )) {
                if (deltaR(aLeptons_eta[ial], aLeptons_phi[ial], hJet_eta[0], hJet_phi[0]) > 0.3 && deltaR(aLepton_eta[ial], aLepton_phi[ial], hJet_eta[1], hJet_phi[1]) > 0.3) {
                    nalep_pt5_Znn++;
                }
            }
        }
..........................................
	*/
	outTree->Fill();  // fill it!
    
	}

  // end loop over TTree entries
    
    /// Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: ";
    sw.Print();

    output->cd();
    outTree->Write();
    output->Close();
    input->Close();

    delete input;
    delete output;
    /*
    delete reader;
    delete readerFJ;
    */
    

    std::cout << "==> GrowTree is done!" << std::endl << std::endl;
    return;
}
