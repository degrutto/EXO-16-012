#include <iostream>
#include "TFile.h"
#include "TFileInfo.h"
#include "TTree.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TString.h"
#include "TCut.h"
#include "TH1.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFileCollection.h"
//#include "HelperNtuples.h"


//#include "HelperFunctions.h"

//#define XROOTD

void Skim(TString fname="ZnnH125", TString outputname = "")
{
  //    gROOT->LoadMacro("HelperFunctions.h" );  // make functions visible to TTreeFormula
    gROOT->SetBatch(1);

    TChain * chain  = new TChain("tree");
    TString dijet   = "";
    TString outdir  = "/afs/cern.ch/work/d/degrutto/public/MiniAOD/ZnnHbb_Spring15_PU20bx25/skimV12_v2/";
    //       TString outdir  = "";
    TString prefix  = "skim_";
    // TString prefix  = "";
    TString suffix  = "tree*.root";
    


    //    TCut selection = baseline.c_str();
    // Different baseline for dimuons
    // MET filters
    // selection += metfilter.c_str();

    // JSON & trigger
    /*   
 if (fname.Contains("Data")) {
        TCut trigger = mettrigger.c_str();
        selection += trigger;
        //selection.Print();
      }
    */

      /*
    } else if (process == "Data_METBTag_R" || process == "Data_METBTag_P") {
        TCut trigger = metbtagtrigger.c_str();
        selection += trigger;
        //selection.Print();
    } else if (process == "Data_SingleMu_R" || process == "Data_SingleMu_P") {
        TCut trigger = singlemutrigger.c_str();
        selection += trigger;
        //selection.Print();
    } else if (process == "Data_SingleEl_R" || process == "Data_SingleEl_P") {
        TCut trigger = singleeltrigger.c_str();
        selection += trigger;
        //selection.Print();
    }
      */


    chain->Add(fname);

    // Sum Count, CountWithPU, CountWithPUP, CountWithPUM
    TObjArray * files = chain->GetListOfFiles();
    TIter next(files);
    TChainElement * chainElem = 0;
    TFile * f2 = 0;
    TH1D * h1 = new TH1D("Count", ";Counts", 16, 0, 16);
    TH1F * htemp = 0;

    while ((chainElem = (TChainElement*) next())) {
      //#ifndef XROOTD
      //      f2 = TFile::Open("dcache:" + TString(chainElem->GetTitle()));
      //#else
      std::cout << "chainElem->GetTitle() " << chainElem->GetTitle() << std::endl; 
        f2 = TFile::Open( TString(chainElem->GetTitle()));
	//#endif
        htemp = (TH1F*) f2->Get("Count");
        h1->SetBinContent(1, h1->GetBinContent(1)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountPosWeight");
        h1->SetBinContent(2, h1->GetBinContent(2)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountNegWeight");
        h1->SetBinContent(3, h1->GetBinContent(3)+htemp->GetBinContent(1));

	/*
        h1->SetBinContent(4, h1->GetBinContent(4)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithMCProd");
        h1->SetBinContent(5, h1->GetBinContent(5)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("CountWithPUMCProd");
        h1->SetBinContent(6, h1->GetBinContent(6)+htemp->GetBinContent(1));
        htemp = (TH1F*) f2->Get("countWithSignalQCDcorrections");
        h1->SetBinContent(7, h1->GetBinContent(7)+htemp->GetBinContent(1));
        */
        std::clog << fname << ": skimmed from " << chainElem->GetTitle() << std::endl;
    }
    
    // LHE Count

    /*    
    TH1D * h2 = new TH1D("LHECount", ";LHE Counts", 16, 0, 16);
    TString process_lhe = fname;
    if (process_lhe.BeginsWith("WJets"))
        process_lhe = "WJets";
    else if (process_lhe.BeginsWith("ZJets"))
        process_lhe = "ZJets";
    else 
        process_lhe = "";
    const std::vector<std::string>& lhecuts = GetLHECuts(process_lhe.Data());
    for (unsigned int i=0; i < lhecuts.size(); i++) {
        TCut cut2 = lhecuts.at(i).c_str();
        h2->SetBinContent(i+1, chain->GetEntries(cut2));
    }
    
    */

    // Make output directory if it doesn't exist
    if (gSystem->AccessPathName(outdir))
        gSystem->mkdir(outdir);
    TString outname = outdir + prefix + Form("%s", outputname.Data());
    std::cout << "outname is " << outname << std::endl;

    TFile* f1 = TFile::Open(outname, "RECREATE");
    //    TTree* t1 = (TTree*) chain->CopyTree(selection);
    // the selection is bypassed, loop instead on the events to be faster
   
  TTree* t1 = (TTree*) chain->CloneTree(0);
  

    Long64_t nentries = chain->GetEntries();
    float met   = 0;
    chain->SetBranchAddress("met_pt",&met);

   
    for (Long64_t i=0;i<nentries; i++) {
      chain->GetEntry(i);
      if (met > 150.) t1->Fill();
      
    }


    std::clog << fname << ": skimmed from " << chain->GetEntriesFast() << " to " << t1->GetEntriesFast() << " entries." << std::endl;
    //    t1->Print();    
    t1->AutoSave();    
    t1->Write();
    h1->Write();
    // h2->Write();
    f1->Close();

 
    delete f1;  
    delete f2;  
    return;
}

