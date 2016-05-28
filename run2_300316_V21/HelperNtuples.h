const std::string tagMC        = "root://eoscms//eos/cms//store/user/capalmer/VHBBHeppyNtuples/V7/";
const std::string tagData      = "root://eoscms//eos/cms//store/user/capalmer/VHBBHeppyNtuples/V7/";
//const std::string baseline     = "(Vtype==4||Vtype==3||Vtype==2) && met_pt>100 && Sum$(abs(Jet_eta)<2.5 && Jet_pt>20)>1 && Sum$(Jet_btagCSV>0.7)>0";
//const std::string baseline     = "met_pt>200 &&   (Jet_btagCSV[hJidx[0]]>0.9 ||  Jet_btagCSV[hJidx[1]]>0.9) && (Vtype==4) && Jet_pt[hJidx[0]]>30  && Jet_pt[hJidx[1]] > 30  && abs(Jet_eta[hJidx[0]])<2.4 &&  abs(Jet_eta[hJidx[1]])<2.4" ;
const std::string baseline     = "(Vtype==4||Vtype==3||Vtype==2) && met_pt>150" ;



//const std::string regression   = "(Vtype==4||Vtype==3||Vtype==2) && met_pt>120 && abs(hJet_eta[0])<2.5 && abs(hJet_eta[1])<2.5 && hJet_id[0]==1 && hJet_id[1]==1 && hJet_genPt[0]>10 && hJet_genPt[1]>10 && abs(hJet_flavour[0])==5 && hJet_pt[0]>20 && hJet_pt[1]>20 && hJet_csv[0]>0 && hJet_csv[1]>0";
const std::string regression   = baseline + "&& ( abs(Jet_mcFlavour[hJCidx[0]])==5  && abs(Jet_mcFlavour[hJCidx[1]])==5 )";
const std::string fjregression = "(Vtype==4||Vtype==3||Vtype==2) && METtype1corr.et>80 && FatH.FatHiggsFlag==1 && nfathFilterJets>0 && abs(fathFilterJets_eta[0])<2.5 && abs(fathFilterJets_eta[1])<2.5 && fathFilterJets_genPt[0]>10 && fathFilterJets_genPt[1]>10 && abs(fathFilterJets_flavour[0])==5 && fathFilterJets_pt[0]>15 && fathFilterJets_pt[1]>15 && fathFilterJets_csv[0]>0 && fathFilterJets_csv[1]>0";
const std::string mettrigger   = "1";
const std::string metfilter    = "1";
const std::string metbtagtrigger = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[49]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[50]==1)) )";
const std::string singlemutrigger = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=208686 && (triggerFlags[23]==1)) )";
const std::string singleeltrigger = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=208686 && (triggerFlags[44]==1)) )";


const float lumiToScale = 3000.0;


//const float lumi = 20000.0;
std::map<std::string, float> GetXsec() {
    std::map<std::string, float> values;
    values["ZnnH125"         ] =      (0.8696 -  0.1057) * 0.577 * 0.2 ;
    values["WlnH125"         ] =        1.380 * 0.577 * 0.1080 * 3. ;
    values["monoH_600"           ] =         0.026 ;
    values["monoH_800"           ] =        0.0288 ;
    values["monoH_1000"           ] =       0.02337 ;
    values["ggZH125"         ] =      2 * 0.1057 * 0.577 * 0.2  ;
    values["WJetsHT100"      ] =     1.23 *  1.29e3 ;
    values["WJetsHT200"      ] =     1.23 *  3.86e2 ;
    values["WJetsHT400"      ] =      1.23 * 47.9 ;
    values["WJetsHT600"      ] =       1.23 * 19.9 ;
    values["WJetsIncl"       ] =    61623.000000 ;
    values["ZJetsHT100"      ] =      409.860000 ;
    values["ZJetsHT200"      ] =      110.880000 ;
    values["ZJetsHT400"      ] =         13.189  ;
    values["ZJetsHT600"      ] =        4.524300 ;
    values["TTPow"       ] =      809.100000 ;
    values["TTMad"          ] =      809.100000 ;
    // W-> l nu = 0.1080 * 3 = 0.324, W -> qq = 0.676
    values["TTMadDiLep"          ] =  0.324 * 0.324 *    809.100000 ;
    values["TTMadSLepT"          ] =  0.5 * 0.324 * 0.676   *  809.100000 ;
    values["TTMadSLepTbar"          ] = 0.5 * 0.324 * 0.676 *  809.100000 ;
    // for singletop https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SingleTopRefXsec
    values["T_s"             ] =        6.35 ;
    values["T_t"             ] =      136.02 ;
    values["T_tW"            ] =       71.7 ;
    values["Tbar_s"          ] =        3.97 ;
    values["Tbar_t"          ] =      80.95;
    values["Tbar_tW"         ] =       23.310000 ;

    values["T_s_comb_lep"             ] =       10.32 *  0.324 ;
    values["T_t_comb_lep"             ] =       (118.44+64.470000) *  0.324 ;
    values["T_t_lep"             ] =      118.440000 * 0.324;
    values["Tbar_t_lep"          ] =       64.470000 * 0.324;



    values["QCDHT100"        ] =   2.75e7;
    values["QCDHT200"        ] =   1.74e6 ;
    values["QCDHT300"        ] =    3.67e5 ;
    values["QCDHT500"       ] =     2.94e4  ;
    values["QCDHT700"       ] =     6.524e3  ;
    values["QCDHT1000"       ] =      1.064e03 ;
    values["QCDHT1500"       ] =      121.5 ;
    values["QCDHT2000"       ] =      2.542e+01 ;
    values["TTMad"           ] =      809.100000 ;
    values["WW"           ] =      118.7 ;
    values["ZZ"           ] =      2.09 * 8.297 ;
    values["ZqqZnunu"           ] =      2.09 * 8.297 * 0.6991 * 0.20 ;
    values["WZ"           ] =      2.09 * 33.85 ;
    return values;


}

/*
for disobon
;-- sigma_VV taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV r9
; also check WW https://cdsweb.cern.ch/record/1460099/files/SMP-12-013-pas.pdf
; also check ZZ https://cdsweb.cern.ch/record/1460099/files/SMP-12-014-pas.pdf
; also check MCFM arXiv:1105.0020

;WW              = 56.7532            , 10000420.0 , 10007584.0000 , "#cccccc"
   ;WZ              = 33.85              , 240665 , 240665 , "#cccccc"
      ;ZZ              = 8.297              ,  9799897.0 ,  9809065.0000 , "#cccccc"

at 13TeV
https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
so WW  118.7 
118.7 / 56.7532 = 2.09

and take this ration for ZZ and WZ

*/


std::vector<std::string> GetLHECuts(const std::string id) {
 
    std::vector<std::string> values;
 
         if (id == "WJets") {
           values.resize(5, "");
        values[0] = "lheHT<100";
        values[1] = "100<=lheHT && lheHT<200";
        values[2] = "200<=lheHT && lheHT<400";
        values[3] = "400<=lheHT && lheHT<600";
        values[4] = "600<=lheHT";
        return values;
    }
         else if (id == "ZJets") {
   
        values.resize(4, "");
        values[0] = "100<=lheHT && lheHT<200";
        values[1] = "200<=lheHT && lheHT<400";
        values[2] = "400<=lheHT && lheHT<600";
        values[3] = "600<=lheHT";
        return values;
    }
    
    return values;
}

std::map<std::string, std::string> GetLHEWeights() {
    std::map<std::string, std::string> values;
    
    values["WJets"           ] = "((lheHT<100) * 5000.0 * 61623.000000 / 10017431.000000 * 1.000000) + ((100<=lheHT && lheHT<200) * 5000.0 * 2234.910000 / 5262249.000000 * 0.322231) + ((200<=lheHT && lheHT<400) * 5000.0 * 580.068000 / 4936055.000000 * 0.393106) + ((400<=lheHT && lheHT<600) * 5000.0 * 68.400300 / 4640575.000000 * 0.820838) + ((600<=lheHT) * 5000.0 * 23.136300 / 4643671.000000 * 0.926380)";
    values["ZJets"           ] = "((100<=lheHT && lheHT<200) * 5000.0 * 409.860000 / 4986410.000000 * 1.000000) + ((200<=lheHT && lheHT<400) * 5000.0 * 110.880000 / 4546455.000000 * 1.000000) + ((400<=lheHT && lheHT<600) * 5000.0 * 13.189000 / 4645939.000000 * 1.000000) + ((600<=lheHT) * 5000.0 * 4.524300 / 4463773.000000 * 1.000000)";
    
    return values;
}

