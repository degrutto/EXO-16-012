/// TMVA Regression
std::vector<std::string> GetInputExpressionsReg() {
    std::vector<std::string> values;
    values.resize(17, "");
    values[0] = "breg_rawptJER := smear_pt_res(Jet_rawPt[hJCidx], Jet_mcPt[hJCidx], Jet_eta[hJCidx])";
    values[1] = "breg_pt := Jet_pt[hJCidx]";
    values[2] = "breg_eta := Jet_eta[hJCidx]";
    values[3] = "breg_leadtrackpt := max(0,Jet_leadTrackPt[hJCidx])";
    values[4] = "breg_vtx3dL := max(0,Jet_vtx3DVal[hJCidx])";
    values[5] = "breg_vtx3deL := max(0,Jet_vtx3DSig[hJCidx])";
    values[6] = "breg_vtxMass := max(0,Jet_vtxMass[hJCidx])";
    values[7] = "breg_vtxPt := max(0,Jet_vtxPt[hJCidx])";
    values[8] = "breg_chf := Jet_chHEF[hJCidx]";
    values[9] = "breg_cef := Jet_chEmEF[hJCidx]";
    values[10] = "breg_nhf := Jet_neHEF[hJCidx]";
    values[11] = "breg_nef := Jet_neHEF[hJCidx]";
    values[12] = "breg_nch := Jet_chMult[hJCidx]";
    values[13] = "breg_ntot := Jet_mult[hJCidx]";
    values[14] = "breg_softlepptrel := max(0,Jet_leptonPtRel[hJCidx])";
    values[15] = "breg_softleppt := max(0,Jet_leptonPt[hJCidx])";
    values[16] = "breg_softlepdR := max(0,Jet_leptonDeltaR[hJCidx])";
    return values;
}

std::vector<std::string> GetInputExpressionLabelsReg() {
    std::vector<std::string> values;
    values.resize(17, "");
    values[0] = "raw (JER corr.) p_{T}(j);GeV;F";
    values[1] = "p_{T}(j);GeV;F";
    values[2] = "#eta(j);;F";
    values[3] = "leading track p_{T}(j);GeV;F";
    values[4] = "sec vtx 3-d flight length(j);cm;F";
    values[5] = "sec vtx 3-d flight length error(j);cm;F";
    values[6] = "sec vtx mass(j);GeV;F";
    values[7] = "sec vtx p_{T}(j);GeV;F";
    values[8] = "charged hadronic frac(j);;F";
    values[9] = "charged electromagnetic frac(j);;F";
    values[10] = "neutral hadronic frac(j);;F";
    values[11] = "neutral electromagnetic frac(j);;F";
    values[12] = "# charged particles(j);;I";
    values[13] = "# particles(j);;I";
    values[14] = "p_{T, rel}(lep-in-j);GeV;F";
    values[15] = "p_{T}(lep-in-j);GeV;F";
    values[16] = "#DeltaR(lep-in-j);;F";
    return values;
}

std::vector<std::string> GetInputExpressionsReg0() {
    std::vector<std::string> values;
    values.resize(17, "");
    values[0] = "smear_pt_res(Jet_rawPt[hJCidx[0]], Jet_mcPt[hJCidx[0]], Jet_eta[hJCidx[0]])";
    values[1] = "Jet_pt[hJCidx[0]]";
    values[2] = "Jet_eta[hJCidx[0]]";
    values[3] = "max(0,Jet_leadTrackPt[hJCidx[0]])";
    values[4] = "max(0,Jet_vtx3DVal[hJCidx[0]])";
    values[5] = "max(0,Jet_vtx3DSig[hJCidx[0]])";
    values[6] = "max(0,Jet_vtxMass[hJCidx[0]])";
    values[7] = "max(0,Jet_vtxPt[hJCidx[0]])";
    values[8] = "Jet_chHEF[hJCidx[0]]";
    values[9] = "Jet_chEmEF[hJCidx[0]]";
    values[10] = "Jet_neHEF[hJCidx[0]]";
    values[11] = "Jet_neHEF[hJCidx[0]]";
    values[12] = "Jet_chMult[hJCidx[0]]";
    values[13] = "Jet_mult[hJCidx[0]]";
    values[14] = "max(0,Jet_leptonPtRel[hJCidx[0]])";
    values[15] = "max(0,Jet_leptonPt[hJCidx[0]])";
    values[16] = "max(0,Jet_leptonDeltaR[hJCidx[0]])";
    return values;
}

std::vector<std::string> GetInputExpressionsReg1() {
    std::vector<std::string> values;
    values.resize(17, "");
    values[0] = "smear_pt_res(Jet_rawPt[hJCidx[1]], Jet_mcPt[hJCidx[1]], Jet_eta[hJCidx[1]])";
    values[1] = "Jet_pt[hJCidx[1]]";
    values[2] = "Jet_eta[hJCidx[1]]";
    values[3] = "max(0,Jet_leadTrackPt[hJCidx[1]])";
    values[4] = "max(0,Jet_vtx3DVal[hJCidx[1]])";
    values[5] = "max(0,Jet_vtx3DSig[hJCidx[1]])";
    values[6] = "max(0,Jet_vtxMass[hJCidx[1]])";
    values[7] = "max(0,Jet_vtxPt[hJCidx[1]])";
    values[8] = "Jet_chHEF[hJCidx[1]]";
    values[9] = "Jet_chEmEF[hJCidx[1]]";
    values[10] = "Jet_neHEF[hJCidx[1]]";
    values[11] = "Jet_neHEF[hJCidx[1]]";
    values[12] = "Jet_chMult[hJCidx[1]]";
    values[13] = "Jet_mult[hJCidx[1]]";
    values[14] = "max(0,Jet_leptonPtRel[hJCidx[1]])";
    values[15] = "max(0,Jet_leptonPt[hJCidx[1]])";
    values[16] = "max(0,Jet_leptonDeltaR[hJCidx[1]])";
    return values;
}

/// TMVA Regression (filter jets)
std::vector<std::string> GetInputExpressionsFJReg() {
    std::vector<std::string> values;
    values.resize(12, "");
    values[0] = "breg_fjrawptJER := smear_pt_res(fathFilterJets_ptRaw, fathFilterJets_genPt, fathFilterJets_eta)";
    values[1] = "breg_fjpt := fathFilterJets_pt";
    values[2] = "breg_fjet := evalEt(fathFilterJets_pt, fathFilterJets_eta, fathFilterJets_phi, fathFilterJets_e)";
    values[3] = "breg_fjmt := evalMt(fathFilterJets_pt, fathFilterJets_eta, fathFilterJets_phi, fathFilterJets_e)";
    values[4] = "breg_fjleadtrackpt := max(0,fathFilterJets_ptLeadTrack)";
    values[5] = "breg_fjvtx3dL := max(0,fathFilterJets_vtx3dL)";
    values[6] = "breg_fjvtx3deL := max(0,fathFilterJets_vtx3deL)";
    values[7] = "breg_fjvtxMass := max(0,fathFilterJets_vtxMass)";
    values[8] = "breg_fjvtxPt := max(0,fathFilterJets_vtxPt)";
    values[9] = "breg_fjntot := fathFilterJets_nconstituents";
    values[10] = "breg_fjeJEC := fathFilterJets_JECUnc";
    values[11] = "breg_fjarea := fathFilterJets_jetArea";
    return values;
}

std::vector<std::string> GetInputExpressionLabelsFJReg() {
    std::vector<std::string> values;
    values.resize(12, "");
    values[0] = "raw (JER corr.) p_{T}(fj);GeV;F";
    values[1] = "p_{T}(fj);GeV;F";
    values[2] = "E_{T}(fj);GeV;F";
    values[3] = "M_{T}(fj);GeV;F";
    values[4] = "leading track p_{T}(fj);GeV;F";
    values[5] = "sec vtx 3-d flight length(fj);cm;F";
    values[6] = "sec vtx 3-d flight length error(fj);cm;F";
    values[7] = "sec vtx mass(fj);GeV;F";
    values[8] = "sec vtx p_{T}(fj);GeV;F";
    values[9] = "# particles(fj);;I";
    values[10] = "JEC uncertainty(fj);GeV;F";
    values[11] = "jet area(fj);;F";
    return values;
}

std::vector<std::string> GetInputExpressionsFJReg0() {
    std::vector<std::string> values;
    values.resize(12, "");
    values[0] = "smear_pt_res(fathFilterJets_ptRaw[0], fathFilterJets_genPt[0], fathFilterJets_eta[0])";
    values[1] = "fathFilterJets_pt[0]";
    values[2] = "evalEt(fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0])";
    values[3] = "evalMt(fathFilterJets_pt[0], fathFilterJets_eta[0], fathFilterJets_phi[0], fathFilterJets_e[0])";
    values[4] = "max(0,fathFilterJets_ptLeadTrack[0])";
    values[5] = "max(0,fathFilterJets_vtx3dL[0])";
    values[6] = "max(0,fathFilterJets_vtx3deL[0])";
    values[7] = "max(0,fathFilterJets_vtxMass[0])";
    values[8] = "max(0,fathFilterJets_vtxPt[0])";
    values[9] = "fathFilterJets_nconstituents[0]";
    values[10] = "fathFilterJets_JECUnc[0]";
    values[11] = "fathFilterJets_jetArea[0]";
    return values;
}

std::vector<std::string> GetInputExpressionsFJReg1() {
    std::vector<std::string> values;
    values.resize(12, "");
    values[0] = "smear_pt_res(fathFilterJets_ptRaw[1], fathFilterJets_genPt[1], fathFilterJets_eta[1])";
    values[1] = "fathFilterJets_pt[1]";
    values[2] = "evalEt(fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1])";
    values[3] = "evalMt(fathFilterJets_pt[1], fathFilterJets_eta[1], fathFilterJets_phi[1], fathFilterJets_e[1])";
    values[4] = "max(0,fathFilterJets_ptLeadTrack[1])";
    values[5] = "max(0,fathFilterJets_vtx3dL[1])";
    values[6] = "max(0,fathFilterJets_vtx3deL[1])";
    values[7] = "max(0,fathFilterJets_vtxMass[1])";
    values[8] = "max(0,fathFilterJets_vtxPt[1])";
    values[9] = "fathFilterJets_nconstituents[1]";
    values[10] = "fathFilterJets_JECUnc[1]";
    values[11] = "fathFilterJets_jetArea[1]";
    return values;
}

std::vector<std::string> GetInputExpressionsFJReg2() {
    std::vector<std::string> values;
    values.resize(12, "");
    values[0] = "smear_pt_res(fathFilterJets_ptRaw[2], fathFilterJets_genPt[2], fathFilterJets_eta[2])";
    values[1] = "fathFilterJets_pt[2]";
    values[2] = "evalEt(fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2])";
    values[3] = "evalMt(fathFilterJets_pt[2], fathFilterJets_eta[2], fathFilterJets_phi[2], fathFilterJets_e[2])";
    values[4] = "max(0,fathFilterJets_ptLeadTrack[2])";
    values[5] = "max(0,fathFilterJets_vtx3dL[2])";
    values[6] = "max(0,fathFilterJets_vtx3deL[2])";
    values[7] = "max(0,fathFilterJets_vtxMass[2])";
    values[8] = "max(0,fathFilterJets_vtxPt[2])";
    values[9] = "fathFilterJets_nconstituents[2]";
    values[10] = "fathFilterJets_JECUnc[2]";
    values[11] = "fathFilterJets_jetArea[2]";
    return values;
}

/// TMVA Classification
//! BDT regression is applied
std::map<std::string, std::string> GetWeightExpressions() {
    std::map<std::string, std::string> values;
    //    values["ZnunuHighPt"     ] = "efflumi * PUweight * triggercorr2012ABCD( (triggerFlags[42]==1 || triggerFlags[39]==1), (triggerFlags[41]==1), METtype1corr.et, max(hJet_csv_nominal[0],hJet_csv_nominal[1]) )";
    values["ZnunuHighPt"     ] = "efflumi * weightTrig * sign(genWeight)";
    values["ZnunuMedPt"      ] = "efflumi * PUweight * triggercorr2012ABCD( (triggerFlags[42]==1 || triggerFlags[39]==1), (triggerFlags[41]==1), METtype1corr.et, max(hJet_csv_nominal[0],hJet_csv_nominal[1]) )";
    values["ZnunuLowPt"      ] = "efflumi * PUweight * triggercorr2012ABCD( (triggerFlags[42]==1 || triggerFlags[39]==1), (triggerFlags[41]==1), METtype1corr.et, max(hJet_csv_nominal[0],hJet_csv_nominal[1]) )";
    values["ZnunuLowCSV"     ] = "efflumi * PUweight * triggercorr2012ABCD( (triggerFlags[42]==1 || triggerFlags[39]==1), (triggerFlags[41]==1), METtype1corr.et, max(hJet_csv_nominal[0],hJet_csv_nominal[1]) )";
    return values;
}

std::map<std::string, std::string> GetTriggerExpressions() {
    std::map<std::string, std::string> values;
    //    values["ZnunuHighPt"     ] = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[42]==1 || triggerFlags[49]==1 || triggerFlags[40]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1)) )";
    values["ZnunuHighPt"     ] = "";
    values["ZnunuMedPt"      ] = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[42]==1 || triggerFlags[49]==1 || triggerFlags[40]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1)) )";
    values["ZnunuLowPt"      ] = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[42]==1 || triggerFlags[49]==1 || triggerFlags[40]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1)) )";
    values["ZnunuLowCSV"     ] = "EVENT.json && ( (190456<=EVENT.run && EVENT.run<=193752 && (triggerFlags[42]==1 || triggerFlags[49]==1 || triggerFlags[40]==1)) || (193752<=EVENT.run && EVENT.run<=208686 && (triggerFlags[42]==1 || triggerFlags[39]==1 || triggerFlags[41]==1)) )";
    return values;
}

std::vector<std::string> GetSelExpressions(const std::string id) {
    std::vector<std::string> values;
    if (id == "ZnunuHighPt") {
        values.resize(6, "");
	//        values[0] = "ZnunuHighPt_VH := (Vtype==4) && METtype1corr.et>170 && HptReg>130 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && mindPhiMETJet_dPhi>0.5 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_pt5_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";

			values[0] = "ZnunuHighPt_VH := (Vtype==4) &&   Jet_btagCSV[hJCidx[0]]>0.97 && Jet_btagCSV[hJCidx[1]]>0.423   &&  met_pt>170 && HCSV_pt>130 && ( (Jet_pt[hJCidx[0]]>80 && Jet_pt[hJCidx[1]]>30) || (Jet_pt[hJCidx[1]]>80 && Jet_pt[hJCidx[0]]>30 ))   &&  MinIf$(abs(deltaPhi(met_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>0.5 && abs(deltaPhi(HCSV_phi, met_phi))>2.5"; // <-- remove ZnunuHighPt_VH : = for now for BDT skim classidication
			//	values[0] = "(Vtype==4) &&   Jet_btagCSV[hJCidx[0]]>0.97 && Jet_btagCSV[hJCidx[1]]>0.423   &&  met_pt>170 && HCSV_pt>130 && ( (Jet_pt[hJCidx[0]]>80 && Jet_pt[hJCidx[1]]>30) || (Jet_pt[hJCidx[1]]>80 && Jet_pt[hJCidx[0]]>30 ))   &&  MinIf$(abs(deltaPhi(met_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>0.5 && abs(deltaPhi(HCSV_phi, met_phi))>2.5"; // <-- remove ZnunuHighPt_VH : = for now for BDT skim classidication

        values[1] = "ZnunuHighPt_ZjLF := (Vtype==4) && (Vtype==4) &&  met_pt>170 &&  max(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])>0. && min(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])>0. && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId>0)==2 && MinIf$(abs(deltaPhi(met_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>0.5";
        values[2] = "ZnunuHighPt_ZjHF := (Vtype==4) &&   met_pt>170 && max(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])>0.814 && min(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])>0.423 && (HCSV_mass<90 || HCSV_mass>150) && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId>0)==2 && MinIf$(abs(deltaPhi(met_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>0.5";
        values[3] = "ZnunuHighPt_WjLF := (Vtype==2||Vtype==3)  &&  met_pt>170 && max(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])>0. && min(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])>0. && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId>0)==2 &&  MinIf$(abs(deltaPhi(met_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>0.5";
        values[4] = "ZnunuHighPt_WjHF := (Vtype==2||Vtype==3) && max(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])>0.814 && min(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])>0.423 && (HCSV_mass<90 || HCSV_mass>150) && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId>0)==2 &&  MinIf$(abs(deltaPhi(met_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>0.5 && met_pt>170";
        values[5] = "ZnunuHighPt_TT := (Vtype==2||Vtype==3) && max(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])>0.814 && min(Jet_btagCSV[hJCidx[0]],Jet_btagCSV[hJCidx[1]])>0.423 && Sum$(Jet_pt>30 && abs(Jet_eta)<4.5 && Jet_puId>0)>2 &&  MinIf$(abs(deltaPhi(met_phi,Jet_phi)) , Jet_pt>25 && abs(Jet_eta)<4.5 && Jet_puId==1)>0.5 && met_pt>170";
        return values;
    }
    else if (id == "ZnunuMedPt") {
        values.resize(6, "");
        values[0] = "ZnunuMedPt_VH := (Vtype==4) && 130<METtype1corr.et && METtype1corr.et<=170 && HptReg>130 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_pt5_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[1] = "ZnunuMedPt_ZjLF := (Vtype==4) && 130<METtype1corr.et && METtype1corr.et<=170 && HptReg>130 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])<=0.898 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && hJet_nhf[0]<0.90 && hJet_nhf[1]<0.90 && hJet_nef[0]<0.90 && hJet_nef[1]<0.90 && nalep_pt5_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[2] = "ZnunuMedPt_ZjHF := (Vtype==4) && 130<METtype1corr.et && METtype1corr.et<=170 && HptReg>130 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && (HmassReg<100 || 140<HmassReg) && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_pt5_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[3] = "ZnunuMedPt_WjLF := (Vtype==2||Vtype==3) && 130<METtype1corr.et && METtype1corr.et<=170 && HptReg>130 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])<=0.898 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && naJets_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[4] = "ZnunuMedPt_WjHF := (Vtype==2||Vtype==3) && 130<METtype1corr.et && METtype1corr.et<=170 && HptReg>130 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && (HmassReg<100 || 140<HmassReg) && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && naJets_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[5] = "ZnunuMedPt_TT := (Vtype==2||Vtype==3) && 130<METtype1corr.et && METtype1corr.et<=170 && HptReg>130 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && (HmassReg<100 || 140<HmassReg) && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && naJets_Znn>=1 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        return values;
    }
    else if (id == "ZnunuLowPt") {
        values.resize(6, "");
        values[0] = "ZnunuLowPt_VH := (Vtype==4) && 100<METtype1corr.et && METtype1corr.et<=130 && HptReg>100 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && METtype1corr.et/sqrt(METtype1corr.sumet)>3 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && naJets_Znn<=1 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_pt5_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[1] = "ZnunuLowPt_ZjLF := (Vtype==4) && 100<METtype1corr.et && METtype1corr.et<=130 && HptReg>100 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])<=0.898 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && METtype1corr.et/sqrt(METtype1corr.sumet)>3 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && naJets_Znn<=1 && hJet_nhf[0]<0.90 && hJet_nhf[1]<0.90 && hJet_nef[0]<0.90 && hJet_nef[1]<0.90 && nalep_pt5_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[2] = "ZnunuLowPt_ZjHF := (Vtype==4) && 100<METtype1corr.et && METtype1corr.et<=130 && HptReg>100 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && (HmassReg<100 || 140<HmassReg) && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && METtype1corr.et/sqrt(METtype1corr.sumet)>3 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && naJets_Znn<=1 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_pt5_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[3] = "ZnunuLowPt_WjLF := (Vtype==2||Vtype==3) && 100<METtype1corr.et && METtype1corr.et<=130 && HptReg>100 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])<=0.898 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && METtype1corr.et/sqrt(METtype1corr.sumet)>3 && naJets_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[4] = "ZnunuLowPt_WjHF := (Vtype==2||Vtype==3) && 100<METtype1corr.et && METtype1corr.et<=130 && HptReg>100 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.244 && (HmassReg<100 || 140<HmassReg) && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && METtype1corr.et/sqrt(METtype1corr.sumet)>3 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && naJets_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[5] = "ZnunuLowPt_TT := (Vtype==2||Vtype==3) && 100<METtype1corr.et && METtype1corr.et<=130 && HptReg>100 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && (HmassReg<100 || 140<HmassReg) && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && METtype1corr.et/sqrt(METtype1corr.sumet)>3 && naJets_Znn>=1 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        return values;
    }
    else if (id == "ZnunuLowCSV") {
        values.resize(3, "");
        values[0] = "ZnunuLowCSV_VH := (Vtype==4) && METtype1corr.et>170 && HptReg>130 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])<=0.244 && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_pt5_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[1] = "ZnunuLowCSV_ZjHF := (Vtype==4) && METtype1corr.et>170 && HptReg>130 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])<=0.244 && (HmassReg<100 || 140<HmassReg) && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && abs(deltaPhi(METtype1corr.phi,METnoPUCh.phi))<0.5 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && nalep_pt5_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        values[2] = "ZnunuLowCSV_WjHF := (Vtype==2||Vtype==3) && METtype1corr.et>170 && HptReg>130 && max(hJet_pt[0],hJet_pt[1])>60 && min(hJet_pt[0],hJet_pt[1])>30 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.679 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])<=0.244 && (HmassReg<100 || 140<HmassReg) && HmassReg<250 && mindPhiMETJet_dPhi>0.7 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.0 && naJets_Znn==0 && hbhe && ecalFlag && cschaloFlag && hcallaserFlag && trackingfailureFlag && eebadscFlag && !isBadHcalEvent && hJet_puJetIdL[0]>0 && hJet_puJetIdL[1]>0 && HmassReg>0";
        return values;
    }
    return values;
}

std::vector<std::string> GetSelMjjExpressions(const std::string id) {
    std::vector<std::string> values;
    if (id == "ZnunuHighPt") {
        values.resize(1, "");
        values[0] = "ZnunuHighPt_VH := HptReg>160 && hJet_pt[0]>80 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.5 && naJets_Znn==0 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.95 && mindPhiMETJet_dPhi>1.5 && HmassReg>0";
        return values;
    }
    else if (id == "ZnunuMedPt") {
        values.resize(1, "");
        values[0] = "ZnunuMedPt_VH := HptReg>130 && hJet_pt[0]>60 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.5 && naJets_Znn==0 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.95 && mindPhiMETJet_dPhi>1.5 && HmassReg>0";
        return values;
    }
    else if (id == "ZnunuLowPt") {
        values.resize(1, "");
        values[0] = "ZnunuLowPt_VH := HptReg>110 && hJet_pt[0]>60 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.5 && naJets_Znn==0 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.95 && mindPhiMETJet_dPhi>1.5 && HmassReg>0";
        return values;
    }
    else if (id == "ZnunuLowCSV") {
        values.resize(1, "");
        values[0] = "ZnunuLowCSV_VH := HptReg>160 && hJet_pt[0]>80 && max(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.898 && min(hJet_csv_nominal[0],hJet_csv_nominal[1])>0.5 && naJets_Znn==0 && abs(deltaPhi(H.phi,METtype1corr.phi))>2.95 && mindPhiMETJet_dPhi>1.5 && HmassReg>0";
        return values;
    }
    return values;
}

TObjString selectFlags_id = " \n//ZnunuHighPt \nselectFlags[0][isyst] = \"ZnunuHighPt_VH\"; \nselectFlags[1][isyst] = \"ZnunuHighPt_ZjLF\"; \nselectFlags[2][isyst] = \"ZnunuHighPt_ZjHF\"; \nselectFlags[3][isyst] = \"ZnunuHighPt_WjLF\"; \nselectFlags[4][isyst] = \"ZnunuHighPt_WjHF\"; \nselectFlags[5][isyst] = \"ZnunuHighPt_TT\"; \n//ZnunuMedPt \nselectFlags[0][isyst] = \"ZnunuMedPt_VH\"; \nselectFlags[1][isyst] = \"ZnunuMedPt_ZjLF\"; \nselectFlags[2][isyst] = \"ZnunuMedPt_ZjHF\"; \nselectFlags[3][isyst] = \"ZnunuMedPt_WjLF\"; \nselectFlags[4][isyst] = \"ZnunuMedPt_WjHF\"; \nselectFlags[5][isyst] = \"ZnunuMedPt_TT\"; \n//ZnunuLowPt \nselectFlags[0][isyst] = \"ZnunuLowPt_VH\"; \nselectFlags[1][isyst] = \"ZnunuLowPt_ZjLF\"; \nselectFlags[2][isyst] = \"ZnunuLowPt_ZjHF\"; \nselectFlags[3][isyst] = \"ZnunuLowPt_WjLF\"; \nselectFlags[4][isyst] = \"ZnunuLowPt_WjHF\"; \nselectFlags[5][isyst] = \"ZnunuLowPt_TT\"; \n//ZnunuLowCSV \nselectFlags[0][isyst] = \"ZnunuLowCSV_VH\"; \nselectFlags[1][isyst] = \"ZnunuLowCSV_ZjHF\"; \nselectFlags[2][isyst] = \"ZnunuLowCSV_WjHF\"; ";


std::vector<std::vector<std::pair<std::string, std::string> > > GetSystExpressions() {                                                                                                               
  std::vector<std::vector<std::pair<std::string, std::string> > > values;                                                                                                                          
  values.resize(1);
  values[0].resize(1); // NONE

  return values;

 }
/*
std::vector<std::vector<std::pair<std::string, std::string> > > GetSystExpressions() {
    std::vector<std::vector<std::pair<std::string, std::string> > > values;
    values.resize(15);
    values[0].resize(1); // NONE
    values[0][0] = make_pair("hJet_ptReg", "hJet_ptReg");
    values[1].resize(10); // CMS_vhbb_res_j_Znn_8TeVUp
    values[1][0] = make_pair("hJet_ptReg", "hJet_ptReg_res_j_up");
    values[1][1] = make_pair("HptReg", "HptReg_res_j_up");
    values[1][2] = make_pair("HmassReg", "HmassReg_res_j_up");
    values[1][3] = make_pair("aJet_pt", "aJet_pt_res_j_up");
    values[1][4] = make_pair("METtype1corr.et", "MET_res_j_up");
    values[1][5] = make_pair("METtype1corr.phi", "METphi_res_j_up");
    values[1][6] = make_pair("naJets_Znn", "naJets_Znn_res_j_up");
    values[1][7] = make_pair("mindPhiMETJet_dPhi", "mindPhiMETJet_dPhi_res_j_up");
    values[1][8] = make_pair("mindRHJet_dR", "mindRHJet_dR_res_j_up");
    values[1][9] = make_pair("FatHmassReg", "FatHmassReg_res_j_up");
    values[2].resize(10); // CMS_vhbb_res_j_Znn_8TeVDown
    values[2][0] = make_pair("hJet_ptReg", "hJet_ptReg_res_j_down");
    values[2][1] = make_pair("HptReg", "HptReg_res_j_down");
    values[2][2] = make_pair("HmassReg", "HmassReg_res_j_down");
    values[2][3] = make_pair("aJet_pt", "aJet_pt_res_j_down");
    values[2][4] = make_pair("METtype1corr.et", "MET_res_j_down");
    values[2][5] = make_pair("METtype1corr.phi", "METphi_res_j_down");
    values[2][6] = make_pair("naJets_Znn", "naJets_Znn_res_j_down");
    values[2][7] = make_pair("mindPhiMETJet_dPhi", "mindPhiMETJet_dPhi_res_j_down");
    values[2][8] = make_pair("mindRHJet_dR", "mindRHJet_dR_res_j_down");
    values[2][9] = make_pair("FatHmassReg", "FatHmassReg_res_j_down");
    values[3].resize(10); // CMS_vhbb_scale_j_Znn_8TeVUp
    values[3][0] = make_pair("hJet_ptReg", "hJet_ptReg_scale_j_up");
    values[3][1] = make_pair("HptReg", "HptReg_scale_j_up");
    values[3][2] = make_pair("HmassReg", "HmassReg_scale_j_up");
    values[3][3] = make_pair("aJet_pt", "aJet_pt_scale_j_up");
    values[3][4] = make_pair("METtype1corr.et", "MET_scale_j_up");
    values[3][5] = make_pair("METtype1corr.phi", "METphi_scale_j_up");
    values[3][6] = make_pair("naJets_Znn", "naJets_Znn_scale_j_up");
    values[3][7] = make_pair("mindPhiMETJet_dPhi", "mindPhiMETJet_dPhi_scale_j_up");
    values[3][8] = make_pair("mindRHJet_dR", "mindRHJet_dR_scale_j_up");
    values[3][9] = make_pair("FatHmassReg", "FatHmassReg_scale_j_up");
    values[4].resize(10); // CMS_vhbb_scale_j_Znn_8TeVDown
    values[4][0] = make_pair("hJet_ptReg", "hJet_ptReg_scale_j_down");
    values[4][1] = make_pair("HptReg", "HptReg_scale_j_down");
    values[4][2] = make_pair("HmassReg", "HmassReg_scale_j_down");
    values[4][3] = make_pair("aJet_pt", "aJet_pt_scale_j_down");
    values[4][4] = make_pair("METtype1corr.et", "MET_scale_j_down");
    values[4][5] = make_pair("METtype1corr.phi", "METphi_scale_j_down");
    values[4][6] = make_pair("naJets_Znn", "naJets_Znn_scale_j_down");
    values[4][7] = make_pair("mindPhiMETJet_dPhi", "mindPhiMETJet_dPhi_scale_j_down");
    values[4][8] = make_pair("mindRHJet_dR", "mindRHJet_dR_scale_j_down");
    values[4][9] = make_pair("FatHmassReg", "FatHmassReg_scale_j_down");
    values[5].resize(2); // CMS_vhbb_eff_bUp
    values[5][0] = make_pair("hJet_csv_nominal", "hJet_csv_upBC");
    values[5][1] = make_pair("aJet_csv_nominal", "aJet_csv_upBC");
    values[6].resize(2); // CMS_vhbb_eff_bDown
    values[6][0] = make_pair("hJet_csv_nominal", "hJet_csv_downBC");
    values[6][1] = make_pair("aJet_csv_nominal", "aJet_csv_downBC");
    values[7].resize(2); // CMS_vhbb_fake_b_8TeVUp
    values[7][0] = make_pair("hJet_csv_nominal", "hJet_csv_upL");
    values[7][1] = make_pair("aJet_csv_nominal", "aJet_csv_upL");
    values[8].resize(2); // CMS_vhbb_fake_b_8TeVDown
    values[8][0] = make_pair("hJet_csv_nominal", "hJet_csv_downL");
    values[8][1] = make_pair("aJet_csv_nominal", "aJet_csv_downL");
    values[9].resize(2); // UEPSUp
    values[9][0] = make_pair("PUweight", "PUweightP");
    values[9][1] = make_pair("efflumi", "efflumi_UEPS_up");
    values[10].resize(2); // UEPSDown
    values[10][0] = make_pair("PUweight", "PUweightM");
    values[10][1] = make_pair("efflumi", "efflumi_UEPS_down");
    values[11].resize(2); // CMS_vhbb_trigger_MET_Znn_8TeVUp
    values[11][0] = make_pair("triggercorr2012ABCD", "triggercorr2012ABCD_MET_up");
    values[11][1] = make_pair("triggerweight2012ABCD", "triggerweight2012ABCD_up");
    values[12].resize(2); // CMS_vhbb_trigger_MET_Znn_8TeVDown
    values[12][0] = make_pair("triggercorr2012ABCD", "triggercorr2012ABCD_MET_down");
    values[12][1] = make_pair("triggerweight2012ABCD", "triggerweight2012ABCD_down");
    values[13].resize(1); // CMS_vhbb_trigger_METCSV_Znn_8TeVUp
    values[13][0] = make_pair("triggercorr2012ABCD", "triggercorr2012ABCD_CSV_up");
    values[14].resize(1); // CMS_vhbb_trigger_METCSV_Znn_8TeVDown
    values[14][0] = make_pair("triggercorr2012ABCD", "triggercorr2012ABCD_CSV_down");
    return values;
}
*/

TObjString systFlags_id = " \nsystFlags[0] = \"NONE\"; \nsystFlags[1] = \"CMS_vhbb_res_j_Znn_8TeVUp\"; \nsystFlags[2] = \"CMS_vhbb_res_j_Znn_8TeVDown\"; \nsystFlags[3] = \"CMS_vhbb_scale_j_Znn_8TeVUp\"; \nsystFlags[4] = \"CMS_vhbb_scale_j_Znn_8TeVDown\"; \nsystFlags[5] = \"CMS_vhbb_eff_bUp\"; \nsystFlags[6] = \"CMS_vhbb_eff_bDown\"; \nsystFlags[7] = \"CMS_vhbb_fake_b_8TeVUp\"; \nsystFlags[8] = \"CMS_vhbb_fake_b_8TeVDown\"; \nsystFlags[9] = \"UEPSUp\"; \nsystFlags[10] = \"UEPSDown\"; \nsystFlags[11] = \"CMS_vhbb_trigger_MET_Znn_8TeVUp\"; \nsystFlags[12] = \"CMS_vhbb_trigger_MET_Znn_8TeVDown\"; \nsystFlags[13] = \"CMS_vhbb_trigger_METCSV_Znn_8TeVUp\"; \nsystFlags[14] = \"CMS_vhbb_trigger_METCSV_Znn_8TeVDown\"; ";

