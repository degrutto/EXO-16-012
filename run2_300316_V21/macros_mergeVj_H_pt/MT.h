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



inline double evalHMETMassiveMt(double m0, double pt0, double phi0, double eta0,double m1, double pt1, double eta1, double phi1)
{

  return TMath::Sqrt(m0*m0 + m1*m1 + 2.0 * (TMath::Sqrt((m0*m0+pt0*pt0)*(m1*m1+pt1*pt1)) - pt0*pt1*TMath::Cos(deltaPhi(phi0, phi1)) ));

  //TLorentzVector v1,v2;                                                                                                
  // v1.SetPtEtaPhiM(pt0,eta0, phi0, m0);                                                                                
  //v2.SetPtEtaPhiM(pt1,eta1, phi1, m1);                                                                                 
  // return (v1+v2).Mt();                                                                                                


}


double ptWeightQCD(int nGenVbosons, double lheHT, int GenVbosons_pdgId){
  double SF = 1.;
  if (lheHT>100 && nGenVbosons==1){
    if (GenVbosons_pdgId == 23){ // Z
      SF =   ((lheHT>100 && lheHT<200)*1.588 * ( 280.35 / (409.860000) ) + (lheHT>200 && lheHT<400)*1.438 * ( 77.67 / ( 110.880000 )) + (lheHT>400 && lheHT<600)*1.494 * (10.73 / (13.189 )) + (lheHT>600)*1.139 * ( 4.116 / (4.524300) ));
    }
    if (abs(GenVbosons_pdgId) == 24){
      SF =   ((lheHT>100 && lheHT<200)*1.588 * ( 1345 / (1.23 *  1.29e3) ) + (lheHT>200 && lheHT<400)*1.438 * ( 359.7 / ( 1.23 *  3.86e2)) + (lheHT>400 && lheHT<600)*1.494 * (48.91 / (1.23 * 47.9 )) + (lheHT>600)*1.139 * ( 18.77 / (1.23 * 19.9) ));
    }
  }
  return SF>0?SF:0;
}

// weights correction for EWK NLO correction
double ptWeightEWK(int nGenVbosons,double GenVbosons_pt,int VtypeSim,int GenVbosons_pdgId){
  double SF = 1.;
  if (nGenVbosons ==1)
    {
      if (VtypeSim == 0 || VtypeSim == 1 || VtypeSim == 4 || VtypeSim == 5)
	{
	  if (GenVbosons_pdgId == 23)
	    {
	      //for Z options
	      if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.1808051+6.04146*(TMath::Power((GenVbosons_pt+759.098),-0.242556));
	    }
	}
      else if (VtypeSim == 2 || VtypeSim == 3)
	{
	  //for W options
	  if (GenVbosons_pdgId == 24 || GenVbosons_pdgId == -24)
	    {
	      if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.830041+7.93714*(TMath::Power((GenVbosons_pt+877.978),-0.213831));
	    }
	}
    }
  return SF>0?SF:0;
}
