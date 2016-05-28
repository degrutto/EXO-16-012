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




