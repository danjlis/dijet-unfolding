#ifndef _DIJETFINDER_H_
#define _DIJETFINDER_H_

#include <iostream>
#include <string>

using std::cout;
using std::endl;

#include <vector>
#include "TMath.h"


struct jet
{
  int id;
  float pt = 0;
  float pt_uncalib = 0;
  float eta = 0;
  float eta_det = 0;
  float phi = 0;
  float t = 0;
  float emcal = 0;
  int matched = 0;
  float dR = 1;
  int istruth = 0;
  
  void print()
  {
    std::cout << "Jet " << id << std::endl;
    std::cout << "    pt/eta/phi : " << pt << " / " << eta << " / " << phi << std::endl;
    if (!istruth) std::cout << "    uncalib/em/det/t   : " << pt_uncalib << " / " << emcal << " / " << eta_det << " / " << t << std::endl;
  };
};

class dijetfinder
{
 public:
  dijetfinder(int cone_size)
    {
      m_cone_size = cone_size;
      m_eta_cut = 1.1 - ( ((float) cone_size ) * 0.1 );
      m_dR_cut = 0.75 * ( ((float) cone_size ) * 0.1 );
    }
  
  ~dijetfinder(){}

  void SetVerbosity( int v ) { m_verbosity = v; }
  
  float getDPHI(float phi1, float phi2);
  float getDR(struct jet j1, struct jet j2);
  std::vector<std::pair<struct jet, struct jet>>  match_dijets(std::vector<struct jet> myrecojets, std::vector<struct jet> mytruthjets);
  bool check_dijet_reco(std::vector<struct jet> myrecojets);
  bool check_dijet_truth(std::vector<struct jet> mytruthjets);

  void setTruthCuts(float pt1, float pt2)
  {
    m_truth_leading_cut = pt1;
    m_truth_subleading_cut = pt2;
  }
  void setRecoCuts(float pt1, float pt2)
  {
    m_reco_leading_cut = pt1;
    m_reco_subleading_cut = pt2;
  }
  
 private:
  int m_verbosity = 0;
  float m_eta_cut = 0.7;
  int m_cone_size = 4;
  float m_dR_cut = 3.;
  float m_dphicut = 3*TMath::Pi()/4.;
  float m_truth_leading_cut = 10;
  float m_truth_subleading_cut = 5;
  float m_reco_leading_cut = 10;
  float m_reco_subleading_cut = 5;
};
#endif
