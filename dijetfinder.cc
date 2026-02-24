#include "dijetfinder.h"

float dijetfinder::getDPHI(float phi1, float phi2)
{
  float dphitr = phi1 - phi2;
  if (dphitr > TMath::Pi())
    {
      dphitr -= 2*TMath::Pi();
    }
  if (dphitr < -1*TMath::Pi())
    {
      dphitr += 2*TMath::Pi();
    }
  return dphitr;
  
}
float dijetfinder::getDR(struct jet j1, struct jet j2)
{ 
  double dphi = getDPHI(j1.phi,j2.phi);

  double dR = sqrt(TMath::Power(j1.eta - j2.eta, 2) + TMath::Power(dphi, 2));
  return dR;
}

std::vector<std::pair<struct jet, struct jet>>  dijetfinder::match_dijets(std::vector<struct jet> myrecojets, std::vector<struct jet> mytruthjets)
{
  std::vector<std::pair<struct jet, struct jet>> matched_dijets = {};
  for (auto jet : myrecojets)
    {
      for (auto tjet : mytruthjets)
	{
	  //if (std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == jet.id;}) != matched_dijets.end()) continue;
	  float dR = fabs(getDPHI(jet.phi, tjet.phi));
	      
	  if (dR < m_dR_cut)
	    {
	      if (std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id;}) != matched_dijets.end()) continue;


	      tjet.matched = 1;
	      tjet.dR = dR;
	      jet.matched = 1;
	      jet.dR = dR;

	      matched_dijets.push_back(std::make_pair(tjet, jet));	
	      break;
	    }	  
	}
    }

  return matched_dijets;
}

bool dijetfinder::check_dijet_reco(std::vector<struct jet> myrecojets)
{
  if (m_verbosity > 5)
    {
      std::cout << " Checking Reco Jets " << std::endl;
    }

  if (myrecojets.size() < 2) return false;
  auto leading_iter = myrecojets.begin();
  
  auto subleading_iter = myrecojets.begin() + 1;

  if (m_verbosity > 5)
    {
      std::cout << " Checking eta w/ " << m_eta_cut << std::endl;
    }

  if (fabs(leading_iter->eta_det) > m_eta_cut || fabs(subleading_iter->eta_det) > m_eta_cut) return false;
  if (fabs(leading_iter->eta) > m_eta_cut || fabs(subleading_iter->eta) > m_eta_cut) return false;

  if (m_verbosity > 5)
    {
      std::cout << " Checking eta_det w/ " << m_eta_cut << std::endl;
    }

  //  if (fabs(leading_iter->eta_det) > m_eta_cut || fabs(subleading_iter->eta_det) > m_eta_cut) return false;	  

  float dphir = fabs(getDPHI(leading_iter->phi, subleading_iter->phi));

  if (m_verbosity > 5)
    {
      std::cout << " Checking dijet w/ " << m_reco_leading_cut << " / " << m_reco_subleading_cut << " / " << m_dphicut << " --> " << leading_iter->pt << " / " << subleading_iter->pt << " / " << dphir << std::endl;
    }

  if (!(leading_iter->pt >= m_reco_leading_cut && subleading_iter->pt >= m_reco_subleading_cut && dphir >= m_dphicut)) return false;

  double jetdeltatime = 17.6*(leading_iter->t - subleading_iter->t);
  double jetleadtime = 17.6*(leading_iter->t);
  bool passleadtime = ( TMath::Abs(jetleadtime +2.0) < 6.0 );
  bool passdijettime = (TMath::Abs(jetdeltatime) < 3.0);	  
  
  bool passbothtime = (passdijettime) && (passleadtime);
  
  if (!passbothtime) return false;

  if (m_verbosity > 5) std::cout << " LOOKS GOOD " << std::endl;
  
  return true;
  
}
bool dijetfinder::check_dijet_truth(std::vector<struct jet> mytruthjets)
{
  if (m_verbosity > 5)
    {
      std::cout << " Checking Truth Jets " << std::endl;
    }

  if (mytruthjets.size() < 2) return false;
    
  auto leading_iter = mytruthjets.begin();
  
  auto subleading_iter = mytruthjets.begin() + 1;

  if (m_verbosity > 5)
    {
      std::cout << " Checking eta w/ " << m_eta_cut << std::endl;
    }


  if (fabs(leading_iter->eta) > m_eta_cut || fabs(subleading_iter->eta) > m_eta_cut) return false;
	  
  float dphir = fabs(getDPHI(leading_iter->phi, subleading_iter->phi));

  if (m_verbosity > 5)
    {
      std::cout << " Checking dijet w/ " << m_truth_leading_cut << " / " << m_truth_subleading_cut << " / " << m_dphicut << " --> " << leading_iter->pt << " / " << subleading_iter->pt << " / " << dphir << std::endl;
    }

  if (!(leading_iter->pt >= m_truth_leading_cut && subleading_iter->pt >= m_truth_subleading_cut && dphir >= m_dphicut)) return false;

  if (m_verbosity > 5) std::cout << " LOOKS GOOD " << std::endl;  
  return true;
}
