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

std::vector<std::pair<struct jet, struct jet>> dijetfinder::match_dijets_smear(std::vector<struct jet> myrecojets, std::vector<struct jet> mytruthjets)
{
  std::vector<std::pair<struct jet, struct jet>> matched_dijets = {};
  for (auto tjet : mytruthjets)
    {
      auto truth_iter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id;});
      float truth_closest_jet = 1.0;
      bool already_matched = false;

      if (truth_iter != matched_dijets.end())
	{
	  already_matched = true;
	  truth_closest_jet = truth_iter->first.dR;
	}

      
      for (auto jet : myrecojets)
	{
	  auto reco_iter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == jet.id;});
	  float reco_closest_jet = 1.0;
	  bool reco_already_matched = false;

	  if (reco_iter != matched_dijets.end())
	    {
	      reco_already_matched = true;
	      reco_closest_jet = truth_iter->first.dR;
	    }
	  
	  float dR = fabs(getDR(jet, tjet));
	  //float dphi = fabs(getDPHI(jet.phi, tjet.phi));
	  if (dR < m_wide_dR_cut && (reco_closest_jet > dR || truth_closest_jet > dR))
	    {
	      if (already_matched)
		{
		  auto truth_iter2 = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id;});
		  matched_dijets.erase(truth_iter2);
		}
	      if (reco_already_matched)
		{
		  auto reco_iter2 = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == jet.id;});
		  matched_dijets.erase(reco_iter2);
		}

	      tjet.matched = 1;
	      tjet.dR = dR;
	      jet.matched = 1;
	      jet.dR = dR;
	      already_matched = true;
	      truth_closest_jet = dR;
	      matched_dijets.push_back(std::make_pair(tjet, jet));	
	      break;
	    }	  
	}
    }

  // collect the rest
  for (auto tjet : mytruthjets)
    {
      auto truth_iter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto b) { return b.first.id == tjet.id;});

      if (truth_iter != matched_dijets.end())
	{
	  continue;
	}

      
      for (auto jet : myrecojets)
	{
	  auto reco_iter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto b) { return b.second.id == jet.id;});

	  if (reco_iter != matched_dijets.end())
	    {
	      continue;
	    }
	  
	  float dR = fabs(getDR(jet, tjet));
	  //float dphi = fabs(getDPHI(jet.phi, tjet.phi));
	  if (dR < m_wide_dR_cut)
	    {
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
std::vector<std::pair<struct jet, struct jet>> dijetfinder::match_dijets_response(std::vector<struct jet> myrecojets, std::vector<struct jet> mytruthjets)
{
  std::vector<std::pair<struct jet, struct jet>> matched_dijets = {};
  for (auto tjet : mytruthjets)
    {
      auto truth_iter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id;});
      float truth_closest_jet = 1.0;
      bool already_matched = false;

      if (truth_iter != matched_dijets.end())
	{
	  already_matched = true;
	  truth_closest_jet = truth_iter->first.dR;
	}

      
      for (auto jet : myrecojets)
	{
	  auto reco_iter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == jet.id;});
	  float reco_closest_jet = 1.0;
	  bool reco_already_matched = false;

	  if (reco_iter != matched_dijets.end())
	    {
	      reco_already_matched = true;
	      reco_closest_jet = truth_iter->first.dR;
	    }
	  
	  float dR = fabs(getDR(jet, tjet));
	  float dphi = fabs(getDPHI(jet.phi, tjet.phi));
	  if (dR < m_dR_cut && (reco_closest_jet > dR || truth_closest_jet > dR))
	    {
	      if (already_matched)
		{
		  auto truth_iter2 = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id;});
		  matched_dijets.erase(truth_iter2);
		}
	      if (reco_already_matched)
		{
		  auto reco_iter2 = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == jet.id;});
		  matched_dijets.erase(reco_iter2);
		}

	      tjet.matched = 1;
	      tjet.dR = dR;
	      jet.matched = 1;
	      jet.dR = dR;
	      already_matched = true;
	      truth_closest_jet = dR;
	      matched_dijets.push_back(std::make_pair(tjet, jet));	
	      break;
	    }	  
	}
    }

  // collect the rest
  for (auto tjet : mytruthjets)
    {
      auto truth_iter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto b) { return b.first.id == tjet.id;});

      if (truth_iter != matched_dijets.end())
	{
	  continue;
	}

      
      for (auto jet : myrecojets)
	{
	  auto reco_iter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto b) { return b.second.id == jet.id;});

	  if (reco_iter != matched_dijets.end())
	    {
	      continue;
	    }
	  
	  float dR = fabs(getDR(jet, tjet));
	  float dphi = fabs(getDPHI(jet.phi, tjet.phi));
	  if (dR < m_dR_cut)
	    {
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

std::vector<std::pair<struct jet, struct jet>>  dijetfinder::match_dijets(std::vector<struct jet> myrecojets, std::vector<struct jet> mytruthjets)
{
  std::vector<std::pair<struct jet, struct jet>> matched_dijets = {};
  for (auto tjet : mytruthjets)
    {
      auto truth_iter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id;});
      float truth_closest_jet = 1.0;
      bool already_matched = false;

      if (truth_iter != matched_dijets.end())
	{
	  already_matched = true;
	  truth_closest_jet = truth_iter->first.dR;
	}
      
      for (auto jet : myrecojets)
	{
	  // auto reco_iter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == jet.id;});
	  //if (reco_iter != matched_dijets.end()) continue;
	  
	  float dR = fabs(getDR(jet, tjet));
	  float dphi = fabs(getDPHI(jet.phi, tjet.phi));
	  if (dphi < m_dR_cut && truth_closest_jet > dR)
	    {
	      if (already_matched)
		{
		  auto truth_iter2 = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id;});
		  matched_dijets.erase(truth_iter2);
		}
	      tjet.matched = 1;
	      tjet.dR = dR;
	      jet.matched = 1;
	      jet.dR = dR;
	      already_matched = true;
	      truth_closest_jet = dR;
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

  if (m_philoc)
    {
      float philead = fabs(leading_iter->phi);
      if (philead > TMath::Pi())
	{
	  philead = 2*TMath::Pi() - philead;
	}
      float phisublead = fabs(subleading_iter->phi);
      if (phisublead > TMath::Pi())
	{
	  phisublead = 2*TMath::Pi() - phisublead;
	}
      if (m_philoc == 1 && !(philead >= TMath::Pi()/4. && philead < 3*TMath::Pi()/4.)) return false;
      if (m_philoc == 2 && !(philead < TMath::Pi()/4. || philead >= 3*TMath::Pi()/4.)) return false;

      if (m_philoc == 1 && !(phisublead >= TMath::Pi()/4. && phisublead < 3*TMath::Pi()/4.)) return false;
      if (m_philoc == 2 && !(phisublead < TMath::Pi()/4. || phisublead >= 3*TMath::Pi()/4.)) return false;
    }
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
  //if (subleading_iter->e/leading_iter->e < 0.3) return false;
  // double jetdeltatime = 17.6*(leading_iter->t - subleading_iter->t);
  // double jetleadtime = 17.6*(leading_iter->t);
  // bool passleadtime = ( TMath::Abs(jetleadtime +2.0) < 6.0 );
  // bool passdijettime = (TMath::Abs(jetdeltatime) < 3.0);	  
  
  // bool passbothtime = (passdijettime) && (passleadtime);
  
  //if (!passbothtime) return false;

  if (m_verbosity > 5) std::cout << " LOOKS GOOD " << std::endl;
  
  return true;
  
}
bool dijetfinder::passes_time_cut(double calib_lead_time, double calib_delta_time)
{
  double x = calib_lead_time;
  double y = x - calib_delta_time;
  if (fabs(x) < 6 && fabs(y) < 3) return true;
  return false;
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

  if (m_philoc)
    {
      float philead = fabs(leading_iter->phi);
      if (philead > TMath::Pi())
	{
	  philead = 2*TMath::Pi() - philead;
	}
      if (m_philoc == 1 && !(philead >= TMath::Pi()/4. && philead < 3*TMath::Pi()/4.)) return false;
      if (m_philoc == 2 && !(philead < TMath::Pi()/4. || philead >= 3*TMath::Pi()/4.)) return false;
      float phisublead = fabs(subleading_iter->phi);
      if (phisublead > TMath::Pi())
	{
	  phisublead = 2*TMath::Pi() - phisublead;
	}
      if (m_philoc == 1 && !(phisublead >= TMath::Pi()/4. && phisublead < 3*TMath::Pi()/4.)) return false;
      if (m_philoc == 2 && !(phisublead < TMath::Pi()/4. || phisublead >= 3*TMath::Pi()/4.)) return false;
    }

  float dphir = fabs(getDPHI(leading_iter->phi, subleading_iter->phi));

  if (m_verbosity > 5)
    {
      std::cout << " Checking dijet w/ " << m_truth_leading_cut << " / " << m_truth_subleading_cut << " / " << m_dphicut << " --> " << leading_iter->pt << " / " << subleading_iter->pt << " / " << dphir << std::endl;
    }

  if (!(leading_iter->pt >= m_truth_leading_cut && subleading_iter->pt >= m_truth_subleading_cut && dphir >= m_dphicut)) return false;

  if (m_verbosity > 5) std::cout << " LOOKS GOOD " << std::endl;  
  return true;
}
