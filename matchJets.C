struct jet
{
  int id;
  float pt = 0;
  float eta = 0;
  float phi = 0;
  int matched = 0;
  int dR = 1;
};

double getDR(struct jet j1, struct jet j2)
{ 
  double dphi = fabs(j1.phi - j2.phi);
  if (dphi > TMath::Pi())
    {
      dphi = 2*TMath::Pi() - dphi;
    }

  double dR = sqrt(TMath::Power(j1.eta - j2.eta, 2) + TMath::Power(dphi, 2));
  return dR;
}
void matchJets()
{

  TFile *f = new TFile("sim/TREE_JET_v1_2_new_ProdA_2024-type12-00000022.root","r");

  TTree *t = (TTree*) f->Get("ttree");

  std::vector<float> *truth_jet_pt_4 = 0;
  std::vector<float> *truth_jet_eta_4 = 0;
  std::vector<float> *truth_jet_phi_4 = 0;
  std::vector<float> *reco_jet_pt_4 = 0;
  std::vector<float> *reco_jet_eta_4 = 0;
  std::vector<float> *reco_jet_phi_4 = 0;

  t->SetBranchAddress("truth_jet_pt_4", &truth_jet_pt_4);
  t->SetBranchAddress("truth_jet_eta_4", &truth_jet_eta_4);
  t->SetBranchAddress("truth_jet_phi_4", &truth_jet_phi_4);
  t->SetBranchAddress("jet_pt_4", &reco_jet_pt_4);
  t->SetBranchAddress("jet_eta_4", &reco_jet_eta_4);
  t->SetBranchAddress("jet_phi_4", &reco_jet_phi_4);


  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> myrecojets;
  std::vector<std::pair<struct jet, struct jet>> matched_dijets;
  int entries = t->GetEntries();

  TH1D *h_truth_aj = new TH1D("h_truth_aj","", 10, 0, 1);
  TH1D *h_truth_dphi = new TH1D("h_truth_dphi","", 64, 0, TMath::Pi());
  TH1D *h_truth_leading_jet = new TH1D("h_truth_leading_jet","", 100, 0, 100);
  TH1D *h_truth_subleading_jet = new TH1D("h_truth_subleading_jet","", 100, 0, 100);
  TH1D *h_matched_reco_aj = new TH1D("h_matched_reco_aj","", 10, 0, 1);
  TH1D *h_matched_reco_dphi = new TH1D("h_matched_reco_dphi","", 64, 0, TMath::Pi());
  TH1D *h_matched_reco_leading_jet = new TH1D("h_matched_reco_leading_jet","", 100, 0, 100);
  TH1D *h_matched_reco_subleading_jet = new TH1D("h_matched_reco_subleading_jet","", 100, 0, 100);
  TH1D *h_matched_truth_aj = new TH1D("h_matched_truth_aj","", 10, 0, 1);
  TH1D *h_matched_truth_dphi = new TH1D("h_matched_truth_dphi","", 64, 0, TMath::Pi());
  TH1D *h_matched_truth_leading_jet = new TH1D("h_matched_truth_leading_jet","", 100, 0, 100);
  TH1D *h_matched_truth_subleading_jet = new TH1D("h_matched_truth_subleading_jet","", 100, 0, 100);
  TH1D *h_matched_dR_leading = new TH1D("h_matched_dR_leading","", 200, 0, 1.0);
  TH1D *h_matched_dR_subleading = new TH1D("h_matched_dR_subleading","", 200, 0, 1.0);
  TH2D *h_matched_dR_leading_truth_pt = new TH2D("h_matched_dR_leading_truth_pt","", 200, 0, 1.0, 100, 0, 100);
  TH2D *h_matched_dR_subleading_truth_pt = new TH2D("h_matched_dR_subleading_truth_pt","", 200, 0, 1.0, 100, 0, 100);
  
  float dRcut = 1.0;
  float leading_cut = 20;
  float subleading_cut = 8;
  float ptcut = subleading_cut;
  float dphicut = 7*TMath::Pi()/8.;
  float etacut = 0.7;
  for (int i = 0; i < entries; i++)
    {
      t->GetEntry(i);

      std::cout << i << " \r" << std::flush;
      int ntruthjets = truth_jet_pt_4->size();
      int nrecojets = reco_jet_pt_4->size();

      matched_dijets.clear();
      mytruthjets.clear();
      myrecojets.clear();
      for (int j = 0; j < ntruthjets;j++)
	{
	  if (truth_jet_pt_4->at(j) > ptcut && fabs(truth_jet_eta_4->at(j)) < etacut)
	    {
	      struct jet tempjet;
	      tempjet.pt = truth_jet_pt_4->at(j);
	      tempjet.eta = truth_jet_eta_4->at(j);
	      tempjet.phi = truth_jet_phi_4->at(j);
	      tempjet.id = j;
	      mytruthjets.push_back(tempjet);
	    }
	}
      
      if (mytruthjets.size() < 2) continue;

      std::sort(mytruthjets.begin(), mytruthjets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

      auto leading_iter = mytruthjets.begin();

      if (leading_iter->pt <= leading_cut) continue;

      auto subleading_iter = std::find_if(leading_iter, mytruthjets.end(), [=] (auto a) {
	  double dphi = fabs(leading_iter->phi - a.phi);
	  if (dphi > TMath::Pi())
	    {
	      dphi = 2*TMath::Pi() - dphi;
	    }

	  if (dphi < dphicut) return false;
	  return true;
	});
      
      if (subleading_iter == mytruthjets.end()) continue;
      float aj = (leading_iter->pt - subleading_iter->pt)/(leading_iter->pt + subleading_iter->pt);
      double dphit = fabs(leading_iter->phi - subleading_iter->phi);
      if (dphit > TMath::Pi())
	{
	  dphit = 2*TMath::Pi() - dphit;
	}


      h_truth_aj->Fill(aj);
      h_truth_dphi->Fill(dphit);
      h_truth_leading_jet->Fill(leading_iter->pt);
      h_truth_subleading_jet->Fill(subleading_iter->pt);
      for (int j = 0; j < nrecojets;j++)
	{
	  if (reco_jet_pt_4->at(j) > ptcut && fabs(reco_jet_eta_4->at(j)) < etacut)
	    {
	      struct jet tempjet;
	      tempjet.pt = reco_jet_pt_4->at(j);
	      tempjet.eta = reco_jet_eta_4->at(j);
	      tempjet.phi = reco_jet_phi_4->at(j);
	      tempjet.id = j;
	      myrecojets.push_back(tempjet);
	    }
	}

      if (myrecojets.size() < 2)
	{
	  continue;
	}
      
      std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

      
      for (auto jet : myrecojets)
	{
	  for (auto tjet : mytruthjets)
	    {
	      float dR = getDR(jet, tjet);
	      
	      if (dR < dRcut && dR < jet.dR)
		{
		  jet.matched = 1;
		  jet.dR = dR;
		  
		  auto riter = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.second.id == jet.id; });
		  auto titer = std::find_if(matched_dijets.begin(), matched_dijets.end(), [=] (auto a) { return a.first.id == tjet.id; });

		  if (titer == matched_dijets.end() && riter == matched_dijets.end())
		    {
		      matched_dijets.push_back(std::make_pair(tjet, jet));
		    }
		  else if (titer == matched_dijets.end())
		    {
		      matched_dijets.erase(riter);
		      matched_dijets.push_back(std::make_pair(tjet, jet));
		    }
		  else if (riter == matched_dijets.end())
		    {
		      if (dR < titer->second.dR)
			{
			  matched_dijets.erase(titer);
			  matched_dijets.push_back(std::make_pair(tjet, jet));
			}
		    }
		  else
		    {
		      if (dR < titer->second.dR)
			{
			  matched_dijets.erase(titer);
			  matched_dijets.push_back(std::make_pair(tjet, jet));
			}
		      matched_dijets.erase(riter);
		    }

		    
		}	  
	    }
	}

      std::sort(matched_dijets.begin(), matched_dijets.end(), [] (auto a, auto b) { return a.second.pt > b.second.pt; });
      bool hit = false;
      double reco_aj = 2;
      double reco_dphi = -99;
      double leading_pt = 0;
      double subleading_pt = 0;
      double leading_dR = 0;
      double subleading_dR = 0;
      double leading_pt_truth = 0;
      double subleading_pt_truth = 0;
      if (matched_dijets.size() >= 2)
	{
	  struct jet leader = matched_dijets.at(0).second;
	  struct jet subleader = matched_dijets.at(1).second;

	  double dphi = fabs(leader.phi - subleader.phi);
	  if (dphi > TMath::Pi())
	    {
	      dphi = 2*TMath::Pi() - dphi;
	    }

	  if (leader.pt > leading_cut && subleader.pt > subleading_cut && dphi > dphicut)
	    {
	      hit = true;
	      reco_dphi = dphi;
	      reco_aj = (leader.pt - subleader.pt)/(leader.pt + subleader.pt);
	      leading_pt = leader.pt;
	      subleading_pt = subleader.pt;
	      leading_pt_truth = matched_dijets.at(0).first.pt;
	      subleading_pt_truth = matched_dijets.at(1).first.pt;
	      leading_dR = getDR(leader, matched_dijets.at(0).first);
	      subleading_dR = getDR(subleader, matched_dijets.at(1).first);
	    }

	}

      if (hit)
	{
	  std::cout << "Truth to Reco = "<< aj << " --> "<<reco_aj<<std::endl;
	  h_matched_truth_aj->Fill(aj);
	  h_matched_truth_dphi->Fill(dphit);
	  h_matched_truth_leading_jet->Fill(leading_iter->pt);
	  h_matched_truth_subleading_jet->Fill(subleading_iter->pt);
	  h_matched_reco_aj->Fill(reco_aj);
	  h_matched_reco_dphi->Fill(reco_dphi);
	  h_matched_reco_leading_jet->Fill(leading_pt);
	  h_matched_reco_subleading_jet->Fill(subleading_pt);
	  h_matched_dR_leading_truth_pt->Fill(leading_dR, leading_pt_truth);
	  h_matched_dR_subleading_truth_pt->Fill(subleading_dR, subleading_pt_truth);
	  h_matched_dR_leading->Fill(leading_dR);
	  h_matched_dR_subleading->Fill(subleading_dR);

	}

    }


    TCanvas *c = new TCanvas("c","c", 500, 500);
  h_matched_truth_aj->SetLineColor(kBlue);
  h_matched_reco_aj->SetLineColor(kRed);
  h_truth_aj->SetLineColor(kBlack);
  h_matched_truth_aj->Scale(1./h_matched_truth_aj->Integral(0, -1, "width"));
  h_matched_reco_aj->Scale(1./h_matched_reco_aj->Integral(0, -1, "width"));
  h_truth_aj->Scale(1./h_truth_aj->Integral(0, -1, "width"));

  h_matched_truth_aj->Draw();
  h_truth_aj->Draw("same");
  h_matched_reco_aj->Draw("same");

  TCanvas *c2 = new TCanvas("c2","c2", 500, 500);
  h_matched_truth_dphi->SetLineColor(kBlue);
  h_matched_reco_dphi->SetLineColor(kRed);
  h_truth_dphi->SetLineColor(kBlack);
  h_matched_truth_dphi->Scale(1./h_matched_truth_dphi->Integral(0, -1, "width"));
  h_matched_reco_dphi->Scale(1./h_matched_reco_dphi->Integral(0, -1, "width"));
  h_truth_dphi->Scale(1./h_truth_dphi->Integral(0, -1, "width"));

  h_matched_truth_dphi->Draw();
  h_truth_dphi->Draw("same");
  h_matched_reco_dphi->Draw("same");

  TCanvas *c3 = new TCanvas("c3","c3", 500, 500);
  h_matched_dR_leading->SetLineColor(kRed);
  h_matched_dR_subleading->SetLineColor(kBlue);
  h_matched_dR_leading->Draw();
  h_matched_dR_subleading->Draw("same");

  TCanvas *c1 = new TCanvas("c1","c1", 500, 500);
  c1->Divide(3, 1);
  c1->cd(1);
  gPad->SetLogy();
  h_truth_leading_jet->SetLineColor(kRed);
  h_truth_subleading_jet->SetLineColor(kBlue);
  h_truth_leading_jet->Draw();
  h_truth_subleading_jet->Draw("same");
  c1->cd(2);
  gPad->SetLogy();
  h_matched_truth_leading_jet->SetLineColor(kRed);
  h_matched_truth_subleading_jet->SetLineColor(kBlue);
  h_matched_truth_leading_jet->Draw();
  h_matched_truth_subleading_jet->Draw("same");
  c1->cd(3);
  gPad->SetLogy();
  h_matched_reco_leading_jet->SetLineColor(kRed);
  h_matched_reco_subleading_jet->SetLineColor(kBlue);
  h_matched_reco_leading_jet->Draw();
  h_matched_reco_subleading_jet->Draw("same");
  
  
  return;
}
