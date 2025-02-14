#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;



#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#endif
#include "../macros/dlUtility.h"

struct jet
{
  int id;
  float pt = 0;
  float eta = 0;
  float phi = 0;
  int matched = 0;
  float dR = 1;
};

std::vector<float> dividebin(TH2* h2, int bin, int ndivs);

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
void unfoldTruth()
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

  float max_pt = 100;
  float binwidth = 5;
  int nbins = max_pt/binwidth;
  int nbins_flat = nbins*nbins;
  std::vector<struct jet> mytruthjets;
  std::vector<struct jet> myrecojets;
  std::vector<std::pair<struct jet, struct jet>> matched_dijets;
  int entries = t->GetEntries();

  TH1D *h_reco_xj = new TH1D("h_reco_xj","", 10, 0, 1);
  TH1D *h_reco_dphi = new TH1D("h_reco_dphi","", 64, 0, TMath::Pi());
  TH1D *h_reco_leading_jet = new TH1D("h_reco_leading_jet","", 100, 0, 100);
  TH1D *h_reco_subleading_jet = new TH1D("h_reco_subleading_jet","", 100, 0, 100);
  TH1D *h_matched_reco_xj = new TH1D("h_matched_reco_xj","", 10, 0, 1);
  TH1D *h_matched_reco_dphi = new TH1D("h_matched_reco_dphi","", 64, 0, TMath::Pi());
  TH1D *h_matched_reco_leading_jet = new TH1D("h_matched_reco_leading_jet","", 100, 0, 100);
  TH1D *h_matched_reco_subleading_jet = new TH1D("h_matched_reco_subleading_jet","", 100, 0, 100);
  TH1D *h_matched_truth_xj = new TH1D("h_matched_truth_xj","", 10, 0, 1);
  TH1D *h_matched_truth_dphi = new TH1D("h_matched_truth_dphi","", 64, 0, TMath::Pi());
  TH1D *h_matched_truth_leading_jet = new TH1D("h_matched_truth_leading_jet","", 100, 0, 100);
  TH1D *h_matched_truth_subleading_jet = new TH1D("h_matched_truth_subleading_jet","", 100, 0, 100);
  TH1D *h_matched_dR_leading = new TH1D("h_matched_dR_leading","", 200, 0, 1.0);
  TH1D *h_matched_dR_subleading = new TH1D("h_matched_dR_subleading","", 200, 0, 1.0);
  TH2D *h_matched_dR_leading_truth_pt = new TH2D("h_matched_dR_leading_truth_pt","", 200, 0, 1.0, 100, 0, 100);
  TH2D *h_matched_dR_subleading_truth_pt = new TH2D("h_matched_dR_subleading_truth_pt","", 200, 0, 1.0, 100, 0, 100);

  RooUnfoldResponse rooResponse((int) nbins_flat, -0.5, nbins_flat - 0.5); 

  TH1D *h_flat_matched_truth = new TH1D("h_matched_truth_flat",";p_{T,1, smear} + p_{T,2, smear}", (int) nbins_flat, -0.5, nbins_flat - 0.5); 
  TH1D *h_flat_matched_reco = new TH1D("h_matched_reco_flat",";p_{T,1, smear} + p_{T,2, smear}", (int) nbins_flat, -0.5, nbins_flat - 0.5); //100, -0.5, 99.5);
  TH1D *h_flat_reco = new TH1D("h_reco_flat",";p_{T,1, smear} + p_{T,2, smear}", (int) nbins_flat, -0.5, nbins_flat - 0.5); 

  TH1D *h_matched_leading_scale = new TH1D("h_matched_leading_scale",";p_{T,reco}/p_{T,truth};", 100, 0, 2);
  TH1D *h_matched_subleading_scale = new TH1D("h_matched_subleading_scale",";p_{T,reco}/p_{T,truth};", 100, 0, 2);

  TH2D *h2_matched_leading_scale = new TH2D("h2_matched_leading_scale",";p_{T,truth};;p_{T,reco}/p_{T,truth};", 15, 20, 50, 100, 0, 2);
  TH2D *h2_matched_subleading_scale = new TH2D("h2_matched_subleading_scale",";p_{T,truth};;p_{T,reco}/p_{T,truth};", 15, 20, 50, 100, 0, 2);

  float dRcut = 1.0;
  float leading_cut = 20;
  float subleading_cut = 8;
  float ptcut = subleading_cut;
  float dphicut = 7*TMath::Pi()/8.;
  float etacut = 0.7;
  int ndijets =0 ;
  for (int i = 0; i < entries; i++)
    {
      t->GetEntry(i);

      std::cout << i << " \r" << std::flush;
      int ntruthjets = truth_jet_pt_4->size();
      int nrecojets = reco_jet_pt_4->size();
      std::pair<int,int> dijet_ids = std::make_pair(-1,-1);
      matched_dijets.clear();
      mytruthjets.clear();
      myrecojets.clear();
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
      
      if (myrecojets.size() < 2) continue;

      std::sort(myrecojets.begin(), myrecojets.end(), [] (auto a, auto b) { return a.pt > b.pt; });

      auto leading_iter = myrecojets.begin();

      if (leading_iter->pt <= leading_cut) continue;

      auto subleading_iter = std::find_if(leading_iter, myrecojets.end(), [=] (auto a) {
	  double dphi = fabs(leading_iter->phi - a.phi);
	  if (dphi > TMath::Pi())
	    {
	      dphi = 2*TMath::Pi() - dphi;
	    }

	  if (dphi < dphicut) return false;
	  return true;
	});
      
      if (subleading_iter == myrecojets.end()) continue;
      float xj = (subleading_iter->pt/leading_iter->pt);
      float aj = (leading_iter->pt - subleading_iter->pt)/(leading_iter->pt + subleading_iter->pt);
      double dphir = fabs(leading_iter->phi - subleading_iter->phi);
      if (dphir > TMath::Pi())
	{
	  dphir = 2*TMath::Pi() - dphir;
	}

      dijet_ids = std::make_pair(leading_iter->id, subleading_iter->id);

      int pt1bin = leading_iter->pt/binwidth;
      int pt2bin = subleading_iter->pt/binwidth;

      for (int j = 0; j < ntruthjets;j++)
	{
	  struct jet tempjet;
	  tempjet.pt = truth_jet_pt_4->at(j);
	  tempjet.eta = truth_jet_eta_4->at(j);
	  tempjet.phi = truth_jet_phi_4->at(j);
	  tempjet.id = j;
	  mytruthjets.push_back(tempjet);	  
	}

      if (mytruthjets.size() < 2) continue;
      std::sort(mytruthjets.begin(), mytruthjets.end(), [] (auto a, auto b) { return a.pt > b.pt; });


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

      int index_leading = -1;
      int index_subleading = -1;

      for (int id = 0; id < matched_dijets.size(); id++)
	{
	  // std::cout << " Matched Jet: " << std::endl;
	  // std::cout << "    Reco: " << matched_dijets.at(id).first.id << "  " << matched_dijets.at(id).first.pt << " / " <<matched_dijets.at(id).first.eta << " / " <<  matched_dijets.at(id).first.phi << std::endl;
	  // std::cout << "    Truth: " << matched_dijets.at(id).second.id << "  " << matched_dijets.at(id).second.pt << " / " <<matched_dijets.at(id).second.eta << " / " <<  matched_dijets.at(id).second.phi << std::endl;
	  // std::cout << "    Match dR = " << matched_dijets.at(id).second.dR << std::endl;
	  if (matched_dijets.at(id).first.id == dijet_ids.first)
	    {
	      index_leading = id;
	    }
	  if (matched_dijets.at(id).first.id == dijet_ids.second)
	    {
	      index_subleading = id;
	    }
	}

      if (index_leading >= 0 && index_subleading >= 0)
	{
s
	  auto truth_leading = matched_dijets.at(index_leading).first;
	  auto truth_subleading = matched_dijets.at(index_subleading).first;
	  
	  h_matched_dR_leading->Fill(matched_dijets.at(index_leading).second.dR);
	  h_matched_dR_subleading->Fill(matched_dijets.at(index_subleading).second.dR);
	  float xjt = (min(truth_subleading.pt, truth_leading.pt)) / max(truth_subleading.pt, truth_leading.pt); 
	  float ajt = (truth_leading.pt - truth_subleading.pt) / (truth_leading.pt + truth_subleading.pt) ;
	  double dphit = fabs(truth_leading.phi - truth_subleading.phi);
	  if (dphit > TMath::Pi())
	    {
	      dphit = 2*TMath::Pi() - dphit;
	    }

	  ndijets++;	  
	  h_matched_reco_xj->Fill(xj);
	  h_matched_reco_dphi->Fill(dphir);
	  h_matched_truth_xj->Fill(xjt);
	  h_matched_truth_dphi->Fill(dphit);

	  int pt1_reco = matched_dijets.at(index_leading).second.pt/binwidth;
	  int pt2_reco = matched_dijets.at(index_subleading).second.pt/binwidth;
	  int pt1_truth = truth_leading.pt/binwidth;
	  int pt2_truth = truth_subleading.pt/binwidth;

	  h2_matched_leading_scale->Fill(truth_leading.pt, matched_dijets.at(index_leading).second.pt/truth_leading.pt);
	  h2_matched_subleading_scale->Fill(truth_subleading.pt, matched_dijets.at(index_subleading).second.pt/truth_subleading.pt);

	  h_matched_leading_scale->Fill(matched_dijets.at(index_leading).second.pt/truth_leading.pt);
	  h_matched_subleading_scale->Fill(matched_dijets.at(index_subleading).second.pt/truth_subleading.pt);

	  h_flat_matched_truth->Fill(pt1_truth*nbins + pt2_truth);
	  h_flat_matched_truth->Fill(pt2_truth*nbins + pt1_truth);
	  h_flat_matched_reco->Fill(pt1_reco*nbins + pt2_reco);
	  h_flat_matched_reco->Fill(pt2_reco*nbins + pt1_reco);

	  rooResponse.Fill(pt1_reco*nbins + pt2_reco, pt1_truth*nbins + pt2_truth);
	  rooResponse.Fill(pt2_reco*nbins + pt1_reco, pt2_truth*nbins + pt1_truth);
	}
      
    }
  
  std::cout << "Number of Dijets: " << ndijets << std::endl;
  TCanvas *c3 = new TCanvas("c3","c3", 500, 500);
  h_matched_dR_leading->SetLineColor(kRed);
  h_matched_dR_subleading->SetLineColor(kBlue);
  h_matched_dR_leading->Draw();
  h_matched_dR_subleading->Draw("same");
  h_flat_matched_reco->Scale(0.5);
  h_flat_matched_truth->Scale(0.5);  
  RooUnfoldBayes   unfold (&rooResponse, h_flat_matched_reco, 4);    // OR

  TH1D* hUnfold= (TH1D*) unfold.Hunfold();

  TCanvas *c4 = new TCanvas("c4","c4", 500, 500);
  h_flat_matched_truth->SetLineColor(kRed);
  h_flat_matched_reco->SetLineColor(kBlue);
  hUnfold->SetMarkerColor(kBlack);
  hUnfold->SetMarkerSize(1);
  hUnfold->SetMarkerStyle(8);
  h_flat_matched_truth->Draw("hist");
  h_flat_matched_reco->Draw("same hist");
  hUnfold->Draw("same p");

  TCanvas *c5 = new TCanvas("c5","c5", 800, 800);
  c5->Divide(2,2);
  c5->cd(1);
  h_matched_leading_scale->Draw();
  c5->cd(2);
  h2_matched_leading_scale->Draw("colz");
  c5->cd(3);
  h_matched_subleading_scale->Draw();
  c5->cd(4);
  h2_matched_subleading_scale->Draw("colz");

  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, 0, max_pt, nbins, 0, max_pt);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, 0, max_pt, nbins, 0, max_pt);
  TH2D *h_pt1pt2_unfold = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}", nbins, 0, max_pt, nbins, 0, max_pt);

  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", 10, 0, 1);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};", 10, 0, 1);
  TH1D *h_xj_unfold = new TH1D("h_xj_unfold", ";x_{J};", 10, 0, 1);

  for (int ib = 0; ib < h_flat_matched_reco->GetNbinsX(); ib++)
    {
      int xbin = ib/nbins;
      int ybin = ib%nbins;
      
      int b = h_pt1pt2_reco->GetBin(xbin+1, ybin+1);

      h_pt1pt2_reco->SetBinContent(b, h_flat_matched_reco->GetBinContent(ib+1));
      h_pt1pt2_truth->SetBinContent(b, h_flat_matched_truth->GetBinContent(ib+1));
      h_pt1pt2_unfold->SetBinContent(b, hUnfold->GetBinContent(ib+1));
    }
  for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy < nbins; iy++)
	{
	  int bin = h_pt1pt2_reco->GetBin(ix+1, iy+1);

	  if (ix > iy)
	    {

	      h_pt1pt2_reco->SetBinContent(bin, h_pt1pt2_reco->GetBinContent(bin)*2.);
	      h_pt1pt2_truth->SetBinContent(bin, h_pt1pt2_truth->GetBinContent(bin)*2.);
	      h_pt1pt2_unfold->SetBinContent(bin, h_pt1pt2_unfold->GetBinContent(bin)*2.);
	    }
	  else if (ix < iy)
	    {	      
	      h_pt1pt2_reco->SetBinContent(bin, 0);
	      h_pt1pt2_truth->SetBinContent(bin, 0);
	      h_pt1pt2_unfold->SetBinContent(bin, 0);
	    }

	}
    }
  for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy <= ix; iy++)
	{
	  int bin = h_pt1pt2_reco->GetBin(ix+1, iy+1);
  
	  std::vector<float> divisions_reco = dividebin(h_pt1pt2_reco, bin, 10);
	  std::vector<float> divisions_truth = dividebin(h_pt1pt2_truth, bin, 10);
	  std::vector<float> divisions_unfold = dividebin(h_pt1pt2_unfold, bin, 10);
	  std::cout << divisions_reco.size() << " -- " << ix << " -- " << iy << std::endl;
	  float sumreco = 0;
	  for (int id = 0; id < divisions_reco.size(); id++)
	    {
	      sumreco += h_pt1pt2_reco->GetBinContent(bin)*divisions_reco.at(id);
	      if (id == 0) continue;
	      
	      h_xj_reco->Fill((id - 1) * 0.1 + 0.05, h_pt1pt2_reco->GetBinContent(bin)*divisions_reco.at(id));
	    }
	  if (sumreco)
	    std::cout << " == " << sumreco << " / " << h_pt1pt2_reco->GetBinContent(bin) << " = " << sumreco/h_pt1pt2_reco->GetBinContent(bin)<<std::endl;

	  for (int id = 0; id < divisions_truth.size(); id++)
	    {
	      if (id == 0) continue;
	      h_xj_truth->Fill((id - 1) *0.1 + 0.05, h_pt1pt2_truth->GetBinContent(bin)*divisions_truth.at(id));
	    }
	  for (int id = 0; id < divisions_unfold.size(); id++)
	    {
	      if (id == 0) continue;
	      h_xj_unfold->Fill((id - 1) *0.1 + 0.05, h_pt1pt2_unfold->GetBinContent(bin)*divisions_unfold.at(id));
	    }
	}
    }

  std::cout << "RECO FLAT  : " <<  h_flat_matched_reco->Integral() << std::endl;
  std::cout << "RECO PT1PT2: " <<  h_pt1pt2_reco->Integral() << std::endl;
  std::cout << "RECO XJ PRJ: " <<  h_xj_reco->Integral() << std::endl;
  std::cout << "RECO XJ    : " <<  h_matched_reco_xj->Integral() << std::endl;

  std::cout << "TRUTH FLAT  : " <<  h_flat_matched_truth->Integral() << std::endl;
  std::cout << "TRUTH PT1PT2: " <<  h_pt1pt2_truth->Integral() << std::endl;
  std::cout << "TRUTH XJ PRJ: " <<  h_xj_truth->Integral() << std::endl;
  std::cout << "TRUTH XJ    : " <<  h_matched_truth_xj->Integral() << std::endl;

  std::cout << "UNFOLD PT1PT2: " <<  h_pt1pt2_unfold->Integral() << std::endl;
  std::cout << "UNFOLD XJ PRJ: " <<  h_xj_unfold->Integral() << std::endl;
	      
  TCanvas *c6 = new TCanvas("c6","c6", 700 ,300);
  c6->Divide(3,1);
  c6->cd(1);
  h_pt1pt2_reco->Draw("colz");
  c6->cd(2);
  h_pt1pt2_truth->Draw("colz");
  c6->cd(3);
  h_pt1pt2_unfold->Draw("colz");
  TCanvas *cunfold = new TCanvas("cunfold","cunfold", 500, 500);
  h_xj_reco->Draw();
  h_xj_unfold->SetLineColor(kBlue);
  h_xj_unfold->Draw("same");
  TCanvas *creco = new TCanvas("creco","creco",500, 700);
  dlutility::ratioPanelCanvas(creco);
  creco->cd(1);
  dlutility::SetLineAtt(h_xj_reco, kBlack,1 ,1 );
  dlutility::SetLineAtt(h_matched_reco_xj, kRed,1 ,1 );
  h_xj_reco->Draw();
  h_matched_reco_xj->Draw("same");
  creco->cd(2);
  TH1D *hd_reco = (TH1D*) h_xj_reco->Clone();
  hd_reco->Divide(h_matched_reco_xj);
  hd_reco->Draw();
  TCanvas *ctruth = new TCanvas("ctruth","ctruth",500, 700);
  dlutility::ratioPanelCanvas(ctruth);
  ctruth->cd(1);
  dlutility::SetLineAtt(h_xj_truth, kBlack,1 ,1 );
  dlutility::SetLineAtt(h_matched_truth_xj, kRed,1 ,1 );
  h_xj_truth->Draw();
  h_matched_truth_xj->Draw("same");
  ctruth->cd(2);
  TH1D *hd_truth = (TH1D*) h_xj_truth->Clone();
  hd_truth->Divide(h_matched_truth_xj);
  hd_truth->Draw();
  

  
  // TH1D *h_truth = new TH1D("h_truth","", 10, 0, 1);
  // TH1D *h_reco = new TH1D("h_reco","", 10, 0, 1);
  // TH1D *h_unfold = new TH1D("h_unfold","", 10, 0, 1);

  // for (int i = 0; i < 100; i++)
  //   {
  //     float v = (i%10)*0.1 + 0.05;
  //     h_unfold->Fill(v, hUnfold->GetBinContent(i+1));
  //     h_matched_unfold->Fill(v, hUnfoldm->GetBinContent(i+1));
  //     h_truth->Fill(v, h_flat_truth_->GetBinContent(i+1));
  //     h_matched_truth->Fill(v, h_flat_matched_truth_->GetBinContent(i+1));
  //     h_matched_reco->Fill(v, h_flat_matched_reco_->GetBinContent(i+1));
  //   }

  // for (int i = 0; i < 10; i++)
  //   {

  //     h_unfold->SetBinError(i+1, sqrt(h_unfold->GetBinContent(i+1)));
  //     h_matched_unfold->SetBinError(i+1, sqrt(h_matched_unfold->GetBinContent(i+1)));
  //     h_truth->SetBinError(i+1, sqrt(h_truth->GetBinContent(i+1)));
  //     h_matched_truth->SetBinError(i+1, sqrt(h_matched_truth->GetBinContent(i+1)));
  //     h_matched_reco->SetBinError(i+1, sqrt(h_matched_reco->GetBinContent(i+1)));
  //   }

  // dlutility::SetLineAtt(h_matched_reco, kBlue, 1, 1);
  // dlutility::SetLineAtt(h_matched_truth, kRed, 1, 1);
  // dlutility::SetLineAtt(h_truth, kBlack, 1, 1);
  // dlutility::SetLineAtt(h_matched_unfold, kRed, 1, 1);
  // dlutility::SetLineAtt(h_unfold, kBlack, 1, 1);
  // dlutility::SetMarkerAtt(h_matched_unfold, kRed, 1, 8);
  // dlutility::SetMarkerAtt(h_unfold, kBlack, 1, 8);

  // h_truth->Scale(1./h_truth->Integral(0, -1, "width"));

  // h_unfold->Scale(1./h_unfold->Integral(0, -1, "width"));
  // h_matched_truth->Scale(1./h_matched_truth->Integral(0, -1, "width"));
  // h_matched_reco->Scale(1./h_matched_reco->Integral(0, -1, "width"));
  // h_matched_unfold->Scale(1./h_matched_unfold->Integral(0, -1, "width"));

  // TCanvas *c5 = new TCanvas("c5","c5", 500, 500);
  // h_truth->Draw();
  // h_matched_reco->Draw("same");
  // h_unfold->Draw("same");

  // TCanvas *c6 = new TCanvas("c6","c6", 500, 500);
  // h_matched_truth->Draw();
  // h_matched_reco->Draw("same");
  // h_matched_unfold->Draw("same");
  return;
}
std::vector<float> dividebin(TH2* h2, int bin, int ndivs)
{
  std::vector<float> divisions;
  int binx, biny, binz;
  h2->GetBinXYZ(bin, binx, biny, binz);
  float x1 = h2->GetXaxis()->GetBinLowEdge(binx);
  float y1 = h2->GetYaxis()->GetBinLowEdge(biny);
  float x2 = x1 +  h2->GetXaxis()->GetBinWidth(binx);
  float y2 = y1 +  h2->GetYaxis()->GetBinWidth(biny);
  float total_area = (x2 - x1) * (y2 - y1);
  if (x1 == y1)
    {
      total_area /= 2.;
    }
  // std::cout << x1 << " , " << y1 << " -->  " << x2 << " , " << y2 << " = " << total_area << std::endl;
  float area_total = 0;
  float dx = 1./(float)ndivs;
  // std::cout << dx << std::endl;
  for (int i = 0; i <= ndivs; i++)
    {
      float slope = ((float)i)*dx;
      // std::cout << "slope =  " << slope << std::endl;
      float yb1 = slope*x1;
      float yb2 = slope*x2;
      // std::cout << "YBin : " << yb1 << " -- " << yb2 << std::endl;
      if (yb1 <= y1 && yb2 <= y1)
	{
	  divisions.push_back(0);
	  continue;
	}
      if (yb1 >= y2 && yb2 > y2)
	{
	  divisions.push_back((total_area - area_total)/total_area);
	  area_total = 25;
	  continue;
	}

      float xb1 = y1/slope;
      float xb2 = y2/slope;

      // std::cout << "XBin : " << xb1 << " -- " << xb2 << std::endl;

      if (xb1 < x1) xb1 = x1;
      if (xb2 > x2) xb2 = x2;

      if (yb1 < y1) yb1 = y1;
      if (yb2 > y2) yb2 = y2;
      // std::cout << "YBin : " << yb1 << " -- " << yb2 << std::endl;
      // std::cout << "XBin : " << xb1 << " -- " << xb2 << std::endl;
      

      float area = 0.5*(xb2 - xb1)*(yb2 - yb1);
      if (yb1 > y1)
	{
	  area += (yb1 - y1)*(x2 - x1);
	}
      if (xb2 < x2)
	{
	  area += (y2 - yb1)*(x2 - xb2);
	}

      area -= area_total;
      area_total += area;
      
      // std::cout << area << std::endl;
      divisions.push_back(area/total_area);
    }
  //  divisions.push_back(total_area - area_total);
  return divisions;  
}
