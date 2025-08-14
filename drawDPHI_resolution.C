#include "dlUtility.h"
#include "read_binning.h"

const int colorsys[5] = {kBlack, kMagenta + 2, kPink + 2, kGreen - 2, kAzure + 4};
const int markersys[4] = {8, 21, 22, 23};
void drawDPHI_resolution(const std::string configfile = "binning.config")
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb(configfile.c_str());

  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();
  Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;

  const int nbinsdphi  = 32;
  float min_dphi = 3*TMath::Pi()/4.;
  float max_dphi = TMath::Pi();
  float stepdphi = (max_dphi - min_dphi)/(float)nbinsdphi;

  float idphi_bins[nbinsdphi+1];
  for (int i = 0; i < nbinsdphi+1; i++)
    {
      idphi_bins[i] = min_dphi + i*stepdphi;
    }
  const int nbins = read_nbins;

  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
    }

  float truth_leading_cut = rb.get_truth_leading_cut();
  float truth_subleading_cut = rb.get_truth_subleading_cut();

  float reco_leading_cut = rb.get_reco_leading_cut();
  float reco_subleading_cut = rb.get_reco_subleading_cut();

  float measure_leading_cut = rb.get_measure_leading_cut();
  float measure_subleading_cut = rb.get_measure_subleading_cut();
    
  int truth_leading_bin = rb.get_truth_leading_bin();
  int truth_subleading_bin = rb.get_truth_subleading_bin();

  int reco_leading_bin = rb.get_reco_leading_bin();
  int reco_subleading_bin = rb.get_reco_subleading_bin();

  int measure_leading_bin = rb.get_measure_leading_bin();
  int measure_subleading_bin = rb.get_measure_subleading_bin();
  const int mbins = rb.get_measure_bins();
  float sample_boundary[4] = {0};
  int binranges[10] = {0};
  int binrangesmin[10] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ir = 0; ir < mbins+1; ir++)
    {
      binranges[ir] = rb.get_measure_region(ir);
      binrangesmin[ir] = rb.get_subleading_measure_region(ir);
      std::cout << binranges[ir] << " -- " <<  binrangesmin[ir] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;


  
  TString filename = "dphi_hists.root";
  TFile *fin = new TFile(filename, "r");
  TProfile2D *h2_sim_match_dphires = (TProfile2D*) fin->Get("h2_sim_match_ddphi");
  const int n2bins = 4;
  TH1D *h_pt1_resolution[n2bins];

  for (int i = 0; i < 4; i++)
    {
      h_pt1_resolution[i] = new TH1D(Form("h_pt1_resolution_%d", i), ";p_{T2} [GeV]; d#Delta#phi;",10, 10, 60);
      for (int j = 0; j < 10; j++)
	{
	  int xbin = i+1;
	  int ybin = j+1;
	  int bbin = h2_sim_match_dphires->GetBin(xbin, ybin);
	  h_pt1_resolution[i]->SetBinContent(ybin, h2_sim_match_dphires->GetBinError(bbin));
	}
    }

  for (int i = 0; i < 4; i++)
    {
      dlutility::SetLineAtt(h_pt1_resolution[i], colorsys[i+1], 1, 1);
      dlutility::SetMarkerAtt(h_pt1_resolution[i], colorsys[i+1], 1+0.1*i, markersys[i]);
    }

  TCanvas *c = new TCanvas("c","c", 500, 500);
  h_pt1_resolution[3]->SetMaximum(0.1);
  h_pt1_resolution[3]->SetMinimum(0.0);
  h_pt1_resolution[3]->Draw("p");
  h_pt1_resolution[2]->Draw("p same");
  h_pt1_resolution[1]->Draw("p same");
  h_pt1_resolution[0]->Draw("p same");


  TLegend *l = new TLegend(0.56, 0.6, 0.87, 0.78);
  l->SetLineWidth(0);
  l->AddEntry(h_pt1_resolution[0], "20 #geq p_{T1} < 30 GeV","p");
  l->AddEntry(h_pt1_resolution[1], "30 #geq p_{T1} < 40 GeV","p");
  l->AddEntry(h_pt1_resolution[2], "40 #geq p_{T1} < 50 GeV","p");
  l->AddEntry(h_pt1_resolution[3], "50 #geq p_{T1} < 60 GeV","p");
  l->Draw("same");

  
  dlutility::DrawSPHENIXpp(0.6, 0.87);
  TProfile *h_sim_match_ddphi = (TProfile*) fin->Get("h_sim_match_ddphi");

  TH1D *h_res = new TH1D("h_res",";<p_{T}>; d#Delta#phi;",25, 10, 60);
  for (int i = 0; i < 25; i++)
    {
      h_res->SetBinContent(i+1, h_sim_match_ddphi->GetBinError(i+1));
    }

  dlutility::SetLineAtt(h_res, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_res, kBlack, 1, 8);
  TCanvas *c1 = new TCanvas("c1","c1", 500, 500);
  h_res->SetTitle(";< p_{T} > [GeV]; #sigma(d#Delta#phi)");
  h_res->SetMaximum(0.17);
  h_res->SetMinimum(0.0);
  h_res->Draw("p");
  dlutility::DrawSPHENIXpp(0.6, 0.87);
  dlutility::drawText("PYTHIA-8", 0.6, 0.77);
  c1->Print("dphi_res.pdf");
}
