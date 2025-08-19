#include <iostream>

using std::cout;
using std::endl;


#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

int makeMCDataCompare(const std::string configfile = "binning_AA.config")
{
  int color_data = kBlack;
  int color_sim = kRed;
  int cone_size = 3;
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb(configfile.c_str());


  std::string j10_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v14_10_new_ProdA_2024-00000030.root";
  std::string j20_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v14_20_new_ProdA_2024-00000030.root";
  std::string j30_file = rb.get_tntuple_location() + "/TREE_MATCH_r0" + std::to_string(cone_size) + "_v14_30_new_ProdA_2024-00000030.root";

  float maxpttruth[3];
  float pt1_truth[3];
  float pt2_truth[3];
  float dphi_truth[3];
  float pt1_reco[3];
  float pt2_reco[3];
  float nrecojets[3];
  float dphi_reco[3];
  float match[3];
  float mbd_vertex[3];
  float centrality[3];

  float n_events[3];
  float b_n_events = 0;
  
  TFile *fin[3];
  fin[0] = new TFile(j10_file.c_str(), "r");
  fin[1] = new TFile(j20_file.c_str(), "r");
  fin[2] = new TFile(j30_file.c_str(), "r");

 
  TNtuple *tn[3];
  for (int i = 0 ; i < 3; i++)
    {
      if (!fin[i])
	{
	  std::cout << "no file: " << i << std::endl;
	  return 0;

	}

      tn[i] = (TNtuple*) fin[i]->Get("tn_match");
      tn[i]->SetBranchAddress("maxpttruth", &maxpttruth[i]);
      tn[i]->SetBranchAddress("pt1_truth", &pt1_truth[i]);
      tn[i]->SetBranchAddress("pt2_truth", &pt2_truth[i]);
      tn[i]->SetBranchAddress("dphi_truth", &dphi_truth[i]);
      tn[i]->SetBranchAddress("pt1_reco", &pt1_reco[i]);
      tn[i]->SetBranchAddress("pt2_reco", &pt2_reco[i]);
      tn[i]->SetBranchAddress("dphi_reco", &dphi_reco[i]);
      tn[i]->SetBranchAddress("nrecojets", &nrecojets[i]);
      tn[i]->SetBranchAddress("matched", &match[i]);
      tn[i]->SetBranchAddress("mbd_vertex", &mbd_vertex[i]);
      tn[i]->SetBranchAddress("centrality", &centrality[i]);

      TNtuple *tn_stats = (TNtuple*) fin[i]->Get("tn_stats");
      tn_stats->SetBranchAddress("nevents", &b_n_events);
      tn_stats->GetEntry(0);
      n_events[i] = b_n_events;

    }
  float cs_10 = (2.889e-6);
  float cs_20 = 5.4067742e-8;
  float cs_30 = (2.505e-9);
  
  float scale_factor[3];
  scale_factor[0] = (n_events[2]/n_events[0]) * cs_10/cs_30;
  scale_factor[1] = (n_events[2]/n_events[1]) * cs_20/cs_30; 
  scale_factor[2] = 1;
  
  Int_t minentries = rb.get_minentries();
  Int_t read_nbins = rb.get_nbins();
  //Int_t primer = rb.get_primer();

  Double_t dphicut = rb.get_dphicut();
  Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;

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

  float sample_boundary[4] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  TH1D *h_truth_lead_sample[3];
  for (int i = 0; i < 3; i++)
    {
      h_truth_lead_sample[i] = new TH1D(Form("h_truth_lead_%d", i), " ; Leading Jet p_{T} [GeV]; counts", 100, 0, 100);
    }
  
  TH1D *h_centrality = new TH1D("h_centrality", ";Centrality; counts", 20, 0, 100);
  TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);

  for (int isample = 0; isample < 3; isample++)
    {
      std::cout << "Sample " << isample << std::endl;
      int entries0 = 0;
      int entries1 = tn[isample]->GetEntries()/2;    
      int entries2 = tn[isample]->GetEntries();

      for (int i = entries0; i < entries2; i++)
	{
	  tn[isample]->GetEntry(i);
	  int inrecojets = nrecojets[isample];
	  double event_scale = scale_factor[isample];
	  
	  float max_truth = 0;
	  float max_reco = 0;
	  float min_truth = 0;
	  float min_reco = 0;

	  if (pt1_truth[isample] >= pt2_truth[isample])
	    {
	      max_truth = pt1_truth[isample];
	      max_reco = pt1_reco[isample];
	      min_truth = pt2_truth[isample];
	      min_reco = pt2_reco[isample];
	    }
	  else
	    {
	      max_truth = pt2_truth[isample];
	      max_reco = pt2_reco[isample];
	      min_truth = pt1_truth[isample];
	      min_reco = pt1_reco[isample];
	    }
	    
	  float pt1_truth_bin = nbins;
	  float pt2_truth_bin = nbins;
	  float pt1_reco_bin = nbins;
	  float pt2_reco_bin = nbins;

	  float e1 = max_truth;
	  float e2 = min_truth;
	  float es1 = max_reco;
	  float es2 = min_reco;

	  h_truth_lead_sample[isample]->Fill(e1, event_scale);

	  if (maxpttruth[isample] < sample_boundary[isample] || maxpttruth[isample] >= sample_boundary[isample+1]) continue;


	  // double smear1 = fgaus->GetRandom();
	  // double smear2 = fgaus->GetRandom();

	  // es1 = es1 + smear1*e1;
	  // es2 = es2 + smear2*e2; 	  
	  
	  float maxi = std::max(es1, es2);
	  float mini = std::min(es1, es2);


	  for (int ib = 0; ib < nbins; ib++)
	    {
	      if ( e1 < ipt_bins[ib+1] && e1 >= ipt_bins[ib])
		{
		  pt1_truth_bin = ib;
		}
	      if ( e2 < ipt_bins[ib+1] && e2 >= ipt_bins[ib])
		{
		  pt2_truth_bin = ib;
		}
	      if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
		{
		  pt1_reco_bin = ib;
		}
	      if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
		{
		  pt2_reco_bin = ib;
		}
	    }
	  
	  bool truth_good = (e1 >= truth_leading_cut && e2 >= truth_subleading_cut && dphi_truth[isample] > dphicuttruth);

	  bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco[isample] > dphicut);

	  if (!truth_good && !reco_good) continue;

	  h_mbd_vertex->Fill(mbd_vertex[isample], event_scale);
	  h_centrality->Fill(centrality[isample], event_scale);	      
	}
    }


  std::string data_file = rb.get_tntuple_location() + "/TNTUPLE_DIJET_r0" + std::to_string(cone_size) + "_v10_1_492_2024p020_v007_gl10-all.root";
  TFile *find = new TFile(data_file.c_str(), "r");
  TNtuple *tnd  = (TNtuple*) find->Get("tn_dijet");;
  if (!tnd)
    {
      std::cout << " no data "<< std::endl;
    }
  float d_mbd_vertex;
  float d_pt1_reco;
  float d_pt2_reco;
  float d_dphi_reco;
  float d_nrecojets;
  float d_trigger;
  float d_centrality;

  tnd->SetBranchAddress("pt1_reco", &d_pt1_reco);
  tnd->SetBranchAddress("pt2_reco", &d_pt2_reco);
  tnd->SetBranchAddress("dphi_reco", &d_dphi_reco);
  tnd->SetBranchAddress("njets", &d_nrecojets);
  tnd->SetBranchAddress("centrality", &d_centrality);
  tnd->SetBranchAddress("trigger", &d_trigger);
  tnd->SetBranchAddress("mbd_vertex", &d_mbd_vertex);

  int d_entries = tnd->GetEntries();

  TH1D *h_d_mbd_vertex = new TH1D("h_d_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
  TH1D *h_d_centrality = new TH1D("h_d_centrality", ";Centrality; counts", 20, 0, 100);

  for (int i = 0; i < d_entries; i++)
    {
      tnd->GetEntry(i);


      float maxi = std::max(d_pt1_reco, d_pt2_reco);
      float mini = std::min(d_pt1_reco, d_pt2_reco);

      float pt1_reco_bin = nbins;
      float pt2_reco_bin = nbins;

      float es1 = d_pt1_reco;
      float es2 = d_pt2_reco;
      if (es1 > ipt_bins[nbins] ) continue;
      for (int ib = 0; ib < nbins; ib++)
	{

	  if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
	    {
	      pt1_reco_bin = ib;
	    }
	  if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
	    {
	      pt2_reco_bin = ib;
	    }
	}
	  
	 
      bool fill_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut);

      if (fill_good)
	{
	  h_d_mbd_vertex->Fill(d_mbd_vertex);
	  h_d_centrality->Fill(d_centrality);
	}
    }

  h_mbd_vertex->Scale(1./h_mbd_vertex->Integral(), "width");
  h_d_mbd_vertex->Scale(1./h_d_mbd_vertex->Integral(), "width");
  TCanvas *c5 = new TCanvas("c5","c5", 500, 500);
  dlutility::ratioPanelCanvas(c5, 0.4);
  c5->cd(1);
  gPad->SetBottomMargin(0.1);
  dlutility::SetLineAtt(h_mbd_vertex, color_sim, 1, 1);
  dlutility::SetLineAtt(h_d_mbd_vertex, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_mbd_vertex, color_sim, 0.5, 8);
  dlutility::SetMarkerAtt(h_d_mbd_vertex, color_data, 0.5, 8);
  dlutility::SetFont(h_d_mbd_vertex, 42, 0.05);
  h_d_mbd_vertex->SetMinimum(0);

  h_d_mbd_vertex->SetTitle("; z_{vtx}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dz_{vtx}} ");
  h_d_mbd_vertex->Draw("p");
  h_mbd_vertex->Draw("same p");

  dlutility::DrawSPHENIX(0.7, 0.8);
  TLegend *leg1 = new TLegend(0.22, 0.65, 0.4, 0.85);
  leg1->SetLineWidth(0);
  leg1->SetTextSize(0.04);
  leg1->SetTextFont(42);
  leg1->AddEntry(h_d_mbd_vertex,"Data","p");
  leg1->AddEntry(h_mbd_vertex,"Sim Reco - HIJING","p");
  leg1->Draw("same");  

  c5->cd(2);
  gPad->SetTopMargin(0.05);

  TH1D *h_compare = (TH1D*) h_d_mbd_vertex->Clone();
  h_compare->SetName("h_mbd_reweight");
  h_compare->SetTitle("; z_{vtx} [cm]; Data/MC");
  h_compare->Divide(h_mbd_vertex);

  dlutility::SetLineAtt(h_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_compare, kBlack, 1, 8);

  dlutility::SetFont(h_compare, 42, 0.07);

  h_compare->SetMinimum(0.0);
  h_compare->SetMaximum(2);

  h_compare->Draw();
  TLine *linl = new TLine(-60, 1, 60 ,1);
  linl->SetLineColor(kBlack);
  linl->SetLineWidth(2);
  linl->SetLineStyle(4);
  linl->Draw("same");
  /* TLine *linep = new TLine(20, 0.2, 40, 0.2); */
  /* dlutility::SetLineAtt(linep, kBlack, 3, 4); */
  /* linep->Draw("same"); */
  /* TLine *linen = new TLine(20, -0.2, 40, -0.2); */
  /* dlutility::SetLineAtt(linen, kBlack, 3, 4); */
  /* linen->Draw("same"); */
  /* TLine *lin0 = new TLine(20, 0, 40, 0); */
  /* dlutility::SetLineAtt(lin0, kBlack, 2, 1); */
  /* lin0->Draw("same"); */

  c5->Print(Form("%s/unfolding_plots/datasim_mbd_compare_r%02d.png", rb.get_code_location().c_str(), cone_size));
  c5->Print(Form("%s/unfolding_plots/datasim_mbd_compare_r%02d.pdf", rb.get_code_location().c_str(), cone_size));

  h_centrality->Scale(1./h_centrality->Integral(), "width");
  h_d_centrality->Scale(1./h_d_centrality->Integral(), "width");
  c5->cd(1);
  gPad->SetBottomMargin(0.1);
  dlutility::SetLineAtt(h_centrality, color_sim, 1, 1);
  dlutility::SetLineAtt(h_d_centrality, color_data, 1, 1);
  dlutility::SetMarkerAtt(h_centrality, color_sim, 0.5, 8);
  dlutility::SetMarkerAtt(h_d_centrality, color_data, 0.5, 8);
  dlutility::SetFont(h_d_centrality, 42, 0.05);
  h_d_centrality->SetMinimum(0);

  h_d_centrality->SetTitle("; Centrality; #frac{1}{N_{pair}}#frac{dN_{pair}}{dCentrality} ");
  h_d_centrality->Draw("p");
  h_centrality->Draw("same p");

  dlutility::DrawSPHENIX(0.7, 0.8);
  TLegend *leg2 = new TLegend(0.22, 0.65, 0.4, 0.85);
  leg2->SetLineWidth(0);
  leg2->SetTextSize(0.04);
  leg2->SetTextFont(42);
  leg2->AddEntry(h_d_centrality,"Data","p");
  leg2->AddEntry(h_centrality,"Sim Reco - HIJING","p");
  leg2->Draw("same");  

  c5->cd(2);

  TH1D *h_c_compare = (TH1D*) h_d_centrality->Clone();
  h_c_compare->SetName("h_cent_reweight");
  h_c_compare->SetTitle("; Centrality; Data/MC");
  h_c_compare->Divide(h_centrality);

  dlutility::SetLineAtt(h_c_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_c_compare, kBlack, 1, 8);

  dlutility::SetFont(h_c_compare, 42, 0.07);

  h_c_compare->SetMinimum(0.0);
  //h_c_compare->SetMaximum(2);

  h_c_compare->Draw();
  TLine *lin2 = new TLine(0, 1, 100 ,1);
  lin2->SetLineColor(kBlack);
  lin2->SetLineWidth(2);
  lin2->SetLineStyle(4);
  lin2->Draw("same");
  /* TLine *linep = new TLine(20, 0.2, 40, 0.2); */
  /* dlutility::SetLineAtt(linep, kBlack, 3, 4); */
  /* linep->Draw("same"); */
  /* TLine *linen = new TLine(20, -0.2, 40, -0.2); */
  /* dlutility::SetLineAtt(linen, kBlack, 3, 4); */
  /* linen->Draw("same"); */
  /* TLine *lin0 = new TLine(20, 0, 40, 0); */
  /* dlutility::SetLineAtt(lin0, kBlack, 2, 1); */
  /* lin0->Draw("same"); */

  c5->Print(Form("%s/unfolding_plots/datasim_centrality_compare_r%02d.png", rb.get_code_location().c_str(), cone_size));
  c5->Print(Form("%s/unfolding_plots/datasim_centrality_compare_r%02d.pdf", rb.get_code_location().c_str(), cone_size));

  
  return 0;
}

