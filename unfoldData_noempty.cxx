#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"

#endif
#include "../macros/dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

int unfoldData_noempty(const std::string configfile = "binning.config", const int niterations = 20, const int cone_size = 4)
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  std::string data_file = "../tntuples/TNTUPLE_DIJET_r0" + std::to_string(cone_size) + "_v6_4_ana462_2024p010_v001_gl10-all.root";
  //std::string data_file = "../tntuples/TNTUPLE_DIJET_v6_2_ana462_2024p010_v001_gl10-00047289-00048291.root";//TNTUPLE_DIJET_v6_1_ana462_2024p010_v001_gl10-00047352-00047733.root";//TNTUPLE_DIJET_v4_2_ana450_2024p009_gl10-00047352-00047733.root";

  float mbd_vertex;
  float pt1_reco;
  float pt2_reco;
  float dphi_reco;
  float nrecojets;
  float trigger;
  TFile *fin = new TFile(data_file.c_str(), "r");
  TNtuple *tn  = (TNtuple*) fin->Get("tn_dijet");;
  if (!tn)
    {
      std::cout << " no data "<< std::endl;
    }
  tn->SetBranchAddress("pt1_reco", &pt1_reco);
  tn->SetBranchAddress("pt2_reco", &pt2_reco);
  tn->SetBranchAddress("dphi_reco", &dphi_reco);
  tn->SetBranchAddress("njets", &nrecojets);
  tn->SetBranchAddress("trigger", &trigger);
  tn->SetBranchAddress("mbd_vertex", &mbd_vertex);

  read_binning rb(configfile.c_str());

  Int_t read_nbins = rb.get_nbins();
  Int_t primer = rb.get_primer();

  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;

  Int_t njet_sys = rb.get_njet_sys();
  Int_t prior_sys = rb.get_prior_sys();
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;

  if (JER_sys != 0)
    {
      std::cout << "Calculating JER extra = " << JER_sys  << std::endl;
    }
  if (JES_sys != 0)
    {
      std::cout << "Calculating JES extra = " << JES_sys  << std::endl;
    }
  if (prior_sys != 0)
    {
      std::cout << "Calculating prior extra = " << prior_sys  << std::endl;
    }

  
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
  float low_trigger[3] = {0};
  
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ib = 0; ib < 3; ib++)
    {
      low_trigger[ib] = rb.get_low_trigger(ib);
      std::cout <<  low_trigger[ib] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  TString responsepath = "response_matrices/response_matrix_r0" + std::to_string(cone_size);
  if (JES_sys > 0)
    {
      responsepath += "_posJES";
    }
  else if (JES_sys < 0)
    {
      responsepath += "_negJES";
    }
  else if (JER_sys > 0)
    {
      responsepath += "_posJER";
    }
  else if (JER_sys < 0)
    {
      responsepath += "_negJER";
    }
  else if (njet_sys > 0)
    {
      responsepath += "_NJET";
    }
  else if (prior_sys)
    {
      responsepath += "_PRIOR";
    }
  if (primer)
    {
      responsepath = "response_matrices/response_matrix_r0" + std::to_string(cone_size) + "_PRIMER" + std::to_string(primer);
    }
  responsepath += ".root";

  
  TFile *fresponse = new TFile(responsepath.Data(),"r");
  
  RooUnfoldResponse *rooResponse = (RooUnfoldResponse*) fresponse->Get("response_noempty");
  if (!rooResponse)
    {
      std::cout << "no repsonse" << std::endl;
      return 1;
    }
  TH1D *h_flat_truth_pt1pt2 = (TH1D*) fresponse->Get("h_truth_flat_pt1pt2"); 
  if (!h_flat_truth_pt1pt2)
    {
      std::cout << "no truth" << std::endl;
      return 1;
    }
  TH1D *h_flat_reco_pt1pt2 = (TH1D*) fresponse->Get("h_reco_flat_pt1pt2"); 
  if (!h_flat_reco_pt1pt2)
    {
      std::cout << "no reco" << std::endl;
      return 1;
    }
  TH1D *h_flat_truth_mapping = (TH1D*) fresponse->Get("h_flat_truth_mapping"); 
  if (!h_flat_truth_mapping)
    {
      std::cout << "no truth" << std::endl;
      return 1;
    }
  TH1D *h_flat_reco_mapping = (TH1D*) fresponse->Get("h_flat_reco_mapping"); 
  if (!h_flat_reco_mapping)
    {
      std::cout << "no reco" << std::endl;
      return 1;
    }
  TH1D *h_flat_truth_skim = (TH1D*) fresponse->Get("h_flat_truth_skim"); 
  if (!h_flat_truth_skim)
    {
      std::cout << "no truth" << std::endl;
      return 1;
    }
  TH1D *h_flat_reco_skim = (TH1D*) fresponse->Get("h_flat_reco_skim"); 
  if (!h_flat_reco_skim)
    {
      std::cout << "no reco" << std::endl;
      return 1;
    }

  TH1D *h_mbd_vertex = new TH1D("h_mbd_vertex", ";z_{vtx}; counts", 120, -60, 60);
  TH1D *h_njets = new TH1D("h_njets", ";N_{Jet}; counts", 51, -0.5, 50.5);
  TH1D *h_reco_xj = new TH1D("h_reco_xj",";A_{J};1/N", nbins, ixj_bins);

  TH2D *h_pt1pt2 = new TH2D("h_pt1pt2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_flat_data_pt1pt2 = new TH1D("h_data_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);

  int nbin_response = nbins*nbins;
  
  int entries = tn->GetEntries();
  for (int i = 0; i < entries; i++)
    {
      tn->GetEntry(i);

      float maxi = std::max(pt1_reco, pt2_reco);
      float mini = std::min(pt1_reco, pt2_reco);

      if (trigger < 17)
	{
	  continue;
	}
      if (nrecojets > 9) continue;
      float pt1_reco_bin = nbins;
      float pt2_reco_bin = nbins;

      float es1 = pt1_reco;
      float es2 = pt2_reco;
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
	  
	 
      bool reco_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_reco > dphicut);

      if (!reco_good) continue;
      
      h_pt1pt2->Fill(es1, es2);
      h_pt1pt2->Fill(es2, es1);
      h_flat_data_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
      h_flat_data_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);
      h_reco_xj->Fill(mini/maxi);
      h_mbd_vertex->Fill(mbd_vertex);
      h_njets->Fill(nrecojets);
    }

  h_flat_data_pt1pt2->Scale(.5);
  TH1D *h_flat_data_skim = (TH1D*) h_flat_reco_skim->Clone();
  h_flat_data_skim->SetName("h_flat_data_skim");
  h_flat_data_skim->Reset();
  histo_opps::skim_down_histo(h_flat_data_skim, h_flat_data_pt1pt2, h_flat_reco_mapping); 
  
  TH1D* h_flat_unfold_skim[niterations];
  TH1D* h_flat_unfold_pt1pt2[niterations];
  int niter = 1;
  for (int iter = 0; iter < niterations; iter++ )
    {
      
      RooUnfoldBayes   unfold (rooResponse, h_flat_data_skim, iter + 1);    // OR
      h_flat_unfold_skim[iter] = (TH1D*) unfold.Hunfold();
      std::cout <<" Nbins skim reco = "<<h_flat_unfold_skim[iter]->GetNbinsX()<<std::endl;
      h_flat_unfold_pt1pt2[iter] = (TH1D*) h_flat_truth_pt1pt2->Clone();
      h_flat_unfold_pt1pt2[iter]->Reset();
      h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d",iter));
      histo_opps::fill_up_histo(h_flat_unfold_skim[iter], h_flat_unfold_pt1pt2[iter], h_flat_truth_mapping);
    }

  
  TCanvas *c = new TCanvas("c","c", 500, 500);

  for (int iter = 0; iter < niterations; iter++)
    {
      dlutility::SetLineAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_flat_unfold_pt1pt2[iter], kBlack, 1, 8);
    }
  
  dlutility::SetLineAtt(h_flat_truth_pt1pt2, kRed, 2, 1);
  dlutility::SetLineAtt(h_flat_data_pt1pt2, kBlue, 2, 1);

  h_flat_truth_pt1pt2->Draw("hist");
  h_flat_data_pt1pt2->Draw("hist same");
  h_flat_unfold_pt1pt2[niter]->Draw("same p");

  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";p_{T1};p_{T2}",nbins, ipt_bins, nbins, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }

  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_truth = new TH1D("h_xj_truth", ";x_{J};",nbins, ixj_bins);
  TH1D *h_xj_unfold[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d", iter), ";x_{J};",nbins, ixj_bins);
    }

  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins);
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins);
    }

  h_pt1pt2_data->SetTitle(";Data p_{T, 1} [GeV]; Data p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_truth->SetTitle(";Truth p_{T, 1} [GeV]; Truth p_{T, 2} [GeV]; Counts * lumi scale ");
  h_pt1pt2_unfold[niter]->SetTitle(";Unfold p_{T, 1} [GeV]; Unfold p_{T, 2} [GeV]; Counts * lumi scale ");
      
  TCanvas *cpt1pt2 = new TCanvas("cpt1pt2","cpt1pt2", 800, 300);
  cpt1pt2->Divide(3, 1);
  cpt1pt2->cd(1);
  gPad->SetLogz();
  gPad->SetRightMargin(0.2);
  h_pt1pt2_data->Draw("colz");
  dlutility::DrawSPHENIXpp(0.22, 0.85);
  
  cpt1pt2->cd(2);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_pt1pt2_truth->Draw("colz");
  cpt1pt2->cd(3);
  gPad->SetRightMargin(0.2);
  gPad->SetLogz();
  h_pt1pt2_unfold[niter]->Draw("colz");


  cpt1pt2->Print(Form("pt1pt2_r%02d.pdf", cone_size));

  histo_opps::project_xj(h_pt1pt2_data, h_xj_data, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
    }


  histo_opps::normalize_histo(h_xj_truth, nbins);
  histo_opps::normalize_histo(h_xj_data, nbins);
  histo_opps::normalize_histo(h_reco_xj, nbins);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::normalize_histo(h_xj_unfold[iter], nbins);
    }

  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
  dlutility::ratioPanelCanvas(cxj);
  cxj->cd(1);
  dlutility::SetLineAtt(h_xj_unfold[niter], kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_xj_unfold[niter], kBlack, 1, 8);

  dlutility::SetLineAtt(h_xj_truth, kRed, 2, 1);
  dlutility::SetMarkerAtt(h_xj_truth, kRed, 1, 8);

  dlutility::SetLineAtt(h_xj_data, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_data, kBlue, 1, 8);

  dlutility::SetLineAtt(h_reco_xj, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_reco_xj, kBlue, 1, 24);


  dlutility::SetFont(h_xj_unfold[niter], 42, 0.04);
  //h_xj_truth->GetXaxis()->SetRangeUser(0.3, 1.001);
  h_xj_truth->SetTitle(";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");
  h_xj_truth->SetMaximum(6);
  h_xj_truth->SetMinimum(0);
  h_xj_truth->Draw("p");
  h_xj_unfold[niter]->Draw("same p");
  h_xj_data->Draw("hist same");
  h_reco_xj->Draw("hist same");
  h_xj_truth->Draw("same hist");
  h_xj_data->Draw("p same");
  h_xj_unfold[niter]->Draw("same hist");
  h_xj_unfold[niter]->Draw("same p");
  dlutility::DrawSPHENIXpp(0.22, 0.84);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
  dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
  dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
  dlutility::drawText("#Delta#phi #geq 3#pi/3", 0.22, 0.59);
  dlutility::drawText("\\mathscr{L} = 25.7 pb^{-1}", 0.22, 0.54);

  TLegend *leg = new TLegend(0.22, 0.35, 0.4, 0.50);
  leg->SetLineWidth(0);
  leg->AddEntry(h_xj_data, "Data");
  leg->AddEntry(h_xj_truth, "Pythia8");
  leg->AddEntry(h_xj_unfold[niter], "Unfolded");
  leg->Draw("same");
  dlutility::drawText(Form("Niter = %d", niter + 1), 0.22, 0.25);
    
  cxj->cd(2);

  TH1D *h_data_compare = (TH1D*) h_xj_unfold[niter]->Clone();
  h_data_compare->Divide(h_xj_truth);
  //h_data_compare->GetXaxis()->SetRangeUser(0.3, 1.001);
  h_data_compare->SetTitle(";x_{J}; Unfold / Truth");
  dlutility::SetFont(h_data_compare, 42, 0.06);
  dlutility::SetLineAtt(h_data_compare, kBlack, 1,1);
  dlutility::SetMarkerAtt(h_data_compare, kBlack, 1,8);

 
  h_data_compare->Draw("p");
  TLine *line = new TLine(0.3, 1, 1, 1);
  line->SetLineStyle(4);
  line->SetLineColor(kRed + 3);
  line->SetLineWidth(2);
  line->Draw("same");

  TCanvas *c_iter = new TCanvas("c_iter", "c_iter");
  TH1D *h_closure[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_closure[iter] = (TH1D*) h_xj_unfold[iter]->Clone();
      h_closure[iter]->SetName(Form("h_closure_%d", iter));
		
      h_closure[iter]->Divide(h_xj_truth);
      h_closure[iter]->SetTitle(";x_{J}; Unfold / Truth");
      dlutility::SetFont(h_closure[iter], 42, 0.06);
    }
  int colors[5] = {kBlue, kBlue - 7, kBlue - 9, kYellow - 7, kYellow +1}; 

  for (int iter = 0; iter < 5; iter++)
    {
      dlutility::SetLineAtt(h_closure[1 + 2*iter], colors[iter], 1,1);
      dlutility::SetMarkerAtt(h_closure[1 + 2*iter], colors[iter], 1,8);
    }
  h_closure[1]->Draw("p");
  h_closure[3]->Draw("p same");
  h_closure[5]->Draw("p same");
  h_closure[7]->Draw("p same");
  h_closure[9]->Draw("p same");
  TLine *line2 = new TLine(0.1, 1, 1, 1);
  line2->SetLineStyle(4);
  line2->SetLineColor(kRed + 3);
  line2->SetLineWidth(2);
  line2->Draw("same");



  dlutility::SetLineAtt(h_xj_data, kBlue, 2, 1);
  dlutility::SetMarkerAtt(h_xj_data, kBlue, 1, 8);

  dlutility::SetLineAtt(h_reco_xj, kBlack, 2, 1);
  dlutility::SetMarkerAtt(h_reco_xj, kBlack, 1, 24);
  
  TCanvas *cproj = new TCanvas("cproj","cproj", 500, 700);
  dlutility::ratioPanelCanvas(cproj);

  cproj->cd(1);
  h_xj_data->SetMaximum(3);
  h_xj_data->SetTitle(";x_{J}; #frac{1}{N_{pairs}}#frac{dN_{pair}}{dx_{J}}");
  dlutility::SetFont(h_xj_data, 42, 0.04);
  
  h_xj_data->Draw();
  h_reco_xj->Draw("same");
  dlutility::DrawSPHENIXpp(0.22, 0.84);
  dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
  dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[measure_leading_bin], ipt_bins[nbins - 1]), 0.22, 0.69);
  dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
  TLegend *legp = new TLegend(0.22, 0.45, 0.4, 0.55);
  legp->SetLineWidth(0);
  legp->SetTextFont(42);
  legp->SetTextSize(0.04);
  legp->AddEntry(h_reco_xj, "Data Filled");
  legp->AddEntry(h_xj_data, "Data Projected");
  legp->Draw("same");

  cproj->cd(2);
  TH1D *h_fillproj_compare = (TH1D*) h_xj_data->Clone();

  dlutility::SetFont(h_fillproj_compare, 42, 0.06);
  h_fillproj_compare->Divide(h_reco_xj);

  h_fillproj_compare->SetMaximum(1.2);
  
  h_fillproj_compare->SetMinimum(0.8);
  h_fillproj_compare->SetTitle(";x_{J};Projected / Filled");
  dlutility::SetLineAtt(h_fillproj_compare, kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_fillproj_compare, kBlack, 1, 8);

  h_fillproj_compare->Draw();
  TLine *line3 = new TLine(0.1, 1, 1, 1);
  line3->SetLineStyle(4);
  line3->SetLineColor(kRed + 3);
  line3->SetLineWidth(2);
  line3->Draw("same");
  cproj->Print(Form("proj_compar_r%02d.pdf", cone_size));
  
  TString unfoldpath = "unfolded_hists/unfolded_hists_r0" + std::to_string(cone_size);
  if (JES_sys > 0)
    {
      unfoldpath += "_posJES";
    }
  if (JES_sys < 0)
    {
      unfoldpath += "_negJES";
    }

  else if (JER_sys > 0)
    {
      unfoldpath += "_posJER";
    }
  else if (JER_sys < 0)
    {
      unfoldpath += "_negJER";
    }
  else if (njet_sys > 0)
    {
      unfoldpath += "_NJET";
    }
  else if (prior_sys)
    {
      unfoldpath += "_PRIOR";
    }

  if (primer)
    {
      unfoldpath = "unfolded_hists/unfolded_hists_r0" + std::to_string(cone_size) + "_PRIMER" + std::to_string(primer);
    }

  unfoldpath += ".root";
  
  TFile *fout = new TFile(unfoldpath.Data(),"recreate");
  h_flat_data_skim->Write();
  h_flat_data_pt1pt2->Write();
  h_flat_reco_pt1pt2->Write();
  h_flat_truth_pt1pt2->Write();
  h_mbd_vertex->Write();
  h_njets->Write();
  
  for (int iter = 0; iter < niterations; ++iter)
    {
      h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d", iter));
      h_flat_unfold_pt1pt2[iter]->Write();
    }
  fout->Close();
  

  return 0;
  
}
