
#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"
const bool NUCLEAR = true;

const int color_unfold_fill = kAzure - 4;
const int color_unfold = kAzure - 5;
const float marker_unfold = 20;
const float msize_unfold = 0.9;
const float lsize_unfold = 1.1;

const int color_pythia = kRed;
const float msize_pythia = 0.9;
const float marker_pythia = 20;
const float lsize_pythia = 1.1;

const int color_herwig = kViolet;
const float msize_herwig = 0.9;
const float marker_herwig = 20;
const float lsize_herwig = 1.1;

const int color_reco = kRed;
const float marker_reco = 24;
const float msize_reco = 0.9;
const float lsize_reco = 1.1;
const int color_data = kAzure - 6;
const float marker_data = 24;
const float msize_data = 0.9;
const float lsize_data = 1.1;
void drawHalfClosure_pp(const int cone_size = 4, const int primer = 0)
{

  std::string sys_string = "HALF_PuP";
  if (primer)
    {
      sys_string = "PRIMER" + std::to_string(primer) + "_HALF_PuP";
    }

  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  const int niterations = 10;

  read_binning rb("binning.config");

  Double_t first_xj = rb.get_first_xj();
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();


  const int nbins = read_nbins;
  const int nbins_pt = read_nbins + 1;
  int first_bin = 0;
  float ipt_bins[nbins_pt+1];
  double  dxj_bins[nbins+1];

  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  for (int i = 0 ; i < nbins + 1; i++)
    {
      std::cout << ipt_bins[i] << " -- " << ixj_bins[i] << std::endl;
      dxj_bins[i] = ixj_bins[i];
      if (dxj_bins[i] > 0.3 && first_bin == 0) first_bin = i;

    }
  ipt_bins[nbins_pt] = 100;
  
  Int_t max_reco_bin = rb.get_maximum_reco_bin();
  
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
  int measure_bins[10] = {0};
  int subleading_measure_bins[10] = {0};
  for (int ib = 0; ib < 4; ib++)
    {
      sample_boundary[ib] = rb.get_sample_boundary(ib);
      std::cout <<  sample_boundary[ib] << std::endl;
    }
  for (int ir = 0; ir < mbins+1; ir++)
    {
      measure_bins[ir] = rb.get_measure_region(ir);
      subleading_measure_bins[ir] = rb.get_subleading_measure_region(ir);
      std::cout << ipt_bins[measure_bins[ir]] << " -- " <<  ipt_bins[subleading_measure_bins[ir]] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;
  

  TFile *fin = new TFile(Form("%s/response_matrices/response_matrix_pp_r%02d_%s.root ", rb.get_code_location().c_str(), cone_size, sys_string.c_str()),"r");
  if (!fin)
    {
      std::cout << " no file " << std::endl;
      return;
    }
  std::cout << __LINE__ << std::endl;
  TH1D *h_flat_reco_pt1pt2 = (TH1D*) fin->Get("h_reco_flat_pt1pt2");
  TH1D *h_flat_truth_pt1pt2 = (TH1D*) fin->Get("h_truth_flat_pt1pt2");
  TH1D *h_flat_unfold_pt1pt2[niterations];
  for (int iter = 0; iter < niterations; iter++)
    {
      h_flat_unfold_pt1pt2[iter] = (TH1D*) fin->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
    }
  std::cout << __LINE__ << std::endl;
  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";#it{p}_{T,1};#it{p}_{T,2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_truth = new TH2D("h_pt1pt2_truth", ";#it{p}_{T,1};#it{p}_{T,2}", nbins_pt, ipt_bins, nbins_pt, ipt_bins);
  TH2D *h_pt1pt2_unfold[niterations];

  for (int iter = 0; iter < niterations; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D("h_pt1pt2_unfold", ";#it{p}_{T,1};#it{p}_{T,2}",nbins_pt, ipt_bins, nbins_pt, ipt_bins);
      h_pt1pt2_unfold[iter]->SetName(Form("h_pt1pt2_unfold_iter%d", iter));
    }

  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, nbins_pt);
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, nbins_pt);
  for (int iter = 0; iter < niterations; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], nbins_pt);
    }

  std::cout << __LINE__ << std::endl;
  TH1D *h_xj_reco_range[mbins];
  TH1D *h_xj_truth_range[mbins];
  TH1D *h_xj_unfold_range[mbins][niterations];
  for (int irange = 0; irange < mbins; irange++)
    {
      h_xj_reco_range[irange] = new TH1D(Form("h_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_xj_truth_range[irange] = new TH1D(Form("h_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_xj_unfold_range[irange][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  std::cout << __LINE__ << std::endl;
  TH1D *h_final_xj_reco_range[mbins];
  TH1D *h_final_xj_truth_range[mbins];
  TH1D *h_final_xj_unfold_range[mbins][niterations];
  std::cout << __LINE__ << std::endl;
  for (int irange = 0; irange < mbins; irange++)
    {
      h_final_xj_reco_range[irange] = new TH1D(Form("h_final_xj_reco_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      h_final_xj_truth_range[irange] = new TH1D(Form("h_final_xj_truth_range_%d", irange), ";x_{J};", nbins, ixj_bins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_final_xj_unfold_%d_iter%d", irange, iter), ";x_{J};",nbins, ixj_bins);
	}
    }
  std::cout << __LINE__ << std::endl;
  for (int irange = 0; irange < mbins; irange++)
    {
      histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco_range[irange], nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth_range[irange], nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold_range[irange][iter], nbins_pt, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}




      histo_opps::finalize_xj(h_xj_truth_range[irange], h_final_xj_truth_range[irange], nbins, first_xj);
      histo_opps::finalize_xj(h_xj_reco_range[irange], h_final_xj_reco_range[irange], nbins, first_xj);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::finalize_xj(h_xj_unfold_range[irange][iter], h_final_xj_unfold_range[irange][iter], nbins, first_xj);
	}

      histo_opps::normalize_histo(h_final_xj_truth_range[irange], nbins);
      histo_opps::normalize_histo(h_final_xj_reco_range[irange], nbins);
      for (int iter = 0; iter < niterations; iter++)
	{
	  histo_opps::normalize_histo(h_final_xj_unfold_range[irange][iter], nbins);
	}

    }
  std::cout << __LINE__ << std::endl;
  
  TCanvas *cxj = new TCanvas("cxj","cxj", 500, 700);
  dlutility::ratioPanelCanvas(cxj);

  TH1D *h_closure_test[mbins][niterations];
  
  for (int niter = 0; niter < niterations; niter++)
    {
      for (int irange = 0; irange < mbins; irange++)
      {
	std::cout << __LINE__ << std::endl;
	cxj->cd(1);
	dlutility::SetLineAtt(h_final_xj_unfold_range[irange][niter], color_unfold, lsize_unfold, 1);
	dlutility::SetMarkerAtt(h_final_xj_unfold_range[irange][niter], color_unfold, msize_unfold, marker_unfold);

	dlutility::SetLineAtt(h_final_xj_truth_range[irange], color_pythia, lsize_pythia, 1);
	dlutility::SetMarkerAtt(h_final_xj_truth_range[irange], color_pythia, msize_pythia, marker_pythia);

	dlutility::SetFont(h_final_xj_truth_range[irange], 42, 0.05);

	h_final_xj_truth_range[irange]->SetMaximum(5);
	h_final_xj_truth_range[irange]->SetMinimum(0);
	h_final_xj_truth_range[irange]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

	TH1D *ht = (TH1D*) h_final_xj_truth_range[irange]->Clone();
	TH1D *hu = (TH1D*) h_final_xj_unfold_range[irange][niter]->Clone();

	ht->Draw("p E1");
	hu->Draw("same p E1");


	dlutility::DrawSPHENIX(0.22, 0.84);
	dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", cone_size*0.1), 0.22, 0.74);
	dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[1]], ipt_bins[measure_bins[2]]), 0.22, 0.69);
	dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
	dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);

	TLegend *leg = new TLegend(0.2, 0.4, 0.4, 0.56);
	leg->SetLineWidth(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);
	leg->AddEntry(h_final_xj_truth_range[irange], "Training Truth","p");
	leg->AddEntry(h_final_xj_unfold_range[irange][niter], "Testing Unfold","p");
	leg->Draw("same");

	cxj->cd(2);

	TH1D *h_data_compare = (TH1D*) h_final_xj_truth_range[irange]->Clone();

	h_data_compare->Reset();
	for (int i = 0; i < h_data_compare->GetNbinsX(); i++)
	  {
	    if (h_final_xj_truth_range[irange]->GetBinContent(i+1) == 0) continue;


	    h_data_compare->SetBinContent(i+1, (h_final_xj_unfold_range[irange][niter]->GetBinContent(i+1) - h_final_xj_truth_range[irange]->GetBinContent(i+1))/(h_final_xj_truth_range[irange]->GetBinContent(i+1)));
	    float err_unfold = h_final_xj_unfold_range[irange][niter]->GetBinError(i+1);
	    float err_truth = h_final_xj_truth_range[irange]->GetBinError(i+1);

	    float con_unfold = h_final_xj_unfold_range[irange][niter]->GetBinContent(i+1);
	    float con_truth = h_final_xj_truth_range[irange]->GetBinContent(i+1);

	    float err_total = fabs(con_unfold/con_truth) * sqrt(TMath::Power(err_unfold/con_unfold, 2) + TMath::Power(err_truth/con_truth, 2));

	    h_data_compare->SetBinError(i+1, err_total);
	  }
	h_data_compare->SetTitle(";x_{J}; Truth - Unfold / Truth");
	dlutility::SetFont(h_data_compare, 42, 0.1, 0.07, 0.07, 0.07);
	dlutility::SetLineAtt(h_data_compare, kBlack, 1,1);
	dlutility::SetMarkerAtt(h_data_compare, kBlack, 1,8);
	h_data_compare->SetMaximum(0.2);
	h_data_compare->SetMinimum(-0.2);

	TH1D *hd = (TH1D*) h_data_compare->Clone();

	hd->Draw("p");

	h_closure_test[irange][niter] = (TH1D*) hd->Clone();
	h_closure_test[irange][niter]->SetName(Form("h_closure_test_%d_%d", irange, niter));
	  
	TLine *line = new TLine(hd->GetBinLowEdge(1), 1, 1, 1);
	line->SetLineStyle(4);
	line->SetLineColor(kRed + 3);
	line->SetLineWidth(2);
	line->Draw("same");
	cxj->Print(Form("%s/unfolding_plots/h_xj_half_closure_pp_r%02d_range_%d_iter_%d_%s.png", rb.get_code_location().c_str(), cone_size, irange, niter, sys_string.c_str()));
	cxj->Print(Form("%s/unfolding_plots/h_xj_half_closure_pp_r%02d_range_%d_iter_%d_%s.pdf", rb.get_code_location().c_str(), cone_size, irange, niter, sys_string.c_str()));
      }
    }

  for (int irange = 0; irange < mbins; irange++)
    {
      TCanvas *cxjclos = new TCanvas("cxjclos","cxjclos", 700, 500);
      dlutility::createCutCanvas(cxjclos);
      cxjclos->cd(1);

      dlutility::SetLineAtt(h_closure_test[irange][0], kBlue, 1, 1);
      dlutility::SetLineAtt(h_closure_test[irange][1], kCyan, 1, 1);
      dlutility::SetLineAtt(h_closure_test[irange][3], kGreen, 1, 1);
      dlutility::SetLineAtt(h_closure_test[irange][5], kViolet, 1, 1);
      dlutility::SetLineAtt(h_closure_test[irange][7], kRed, 1, 1);

      dlutility::SetMarkerAtt(h_closure_test[irange][0], kBlue, 1, 8);
      dlutility::SetMarkerAtt(h_closure_test[irange][1], kCyan, 1, 8);
      dlutility::SetMarkerAtt(h_closure_test[irange][3], kGreen, 1, 8);
      dlutility::SetMarkerAtt(h_closure_test[irange][5], kViolet, 1, 8);
      dlutility::SetMarkerAtt(h_closure_test[irange][7], kRed, 1, 8);
      dlutility::SetFont(h_closure_test[irange][0], 42, 0.04);
      h_closure_test[irange][0]->SetMinimum(-0.2);
      h_closure_test[irange][0]->SetMaximum(0.2);
      h_closure_test[irange][0]->Draw("p E1");
      //h_closure_test[irange][0]->SetTitleOffset(0.9);
      h_closure_test[irange][0]->SetTitle(";x_{J}; Half Truth - Half Unfold / Half Truth");
      h_closure_test[irange][1]->Draw("p E1 same");
      h_closure_test[irange][3]->Draw("p E1 same");
      h_closure_test[irange][5]->Draw("p E1 same");
      h_closure_test[irange][7]->Draw("p E1 same");    

      cxjclos->cd(2);

      dlutility::DrawSPHENIXcut(0.02, 0.9, 0, 0.08);
      dlutility::drawText(Form("anti-#it{k}_{t} #kern[-0.3]{#it{R}} = %0.1f", 0.1*cone_size), 0.02, 0.78, 0, kBlack, 0.08);
      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.02, 0.73, 0, kBlack, 0.08);
      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.02,0.68, 0, kBlack, 0.08);
      dlutility::drawText("#Delta#phi #geq 7#pi/8", 0.02,0.63, 0, kBlack, 0.08);

      dlutility::drawText("Closure Test", 0.02,0.53, 0, kBlack, 0.08);

      TLegend *leg4 = new TLegend(0.02, 0.25, 0.9, 0.48);
      leg4->SetTextSize(0.08);
      leg4->SetTextFont(42);
      leg4->SetLineWidth(0);
      leg4->AddEntry(h_closure_test[irange][0],"1 iteration","l");
      leg4->AddEntry(h_closure_test[irange][1],"2 iteration","l");
      leg4->AddEntry(h_closure_test[irange][3],"4 iteration","l");
      leg4->AddEntry(h_closure_test[irange][5],"6 iteration","l");
      leg4->AddEntry(h_closure_test[irange][7],"8 iteration","l");
      leg4->Draw();
      cxjclos->Print(Form("%s/unfolding_plots/h_xj_half_closure_pp_r%02d_range_%d_iter_all_%s.pdf", rb.get_code_location().c_str(), cone_size, irange, sys_string.c_str()));
      cxjclos->Print(Form("%s/unfolding_plots/h_xj_half_closure_pp_r%02d_range_%d_iter_all_%s.png", rb.get_code_location().c_str(), cone_size, irange, sys_string.c_str()));

      TCanvas *cclos = new TCanvas("cxjclos","cxjclos", 1000, 800);
      cclos->Divide(3,2);

      int nbins1 = h_closure_test[irange][0]->GetNbinsX();
  
      for (int i = 0; i < 6; i++)
	{
	  cclos->cd(i+1);
     
	  TH1D *h_niter_closure = new TH1D("h_niter_closure",";N_{iterations}; Rel. Error", niterations, 0, niterations);
	  TH1D *h_niter_stat = new TH1D("h_niter_stat",";N_{iterations}; Rel. Error", niterations, 0, niterations);
	  TH1D *h_niter_quad = new TH1D("h_niter_quad",";N_{iterations}; Rel. Error", niterations, 0, niterations);
		  
	  for (int iter = 0; iter < niterations; iter++)
	    {
	      std::cout << fabs(h_closure_test[irange][iter]->GetBinCenter(nbins1 - i)) << std::endl;
	      float closure = fabs(h_closure_test[irange][iter]->GetBinContent(nbins1 - i));
	      h_niter_closure->SetBinContent(iter + 1, fabs(h_closure_test[irange][iter]->GetBinContent(nbins1 - i)));
	      float stat=fabs(h_final_xj_unfold_range[irange][iter]->GetBinError(nbins1 - i)/h_final_xj_truth_range[irange]->GetBinContent(nbins1 - i));
	      h_niter_stat->SetBinContent(iter + 1, stat);
	      h_niter_quad->SetBinContent(iter + 1, sqrt(stat*stat + closure*closure));
	    }

	  h_niter_closure->SetLineColor(kBlack);
	  //h_niter_stat->SetLineColor(kBlue);
	  //h_niter_closure->SetLineColor(kBlack);
	  h_niter_closure->SetMaximum(0.15);
	  h_niter_closure->Draw("hist");

	  dlutility::drawText(Form("x_{J} = %0.2f", h_closure_test[irange][0]->GetBinCenter(nbins1 - i)), 0.22, 0.8);


	}

      cclos->Print(Form("%s/unfolding_plots/h_xj_half_closure_pp_r%02d_range_%d_iter_xjbin_%s.pdf", rb.get_code_location().c_str(), cone_size, irange, sys_string.c_str()));
    }      

  TFile *f_out = new TFile(Form("%s/uncertainties/half_closure_pp_r%02d_%s.root", rb.get_code_location().c_str(), cone_size, sys_string.c_str()),"recreate");
  for (int i = 0; i < 10; i++)
    {
      for (int irange = 0; irange < mbins; irange++)
	{
	  h_closure_test[irange][i]->Write();
	}
    }
  f_out->Close();  
  return;
}
