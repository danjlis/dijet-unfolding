#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"
const bool NUCLEAR = true;

int cone_sizes[7] = {2, 3, 4, 5, 6, 7, 8};
int color_unfold[7] = {kSpring + 2, kOrange - 2, kBlack, kRed + 2, kViolet + 2, kPink + 3, kAzure - 5};
int color_fill_unfold[7] = {kSpring + 3, kOrange -  1, kAzure - 4, kRed + 3, kViolet + 3};
const float marker_unfold = 20;
const float msize_unfold = 0.9;
const float lsize_unfold = 1.1;

void drawRadialScan()
{
  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb("binning.config");

  Double_t first_xj = rb.get_first_xj();
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;
  int first_bin = 0;
  float ipt_bins[nbins+1];
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
   
  const int niterations = 10;
  const int niter = 2;
  TProfile *h_xj_rms[7][mbins][niterations];
  TFile *finu[7];
  TFile *fins[7];
  TFile *fin[7];  
  TH1D *h_total_sys_range[7][mbins][niterations];
  TH1D *h_total_sys_neg_range[7][mbins][niterations];

  TH1D *h_flat_unfold_pt1pt2[7][niterations];
  TH2D *h_pt1pt2_unfold[7][niterations];

  TH1D *h_xj_unfold[7][niterations];
  TH1D *h_xj_unfold_range[7][mbins][niterations];

  TH1D *h_xjunc_unfold_range[7][mbins][niterations];
  TH1D *h_xjunc_unfold[7][niterations];

  TH1D *h_final_xj_unfold_range[7][mbins][niterations];
  TH1D *h_final_xj_systematics[7][mbins][niterations];
  TGraphAsymmErrors *g_final_xj_systematics[7][mbins][niterations];

  for (int ic = 0; ic < 7; ic++)
    {
      finu[ic] = new TFile(Form("uncertainties/uncertainties_pp_r%02d_nominal.root", cone_sizes[ic]),"r");
      if (!finu[ic])
	{
	  std::cout << " no unc " << std::endl;
	  return;
	}

      for (int irange = 0; irange < mbins; irange++)
	{
	  for (int iter = 0; iter < niterations; ++iter)
	    {
	      h_xj_rms[ic][irange][iter] = (TProfile*) finu[ic]->Get(Form("hp_xj_range_%d_%d", irange, iter));
	      h_xj_rms[ic][irange][iter]->SetName(Form("hp_xj_range_%d_%d_r%02d", irange, iter, cone_sizes[ic]));
	    }
	}

       fins[ic] = new TFile(Form("uncertainties/systematics_pp_r%02d.root", cone_sizes[ic]),"r");
       if (!fins[ic])
	 {
	   std::cout << " no sys " << std::endl;
	   return;
	 }
       for (int irange = 0; irange < mbins; irange++)
	 {
	   for (int iter = 0; iter < niterations; ++iter)
	     {
	       h_total_sys_range[ic][irange][iter] = (TH1D*) fins[ic]->Get(Form("h_total_sys_range_%d_iter_%d", irange, iter));
	       h_total_sys_neg_range[ic][irange][iter] = (TH1D*) fins[ic]->Get(Form("h_total_sys_neg_range_%d_iter_%d", irange, iter));
	       h_total_sys_range[ic][irange][iter]->SetName(Form("h_total_sys_range_%d_iter_%d_r%02d", irange, iter, cone_sizes[ic]));
	       h_total_sys_neg_range[ic][irange][iter]->SetName(Form("h_total_sys_neg_range_%d_iter_%d_r%02d", irange, iter, cone_sizes[ic]));
	     }
	 }

       fin[ic] = new TFile(Form("unfolding_hists/unfolding_hists_pp_r%02d_nominal.root", cone_sizes[ic]),"r");
       if (!fin[ic])
	 {
	   std::cout << " no file " << std::endl;
	   return;
	 }

       for (int iter = 0; iter < niterations; iter++)
	 {
	   h_flat_unfold_pt1pt2[ic][iter] = (TH1D*) fin[ic]->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
	   h_flat_unfold_pt1pt2[ic][iter]->SetName(Form("h_flat_unfold_pt1pt2_%d_r%02d", iter, cone_sizes[ic]));
	 }


       for (int iter = 0; iter < niterations; iter++)
	 {
	   h_pt1pt2_unfold[ic][iter] = new TH2D(Form("h_pt1pt2_unfold_iter%d_r%02d", iter, cone_sizes[ic]), ";#it{p}_{T,1};#it{p}_{T,2}",nbins, ipt_bins, nbins, ipt_bins);
	 }

       for (int iter = 0; iter < niterations; iter++)
	 {
	   h_xj_unfold[ic][iter] = new TH1D(Form("h_xj_unfold_iter%d_r%02d", iter, cone_sizes[ic]), ";x_{J};",nbins, ixj_bins);
	 }


       for (int irange = 0; irange < mbins; irange++)
	 {
	   for (int iter = 0; iter < niterations; iter++)
	     {
	       h_xj_unfold_range[ic][irange][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d_r%02d", irange, iter, cone_sizes[ic]), ";x_{J};",nbins, ixj_bins);
	     }
	 }
       
       for (int iter = 0; iter < niterations; iter++)
	 {
	   h_xjunc_unfold[ic][iter] = new TH1D(Form("h_xjunc_unfold_iter%d_r%02d", iter, cone_sizes[ic]), ";x_{J};",nbins, ixj_bins);
	 }


       for (int irange = 0; irange < mbins; irange++)
	 {
	   for (int iter = 0; iter < niterations; iter++)
	     {
	       h_xjunc_unfold_range[ic][irange][iter] = new TH1D(Form("h_xjunc_unfold_%d_iter%d_r%02d", irange, iter, cone_sizes[ic]), ";x_{J};",nbins, ixj_bins);
	     }
	 }

       for (int iter = 0; iter < niterations; iter++)
	 {
	   histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[ic][iter], h_pt1pt2_unfold[ic][iter], nbins);
	 }
       for (int iter = 0; iter < niterations; iter++)
	 {
	   histo_opps::project_xj(h_pt1pt2_unfold[ic][iter], h_xj_unfold[ic][iter], nbins, measure_leading_bin, nbins - 2, measure_subleading_bin, nbins - 2);
	 }
  
       for (int iter = 0; iter < niterations; iter++)
	 {
	   histo_opps::normalize_histo(h_xj_unfold[ic][iter], nbins);
	 }
       for (int irange = 0; irange < mbins; irange++)
	 {
	   for (int iter = 0; iter < niterations; iter++)
	     {
	       h_final_xj_unfold_range[ic][irange][iter] = new TH1D(Form("h_final_xj_unfold_%d_iter%d_r%02d", irange, iter, cone_sizes[ic]), ";x_{J};",nbins, ixj_bins);
	       h_final_xj_systematics[ic][irange][iter] = new TH1D(Form("h_final_xj_systematics_%d_%d_r%02d", irange, iter, cone_sizes[ic]), ";x_{J};",nbins, ixj_bins);
	     }
	 }

       for (int irange = 0; irange < mbins; irange++)
	 {
	   for (int iter = 0; iter < niterations; iter++)
	     {
	       histo_opps::project_xj(h_pt1pt2_unfold[ic][iter], h_xj_unfold_range[ic][irange][iter], nbins, measure_bins[irange], measure_bins[irange+1], measure_subleading_bin, nbins - 2);
	}

	   for (int iter = 0; iter < niterations; iter++)
	     {
	       histo_opps::normalize_histo(h_xj_unfold_range[ic][irange][iter], nbins);
	       histo_opps::normalize_histo(h_xj_rms[ic][irange][iter], nbins);
	     }
	   for (int iter = 0; iter < niterations; iter++)
	     {
	       histo_opps::finalize_xj(h_xj_unfold_range[ic][irange][iter], h_final_xj_unfold_range[ic][irange][iter], nbins, first_xj);
	       histo_opps::set_xj_errors(h_final_xj_unfold_range[ic][irange][iter], h_xj_rms[ic][irange][iter], nbins);
	       h_final_xj_systematics[ic][irange][iter] = (TH1D*) h_final_xj_unfold_range[ic][irange][iter]->Clone();
	       h_final_xj_systematics[ic][irange][iter]->SetName(Form("h_final_xj_systematics_%d_%d_%d", irange, iter, cone_sizes[ic]));
	       g_final_xj_systematics[ic][irange][iter] = histo_opps::get_xj_systematics(h_final_xj_systematics[ic][irange][iter], h_total_sys_neg_range[ic][irange][iter], h_total_sys_range[ic][irange][iter], nbins);
	     }
	 }
    }
  
  TCanvas *cxj_money = new TCanvas("cxj_money","cxj_money", 500, 500);

  for (int irange = 0; irange < mbins; irange++)
    {
      for (int ic = 0; ic < 7; ic++)
	{
	  dlutility::SetLineAtt(h_final_xj_unfold_range[ic][irange][niter], color_unfold[ic], lsize_unfold, 1);
	  dlutility::SetMarkerAtt(h_final_xj_unfold_range[ic][irange][niter], color_unfold[ic], msize_unfold, marker_unfold);

	  dlutility::SetLineAtt(g_final_xj_systematics[ic][irange][niter], color_unfold[ic], lsize_unfold, 1);
	  dlutility::SetMarkerAtt(g_final_xj_systematics[ic][irange][niter], color_unfold[ic], msize_unfold, marker_unfold);
	  g_final_xj_systematics[ic][irange][niter]->SetFillColorAlpha(color_fill_unfold[ic], 0.3); 

	  h_final_xj_unfold_range[ic][irange][niter] = (TH1D*) h_final_xj_unfold_range[ic][irange][niter]->Rebin(nbins - first_bin, Form("h_rebin_unf_%d_%d", irange, ic), &dxj_bins[first_bin]);
	}

      gPad->SetTopMargin(0.05);
      gPad->SetRightMargin(0.05);
      gPad->SetLeftMargin(0.17);
      gPad->SetBottomMargin(0.17);
      h_final_xj_unfold_range[4][irange][niter]->GetYaxis()->SetTitleOffset(1.8);
      dlutility::SetFont(h_final_xj_unfold_range[4][irange][niter], 42, 0.06, 0.04, 0.05, 0.05);


      h_final_xj_unfold_range[4][irange][niter]->SetMaximum(6);
      h_final_xj_unfold_range[4][irange][niter]->SetMinimum(0);
      h_final_xj_unfold_range[4][irange][niter]->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

      //ht->Draw("E4 same");
      h_final_xj_unfold_range[4][irange][niter]->Draw("p E1");
      //g_final_xj_systematics[4][irange][niter]->Draw("same p E2");
      for (int ic = 0; ic < 7; ic++)
	{
	  h_final_xj_unfold_range[ic][irange][niter]->Draw("same p E1");
	  //g_final_xj_systematics[ic][irange][niter]->Draw("same p E2");
	}

      float top = 0.88;
      float ss = 0.05;
      dlutility::DrawSPHENIXpp(0.22, top);

      dlutility::drawText("anti-#it{k}_{t}", 0.22, top - 2*ss);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq #it{p}_{T,1} < %2.1f GeV} ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.22, top - 3*ss);
      dlutility::drawText(Form("#it{p}_{T,2} #kern[-0.07]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, top - 4*ss);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, top - 5*ss);

      
      TLegend *leg = new TLegend(0.65, top - 5.5*ss, 0.9, 0.92);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_final_xj_unfold_range[0][irange][niter], "R = 0.2");
      leg->AddEntry(h_final_xj_unfold_range[1][irange][niter], "R = 0.3");
      leg->AddEntry(h_final_xj_unfold_range[2][irange][niter], "R = 0.4");
      leg->AddEntry(h_final_xj_unfold_range[3][irange][niter], "R = 0.5");
      leg->AddEntry(h_final_xj_unfold_range[4][irange][niter], "R = 0.6");
      leg->AddEntry(h_final_xj_unfold_range[5][irange][niter], "R = 0.7");
      leg->AddEntry(h_final_xj_unfold_range[6][irange][niter], "R = 0.8");
      
      leg->Draw("same");

      cxj_money->Print(Form("h_final_xj_unfolded_radial_range_%d.png", irange));
      cxj_money->Print(Form("h_final_xj_unfolded_radial_range_%d.pdf", irange));
    }

  return;
}
