
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

void drawRadialScan_pp()
{
  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb("binning.config");

  Double_t first_xj = rb.get_first_xj();
  
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;
  const int nbins_pt = read_nbins+1;
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

  TFile *finalsys[7];
  TH2D *h_average_xj[7][niterations];
  for (int ic = 1; ic < 7; ic++)
    {
      int cone_size = ic+2;

      finalsys[ic] = new TFile(Form("%s/uncertainties/systematics_pp_r0%d.root",   rb.get_code_location().c_str(), cone_size),"r");

      if (!finalsys[ic])
	{
	  std::cout << "no file " << ic+2 << std::endl;
	}
      for (int iter = 0; iter < niterations; iter++)
	{
	  h_average_xj[ic][iter] = (TH2D*) finalsys[ic]->Get(Form("h_average_xj_%d", iter));
	  h_average_xj[ic][iter]->SetName(Form("h_average_xj_r%02d_%d", cone_size, iter));
	}
    }
  
  TFile *finalout[7];

  TGraphAsymmErrors *g_final_shift_xj_truth_range[7][3];
  TGraphAsymmErrors *g_final_shift_xj_sys_range[7][3];
  TGraphAsymmErrors *g_final_shift_xj_unfold_range[7][3];
  TH1D *h_final_xj_data_range[7][3];
  TH1D *h_average_xj_truth[7];
  TH1D *h_average_xj_herwig[7];
  
  

  for (int ic = 1; ic < 7; ic++)
    {
      int cone_size = ic+2;

      finalout[ic] = new TFile(Form("%s/final_plots/final_plots_pp_r0%d.root",   rb.get_code_location().c_str(), cone_size),"r");

      if (!finalout[ic])
	{
	  std::cout << "no file " << ic+2 << std::endl;
	}
      h_average_xj_truth[ic] = (TH1D*) finalout[ic]->Get("h_average_xj_truth");
      h_average_xj_truth[ic]->SetName(Form("h_average_xj_truth_r%02d", cone_size));
      h_average_xj_herwig[ic] = (TH1D*) finalout[ic]->Get("h_average_xj_herwig");
      h_average_xj_herwig[ic]->SetName(Form("h_average_xj_herwig_r%02d", cone_size));
      
      for (int irange = 0; irange < mbins; irange++)
	{
	  g_final_shift_xj_unfold_range[ic][irange] = (TGraphAsymmErrors*) finalout[ic]->Get(Form("g_final_shift_xj_unfold_range_%d", irange));
	  if (!g_final_shift_xj_unfold_range[ic][irange])
	    {
	      std::cout << "nohist 1 : " << ic+2 << std::endl;
	    }

	  g_final_shift_xj_sys_range[ic][irange] = (TGraphAsymmErrors*) finalout[ic]->Get(Form("g_final_shift_xj_sys_range_%d", irange));
	  if (!g_final_shift_xj_sys_range[ic][irange])
	    {
	      std::cout << "nohist 1 : " << ic+2 << std::endl;
	    }

	  h_final_xj_data_range[ic][irange] = (TH1D*) finalout[ic]->Get(Form("h_final_xj_data_range_%d", irange));
	  h_final_xj_data_range[ic][irange]->SetName(Form("h_final_xj_data_range_%d_r%02d", irange, cone_size));
	  g_final_shift_xj_sys_range[ic][irange]->SetName(Form("g_final_shift_xj_sys_range_%d_r%02d", irange, cone_size));
	  g_final_shift_xj_unfold_range[ic][irange]->SetName(Form("g_final_shift_xj_unfold_range_%d_r%02d", irange, cone_size));
	}
    }
  TCanvas *cxj_money_all = new TCanvas("cxj_all","cxj_all", 1500, 400);
  dlutility::systematic_split_canvas(cxj_money_all, 3, 0);

  TCanvas *cxj_money = new TCanvas("cxj_money","cxj_money", 1000, 1000);
  for (int irange = 0; irange < mbins; irange++)
    {
      cxj_money->cd();
      for (int ic = 1; ic < 7; ic++)
	{
	  double binwidth = ixj_bins[first_bin + 1] - ixj_bins[first_bin];
	  double fix_shift = binwidth / 6.;

	  int point = 0;
	  for (int b=0; b < (nbins - first_bin); b++)
	    {
	      double xw = ixj_bins[first_bin] - ixj_bins[first_bin]; 
	      double x = ixj_bins[first_bin + b] + 0.5*xw;
	      double ex = fix_shift/2.;
	      double sx= x + (ic-1)*fix_shift;


	      g_final_shift_xj_unfold_range[ic][irange]->SetPointX(b, sx);
	      g_final_shift_xj_sys_range[ic][irange]->SetPointEXhigh(b, ex);
	      g_final_shift_xj_sys_range[ic][irange]->SetPointEXlow(b, ex);
	      g_final_shift_xj_sys_range[ic][irange]->SetPointX(b, sx);

	    }

	}

      int marker_range = 20;
      float msize_range = 1.5;
      for (int ic = 1; ic < 7; ic++)
	{
	  dlutility::SetLineAtt(g_final_shift_xj_sys_range[ic][irange], color_unfold[ic], lsize_unfold, 1);
	  dlutility::SetMarkerAtt(g_final_shift_xj_sys_range[ic][irange], color_unfold[ic], msize_range, marker_range);
	  g_final_shift_xj_sys_range[ic][irange]->SetFillColorAlpha(color_unfold[ic], 0.3); 
	  
	  dlutility::SetLineAtt(g_final_shift_xj_unfold_range[ic][irange], color_unfold[ic], lsize_unfold, 1);
	  dlutility::SetMarkerAtt(g_final_shift_xj_unfold_range[ic][irange], color_unfold[ic], msize_range, marker_range);
	  
	}
	

      
      TH1D *hblank1 = new TH1D("hb","", 1, 0.25, 1);
      
      hblank1->GetYaxis()->SetTitleOffset(1.8);
      dlutility::SetFont(hblank1, 42, 0.05, 0.05, 0.05, 0.05);


      hblank1->SetMaximum(6.0);
      hblank1->SetMinimum(0);
      hblank1->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;

      //ht->Draw("E4 same");

      gPad->SetTopMargin(0.05);
      gPad->SetRightMargin(0.05);
      gPad->SetLeftMargin(0.17);
      gPad->SetBottomMargin(0.17);
      
      hblank1->Draw("p E1");

      for (int ic = 1; ic < 7; ic++)
	{
	  g_final_shift_xj_unfold_range[ic][irange]->Draw("same p E1");

	  g_final_shift_xj_sys_range[ic][irange]->Draw("same p E2");
	}

      float top = 0.88;
      float ss = 0.05;
      float xtitle = 0.24;
      dlutility::DrawSPHENIXpp(xtitle, top);

      dlutility::drawText("anti-#it{k}_{t}", xtitle, top - 2*ss);
      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), xtitle, top - 3*ss);
      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), xtitle, top - 4*ss);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", xtitle, top - 5*ss);

      
      TLegend *leg = new TLegend(0.22, top - 10.5*ss, xtitle+0.2, top - 5.5*ss);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      //leg->AddEntry(g_final_shift_xj_unfold_range[0][irange], "R = 0.2");
      leg->AddEntry(g_final_shift_xj_unfold_range[1][irange], "R = 0.3");
      leg->AddEntry(g_final_shift_xj_unfold_range[2][irange], "R = 0.4");
      leg->AddEntry(g_final_shift_xj_unfold_range[3][irange], "R = 0.5");
      leg->AddEntry(g_final_shift_xj_unfold_range[4][irange], "R = 0.6");
      leg->AddEntry(g_final_shift_xj_unfold_range[5][irange], "R = 0.7");
      leg->AddEntry(g_final_shift_xj_unfold_range[6][irange], "R = 0.8");
      
      
      leg->Draw("same");
      cxj_money->Print(Form("h_final_xj_unfolded_radial_range_%d.png", irange));
      cxj_money->Print(Form("h_final_xj_unfolded_radial_range_%d.pdf", irange));

      cxj_money_all->cd(irange + 1);
      float xx = 0.08;
      if (irange == 0)
	{
	  xx = 0.38;
	}
      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), xx, 0.8);
    }

  cxj_money_all->cd(4);
  
  dlutility::DrawSPHENIXppsize(0.02, 0.9, 0.06);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.02,0.74, 0, kBlack, 0.06);
  dlutility::drawText("|#eta| < 1.1 - R", 0.02,0.63, 0, kBlack, 0.06);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.02, 0.52, 0, kBlack, 0.06);

  TLegend *leg4 = new TLegend(0.02, 0.2, 0.9, 0.58);
  leg4->SetTextSize(0.08);
  leg4->SetTextFont(42);
  leg4->SetLineWidth(0);
  //leg->AddEntry(g_final_shift_xj_unfold_range[0][irange], "R = 0.2");
  leg4->AddEntry(g_final_shift_xj_unfold_range[1][0], "R = 0.3");
  leg4->AddEntry(g_final_shift_xj_unfold_range[2][0], "R = 0.4");
  leg4->AddEntry(g_final_shift_xj_unfold_range[3][0], "R = 0.5");
  leg4->AddEntry(g_final_shift_xj_unfold_range[4][0], "R = 0.6");
  leg4->AddEntry(g_final_shift_xj_unfold_range[5][0], "R = 0.7");
  leg4->AddEntry(g_final_shift_xj_unfold_range[6][0], "R = 0.8");
  leg4->Draw("same");

  cxj_money_all->Print("h_final_xj_unfolded_radial_range_all.pdf");

  TCanvas *ca = new TCanvas("ca","ca", 500, 500);
  TGraphAsymmErrors *g_average_xj[3];
  TGraph *g_average_xj_truth[3];
  TGraph *g_average_xj_herwig[3];
  TGraphAsymmErrors *g_scale_average_xj[3];
  TGraph *g_scale_average_xj_truth[3];
  TGraph *g_scale_average_xj_herwig[3];

  for (int i = 0; i < 3; i++)
    {
      std::cout << "range " << i << std::endl;
      g_scale_average_xj[i] = new TGraphAsymmErrors(6);
      g_scale_average_xj_truth[i] = new TGraph(6);
      g_scale_average_xj_herwig[i] = new TGraph(6);

      g_average_xj[i] = new TGraphAsymmErrors(6);
      g_average_xj_truth[i] = new TGraph(6);
      g_average_xj_herwig[i] = new TGraph(6);

      float scale = 1 - (2 - i)*0.25;
      
      for (int j = 0; j < 6; j++)
	{
	  g_average_xj_truth[i]->SetPoint(j, 0.3 + 0.1*j, h_average_xj_truth[j+1]->GetBinContent(i+1));
	  g_average_xj_herwig[i]->SetPoint(j, 0.3 + 0.1*j, h_average_xj_herwig[j+1]->GetBinContent(i+1));
	  g_average_xj[i]->SetPoint(j, 0.3 + 0.1*(j), h_average_xj[j+1][niter]->GetBinContent(1+i, 1));
	  g_average_xj[i]->SetPointError(j, 0.05, 0.05, h_average_xj[j+1][niter]->GetBinContent(1+i, 3), h_average_xj[j+1][niter]->GetBinContent(1+i, 2));	  

	  g_scale_average_xj_truth[i]->SetPoint(j, 0.3 + 0.1*j, h_average_xj_truth[j+1]->GetBinContent(i+1)*scale);
	  g_scale_average_xj_herwig[i]->SetPoint(j, 0.3 + 0.1*j, h_average_xj_herwig[j+1]->GetBinContent(i+1)*scale);
	  g_scale_average_xj[i]->SetPoint(j, 0.3 + 0.1*(j), h_average_xj[j+1][niter]->GetBinContent(1+i, 1)*scale);
	  g_scale_average_xj[i]->SetPointError(j, 0.05, 0.05, h_average_xj[j+1][niter]->GetBinContent(1+i, 3)*scale, h_average_xj[j+1][niter]->GetBinContent(1+i, 2)*scale);	  

	}
    }
  int kcolor[3] = {kAzure-2, kRed, kGreen+3};
  dlutility::SetMarkerAtt(g_average_xj[1], kcolor[1], 1, 20);
  dlutility::SetMarkerAtt(g_average_xj[2], kcolor[2], 1, 21);
  dlutility::SetMarkerAtt(g_average_xj[0], kcolor[0], 1, 22);
  dlutility::SetLineAtt(g_average_xj[1], kcolor[1], 1, 1);
  dlutility::SetLineAtt(g_average_xj[2], kcolor[2], 1, 1);
  dlutility::SetLineAtt(g_average_xj[0], kcolor[0], 1, 1);

  dlutility::SetLineAtt(g_average_xj_truth[1], kcolor[1], 3, 1);
  dlutility::SetLineAtt(g_average_xj_truth[2], kcolor[2], 3, 1);
  dlutility::SetLineAtt(g_average_xj_truth[0], kcolor[0], 3, 1);

  dlutility::SetLineAtt(g_average_xj_herwig[1], kcolor[1], 3, 8);
  dlutility::SetLineAtt(g_average_xj_herwig[2], kcolor[2], 3, 8);
  dlutility::SetLineAtt(g_average_xj_herwig[0], kcolor[0], 3, 8);

  g_average_xj[0]->SetFillColorAlpha(kcolor[0], 0.3);
  g_average_xj[1]->SetFillColorAlpha(kcolor[1], 0.3);
  g_average_xj[2]->SetFillColorAlpha(kcolor[2], 0.3);

  dlutility::SetMarkerAtt(g_scale_average_xj[1], kcolor[1], 1, 20);
  dlutility::SetMarkerAtt(g_scale_average_xj[2], kcolor[2], 1, 21);
  dlutility::SetMarkerAtt(g_scale_average_xj[0], kcolor[0], 1, 22);
  dlutility::SetLineAtt(g_scale_average_xj[1], kcolor[1], 1, 1);
  dlutility::SetLineAtt(g_scale_average_xj[2], kcolor[2], 1, 1);
  dlutility::SetLineAtt(g_scale_average_xj[0], kcolor[0], 1, 1);

  dlutility::SetLineAtt(g_scale_average_xj_truth[1], kcolor[1], 3, 1);
  dlutility::SetLineAtt(g_scale_average_xj_truth[2], kcolor[2], 3, 1);
  dlutility::SetLineAtt(g_scale_average_xj_truth[0], kcolor[0], 3, 1);

  dlutility::SetLineAtt(g_scale_average_xj_herwig[1], kcolor[1], 3, 8);
  dlutility::SetLineAtt(g_scale_average_xj_herwig[2], kcolor[2], 3, 8);
  dlutility::SetLineAtt(g_scale_average_xj_herwig[0], kcolor[0], 3, 8);

  g_scale_average_xj[0]->SetFillColorAlpha(kcolor[0], 0.3);
  g_scale_average_xj[1]->SetFillColorAlpha(kcolor[1], 0.3);
  g_scale_average_xj[2]->SetFillColorAlpha(kcolor[2], 0.3);

  TH1D *h_bb = new TH1D("h_bb", "; Jet Radius (#it{R}); #LT x_{J} #GT;", 6, 0.25, 0.85);
  h_bb->SetMinimum(0.6);
  h_bb->SetMaximum(1.15);
  h_bb->GetYaxis()->SetTitleOffset(1.5);
  dlutility::SetFont(h_bb, 42, 0.05, 0.05, 0.05, 0.05);

  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.17);
  gPad->SetBottomMargin(0.17);

  h_bb->Draw();

  g_average_xj_truth[2]->Draw("same l");
  g_average_xj_truth[1]->Draw("same l");
  g_average_xj_truth[0]->Draw("same l");

  g_average_xj_herwig[2]->Draw("same l");
  g_average_xj_herwig[1]->Draw("same l");
  g_average_xj_herwig[0]->Draw("same l");

  g_average_xj[2]->Draw("same p E2");
  g_average_xj[1]->Draw("same p E2");
  g_average_xj[0]->Draw("same p E2");

  dlutility::DrawSPHENIXpp(0.23, 0.85);

  dlutility::drawText("anti-#it{k}_{t}", 0.23, 0.85 - 2*0.05);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.23, 0.85 - 3*0.05);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.23, 0.85 - 4*0.05);

  
  TLegend *legrad = new TLegend(0.56, 0.75, 0.75, 0.9);
  legrad->SetHeader("Data");
  legrad->SetLineWidth(0);
  legrad->SetTextSize(0.03);
  legrad->SetTextFont(42);
  //legrad->AddEntry(g_final_shift_xj_unfold_range[0][irange], "R = 0.2");
  int iir = 0;
  legrad->AddEntry(g_average_xj[0],Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV", ipt_bins[measure_bins[iir]], ipt_bins[measure_bins[iir+1]]));
  iir = 1;
  legrad->AddEntry(g_average_xj[1],Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV",ipt_bins[measure_bins[iir]], ipt_bins[measure_bins[iir+1]] ));
  iir = 2;
  legrad->AddEntry(g_average_xj[2],Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV", ipt_bins[measure_bins[iir]], ipt_bins[measure_bins[iir+1]]));
  legrad->Draw("same");

  TLegend *legrad2 = new TLegend(0.56, 0.62, 0.75, 0.72);
  legrad2->SetHeader("Simulation");
  legrad2->SetLineWidth(0);
  legrad2->SetTextSize(0.03);
  legrad2->SetTextFont(42);
  //legrad2->AddEntry(g_final_shift_xj_unfold_range[0][irange], "R = 0.2");
  legrad2->AddEntry(g_average_xj_truth[0],"PYTHIA-8");
  legrad2->AddEntry(g_average_xj_herwig[0],"HERWIG 7.3");
  legrad2->Draw("same");

  ca->Print("h_final_xj_unfolded_radial_average.pdf");

  h_bb->SetMinimum(0.2);
  h_bb->SetMaximum(1.5);
  h_bb->Draw();

  g_scale_average_xj_truth[2]->Draw("same l");
  g_scale_average_xj_truth[1]->Draw("same l");
  g_scale_average_xj_truth[0]->Draw("same l");

  g_scale_average_xj_herwig[2]->Draw("same l");
  g_scale_average_xj_herwig[1]->Draw("same l");
  g_scale_average_xj_herwig[0]->Draw("same l");

  g_scale_average_xj[2]->Draw("same p E2");
  g_scale_average_xj[1]->Draw("same p E2");
  g_scale_average_xj[0]->Draw("same p E2");

  dlutility::DrawSPHENIXpp(0.2, 0.87);

  dlutility::drawText("anti-#it{k}_{t}", 0.2, 0.87 - 2*0.05);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.2, 0.87 - 3*0.05);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.2, 0.87 - 4*0.05);

  legrad = new TLegend(0.52, 0.72, 0.71, 0.92);
  legrad->SetHeader("Data");
  legrad->SetLineWidth(0);
  legrad->SetTextSize(0.03);
  legrad->SetTextFont(42);
  //legrad->AddEntry(g_final_shift_xj_unfold_range[0][irange], "R = 0.2");
  iir = 0;
  legrad->AddEntry(g_average_xj[0],Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV #times 1/2", ipt_bins[measure_bins[iir]], ipt_bins[measure_bins[iir+1]]));
  iir = 1;
  legrad->AddEntry(g_average_xj[1],Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV #times 3/4",ipt_bins[measure_bins[iir]], ipt_bins[measure_bins[iir+1]] ));
  iir = 2;
  legrad->AddEntry(g_average_xj[2],Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV #times 1", ipt_bins[measure_bins[iir]], ipt_bins[measure_bins[iir+1]]));
  legrad->Draw("same");

  legrad2 = new TLegend(0.52, 0.58, 0.71, 0.72);
  legrad2->SetHeader("MC");
  legrad2->SetLineWidth(0);
  legrad2->SetTextSize(0.03);
  legrad2->SetTextFont(42);
  //legrad2->AddEntry(g_final_shift_xj_unfold_range[0][irange], "R = 0.2");
  legrad2->AddEntry(g_average_xj_truth[0],"PYTHIA-8");
  legrad2->AddEntry(g_average_xj_herwig[0],"HERWIG 7.3");
  legrad2->Draw("same");

  ca->Print("h_final_xj_unfolded_radial_average_scale.pdf");

  return;
}
