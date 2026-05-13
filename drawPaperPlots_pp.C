
#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"
const bool NUCLEAR = true;
int color_pythia = kRed;//kOrange+7;
int color_herwig = kGreen+2;//kSpring-6;
int cone_sizes[7] = {2, 3, 4, 5, 6, 7, 8};
int color_unfold[7] = {kSpring + 2, kOrange - 2, kBlack, kRed + 2, kViolet + 2, kPink + 3, kAzure - 5};
int color_fill_unfold[7] = {kSpring + 3, kOrange -  1, kAzure - 4, kRed + 3, kViolet + 3};
const int color_unfold_fill = kAzure - 8;
const float marker_unfold = 20;
const float msize_unfold = 0.9;
const float lsize_unfold = 1.1;

void drawPaperPlots_pp(int cone1 = 4, int cone2 = 6)
{
  std::cout << __LINE__ << std::endl;
  gStyle->SetCanvasPreferGL(0);
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();
  std::cout << __LINE__ << std::endl;
  read_binning rb("binning.config");
  std::cout << __LINE__ << std::endl;
  Double_t first_xj = rb.get_first_xj();
  std::cout << __LINE__ << std::endl;
  Int_t read_nbins = rb.get_nbins();
  std::cout << __LINE__ << std::endl;
  Double_t dphicut = rb.get_dphicut();
  std::cout << __LINE__ << std::endl;
  const int nbins = read_nbins;
  const int nbins_pt = read_nbins+1;
  int first_bin = 0;
  float ipt_bins[nbins_pt+1];
  double  dxj_bins[nbins+1];

  float ixj_bins[nbins+1];

  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);
  std::cout << __LINE__ << std::endl;
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
  std::cout << __LINE__ << std::endl;
  const int niterations = 10;
  //const int niter = 2;
  std::cout << __LINE__ << std::endl;
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
    }
  
  TFile *finalout[7];

  TGraphAsymmErrors *g_final_shift_xj_truth_range[7][3];
  TGraphAsymmErrors *g_final_shift_xj_herwig_range[7][3];
  TGraphAsymmErrors *g_final_shift_xj_sys_range[7][3];
  TGraphAsymmErrors *g_final_shift_xj_unfold_range[7][3];
  TH1D *h_final_xj_data_range[7][3];
  TH1D *h_average_xj_truth[7];
  TH1D *h_average_xj_herwig[7];
  
  
  int niter = 1;
  for (int ic = 1; ic < 7; ic++)
    {
      int cone_size = ic+2;
      std::cout << cone_size << std::endl;
      finalout[ic] = new TFile(Form("%s/final_plots/final_plots_pp_r0%d.root",   rb.get_code_location().c_str(), cone_size),"r");

      if (!finalout[ic])
	{
	  std::cout << "no file " << ic+2 << std::endl;
	}
      
      for (int irange = 0; irange < mbins; irange++)
	{
	  g_final_shift_xj_unfold_range[ic][irange] = (TGraphAsymmErrors*) finalout[ic]->Get(Form("g_final_xj_statistics_%d_%d", irange, niter));
	  if (!g_final_shift_xj_unfold_range[ic][irange])
	    {
	      std::cout << "nohist 1 : " << ic+2 << std::endl;
	    }
	  g_final_shift_xj_herwig_range[ic][irange] = (TGraphAsymmErrors*) finalout[ic]->Get(Form("g_final_xj_herwig_%d", irange));
	  if (!g_final_shift_xj_herwig_range[ic][irange])
	    {
	      std::cout << "nohist 1 : " << ic+2 << std::endl;
	    }

	  g_final_shift_xj_truth_range[ic][irange] = (TGraphAsymmErrors*) finalout[ic]->Get(Form("g_final_xj_truth_%d", irange));
	  if (!g_final_shift_xj_truth_range[ic][irange])
	    {
	      std::cout << "nohist 1 : " << ic+2 << std::endl;
	    }
	  g_final_shift_xj_sys_range[ic][irange] = (TGraphAsymmErrors*) finalout[ic]->Get(Form("g_final_xj_systematics_%d_%d", irange, niter));
	  if (!g_final_shift_xj_sys_range[ic][irange])
	    {
	      std::cout << "nohist 1 : " << ic+2 << std::endl;
	    }
	  dlutility::SetLineAtt(g_final_shift_xj_herwig_range[ic][irange], color_herwig, 2, 1);
	  dlutility::SetLineAtt(g_final_shift_xj_truth_range[ic][irange], color_pythia, 2, 1);
	  g_final_shift_xj_sys_range[ic][irange]->SetName(Form("g_final_shift_xj_sys_range_%d_r%02d", irange, cone_size));
	  g_final_shift_xj_sys_range[ic][irange]->SetFillColorAlpha(color_unfold_fill, 0.3); 
	  g_final_shift_xj_unfold_range[ic][irange]->SetName(Form("g_final_shift_xj_unfold_range_%d_r%02d", irange, cone_size));
	  g_final_shift_xj_sys_range[ic][irange]->SetName(Form("g_final_shift_xj_sys_range_%d_r%02d", irange, cone_size));
	  g_final_shift_xj_truth_range[ic][irange]->SetName(Form("g_final_shift_xj_truth_range_%d_r%02d", irange, cone_size));

	}
    }
  TCanvas *c1 = new TCanvas("c1", "Multi-Pad Layout", 1200, 600);
    
  // Define the split between the plots (left) and the text (right)
  double splitX = 0.8; 
    
  // Main pad for the 2x3 grid
  TPad *gridContainer = new TPad("gridContainer", "", 0, 0, splitX, 1);
  gridContainer->SetMargin(0, 0, 0, 0); // Remove container margins
  gridContainer->Draw();
    
  // Pad for the text on the far right
  c1->cd();
  TPad *textPad = new TPad("textPad", "Text Pad", splitX, 0, 1, 1);
  textPad->SetFillColor(0);
  textPad->Draw();
  textPad->cd();
  float size = 0.1;
  float x = 0.01;
  dlutility::DrawSPHENIXpp(x, 0.84, size);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", measure_subleading_cut), x, 0.73, 0, kBlack, size);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", x, 0.65, 0, kBlack, size);
  dlutility::drawText("|#eta| #leq 1.1 - #it{R}", x, 0.57, 0, kBlack, size);
  TLegend *leg = new TLegend(x, 0.2, 0.5, 0.5);
  leg->SetLineWidth(0);
  leg->SetTextFont(43);
  leg->SetTextSize(20);
  leg->AddEntry(g_final_shift_xj_sys_range[2][1], "Data");
  leg->AddEntry(g_final_shift_xj_truth_range[2][1], "PYTHIA-8","l");
  leg->AddEntry(g_final_shift_xj_herwig_range[2][1], "HERWIG 7.3","l");
  leg->Draw("same");
  // Create the 2x3 grid inside the containero
  gridContainer->cd();
  TPad *p[2][3];
    
  double xStep = 1.0 / 3.0;
  double yStep = 1.0 / 2.0;

  double marginleft = 0.1;
  double marginright = 0.01;
  double marginbottom = 0.15;
  double margintop = 0.05;
  double y_div = (1 - margintop - marginbottom)/2.0;
  double x_div = (1 - marginright - marginleft)/3.0;
  double x1[2][3] = {0};
  double y1[2][3] = {0};
  double x2[2][3] = {0};
  double y2[2][3] = {0};
    
  for (int i = 0; i < 2; i++)
    {
      for ( int j  = 0; j < 3; j++)
	{
	  double xx1 = 0;
	  double yy1 = 0;
	  double xx2 = 1;
	  double yy2 = 1;

	  if (j != 2)
	    {
	      xx2 = marginleft + (j+1)*x_div;
	    }
	  if (j > 0)
	    {
	      xx1 = marginleft + j*x_div;
	    }
	  if (i != 1)
	    {
	      yy1 = marginbottom + (1 - i)*y_div;
	    }
	  if (i > 0)
	    {
	      yy2 = marginbottom + i*y_div;
	    }
	  x1[i][j] = xx1;
	  y1[i][j] = yy1;
	  x2[i][j] = xx2;
	  y2[i][j] = yy2;
	}
    }
  for (int i = 0; i < 2; i++) { // Rows (top to bottom)
    for (int j = 0; j < 3; j++) { // Columns (left to right)
            

      TString name = Form("p_%d_%d", i, j);
      p[i][j] = new TPad(name, name, x1[i][j], y1[i][j], x2[i][j], y2[i][j]);
            
      // --- The "Merging" Logic ---
      // Remove right margins for columns 0 and 1
      if (j < 2) p[i][j]->SetRightMargin(0);
      else
	{
	  p[i][j]->SetRightMargin(marginright/(x2[i][j] - x1[i][j]));
	}
      // Remove left margins for columns 1 and 2
      if (j > 0) p[i][j]->SetLeftMargin(0);
      else
	{
	  p[i][j]->SetLeftMargin(marginleft/(x2[i][j] - x1[i][j]));
	}
      // Remove top margins for the bottom row
      if (i == 1) p[i][j]->SetTopMargin(0);
      else
	{
	  p[i][j]->SetTopMargin(margintop/(y2[i][j] - y1[i][j]));
	}
      // Remove bottom margins for the top row
      if (i == 0) p[i][j]->SetBottomMargin(0);
      else
	{
	  p[i][j]->SetBottomMargin(marginbottom/(y2[i][j] - y1[i][j]));
	}

      p[i][j]->Draw();
      p[i][j]->cd();
	    
      TH1D *h = new TH1D(Form("h%d%d",i,j), ";x_{J};#drac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}", 10, first_xj, 1);
      h->SetMaximum(5);
      float tfont = 20;//;// * (x2[0][0] - x1[0][0])/(x2[i][j] - x1[i][j]);
      float infont = 20;;//0.07;// * (y2[0][0] - y1[0][0])/(y2[i][j] - y1[i][j]);

      dlutility::SetFont(h, 43, tfont);
      h->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;
      h->GetYaxis()->SetTitleOffset(1.7);
      h->Draw();
      // Dummy histogram to visualize the merge
      int index = 0;
      if (i == 0)
	{
	  index = cone1 - 2;
	}
      else
	{
	  index = cone2 - 2;
	}
      g_final_shift_xj_truth_range[index][j]->Draw("l same");
      g_final_shift_xj_herwig_range[index][j]->Draw("l same");

      g_final_shift_xj_sys_range[index][j]->Draw("same p E2");
      g_final_shift_xj_unfold_range[index][j]->Draw("same p E1");

      double xtext = 0.1;
      double ytext = 0.85;
      double ytext2 = 0.75;
      if (j == 0){
	xtext = marginleft/(x2[i][j] - x1[i][j]) + xtext*x_div/(x2[i][j] - x1[i][j]);
      }
      else if (j == 2)
	{
	  xtext = xtext*x_div/(x2[i][j] - x1[i][j]);
	}
	    
      if (i == 0)
	{
	  ytext = ytext*y_div/(y2[i][j] - y1[i][j]);
	  ytext2 = ytext2*y_div/(y2[i][j] - y1[i][j]);
	}
      else if (i == 1)
	{
	  ytext = marginbottom /(y2[i][j] - y1[i][j]) + ytext*y_div/(y2[i][j] - y1[i][j]);
	  ytext2 = marginbottom /(y2[i][j] - y1[i][j]) + ytext2*y_div/(y2[i][j] - y1[i][j]);
	}
	    
      if (j == 0)	    dlutility::drawText(Form("#it{R} = %0.1f", (index+2)*0.1), xtext, ytext, 0, kBlack, infont, 43);
      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[j]], ipt_bins[measure_bins[j+1]]), xtext, ytext2, 0, kBlack, infont, 43);


      // Only show Y axis on the leftmost pads
      if (j != 0) {
	h->GetYaxis()->SetLabelSize(0);
	h->GetYaxis()->SetTitleSize(0);
      }
      // Only show X axis on the bottom row
      if (i != 1) {
	h->GetXaxis()->SetLabelSize(0);
	h->GetXaxis()->SetTitleSize(0);
      }
            
      gridContainer->cd();
    }
  }

  c1->Print(Form("xj_paper_r0%d_r0%d.pdf", cone1, cone2));
  c1->Print(Form("xj_paper_r0%d_r0%d.png", cone1, cone2));
  TFile *fin[7];
  TH1D *h_data_dphi_range[7][3];
  TGraphAsymmErrors *h_data_dphi_range_sys[7][3];
  TGraphErrors *h_data_dphi_range_stat[7][3];
  TH1D *h_truth_match_dphi_range[7][3];
  TH1D *h_herwig_dphi_range[7][3];

  for (int i = 1; i < 7;i++)
    {
      
      fin[i] = new TFile(Form("dphihists/dphi_summary_pp_r0%d.root", i+2),"r");
      for (int irange = 0; irange < 3; irange++)
	{
	  std::cout << i << " / " << irange << std::endl;
	  h_data_dphi_range[i][irange] = (TH1D*) fin[i]->Get(Form("h_data_dphi_range_0_%d", irange));
	  h_data_dphi_range_sys[i][irange] = (TGraphAsymmErrors*) fin[i]->Get(Form("h_data_dphi_range_sys_%d", irange));
	  h_truth_match_dphi_range[i][irange] = (TH1D*) fin[i]->Get(Form("h_truth_match_dphi_range_0_%d", irange));
	  h_herwig_dphi_range[i][irange] = (TH1D*) fin[i]->Get(Form("h_herwig_dphi_range_%d", irange));
	  dlutility::SetLineAtt(h_herwig_dphi_range[i][irange], color_herwig, 2, 1);
	  dlutility::SetLineAtt(h_truth_match_dphi_range[i][irange], color_pythia, 2, 1);
	  std::cout << "done" << std::endl;
	  std::cout << __LINE__<< std::endl;

	  h_data_dphi_range[i][irange]->SetName(Form("h_data_dphi_range_%d_r%02d", irange, i + 2));
	  h_data_dphi_range_stat[i][irange] = new TGraphErrors(h_data_dphi_range[i][irange]);
	  h_data_dphi_range_stat[i][irange]->SetName(Form("h_data_dphi_range_stat_%d_r%02d", irange, i + 2));
	  int nn = h_data_dphi_range_stat[i][irange]->GetN();
	  std::cout << "nn " << nn << std::endl;
	  for (int in = 0; in < nn; in++)
	    {
	      std::cout << in << std::endl;	      
	      double yerr = h_data_dphi_range_stat[i][irange]->GetErrorY(in);
	      h_data_dphi_range_stat[i][irange]->SetPointError(in, 0, yerr);
	    }

	  h_data_dphi_range_stat[i][irange]->SetMarkerSize(1);
	  std::cout << __LINE__<< std::endl;
	  h_truth_match_dphi_range[i][irange]->SetName(Form("h_truth_match_dphi_range_%d_r%02d", irange, i + 2));
	  std::cout << __LINE__<< std::endl;
	  h_herwig_dphi_range[i][irange]->SetName(Form("h_herwig_dphi_range_%d_r%02d", irange, i + 2));
	  std::cout << __LINE__<< std::endl;
	  h_data_dphi_range_sys[i][irange]->SetName(Form("h_data_dphi_range_sys_%d_r%02d", irange, i + 2));
	  h_data_dphi_range_sys[i][irange]->SetFillColorAlpha(color_unfold_fill, 0.3);
	}
    }

  TCanvas *c2 = new TCanvas("c2", "Multi-Pad Layout", 1200, 600);
    
  // Define the split between the plots (left) and the text (right)
    
  // Main pad for the 2x3 grid
  TPad *gridContainer2 = new TPad("gridContainer2", "", 0, 0, splitX, 1);
  gridContainer2->SetMargin(0, 0, 0, 0); // Remove container margins
  gridContainer2->Draw();
    
  // Pad for the text on the far right
  c2->cd();
  TPad *textPad2 = new TPad("textPad2", "Text Pad2", splitX, 0, 1, 1);
  textPad2->SetFillColor(0);
  textPad2->Draw();
  textPad2->cd();
  /* float size = 0.1; */
  /* float x = 0.01; */
  dlutility::DrawSPHENIXpp(x, 0.84, size);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", measure_subleading_cut), x, 0.73, 0, kBlack, size);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", x, 0.65, 0, kBlack, size);
  dlutility::drawText("|#eta| #leq 1.1 - #it{R}", x, 0.57, 0, kBlack, size);
  TLegend *leg2 = new TLegend(x, 0.2, 0.5, 0.5);
  leg2->SetLineWidth(0);
  leg2->SetTextFont(43);
  leg2->SetTextSize(20);
  leg2->AddEntry(h_data_dphi_range[2][1], "Data");
  leg2->AddEntry(h_truth_match_dphi_range[2][1], "PYTHIA-8","l");
  leg2->AddEntry(h_herwig_dphi_range[2][1], "HERWIG 7.3","l");
  leg2->Draw("same");
  // Create the 2x3 grid inside the containero
  gridContainer2->cd();
  TPad *p2[2][3];
    
    
  for (int i = 0; i < 2; i++) { // Rows (top to bottom)
    for (int j = 0; j < 3; j++) { // Columns (left to right)
            

      TString name = Form("p_%d_%d", i, j);
      p2[i][j] = new TPad(name, name, x1[i][j], y1[i][j], x2[i][j], y2[i][j]);
            
      // --- The "Merging" Logic ---
      // Remove right margins for columns 0 and 1
      if (j < 2) p2[i][j]->SetRightMargin(0);
      else
	{
	  p2[i][j]->SetRightMargin(marginright/(x2[i][j] - x1[i][j]));
	}
      // Remove left margins for columns 1 and 2
      if (j > 0) p2[i][j]->SetLeftMargin(0);
      else
	{
	  p2[i][j]->SetLeftMargin(marginleft/(x2[i][j] - x1[i][j]));
	}
      // Remove top margins for the bottom row
      if (i == 1) p2[i][j]->SetTopMargin(0);
      else
	{
	  p2[i][j]->SetTopMargin(margintop/(y2[i][j] - y1[i][j]));
	}
      // Remove bottom margins for the top row
      if (i == 0) p2[i][j]->SetBottomMargin(0);
      else
	{
	  p2[i][j]->SetBottomMargin(marginbottom/(y2[i][j] - y1[i][j]));
	}

      p2[i][j]->Draw();
      p2[i][j]->cd();
	    
      TH1D *h = new TH1D(Form("h%d%d",i,j), ";x_{J};#drac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}", 10, TMath::Pi()*3./4., TMath::Pi());
      h->SetMaximum(10);
      float tfont = 20;//;// * (x2[0][0] - x1[0][0])/(x2[i][j] - x1[i][j]);
      float infont = 20;;//0.07;// * (y2[0][0] - y1[0][0])/(y2[i][j] - y1[i][j]);

      dlutility::SetFont(h, 43, tfont);
      h->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}");;
      h->GetYaxis()->SetTitleOffset(1.7);
      h->Draw();
      // Dummy histogram to visualize the merge
      int index = 0;
      if (i == 0)
	{
	  index = cone1 - 2;
	}
      else
	{
	  index = cone2 - 2;
	}
      h_truth_match_dphi_range[index][j]->Draw("l hist same");
      h_herwig_dphi_range[index][j]->Draw("l hist same");

      h_data_dphi_range_sys[index][j]->Draw("same p E2");
      h_data_dphi_range_stat[index][j]->Draw("same p E0");

      double xtext = 0.1;
      double ytext = 0.85;
      double ytext2 = 0.75;
      if (j == 0){
	xtext = marginleft/(x2[i][j] - x1[i][j]) + xtext*x_div/(x2[i][j] - x1[i][j]);
      }
      else if (j == 2)
	{
	  xtext = xtext*x_div/(x2[i][j] - x1[i][j]);
	}
	    
      if (i == 0)
	{
	  ytext = ytext*y_div/(y2[i][j] - y1[i][j]);
	  ytext2 = ytext2*y_div/(y2[i][j] - y1[i][j]);
	}
      else if (i == 1)
	{
	  ytext = marginbottom /(y2[i][j] - y1[i][j]) + ytext*y_div/(y2[i][j] - y1[i][j]);
	  ytext2 = marginbottom /(y2[i][j] - y1[i][j]) + ytext2*y_div/(y2[i][j] - y1[i][j]);
	}
	    
      if (j == 0)	    dlutility::drawText(Form("#it{R} = %0.1f", (index+2)*0.1), xtext, ytext, 0, kBlack, infont, 43);
      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[j]], ipt_bins[measure_bins[j+1]]), xtext, ytext2, 0, kBlack, infont, 43);


      // Only show Y axis on the leftmost pads
      if (j != 0) {
	h->GetYaxis()->SetLabelSize(0);
	h->GetYaxis()->SetTitleSize(0);
      }
      // Only show X axis on the bottom row
      if (i != 1) {
	h->GetXaxis()->SetLabelSize(0);
	h->GetXaxis()->SetTitleSize(0);
      }
            
      gridContainer2->cd();
    }
  }
  c2->Print(Form("dphi_paper_r0%d_r0%d.pdf", cone1, cone2));
  c2->Print(Form("dphi_paper_r0%d_r0%d.png", cone1, cone2));
  TFile *fpythia = new TFile("truth_hists/PythiaTuneVariationsForDan.root","r");

  TH1D *h1_xJ[6];
  int colors[6] = {kRed + 2, kViolet + 2, kBlue + 3, kCyan - 2, kGreen + 2, kOrange+2};
  for (int i = 0; i < 6; i++)
    {
      std::cout << i << std::endl;
      h1_xJ[i] = (TH1D*) fpythia->Get(Form("h1_xJ_%d", i));
      dlutility::SetLineAtt(h1_xJ[i], colors[i], 3, 1);
    }

  TCanvas *c_money = new TCanvas("c","c", 700, 700);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.03);
  TH1D *hb = new TH1D("hb", ";x_{J};#frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}", 10, first_xj, 1);
  dlutility::SetFont(hb, 42, 0.05);   
  hb->SetMaximum(4);
  hb->GetYaxis()->SetTitleOffset(1.5);
  hb->Draw();
  g_final_shift_xj_sys_range[2][1]->Draw("same p E2");
  g_final_shift_xj_unfold_range[2][1]->Draw("same p E1");
  for (int i = 0; i < 6; i++)
    {
      h1_xJ[i]->Draw("same hist l");
    }
    
  float side = 0.24;
  float side2= 0.66;
  float top = 0.9;
  float top2 = 0.5;
  float fonty = 0.04;
  float ss = 0.05;
  float fonty2 = 0.03;
  float ss2 = 0.04;
  dlutility::DrawSPHENIXpp(side, top, fonty);
  dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", 0.4), side2, top2 - 2*ss2, 0, kBlack, fonty2);
  dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[1]], ipt_bins[measure_bins[1+1]]), side2, top2 - 3*ss2, 0, kBlack, fonty2);
  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", measure_subleading_cut), side2, top2 - 4*ss2, 0, kBlack, fonty2);
  dlutility::drawText("#Delta#phi #geq 3#pi/4", side2, top2 - 5*ss2, 0, kBlack, fonty2);
  dlutility::drawText("|#eta| #leq 0.7", side2, top2 - 6*ss2, 0, kBlack, fonty2);
  TLegend *legg = new TLegend(side, 0.58, 0.38, top - 2*ss2);
  legg->SetLineWidth(0);
  legg->SetTextFont(42);
  legg->SetTextSize(0.03);
  legg->AddEntry(g_final_shift_xj_sys_range[2][1], "Data");
  legg->AddEntry(h1_xJ[0],"FSR #alpha_{S} = 0.10","l");
  legg->AddEntry(h1_xJ[1],"FSR #alpha_{S} = 0.15","l");
  legg->AddEntry(h1_xJ[2],"ISR #alpha_{S} = 0.10","l");
  legg->AddEntry(h1_xJ[3],"ISR #alpha_{S} = 0.15","l");
  legg->AddEntry(h1_xJ[4],"Max ISR Evolution Scale = 0.5","l");
  legg->AddEntry(h1_xJ[5],"Max ISR Evolution Scale = 1.5","l");
  legg->Draw("same");
  c_money->Print("xj_paper_r04_pythia.pdf");
  c_money->Print("xj_paper_r04_pythia.png");
  return;
}
