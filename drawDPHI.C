#include "../macros/dlUtility.h"
#include "read_binning.h"

const int color_herwig = kViolet;
const float msize_herwig = 0.9;
const float marker_herwig = 20;
const float lsize_herwig = 1.1;

const int color_unfold_fill = kAzure - 4;
const int color_unfold = kAzure - 5;
const float marker_unfold = 20;
const float msize_unfold = 0.9;
const float lsize_unfold = 1.1;
const int color_pythia = kRed;
const float msize_pythia = 0.9;
const float marker_pythia = 20;
const float lsize_pythia = 1.1;
const int color_reco = kRed;
const float marker_reco = 24;
const float msize_reco = 0.9;
const float lsize_reco = 1.1;
const int color_data = kAzure - 6;
const float marker_data = 24;
const float msize_data = 0.9;
const float lsize_data = 1.1;


const int colorsys[5] = {kBlack, kMagenta + 2, kGreen - 2, kCyan + 1, kPink + 2};

std::pair<float, float> get_ordered_pair(float e1, float e2);

void drawDPHI(const std::string configfile = "binning.config")
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


  TFile *fher = new TFile("herwig_hist.root","r");
  TH1D *h_herwig_leading_dphi[3];
  TH1D *h_herwig_subleading_dphi[3];
  for (int i = 0; i < 3; i++)
    {
      h_herwig_leading_dphi[i] = (TH1D*) fher->Get(Form("h_truth_leading_dphi_%d", i));
      h_herwig_leading_dphi[i]->Rebin(2);
      h_herwig_leading_dphi[i]->Scale(1./h_herwig_leading_dphi[i]->Integral(0, -1, "width"));

      h_herwig_subleading_dphi[i] = (TH1D*) fher->Get(Form("h_truth_subleading_dphi_%d", i));
      h_herwig_subleading_dphi[i]->Rebin(2);
      h_herwig_subleading_dphi[i]->Scale(1./h_herwig_subleading_dphi[i]->Integral(0, -1, "width"));
      
    }

  float scale_pt1_reco[5][3] = {{0}};
  float scale_pt2_reco[5][3] = {{0}};
  float scale_pt1_truth[5][3] = {{0}};
  float scale_pt2_truth[5][3] = {{0}};
  TString filenames[5];
  filenames[0] = "dphi_hists.root";
  filenames[1] = "dphi_hists_negJER.root";
  filenames[2] = "dphi_hists_posJER.root";
  filenames[3] = "dphi_hists_negJES.root";
  filenames[4] = "dphi_hists_posJES.root";

  TFile *fin[5];
  
  TH3D *h_data_pt1pt2dphi[5];//new TH3D("h_truth_pt1pt2",";p_{T,1, smear};p_{T2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  
  TH3D *h_truth_pt1pt2dphi[5];//new TH3D("h_reco_match_pt1pt2dphi",";p_{T,1, smear};p_{T2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_truth_match_pt1pt2dphi[5];//new TH3D("h_reco_match_pt1pt2dphi",";p_{T,1, smear};p_{T2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphi[5];//new TH3D("h_reco_match_pt1pt2dphi",";p_{T,1, smear};p_{T2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphitruth[5];//new TH3D("h_reco_match_pt1pt2dphitruth",";p_{T,1, smear};p_{T2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);

  TH1D *h_all_reco_dphi_range[7][3];
  TH1D *h_all_truth_dphi_range[7][3];
  TH1D *h_corr_reco_dphi_range[7][3];
  TH1D *h_corr_truth_dphi_range[7][3];

  TH1D *h_all_reco_dphi_rangemin[7][3];
  TH1D *h_all_truth_dphi_rangemin[7][3];
  TH1D *h_corr_reco_dphi_rangemin[7][3];
  TH1D *h_corr_truth_dphi_rangemin[7][3];
  
  TH1D *h_data_dphi_range[7][3];
  TH1D *h_truth_dphi_range[5][3];
  TH1D *h_truth_match_dphi_range[5][3];
  TH1D *h_reco_dphi_range[5][3];

  TH1D *h_data_dphi_rangemin[7][3];
  TH1D *h_truth_dphi_rangemin[5][3];
  TH1D *h_truth_match_dphi_rangemin[5][3];
  TH1D *h_reco_dphi_rangemin[5][3];

  TH1D *h_correction_factors_range[5][3];
  TH1D *h_correction_factors_rangemin[5][3];

  TH1D *h_efficiency_factors_range[5][3];
  TH1D *h_efficiency_factors_rangemin[5][3];

  // Get histograms
  for (int isys = 0; isys < 5; isys++)
    {
      std::cout << filenames[isys] << std::endl;
      fin[isys] = new TFile(filenames[isys].Data(), "r");
      h_data_pt1pt2dphi[isys] = (TH3D*) fin[isys]->Get("h_data_pt1pt2dphi");
      h_reco_match_pt1pt2dphi[isys] = (TH3D*) fin[isys]->Get("h_reco_pt1pt2dphi");;
      h_truth_match_pt1pt2dphi[isys] = (TH3D*) fin[isys]->Get("h_truth_pt1pt2dphi");;
      h_truth_pt1pt2dphi[isys] = (TH3D*) fin[isys]->Get("h_truth_match_pt1pt2dphi");;
      h_reco_match_pt1pt2dphitruth[isys] = (TH3D*) fin[isys]->Get("h_reco_match_pt1pt2dphitruth");;

      
      h_data_pt1pt2dphi[isys]->SetName(Form("h_data_pt1pt2dphi_%d", isys));
      
      h_reco_match_pt1pt2dphi[isys]->SetName(Form("h_reco_pt1pt2dphi_%d", isys));
      
      h_truth_match_pt1pt2dphi[isys]->SetName(Form("h_truth_match_pt1pt2dphi_%d", isys));
      
      h_truth_pt1pt2dphi[isys]->SetName(Form("h_truth_pt1pt2dphi_%d", isys));
      
      h_reco_match_pt1pt2dphitruth[isys]->SetName(Form("h_reco_pt1pt2dphitruth_%d", isys));
      for (int i = 0; i < 3; i++)
	{
	  h_all_reco_dphi_range[isys][i] = (TH1D*) fin[isys]->Get(Form("h_all_counts_leading_reco_%d", i));
	  h_all_truth_dphi_range[isys][i] = (TH1D*) fin[isys]->Get(Form("h_all_counts_leading_truth_%d", i));
	  h_corr_reco_dphi_range[isys][i] = (TH1D*) fin[isys]->Get(Form("h_correlated_counts_leading_reco_%d", i));
	  h_corr_truth_dphi_range[isys][i] = (TH1D*) fin[isys]->Get(Form("h_correlated_counts_leading_truth_%d", i));

	  h_all_reco_dphi_rangemin[isys][i] = (TH1D*) fin[isys]->Get(Form("h_all_counts_subleading_reco_%d", i));
	  h_all_truth_dphi_rangemin[isys][i] = (TH1D*) fin[isys]->Get(Form("h_all_counts_subleading_truth_%d", i));
	  h_corr_reco_dphi_rangemin[isys][i] = (TH1D*) fin[isys]->Get(Form("h_correlated_counts_subleading_reco_%d", i));
	  h_corr_truth_dphi_rangemin[isys][i] = (TH1D*) fin[isys]->Get(Form("h_correlated_counts_subleading_truth_%d", i));

	  h_all_reco_dphi_range[isys][i]->SetName(Form("h_all_counts_leading_reco_%d_%d", i, isys));
	  h_all_truth_dphi_range[isys][i]->SetName(Form("h_all_counts_leading_truth_%d_%d", i, isys));
	  h_corr_reco_dphi_range[isys][i]->SetName(Form("h_correlated_counts_leading_reco_%d_%d", i, isys));
	  h_corr_truth_dphi_range[isys][i]->SetName(Form("h_correlated_counts_leading_truth_%d_%d", i, isys));

	  h_all_reco_dphi_rangemin[isys][i]->SetName(Form("h_all_counts_subleading_reco_%d_%d", i, isys));
	  h_all_truth_dphi_rangemin[isys][i]->SetName(Form("h_all_counts_subleading_truth_%d_%d", i, isys));
	  h_corr_reco_dphi_rangemin[isys][i]->SetName(Form("h_correlated_counts_subleading_reco_%d_%d", i, isys));
	  h_corr_truth_dphi_rangemin[isys][i]->SetName(Form("h_correlated_counts_subleading_truth_%d_%d", i, isys));

	  h_all_reco_dphi_range[isys][i]->SetDirectory(0);
	  h_all_truth_dphi_range[isys][i]->SetDirectory(0);
	  h_corr_reco_dphi_range[isys][i]->SetDirectory(0);
	  h_corr_truth_dphi_range[isys][i]->SetDirectory(0);

	  h_all_reco_dphi_rangemin[isys][i]->SetDirectory(0);
	  h_all_truth_dphi_rangemin[isys][i]->SetDirectory(0);
	  h_corr_reco_dphi_rangemin[isys][i]->SetDirectory(0);
	  h_corr_truth_dphi_rangemin[isys][i]->SetDirectory(0);
	}
      
      h_data_pt1pt2dphi[isys]->SetDirectory(0);//(Form("h_data_pt1pt2dphi_%d", isys));
      h_reco_match_pt1pt2dphi[isys]->SetDirectory(0);//(Form("h_reco_pt1pt2dphi_%d", isys));
      h_truth_match_pt1pt2dphi[isys]->SetDirectory(0);//(Form("h_reco_pt1pt2dphi_%d", isys));
      h_truth_pt1pt2dphi[isys]->SetDirectory(0);//(Form("h_reco_pt1pt2dphi_%d", isys));
      h_reco_match_pt1pt2dphitruth[isys]->SetDirectory(0);//(Form("h_reco_pt1pt2dphitruth_%d", isys));
            
      // corrections are truth / reco
      for (int i = 0; i< 3; i++)
	{
	        
	  h_data_dphi_range[isys][i] = h_data_pt1pt2dphi[isys]->ProjectionZ(Form("h_data_dphi_range_%d_%d", isys, i), binranges[i], binranges[i+1], measure_subleading_bin, -1);

	  h_truth_dphi_range[isys][i] = h_truth_pt1pt2dphi[isys]->ProjectionZ(Form("h_truth_dphi_range_%d_%d", isys, i), binranges[i], binranges[i+1], measure_subleading_bin, -1);

	  h_truth_match_dphi_range[isys][i] = h_truth_match_pt1pt2dphi[isys]->ProjectionZ(Form("h_truth_match_dphi_range_%d_%d", isys, i), binranges[i], binranges[i+1], measure_subleading_bin, -1);

	  h_reco_dphi_range[isys][i] = h_reco_match_pt1pt2dphi[isys]->ProjectionZ(Form("h_reco_dphi_range_%d_%d", isys, i), binranges[i], binranges[i+1], measure_subleading_bin, -1);
	        
	  h_correction_factors_range[isys][i] = (TH1D*) h_truth_match_dphi_range[isys][i]->Clone();
	  h_correction_factors_range[isys][i]->SetName(Form("h_correction_factors_range_%d_%d", isys, i));
	  h_correction_factors_range[isys][i]->Divide( h_reco_dphi_range[isys][i] );

	  for (int ib = 0; ib < h_correction_factors_range[isys][i]->GetNbinsX(); ib++)
	    {
	      float cov = h_corr_reco_dphi_range[isys][i]->GetBinContent(ib+1);
	      float err1 = 0;
	      err1 += TMath::Power(h_reco_dphi_range[isys][i]->GetBinError(ib+1)/h_reco_dphi_range[isys][i]->GetBinContent(ib+1), 2);
	      err1 += TMath::Power(h_truth_match_dphi_range[isys][i]->GetBinError(ib+1)/h_truth_match_dphi_range[isys][i]->GetBinContent(ib+1), 2);
	      err1 -= 2.*cov/(h_reco_dphi_range[isys][i]->GetBinContent(ib+1)*h_truth_match_dphi_range[isys][i]->GetBinContent(ib+1));
	      float errr = h_correction_factors_range[isys][i]->GetBinContent(ib+1)*err1;

	      h_correction_factors_range[isys][i]->SetBinError(ib+1, errr);
	    }
	  
	  
	  h_efficiency_factors_range[isys][i] = (TH1D*) h_truth_dphi_range[isys][i]->Clone();
	  h_efficiency_factors_range[isys][i]->SetName(Form("h_efficiency_factors_range_%d_%d", isys, i));
	  h_efficiency_factors_range[isys][i]->Divide( h_truth_match_dphi_range[isys][i] );

	  if (isys == 0)
	    {
	      h_data_dphi_range[5][i] = h_data_pt1pt2dphi[isys]->ProjectionZ(Form("h_data_dphi_range_sys_half_%d_%d", isys, i), binranges[i], binranges[i+1], measure_subleading_bin, -1);
	      h_data_dphi_range[6][i] = h_data_pt1pt2dphi[isys]->ProjectionZ(Form("h_data_dphi_range_sys_double_%d_%d", isys, i), binranges[i], binranges[i+1], measure_subleading_bin, -1);
	    }
	        
	  for (int ib = 0; ib < h_truth_dphi_range[isys][i]->GetNbinsX(); ib++)
	    {
	      float va = h_data_dphi_range[isys][i]->GetBinContent(ib+1);
	      float vb = h_correction_factors_range[isys][i]->GetBinContent(ib+1);
	      float ea = h_data_dphi_range[isys][i]->GetBinError(ib+1);
	      float eb = h_correction_factors_range[isys][i]->GetBinError(ib+1);
	      float vv = va*vb;
	      float new_err = vv*sqrt( pow(ea/va, 2) + pow(eb/vb, 2));
	      h_data_dphi_range[isys][i]->SetBinContent(ib+1, vv);
	      h_data_dphi_range[isys][i]->SetBinError(ib+1, new_err);
	      if (isys == 0)
		{
		  float half_fac = ( h_correction_factors_range[isys][i]->GetBinContent(ib+1) - 1)*0.5 + 1;
		  float double_fac = ( h_correction_factors_range[isys][i]->GetBinContent(ib+1) - 1)*1.5 + 1;

		  float va = h_data_dphi_range[5][i]->GetBinContent(ib+1);
		  float ea = h_data_dphi_range[isys][i]->GetBinError(ib+1);
		  float eb = (h_correction_factors_range[isys][i]->GetBinError(ib+1)*0.5);
		  float ec = (h_correction_factors_range[isys][i]->GetBinError(ib+1)*1.5);
		  float vvh = va*half_fac;
		  float vvd = va*double_fac;
		  float new_err_h = vvh*sqrt( pow(ea/va, 2) + pow(eb/half_fac, 2));
		  float new_err_d = vvd*sqrt( pow(ea/va, 2) + pow(eb/double_fac, 2));

		  h_data_dphi_range[5][i]->SetBinContent(ib+1, vvh);
		  h_data_dphi_range[5][i]->SetBinError(ib+1, new_err_h);
		  h_data_dphi_range[6][i]->SetBinContent(ib+1, vvd);
		  h_data_dphi_range[6][i]->SetBinError(ib+1, new_err_d);
		}
	    }
	        
	  h_data_dphi_range[isys][i]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
	  h_truth_dphi_range[isys][i]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
	  h_truth_match_dphi_range[isys][i]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
	  h_reco_dphi_range[isys][i]->SetDirectory(0);//(Form("h_reco_match_dphi_%d", isys));
	  //h_data_dphi_range[isys][i]->Rebin(2);
	  //h_truth_dphi_range[isys][i]->Rebin(2);
	  //h_truth_match_dphi_range[isys][i]->Rebin(2);
	  //h_reco_dphi_range[isys][i]->Rebin(2);
	  if (isys == 0)
	    {
	      //h_data_dphi_range[5][i]->Rebin(2);
	      h_data_dphi_range[5][i]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
	      //h_data_dphi_range[6][i]->Rebin(2);
	      h_data_dphi_range[6][i]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
	    }
      

	  if (i == 1)
	    {
	      for (int j = 0; j < 3; j++)
		{
		        
		  h_data_dphi_rangemin[isys][j] = h_data_pt1pt2dphi[isys]->ProjectionZ(Form("h_data_dphi_rangemin_%d_%d", isys, j), binranges[i], binranges[i+1], binrangesmin[j], binrangesmin[j+1]);
		  h_truth_dphi_rangemin[isys][j] = h_truth_pt1pt2dphi[isys]->ProjectionZ(Form("h_truth_dphi_rangemin_%d_%d", isys, j), binranges[i], binranges[i+1], binrangesmin[j], binrangesmin[j+1]);
		  h_truth_match_dphi_rangemin[isys][j] = h_truth_match_pt1pt2dphi[isys]->ProjectionZ(Form("h_truth_match_dphi_rangemin_%d_%d", isys, j), binranges[i], binranges[i+1], binrangesmin[j], binrangesmin[j+1]);
		  h_reco_dphi_rangemin[isys][j] = h_reco_match_pt1pt2dphi[isys]->ProjectionZ(Form("h_reco_dphi_rangemin_%d_%d", isys, j), binranges[i], binranges[i+1], binrangesmin[j], binrangesmin[j+1]);

		  h_correction_factors_rangemin[isys][j] = (TH1D*) h_truth_match_dphi_rangemin[isys][j]->Clone();
		  h_correction_factors_rangemin[isys][j]->SetName(Form("h_correction_factors_rangemin_%d_%d", isys, j));
		  h_correction_factors_rangemin[isys][j]->Divide( h_reco_dphi_rangemin[isys][j] );

		  for (int ib = 0; ib < h_correction_factors_rangemin[isys][j]->GetNbinsX(); ib++)
		    {
		      float cov = h_corr_reco_dphi_rangemin[isys][j]->GetBinContent(ib+1);
		      float err1 = 0;
		      err1 += TMath::Power(h_reco_dphi_rangemin[isys][j]->GetBinError(ib+1)/h_reco_dphi_rangemin[isys][j]->GetBinContent(ib+1), 2);
		      err1 += TMath::Power(h_truth_match_dphi_rangemin[isys][j]->GetBinError(ib+1)/h_truth_match_dphi_rangemin[isys][j]->GetBinContent(ib+1), 2);
		      err1 -= 2.*cov/(h_reco_dphi_rangemin[isys][j]->GetBinContent(ib+1)*h_truth_match_dphi_rangemin[isys][j]->GetBinContent(ib+1));
		      float errr = h_correction_factors_rangemin[isys][j]->GetBinContent(ib+1)*err1;

		      h_correction_factors_rangemin[isys][j]->SetBinError(ib+1, errr);
		    }

		  h_efficiency_factors_rangemin[isys][j] = (TH1D*) h_truth_dphi_rangemin[isys][j]->Clone();
		  h_efficiency_factors_rangemin[isys][j]->SetName(Form("h_efficiency_factors_rangemin_%d_%d", isys, j));
		  h_efficiency_factors_rangemin[isys][j]->Divide( h_truth_match_dphi_rangemin[isys][j] );

		  if (isys == 0)
		    {
		      h_data_dphi_rangemin[5][j] = h_data_pt1pt2dphi[isys]->ProjectionZ(Form("h_data_dphi_rangemin_sys_half_%d_%d", isys, j), binranges[i], binranges[i+1],  binrangesmin[j], binrangesmin[j+1]);//measure_subleading_bin, -1);
		      h_data_dphi_rangemin[6][j] = h_data_pt1pt2dphi[isys]->ProjectionZ(Form("h_data_dphi_rangemin_sys_double_%d_%d", isys, j), binranges[i], binranges[i+1],  binrangesmin[j], binrangesmin[j+1]);
		    }

		  for (int ib = 0; ib < h_truth_dphi_rangemin[isys][j]->GetNbinsX(); ib++)
		    {
		      float va = h_data_dphi_rangemin[isys][j]->GetBinContent(ib+1);
		      float vb = h_correction_factors_rangemin[isys][j]->GetBinContent(ib+1);
		      float ea = h_data_dphi_rangemin[isys][j]->GetBinError(ib+1);
		      float eb = h_correction_factors_rangemin[isys][j]->GetBinError(ib+1);
		      float vv = va*vb;
		      float new_err = vv*sqrt( pow(ea/va, 2) + pow(eb/vb, 2));
		      h_data_dphi_rangemin[isys][j]->SetBinContent(ib+1, vv);
		      h_data_dphi_rangemin[isys][j]->SetBinError(ib+1, new_err);
		      if (isys == 0)
			{
			  float half_fac = ( h_correction_factors_rangemin[isys][j]->GetBinContent(ib+1) - 1)*0.5 + 1;
			  float double_fac = ( h_correction_factors_rangemin[isys][j]->GetBinContent(ib+1) - 1)*1.5 + 1;

			  float va = h_data_dphi_rangemin[5][j]->GetBinContent(ib+1);
			  float ea = h_data_dphi_rangemin[isys][j]->GetBinError(ib+1);
			  float eb = (h_correction_factors_rangemin[isys][j]->GetBinError(ib+1)*0.5);
			  float ec = (h_correction_factors_rangemin[isys][j]->GetBinError(ib+1)*1.5);
			  float vvh = va*half_fac;
			  float vvd = va*double_fac;
			  float new_err_h = vvh*sqrt( pow(ea/va, 2) + pow(eb/half_fac, 2));
			  float new_err_d = vvd*sqrt( pow(ea/va, 2) + pow(eb/double_fac, 2));

			  h_data_dphi_rangemin[5][j]->SetBinContent(ib+1, vvh);
			  h_data_dphi_rangemin[5][j]->SetBinError(ib+1, new_err_h);
			  h_data_dphi_rangemin[6][j]->SetBinContent(ib+1, vvd);
			  h_data_dphi_rangemin[6][j]->SetBinError(ib+1, new_err_d);
			}
		    }
		  //h_data_dphi_rangemin[isys][j]->Rebin(2);//(1./h_data_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  //h_truth_dphi_rangemin[isys][j]->Rebin(2);//(1./h_truth_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  //h_truth_match_dphi_rangemin[isys][j]->Rebin(2);//(1./h_truth_match_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  //h_reco_dphi_rangemin[isys][j]->Rebin(2);//(1./h_reco_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  
		  scale_pt2_truth[isys][j] = 1./h_truth_match_dphi_rangemin[isys][j]->Integral(0, -1, "width");
		  scale_pt2_reco[isys][j] = 1./h_reco_dphi_rangemin[isys][j]->Integral(0, -1, "width");
		  
		  h_data_dphi_rangemin[isys][j]->Scale(1./h_data_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  h_truth_dphi_rangemin[isys][j]->Scale(1./h_truth_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  h_truth_match_dphi_rangemin[isys][j]->Scale(1./h_truth_match_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  h_reco_dphi_rangemin[isys][j]->Scale(1./h_reco_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  
		  h_data_dphi_rangemin[isys][j]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
		  h_truth_dphi_rangemin[isys][j]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
		  h_truth_match_dphi_rangemin[isys][j]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
		  h_reco_dphi_rangemin[isys][j]->SetDirectory(0);//(Form("h_reco_match_dphi_%d", isys));

		  if (isys == 0)
		    {
		      //h_data_dphi_rangemin[5][j]->Rebin(2);
		      h_data_dphi_rangemin[5][j]->Scale(1./h_data_dphi_rangemin[5][j]->Integral(0, -1, "width"));
		      h_data_dphi_rangemin[5][j]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
		      //h_data_dphi_rangemin[6][j]->Rebin(2);
		      h_data_dphi_rangemin[6][j]->Scale(1./h_data_dphi_rangemin[6][j]->Integral(0, -1, "width"));
		      h_data_dphi_rangemin[6][j]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
		    }
		}
	    }

	  scale_pt1_truth[isys][i] = 1./h_truth_match_dphi_range[isys][i]->Integral(0, -1, "width");
	  scale_pt1_reco[isys][i] = 1./h_reco_dphi_range[isys][i]->Integral(0, -1, "width");

	  h_data_dphi_range[isys][i]->Scale(1./h_data_dphi_range[isys][i]->Integral(0, -1, "width"));
	  h_truth_dphi_range[isys][i]->Scale(1./h_truth_dphi_range[isys][i]->Integral(0, -1, "width"));
	  h_truth_match_dphi_range[isys][i]->Scale(1./h_truth_match_dphi_range[isys][i]->Integral(0, -1, "width"));
	  h_reco_dphi_range[isys][i]->Scale(1./h_reco_dphi_range[isys][i]->Integral(0, -1, "width"));
	  if (isys == 0)
	    {
	      h_data_dphi_range[5][i]->Scale(1./h_data_dphi_range[5][i]->Integral(0, -1, "width"));
	      h_data_dphi_range[6][i]->Scale(1./h_data_dphi_range[6][i]->Integral(0, -1, "width"));
	    }
      
	  std::cout << isys << " " << i << " Mean  "<< h_data_dphi_range[isys][i]->GetMean() << std::endl;
	}
      
    }


  TH1D *h_jer_sys_range[2][3];
  TH1D *h_jes_sys_range[2][3];
  TH1D *h_corr_sys_range[2][3];
  TH1D *h_all_sys_range[2][3];
  
  TH1D *h_jer_sys_rangemin[2][3];
  TH1D *h_jes_sys_rangemin[2][3];
  TH1D *h_corr_sys_rangemin[2][3];
  TH1D *h_all_sys_rangemin[2][3];
  

  // make two copies, on with sat and one with sys
  TGraphAsymmErrors *h_data_dphi_range_sys[3];
  TGraphAsymmErrors *h_data_dphi_rangemin_sys[3];

  TGraphAsymmErrors *h_ratio_dphi_range_sys[3];
  TGraphAsymmErrors *h_ratio_dphi_rangemin_sys[3];

  for (int i = 0 ; i < 3; i++)
    {
      h_data_dphi_range_sys[i] = new TGraphAsymmErrors(h_data_dphi_range[0][i]);
      h_data_dphi_range_sys[i]->SetName(Form("h_data_dphi_range_sys_%d", i));
      h_ratio_dphi_range_sys[i] = new TGraphAsymmErrors(h_data_dphi_range[0][i]);
      h_ratio_dphi_range_sys[i]->SetName(Form("h_ratio_dphi_range_sys_%d", i));

      h_jer_sys_range[0][i] = (TH1D*) h_data_dphi_range[0][i]->Clone();
      h_jer_sys_range[0][i]->SetName(Form("h_jer_sys_range_pos_%d", i));
      h_jer_sys_range[0][i]->Reset();//(Form("h_jer_sys_range_pos_%d", i));
      h_jes_sys_range[0][i] = (TH1D*) h_data_dphi_range[0][i]->Clone();
      h_jes_sys_range[0][i]->SetName(Form("h_jes_sys_range_pos_%d", i));
      h_jes_sys_range[0][i]->Reset();//(Form("h_jes_sys_range_pos_%d", i));
      h_corr_sys_range[0][i] = (TH1D*) h_data_dphi_range[0][i]->Clone();
      h_corr_sys_range[0][i]->SetName(Form("h_corr_sys_range_pos_%d", i));
      h_corr_sys_range[0][i]->Reset();//(Form("h_corr_sys_range_pos_%d", i));
      h_all_sys_range[0][i] = (TH1D*) h_data_dphi_range[0][i]->Clone();
      h_all_sys_range[0][i]->SetName(Form("h_all_sys_range_pos_%d", i));
      h_all_sys_range[0][i]->Reset();//(Form("h_all_sys_range_pos_%d", i));

      h_jer_sys_range[1][i] = (TH1D*) h_data_dphi_range[0][i]->Clone();
      h_jer_sys_range[1][i]->SetName(Form("h_jer_sys_range_neg_%d", i));
      h_jer_sys_range[1][i]->Reset();//(Form("h_jer_sys_range_neg_%d", i));
      h_jes_sys_range[1][i] = (TH1D*) h_data_dphi_range[0][i]->Clone();
      h_jes_sys_range[1][i]->SetName(Form("h_jes_sys_range_neg_%d", i));
      h_jes_sys_range[1][i]->Reset();//(Form("h_jes_sys_range_neg_%d", i));
      h_corr_sys_range[1][i] = (TH1D*) h_data_dphi_range[0][i]->Clone();
      h_corr_sys_range[1][i]->SetName(Form("h_corr_sys_range_neg_%d", i));
      h_corr_sys_range[1][i]->Reset();//(Form("h_corr_sys_range_neg_%d", i));

      h_all_sys_range[1][i] = (TH1D*) h_data_dphi_range[0][i]->Clone();
      h_all_sys_range[1][i]->SetName(Form("h_all_sys_range_neg_%d", i));
      h_all_sys_range[1][i]->Reset();//(Form("h_all_sys_range_neg_%d", i));

      std::cout << "range: " << i << std::endl;
      for (int ib = 0; ib < h_data_dphi_range[0][i]->GetNbinsX(); ib++)
	{
	  if (h_data_dphi_range[0][i]->GetBinContent(ib+1) == 0) continue;
	  // Correction
	  float err1 = h_data_dphi_range[0][i]->GetBinContent(ib+1) - h_data_dphi_range[5][i]->GetBinContent(ib+1);
	  float err2 = h_data_dphi_range[0][i]->GetBinContent(ib+1) - h_data_dphi_range[6][i]->GetBinContent(ib+1); 
	  std::pair corr_pair = get_ordered_pair(err1, err2);

	  h_corr_sys_range[0][i]->SetBinContent(ib+1, corr_pair.first/h_data_dphi_range[0][i]->GetBinContent(ib+1));
	  h_corr_sys_range[1][i]->SetBinContent(ib+1, -1*corr_pair.second/h_data_dphi_range[0][i]->GetBinContent(ib+1));
	  // JES
	  err1 = h_data_dphi_range[0][i]->GetBinContent(ib+1) - h_data_dphi_range[1][i]->GetBinContent(ib+1);
	  err2 = h_data_dphi_range[0][i]->GetBinContent(ib+1) - h_data_dphi_range[2][i]->GetBinContent(ib+1); 
	  std::pair jes_pair = get_ordered_pair(err1, err2);
	  h_jes_sys_range[0][i]->SetBinContent(ib+1, jes_pair.first/h_data_dphi_range[0][i]->GetBinContent(ib+1));
	  h_jes_sys_range[1][i]->SetBinContent(ib+1, -1*jes_pair.second/h_data_dphi_range[0][i]->GetBinContent(ib+1));

	  // JER
	  err1 = h_data_dphi_range[0][i]->GetBinContent(ib+1) - h_data_dphi_range[3][i]->GetBinContent(ib+1);
	  err2 = h_data_dphi_range[0][i]->GetBinContent(ib+1) - h_data_dphi_range[4][i]->GetBinContent(ib+1); 
	  
	  std::pair jer_pair = get_ordered_pair(err1, err2);
	  h_jer_sys_range[0][i]->SetBinContent(ib+1, jer_pair.first/h_data_dphi_range[0][i]->GetBinContent(ib+1));
	  h_jer_sys_range[1][i]->SetBinContent(ib+1, -1*jer_pair.second/h_data_dphi_range[0][i]->GetBinContent(ib+1));


	  float sys_pos = sqrt(corr_pair.first*corr_pair.first +
			       jes_pair.first*jes_pair.first +
			       jer_pair.first*jer_pair.first);
	  float sys_neg = sqrt(corr_pair.second*corr_pair.second +
			       jes_pair.second*jes_pair.second +
			       jer_pair.second*jer_pair.second);
	  
	  h_all_sys_range[0][i]->SetBinContent(ib+1, sys_pos/h_data_dphi_range[0][i]->GetBinContent(ib+1));
	  h_all_sys_range[1][i]->SetBinContent(ib+1, -1*sys_neg/h_data_dphi_range[0][i]->GetBinContent(ib+1));

	  
	  h_data_dphi_range_sys[i]->SetPointError(ib, h_data_dphi_range[0][i]->GetBinWidth(ib+1)/2., h_data_dphi_range[0][i]->GetBinWidth(ib+1)/2., fabs(sys_neg), fabs(sys_pos));
	  h_ratio_dphi_range_sys[i]->SetPointY(ib, 1);//h_data_dphi_range[0][i]->GetBinContent(ib+1)/h_truth_match_dphi_range[0][i]->GetBinContent(ib+1));
	  h_ratio_dphi_range_sys[i]->SetPointError(ib, h_data_dphi_range[0][i]->GetBinWidth(ib+1)/2., h_data_dphi_range[0][i]->GetBinWidth(ib+1)/2., fabs(sys_neg)/h_data_dphi_range[0][i]->GetBinContent(ib+1), fabs(sys_pos)/h_data_dphi_range[0][i]->GetBinContent(ib+1));

	}

      h_data_dphi_rangemin_sys[i] = new TGraphAsymmErrors(h_data_dphi_rangemin[0][i]);
      h_data_dphi_rangemin_sys[i]->SetName(Form("h_data_dphi_rangemin_sys_%d", i));
      h_ratio_dphi_rangemin_sys[i] = new TGraphAsymmErrors(h_data_dphi_rangemin[0][i]);
      h_ratio_dphi_rangemin_sys[i]->SetName(Form("h_ratio_dphi_rangemin_sys_%d", i));

      h_jer_sys_rangemin[0][i] = (TH1D*) h_data_dphi_rangemin[0][i]->Clone();
      h_jer_sys_rangemin[0][i]->SetName(Form("h_jer_sys_rangemin_pos_%d", i));
      h_jer_sys_rangemin[0][i]->Reset();//(Form("h_jer_sys_rangemin_pos_%d", i));
      h_jes_sys_rangemin[0][i] = (TH1D*) h_data_dphi_rangemin[0][i]->Clone();
      h_jes_sys_rangemin[0][i]->SetName(Form("h_jes_sys_rangemin_pos_%d", i));
      h_jes_sys_rangemin[0][i]->Reset();//(Form("h_jes_sys_rangemin_pos_%d", i));
      h_corr_sys_rangemin[0][i] = (TH1D*) h_data_dphi_rangemin[0][i]->Clone();
      h_corr_sys_rangemin[0][i]->SetName(Form("h_corr_sys_rangemin_pos_%d", i));
      h_corr_sys_rangemin[0][i]->Reset();//(Form("h_corr_sys_rangemin_pos_%d", i));

      h_all_sys_rangemin[0][i] = (TH1D*) h_data_dphi_rangemin[0][i]->Clone();
      h_all_sys_rangemin[0][i]->SetName(Form("h_all_sys_rangemin_pos_%d", i));
      h_all_sys_rangemin[0][i]->Reset();//(Form("h_all_sys_rangemin_pos_%d", i));

      h_jer_sys_rangemin[1][i] = (TH1D*) h_data_dphi_rangemin[0][i]->Clone();
      h_jer_sys_rangemin[1][i]->SetName(Form("h_jer_sys_rangemin_neg_%d", i));
      h_jer_sys_rangemin[1][i]->Reset();//(Form("h_jer_sys_rangemin_neg_%d", i));
      h_jes_sys_rangemin[1][i] = (TH1D*) h_data_dphi_rangemin[0][i]->Clone();
      h_jes_sys_rangemin[1][i]->SetName(Form("h_jes_sys_rangemin_neg_%d", i));
      h_jes_sys_rangemin[1][i]->Reset();//(Form("h_jes_sys_rangemin_neg_%d", i));
      h_corr_sys_rangemin[1][i] = (TH1D*) h_data_dphi_rangemin[0][i]->Clone();
      h_corr_sys_rangemin[1][i]->SetName(Form("h_corr_sys_rangemin_neg_%d", i));
      h_corr_sys_rangemin[1][i]->Reset();//(Form("h_corr_sys_rangemin_neg_%d", i));

      h_all_sys_rangemin[1][i] = (TH1D*) h_data_dphi_rangemin[0][i]->Clone();
      h_all_sys_rangemin[1][i]->SetName(Form("h_all_sys_rangemin_neg_%d", i));
      h_all_sys_rangemin[1][i]->Reset();//(Form("h_all_sys_rangemin_neg_%d", i));

      std::cout << "rangemin: " << i << std::endl;
      for (int ib = 0; ib < h_data_dphi_rangemin[0][i]->GetNbinsX(); ib++)
	{

	  if (h_data_dphi_rangemin[0][i]->GetBinContent(ib+1) == 0) continue;
	  std::cout << ib << std::endl;
	  // Correction
	  float err1 = h_data_dphi_rangemin[0][i]->GetBinContent(ib+1) - h_data_dphi_rangemin[5][i]->GetBinContent(ib+1);
	  float err2 = h_data_dphi_rangemin[0][i]->GetBinContent(ib+1) - h_data_dphi_rangemin[6][i]->GetBinContent(ib+1); 

	  std::pair corr_pair = get_ordered_pair(err1, err2);

	  h_corr_sys_rangemin[0][i]->SetBinContent(ib+1, corr_pair.first/h_data_dphi_rangemin[0][i]->GetBinContent(ib+1));
	  h_corr_sys_rangemin[1][i]->SetBinContent(ib+1, -1*corr_pair.second/h_data_dphi_rangemin[0][i]->GetBinContent(ib+1));

	  // JES
	  err1 = h_data_dphi_rangemin[0][i]->GetBinContent(ib+1) - h_data_dphi_rangemin[1][i]->GetBinContent(ib+1);
	  err2 = h_data_dphi_rangemin[0][i]->GetBinContent(ib+1) - h_data_dphi_rangemin[2][i]->GetBinContent(ib+1); 
	  std::pair jes_pair = get_ordered_pair(err1, err2);
	  h_jes_sys_rangemin[0][i]->SetBinContent(ib+1, jes_pair.first/h_data_dphi_rangemin[0][i]->GetBinContent(ib+1));
	  h_jes_sys_rangemin[1][i]->SetBinContent(ib+1, -1*jes_pair.second/h_data_dphi_rangemin[0][i]->GetBinContent(ib+1));

	  // JER
	  err1 = h_data_dphi_rangemin[0][i]->GetBinContent(ib+1) - h_data_dphi_rangemin[3][i]->GetBinContent(ib+1);
	  err2 = h_data_dphi_rangemin[0][i]->GetBinContent(ib+1) - h_data_dphi_rangemin[4][i]->GetBinContent(ib+1); 
	  
	  std::pair jer_pair = get_ordered_pair(err1, err2);
	  h_jer_sys_rangemin[0][i]->SetBinContent(ib+1, jer_pair.first/h_data_dphi_rangemin[0][i]->GetBinContent(ib+1));
	  h_jer_sys_rangemin[1][i]->SetBinContent(ib+1, -1*jer_pair.second/h_data_dphi_rangemin[0][i]->GetBinContent(ib+1));

	  std::cout << corr_pair.first << " - " <<corr_pair.second << std::endl;
	  std::cout << jes_pair.first << " - " <<jes_pair.second << std::endl;
	  std::cout << jer_pair.first << " - " <<jer_pair.second << std::endl;
	  float sys_pos = sqrt(corr_pair.first*corr_pair.first +
			       jes_pair.first*jes_pair.first +
			       jer_pair.first*jer_pair.first);
	  float sys_neg = sqrt(corr_pair.second*corr_pair.second +
			       jes_pair.second*jes_pair.second +
			       jer_pair.second*jer_pair.second);
	  std::cout << sys_pos << " / " << sys_neg << std::endl;

	  h_all_sys_rangemin[0][i]->SetBinContent(ib+1, sys_pos/h_data_dphi_rangemin[0][i]->GetBinContent(ib+1));
	  h_all_sys_rangemin[1][i]->SetBinContent(ib+1, -1*sys_neg/h_data_dphi_rangemin[0][i]->GetBinContent(ib+1));

	  h_data_dphi_rangemin_sys[i]->SetPointError(ib, h_data_dphi_rangemin[0][i]->GetBinWidth(ib+1)/2., h_data_dphi_rangemin[0][i]->GetBinWidth(ib+1)/2., fabs(sys_neg), fabs(sys_pos));
	  h_ratio_dphi_rangemin_sys[i]->SetPointY(ib, 1);//h_data_dphi_rangemin[0][i]->GetBinContent(ib+1)/h_truth_match_dphi_rangemin[0][i]->GetBinContent(ib+1));
	  h_ratio_dphi_rangemin_sys[i]->SetPointError(ib, h_data_dphi_rangemin[0][i]->GetBinWidth(ib+1)/2., h_data_dphi_rangemin[0][i]->GetBinWidth(ib+1)/2., fabs(sys_neg)/h_data_dphi_rangemin[0][i]->GetBinContent(ib+1), fabs(sys_pos)/h_data_dphi_rangemin[0][i]->GetBinContent(ib+1));

	}
    }
  std::cout << __LINE__ << std::endl;
  TCanvas *c_sysdiv = new TCanvas("c_sysdiv", "c_sysdiv", 500, 500);
  for (int i = 0; i < 3; i++)
    {

      for (int j = 0 ; j < 2 ; j++)
	{
	  dlutility::SetLineAtt(h_all_sys_range[j][i], kBlack, 3, 1);
	  dlutility::SetLineAtt(h_corr_sys_range[j][i], colorsys[1], 3, 1);
	  dlutility::SetLineAtt(h_jer_sys_range[j][i], colorsys[2], 3, 1);
	  dlutility::SetLineAtt(h_jes_sys_range[j][i], colorsys[3], 3, 1);
	  
	  dlutility::SetMarkerAtt(h_all_sys_range[j][i], kBlack, 1, 8);
	}
      h_all_sys_range[0][i]->SetTitle(";#Delta#phi; Systematic Error;");
      h_all_sys_range[0][i]->SetMaximum(1.2);;
      h_all_sys_range[0][i]->SetMinimum(-0.6);;
      dlutility::SetFont(h_all_sys_range[0][i], 42, 0.06, 0.05, 0.05, 0.05);
      h_all_sys_range[0][i]->GetYaxis()->SetTitleOffset(1.5);
      h_all_sys_range[0][i]->Draw("p");
      h_all_sys_range[1][i]->Draw("p same");

      h_corr_sys_range[0][i]->Draw("hist same");
      h_corr_sys_range[1][i]->Draw("hist same");
      h_jer_sys_range[0][i]->Draw("hist same");
      h_jer_sys_range[1][i]->Draw("hist same");
      h_jes_sys_range[0][i]->Draw("hist same");
      h_jes_sys_range[1][i]->Draw("hist same");
      dlutility::DrawSPHENIXppPrelim(0.22, 0.86);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.2]{#it{R} = 0.4}", 0.22, 0.76);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,1} < %2.1f GeV}", ipt_bins[binranges[i]], ipt_bins[binranges[i+1]]), 0.22, 0.71);
      dlutility::drawText(Form("p_{T,2} #kern[-0.08]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, 0.66);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, 0.59);


      TLegend *leg = new TLegend(0.61, 0.72, 0.75, 0.88);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.037);
      leg->SetTextFont(42);
      leg->AddEntry(h_all_sys_range[0][i], "All Systematics","p");
      leg->AddEntry(h_jer_sys_range[0][i], "JER");
      leg->AddEntry(h_jes_sys_range[0][i], "JES");
      leg->AddEntry(h_corr_sys_range[0][i], "Reco-correction");
	    

      leg->Draw("same");

      c_sysdiv->Print(Form("h_sys_all_range_%d.pdf", i));
      c_sysdiv->Print(Form("h_sys_all_range_%d.png", i));
    }
  TCanvas *c_sysdiv2 = new TCanvas("c_sysdiv2", "c_sysdiv2", 500, 500);
  for (int i = 0; i < 3; i++)
    {
      for (int j = 0 ; j < 2 ; j++)
	{
	  dlutility::SetLineAtt(h_all_sys_rangemin[j][i], kBlack, 3, 1);
	  dlutility::SetLineAtt(h_corr_sys_rangemin[j][i], colorsys[1], 3, 1);
	  dlutility::SetLineAtt(h_jer_sys_rangemin[j][i], colorsys[2], 3, 1);
	  dlutility::SetLineAtt(h_jes_sys_rangemin[j][i], colorsys[3], 3, 1);
	  
	  dlutility::SetMarkerAtt(h_all_sys_rangemin[j][i], kBlack, 1, 8);
	}
      h_all_sys_rangemin[0][i]->SetTitle(";#Delta#phi; Systematic Error;");
      h_all_sys_rangemin[0][i]->SetMaximum(1.2);;
      h_all_sys_rangemin[0][i]->SetMinimum(-0.6);;
      dlutility::SetFont(h_all_sys_rangemin[0][i], 42, 0.06, 0.05, 0.05, 0.05);
      h_all_sys_rangemin[0][i]->GetYaxis()->SetTitleOffset(1.5);

      std::cout << "bin 10: " << h_all_sys_rangemin[0][i]->GetBinContent(10) << std::endl;;//SetMaximum(2.0);;
      std::cout << "bin 10: " << h_all_sys_rangemin[1][i]->GetBinContent(10) << std::endl;;//SetMaximum(2.0);;
      h_all_sys_rangemin[0][i]->Draw("p");
      h_all_sys_rangemin[1][i]->Draw("p same");

      h_corr_sys_rangemin[0][i]->Draw("hist same");
      h_corr_sys_rangemin[1][i]->Draw("hist same");
      h_jer_sys_rangemin[0][i]->Draw("hist same");
      h_jer_sys_rangemin[1][i]->Draw("hist same");
      h_jes_sys_rangemin[0][i]->Draw("hist same");
      h_jes_sys_rangemin[1][i]->Draw("hist same");

      dlutility::DrawSPHENIXppPrelim(0.22, 0.86);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.2]{#it{R} = 0.4}", 0.22, 0.76);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,1} < %2.1f GeV}", ipt_bins[binranges[1]], ipt_bins[binranges[1+1]]), 0.22, 0.71);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,2} < %2.1f GeV}", ipt_bins[binrangesmin[i]], ipt_bins[binrangesmin[i+1]]), 0.22, 0.66);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, 0.61);


      TLegend *leg = new TLegend(0.61, 0.72, 0.75, 0.88);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_all_sys_rangemin[0][i], "All Systematics","p");
      leg->AddEntry(h_jer_sys_rangemin[0][i], "JER");
      leg->AddEntry(h_jes_sys_rangemin[0][i], "JES");
      leg->AddEntry(h_corr_sys_rangemin[0][i], "Reco-correction");
	    

      leg->Draw("same");

      c_sysdiv2->Print(Form("h_sys_all_rangemin_%d.pdf", i));
      c_sysdiv2->Print(Form("h_sys_all_rangemin_%d.png", i));
    }
  
  TCanvas *c_cor = new TCanvas("c_cor","c_cor", 500, 700);
  dlutility::ratioPanelCanvas(c_cor, 0.4);
  for (int i = 0; i < 3; i++)
    {

      c_cor->cd(1);
      gPad->SetBottomMargin(0.1);
      dlutility::SetLineAtt(h_truth_match_dphi_range[0][i], kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_truth_match_dphi_range[0][i], kBlack, 1, 8);
      dlutility::SetLineAtt(h_reco_dphi_range[0][i], kRed, 1, 1);
      dlutility::SetMarkerAtt(h_reco_dphi_range[0][i], kRed, 1, 8);

      float top = 0.8;
      float ss = 0.05;

      dlutility::SetFont(h_truth_match_dphi_range[0][i], 42, 0.05, 0.05, 0.05, 0.05);
      //h_truth_dphi_range[0][i]->Scale(1./ h_truth_dphi_range[0][i]->Integral(), "width");
      //h_reco_dphi_range[0][i]->Scale(1./ h_reco_dphi_range[0][i]->Integral(), "width");
      h_truth_match_dphi_range[0][i]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}");
      h_truth_match_dphi_range[0][i]->SetMaximum(10);
      h_truth_match_dphi_range[0][i]->SetMinimum(0);
      h_truth_match_dphi_range[0][i]->Draw("p");
      h_reco_dphi_range[0][i]->Draw("p same");

      dlutility::DrawSPHENIXppPrelimSpace(0.22, top , 0.1, 0, 1, 0, 1, "PYTHIA-8");
      
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R} = 0.4}", 0.22, top - 3*ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,1} < %2.1f GeV}", ipt_bins[binranges[i]], ipt_bins[binranges[i+1]]), 0.22, top - 4*ss);
      dlutility::drawText(Form("p_{T,2} #kern[-0.08]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, top - 5*ss);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, top - 6*ss);

      TLegend *leg1 = new TLegend(0.20, top - 9*ss, 0.63, top - 7*ss);
      leg1->SetLineWidth(0);
      leg1->SetTextSize(0.04);
      leg1->SetTextFont(42);
      leg1->AddEntry(h_truth_match_dphi_range[0][i], "Truth-particle","p");
      leg1->AddEntry(h_reco_dphi_range[0][i], "Reconstructed","p");
      leg1->Draw("same");  

      c_cor->cd(2);
      gPad->SetTopMargin(0.05);

      TH1D *hd = (TH1D*) h_truth_match_dphi_range[0][i]->Clone();
      hd->Divide(h_reco_dphi_range[0][i]);
      for (int ib = 0; ib < hd->GetNbinsX(); ib++)
	{
	  float cov = pow ( scale_pt1_reco[0][i]*scale_pt2_truth[0][i]*h_reco_dphi_range[0][i]->GetBinError(ib+1)*h_truth_match_dphi_range[0][i]->GetBinError(ib+1) , 2)/(h_corr_reco_dphi_range[0][i]->GetBinContent(ib+1)*scale_pt1_reco[0][i]*scale_pt2_truth[0][i]);

	    //h_all_truth_dphi_range[0][i]->GetBinContent(ib+1)*h_all_reco_dphi_range[0][i]->GetBinContent(ib+1)/h_corr_reco_dphi_range[0][i]->GetBinContent(ib+1);
	  float err1 = 0;
	  std::cout << cov << std::endl;
	  err1 += TMath::Power(h_reco_dphi_range[0][i]->GetBinError(ib+1)/h_reco_dphi_range[0][i]->GetBinContent(ib+1), 2);
	  err1 += TMath::Power(h_truth_match_dphi_range[0][i]->GetBinError(ib+1)/h_truth_match_dphi_range[0][i]->GetBinContent(ib+1), 2);
	  err1 -= 2.*cov/(h_reco_dphi_range[0][i]->GetBinContent(ib+1)*h_truth_match_dphi_range[0][i]->GetBinContent(ib+1));
	  err1 = sqrt(err1);
	  float errr = hd->GetBinContent(ib+1)*err1;
	  std::cout << ib << " : " << hd->GetBinContent(ib+1) << "  " << errr << std::endl;
	  hd->SetBinError(ib+1, errr);
	}
      hd->SetTitle(";#Delta#phi; Truth MC / Reco MC");
      dlutility::SetFont(hd, 42, 0.09, 0.08, 0.075, 0.075);
      hd->GetYaxis()->SetTitleOffset(0.8);
      hd->SetMaximum(2);
      hd->SetMinimum(0);
      hd->Draw("p");
      TLine *li1 = new TLine(TMath::Pi()*3./4., 1, TMath::Pi(), 1);
      dlutility::SetLineAtt(li1, kBlack, 1, 4);
      li1->Draw("same");

      c_cor->Print(Form("correction_range_%d.pdf", i));
      c_cor->Print(Form("correction_range_%d.png", i));
    }
  for (int i = 0; i < 3; i++)
    {
      c_cor->cd(1);
      gPad->SetBottomMargin(0.1);
      dlutility::SetLineAtt(h_truth_match_dphi_rangemin[0][i], kBlack, 1, 1);
      dlutility::SetMarkerAtt(h_truth_match_dphi_rangemin[0][i], kBlack, 1, 8);
      dlutility::SetLineAtt(h_reco_dphi_rangemin[0][i], kRed, 1, 1);
      dlutility::SetMarkerAtt(h_reco_dphi_rangemin[0][i], kRed, 1, 8);

      dlutility::SetFont(h_truth_match_dphi_rangemin[0][i], 42, 0.05);
      //h_truth_dphi_rangemin[0][i]->Scale(1./ h_truth_dphi_rangemin[0][i]->Integral(), "width");
      //h_reco_dphi_rangemin[0][i]->Scale(1./ h_reco_dphi_rangemin[0][i]->Integral(), "width");
      h_truth_match_dphi_rangemin[0][i]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}");
  
      h_truth_match_dphi_rangemin[0][i]->SetMaximum(10);
      h_truth_match_dphi_rangemin[0][i]->SetMinimum(0);
      h_truth_match_dphi_rangemin[0][i]->Draw("p");
      h_reco_dphi_rangemin[0][i]->Draw("p same");


      float top = 0.8;
      float ss = 0.05;


      dlutility::DrawSPHENIXppPrelimSpace(0.22, top , 0.1, 0, 1, 0, 1, "PYTHIA-8");
      
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R} = 0.4}", 0.22, top - 3*ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,1} < %2.1f GeV}", ipt_bins[binranges[i]], ipt_bins[binranges[i+1]]), 0.22, top - 4*ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,2} < %2.1f GeV}", ipt_bins[binrangesmin[i]], ipt_bins[binrangesmin[i+1]]), 0.22, top - 5*ss);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, top - 6*ss);

      TLegend *leg1 = new TLegend(0.20, top - 9*ss, 0.63, top - 7*ss);
      leg1->SetLineWidth(0);
      leg1->SetTextSize(0.04);
      leg1->SetTextFont(42);
      leg1->AddEntry(h_truth_match_dphi_range[0][i], "Truth-particle","p");
      leg1->AddEntry(h_reco_dphi_range[0][i], "Reconstructed","p");
      leg1->Draw("same");  


      c_cor->cd(2);
      gPad->SetTopMargin(0.05);

      TH1D *hd = (TH1D*) h_truth_match_dphi_rangemin[0][i]->Clone();
      hd->Divide(h_reco_dphi_rangemin[0][i]);
      for (int ib = 0; ib < hd->GetNbinsX(); ib++)
	{
	  float cov = pow ( scale_pt1_reco[0][i]*scale_pt2_truth[0][i]*h_reco_dphi_rangemin[0][i]->GetBinError(ib+1)*h_truth_match_dphi_rangemin[0][i]->GetBinError(ib+1) , 2)/(h_corr_reco_dphi_rangemin[0][i]->GetBinContent(ib+1)*scale_pt1_reco[0][i]*scale_pt2_truth[0][i]);

	    //h_all_truth_dphi_rangemin[0][i]->GetBinContent(ib+1)*h_all_reco_dphi_rangemin[0][i]->GetBinContent(ib+1)/h_corr_reco_dphi_rangemin[0][i]->GetBinContent(ib+1);
	  float err1 = 0;
	  std::cout << cov << std::endl;
	  err1 += TMath::Power(h_reco_dphi_rangemin[0][i]->GetBinError(ib+1)/h_reco_dphi_rangemin[0][i]->GetBinContent(ib+1), 2);
	  err1 += TMath::Power(h_truth_match_dphi_rangemin[0][i]->GetBinError(ib+1)/h_truth_match_dphi_rangemin[0][i]->GetBinContent(ib+1), 2);
	  err1 -= 2.*cov/(h_reco_dphi_rangemin[0][i]->GetBinContent(ib+1)*h_truth_match_dphi_rangemin[0][i]->GetBinContent(ib+1));
	  err1 = sqrt(err1);
	  float errr = hd->GetBinContent(ib+1)*err1;
	  std::cout << ib << " : " << hd->GetBinContent(ib+1) << "  " << errr << std::endl;
	  hd->SetBinError(ib+1, errr);
	}
      /* hd->Divide(h_reco_dphi_rangeminmin[0][i]); */
      /* for (int ib = 0; ib < hd->GetNbinsX(); ib++) */
      /* 	{ */
      /* 	  float cov = h_corr_reco_dphi_rangemin[0][i]->GetBinContent(ib+1); */
      /* 	  float err1 = 0; */
      /* 	  err1 += TMath::Power(h_reco_dphi_rangemin[0][i]->GetBinError(ib+1)/h_reco_dphi_rangemin[0][i]->GetBinContent(ib+1), 2); */
      /* 	  err1 += TMath::Power(h_truth_match_dphi_rangemin[0][i]->GetBinError(ib+1)/h_truth_match_dphi_rangemin[0][i]->GetBinContent(ib+1), 2); */
      /* 	  err1 -= 2.*cov/(h_reco_dphi_rangemin[0][i]->GetBinContent(ib+1)*h_truth_match_dphi_rangemin[0][i]->GetBinContent(ib+1)); */
      /* 	  err1 = sqrt(err1); */
      /* 	  float errr = hd->GetBinContent(ib+1)*err1; */

      /* 	  hd->SetBinError(ib+1, errr); */
      /* 	} */

      hd->SetTitle(";#Delta#phi; Truth MC / Reco MC");
      dlutility::SetFont(hd, 42, 0.09, 0.08, 0.075, 0.075);
      hd->GetYaxis()->SetTitleOffset(0.8);
      hd->SetMaximum(2);
      hd->SetMinimum(0);

      hd->Draw("p");
      TLine *li1 = new TLine(TMath::Pi()*3./4., 1, TMath::Pi(), 1);
      dlutility::SetLineAtt(li1, kBlack, 1, 4);
      li1->Draw("same");

      c_cor->Print(Form("correction_rangemin_%d.pdf", i));
      c_cor->Print(Form("correction_rangemin_%d.png", i));
    }
  
  TCanvas *c_sys = new TCanvas("c_sys","c_sys", 1000, 500);
  c_sys->Divide(3, 1);
  for (int j = 0; j < 3; j++)
    {
      c_sys->cd(j+1);
      for (int i = 0; i < 5; i++)
	{
	  dlutility::SetLineAtt(h_data_dphi_range[i][j], colorsys[i], 1, 1);
	  dlutility::SetMarkerAtt(h_data_dphi_range[i][j], colorsys[i], 0.5, 8);
	}
      dlutility::SetLineAtt(h_data_dphi_range[5][j], colorsys[1], 1, 1);
      dlutility::SetMarkerAtt(h_data_dphi_range[5][j], colorsys[1], 0.5, 8);
      dlutility::SetLineAtt(h_data_dphi_range[6][j], colorsys[2], 1, 1);
      dlutility::SetMarkerAtt(h_data_dphi_range[6][j], colorsys[2], 0.5, 8);

      std::cout << h_data_dphi_range[0][j]->GetMean() << std::endl;
      h_data_dphi_range[0][j]->Draw("p");
      h_data_dphi_range[5][j]->Draw("p same");
      h_data_dphi_range[6][j]->Draw("p same");

    }
  


  TCanvas *c_sys1 = new TCanvas("c_sys1","c_sys1", 500, 600);
  dlutility::ratioPanelCanvas(c_sys1, 0.3);

  for (int irange = 0; irange < 3; irange++)
    {
      c_sys1->cd(1);
      dlutility::SetLineAtt(h_data_dphi_range[0][irange], color_unfold, 1, 1);
      dlutility::SetMarkerAtt(h_data_dphi_range[0][irange], color_unfold, 1, 8);
      dlutility::SetLineAtt(h_data_dphi_range_sys[irange], color_unfold, 1, 1);
      dlutility::SetMarkerAtt(h_data_dphi_range_sys[irange], color_unfold, 1, 8);

      h_data_dphi_range_sys[irange]->SetFillColorAlpha(color_unfold_fill, 0.3);
      h_data_dphi_range[0][irange]->SetMaximum(10);

      dlutility::SetFont(h_data_dphi_range[0][irange], 42, 0.06);
      dlutility::SetLineAtt(h_truth_match_dphi_range[0][irange], kRed, 3, 1);
      h_data_dphi_range[0][irange]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}");

      dlutility::SetLineAtt(h_herwig_leading_dphi[irange], color_herwig, 3, 1);

      h_data_dphi_range[0][irange]->Draw("p E1");
      h_data_dphi_range[0][irange]->Draw("p E1 same");
      h_data_dphi_range_sys[irange]->Draw("E2 p same");
      float ss = 0.05;
      float topp = 0.8;
      dlutility::DrawSPHENIXppPrelimsize(0.22, topp, ss);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R} = 0.4}", 0.22, topp - 2*(ss+0.01), 0, kBlack,ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,1} < %2.1f GeV}", ipt_bins[binranges[irange]], ipt_bins[binranges[irange+1]]), 0.22, topp - 3*(ss+0.01),0, kBlack,ss);
      dlutility::drawText(Form("p_{T,2} #kern[-0.08]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, topp - 4*(ss+0.01), 0, kBlack,ss);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, topp - 5*(ss+0.01),0, kBlack,ss);
      TLegend *leg = new TLegend(0.22, topp - 9*(ss+0.01), 0.4, topp - 6*(ss+0.01));
      leg->SetLineWidth(0);
      leg->SetTextSize(0.05);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_range[0][irange], "Data","lp");
      leg->Draw("same");
  
      c_sys1->cd(2);
      TH1D *hratstatpyth = (TH1D*) h_truth_match_dphi_range[0][irange]->Clone();
      TH1D *hratstatherwig = (TH1D*) h_herwig_leading_dphi[irange]->Clone();
      dlutility::SetFont(hratstatpyth, 42, 0.16, 0.15, 0.14, 0.14);
      hratstatpyth->GetYaxis()->SetTitleOffset(0.5);
      hratstatpyth->GetXaxis()->SetTitleOffset(0.8);
      dlutility::SetLineAtt(hratstatpyth, color_pythia, 1, 1);
      dlutility::SetMarkerAtt(hratstatpyth, color_pythia, 1, 8);
      dlutility::SetLineAtt(hratstatherwig, color_herwig, 1, 1);
      dlutility::SetMarkerAtt(hratstatherwig, color_herwig, 1, 8);

      dlutility::SetLineAtt(  h_ratio_dphi_range_sys[irange], kBlack, 1, 1);
      dlutility::SetMarkerAtt(  h_ratio_dphi_range_sys[irange], kBlack, 1, 1);
      h_ratio_dphi_range_sys[irange]->SetFillColorAlpha(kBlack, 0.3);
  
      hratstatpyth->Divide(h_data_dphi_range[0][irange]);
      hratstatherwig->Divide(h_data_dphi_range[0][irange]);

      TH1D *hcopy = (TH1D*) hratstatpyth->Clone();
      hcopy->Reset();
      
      dlutility::SetFont(hcopy, 42, 0.16, 0.15, 0.14, 0.14);
      hcopy->SetMaximum(3);
      hcopy->SetMinimum(0);
      hcopy->SetTitle(";#Delta#phi; MC/Data");
      hcopy->Draw("p E2");
      h_ratio_dphi_range_sys[irange]->Draw("p E2 same");
      
      //hratstatpyth->Draw("p E1");
      //hratstatherwig->Draw("same p E1");

      TLine *li = new TLine(TMath::Pi()*3/4., 1, TMath::Pi(), 1);
      dlutility::SetLineAtt(li, kBlack, 1, 4);
      li->Draw("same");

      c_sys1->Print(Form("dphi_data_range_%d.pdf", irange));
      c_sys1->Print(Form("dphi_data_range_%d.png", irange));

      c_sys1->cd(1);

      h_data_dphi_range[0][irange]->Draw("p E1");
      h_truth_match_dphi_range[0][irange]->Draw("same hist C");
      h_data_dphi_range[0][irange]->Draw("p E1 same");
      h_data_dphi_range_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXppPrelimsize(0.22, topp, ss);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R} = 0.4}", 0.22, topp - 2*(ss+0.01), 0, kBlack,ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,1} < %2.1f GeV}", ipt_bins[binranges[irange]], ipt_bins[binranges[irange+1]]), 0.22, topp - 3*(ss+0.01),0, kBlack,ss);
      dlutility::drawText(Form("p_{T,2} #kern[-0.08]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, topp - 4*(ss+0.01), 0, kBlack,ss);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, topp - 5*(ss+0.01),0, kBlack,ss);
      leg = new TLegend(0.22, topp - 9*(ss+0.01), 0.4, topp - 6*(ss+0.01));
      leg->SetLineWidth(0);
      leg->SetTextSize(0.05);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_range[0][irange], "Data","lp");
      leg->AddEntry(h_truth_match_dphi_range[0][irange], "PYTHIA-8","l");
      leg->Draw("same");
  
      c_sys1->cd(2);
      hratstatpyth->SetMaximum(3.0);
      hratstatpyth->SetMinimum(0);
      hratstatpyth->SetTitle(";#Delta#phi; MC/Data");

      hratstatpyth->Draw("p E1");
      h_ratio_dphi_range_sys[irange]->Draw("p E2 same");
      li->Draw("same");

      c_sys1->Print(Form("dphi_data_pythia_range_%d.pdf", irange));
      c_sys1->Print(Form("dphi_data_pythia_range_%d.png", irange));

      c_sys1->cd(1);
      h_data_dphi_range[0][irange]->Draw("p E1");
      h_truth_match_dphi_range[0][irange]->Draw("same hist C");
      h_herwig_leading_dphi[irange]->Draw("same hist C");
      h_data_dphi_range[0][irange]->Draw("p E1 same");
      h_data_dphi_range_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXppPrelimsize(0.22, topp, ss);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R} = 0.4}", 0.22, topp - 2*(ss+0.01), 0, kBlack,ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,1} < %2.1f GeV}", ipt_bins[binranges[irange]], ipt_bins[binranges[irange+1]]), 0.22, topp - 3*(ss+0.01),0, kBlack,ss);
      dlutility::drawText(Form("p_{T,2} #kern[-0.08]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, topp - 4*(ss+0.01), 0, kBlack,ss);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, topp - 5*(ss+0.01),0, kBlack,ss);
      leg = new TLegend(0.22, topp - 9*(ss+0.01), 0.4, topp - 6*(ss+0.01));
      leg->SetLineWidth(0);
      leg->SetTextSize(0.05);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_range[0][irange], "Data","lp");
      leg->AddEntry(h_truth_match_dphi_range[0][irange], "PYTHIA-8","l");
      leg->AddEntry(h_herwig_leading_dphi[irange], "Herwig 7.3","l");
      leg->Draw("same");
  
      c_sys1->cd(2);
      hratstatpyth->Draw("p E1");
      hratstatherwig->Draw("same p E1");
      h_ratio_dphi_range_sys[irange]->Draw("p E2 same");
      
      li->Draw("same");
      c_sys1->Print(Form("dphi_range_%d.pdf", irange));
      c_sys1->Print(Form("dphi_range_%d.png", irange));

    }

  for (int irange = 0; irange < 3; irange++)
    {
      c_sys1->cd(1);

      dlutility::SetLineAtt(h_data_dphi_rangemin[0][irange], color_unfold, 1, 1);
      dlutility::SetMarkerAtt(h_data_dphi_rangemin[0][irange], color_unfold, 1, 8);
      dlutility::SetLineAtt(h_data_dphi_rangemin_sys[irange], color_unfold, 1, 1);
      dlutility::SetMarkerAtt(h_data_dphi_rangemin_sys[irange], color_unfold, 1, 8);
      h_data_dphi_rangemin_sys[irange]->SetFillColorAlpha(color_unfold_fill, 0.3);

      /* h_truth_match_dphi_rangemin[0][irange]->Scale(1./ h_truth_match_dphi_rangemin[0][irange]->Integral(), "width"); */
      /* h_data_dphi_rangemin[0][irange]->Scale(1./ h_data_dphi_rangemin[0][irange]->Integral(), "width"); */
      /* h_data_dphi_rangemin_sys[irange]->Scale(1./ h_data_dphi_rangemin_sys[irange]->Integral(), "width"); */

      h_data_dphi_rangemin[0][irange]->SetMaximum(10);
      dlutility::SetLineAtt(h_herwig_subleading_dphi[irange], color_herwig, 3, 1);
      dlutility::SetFont(h_data_dphi_rangemin[0][irange], 42, 0.06);
      dlutility::SetLineAtt(h_truth_match_dphi_rangemin[0][irange], kRed, 3, 1);
      h_data_dphi_rangemin[0][irange]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}");
      h_data_dphi_rangemin[0][irange]->Draw("p E1");
      h_data_dphi_rangemin[0][irange]->Draw("p E1 same");
      h_data_dphi_rangemin_sys[irange]->Draw("E2 p same");
      //hcopy
      float ss = 0.05;
      float topp = 0.8;
      dlutility::DrawSPHENIXppPrelimsize(0.22, topp, ss);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R} = 0.4}", 0.22, topp - 2*(ss+0.01), 0, kBlack,ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,1} < %2.1f GeV}", ipt_bins[binranges[1]], ipt_bins[binranges[1+1]]), 0.22, topp - 3*(ss+0.01),0, kBlack,ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,2} < %2.1f GeV} ", ipt_bins[binrangesmin[irange]], ipt_bins[binrangesmin[irange+1]]), 0.22, topp - 4*(ss+0.01), 0, kBlack, ss);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, topp - 5*(ss+0.01),0, kBlack,ss);
      TLegend *leg = new TLegend(0.22, topp - 9*(ss+0.01), 0.4, topp - 6*(ss+0.01));
      leg->SetLineWidth(0);
      leg->SetTextSize(0.05);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_rangemin[0][irange], "Data","lp");
      leg->AddEntry(h_truth_match_dphi_rangemin[0][irange], "PYTHIA-8","l");
      leg->AddEntry(h_herwig_subleading_dphi[irange], "Herwig 7.3","l");
      leg->Draw("same");
      c_sys1->cd(2);
      TH1D *hratstatpyth = (TH1D*) h_truth_match_dphi_rangemin[0][irange]->Clone();
      TH1D *hratstatherwig = (TH1D*) h_herwig_subleading_dphi[irange]->Clone();
      TH1D *hcopy = (TH1D*) hratstatpyth->Clone();
      hcopy->Reset();      
      dlutility::SetFont(hcopy, 42, 0.16, 0.15, 0.14, 0.14);
      hcopy->SetMaximum(3);
      hcopy->SetMinimum(0);
      hcopy->SetTitle(";#Delta#phi; MC/Data");
      hcopy->Draw("p E2");
      h_ratio_dphi_rangemin_sys[irange]->Draw("p E2 same");

      
      dlutility::SetFont(hratstatpyth, 42, 0.16, 0.15, 0.14, 0.14);
      hratstatpyth->GetYaxis()->SetTitleOffset(0.5);
      hratstatpyth->GetXaxis()->SetTitleOffset(0.8);
      dlutility::SetLineAtt(hratstatpyth, color_pythia, 1, 1);
      dlutility::SetMarkerAtt(hratstatpyth, color_pythia, 1, 8);
      dlutility::SetLineAtt(hratstatherwig, color_herwig, 1, 1);
      dlutility::SetMarkerAtt(hratstatherwig, color_herwig, 1, 8);

      dlutility::SetLineAtt(  h_ratio_dphi_rangemin_sys[irange], kBlack, 1, 1);
      dlutility::SetMarkerAtt(  h_ratio_dphi_rangemin_sys[irange], kBlack, 1, 1);
      h_ratio_dphi_rangemin_sys[irange]->SetFillColorAlpha(kBlack, 0.3);
  
      hratstatpyth->Divide(h_data_dphi_rangemin[0][irange]);
      hratstatherwig->Divide(h_data_dphi_rangemin[0][irange]);

      
      hratstatpyth->SetMaximum(3.0);
      hratstatpyth->SetMinimum(0);

      hratstatpyth->SetTitle(";#Delta#phi; MC/Data");
 
      TLine *li = new TLine(TMath::Pi()*3/4., 1, TMath::Pi(), 1);
      dlutility::SetLineAtt(li, kBlack, 1, 4);
      li->Draw("same");

      c_sys1->Print(Form("dphi_data_rangemin_%d.pdf", irange));
      c_sys1->Print(Form("dphi_data_rangemin_%d.png", irange));

      c_sys1->cd(1);

      
      h_data_dphi_rangemin[0][irange]->Draw("p E1");
      h_truth_match_dphi_rangemin[0][irange]->Draw("same hist C");
      h_data_dphi_rangemin[0][irange]->Draw("p E1 same");
      h_data_dphi_rangemin_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXppPrelimsize(0.22, topp, ss);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R} = 0.4}", 0.22, topp - 2*(ss+0.01), 0, kBlack,ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,1} < %2.1f GeV}", ipt_bins[binranges[1]], ipt_bins[binranges[1+1]]), 0.22, topp - 3*(ss+0.01),0, kBlack,ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,2} < %2.1f GeV} ", ipt_bins[binrangesmin[irange]], ipt_bins[binrangesmin[irange+1]]), 0.22, topp - 4*(ss+0.01), 0, kBlack, ss);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, topp - 5*(ss+0.01),0, kBlack,ss);
      leg = new TLegend(0.22, topp - 9*(ss+0.01), 0.4, topp - 6*(ss+0.01));
      leg->SetLineWidth(0);
      leg->SetTextSize(0.05);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_rangemin[0][irange], "Data","lp");
      leg->AddEntry(h_truth_match_dphi_rangemin[0][irange], "PYTHIA-8","l");
      leg->Draw("same");
      c_sys1->cd(2);

      hratstatpyth->SetTitle(";#Delta#phi; MC/Data");
      hratstatpyth->Draw("p E1");
      h_ratio_dphi_rangemin_sys[irange]->Draw("p E2 same");
      hratstatpyth->Draw("p E1 same");

      li->Draw("same");
      c_sys1->Print(Form("dphi_data_pythia_rangemin_%d.pdf", irange));
      c_sys1->Print(Form("dphi_data_pythia_rangemin_%d.png", irange));

      c_sys1->cd(1);

      h_data_dphi_rangemin[0][irange]->Draw("p E1");
      h_truth_match_dphi_rangemin[0][irange]->Draw("same hist C");
      h_herwig_subleading_dphi[irange]->Draw("same hist C");
      h_data_dphi_rangemin[0][irange]->Draw("p E1 same");
      h_data_dphi_rangemin_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXppPrelimsize(0.22, topp, ss);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R} = 0.4}", 0.22, topp - 2*(ss+0.01), 0, kBlack,ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,1} < %2.1f GeV}", ipt_bins[binranges[1]], ipt_bins[binranges[1+1]]), 0.22, topp - 3*(ss+0.01),0, kBlack,ss);
      dlutility::drawText(Form("%2.1f #kern[-0.08]{#leq p_{T,2} < %2.1f GeV} ", ipt_bins[binrangesmin[irange]], ipt_bins[binrangesmin[irange+1]]), 0.22, topp - 4*(ss+0.01), 0, kBlack, ss);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, topp - 5*(ss+0.01),0, kBlack,ss);
      leg = new TLegend(0.22, topp - 9*(ss+0.01), 0.4, topp - 6*(ss+0.01));
      leg->SetLineWidth(0);
      leg->SetTextSize(0.05);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_rangemin[0][irange], "Data","lp");
      leg->AddEntry(h_truth_match_dphi_rangemin[0][irange], "PYTHIA-8","l");
      leg->AddEntry(h_herwig_subleading_dphi[irange], "Herwig 7.3","l");
      leg->Draw("same");
      c_sys1->cd(2);

      hratstatpyth->Draw("p E1");
      hratstatherwig->Draw("same p E1");
      h_ratio_dphi_rangemin_sys[irange]->Draw("p E2 same");
      hratstatpyth->Draw("p E1 same");
      hratstatherwig->Draw("same p E1 same ");

      li->Draw("same");
      c_sys1->Print(Form("dphi_rangemin_%d.pdf", irange));
      c_sys1->Print(Form("dphi_rangemin_%d.png", irange));

    }

  TCanvas *c_money = new TCanvas("c_money", "c_money", 500, 500);
  for (int irange = 0; irange < 3; irange++)
    {
      dlutility::SetLineAtt(h_data_dphi_rangemin[0][irange], color_unfold, 1, 1);
      dlutility::SetMarkerAtt(h_data_dphi_rangemin[0][irange], color_unfold, 1, 8);
      dlutility::SetLineAtt(h_data_dphi_rangemin_sys[irange], color_unfold, 1, 1);
      dlutility::SetMarkerAtt(h_data_dphi_rangemin_sys[irange], color_unfold, 1, 8);
      h_data_dphi_rangemin_sys[irange]->SetFillColorAlpha(color_unfold_fill, 0.3);

      /* h_truth_match_dphi_rangemin[0][irange]->Scale(1./ h_truth_match_dphi_rangemin[0][irange]->Integral(), "width"); */
      /* h_data_dphi_rangemin[0][irange]->Scale(1./ h_data_dphi_rangemin[0][irange]->Integral(), "width"); */
      /* h_data_dphi_rangemin_sys[irange]->Scale(1./ h_data_dphi_rangemin_sys[irange]->Integral(), "width"); */

      h_data_dphi_rangemin[0][irange]->SetMaximum(10);
      dlutility::SetLineAtt(h_herwig_subleading_dphi[irange], color_herwig, 3, 1);
      dlutility::SetFont(h_data_dphi_rangemin[0][irange], 42, 0.06, 0.04, 0.05, 0.05 );
      dlutility::SetLineAtt(h_truth_match_dphi_rangemin[0][irange], kRed, 3, 1);
      h_data_dphi_rangemin[0][irange]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}");
      h_data_dphi_rangemin[0][irange]->Draw("p E1");
      h_truth_match_dphi_rangemin[0][irange]->Draw("same hist C");
      h_herwig_subleading_dphi[irange]->Draw("same hist C");
      h_data_dphi_rangemin[0][irange]->Draw("p E1 same");
      h_data_dphi_rangemin_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXppPrelim(0.22, 0.84);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R = 0.4}}", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq p_{T,1} < %2.1f GeV} ", ipt_bins[binranges[1]], ipt_bins[binranges[1+1]]), 0.22, 0.69);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq p_{T,2} < %2.1f GeV}", ipt_bins[binrangesmin[irange]], ipt_bins[binrangesmin[irange+1]]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, 0.59);

      TLegend *leg = new TLegend(0.22, 0.39, 0.4, 0.54);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_rangemin[0][irange], "Data","lp");
      leg->AddEntry(h_truth_match_dphi_rangemin[0][irange], "PYTHIA-8", "l");
      leg->AddEntry(h_herwig_subleading_dphi[irange], "Herwig 7.3","l");
      leg->Draw("same");

      c_money->Print(Form("h_final_dphi_rangemin_%d.png", irange));
      c_money->Print(Form("h_final_dphi_rangemin_%d.pdf", irange));

      h_data_dphi_rangemin[0][irange]->Draw("p E1");
      h_data_dphi_rangemin_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXppPrelim(0.22, 0.84);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R = 0.4}}", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq p_{T,1} < %2.1f GeV} ", ipt_bins[binranges[1]], ipt_bins[binranges[1+1]]), 0.22, 0.69);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq p_{T,2} < %2.1f GeV}", ipt_bins[binrangesmin[irange]], ipt_bins[binrangesmin[irange+1]]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, 0.59);

      leg = new TLegend(0.22, 0.39, 0.4, 0.54);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_rangemin[0][irange], "Data","lp");
      leg->Draw("same");

      c_money->Print(Form("h_final_data_dphi_rangemin_%d.png", irange));
      c_money->Print(Form("h_final_data_dphi_rangemin_%d.pdf", irange));

      h_data_dphi_rangemin[0][irange]->Draw("p E1");
      h_truth_match_dphi_rangemin[0][irange]->Draw("same hist C");
      h_data_dphi_rangemin[0][irange]->Draw("p E1 same");
      h_data_dphi_rangemin_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXppPrelim(0.22, 0.84);
      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R = 0.4}}", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq p_{T,1} < %2.1f GeV} ", ipt_bins[binranges[1]], ipt_bins[binranges[1+1]]), 0.22, 0.69);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq p_{T,2} < %2.1f GeV}", ipt_bins[binrangesmin[irange]], ipt_bins[binrangesmin[irange+1]]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, 0.59);

      leg = new TLegend(0.22, 0.39, 0.4, 0.54);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_rangemin[0][irange], "Data","lp");
      leg->AddEntry(h_truth_match_dphi_rangemin[0][irange], "PYTHIA-8", "l");
      leg->Draw("same");

      c_money->Print(Form("h_final_data_pythia_dphi_rangemin_%d.png", irange));
      c_money->Print(Form("h_final_data_pythia_dphi_rangemin_%d.pdf", irange));


      
    }
  for (int irange = 0; irange < 3; irange++)
    {
      dlutility::SetLineAtt(h_data_dphi_range[0][irange], color_unfold, 1, 1);
      dlutility::SetMarkerAtt(h_data_dphi_range[0][irange], color_unfold, 1, 8);
      dlutility::SetLineAtt(h_data_dphi_range_sys[irange], color_unfold, 1, 1);
      dlutility::SetMarkerAtt(h_data_dphi_range_sys[irange], color_unfold, 1, 8);
      h_data_dphi_range_sys[irange]->SetFillColorAlpha(color_unfold_fill, 0.3);

      /* h_truth_match_dphi_range[0][irange]->Scale(1./ h_truth_match_dphi_range[0][irange]->Integral(), "width"); */
      /* h_data_dphi_range[0][irange]->Scale(1./ h_data_dphi_range[0][irange]->Integral(), "width"); */
      /* h_data_dphi_range_sys[irange]->Scale(1./ h_data_dphi_range_sys[irange]->Integral(), "width"); */

      h_data_dphi_range[0][irange]->SetMaximum(10);
      dlutility::SetLineAtt(h_herwig_subleading_dphi[irange], color_herwig, 3, 1);
      dlutility::SetFont(h_data_dphi_range[0][irange], 42, 0.06, 0.04, 0.05, 0.05 );
      dlutility::SetLineAtt(h_truth_match_dphi_range[0][irange], kRed, 3, 1);
      h_data_dphi_range[0][irange]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}");
      h_data_dphi_range[0][irange]->Draw("p E1");
      h_truth_match_dphi_range[0][irange]->Draw("same hist C");
      h_herwig_leading_dphi[irange]->Draw("same hist C");
      h_data_dphi_range[0][irange]->Draw("p E1 same");
      h_data_dphi_range_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXppPrelim(0.22, 0.84);

      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R = 0.4}}", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq p_{T,1} < %2.1f GeV} ", ipt_bins[binranges[irange]], ipt_bins[binranges[irange+1]]), 0.22, 0.69);
      dlutility::drawText(Form("p_{T,2} #kern[-0.07]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, 0.59);


      TLegend *leg = new TLegend(0.22, 0.39, 0.4, 0.54);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_range[0][irange], "Data","lp");
      leg->AddEntry(h_truth_match_dphi_range[0][irange], "PYTHIA-8", "l");
      leg->AddEntry(h_herwig_subleading_dphi[irange], "Herwig 7.3","l");
      leg->Draw("same");
      c_money->Print(Form("h_final_dphi_range_%d.png", irange));
      c_money->Print(Form("h_final_dphi_range_%d.pdf", irange));
      h_data_dphi_range[0][irange]->SetMaximum(10);
      dlutility::SetLineAtt(h_herwig_subleading_dphi[irange], color_herwig, 3, 1);
      dlutility::SetFont(h_data_dphi_range[0][irange], 42, 0.06, 0.04, 0.05, 0.05 );
      dlutility::SetLineAtt(h_truth_match_dphi_range[0][irange], kRed, 3, 1);


      h_data_dphi_range[0][irange]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}");
      h_data_dphi_range[0][irange]->Draw("p E1");
      h_data_dphi_range[0][irange]->Draw("p E1 same");
      h_data_dphi_range_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXppPrelim(0.22, 0.84);

      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R = 0.4}}", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq p_{T,1} < %2.1f GeV} ", ipt_bins[binranges[irange]], ipt_bins[binranges[irange+1]]), 0.22, 0.69);
      dlutility::drawText(Form("p_{T,2} #kern[-0.07]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, 0.59);


      leg = new TLegend(0.22, 0.39, 0.4, 0.54);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_range[0][irange], "Data","lp");
      leg->Draw("same");
      c_money->Print(Form("h_final_data_dphi_range_%d.png", irange));
      c_money->Print(Form("h_final_data_dphi_range_%d.pdf", irange));

      h_data_dphi_range[0][irange]->SetMaximum(10);
      dlutility::SetLineAtt(h_herwig_subleading_dphi[irange], color_herwig, 3, 1);
      dlutility::SetFont(h_data_dphi_range[0][irange], 42, 0.06, 0.04, 0.05, 0.05 );
      dlutility::SetLineAtt(h_truth_match_dphi_range[0][irange], kRed, 3, 1);
      h_data_dphi_range[0][irange]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{dN_{pair}}{d#Delta#phi}");
      h_data_dphi_range[0][irange]->Draw("p E1");
      h_truth_match_dphi_range[0][irange]->Draw("same hist C");
      h_data_dphi_range[0][irange]->Draw("p E1 same");
      h_data_dphi_range_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXppPrelim(0.22, 0.84);

      dlutility::drawText("anti-#it{k_{t}} #kern[-0.1]{#it{R = 0.4}}", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f #kern[-0.07]{#leq p_{T,1} < %2.1f GeV} ", ipt_bins[binranges[irange]], ipt_bins[binranges[irange+1]]), 0.22, 0.69);
      dlutility::drawText(Form("p_{T,2} #kern[-0.07]{#geq %2.1f GeV}", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #kern[-0.15]{#geq 3#pi/4}", 0.22, 0.59);


      leg = new TLegend(0.22, 0.39, 0.4, 0.54);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_range[0][irange], "Data","lp");
      leg->AddEntry(h_truth_match_dphi_range[0][irange], "PYTHIA-8", "l");
      leg->Draw("same");
      c_money->Print(Form("h_final_data_pythia_dphi_range_%d.png", irange));
      c_money->Print(Form("h_final_data_pythia_dphi_range_%d.pdf", irange));

    }

  
  TCanvas *cc = new TCanvas("cc","cc",500, 500);
  dlutility::SetLineAtt(h_correction_factors_range[0][0], kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_correction_factors_range[0][0], kBlack, 1, 8);
  dlutility::SetLineAtt(h_correction_factors_range[0][1], kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_correction_factors_range[0][1], kBlue, 1, 8);
  dlutility::SetLineAtt(h_correction_factors_range[0][2], kRed, 1, 1);
  dlutility::SetMarkerAtt(h_correction_factors_range[0][2], kRed, 1, 8);


  h_correction_factors_range[0][0]->SetMaximum(100);
  h_correction_factors_range[0][0]->SetMinimum(00);
  h_correction_factors_range[0][0]->Draw("p");
  h_correction_factors_range[0][1]->Draw("same p");
  h_correction_factors_range[0][2]->Draw("same p");

  TCanvas *cc2 = new TCanvas("cc2","cc2",500, 500);
  dlutility::SetLineAtt(h_correction_factors_rangemin[0][0], kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_correction_factors_rangemin[0][0], kBlack, 1, 8);
  dlutility::SetLineAtt(h_correction_factors_rangemin[0][1], kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_correction_factors_rangemin[0][1], kBlue, 1, 8);
  dlutility::SetLineAtt(h_correction_factors_rangemin[0][2], kRed, 1, 1);
  dlutility::SetMarkerAtt(h_correction_factors_rangemin[0][2], kRed, 1, 8);


  h_correction_factors_rangemin[0][0]->SetMaximum(100);
  h_correction_factors_rangemin[0][0]->SetMinimum(00);

  h_correction_factors_rangemin[0][0]->Draw("p");
  h_correction_factors_rangemin[0][1]->Draw("same p");
  h_correction_factors_rangemin[0][2]->Draw("same p");

  
  TCanvas *ce = new TCanvas("ce","ce",500, 500);
  dlutility::SetLineAtt(h_efficiency_factors_range[0][0], kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_efficiency_factors_range[0][0], kBlack, 1, 8);
  dlutility::SetLineAtt(h_efficiency_factors_range[0][1], kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_efficiency_factors_range[0][1], kBlue, 1, 8);
  dlutility::SetLineAtt(h_efficiency_factors_range[0][2], kRed, 1, 1);
  dlutility::SetMarkerAtt(h_efficiency_factors_range[0][2], kRed, 1, 8);


  h_efficiency_factors_range[0][0]->SetMaximum(100);
  h_efficiency_factors_range[0][0]->SetMinimum(00);
  h_efficiency_factors_range[0][0]->Draw("p");
  h_efficiency_factors_range[0][1]->Draw("same p");
  h_efficiency_factors_range[0][2]->Draw("same p");

  TCanvas *ce2 = new TCanvas("ce2","ce2",500, 500);
  dlutility::SetLineAtt(h_efficiency_factors_rangemin[0][0], kBlack, 1, 1);
  dlutility::SetMarkerAtt(h_efficiency_factors_rangemin[0][0], kBlack, 1, 8);
  dlutility::SetLineAtt(h_efficiency_factors_rangemin[0][1], kBlue, 1, 1);
  dlutility::SetMarkerAtt(h_efficiency_factors_rangemin[0][1], kBlue, 1, 8);
  dlutility::SetLineAtt(h_efficiency_factors_rangemin[0][2], kRed, 1, 1);
  dlutility::SetMarkerAtt(h_efficiency_factors_rangemin[0][2], kRed, 1, 8);


  h_efficiency_factors_rangemin[0][0]->SetMaximum(100);
  h_efficiency_factors_rangemin[0][0]->SetMinimum(00);

  h_efficiency_factors_rangemin[0][0]->Draw("p");
  h_efficiency_factors_rangemin[0][1]->Draw("same p");
  h_efficiency_factors_rangemin[0][2]->Draw("same p");


}

std::pair<float, float> get_ordered_pair(float e1, float e2)
{
  float e1f = fabs(e1);
  float e2f = fabs(e2);
  if (e1 > 0 && e2 < 0)
    {
      return std::make_pair(e1f, e2f);      
    }
  else if (e2 > 0 && e1 < 0)
    {
      return std::make_pair(e2f, e1f);      
    }
  else if (e1 > 0 && e2 > 0)
    {
      return std::make_pair(max(e1f, e2f), 0);      
    }
  else if (e1 < 0 && e2 < 0)
    {
      return std::make_pair(0, max(e1f, e2f));      
    }
  else
    return std::make_pair(0, 0);      

}
