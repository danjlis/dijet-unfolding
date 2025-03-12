#include "../macros/dlUtility.h"
#include "read_binning.h"

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


  
  TString filenames[5];
  filenames[0] = "dphi_hists.root";
  filenames[1] = "dphi_hists_negJER.root";
  filenames[2] = "dphi_hists_posJER.root";
  filenames[3] = "dphi_hists_negJES.root";
  filenames[4] = "dphi_hists_posJES.root";

  TFile *fin[5];
  
  TH3D *h_data_pt1pt2dphi[5];//new TH3D("h_truth_pt1pt2",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  
  TH3D *h_truth_pt1pt2dphi[5];//new TH3D("h_reco_match_pt1pt2dphi",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_truth_match_pt1pt2dphi[5];//new TH3D("h_reco_match_pt1pt2dphi",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphi[5];//new TH3D("h_reco_match_pt1pt2dphi",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphitruth[5];//new TH3D("h_reco_match_pt1pt2dphitruth",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  
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
      h_truth_match_pt1pt2dphi[isys] = (TH3D*) fin[isys]->Get("h_truth_pt1pt2");;
      h_truth_pt1pt2dphi[isys] = (TH3D*) fin[isys]->Get("h_truth_match_pt1pt2");;
      h_reco_match_pt1pt2dphitruth[isys] = (TH3D*) fin[isys]->Get("h_reco_match_pt1pt2dphitruth");;

      
      h_data_pt1pt2dphi[isys]->SetName(Form("h_data_pt1pt2dphi_%d", isys));
      
      h_reco_match_pt1pt2dphi[isys]->SetName(Form("h_reco_pt1pt2dphi_%d", isys));
      
      h_truth_match_pt1pt2dphi[isys]->SetName(Form("h_truth_match_pt1pt2dphi_%d", isys));
      
      h_truth_pt1pt2dphi[isys]->SetName(Form("h_truth_pt1pt2dphi_%d", isys));
      
      h_reco_match_pt1pt2dphitruth[isys]->SetName(Form("h_reco_pt1pt2dphitruth_%d", isys));

            
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

	      h_data_dphi_range[isys][i]->SetBinContent(ib+1, h_data_dphi_range[isys][i]->GetBinContent(ib+1)*h_correction_factors_range[isys][i]->GetBinContent(ib+1));//*h_efficiency_factors_range[isys][i]->GetBinContent(ib+1));
	      if (isys == 0)
		{
		  float half_fac = ( h_correction_factors_range[isys][i]->GetBinContent(ib+1) - 1)*0.5 + 1;
		  float double_fac = ( h_correction_factors_range[isys][i]->GetBinContent(ib+1) - 1)*1.5 + 1;
		  h_data_dphi_range[5][i]->SetBinContent(ib+1, h_data_dphi_range[5][i]->GetBinContent(ib+1)*half_fac);//*h_efficiency_factors_range[isys][i]->GetBinContent(ib+1));
		  h_data_dphi_range[6][i]->SetBinContent(ib+1, h_data_dphi_range[6][i]->GetBinContent(ib+1)*double_fac);//*h_efficiency_factors_range[isys][i]->GetBinContent(ib+1));
		}
	    }
	        
	  h_data_dphi_range[isys][i]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
	  h_truth_dphi_range[isys][i]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
	  h_truth_match_dphi_range[isys][i]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
	  h_reco_dphi_range[isys][i]->SetDirectory(0);//(Form("h_reco_match_dphi_%d", isys));
	  h_data_dphi_range[isys][i]->Rebin(2);
	  h_truth_dphi_range[isys][i]->Rebin(2);
	  h_truth_match_dphi_range[isys][i]->Rebin(2);
	  h_reco_dphi_range[isys][i]->Rebin(2);
	  if (isys == 0)
	    {
	      h_data_dphi_range[5][i]->Rebin(2);
	      h_data_dphi_range[5][i]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
	      h_data_dphi_range[6][i]->Rebin(2);
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

		      h_data_dphi_rangemin[isys][j]->SetBinContent(ib+1, h_data_dphi_rangemin[isys][j]->GetBinContent(ib+1)*h_correction_factors_rangemin[isys][j]->GetBinContent(ib+1));//*h_efficiency_factors_rangemin[isys][j]->GetBinContent(ib+1));
		      if (isys == 0)
			{
			  float half_fac = ( h_correction_factors_rangemin[isys][j]->GetBinContent(ib+1) - 1)*0.5 + 1;
			  float double_fac = ( h_correction_factors_rangemin[isys][j]->GetBinContent(ib+1) - 1)*1.5 + 1;

			  h_data_dphi_rangemin[5][j]->SetBinContent(ib+1, h_data_dphi_rangemin[5][j]->GetBinContent(ib+1)*half_fac);//*h_efficiency_factors_rangemin[isys][j]->GetBinContent(ib+1));
			  h_data_dphi_rangemin[6][j]->SetBinContent(ib+1, h_data_dphi_rangemin[6][j]->GetBinContent(ib+1)*double_fac);//*h_efficiency_factors_rangemin[isys][j]->GetBinContent(ib+1));
			}
		    }


		  h_data_dphi_rangemin[isys][j]->Rebin(2);//(1./h_data_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  h_truth_dphi_rangemin[isys][j]->Rebin(2);//(1./h_truth_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  h_truth_match_dphi_rangemin[isys][j]->Rebin(2);//(1./h_truth_match_dphi_rangemin[isys][j]->Integral(0, -1, "width"));
		  h_reco_dphi_rangemin[isys][j]->Rebin(2);//(1./h_reco_dphi_rangemin[isys][j]->Integral(0, -1, "width"));

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
		      h_data_dphi_rangemin[5][j]->Rebin(2);
		      h_data_dphi_rangemin[5][j]->Scale(1./h_data_dphi_rangemin[5][j]->Integral(0, -1, "width"));
		      h_data_dphi_rangemin[5][j]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
		      h_data_dphi_rangemin[6][j]->Rebin(2);
		      h_data_dphi_rangemin[6][j]->Scale(1./h_data_dphi_rangemin[6][j]->Integral(0, -1, "width"));
		      h_data_dphi_rangemin[6][j]->SetDirectory(0);//(Form("h_truth_match_dphi_%d", isys));
		    }
		}
	    }
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
	  h_ratio_dphi_range_sys[i]->SetPointY(ib, h_data_dphi_range[0][i]->GetBinContent(ib+1)/h_truth_match_dphi_range[0][i]->GetBinContent(ib+1));
	  h_ratio_dphi_range_sys[i]->SetPointError(ib, h_data_dphi_range[0][i]->GetBinWidth(ib+1)/2., h_data_dphi_range[0][i]->GetBinWidth(ib+1)/2., fabs(sys_neg)/h_truth_match_dphi_range[0][i]->GetBinContent(ib+1), fabs(sys_pos)/h_truth_match_dphi_range[0][i]->GetBinContent(ib+1));

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
	  h_ratio_dphi_rangemin_sys[i]->SetPointY(ib, h_data_dphi_rangemin[0][i]->GetBinContent(ib+1)/h_truth_match_dphi_rangemin[0][i]->GetBinContent(ib+1));
	  h_ratio_dphi_rangemin_sys[i]->SetPointError(ib, h_data_dphi_rangemin[0][i]->GetBinWidth(ib+1)/2., h_data_dphi_rangemin[0][i]->GetBinWidth(ib+1)/2., fabs(sys_neg)/h_truth_match_dphi_rangemin[0][i]->GetBinContent(ib+1), fabs(sys_pos)/h_truth_match_dphi_rangemin[0][i]->GetBinContent(ib+1));

	}
    }
  std::cout << __LINE__ << std::endl;
  TCanvas *c_sysdiv = new TCanvas("c_sysdiv", "c_sysdiv", 700, 500);
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
      h_all_sys_range[0][i]->SetMaximum(1.0);;
      h_all_sys_range[0][i]->SetMinimum(-0.6);;
      h_all_sys_range[0][i]->Draw("p");
      h_all_sys_range[1][i]->Draw("p same");

      h_corr_sys_range[0][i]->Draw("hist same");
      h_corr_sys_range[1][i]->Draw("hist same");
      h_jer_sys_range[0][i]->Draw("hist same");
      h_jer_sys_range[1][i]->Draw("hist same");
      h_jes_sys_range[0][i]->Draw("hist same");
      h_jes_sys_range[1][i]->Draw("hist same");
      dlutility::DrawSPHENIXpp(0.22, 0.84);
      dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[binranges[i]], ipt_bins[binranges[i+1]]), 0.22, 0.69);
      dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
      TLegend *leg = new TLegend(0.61, 0.72, 0.79, 0.88);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_all_sys_range[0][i], "All Systematics");
      leg->AddEntry(h_jer_sys_range[0][i], "JER");
      leg->AddEntry(h_jes_sys_range[0][i], "JES");
      leg->AddEntry(h_corr_sys_range[0][i], "Reco-correction");
	    

      leg->Draw("same");

      c_sysdiv->Print(Form("h_sys_all_range_%d.pdf", i));
      c_sysdiv->Print(Form("h_sys_all_range_%d.png", i));
    }
  TCanvas *c_sysdiv2 = new TCanvas("c_sysdiv2", "c_sysdiv2", 700, 500);
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
      h_all_sys_rangemin[0][i]->SetMaximum(1.0);;
      h_all_sys_rangemin[0][i]->SetMinimum(-0.6);;
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

      dlutility::DrawSPHENIXpp(0.22, 0.84);
      dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[binranges[1]], ipt_bins[binranges[1+1]]), 0.22, 0.69);
      dlutility::drawText(Form("%2.1f GeV #leq p_{T,2} < %2.1f GeV ", ipt_bins[binrangesmin[i]], ipt_bins[binrangesmin[i+1]]), 0.22, 0.64);
      //      dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
      TLegend *leg = new TLegend(0.61, 0.72, 0.79, 0.88);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_all_sys_rangemin[0][i], "All Systematics");
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

      dlutility::SetFont(h_truth_match_dphi_range[0][i], 42, 0.05);
      //h_truth_dphi_range[0][i]->Scale(1./ h_truth_dphi_range[0][i]->Integral(), "width");
      //h_reco_dphi_range[0][i]->Scale(1./ h_reco_dphi_range[0][i]->Integral(), "width");
      h_truth_match_dphi_range[0][i]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{d#Delta#phi}{dN_{pair}}");
      h_truth_match_dphi_range[0][i]->SetMaximum(10);
      h_truth_match_dphi_range[0][i]->SetMinimum(0);
      h_truth_match_dphi_range[0][i]->Draw("p");
      h_reco_dphi_range[0][i]->Draw("p same");

      dlutility::DrawSPHENIXpp(0.22, 0.84);
      dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[binranges[i]], ipt_bins[binranges[i+1]]), 0.22, 0.69);
      dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);

      TLegend *leg1 = new TLegend(0.20, 0.36, 0.63, 0.56);
      leg1->SetLineWidth(0);
      leg1->SetTextSize(0.04);
      leg1->SetTextFont(42);
      leg1->AddEntry(h_truth_match_dphi_range[0][i], "PYTHIA-8 truth","p");
      leg1->AddEntry(h_reco_dphi_range[0][i], "PYTHIA-8 reco","p");
      leg1->Draw("same");  

      c_cor->cd(2);
      gPad->SetTopMargin(0.05);

      TH1D *hd = (TH1D*) h_truth_match_dphi_range[0][i]->Clone();
      hd->Divide(h_reco_dphi_range[0][i]);
      hd->SetTitle(";#Delta#phi; Truth MC / Reco MC");
      dlutility::SetFont(hd, 42, 0.08);
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
      h_truth_match_dphi_rangemin[0][i]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{d#Delta#phi}{dN_{pair}}");
  
      h_truth_match_dphi_rangemin[0][i]->SetMaximum(10);
      h_truth_match_dphi_rangemin[0][i]->SetMinimum(0);
      h_truth_match_dphi_rangemin[0][i]->Draw("p");
      h_reco_dphi_rangemin[0][i]->Draw("p same");

      dlutility::DrawSPHENIXpp(0.22, 0.84);
      dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[binranges[1]], ipt_bins[binranges[1+1]]), 0.22, 0.69);
      dlutility::drawText(Form("%2.1f GeV #leq p_{T,2} < %2.1f GeV ", ipt_bins[binrangesmin[i]], ipt_bins[binrangesmin[i+1]]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
      TLegend *leg1 = new TLegend(0.20, 0.36, 0.63, 0.56);
      leg1->SetLineWidth(0);
      leg1->SetTextSize(0.04);
      leg1->SetTextFont(42);
      leg1->AddEntry(h_truth_match_dphi_rangemin[0][i], "PYTHIA-8 truth","p");
      leg1->AddEntry(h_reco_dphi_rangemin[0][i], "PYTHIA-8 reco","p");
      leg1->Draw("same");  

      c_cor->cd(2);
      gPad->SetTopMargin(0.05);

      TH1D *hd = (TH1D*) h_truth_match_dphi_rangemin[0][i]->Clone();
      hd->Divide(h_reco_dphi_rangemin[0][i]);
      hd->SetTitle(";#Delta#phi; Truth MC / Reco MC");
      dlutility::SetFont(hd, 42, 0.08);
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

      dlutility::SetFont(h_data_dphi_range[0][irange], 42, 0.04);
      dlutility::SetLineAtt(h_truth_match_dphi_range[0][irange], kRed, 3, 1);
      h_data_dphi_range[0][irange]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{d#Delta#phi}{dN_{pair}}");


      h_data_dphi_range[0][irange]->Draw("p E1");
      h_truth_match_dphi_range[0][irange]->Draw("same hist C");
      h_data_dphi_range[0][irange]->Draw("p E1 same");
      h_data_dphi_range_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXpp(0.22, 0.84);
      dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[binranges[irange]], ipt_bins[binranges[irange+1]]), 0.22, 0.69);
      dlutility::drawText(Form("p_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
      TLegend *leg = new TLegend(0.22, 0.34, 0.4, 0.49);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_range[0][irange], "Data");
      leg->AddEntry(h_truth_match_dphi_range[0][irange], "PYTHIA-8");
      leg->Draw("same");
  
      c_sys1->cd(2);
      TH1D *hratstat = (TH1D*) h_data_dphi_range[0][irange]->Clone();
      dlutility::SetFont(hratstat, 42, 0.09);
      dlutility::SetLineAtt(hratstat, kBlack, 1, 1);
      dlutility::SetMarkerAtt(hratstat, kBlack, 1, 8);

      dlutility::SetLineAtt(  h_ratio_dphi_range_sys[irange], kBlack, 1, 1);
      dlutility::SetMarkerAtt(  h_ratio_dphi_range_sys[irange], kBlack, 1, 8);
      h_ratio_dphi_range_sys[irange]->SetFillColorAlpha(kBlack, 0.3);
  
      hratstat->Divide(h_truth_match_dphi_range[0][irange]);
      hratstat->SetMaximum(1.5);
      hratstat->SetMinimum(0.5);

      hratstat->SetTitle(";#Delta#phi; Data/PYTHIA-8");
      hratstat->Draw("p E1");
      h_ratio_dphi_range_sys[irange]->Draw("p E2 same");
      TLine *li = new TLine(TMath::Pi()*3/4., 1, TMath::Pi(), 1);
      dlutility::SetLineAtt(li, kBlack, 1, 4);
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

      dlutility::SetFont(h_data_dphi_rangemin[0][irange], 42, 0.04);
      dlutility::SetLineAtt(h_truth_match_dphi_rangemin[0][irange], kRed, 3, 1);
      h_data_dphi_rangemin[0][irange]->SetTitle(";#Delta#phi; #frac{1}{N_{pair}}#frac{d#Delta#phi}{dN_{pair}}");
      h_data_dphi_rangemin[0][irange]->Draw("p E1");
      h_truth_match_dphi_rangemin[0][irange]->Draw("same hist C");
      h_data_dphi_rangemin[0][irange]->Draw("p E1 same");
      h_data_dphi_rangemin_sys[irange]->Draw("E2 p same");

      dlutility::DrawSPHENIXpp(0.22, 0.84);
      dlutility::drawText("anti-k_{T} R = 0.4", 0.22, 0.74);
      dlutility::drawText(Form("%2.1f GeV #leq p_{T,1} < %2.1f GeV ", ipt_bins[binranges[1]], ipt_bins[binranges[1+1]]), 0.22, 0.69);
      dlutility::drawText(Form("%2.1f GeV #leq p_{T,2} < %2.1f GeV ", ipt_bins[binrangesmin[irange]], ipt_bins[binrangesmin[irange+1]]), 0.22, 0.64);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.22, 0.59);
      TLegend *leg = new TLegend(0.22, 0.34, 0.4, 0.49);
      leg->SetLineWidth(0);
      leg->SetTextSize(0.04);
      leg->SetTextFont(42);
      leg->AddEntry(h_data_dphi_rangemin[0][irange], "Data");
      leg->AddEntry(h_truth_match_dphi_rangemin[0][irange], "PYTHIA-8");
      leg->Draw("same");
  
      c_sys1->cd(2);
      TH1D *hratstat = (TH1D*) h_data_dphi_rangemin[0][irange]->Clone();
      dlutility::SetFont(hratstat, 42, 0.09);
      dlutility::SetLineAtt(hratstat, kBlack, 1, 1);
      dlutility::SetMarkerAtt(hratstat, kBlack, 1, 8);

      dlutility::SetLineAtt(  h_ratio_dphi_rangemin_sys[irange], kBlack, 1, 1);
      dlutility::SetMarkerAtt(  h_ratio_dphi_rangemin_sys[irange], kBlack, 1, 8);
      h_ratio_dphi_rangemin_sys[irange]->SetFillColorAlpha(kBlack, 0.3);
  
      hratstat->Divide(h_truth_match_dphi_rangemin[0][irange]);
      hratstat->SetMaximum(1.5);
      hratstat->SetMinimum(0.5);

      hratstat->SetTitle(";#Delta#phi; Data/PYTHIA-8");
      hratstat->Draw("p E1");
      h_ratio_dphi_rangemin_sys[irange]->Draw("p E2 same");
      TLine *li = new TLine(TMath::Pi()*3/4., 1, TMath::Pi(), 1);
      dlutility::SetLineAtt(li, kBlack, 1, 4);
      li->Draw("same");
      c_sys1->Print(Form("dphi_rangemin_%d.pdf", irange));
      c_sys1->Print(Form("dphi_rangemin_%d.png", irange));
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
