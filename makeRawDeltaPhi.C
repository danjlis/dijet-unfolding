#include "dlUtility.h"
#include "read_binning.h"

void makeRawDeltaPhi(const std::string configfile = "binning.config")
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  std::string data_file = "TNTUPLE_DIJET_v6_2_ana468_2024p012_v001_gl10-00047289-00048291.root";
  std::string j10_file = "TREE_MATCH_v6_10_new_ProdA_2024-00000021.root";
  std::string j20_file = "TREE_MATCH_v6_20_new_ProdA_2024-00000021.root";
  std::string j30_file = "TREE_MATCH_v6_30_new_ProdA_2024-00000021.root";

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

  float pt1_data;
  float pt2_data;
  float dphi_data;

  TFile *fin[3];
  fin[0] = new TFile(j10_file.c_str(), "r");
  fin[1] = new TFile(j20_file.c_str(), "r");
  fin[2] = new TFile(j30_file.c_str(), "r");

  TNtuple *tn[3];
  for (int i = 0 ; i < 3; i++)
    {
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
    }
  float cs_10 = (3.646e-6);
  float cs_20 = 5.4067742e-8;
  float cs_30 = (2.505e-9);
  float scale_factor[3];
  scale_factor[0] = (9987000./9988000)*cs_10/cs_30;
  scale_factor[1] = (9987000./9954000.)*cs_20/cs_30; 
  //scale_factor[1] = (3.646e-6)/(2.505e-9);//4.197e-2;
  scale_factor[2] = 1;//4.197e-2;
  
  read_binning rb(configfile.c_str());

  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();
  Double_t dphicuttruth = dphicut;//TMath::Pi()/2.;
  
  TFile *find = new TFile(data_file.c_str(), "r");
  TNtuple *tnd  = (TNtuple*) find->Get("tn_dijet");;
  if (!tnd)
    {
      std::cout << " no data "<< std::endl;
    }
  tnd->SetBranchAddress("pt1_reco", &pt1_data);
  tnd->SetBranchAddress("pt2_reco", &pt2_data);
  tnd->SetBranchAddress("dphi_reco", &dphi_data);


  const int nbinsdphi  = 16;
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

  TH1D *h_data_dphi = new TH1D("h_data_dphi",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_truth_dphi = new TH1D("h_truth_dphi",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_reco_dphi = new TH1D("h_reco_dphi",";#Delta#phi;1/N",nbinsdphi, idphi_bins);

  TH1D *h_truth_match_dphi = new TH1D("h_truth_match_dphi",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_reco_match_dphi = new TH1D("h_reco_match_dphi",";#Delta#phi;1/N",nbinsdphi, idphi_bins);

  TProfile *h_sim_match_ddphi = new TProfile("h_sim_match_ddphi", ";<p_{T}> [GeV]; d#Delta#phi", 25, 10, 60, "s");
  TProfile2D *h2_sim_match_ddphi = new TProfile2D("h2_sim_match_ddphi", ";p_{T1} [GeV] ;p_{T2} [GeV]; d#Delta#phi", 4, 20, 60, 10, 10, 60, "s");

  TH3D *h_truth_pt1pt2dphi = new TH3D("h_truth_pt1pt2dphi",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_pt1pt2dphi = new TH3D("h_reco_pt1pt2dphi",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);

  TH3D *h_truth_match_pt1pt2dphi = new TH3D("h_truth_match_pt1pt2dphi",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphi = new TH3D("h_reco_match_pt1pt2dphi",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphitruth = new TH3D("h_reco_match_pt1pt2dphitruth",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);

  TH3D *h_data_pt1pt2dphi = new TH3D("h_data_pt1pt2dphi",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  // Data

  TH1D *h_data_dphi_counts = new TH1D("h_data_dphi_counts",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_truth_dphi_counts = new TH1D("h_truth_dphi_counts",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_reco_dphi_counts = new TH1D("h_reco_dphi_counts",";#Delta#phi;1/N",nbinsdphi, idphi_bins);

  TH1D *h_truth_match_dphi_counts = new TH1D("h_truth_match_dphi_counts",";#Delta#phi;1/N",nbinsdphi, idphi_bins);
  TH1D *h_reco_match_dphi_counts = new TH1D("h_reco_match_dphi_counts",";#Delta#phi;1/N",nbinsdphi, idphi_bins);

  TProfile *h_sim_match_ddphi_counts = new TProfile("h_sim_match_ddphi_counts", ";<p_{T}> [GeV]; d#Delta#phi", 25, 10, 60, "s");
  TProfile2D *h2_sim_match_ddphi_counts = new TProfile2D("h2_sim_match_ddphi_counts", ";p_{T1} [GeV] ;p_{T2} [GeV]; d#Delta#phi", 4, 20, 60, 10, 10, 60, "s");

  TH3D *h_truth_pt1pt2dphi_counts = new TH3D("h_truth_pt1pt2dphi_counts",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_pt1pt2dphi_counts = new TH3D("h_reco_pt1pt2dphi_counts",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);

  TH3D *h_truth_match_pt1pt2dphi_counts = new TH3D("h_truth_match_pt1pt2dphi_counts",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphi_counts = new TH3D("h_reco_match_pt1pt2dphi_counts",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  TH3D *h_reco_match_pt1pt2dphi_countstruth = new TH3D("h_reco_match_pt1pt2dphi_countstruth",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);

  TH3D *h_data_pt1pt2dphi_counts = new TH3D("h_data_pt1pt2dphi_counts",";p_{T,1, smear};p_{T,2, smear}", nbins, ipt_bins, nbins, ipt_bins,nbinsdphi, idphi_bins);
  // Data

  TH1D *h_correlated_counts_leading_reco[3];
  TH1D *h_correlated_counts_subleading_reco[3];
  TH1D *h_all_counts_leading_reco[3];
  TH1D *h_all_counts_subleading_reco[3];

    TH1D *h_correlated_counts_leading_truth[3];
  TH1D *h_correlated_counts_subleading_truth[3];
  TH1D *h_all_counts_leading_truth[3];
  TH1D *h_all_counts_subleading_truth[3];

  for (int i = 0; i < 3; i++)
    {
      h_correlated_counts_leading_reco[i] = new TH1D(Form("h_correlated_counts_leading_reco_%d", i), ";;", nbinsdphi, idphi_bins);
      h_correlated_counts_subleading_reco[i] = new TH1D(Form("h_correlated_counts_subleading_reco_%d", i), ";;", nbinsdphi, idphi_bins);
      h_all_counts_leading_reco[i] = new TH1D(Form("h_all_counts_leading_reco_%d", i), ";;", nbinsdphi, idphi_bins);
      h_all_counts_subleading_reco[i] = new TH1D(Form("h_all_counts_subleading_reco_%d", i), ";;", nbinsdphi, idphi_bins);
      h_correlated_counts_leading_truth[i] = new TH1D(Form("h_correlated_counts_leading_truth_%d", i), ";;", nbinsdphi, idphi_bins);
      h_correlated_counts_subleading_truth[i] = new TH1D(Form("h_correlated_counts_subleading_truth_%d", i), ";;", nbinsdphi, idphi_bins);
      h_all_counts_leading_truth[i] = new TH1D(Form("h_all_counts_leading_truth_%d", i), ";;", nbinsdphi, idphi_bins);
      h_all_counts_subleading_truth[i] = new TH1D(Form("h_all_counts_subleading_truth_%d", i), ";;", nbinsdphi, idphi_bins);

    }
  int entriesd = tnd->GetEntries();

  TF1 *fJES = new TF1("fJES", "pol4", 0, 100);
  fJES->SetParameters(0.51, 1.42, -8.57e-3, 1.69e-4, -1.13e-6); 

  TF1 *fgaus = new TF1("fgaus", "gaus");
  fgaus->SetRange(-0.5, 0.5);
  Double_t JES_sys = rb.get_jes_sys();
  Double_t JER_sys = rb.get_jer_sys();
  std::cout << "JES = " << JES_sys << std::endl;
  std::cout << "JER = " << JER_sys << std::endl;
  float width = 0.1 + JER_sys;
  fgaus->SetParameters(1, 0, width);

   // Vertex Reweight
  Int_t vtx_sys = rb.get_vtx_sys();
  TFile *fvtx = new TFile("vertex_reweight.root","r");
  TH1D *h_mbd_reweight = (TH1D*) fvtx->Get("h_mbd_reweight");

  std::vector<std::pair<float, float>> vertex_scales;

  for (int ib = 0; ib < h_mbd_reweight->GetNbinsX(); ib++)
    {
      vertex_scales.push_back(std::make_pair(h_mbd_reweight->GetBinLowEdge(ib+1) + h_mbd_reweight->GetBinWidth(ib+1), h_mbd_reweight->GetBinContent(ib+1)));
    }
 
  for (int i = 0; i < entriesd; i++)
    {
      tnd->GetEntry(i);
      float maxi = std::max(pt1_data, pt2_data);
      float mini = std::min(pt1_data, pt2_data);

      float pt1_data_bin = nbins;
      float pt2_data_bin = nbins;

      float es1 = fJES->Eval(pt1_data);
      float es2 = fJES->Eval(pt2_data);

      if (es1 > ipt_bins[nbins] ) continue;
      for (int ib = 0; ib < nbins; ib++)
	{

	  if ( es1 < ipt_bins[ib+1] && es1 >= ipt_bins[ib])
	    {
	      pt1_data_bin = ib;
	    }
	  if ( es2 < ipt_bins[ib+1] && es2 >= ipt_bins[ib])
	    {
	      pt2_data_bin = ib;
	    }
	}	  
	 
      bool data_good = (maxi >= reco_leading_cut && mini >= reco_subleading_cut && dphi_data > dphicut);

      if (!data_good) continue;
      
      h_data_pt1pt2dphi->Fill(es1, es2, dphi_data);

      h_data_dphi->Fill(dphi_data);

    }

  // Sim
  for (int isample = 0; isample < 3; isample++)
    {
      std::cout << "Sample " << isample << std::endl;
      int entries2 = tn[isample]->GetEntries();
      for (int i = 0; i < entries2; i++)
	{
	  tn[isample]->GetEntry(i);
	  
	  if (maxpttruth[isample] < sample_boundary[isample] || maxpttruth[isample] >= sample_boundary[isample+1]) continue;

	  float event_scale =  scale_factor[isample];
	  // Vertex Rewieghting
	  for (int ib = 0; ib < vertex_scales.size(); ib++)
	    {
	      if (mbd_vertex[isample] < vertex_scales.at(ib).first)
		{
		  event_scale *= vertex_scales.at(ib).second;
		  //if (i < 100) std:cout << "found z = " << vertex_scales.at(ib).first << " " <<vertex_scales.at(ib).second<<std::endl;
		  break;
		}
	    }
	  
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
	  float es1 = fJES->Eval(max_reco);
	  float es2 = fJES->Eval(min_reco);
	  
	  double smear1 = fgaus->GetRandom();
	  double smear2 = fgaus->GetRandom();

	  if (JES_sys != 0)
	    {
	      es1 = es1 + (JES_sys + smear1)*e1;
	      es2 = es2 + (JES_sys + smear2)*e2; 
	    }
	  else
	    {
	      es1 = es1 + smear1*e1;
	      es2 = es2 + smear2*e2; 	  
	    }
	  
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
	  
	  bool truth_good = (e1 >= truth_leading_cut && e2 >= truth_subleading_cut && fabs(dphi_truth[isample]) > dphicuttruth);
	  bool reco_good = (maxi >= measure_leading_cut && mini >= measure_subleading_cut && fabs(dphi_reco[isample]) > dphicut);
	  if (match[isample])
	    {
	      float ddphi = dphi_truth[isample] - dphi_reco[isample];
	      if (ddphi > TMath::Pi()) ddphi -= 2*TMath::Pi();
	      if (ddphi < -1*TMath::Pi()) ddphi += 2*TMath::Pi();
	      h_sim_match_ddphi->Fill((max_truth +  min_truth)/2., ddphi, event_scale);
	    }
	  if (truth_good)
	    {
	      h_truth_dphi->Fill(dphi_truth[isample], event_scale);
	      h_truth_pt1pt2dphi->Fill(max_truth, min_truth, dphi_truth[isample], event_scale);
	    }
	  if (reco_good)
	    {
	      h_reco_dphi->Fill(dphi_reco[isample], event_scale);
	      h_reco_pt1pt2dphi->Fill(maxi, mini, dphi_reco[isample], event_scale);
	    }

	  if (match[isample] && reco_good && truth_good)
	    {
	      h_truth_match_dphi->Fill(dphi_truth[isample], event_scale);
	      h_reco_match_dphi->Fill(dphi_reco[isample], event_scale);
	      h_truth_match_pt1pt2dphi->Fill(max_truth, min_truth, dphi_truth[isample], event_scale);
	      h_reco_match_pt1pt2dphi->Fill(maxi, mini, dphi_reco[isample], event_scale);
	      h_reco_match_pt1pt2dphitruth->Fill(maxi, mini, dphi_truth[isample], event_scale);
	      h_truth_dphi->Fill(dphi_truth[isample], event_scale);
	      h_reco_dphi->Fill(dphi_reco[isample], event_scale);
	      h_truth_match_dphi->Fill(dphi_truth[isample], event_scale);
	      h_reco_match_dphi->Fill(dphi_reco[isample], event_scale);

	      h_truth_match_dphi_counts->Fill(dphi_truth[isample]);
	      h_reco_match_dphi_counts->Fill(dphi_reco[isample]);
	      h_truth_match_pt1pt2dphi_counts->Fill(max_truth, min_truth, dphi_truth[isample]);
	      h_reco_match_pt1pt2dphi_counts->Fill(maxi, mini, dphi_reco[isample]);
	      h_reco_match_pt1pt2dphi_countstruth->Fill(maxi, mini, dphi_truth[isample]);
	      h_truth_dphi_counts->Fill(dphi_truth[isample]);
	      h_reco_dphi_counts->Fill(dphi_reco[isample]);
	      h_truth_match_dphi_counts->Fill(dphi_truth[isample]);
	      h_reco_match_dphi_counts->Fill(dphi_reco[isample]);

	      bool same_dphi = false;
	      for (int ip = 0; ip < nbinsdphi; ip++)
		{
		  
		  if (dphi_truth[isample] >= idphi_bins[ip] && dphi_truth[isample] < idphi_bins[ip+1] && dphi_reco[isample] >= idphi_bins[ip] && dphi_reco[isample] < idphi_bins[ip+1])
		    {
		      same_dphi = true;
		      break;
		    }
		}

	      for (int ilead = 0; ilead < 3; ilead++)
		{
		  int same_leading = 0;
		  if (!(max_truth < ipt_bins[binranges[ilead]] || max_truth >= ipt_bins[binranges[ilead+1]]))
		    {
		      same_leading++;
		      h_all_counts_leading_truth[ilead]->Fill(dphi_truth[isample]);
		    }
		  if (!(maxi < ipt_bins[binranges[ilead]] || maxi >= ipt_bins[binranges[ilead+1]]))
		    {
		      same_leading++;
		      h_all_counts_leading_reco[ilead]->Fill(dphi_reco[isample]);
		    }
		  if (same_dphi && same_leading == 2)
		    {
		      h_correlated_counts_leading_reco[ilead]->Fill(dphi_reco[isample]);
		      h_correlated_counts_leading_truth[ilead]->Fill(dphi_truth[isample]);
		    }

		    for (int isub = 0; isub < 3; isub++)
		    {
		      int same_subleading = 0;
		      if (!(min_truth < ipt_bins[binrangesmin[isub]] || min_truth >= ipt_bins[binrangesmin[isub+1]]))
			{
			  same_subleading++;
			  h_all_counts_subleading_truth[isub]->Fill(dphi_truth[isample]);
			}
		      if (!(mini < ipt_bins[binrangesmin[isub]] || mini >= ipt_bins[binrangesmin[isub+1]]))
			{
			  same_subleading++;
			  h_all_counts_subleading_reco[isub]->Fill(dphi_reco[isample]);
			}
		      if (same_dphi && same_subleading == 2)
			{
			  h_correlated_counts_subleading_reco[isub]->Fill(dphi_reco[isample]);
			  h_correlated_counts_subleading_truth[isub]->Fill(dphi_truth[isample]);
			}
		    }
		}
	      //h_sim_match_ddphi->Fill((max_truth +  min_truth)/2., dphi_truth[isample] - dphi_reco[isample], event_scale);
	      h2_sim_match_ddphi->Fill(max_truth,  min_truth, dphi_truth[isample] - dphi_reco[isample], event_scale);
	    }
	}
    }
  TString outpath = "dphi_hists.root";

  if (JES_sys > 0)
    {
      outpath = "dphi_hists_posJES.root";
    }
  if (JES_sys < 0)
    {
      outpath = "dphi_hists_negJES.root";
    }
  if (JER_sys > 0)
    {
      outpath = "dphi_hists_posJER.root";
    }
  if (JER_sys < 0)
    {
      outpath = "dphi_hists_negJER.root";
    }

  TFile *fout = new TFile(outpath.Data(), "recreate");

  h_sim_match_ddphi->Write();
  h2_sim_match_ddphi->Write();
  h_data_pt1pt2dphi->Write();
  h_reco_match_pt1pt2dphi->Write();
  h_reco_match_pt1pt2dphitruth->Write();
  h_truth_match_pt1pt2dphi->Write();

  h_reco_pt1pt2dphi->Write();
  h_truth_pt1pt2dphi->Write();
  h_reco_dphi->Write();
  h_truth_dphi->Write();

  h_reco_match_dphi->Write();
  h_truth_match_dphi->Write();

  h_sim_match_ddphi_counts->Write();
  h2_sim_match_ddphi_counts->Write();
  h_data_pt1pt2dphi_counts->Write();
  h_reco_match_pt1pt2dphi_counts->Write();
  h_reco_match_pt1pt2dphi_countstruth->Write();
  h_truth_match_pt1pt2dphi_counts->Write();

  h_reco_pt1pt2dphi_counts->Write();
  h_truth_pt1pt2dphi_counts->Write();
  h_reco_dphi_counts->Write();
  h_truth_dphi_counts->Write();

  h_reco_match_dphi_counts->Write();
  h_truth_match_dphi_counts->Write();

    for (int i = 0; i < 3; i++) {
      h_correlated_counts_leading_reco[i]->Write();
      h_correlated_counts_subleading_reco[i]->Write();
      h_correlated_counts_leading_truth[i]->Write();
      h_correlated_counts_subleading_truth[i]->Write();

      h_all_counts_leading_reco[i]->Write();
      h_all_counts_subleading_reco[i]->Write();
      h_all_counts_leading_truth[i]->Write();
      h_all_counts_subleading_truth[i]->Write();

    }
  std::cout << h_reco_dphi->GetMean() <<std::endl;  
  fout->Close();
}
