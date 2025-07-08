#include "dlUtility.h"
#include "read_binning.h"
void makeRawXJ()
{
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  std::string mb_file = "TREE_MATCH_v4_8_new_ProdA_2024-00000021.root";
  std::string j10_file = "TREE_MATCH_v6_10_new_ProdA_2024-00000021.root";
  std::string j30_file = "TREE_MATCH_v6_30_new_ProdA_2024-00000021.root";
  std::string data_file = "TNTUPLE_DIJET_v6_1_ana468_2024p012_v001_gl10-00047352-00047733.root";//TNTUPLE_DIJET_v4_2_ana450_2024p009_gl10-00047352-00047733.root";

  float maxpttruth[3];
  float pt1_truth[3];
  float pt2_truth[3];
  float dphi_truth[3];
  float pt1_reco[3];
  float pt2_reco[3];
  float dphi_reco[3];
  float match[3];
  float pt1_data;
  float pt2_data;
  float dphi_data;

  float sample_boundary_goal[4];
  sample_boundary_goal[0] = 3;
  sample_boundary_goal[1] = 13;
  sample_boundary_goal[2] = 29.99;
  sample_boundary_goal[3] = 100;

  float sample_boundary[4];
  sample_boundary[0] = 0;
  sample_boundary[1] = 0;
  sample_boundary[2] = 0;
  sample_boundary[3] = 100;
  
  TFile *find = new TFile(data_file.c_str(), "r");
  TNtuple *tnd  = (TNtuple*) find->Get("tn_dijet");;
  if (!tnd)
    {
      std::cout << " no data "<< std::endl;
    }
  tnd->SetBranchAddress("pt1_reco", &pt1_data);
  tnd->SetBranchAddress("pt2_reco", &pt2_data);
  tnd->SetBranchAddress("dphi_reco", &dphi_data);


  TFile *fin[3];
  fin[0] = new TFile(mb_file.c_str(), "r");
  fin[1] = new TFile(j10_file.c_str(), "r");
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
      tn[i]->SetBranchAddress("matched", &match[i]);
    }

  float scale_factor[3];
  scale_factor[0] = 1;
  scale_factor[1] = (9987000./10000000.)*(3.646e-6)/(2.505e-9);//4.197e-2;
  scale_factor[2] = 1;//4.197e-2;
  
  read_binning rb("binning.config");

  Int_t read_nbins = rb.get_nbins();
  Int_t read_bbin = rb.get_bbins();
  Double_t read_fixed = rb.get_fixed();
  Double_t read_minimum = rb.get_minimum();

  const int nbins = read_nbins;

  float ipt_bins[nbins+1];
  float ixj_bins[nbins+1];

  float alpha = TMath::Power(read_fixed/read_minimum, 1/(float)read_bbin);
  float bin_min = read_minimum;
  float bin_max = read_minimum*TMath::Power(alpha, nbins );

  float bin_xj10 = 1.0;
  float bin_xj1 = 1.0*(bin_min/bin_max);
  std::cout << "Alpha = " << alpha <<std::endl;

  int measure_leading_bin = 0;
  int measure_subleading_bin = 0;
  int truth_leading_bin = 0;
  int truth_subleading_bin = 0;
  float truth_leading_cut = 0;
  float reco_leading_cut = 0;
  float reco_meas_leading_cut = 0;
  float truth_subleading_cut = 0;
  float reco_subleading_cut = 0;
  float reco_meas_subleading_cut = 0;

  float measure_leading_goal = 20;
  float truth_leading_goal = 13;
  float truth_subleading_goal = 0;

  // Make binning
  for (int i = 0; i < nbins+1; i++)
    {
      float ipt = bin_min*TMath::Power(alpha, (float)i);
      float ixj = bin_xj1*TMath::Power(alpha, (float)i);
      ipt_bins[i] = ipt;
      ixj_bins[i] = ixj;

      if (measure_leading_bin == 0 && ipt >= measure_leading_goal) measure_leading_bin = i;
      if (truth_leading_bin == 0 && ipt >= truth_leading_goal) truth_leading_bin = i;
      if (truth_subleading_bin == 0 && ipt >= truth_subleading_goal) truth_subleading_bin = i;
      std::cout << i << " : " <<  ipt << " -- " << ixj <<  std::endl;
    }

  truth_leading_cut=ipt_bins[truth_leading_bin];
  truth_subleading_cut=ipt_bins[truth_subleading_bin];

  reco_leading_cut=ipt_bins[truth_leading_bin + 4];
  reco_subleading_cut=ipt_bins[truth_subleading_bin + 4];

  reco_meas_leading_cut=ipt_bins[truth_leading_bin + 4];
  reco_meas_subleading_cut=ipt_bins[truth_subleading_bin + 4];
  measure_subleading_bin = truth_subleading_bin + 4;

  float dphicut = 3*TMath::Pi()/4.;
  float dphicuttruth = 3*TMath::Pi()/4.;

  TH1D *h_data_fill_xj = new TH1D("h_data_fill_xj",";x_{J};1/N",nbins, ixj_bins);
  TH1D *h_reco_fill_xj = new TH1D("h_reco_fill_xj",";x_{J};1/N",nbins, ixj_bins);
  TH1D *h_data_proj_xj = new TH1D("h_data_proj_xj",";x_{J};1/N",nbins, ixj_bins);
  TH1D *h_reco_proj_xj = new TH1D("h_reco_proj_xj",";x_{J};1/N",nbins, ixj_bins);

  TH1D *h_flat_data_pt1pt2 = new TH1D("h_data_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  TH1D *h_flat_reco_pt1pt2 = new TH1D("h_reco_flat_pt1pt2",";p_{T,1, smear} + p_{T,2, smear}", nbins*nbins, 0, nbins*nbins);
  
  TH1D *h_data_xj = new TH1D("h_data_xj",";x_{J};1/N",20, 0, 1);
  TH1D *h_reco_xj = new TH1D("h_reco_xj",";x_{J};1/N",20, 0, 1);
  TH1D *h_match_reco_xj = new TH1D("h_match_reco_xj",";x_{J};1/N",20, 0, 1);
  TH1D *h_truth_xj = new TH1D("h_truth_xj",";x_{J};1/N",20, 0, 1);
  TH1D *h_match_truth_xj = new TH1D("h_match_truth_xj",";x_{J};1/N",20, 0, 1);
  // Data
  int entriesd = tnd->GetEntries();
  for (int i = 0; i < entriesd; i++)
    {
      tnd->GetEntry(i);
      float maxi = std::max(pt1_data, pt2_data);
      float mini = std::min(pt1_data, pt2_data);

      int pt1_reco_bin=nbins;
      int pt2_reco_bin=nbins;

      for (int ib = 0; ib < nbins; ib++)
	{

	  if ( maxi < ipt_bins[ib+1] && maxi >= ipt_bins[ib])
	    {
	      pt1_reco_bin = ib;
	    }
	  if ( mini < ipt_bins[ib+1] && mini >= ipt_bins[ib])
	    {
	      pt2_reco_bin = ib;
	    }
	}

      bool data_good = (maxi >= reco_meas_leading_cut && mini >= reco_meas_subleading_cut && dphi_data > dphicut);

      if (!data_good) continue;
      
      h_data_xj->Fill(mini/maxi);
      h_flat_data_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin);
      h_flat_data_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin);
      h_data_fill_xj->Fill(mini/maxi);

    }

  // Sim
  for (int isample = 1; isample < 3; isample++)
    {
      std::cout << "Sample " << isample << std::endl;
      int entries2 = tn[isample]->GetEntries();
      for (int i = 0; i < entries2; i++)
	{
	  tn[isample]->GetEntry(i);
	  
	  if (maxpttruth[isample] < sample_boundary[isample] || maxpttruth[isample] >= sample_boundary[isample+1]) continue;
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
	    
	  float maxi = std::max(pt1_reco[isample], pt2_reco[isample]);
	  float mini = std::min(pt1_reco[isample], pt2_reco[isample]);
	  float e1 = max_truth;
	  float e2 = min_truth;
	  float es1 = max_reco;
	  float es2 = min_reco;
	  int pt1_reco_bin=nbins;
	  int pt2_reco_bin=nbins;
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
	  

	  
	  bool truth_good = (e1 >= truth_leading_cut && e2 >= truth_subleading_cut && dphi_truth[isample] > dphicuttruth);
	  bool reco_good = (maxi >= reco_meas_leading_cut && mini >= reco_meas_subleading_cut && dphi_reco[isample] > dphicut);

	  if (!truth_good && !reco_good) continue;

	  if (truth_good)
	    {
	      h_truth_xj->Fill(e2/e1, scale_factor[isample]);
	    }
	  if (reco_good)
	    {
	      h_reco_xj->Fill(mini/maxi, scale_factor[isample]);	
	    }

	  if (match[isample] && reco_good && truth_good)
	    {
	      h_match_truth_xj->Fill(e2/e1, scale_factor[isample]);
	      h_match_reco_xj->Fill(mini/maxi, scale_factor[isample]);	

	      h_flat_reco_pt1pt2->Fill(pt1_reco_bin + nbins*pt2_reco_bin, scale_factor[isample]);
	      h_flat_reco_pt1pt2->Fill(pt2_reco_bin + nbins*pt1_reco_bin, scale_factor[isample]);
	      h_reco_fill_xj->Fill(mini/maxi, scale_factor[isample]);

	    }
	}
    }
  std::cout << " method versus manual " << std::endl;
  std::cout << h_match_reco_xj->Integral(0, -1, "width") << std::endl;
  std::cout << h_match_reco_xj->Integral(1, h_match_reco_xj->GetNbinsX(), "width") << std::endl;
  std::cout << h_match_reco_xj->Integral(1, h_match_reco_xj->GetNbinsX()) << std::endl;

  double   integral = 0;
  for  (int ibin = 0; ibin < h_match_reco_xj->GetNbinsX() ; ibin ++)
    {
      integral = h_match_reco_xj->GetBinWidth(ibin+1)*h_match_reco_xj->GetBinContent(ibin+1);
    }
  std::cout << integral << std::endl;

  h_match_truth_xj->Scale(1./h_match_truth_xj->Integral(0, -1),"width");
  h_match_reco_xj->Scale(1./h_match_reco_xj->Integral(0, -1),"width");
  h_data_xj->Scale(1./h_data_xj->Integral(0, -1),"width");
  

  
  TCanvas *call = new TCanvas("call","call", 700, 700);

  dlutility::SetMarkerAtt(h_data_xj, kAzure - 6, 1, 8);
  dlutility::SetLineAtt(h_data_xj, kAzure - 6, 1, 1);

  dlutility::SetMarkerAtt(h_match_reco_xj, kRed - 2, 1, 8);
  dlutility::SetLineAtt(h_match_reco_xj, kRed - 2, 1, 1);

  dlutility::SetMarkerAtt(h_match_truth_xj, kBlack, 0.9, 24);
  dlutility::SetLineAtt(h_match_truth_xj, kBlack, 1, 1);
  h_match_truth_xj->SetTitle(";x_{J};#frac{1}{N_{pairs}}#frac{dN_{pairs}}{dx_{J}}");
  h_match_truth_xj->GetYaxis()->SetTitleOffset(1.6);
  h_match_truth_xj->Draw("p");
  h_match_reco_xj->Draw("same p");
  h_data_xj->Draw("same p");

  dlutility::DrawSPHENIXpp(0.22, 0.87);
  dlutility::drawText(Form("p_{T,1} > %2.1f GeV", reco_leading_cut), 0.22, 0.77);
  dlutility::drawText(Form("p_{T,2} > %2.1f GeV", reco_subleading_cut), 0.22, 0.72);
  dlutility::drawText("Truth Matched", 0.22, 0.67);
  TLegend *leg = new TLegend(0.22, 0.4, 0.37, 0.62);
  leg->SetLineWidth(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_match_truth_xj, "Pythia8 Truth","p");
  leg->AddEntry(h_match_reco_xj, "Pythia8 Reco","p");
  leg->AddEntry(h_data_xj, "Data","p");
  leg->Draw("same");
  call->Print("h_xj_all.pdf");
  call->Print("h_xj_all.png");


  // projection

  h_flat_reco_pt1pt2->Scale(.5);
  h_flat_data_pt1pt2->Scale(.5);

  TH2D *h_pt1pt2_reco = new TH2D("h_pt1pt2_reco", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);
  TH2D *h_pt1pt2_data = new TH2D("h_pt1pt2_data", ";p_{T1};p_{T2}", nbins, ipt_bins, nbins, ipt_bins);

  TH1D *h_xj_reco = new TH1D("h_xj_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xj_data = new TH1D("h_xj_data", ";x_{J};",nbins, ixj_bins);

  TH1D *h_xjunc_reco = new TH1D("h_xjunc_reco", ";x_{J};", nbins, ixj_bins);
  TH1D *h_xjunc_data = new TH1D("h_xjunc_data", ";x_{J};",nbins, ixj_bins);

  for (int ib = 0; ib < nbins*nbins; ib++)
    {
      int xbin = ib/nbins;
      int ybin = ib%nbins;
      
      int b = h_pt1pt2_reco->GetBin(xbin+1, ybin+1);

      h_pt1pt2_reco->SetBinContent(b, h_flat_reco_pt1pt2->GetBinContent(ib+1));
      h_pt1pt2_data->SetBinContent(b, h_flat_data_pt1pt2->GetBinContent(ib+1));
      h_pt1pt2_reco->SetBinError(b, h_flat_reco_pt1pt2->GetBinError(ib+1));
      h_pt1pt2_data->SetBinError(b, h_flat_data_pt1pt2->GetBinError(ib+1));
    }

  for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy < nbins; iy++)
	{
	  int bin = h_pt1pt2_reco->GetBin(ix+1, iy+1);

	  if (ix > iy)
	    {

	      h_pt1pt2_reco->SetBinContent(bin, h_pt1pt2_reco->GetBinContent(bin)*2.);
	      h_pt1pt2_data->SetBinContent(bin, h_pt1pt2_data->GetBinContent(bin)*2.);
	      
	      h_pt1pt2_reco->SetBinError(bin, h_pt1pt2_reco->GetBinError(bin)/sqrt(2));
	      h_pt1pt2_data->SetBinError(bin, h_pt1pt2_data->GetBinError(bin)/sqrt(2));

	    }
	  else if (ix < iy)
	    {	      
	      h_pt1pt2_reco->SetBinContent(bin, 0);
	      h_pt1pt2_data->SetBinContent(bin, 0);
	    }

	}
    }

  for (int ix = 0; ix < nbins; ix++)
    {
      for (int iy = 0; iy <= ix; iy++)
	{
	  int low =  iy - ix - 1;
	  int high = iy - ix + 1;

	  int xjbin_low = nbins + low + 1;
	  int xjbin_high = nbins + low + 2;
	  int bin = h_pt1pt2_reco->GetBin(ix+1, iy+1);
	  std::cout << ix << " -- " << iy << " --> " << low << " --> " << xjbin_high << "( " << h_xj_reco->GetBinCenter(xjbin_high)<<" ) " << " + " << xjbin_low << "( " << h_xj_reco->GetBinCenter(xjbin_low)<<" ) "<< std::endl;

	  if (ix < measure_leading_bin) continue;
	  if (iy < measure_subleading_bin) continue;
	  if (ix > nbins - 1) continue;
	  if (iy > nbins - 1) continue;

	  if (ix == iy)
	    {
	      h_xj_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low), h_pt1pt2_reco->GetBinContent(bin));
	      h_xjunc_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_reco->GetBinError(bin), 2));
	    }
	  else
	    {
	      h_xj_reco->Fill(h_xj_reco->GetBinCenter(xjbin_high), h_pt1pt2_reco->GetBinContent(bin)/2.);
	      h_xj_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low), h_pt1pt2_reco->GetBinContent(bin)/2.);
	      h_xjunc_reco->Fill(h_xj_reco->GetBinCenter(xjbin_high),TMath::Power( h_pt1pt2_reco->GetBinError(bin)/sqrt(2), 2));
	      h_xjunc_reco->Fill(h_xj_reco->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_reco->GetBinError(bin)/sqrt(2), 2));
	    }

	  if (ix == iy)
	    {
	      h_xj_data->Fill(h_xj_data->GetBinCenter(xjbin_low), h_pt1pt2_data->GetBinContent(bin));
	      h_xjunc_data->Fill(h_xj_data->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_data->GetBinError(bin), 2));
	    }
	  else
	    {
	      h_xj_data->Fill(h_xj_data->GetBinCenter(xjbin_high), h_pt1pt2_data->GetBinContent(bin)/2.);
	      h_xj_data->Fill(h_xj_data->GetBinCenter(xjbin_low), h_pt1pt2_data->GetBinContent(bin)/2.);
	      h_xjunc_data->Fill(h_xj_data->GetBinCenter(xjbin_high),TMath::Power( h_pt1pt2_data->GetBinError(bin)/sqrt(2), 2));
	      h_xjunc_data->Fill(h_xj_data->GetBinCenter(xjbin_low),TMath::Power( h_pt1pt2_data->GetBinError(bin)/sqrt(2), 2));
	    }	  
	}
    }
  for (int i = 1; i <= nbins; i++)
    {
      h_xj_reco->SetBinError(i, sqrt(h_xjunc_reco->GetBinError(i)));
      h_xj_data->SetBinError(i, sqrt(h_xjunc_data->GetBinError(i)));
    }


  h_xj_reco->Scale(1./h_xj_reco->Integral(0, -1),"width");
  h_xj_data->Scale(1./h_xj_data->Integral(0, -1),"width");

  std::cout << " method versus manual " << std::endl;
  std::cout << h_reco_fill_xj->Integral(0, -1, "width") << std::endl;
  std::cout << h_reco_fill_xj->Integral(1, h_reco_fill_xj->GetNbinsX(), "width") << std::endl;
  std::cout << h_reco_fill_xj->Integral(1, h_reco_fill_xj->GetNbinsX()) << std::endl;

  double   integral2 = 0;
  for  (int ibin = 0; ibin < h_reco_fill_xj->GetNbinsX() ; ibin ++)
    {
      integral2 = h_reco_fill_xj->GetBinWidth(ibin+1)*h_reco_fill_xj->GetBinContent(ibin+1);
    }
  std::cout << integral2 << std::endl;

  h_reco_fill_xj->Scale(1./h_reco_fill_xj->Integral(0, -1),"width");
  h_data_fill_xj->Scale(1./h_data_fill_xj->Integral(0, -1),"width");
  TCanvas *cproj = new TCanvas("cproj","cproj", 600, 600);

  dlutility::SetMarkerAtt(h_xj_reco, kRed - 2, 1, 25);
  dlutility::SetLineAtt(h_xj_reco, kRed - 2, 1, 1);

  dlutility::SetMarkerAtt(h_xj_data, kAzure - 6, 1, 25);
  dlutility::SetLineAtt(h_xj_data, kAzure - 6, 1, 1);

  dlutility::SetMarkerAtt(h_reco_fill_xj, kRed - 2, 1, 21);
  dlutility::SetLineAtt(h_reco_fill_xj, kRed - 2, 1, 1);

  dlutility::SetMarkerAtt(h_data_fill_xj, kAzure - 6, 1, 21);
  dlutility::SetLineAtt(h_data_fill_xj, kAzure - 6, 1, 1);



  h_match_truth_xj->Draw("p");
  h_match_reco_xj->Draw("same p");
  h_data_xj->Draw("same p");
  h_data_fill_xj->Draw("same");
  h_reco_fill_xj->Draw("same");

  h_xj_data->Draw("same");
  h_xj_reco->Draw("same");

}
