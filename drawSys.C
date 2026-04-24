//#include "systematics.h"

#include "dlUtility.h"
#include "read_binning.h"
#include "histo_opps.h"

#include "SystematicInfo.h"



void drawSys(const int cone_size = 4)
{

  gSystem->Load("SystematicInfo_cc.so");
  
  gStyle->SetOptStat(0);
  dlutility::SetyjPadStyle();

  read_binning rb("binning.config");

  Double_t first_xj = rb.get_first_xj();
  int first_bin = 0;
  Int_t read_nbins = rb.get_nbins();

  Double_t dphicut = rb.get_dphicut();

  const int nbins = read_nbins;
  const int nbins_pt = read_nbins+1;

  float ipt_bins[nbins_pt+1];
  float ixj_bins[nbins+1];
  double  dxj_bins[nbins+1];
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
      std::cout << measure_bins[ir] << " -- " <<  subleading_measure_bins[ir] << std::endl;
    }

  std::cout << "Truth1: " << truth_leading_cut << std::endl;
  std::cout << "Reco 1: " <<  reco_leading_cut << std::endl;
  std::cout << "Meas 1: " <<  measure_leading_cut << std::endl;
  std::cout << "Truth2: " <<  truth_subleading_cut << std::endl;
  std::cout << "Reco 2: " <<  reco_subleading_cut << std::endl;
  std::cout << "Meas 2: " <<  measure_subleading_cut << std::endl;

  
  const int niterations = 10;

  TH2D *h_average_xj[niterations];

  std::vector<std::pair<double, int>> *variation_average_xj[niterations][3];

  for (int iter = 0 ; iter < niterations; iter++)
    {
      h_average_xj[iter] = new TH2D(Form("h_average_xj_%d", iter),";range;variation;", 3, 0, 3, 3, 0, 3);
      for (int j = 0 ; j < 3; j++)
	{
	  variation_average_xj[iter][j] = new std::vector<std::pair<double, int>>();
	}
    }

  
  TFile *finu = new TFile(Form("%s/uncertainties/uncertainties_pp_r%02d_nominal.root", rb.get_code_location().c_str(), cone_size),"r");

  TProfile *h_xj_rms[niterations];
  for (int iter = 0; iter < niterations; ++iter)
    {
      h_xj_rms[iter] = (TProfile*) finu->Get(Form("hp_xj_%d", iter));
    }

  std::map<std::string, SystematicInfo> sys_list;
  std::map<std::string, std::string> sys_titles;
  std::map<std::string, std::vector<std::string>> sys_draw;

  std::ifstream infile("systematic_list.txt");
  std::string line;
  int isys = 0;
  int nom = 1;
  while (std::getline(infile, line)) {
    if (line.empty()) continue;

    std::stringstream ss(line);
    std::string name, group_str, val_str, path, config, title,color,marker,msize,lsize;

    // read directly into variables
    std::getline(ss, name, ',');
    std::getline(ss, group_str, ',');
    std::getline(ss, val_str, ',');
    std::getline(ss, path, ',');
    std::getline(ss, config, ',');
    std::getline(ss, title, ',');
    std::getline(ss, color, ',');
    std::getline(ss, marker, ',');
    std::getline(ss, msize, ',');
    std::getline(ss, lsize, ',');

    int value = std::stoi(val_str);
    int vmarker = std::stoi(marker);
    int vmsize = std::stoi(msize);
    int vlsize = std::stoi(lsize);

    int icolor = dlutility::parseRootColor(color);

    std::string filename = Form(path.c_str(), rb.get_code_location().c_str(), cone_size);    
    std::cout << "Name: " << name
	      << ", Int: " << value
	      << ", Path: " << filename
	      << ", Config: " << config
      	      << ", Title: " << title
	      << std::endl;


    SystematicInfo a_sys;

    a_sys.set_name(name.c_str());
    a_sys.set_nominal(nom);
    a_sys.set_reverse(value);
    a_sys.setNIterations(10);
    a_sys.setFile(filename.c_str());
    a_sys.setDrawingPars(icolor, icolor, vmarker, vmsize, vlsize);
    a_sys.setup(config.c_str());
    if(strcmp("HALF", group_str.c_str()) == 0)
      {
	std::cout << "setting half"<< std::endl;
	a_sys.setIsHalf(true);
      }
    a_sys.gethistos();
    nom = 0;
    
    if (!a_sys.isGood()) continue;

    sys_list[name] = a_sys;
    sys_titles[name] = title;
    
    if(strcmp("NA", group_str.c_str()) != 0)
	sys_draw[group_str.c_str()].push_back(name);
    
  }


  // JES/JER systematics
  std::map<std::string, std::array<std::array<TH1D*, 3>, 10>> h_sys;
  std::map<std::string, std::array<std::array<TH1D*, 3>, 10>> h_sys_neg;

  TH1D *h_total_sys_range[mbins][niterations];
  TH1D *h_total_sys_neg_range[mbins][niterations];

  double positive_systematics[nbins][mbins][niterations];
  double negative_systematics[nbins][mbins][niterations];

  for (int irange = 0; irange < mbins; irange++)
    {
      for (int iter = 0 ; iter < niterations; iter++)
	{
	  h_total_sys_range[irange][iter] = new TH1D(Form("h_total_sys_range_%d_iter_%d", irange, iter),";x_{J}; #frac{Var - Nom}{Nom}", nbins - first_bin, &dxj_bins[first_bin]);
	  h_total_sys_neg_range[irange][iter] = new TH1D(Form("h_total_sys_neg_range_%d_iter_%d", irange, iter),";x_{J}; #frac{Var - Nom}{Nom}", nbins - first_bin, &dxj_bins[first_bin]);
	  for (int ib = 0 ;ib < nbins; ib++)
	    {
	      positive_systematics[ib][irange][iter] = 0;
	      negative_systematics[ib][irange][iter] = 0;
	    }
	}
    }
  TCanvas *cxj_split = new TCanvas("cxj_split","cxj_split", 1800, 600);

  TH1D *h[5];
  TH1D *hf[5];
  TH1D *h_compare[5];

  TH1D *h_final;
  SystematicInfo fsys = sys_list["nominal"];
  
  for (auto sys_group : sys_draw)
    {
      std::cout << sys_group.first.c_str() << std::endl;
      if (strcmp(sys_group.first.c_str(), "HALF") == 0) continue;
      std::cout << sys_group.first.c_str() << std::endl;
      for (int niter = 0; niter < niterations; niter++)
	{
	  cxj_split->Clear();
      
	  dlutility::systematic_split_canvas(cxj_split, 3, 0.55);
      
	  for (int irange = 0; irange < mbins; irange++)
	    {
	      std::cout << irange << std::endl;
	      cxj_split->cd(1 + irange*2);
	      for (int isys = 0; isys < sys_group.second.size(); isys++)
		{
		  std::string sys_name = sys_group.second.at(isys);
		  SystematicInfo ssys = sys_list[sys_name.c_str()];
		  h[isys] = ssys.get_xj(irange, niter);
		  dlutility::SetLineAtt(h[isys], ssys.getColor(), ssys.getLineSize(), 1);
		  dlutility::SetMarkerAtt(h[isys], ssys.getColor(), ssys.getMarkerSize(), ssys.getMarker());
		  hf[isys] = (TH1D*) h[isys]->Rebin(nbins - first_bin, Form("h_rebin_%s_%d", ssys.get_name().c_str(), irange), &dxj_bins[first_bin]);
		}

	      h_final = fsys.get_xj(irange, niter);
	      dlutility::SetLineAtt(h_final, fsys.getColor(), fsys.getLineSize(), 1);
	      dlutility::SetMarkerAtt(h_final, fsys.getColor(), fsys.getMarkerSize(), fsys.getMarker());

	      dlutility::SetFont(h_final, 42, 0.10);

	      h_final->SetMaximum(6.0);
	      h_final->SetMinimum(0);
	      h_final->SetTitle(";x_{J}; #frac{1}{N_{pair}}#frac{dN_{pair}}{dx_{J}}");;
	      
	      TH1D *hu = (TH1D*) h_final->Rebin(nbins - first_bin, Form("h_rebin_unf_%d", irange), &dxj_bins[first_bin]);
	      hu->GetYaxis()->ChangeLabel(1, -1, 0);
	      hu->Draw("p");
	      TH1D *huc = (TH1D*) hu->Clone("huc");
	      huc->GetYaxis()->SetNdivisions(508);
	      huc->GetYaxis()->ChangeLabel(-1, -1, 0);
	      huc->GetYaxis()->ChangeLabel(1, -1, 0);
	      huc->Reset();
	      huc->SetMinimum(-0.4);
	      huc->SetMaximum(0.4);
	      huc->SetTitle(";x_{J}; #frac{Var. - Nom.}{Nom}");
	      huc->GetYaxis()->SetTitleOffset(1.1);
	      huc->GetXaxis()->SetTitleOffset(1.2);

	      dlutility::SetFont(huc, 42, 0.09);
	      huc->GetXaxis()->SetTitleSize(0.11);
	      huc->SetMaximum(0.4);
	      huc->SetMinimum(-0.4);
	      

	      for (int isys = 0; isys < sys_group.second.size(); isys++)
		{
		  std::string sys_name = sys_group.second.at(isys);
		  SystematicInfo ssys = sys_list[sys_name.c_str()];

		  hf[isys]->Draw("same p");

		  h_compare[isys] = (TH1D*) hf[isys]->Clone();
		  h_compare[isys]->Divide(hu);
		  dlutility::SetLineAtt(h_compare[isys], ssys.getColor(), 1,1);
		  dlutility::SetMarkerAtt(h_compare[isys], ssys.getColor(), 1,8);
	      
		  for (int ib = 0; ib < h_compare[isys]->GetNbinsX(); ib++)
		    {
		  
		      if (h_compare[isys]->GetBinCenter(ib+1) < 0.3) continue;
		      h_compare[isys]->SetBinContent(ib+1, h_compare[isys]->GetBinContent(ib+1) - 1);
		    }

		  h_sys[ssys.get_name().c_str()][niter][irange] = (TH1D*) h_compare[isys]->Clone();
		  h_sys[ssys.get_name().c_str()][niter][irange]->SetName(Form("h_sys_%s_%d_%d", ssys.get_name().c_str(), niter, irange));
		  h_sys_neg[ssys.get_name().c_str()][niter][irange] = (TH1D*) h_compare[isys]->Clone();
		  h_sys_neg[ssys.get_name().c_str()][niter][irange]->SetName(Form("h_sys_neg_%s_%d_%d", ssys.get_name().c_str(), niter, irange));
		  h_sys_neg[ssys.get_name().c_str()][niter][irange]->Scale(-1);
		  
		  h_compare[isys]->SetFillColorAlpha(ssys.getColor(), 0.5);
	      
		  h_compare[isys]->SetLineColor(kBlack);

		}
	      if (irange == 0)
		{
		  dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.4, 0.8 ,0, kBlack, 0.08);
		  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.4, 0.7, 0, kBlack, 0.08);
		}
	      else
		{
		  dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.1, 0.8, 0, kBlack, 0.08);
		  dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.1, 0.7, 0, kBlack, 0.08);
		}
	      
	      
	      cxj_split->cd(2 + irange*2);
	      huc->Draw();
	      for (int isys = 0; isys < sys_group.second.size(); isys++)
		h_compare[isys]->Draw("hist same");
	      
	  
	      TLine *line = new TLine(huc->GetBinLowEdge(1), 0, 1, 0);
	      line->SetLineStyle(1);
	      line->SetLineColor(kBlack);
	      line->SetLineWidth(1);
	      line->Draw("same");
	      
	    }
	  cxj_split->cd(7);
	  dlutility::DrawSPHENIXpp(0.15, 0.84, 0.08);
	  dlutility::drawText(Form("anti-k_{t} R = %0.1f", cone_size*0.1), 0.15, 0.73, 0, kBlack, 0.08);
	  dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.15, 0.65, 0, kBlack, 0.08);
	  dlutility::drawText("|#eta| #leq 0.7", 0.15, 0.57, 0, kBlack, 0.08);
	  dlutility::drawText(Form("N_{iter} = %d", niter + 1), 0.15, 0.49, 0, kBlack, 0.08);
	  TLegend *leg = new TLegend(0.15, 0.20, 0.4, 0.41);
	  leg->SetLineWidth(0);
	  leg->SetTextSize(0.08);
	  leg->SetTextFont(42);
	  leg->AddEntry(h_final, sys_titles["nominal"].c_str());
	  for (int isys = 0; isys < sys_group.second.size(); isys++)
	    {
	      std::string sys_name = sys_group.second.at(isys);
	      SystematicInfo ssys = sys_list[sys_name.c_str()];

	      leg->AddEntry(hf[isys], sys_titles[sys_name.c_str()].c_str());
	    }
	  leg->Draw("same");

	  cxj_split->SaveAs(Form("%s/systematic_plots/h_%s_xj_unfolded_pp_r%02d_range_all_iter_%d.png",  rb.get_code_location().c_str(), sys_group.first.c_str(), cone_size, niter));
	  cxj_split->SaveAs(Form("%s/systematic_plots/h_%s_xj_unfolded_pp_r%02d_range_all_iter_%d.pdf",  rb.get_code_location().c_str(), sys_group.first.c_str(), cone_size, niter));
	}
    }
	  
  
  
  for (int iter = 0 ; iter < niterations; iter++)
    {
      for (int irange = 0; irange < mbins; irange++)
	{

	  TH1D *hh = sys_list["HALF"].get_half(irange, iter);
	  if (!hh)
	    {
	      std::cout << "no half" << std::endl;
	      return;
	    }
	  h_sys["HALF"][iter][irange] = (TH1D*) hh->Clone();
	  h_sys["HALF"][iter][irange]->SetName(Form("h_sys_HALF_%d_%d", iter, irange));

	  h_sys_neg["HALF"][iter][irange] = (TH1D*) (sys_list["HALF"].get_half(irange, iter))->Clone();
	  h_sys_neg["HALF"][iter][irange]->SetName(Form("h_sys_neg_HALF_%d_%d", iter, irange));
	  h_sys_neg["HALF"][iter][irange]->Scale(-1);

	  for (auto sys_group : sys_draw)
	    {
	      bool use_group = true;
	      int reverse = 0;
	      for (int isys = 0; isys < sys_group.second.size(); isys++)
		{
		  std::string sys_name = sys_group.second.at(isys);
		  int rev = sys_list[sys_name.c_str()].get_reverse();
		  if (rev < 0)
		    {
		      use_group = false;
		    }
		  else if (rev == 1)
		    {
		      reverse = rev;
		    }
		  if (strcmp(sys_group.first.c_str(), "HALF") == 0)
		    {
		      h[isys] = sys_list[sys_name.c_str()].get_half(irange, iter);
		    }
		  else
		    {
		      h[isys] = h_sys[sys_name.c_str()][iter][irange];
		    }
		}

	      for (int i = 0; i < h[isys]->GetNbinsX(); i++)
		{
		  double min_v = 0;
		  double max_v = 0;
		  double drev = 0;
		  for (int isys = 0; isys < sys_group.second.size(); isys++)
		    {
		      double c = h[isys]->GetBinContent(i+1);
		      drev = c;
		      if (c < min_v) min_v = c;
		      if (c > max_v) max_v = c;
		    }
		  if (reverse)
		    {
		      positive_systematics[i][irange][iter] += drev*drev;
		      negative_systematics[i][irange][iter] += drev*drev;
		    }
		  else
		    {
		      positive_systematics[i][irange][iter] += max_v*max_v;
		      negative_systematics[i][irange][iter] += min_v*min_v;
		    }
		}		  
	    }
	}
    }
  for (int irange = 0; irange < mbins; irange++)
    {
      for (int iter = 0 ; iter < niterations; iter++)
	{
	  for (int i = 0; i < h_total_sys_range[irange][iter]->GetNbinsX(); i++)
	    {
	      h_total_sys_range[irange][iter]->SetBinContent(i+1, sqrt(positive_systematics[i][irange][iter]));
	      h_total_sys_neg_range[irange][iter]->SetBinContent(i+1, -1*sqrt(negative_systematics[i][irange][iter]));
	    }
	  dlutility::SetMarkerAtt(h_total_sys_neg_range[irange][iter], kBlack, 1, 8);
	  dlutility::SetMarkerAtt(h_total_sys_range[irange][iter], kBlack, 1, 8);
	}
    }

  TCanvas *cxjsys = new TCanvas("cxjsys","cxjsys", 1500, 400);
  for (int iter = 0; iter < niterations; iter++)
    {
      dlutility::systematic_split_canvas(cxjsys, 3, 0);
      TLegend *leg4 = new TLegend(0.02, 0.30, 0.9, 0.58);
      leg4->SetTextSize(0.08);
      leg4->SetTextFont(42);
      leg4->SetLineWidth(0);

      for (int irange = 0; irange < mbins; irange++)
	{
	  std::cout << irange << std::endl;
	  cxjsys->cd(irange+1);

	  h_total_sys_range[irange][iter]->SetMinimum(-0.5);
	  h_total_sys_range[irange][iter]->SetMaximum(0.5);
	  dlutility::SetFont(h_total_sys_range[irange][iter], 42, 0.08, 0.07, 0.07, 0.07);       	  
	  h_total_sys_range[irange][iter]->GetXaxis()->SetRangeUser(dxj_bins[first_bin], dxj_bins[nbins]);
	  h_total_sys_range[irange][iter]->Draw("p");
	  h_total_sys_neg_range[irange][iter]->Draw("p same");

	  for (auto sys : sys_list)
	    {
	      bool use_group = true;
	      int reverse = 0;
	      std::string sys_name = sys.first;
	      int rev = sys.second.get_reverse();
	      if (rev < 0)
		{
		  continue;
		}
	      h_sys[sys_name.c_str()][iter][irange]->Draw("same hist");
	      if (rev == 1)
		{
		  h_sys_neg[sys_name.c_str()][iter][irange]->Draw("same hist");
		}
	      if (irange == 0)
		{
		  leg4->AddEntry(h_sys[sys_name.c_str()][iter][irange], sys_titles[sys_name.c_str()].c_str());
		}
	    }

	  if (irange == 0)
	    {
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.94, 0.85 ,1, kBlack, 0.055);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.94, 0.76, 1, kBlack, 0.055);
	    }
	  else
	    {
	      dlutility::drawText(Form("%2.1f #leq #it{p}_{T,1} < %2.1f GeV ", ipt_bins[measure_bins[irange]], ipt_bins[measure_bins[irange+1]]), 0.9, 0.85, 1, kBlack, 0.06);
	      dlutility::drawText(Form("#it{p}_{T,2} #geq %2.1f GeV", ipt_bins[measure_subleading_bin]), 0.9, 0.76, 1, kBlack, 0.06);
	    }
	}
      cxjsys->cd(4);
      dlutility::DrawSPHENIXppsize(0.02, 0.9, 0.06);
      dlutility::drawText(Form("anti-#it{k}_{t} #it{R} = %0.1f", 0.1*cone_size), 0.02, 0.78, 0, kBlack, 0.06);
      dlutility::drawText(Form("|#eta| < %0.1f", 1.1 - cone_size*0.1), 0.02,0.69, 0, kBlack, 0.06);
      dlutility::drawText("#Delta#phi #geq 3#pi/4", 0.02,0.63, 0, kBlack, 0.06);
      leg4->Draw();
      cxjsys->SaveAs(Form("%s/systematic_plots/h_SYS_xj_unfolded_pp_r%02d_range_all_iter_%d.png",  rb.get_code_location().c_str(), cone_size, iter));
      cxjsys->SaveAs(Form("%s/systematic_plots/h_SYS_xj_unfolded_pp_r%02d_range_all_iter_%d.pdf",  rb.get_code_location().c_str(), cone_size, iter));
    }
  
  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     for (int irange = 0; irange < 3; irange++) */
  /* 	{ */
  /* 	  double pos_var = 0; */
  /* 	  double neg_var = 0; */
  /* 	  double nominal = h_average_xj[iter]->GetBinContent(1+irange, 1); */
  /* 	  for (auto b : *variation_average_xj[iter][irange]) */
  /* 	    { */
  /* 	      double diff = b.first - nominal; */

  /* 	      if (diff >= 0 || b.second == 1) */
  /* 		{ */
  /* 		  pos_var += diff*diff; */
  /* 		} */
  /* 	      if (diff < 0 || b.second == 1) */
  /* 		{ */
  /* 		  neg_var += diff*diff; */
  /* 		} */
  /* 	    } */
  /* 	  pos_var = sqrt(pos_var); */
  /* 	  neg_var = sqrt(neg_var); */
  /* 	  h_average_xj[iter]->SetBinContent(1 + irange, 2, pos_var); */
  /* 	  h_average_xj[iter]->SetBinContent(1 + irange, 3, neg_var); */
  /* 	} */
  /*   } */
  
  /* TFile *fsys = new TFile(Form("%s/uncertainties/systematics_pp_r%02d.root",  rb.get_code_location().c_str(), cone_size), "recreate"); */

  /* for (int iter = 0; iter < niterations; iter++) */
  /*   { */
  /*     h_average_xj[iter]->Write(); */
      
  /*     for (int irange = 0; irange < mbins; irange++) */
  /* 	{ */
  /* 	  h_sys_half[irange][iter]->Write(); */
  /* 	  h_total_sys_range[irange][iter]->Write(); */
  /* 	  h_total_sys_neg_range[irange][iter]->Write(); */
  /* 	} */
  /*   } */
  /* fsys->Close(); */

  
  return;
}

