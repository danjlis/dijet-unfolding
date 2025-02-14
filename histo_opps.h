#ifndef HISTO_OPPS_H
#define HISTO_OPPS_H

#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TGraph.h"
#include "TLatex.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TColor.h"
namespace histo_opps
{

  void get_xj_systematics(TH1D *h1, TH1D *hs, const int nbins)
  {
    for (int i = 0; i < nbins; i++)
      {
	h1->SetBinError(i+1, hs->GetBinContent(i+1)*h1->GetBinContent(i+1));
      }
  }
  void set_xj_errors(TH1D *h1, TProfile *hp, const int nbins)
  {
    for (int i = 0; i < nbins; i++)
      {
	float old_error = h1->GetBinError(i+1);
	float stat_error = hp->GetBinError(i+1);
	float new_error = sqrt(TMath::Power(old_error, 2) + TMath::Power(stat_error, 2));
	h1->SetBinError(i+1, new_error);
      }
  }
  void finalize_xj(TH1D *h1, TH1D *h2, const int nbins, float first_xj)
  {
    for (int i = 0; i < nbins; i++)
      {
	if (h1->GetBinLowEdge(i+1) >= first_xj)
	  {
	    h2->SetBinContent(i+1, h1->GetBinContent(i+1));
	    h2->SetBinError(i+1, h1->GetBinError(i+1));
	  }
      }
  }

  void finalize_xj(TProfile *h1, TH1D *h2, const int nbins, float first_xj)
  {
    for (int i = 0; i < nbins; i++)
      {
	if (h1->GetBinLowEdge(i+1) >= first_xj)
	  {
	    h2->SetBinContent(i+1, h1->GetBinContent(i+1));
	    h2->SetBinError(i+1, h1->GetBinError(i+1));
	  }
      }
  }

  void normalize_histo(TH1D* h, const int nbins)
  {
    Double_t bin_contents[30] = {0};
    Double_t bin_errors[30] = {0};

    Double_t integral = 0;

    for (int i = 0; i < nbins; i++)
      {

        Double_t v = h->GetBinContent(i+1);
        Double_t w = h->GetBinWidth(i+1);
        Double_t err = h->GetBinError(i+1);

        Double_t vw = v/w;
        Double_t ew = err/w;

        integral += v;

        bin_contents[i] = v/w;
        bin_errors[i] = ew;

      }

    for (int i = 0; i < nbins; i++)
      {
        h->SetBinContent(i+1, bin_contents[i]/integral);
        h->SetBinError(i+1, bin_errors[i]/integral);
      }
    return;
  }

  void tprofile_to_histo(TProfile *hp, TH1D* h, const int nbins)
  {
    for (int ib = 0; ib < nbins; ib++)
      {
	h->SetBinContent(ib+1, hp->GetBinContent(ib+1));
	h->SetBinError(ib+1, hp->GetBinError(ib+1));
      }
  }

  void project_xj(TH2D* hpt1pt2, TH1D* h_xj, const int nbins, const int start_leading_bin, const int end_leading_bin, const int start_subleading_bin, const int end_subleading_bin)
  {

    TH1D *h_unc = (TH1D*) h_xj->Clone();
    h_unc->Reset();

    TH2D *h_asym_pt1pt2 = (TH2D*) hpt1pt2->Clone();
  
    for (int ix = 0; ix < nbins; ix++)
      {
	for (int iy = 0; iy < nbins; iy++)
	  {
	    int bin = h_asym_pt1pt2->GetBin(ix+1, iy+1);

	    if (ix > iy)
	      {
		h_asym_pt1pt2->SetBinContent(bin, h_asym_pt1pt2->GetBinContent(bin)*2.);
		h_asym_pt1pt2->SetBinError(bin, h_asym_pt1pt2->GetBinError(bin)*2);
	      }
	    else if (ix < iy)
	      {
		h_asym_pt1pt2->SetBinContent(bin, 0);
		h_asym_pt1pt2->SetBinError(bin, 0);
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
	    //std::cout << ix << " / " << iy << " -- > " << xjbin_low << "--"<<xjbin_high<<std::endl;
	    int bin = h_asym_pt1pt2->GetBin(ix+1, iy+1);

	    if (ix < start_leading_bin) continue;
	    if (iy < start_subleading_bin) continue;
	    if (ix >= end_leading_bin) continue;
	    if (iy >= end_subleading_bin) continue;

	    if (ix == iy)
	      {

		h_xj->Fill(h_xj->GetBinCenter(xjbin_low), h_asym_pt1pt2->GetBinContent(bin));
		h_unc->Fill(h_xj->GetBinCenter(xjbin_low),TMath::Power( h_asym_pt1pt2->GetBinError(bin), 2));
	      }
	    else
	      {
		h_xj->Fill(h_xj->GetBinCenter(xjbin_high), h_asym_pt1pt2->GetBinContent(bin)/2.);
		h_xj->Fill(h_xj->GetBinCenter(xjbin_low), h_asym_pt1pt2->GetBinContent(bin)/2.);
		h_unc->Fill(h_xj->GetBinCenter(xjbin_high),TMath::Power( h_asym_pt1pt2->GetBinError(bin)/sqrt(2), 2));
		h_unc->Fill(h_xj->GetBinCenter(xjbin_low),TMath::Power( h_asym_pt1pt2->GetBinError(bin)/sqrt(2), 2));
	      }
	  }
      }

    for (int ix = 0; ix < nbins; ix++)
      {
	h_xj->SetBinError(ix+1, sqrt(h_unc->GetBinContent(ix+1)));
      }

    return;
  }

  void make_sym_pt1pt2(TH1D *hflat, TH2D* hpt1pt2, const int nbins)
  {

    for (int ib = 0; ib < nbins*nbins; ib++)
      {
	int xbin = ib/nbins;
	int ybin = ib%nbins;
      
	int b = hpt1pt2->GetBin(xbin+1, ybin+1);

	hpt1pt2->SetBinContent(b, hflat->GetBinContent(ib+1));

	hpt1pt2->SetBinError(b, hflat->GetBinError(ib+1));
      }
    return;
  }

  void skim_down_histo(TH1D *h_skim, TH1D *h_full, TH1D *h_mapping)
  {

    int nbins = h_mapping->GetXaxis()->GetNbins();

    for (int ib = 0; ib < nbins*nbins; ib++)
      {
	if (!(h_mapping->GetBinContent(ib+1) == 0))
	  {
	    h_skim->SetBinContent(h_mapping->GetBinContent(ib+1), h_full->GetBinContent(ib+1));
	    h_skim->SetBinError(h_mapping->GetBinContent(ib+1), h_full->GetBinError(ib+1));
	  }
      }
  }

  void fill_up_histo(TH1D *h_skim, TH1D *h_full, TH1D *h_mapping)
  {
    
    int nbins = h_mapping->GetXaxis()->GetNbins();
 
    for (int ib = 0; ib < nbins*nbins; ib++)
	{
	  int bin = h_mapping->GetBinContent(ib+1);
	  if (bin)
	    {
	      h_full->SetBinContent(ib+1, h_skim->GetBinContent(bin));
	      h_full->SetBinError(ib+1, h_skim->GetBinError(bin));
	    }
	}
  }
};

#endif
