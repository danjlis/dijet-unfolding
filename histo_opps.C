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
/* void normalize_histo(TH1D* h, const int nbins) */
/* { */
/*   Double_t bin_contents[30] = {0}; */
/*   Double_t bin_errors[30] = {0}; */

/*   Double_t integral = 0; */
/*   for (int i = 0; i < nbins; i++) */
/*     { */
/*       Double_t v = h->GetBinContent(i+1); */
/*       Double_t w = h->GetBinWidth(i+1); */
/*       Double_t err = h->GetBinError(i+1); */

/*       Double_t vw = v/w; */
/*       Double_t ew = err/w; */

/*       integral += vw; */

/*       bin_contents[i] = vw; */
/*       bin_errors[i] = vw; */

/*     } */
/*   for (int i = 0; i < nbins; i++) */
/*     { */
/*       h->SetBinContent(i+1, bin_contents[i]/integral); */
/*       h->SetBinError(i+1, bin_errors[i]/integral); */
/*     } */
/*   return; */
/* } */
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
	      h_asym_pt1pt2->SetBinError(bin, h_asym_pt1pt2->GetBinError(bin)*sqrt(2));
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

void histo_opps()

{
  std::cout << "hello world" << std::endl;
}
#endif
