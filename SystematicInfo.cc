#include "SystematicInfo.h"

void SystematicInfo::gethistos_half()
{
  std::cout << "GETTING HALF" << std::endl;
  TFile *fin = new TFile(m_filename.c_str(),"r");
  if (!fin || fin->IsZombie()) {
    m_broken = 1;
    std::cout << "no file " << std::endl;
    return;
  }
  for (int irange = 0; irange < m_mbins; irange++)
    {
      for (int iter = 0; iter < m_niter; iter++)
	{

	  h_closure_test[irange][iter] = (TH1D*) fin->Get(Form("h_closure_test_%d_%d", irange, iter));
	  if (!h_closure_test[irange][iter]) {
	    std::cout << "no hist " << std::endl;
	    m_broken = 1;
	    return;
	  }
	  h_closure_test[irange][iter]->SetDirectory(0);
	  h_closure_test[irange][iter]->SetName(Form("h_closure_test_range_%s_%d_%d", m_name.c_str(), irange, iter));

	  h_sys_half[irange][iter] = (TH1D*) h_closure_test[irange][iter]->Rebin(m_nbins - m_first_bin, Form("h_rebin_half_%s_%d_%d", m_name.c_str(), irange, iter), &dxj_bins[m_first_bin]);
	  h_sys_half[irange][iter]->SetDirectory(0);
	}
    }
  fin->Close();
  delete fin;
  
}
void SystematicInfo::setup(const std::string configname)
{

  read_binning rb(configname.c_str());

  m_first_xj = rb.get_first_xj();
  m_first_bin = 0;

  m_nbins = rb.get_nbins();
  m_nbins_pt = m_nbins+1;
  rb.get_pt_bins(ipt_bins);
  rb.get_xj_bins(ixj_bins);

  for (int i = 0 ; i < m_nbins + 1; i++)
    {
      dxj_bins[i] = ixj_bins[i];
      if (ixj_bins[i] > 0.3 && m_first_bin == 0) m_first_bin = i;
    }

  ipt_bins[m_nbins_pt] = 100;

  m_truth_leading_cut = rb.get_truth_leading_cut();
  m_truth_subleading_cut = rb.get_truth_subleading_cut();

  m_reco_leading_cut = rb.get_reco_leading_cut();
  m_reco_subleading_cut = rb.get_reco_subleading_cut();

  m_measure_leading_cut = rb.get_measure_leading_cut();
  m_measure_subleading_cut = rb.get_measure_subleading_cut();
    
  m_truth_leading_bin = rb.get_truth_leading_bin();
  m_truth_subleading_bin = rb.get_truth_subleading_bin();

  m_reco_leading_bin = rb.get_reco_leading_bin();
  m_reco_subleading_bin = rb.get_reco_subleading_bin();

  m_measure_leading_bin = rb.get_measure_leading_bin();
  m_measure_subleading_bin = rb.get_measure_subleading_bin();

  m_mbins = rb.get_measure_bins();

  for (int ir = 0; ir < m_mbins+1; ir++)
    {

      m_measure_bins[ir] = rb.get_measure_region(ir);

      m_subleading_measure_bins[ir] = rb.get_subleading_measure_region(ir);
    }
}
void SystematicInfo::gethistos()
{
  m_broken = 0;
  if (m_isHalf)
    {
      gethistos_half();
      return;
    }
  
  std::cout << "getting histos from : " << m_filename << std::endl;
  TFile *fin = new TFile(m_filename.c_str(),"r");

  if (!fin || fin->IsZombie()) {
    m_broken = 1;
    return;
  }
  h_flat_data_pt1pt2 = (TH1D*) fin->Get("h_data_flat_pt1pt2");    

  if (!h_flat_data_pt1pt2) {
    m_broken = 1;
    return;
  }
  h_flat_data_pt1pt2->SetDirectory(0);
  h_flat_data_pt1pt2->SetName(Form("h_flat_data_pt1pt2_%s", m_name.c_str()));

  h_flat_reco_pt1pt2 = (TH1D*) fin->Get("h_reco_flat_pt1pt2");
  if (!h_flat_reco_pt1pt2) {
    m_broken = 1;
    return;
  }
  h_flat_reco_pt1pt2->SetDirectory(0);
  h_flat_reco_pt1pt2->SetName(Form("h_flat_reco_pt1pt2_%s", m_name.c_str()));

  h_flat_truth_pt1pt2 = (TH1D*) fin->Get("h_truth_flat_pt1pt2");
  if (!h_flat_truth_pt1pt2) {
    m_broken = 1;
    return;
  }
  h_flat_truth_pt1pt2->SetDirectory(0);
  h_flat_truth_pt1pt2->SetName(Form("h_flat_truth_pt1pt2_%s", m_name.c_str()));
  for (int iter = 0; iter < m_niter; iter++)
    {
      h_flat_unfold_pt1pt2[iter] = (TH1D*) fin->Get(Form("h_flat_unfold_pt1pt2_%d", iter));
      if (!h_flat_unfold_pt1pt2[iter])
	{
	  m_broken = 1;
	  return;
	}
      h_flat_unfold_pt1pt2[iter]->SetDirectory(0);
      h_flat_unfold_pt1pt2[iter]->SetName(Form("h_flat_unfold_pt1pt2_%d_%s", iter, m_name.c_str()));
    }

  fin->Close();
  delete fin;

  h_pt1pt2_reco = new TH2D(Form("h_pt1pt2_reco_%s", m_name.c_str()), ";#it{p}_{T,1};#it{p}_{T,2}", m_nbins_pt, ipt_bins, m_nbins_pt, ipt_bins);
  h_pt1pt2_truth = new TH2D(Form("h_pt1pt2_truth_%s", m_name.c_str()), ";#it{p}_{T,1};#it{p}_{T,2}", m_nbins_pt, ipt_bins, m_nbins_pt, ipt_bins);
  h_pt1pt2_data = new TH2D(Form("h_pt1pt2_data_%s", m_name.c_str()), ";#it{p}_{T,1};#it{p}_{T,2}", m_nbins_pt, ipt_bins, m_nbins_pt, ipt_bins);
  
  for (int iter = 0; iter < m_niter; iter++)
    {
      h_pt1pt2_unfold[iter] = new TH2D(Form("h_pt1pt2_unfold_iter%d_%s", iter, m_name.c_str()), ";#it{p}_{T,1};#it{p}_{T,2}",m_nbins_pt, ipt_bins, m_nbins_pt, ipt_bins);
    }

  h_xj_reco = new TH1D(Form("h_xj_reco_%s", m_name.c_str()), ";x_{J};", m_nbins, ixj_bins);
  h_xj_truth = new TH1D(Form("h_xj_truth_%s", m_name.c_str()), ";x_{J};", m_nbins, ixj_bins);
  h_xj_data = new TH1D(Form("h_xj_data_%s", m_name.c_str()), ";x_{J};", m_nbins, ixj_bins);
  for (int iter = 0; iter < m_niter; iter++)
    {
      h_xj_unfold[iter] = new TH1D(Form("h_xj_unfold_iter%d_%s", iter, m_name.c_str()), ";x_{J};",m_nbins, ixj_bins);
    }
  
  for (int irange = 0; irange < m_mbins; irange++)
    {
      h_xj_reco_range[irange] = new TH1D(Form("h_xj_reco_range_%d_%s", irange, m_name.c_str()), ";x_{J};", m_nbins, ixj_bins);
      h_xj_truth_range[irange] = new TH1D(Form("h_xj_truth_range_%d_%s", irange, m_name.c_str()), ";x_{J};", m_nbins, ixj_bins);
      h_xj_data_range[irange] = new TH1D(Form("h_xj_data_range_%d_%s", irange, m_name.c_str()), ";x_{J};", m_nbins, ixj_bins);
      for (int iter = 0; iter < m_niter; iter++)
	{
	  h_xj_unfold_range[irange][iter] = new TH1D(Form("h_xj_unfold_%d_iter%d_%s", irange, iter, m_name.c_str()), ";x_{J};",m_nbins, ixj_bins);
	}
    }

    
  histo_opps::make_sym_pt1pt2(h_flat_reco_pt1pt2, h_pt1pt2_reco, m_nbins_pt);
  histo_opps::make_sym_pt1pt2(h_flat_truth_pt1pt2, h_pt1pt2_truth, m_nbins_pt);
  histo_opps::make_sym_pt1pt2(h_flat_data_pt1pt2, h_pt1pt2_data, m_nbins_pt);
  for (int iter = 0; iter < m_niter; iter++)
    {
      histo_opps::make_sym_pt1pt2(h_flat_unfold_pt1pt2[iter], h_pt1pt2_unfold[iter], m_nbins_pt);
    }

  histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco, m_nbins_pt, m_measure_leading_bin, m_nbins - 2, m_measure_subleading_bin, m_nbins - 2);
  histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth, m_nbins_pt, m_measure_leading_bin, m_nbins - 2, m_measure_subleading_bin, m_nbins - 2);
  histo_opps::project_xj(h_pt1pt2_data, h_xj_data, m_nbins_pt, m_measure_leading_bin, m_nbins - 2, m_measure_subleading_bin, m_nbins - 2);
  for (int iter = 0; iter < m_niter; iter++)
    {
      histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold[iter], m_nbins_pt, m_measure_leading_bin, m_nbins - 2, m_measure_subleading_bin, m_nbins - 2);
    }
  
  histo_opps::normalize_histo(h_xj_reco, m_nbins);
  histo_opps::normalize_histo(h_xj_truth, m_nbins);
  histo_opps::normalize_histo(h_xj_data, m_nbins);
  for (int iter = 0; iter < m_niter; iter++)
    {
      histo_opps::normalize_histo(h_xj_unfold[iter], m_nbins);
    }
    
  for (int irange = 0; irange < m_mbins; irange++)
    {
      h_final_xj_reco_range[irange] = new TH1D(Form("h_final_xj_reco_range_%d_%s", irange, m_name.c_str()), ";x_{J};", m_nbins, ixj_bins);
      h_final_xj_truth_range[irange] = new TH1D(Form("h_final_xj_truth_range_%d_%s", irange, m_name.c_str()), ";x_{J};", m_nbins, ixj_bins);
      h_final_xj_data_range[irange] = new TH1D(Form("h_final_xj_data_range_%d_%s", irange, m_name.c_str()), ";x_{J};", m_nbins, ixj_bins);
      for (int iter = 0; iter < m_niter; iter++)
	{
	  h_final_xj_unfold_range[irange][iter] = new TH1D(Form("h_final_xj_unfold_%d_iter%d_%s", irange, iter, m_name.c_str()), ";x_{J};", m_nbins, ixj_bins);
	}
	

      histo_opps::project_xj(h_pt1pt2_reco, h_xj_reco_range[irange], m_nbins_pt, m_measure_bins[irange], m_measure_bins[irange+1], m_measure_subleading_bin, m_nbins - 2);
      histo_opps::project_xj(h_pt1pt2_truth, h_xj_truth_range[irange], m_nbins_pt, m_measure_bins[irange], m_measure_bins[irange+1], m_measure_subleading_bin, m_nbins - 2);
      histo_opps::project_xj(h_pt1pt2_data, h_xj_data_range[irange], m_nbins_pt, m_measure_bins[irange], m_measure_bins[irange+1], m_measure_subleading_bin, m_nbins - 2);
      for (int iter = 0; iter < m_niter; iter++)
	{
	  histo_opps::project_xj(h_pt1pt2_unfold[iter], h_xj_unfold_range[irange][iter], m_nbins_pt, m_measure_bins[irange], m_measure_bins[irange+1], m_measure_subleading_bin, m_nbins - 2);
	}
	
	
      histo_opps::finalize_xj(h_xj_reco_range[irange], h_final_xj_reco_range[irange], m_nbins, m_first_xj);
      histo_opps::finalize_xj(h_xj_truth_range[irange], h_final_xj_truth_range[irange], m_nbins, m_first_xj);
      histo_opps::finalize_xj(h_xj_data_range[irange], h_final_xj_data_range[irange], m_nbins, m_first_xj);
      for (int iter = 0; iter < m_niter; iter++)
	{
	  histo_opps::finalize_xj(h_xj_unfold_range[irange][iter], h_final_xj_unfold_range[irange][iter], m_nbins, m_first_xj);
	}
      histo_opps::normalize_histo(h_final_xj_reco_range[irange], m_nbins);
      histo_opps::normalize_histo(h_final_xj_truth_range[irange], m_nbins);
      histo_opps::normalize_histo(h_final_xj_data_range[irange], m_nbins);
      for (int iter = 0; iter < m_niter; iter++)
	{
	  histo_opps::normalize_histo(h_final_xj_unfold_range[irange][iter], m_nbins);
	  average_xj[iter][irange] = histo_opps::get_average_xj(h_final_xj_unfold_range[irange][iter]);
	}
    }    
}
