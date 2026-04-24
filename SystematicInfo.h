#ifndef SYSTEMATICINFO_H
#define SYSTEMATICINFO_H

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <iostream>


#include "read_binning.h"
#include "histo_opps.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
class SystematicInfo
{
 public:
  SystematicInfo(){}
  
  void set_nominal(int nominal){
    m_nominal = nominal;
  }
  void set_name(std::string name){
    m_name = name;
  }
  void set_reverse(int reverse){
    m_reverse = reverse;
  }
  
  std::string get_name()
    {
      return m_name;
    }
  
  int get_reverse()
  {
    return m_reverse;
  }
  
  int isNominal()    
  {
    return m_nominal;
  }

  void setNIterations(int niter)
  {
    m_niter = niter;
  }

  void setFile(const std::string filename)
  {
    m_filename = filename;
  }

  void setDrawingPars(int color, int color_fill, int marker, float msize, float lsize)
  {
    m_color = color;
    m_color_fill = color_fill;
    m_marker = marker;
    m_msize = msize;
    m_lsize = lsize;
  }

  int getColor() { return m_color; };
  int getFillColor() { return m_color_fill; };
  int getMarker() { return m_marker; };
  int getMarkerSize() { return m_msize; };
  int getLineSize() { return m_lsize; };

  void setup(const std::string configname);

  void gethistos();

  void gethistos_half();
  
  bool isGood()
  {
    return !m_broken;
  }
  TH1D *get_xj(int irange, int iter)
  {
    return h_final_xj_unfold_range[irange][iter];
  }

  TH1D *get_half(int irange, int iter)
  {
    return h_sys_half[irange][iter];
  }
  
  double get_average_xj(int irange, int iter)
  {
    return average_xj[iter][irange];
  }
  void setIsHalf(bool half)
  {
    m_isHalf = half;
  }
  bool isHalf()
  {
    return m_isHalf;
  }
  
 private:
  bool m_isHalf{false};
  int m_color = kBlack;
  int m_color_fill = kBlack;
  int m_marker = 1;
  float m_msize = 1;
  float m_lsize = 1;

  int m_broken = 0;
  int m_niter = 10;
  int m_nominal = 0;
  int m_reverse = 0;
  std::string m_name = "nominal";
  std::string m_filename = "";

  float m_first_xj = 0;
  int m_first_bin = 0;
  int m_nbins = 0;
  int m_nbins_pt = 0;
  int m_mbins = 0;
  int m_measure_bins[10] = {0};
  int m_subleading_measure_bins[10] = {0};


  double average_xj[20][3] = {0};
  float ipt_bins[100]={0};
  float ixj_bins[100]={0};
  double dxj_bins[100]={0};

  TH1D *h_flat_data_pt1pt2{nullptr};
  TH1D *h_flat_reco_pt1pt2{nullptr};
  TH1D *h_flat_truth_pt1pt2{nullptr};
  TH1D *h_flat_unfold_pt1pt2[20] = {nullptr};

  TH2D *h_pt1pt2_reco{nullptr};
  TH2D *h_pt1pt2_truth{nullptr};
  TH2D *h_pt1pt2_data{nullptr};
  TH2D *h_pt1pt2_unfold[20] = {nullptr};
  TH1D *h_xj_reco{nullptr};
  TH1D *h_xj_truth{nullptr};
  TH1D *h_xj_data{nullptr};
  TH1D *h_xj_unfold[20] = {nullptr};
  
  TH1D *h_xj_reco_range[3] = {nullptr};
  TH1D *h_xj_truth_range[3] = {nullptr};
  TH1D *h_xj_data_range[3] = {nullptr};
  TH1D *h_xj_unfold_range[3][20] = {nullptr};
  TH1D *h_final_xj_reco_range[3] = {nullptr};
  TH1D *h_final_xj_truth_range[3] = {nullptr};
  TH1D *h_final_xj_data_range[3] = {nullptr};
  TH1D *h_final_xj_unfold_range[3][20] = {nullptr};

  TH1D *h_sys_half[3][20] = {nullptr};
  TH1D *h_closure_test[3][20] = {nullptr};

  float m_truth_leading_cut = 0;//rb.get_truth_leading_cut();
  float m_truth_subleading_cut = 0;//rb.get_truth_subleading_cut();

  float m_reco_leading_cut = 0;//rb.get_reco_leading_cut();
  float m_reco_subleading_cut = 0;//rb.get_reco_subleading_cut();

  float m_measure_leading_cut = 0;//rb.get_measure_leading_cut();
  float m_measure_subleading_cut = 0;//rb.get_measure_subleading_cut();
    
  int m_truth_leading_bin = 0;//rb.get_truth_leading_bin();
  int m_truth_subleading_bin = 0;//rb.get_truth_subleading_bin();

  int m_reco_leading_bin = 0;//rb.get_reco_leading_bin();
  int m_reco_subleading_bin = 0;//rb.get_reco_subleading_bin();

  int m_measure_leading_bin = 0;//rb.get_measure_leading_bin();
  int m_measure_subleading_bin = 0;//rb.get_measure_subleading_bin();

};
#endif
