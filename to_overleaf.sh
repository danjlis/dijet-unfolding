#! /usr/bin/env bash

conesize=$1

if [[ -d to_overleaf/r0${conesize} ]];
then
    rm -rf to_overleaf/r0${conesize}/*
else
    mkdir -p to_overleaf/r0${conesize}/
fi

mkdir -p to_overleaf/r0${conesize}/final_plots/
mkdir -p to_overleaf/r0${conesize}/unfolding_plots/
mkdir -p to_overleaf/r0${conesize}/sys_plots/
mkdir -p to_overleaf/r0${conesize}/truth_plots/
mkdir -p to_overleaf/r0${conesize}/jer_plots/
mkdir -p to_overleaf/r0${conesize}/dphi_plots/

unfolding_dir=to_overleaf/r0${conesize}/unfolding_plots/
sys_dir=to_overleaf/r0${conesize}/sys_plots/
final_dir=to_overleaf/r0${conesize}/final_plots/
truth_dir=to_overleaf/r0${conesize}/truth_plots/
jer_dir=to_overleaf/r0${conesize}/jer_plots/
dphi_dir=to_overleaf/r0${conesize}/dphi_plots/

cp unfolding_plots/combined_sample_pp_r0${conesize}_nominal.pdf ${truth_dir}

cp unfolding_plots/etaphi_comp_pp_r0${conesize}_PRIMER1_nominal.pdf ${unfolding_dir}
cp unfolding_plots/spectra_comp_pp_r0${conesize}_PRIMER1_nominal.pdf ${unfolding_dir}
cp unfolding_plots/datasim_eta_lead_pp_r0${conesize}_nominal.pdf ${unfolding_dir}
cp unfolding_plots/datasim_mbd_pp_r0${conesize}_nominal.pdf ${unfolding_dir}
cp unfolding_plots/datasim_emfrac_pp_r0${conesize}_EMFRAC.pdf ${unfolding_dir}
cp unfolding_plots/h_xj_half_closure_pp_r0${conesize}_range_*_iter_all_nominal.pdf ${unfolding_dir}
cp unfolding_plots/h_xj_half_closure_pp_r0${conesize}_range_*_iter_xjbin_nominal.pdf ${unfolding_dir}
cp unfolding_plots/iteration_tune_pp_r0${conesize}_nominal.pdf ${unfolding_dir}
cp unfolding_plots/iteration_tune_pp_r0${conesize}_PRIMER2_nominal.pdf ${unfolding_dir}
cp unfolding_plots/prior_matrix_pp_r0${conesize}_nominal.pdf ${unfolding_dir}
cp unfolding_plots/response_matrix_pp_r0${conesize}_nominal.pdf ${unfolding_dir}
cp unfolding_plots/h_xj_full_closure_r0${conesize}_range_*_iter_5_nominal.pdf ${unfolding_dir}

cp final_plots/h_final_xj_unfolded_pp_r0${conesize}_range_*.pdf ${final_dir}
cp final_plots/h_xj_unfolded_pp_r0${conesize}_range_*.pdf ${final_dir}

cp systematic_plots//h_JER_xj_unfolded_pp_r0${conesize}_range_all_iter_1.pdf ${sys_dir}
cp systematic_plots//h_JES_xj_unfolded_pp_r0${conesize}_range_all_iter_1.pdf ${sys_dir}
cp systematic_plots//h_PRIOR_xj_unfolded_pp_r0${conesize}_range_all_iter_1.pdf ${sys_dir}
cp systematic_plots//h_TRIGGER_xj_unfolded_pp_r0${conesize}_range_all_iter_1.pdf ${sys_dir}
cp systematic_plots//h_PHI_xj_unfolded_pp_r0${conesize}_range_all_iter_1.pdf ${sys_dir}
cp systematic_plots//h_CA_xj_unfolded_pp_r0${conesize}_range_all_iter_1.pdf ${sys_dir}
cp systematic_plots//h_EMFRAC_xj_unfolded_pp_r0${conesize}_range_all_iter_1.pdf ${sys_dir}
cp systematic_plots//h_SYS_xj_unfolded_pp_r0${conesize}_range_all_iter_1.pdf ${sys_dir}

cp ../../JESR/r0${conesize}/jer/datasim_jer_r0${conesize}_im_v_bi_3jet3_0_closure_PYTHIA.pdf  ${jer_dir}
cp ../../JESR/r0${conesize}/jer/sigp_r0${conesize}_3jet3_pta4_0_closure_PYTHIA.pdf  ${jer_dir}
cp ../../JESR/r0${conesize}/jer/sige_r0${conesize}_3jet3_pta4_0_closure_PYTHIA.pdf  ${jer_dir}
cp ../../JESR/r0${conesize}/jer/aj_r0${conesize}_3jet3_pta4_0_closure_PYTHIA.pdf   ${jer_dir}

cp unfolding_plots/h_xj_half_closure_pp_r0${conesize}_range_0_iter_all_HALF_nominal.pdf ${unfolding_dir}
cp unfolding_plots/h_xj_half_closure_pp_r0${conesize}_range_1_iter_all_HALF_nominal.pdf ${unfolding_dir}
cp unfolding_plots/h_xj_half_closure_pp_r0${conesize}_range_2_iter_all_HALF_nominal.pdf ${unfolding_dir}

cp ../../JESR/r0${conesize}/jer/jer_fits_r0${conesize}_0_closure_PYTHIA.pdf ${jer_dir}
cp ../../JESR/r0${conesize}/jer/sigpsigv_data_r0${conesize}_0_closure_PYTHIA.pdf ${jer_dir}
cp unfolding_plots/response_matrix_pp_r0${conesize}_nominal.pdf $unfolding_dir
cp unfolding_plots/response_matrix_zoom_pp_r0${conesize}_nominal.pdf $unfolding_dir

cp dphi_plots/dphi_data_pythia_pp_r0${conesize}_range_*.pdf $dphi_dir
cp dphi_plots/h_final_data_pythia_dphi_pp_r0${conesize}_range_*.pdf $dphi_dir
cp dphi_plots/h_final_data_pythia_dphi_pp_r0${conesize}_rangemin_*.pdf $dphi_dir
cp dphi_plots/correction_pp_r0${conesize}_rangemin_*.pdf ${dphi_dir}
cp dphi_plots/correction_pp_r0${conesize}_range_*.pdf ${dphi_dir}
cp dphi_plots/h_sys_all_pp_r0${conesize}_range_*.pdf ${dphi_dir}
cp dphi_plots/h_sys_all_pp_r0${conesize}_rangemin_*.pdf ${dphi_dir}
