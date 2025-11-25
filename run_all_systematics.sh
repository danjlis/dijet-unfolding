#!/usr/bin/env bash

for i in $(seq 0 3); do
    echo "cent bin ${i}"
    bash run_all_sys_AA.sh 3 ${i} binning_negJES_AA.config
    bash run_all_sys_AA.sh 3 ${i} binning_posJES_AA.config
    bash run_all_sys_AA.sh 3 ${i} binning_negJER_AA.config
    bash run_all_sys_AA.sh 3 ${i} binning_posJER_AA.config

    bash run_all_sys_AA.sh 3 ${i} binning_ZYAM_AA.config
    bash run_all_sys_AA.sh 3 ${i} binning_prior_AA.config
    root -l -b -q "drawSys_AA(3, ${i})"
done

