#!/usr/bin/env bash

conesize=3

for cent in $(seq 0 3); do

    #root -l -q -b "unfoldDataUncertainties_noempty.cxx(\"binning.config\"${conesize})"

    root -l -q -b "getBackground.C(${conesize}, ${cent}, \"binning_AA.config\")"
    
    root -l -q -b "createResponse_noempty_AA.cxx(\"binning_AA.config\", 0, 10, ${conesize}, ${cent}, 1)"
    
    root -l -q -b "unfoldData_noempty_AA.cxx(\"binning_AA.config\", 10, ${conesize}, ${cent}, 1)"
    
    root -l -q -b "getCentralityReweighting.C(${conesize}, ${cent}, \"binning_AA.config\")"
    
    root -l -q -b "createResponse_noempty_AA.cxx(\"binning_AA.config\", 0, 10, ${conesize}, ${cent}, 2)"
    
    root -l -q -b "unfoldData_noempty_AA.cxx(\"binning_AA.config\", 10, ${conesize}, ${cent}, 2)"
    
    # root -l -q -b "createResponse_noempty_AA.cxx(\"binning.config\", 1, 10, ${conesize}, 0)"
    # root -l -q -b "unfoldHalf_noempty_AA.cxx(\"binning.config\", 10, ${conesize}, 2)"
    
    root -l -q -b "createResponse_noempty_AA.cxx(\"binning_AA.config\", 0, 10, ${conesize}, ${cent})"
    
    root -l -q -b "unfoldData_noempty_AA.cxx(\"binning_AA.config\", 10, ${conesize}, ${cent})"
    
    root -l -q -b "unfoldDataUncertainties_noempty_AA.cxx(10, ${conesize}, ${cent})"


    bash run_all_sys_AA.sh 3 ${cent} binning_negJES_AA.config
    bash run_all_sys_AA.sh 3 ${cent} binning_posJES_AA.config
    bash run_all_sys_AA.sh 3 ${cent} binning_negJER_AA.config
    bash run_all_sys_AA.sh 3 ${cent} binning_posJER_AA.config

    bash run_all_sys_AA.sh 3 ${cent} binning_ZYAM_AA.config
    bash run_all_sys_AA.sh 3 ${cent} binning_prior_AA.config

    root -l -b -q "drawSys_AA.C(3, ${cent})"
    root -l -b -q "drawFinalUnfold_AA.C(3, ${cent})"

done
