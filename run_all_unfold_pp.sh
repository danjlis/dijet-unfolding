#!/usr/bin/env bash


conesize=$1

#root -l -q -b "unfoldDataUncertainties_noempty.cxx(\"binning.config\"${conesize})"


root -l -q -b "createResponse_noempty_pp.cxx(\"binning.config\", 0, 10, ${conesize}, 1)"

root -l -q -b "unfoldData_noempty_pp.cxx(\"binning.config\", 10, ${conesize}, 1)"

root -l -q -b "getVtxReweighting.C(${conesize}, \"binning.config\")"

root -l -q -b "createResponse_noempty_pp.cxx(\"binning.config\", 0, 10, ${conesize}, 2)"

root -l -q -b "unfoldData_noempty_pp.cxx(\"binning.config\", 10, ${conesize}, 2)"

root -l -q -b "createResponse_noempty_pp.cxx(\"binning.config\", 0, 10, ${conesize})"

root -l -q -b "unfoldData_noempty_pp.cxx(\"binning.config\", 10, ${conesize})"

root -l -q -b "unfoldDataUncertainties_noempty_pp.cxx(10, ${conesize})"


bash run_all_sys_pp.sh ${conesize} binning_negJES.config
bash run_all_sys_pp.sh ${conesize} binning_posJES.config
bash run_all_sys_pp.sh ${conesize} binning_negJER.config
bash run_all_sys_pp.sh ${conesize} binning_posJER.config
