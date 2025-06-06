#!/usr/bin/env bash
source ./roounfold/RooUnfold-master/build/setup.sh

root -l -q -b "createResponse_full_noempty.cxx(\"binning.config\")"

root -l -q -b "createResponse_full_noempty.cxx(\"binning_negJER.config\")"

root -l -q -b "createResponse_full_noempty.cxx(\"binning_posJER.config\")"

root -l -q -b "createResponse_full_noempty.cxx(\"binning_negJES.config\")"

root -l -q -b "createResponse_full_noempty.cxx(\"binning_posJES.config\")"

root -l -q -b "createResponse_full_noempty.cxx(\"binning_njet.config\")"

root -l -q -b "createResponse_full_noempty.cxx(\"binning_VTX.config\")"

root -l -q -b "unfoldData_noempty.cxx(\"binning.config\")"

root -l -q -b "unfoldData_noempty.cxx(\"binning_posJER.config\")"

root -l -q -b "unfoldData_noempty.cxx(\"binning_negJER.config\")"

root -l -q -b "unfoldData_noempty.cxx(\"binning_negJES.config\")"

root -l -q -b "unfoldData_noempty.cxx(\"binning_posJES.config\")"

root -l -q -b "unfoldData_noempty.cxx(\"binning_njet.config\")"

root -l -q -b "unfoldData_noempty.cxx(\"binning_vtx.config\")"

