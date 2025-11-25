#!/usr/bin/env bash

conesize=$1

#root -l -q -b "unfoldDataUncertainties_noempty.cxx(\"binning.config\"${conesize})"

root -l -q -b "createResponse_noempty.cxx(\"binning_primer.config\", 0, 10, ${conesize})"

root -l -q -b "unfoldData_noempty.cxx(\"binning_primer.config\", 10, ${conesize})"

root -l -q -b "getVtxReweighting.C(${conesize})"

root -l -q -b "createResponse_noempty.cxx(\"binning_primer2.config\", 0, 10, ${conesize})"

root -l -q -b "unfoldData_noempty.cxx(\"binning_primer2.config\", 10, ${conesize})"

# root -l -q -b "createResponse_noempty.cxx(\"binning.config\", 1, 10, ${conesize})"

# root -l -q -b "unfoldHalf_noempty.cxx(\"binning.config\", 10, ${conesize})"

root -l -q -b "createResponse_noempty.cxx(\"binning.config\", 0, 10, ${conesize})"

root -l -q -b "unfoldData_noempty.cxx(\"binning.config\", 10, ${conesize})"

root -l -q -b "unfoldDataUncertainties_noempty.cxx(10, ${conesize})"

root -l -q -b "createResponse_noempty.cxx(\"binning_prior.config\", 0, 10, ${conesize})"

root -l -q -b "unfoldData_noempty.cxx(\"binning_prior.config\", 10, ${conesize})"

root -l -q -b "createResponse_noempty.cxx(\"binning_negJER.config\", 0, 10, ${conesize})"

root -l -q -b "createResponse_noempty.cxx(\"binning_posJER.config\", 0, 10, ${conesize})"

root -l -q -b "createResponse_noempty.cxx(\"binning_negJES.config\", 0, 10, ${conesize})"

root -l -q -b "createResponse_noempty.cxx(\"binning_posJES.config\", 0, 10, ${conesize})"

root -l -q -b "createResponse_noempty.cxx(\"binning_njet.config\", 0, 10, ${conesize})"

root -l -q -b "unfoldData_noempty.cxx(\"binning_negJER.config\", 10, ${conesize})"

root -l -q -b "unfoldData_noempty.cxx(\"binning_posJER.config\", 10, ${conesize})"

root -l -q -b "unfoldData_noempty.cxx(\"binning_negJES.config\", 10, ${conesize})"

root -l -q -b "unfoldData_noempty.cxx(\"binning_posJES.config\", 10, ${conesize})"

root -l -q -b "unfoldData_noempty.cxx(\"binning_njet.config\", 10, ${conesize})"
