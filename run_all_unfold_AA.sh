#!/usr/bin/env bash


conesize=$1
cent=$2
sys=$3
if [[ -z "$cent" || -z "$conesize" ]]; then
    echo "Usage: $0 conesize centrality"
    exit 1
fi

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
