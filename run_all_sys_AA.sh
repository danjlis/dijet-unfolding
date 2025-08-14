#!/usr/bin/env bash


conesize=$1
cent=$2
sysconfig=$3
if [[ -z "$cent" || -z "$conesize" || -z "$sysconfig" ]]; then
    echo "Usage: $0 conesize centrality"
    exit 1
fi

root -l -q -b "getBackground.C(${conesize}, ${cent}, \"binning_AA.config\")"

root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent}, 1)"

root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent}, 1)"

root -l -q -b "getCentralityReweighting.C(${conesize}, ${cent}, \"${sysconfig}\")"

root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent}, 2)"

root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent}, 2)"

# root -l -q -b "createResponse_noempty_AA.cxx(\"binning.config\", 1, 10, ${conesize}, 0)"
# root -l -q -b "unfoldHalf_noempty_AA.cxx(\"binning.config\", 10, ${conesize}, 2)"

root -l -q -b "createResponse_noempty_AA.cxx(\"${sysconfig}\", 0, 10, ${conesize}, ${cent})"

root -l -q -b "unfoldData_noempty_AA.cxx(\"${sysconfig}\", 10, ${conesize}, ${cent})"


