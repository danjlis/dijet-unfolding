#!/usr/bin/env bash


conesize=$1
sysconfig=$2

root -l -q -b "createResponse_noempty_pp.cxx(\"${sysconfig}\", 0, 10, ${conesize}, 1)"

root -l -q -b "unfoldData_noempty_pp.cxx(\"${sysconfig}\", 10, ${conesize}, 1)"

#root -l -q -b "getCentralityReweighting.C(${conesize}, \"${sysconfig}\")"

root -l -q -b "createResponse_noempty_pp.cxx(\"${sysconfig}\", 0, 10, ${conesize}, 2)"

root -l -q -b "unfoldData_noempty_pp.cxx(\"${sysconfig}\", 10, ${conesize}, 2)"

root -l -q -b "createResponse_noempty_pp.cxx(\"${sysconfig}\", 0, 10, ${conesize})"

root -l -q -b "unfoldData_noempty_pp.cxx(\"${sysconfig}\", 10, ${conesize})"


