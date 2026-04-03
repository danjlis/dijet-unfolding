#!/usr/bin/env bash

source ~/.profile

conesize=$1
sysconfig=$2
usefakes=$3
final=$4

#root -l -q -b "unfoldDataUncertainties_noempty.cxx(\"binning.config\"${conesize})"
if [[ $final -eq 0 ]]; then
    ./createResponse_noempty_pp -c ${sysconfig} -r ${conesize} -n 10 -p 1 -h 0 -f ${usefakes} -e 0 -m 1

    ./createResponse_noempty_pp -c ${sysconfig} -r ${conesize} -n 10 -p 1 -h 0 -f ${usefakes} -e 0
    
    ./unfoldData_noempty_pp -c ${sysconfig} -r ${conesize} -n 10 -p 1
    
    root -l -q -b "getVtxReweighting.C(${conesize}, \"${sysconfig}\")"
    
    ./createResponse_noempty_pp -c ${sysconfig} -r ${conesize} -n 10 -p 2 -h 0 -f ${usefakes} -e 0
    ./unfoldData_noempty_pp -c ${sysconfig} -r ${conesize} -n 10 -p 2
fi
./createResponse_noempty_pp -c ${sysconfig} -r ${conesize} -n 10 -p 0 -h 0 -f ${usefakes} -e 0
./unfoldData_noempty_pp -c ${sysconfig} -r ${conesize} -n 10 -p 0

