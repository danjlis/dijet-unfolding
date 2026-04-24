#!/usr/bin/env bash
source ~/.profile

conesize=$1
fakes=1
config=$2
#./createResponse_noempty_pp -r ${conesize} -c ${config} -n 10 -p 1 -h 0 -f ${fakes} -e 0 -v 2 -m 1

./createResponse_noempty_pp -c ${config} -r ${conesize} -n 10 -p 1 -h 0 -f ${fakes} -e 0

#./unfoldData_noempty_pp -c ${config} -r ${conesize} -n 10 -p 1

#root -l -q -b "drawFullClosure.C(${conesize}, 1)"

