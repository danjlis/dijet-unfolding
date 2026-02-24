#!/usr/bin/env bash
source ~/.profile

conesize=$1
fakes=$2
onlynom=$3
interactive=$4

interaction () {

    if [[ $1 -gt 0 ]]; then
	read -n 1 -s -r -p "Press any key to continue..."
    fi
}

echo "First PRIMER"

echo " FULL" 
./createResponse_noempty_pp -c binning.config -r ${conesize} -n 10 -p 1 -h 0 -f ${fakes} -e 0

interaction $interactive

echo " DATA"
./unfoldData_noempty_pp -c binning.config -r ${conesize} -n 10 -p 1 
interaction $interactive


echo " DRAW FULL"
root -l -q -b "drawFullClosure.C(${conesize}, 1)"

interaction $interactive

echo " DRAW VTX"
root -l -q -b "getVtxReweighting.C(${conesize}, \"binning.config\")"

interaction $interactive

echo "PRIMER 2"

echo " FULL"
./createResponse_noempty_pp -c binning.config -r ${conesize} -n 10 -p 2 -h 0 -f ${fakes} -e 0

interaction $interactive

echo " DATA"
./unfoldData_noempty_pp -c binning.config -r ${conesize} -n 10 -p 2

interaction $interactive


echo " UNCERTAINTIES" 
./unfoldDataUncertainties_noempty_pp -c binning.config -r ${conesize} -n 10 -p 2

interaction $interactive

echo " DRAW FULL"
root -l -q -b "drawFullClosure.C(${conesize}, 2)"

interaction $interactive

echo " DRAW ITERATION"
root -l -q -b "makeIterationPlot_pp.C(${conesize}, 2)"
interaction $interactive

echo "NOMINAL "
echo " FULL"
./createResponse_noempty_pp -c binning.config -r ${conesize} -n 10 -p 0 -h 0 -f ${fakes} -e 0

interaction $interactive

echo " DATA"
./unfoldData_noempty_pp -c binning.config -r ${conesize} -n 10 -p 0

interaction $interactive

echo " UNCERTAINTIES"

./unfoldDataUncertainties_noempty_pp -c binning.config -r ${conesize} -n 10 -p 0 

interaction $interactive

echo " FULL"
root -l -q -b "drawFullClosure.C(${conesize}, 0)"

interaction $interactive

root -l -q -b "makeIterationPlot_pp.C(${conesize}, 0)"

if [ $onlynom -eq 0 ];
then
    echo " SYS"
    bash run_all_sys_pp.sh ${conesize} binning_negJES.config ${fakes}
    bash run_all_sys_pp.sh ${conesize} binning_posJES.config ${fakes}
    bash run_all_sys_pp.sh ${conesize} binning_negJER.config ${fakes}
    bash run_all_sys_pp.sh ${conesize} binning_posJER.config ${fakes}
    bash run_all_sys_pp.sh ${conesize} binning_herwig.config ${fakes}

    root -l -b -q "drawSysJESJER.C(${conesize})"
    root -l -b -q "drawFinalUnfold.C(${conesize})"
fi

