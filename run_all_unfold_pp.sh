#!/usr/bin/env bash
source ~/.profile

conesize=$1
fakes=$2
onlynom=$3
interactive=$4
final=$5
interaction () {

    if [[ $1 -gt 0 ]]; then
	read -n 1 -s -r -p "Press any key to continue..."
    fi
}


if [[ $final -eq 0 ]]; then
    ./createResponse_noempty_pp -r ${conesize} -c "binning.config" -n 10 -p 1 -h 0 -f ${fakes} -e 0 -v 2 -m 1

    ./createResponse_noempty_pp -c binning.config -r ${conesize} -n 10 -p 1 -h 0 -f ${fakes} -e 0

    interaction $interactive

    ./unfoldData_noempty_pp -c binning.config -r ${conesize} -n 10 -p 1 

    interaction $interactive

    root -l -q -b "drawFullClosure.C(${conesize}, 1)"

    interaction $interactive

    root -l -q -b "getVtxReweighting.C(${conesize}, \"binning.config\")"

    interaction $interactive

    ./createResponse_noempty_pp -c binning.config -r ${conesize} -n 10 -p 2 -h 0 -f ${fakes} -e 0

    interaction $interactive

    ./unfoldData_noempty_pp -c binning.config -r ${conesize} -n 10 -p 2

    interaction $interactive

    ./unfoldDataUncertainties_noempty_pp -c binning.config -r ${conesize} -n 10 -p 2

    interaction $interactive

    root -l -q -b "drawFullClosure.C(${conesize}, 2)"

    interaction $interactive

    root -l -q -b "makeIterationPlot_pp.C(${conesize}, 2)"

    interaction $interactive
fi

echo "NOMINAL"
./createResponse_noempty_pp -c binning.config -r ${conesize} -n 10 -p 0 -h 0 -f ${fakes} -e 0
./createResponse_noempty_pp -c binning.config -r ${conesize} -n 10 -p 0 -h 1 -f ${fakes} -e 0

interaction $interactive

./unfoldData_noempty_pp -c binning.config -r ${conesize} -n 10 -p 0

interaction $interactive


if [[ $final -eq 0 ]]; then
    ./unfoldDataUncertainties_noempty_pp -c binning.config -r ${conesize} -n 10 -p 0 
    
    interaction $interactive
    
    root -l -q -b "drawFullClosure.C(${conesize}, 0)"
    root -l -q -b "drawHalfClosure_pp.C(${conesize}, 0)"
    
    interaction $interactive
    
    root -l -q -b "makeIterationPlot_pp.C(${conesize}, 0)"
    if [ $onlynom -eq 0 ];
    then
	echo " SYS"
	bash run_all_sys_pp.sh ${conesize} binning_negJES.config ${fakes} ${final}
	bash run_all_sys_pp.sh ${conesize} binning_posJES.config ${fakes} ${final}
	bash run_all_sys_pp.sh ${conesize} binning_negJER.config ${fakes} ${final}
	bash run_all_sys_pp.sh ${conesize} binning_posJER.config ${fakes} ${final}
	bash run_all_sys_pp.sh ${conesize} binning_herwig.config ${fakes} ${final}
	bash run_all_sys_pp.sh ${conesize} binning_trigger.config ${fakes} ${final}
	bash run_all_sys_pp.sh ${conesize} binning_horizontal.config ${fakes} ${final}
	bash run_all_sys_pp.sh ${conesize} binning_vertical.config ${fakes} ${final}
	bash run_all_sys_pp.sh ${conesize} binning_0.config ${fakes} ${final}
	bash run_all_sys_pp.sh ${conesize} binning_1.5.config ${fakes} ${final}
	root -l -b -q "drawSysJESJER.C(${conesize})"
	root -l -b -q "drawFinalUnfold.C(${conesize})"
    fi
fi

