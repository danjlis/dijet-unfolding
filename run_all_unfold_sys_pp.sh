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

bash run_all_sys_pp.sh ${conesize} binning_negJES.config ${fakes} ${final}
bash run_all_sys_pp.sh ${conesize} binning_posJES.config ${fakes} ${final}
bash run_all_sys_pp.sh ${conesize} binning_negJER.config ${fakes} ${final}
bash run_all_sys_pp.sh ${conesize} binning_posJER.config ${fakes} ${final}
bash run_all_sys_pp.sh ${conesize} binning_herwig.config ${fakes} ${final}

