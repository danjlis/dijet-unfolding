export DIJET_UNFOLDING_PATH=/sphenix/user/dlis/Projects/jet/ppg08-analysis/dijet-unfolding/
export DIJET_TNTUPLE_PATH=/sphenix/user/dlis/Projects/jet/tntuples/

export AUAU_DATA_PATH=/sphenix/tg/tg01/jets/dlis/data/v81/
export AUAU_SIM_PATH=/sphenix/tg/tg01/jets/dlis/hijing/v10/

export PP_SIM_PATH=/sphenix/tg/tg01/jets/dlis/sim/pythia/
export PP_DATA_PATH=/sphenix/tg/tg01/jets/dlis/data/all/

if [[ ! -d final_plots ]];
then 
    mkdir final_plots
    mkdir auau_plots
    mkdir jer_plots
    mkdir njet
    mkdir response_matrices
    mkdir setup_env.sh
    mkdir systematic_plots
    mkdir truth_hists
    mkdir uncertainties
    mkdir unfolding_hists
    mkdir unfolding_plots
    mkdir vertex
fi
