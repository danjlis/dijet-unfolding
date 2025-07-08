export DIJET_UNFOLDING_PATH=/sphenix/user/dlis/Projects/jet/ppg08-analysis/dijet-unfolding/
export DIJET_TNTUPLE_PATH=/sphenix/user/dlis/Projects/jet/tntuples/

if [[ ! -d final_plots ]];
then 
    mkdir final_plots
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
