export DIJET_UNFOLDING_PATH=/Users/daniel/sPHENIX/ppg08/ppg08-analysis/dijet-unfolding/
export DIJET_TNTUPLE_PATH=/Users/daniel/sPHENIX/ppg08/trees/

export AUAU_DATA_PATH=/Users/daniel/sPHENIX/ppg08/trees/
export AUAU_SIM_PATH=/Users/daniel/sPHENIX/ppg08/trees/

export JESR_PATH=/Users/daniel/sPHENIX/ppg08/JESR/

export PP_SIM_PATH=/Users/daniel/sPHENIX/ppg08/trees/
export PP_DATA_PATH=/Users/daniel/sPHENIX/ppg08/trees/

if [[ ! -d final_plots ]];
then 
    mkdir final_plots
    mkdir auau_plots
    mkdir jer_plots
    mkdir jer
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
