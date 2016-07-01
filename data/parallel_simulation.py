# import os
#os.chdir('/data/emily/WF/SL_ABM')
#delete
# os.chdir('C:\\Users\\Emily Burchfield\\Dropbox\\Vanderbilt\\Dissertation\\Method\\ABM\\project')
import jg_agrowell_abm as abm
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import multiprocessing as mp
import sys
#
# If necessary, run
# pip install cython
# pip install feather-format
#
# or
# sudo -H pip install cython
# sudo -H pip install feather-format
#
import feather

n_years = 20
random.seed(1818)

###############################################################################
# SIMULATIONS
###############################################################################

fn = "data/payoff.csv"
# sim_dir = "."
# os.chdir(sim_dir)

payoff = abm.load_payoff(fn)


def one_sim(risk_model, i):
    print("Starting sim %d for risk model %d" % (i + 1, risk_model))
    field_canals = abm.build_field_canals(10, 15, payoff, risk_model=risk_model,
                                          risk_parameters=abm.RiskParameters())
    simulation = abm.Simulation(field_canals)
    simulation.simulate(n_years)
    sid = simulation.iter_data
    sid.loc[:,'simulation'] = pd.Series(i + 1, index = sid.index)
    # iter_data = iter_data.append(sid, ignore_index = True)
    print("Finishing sim %d for risk model %d" % (i + 1, risk_model))
    return sid

def do_simulations(n_sims = 1, risk_model = 0):
#def do_simulations(n_sims = 1, risk_model = 0):
    print("Starting job for risk model %d" % risk_model)
    pool = mp.Pool(50)
    iter_data = pd.DataFrame([])
    simulations = []
    res = [pool.apply_async(one_sim, (risk_model, i)) for i in range(n_sims)]
    try:
        sim_data = [x.get(timeout = 300) for x in res]
    except mp.TimeoutError:
        print("Simulation failed: timeout.")
        del(pool)
        sys.exit(1)
    print("Pool.map finished for risk_model %d" % risk_model)
    for sid in sim_data:
        iter_data = iter_data.append(sid, ignore_index = True)
    print("Finishing job for risk model %d" % risk_model)
    return iter_data

eu_results = do_simulations(100, risk_model = 0)
eu_iter_data = eu_results
# eu_iter_data.to_csv('eu_iter_data_p.csv')
feather.write_dataframe(eu_iter_data, 'eu_iter_data_p.feather')

regret_results = do_simulations(100, risk_model = 1)
regret_iter_data = regret_results
# regret_iter_data.to_csv('regret_iter_data_p.csv')
feather.write_dataframe(regret_iter_data, 'regret_iter_data_p.feather')

prospect_results = do_simulations(100, risk_model = 2)
prospect_iter_data = prospect_results
#prospect_iter_data.to_csv('prospect_iter_data_p.csv')
feather.write_dataframe(prospect_iter_data, 'prospect_iter_data_p.feather')

mixed_results = do_simulations(100, risk_model = 3)
mixed_iter_data = mixed_results
# mixed_iter_data.to_csv('mixed_iter_data_p.csv')
feather.write_dataframe(mixed_iter_data, 'mixed_iter_data_p.feather')

