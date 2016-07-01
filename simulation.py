# import os
#os.chdir('/data/emily/WF/SL_ABM')
#delete
# os.chdir('C:\\Users\\Emily Burchfield\\Dropbox\\Vanderbilt\\Dissertation\\Method\\ABM\\project')
import agrowell_abm as abm
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
# from multiprocessing import Pool
# from functools import partial

n_years = 20
random.seed(1818)

###############################################################################
# SIMULATIONS
###############################################################################

fn = "data/payoff.csv"
# sim_dir = "."
# os.chdir(sim_dir)

payoff = abm.load_payoff(fn)

# pool = Pool(4)

def one_sim(risk_model, i):
    print("Starting sim %d for risk model %d" % (i + 1, risk_model))
    field_canals = abm.build_field_canals(10, 15, payoff, risk_model=risk_model,
                                          risk_parameters=abm.RiskParameters())
    simulation = abm.Simulation(field_canals)
    simulation.simulate(n_years)
    sid = simulation.iter_data
    sid.loc[:,'simulation'] = pd.Series(i + 1, index = sid.index)
    # iter_data = iter_data.append(sid, ignore_index = True)
    return sid

#def do_simulations(n_sims = 1, risk_model = 0, pool = pool):
def do_simulations(n_sims = 1, risk_model = 0):
    print("Starting job for risk model %d" % risk_model)
    iter_data = pd.DataFrame([])
    simulations = []
    #f = partial(one_sim, risk_model)
    #sim_data = pool.map(f, range(n_sims))
    # for sid in sim_data:
    #    iter_data = iter_data.append(sid, ignore_index = True)
    for i in range(n_sims):
        print("Replication %d for risk model %d" % (i + 1, risk_model))
        field_canals = abm.build_field_canals(10, 15, payoff, risk_model=risk_model,
                                              risk_parameters=abm.RiskParameters())
        simulation = abm.Simulation(field_canals)
        simulation.simulate(n_years)
        simulations.append(simulation)
        sid = simulation.iter_data
        sid.loc[:,'simulation'] = pd.Series(i + 1, index = sid.index)
        iter_data = iter_data.append(sid, ignore_index = True)
    return {'simulations':simulations, 'iter_data':iter_data}

eu_results = do_simulations(100, risk_model = 0)
eu_iter_data = eu_results['iter_data']
eu_iter_data.to_csv('eu_iter_data.csv')

regret_results = do_simulations(100, risk_model = 1)
regret_iter_data = regret_results['iter_data']
regret_iter_data.to_csv('regret_iter_data.csv')

prospect_results = do_simulations(100, risk_model = 2)
prospect_iter_data = prospect_results['iter_data']
prospect_iter_data.to_csv('prospect_iter_data.csv')

mixed_results = do_simulations(100, risk_model = 3)
mixed_iter_data = mixed_results['iter_data']
mixed_iter_data.to_csv('mixed_iter_data.csv')
