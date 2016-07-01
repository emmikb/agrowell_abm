import sys
sys.path.insert(0, '/home/emily/SL_ABM')
import agrowell_abm as abm
import numpy as np
import pandas as pd
import random
import multiprocessing as mp

d = '/data/emily/SL_ABM_sim_data/'
n_years = 20
random.seed(1818)

###############################################################################
# SIMULATIONS
###############################################################################

fn = "/home/emily/SL_ABM/data/payoff.csv"
payoff = abm.load_payoff(fn)
     
def one_sim(risk_model, i, parameter_name, new_value):
	print("Starting sim %d for risk model %d" % (i + 1, risk_model))
	risk_parameters = abm.RiskParameters()
	
	nparam = len(parameter_name)
	if nparam == 2:
		setattr(risk_parameters, parameter_name[0], new_value[0])
		setattr(risk_parameters, parameter_name[1], new_value[1])		
	if nparam == 3:
		setattr(risk_parameters, parameter_name[0], new_value[0])
		setattr(risk_parameters, parameter_name[1], new_value[1])
		setattr(risk_parameters, parameter_name[2], new_value[2])
	else:
		setattr(risk_parameters, parameter_name[0], new_value[0])	#new_value[i]
	
	field_canals = abm.build_field_canals(10, 15, payoff, risk_model=risk_model, risk_parameters=risk_parameters)
	simulation = abm.Simulation(field_canals)
	simulation.simulate(n_years)
	sid = simulation.iter_data
	sid.loc[:,'simulation'] = pd.Series(i + 1, index = sid.index)
	# iter_data = iter_data.append(sid, ignore_index = True)
	print("Finishing sim %d for risk model %d" % (i + 1, risk_model))
	return sid

def do_simulations(n_sims, risk_model, parameter_name, new_value):
	print("Starting job for risk model %d" % risk_model)
	pool = mp.Pool(50)
	iter_data = pd.DataFrame([])
	simulations = []
	res = [pool.apply_async(one_sim, (risk_model, i, parameter_name, new_value)) for i in range(n_sims)]

	try:
		sim_data = [x.get(timeout = 300) for x in res]
	except mp.TimeoutError:
		print("Simulation failed: timeout.")
		del(pool)
		sys.exit(1)
	nparam = len(parameter_name)
	
	if nparam == 1:
		print("Pool.map finished for risk_model %d with parameter %s set as %s." % (risk_model, parameter_name, new_value))  
	if nparam == 2:	
		print("Pool.map finished for risk_model %d with parameter %s set as %s and parameter %s set as %s." % (risk_model, parameter_name[0], new_value[0], parameter_name[1], new_value[1]))
	if nparam == 3:	
		print("Pool.map finished for risk_model %d with parameter %s set as %s, parameter %s set as %s, and parameter %s set as %s." % (risk_model, parameter_name[0], new_value[0], parameter_name[1], new_value[1], parameter_name[2], new_value[2]))
	  

	for sid in sim_data:
		iter_data = iter_data.append(sid, ignore_index = True)   
	
	if nparam == 1:
		print("Finished job for risk_model %d with parameter %s set as %s." % (risk_model, parameter_name, new_value))
	if nparam == 2:	
		print("Finished job for risk_model %d with parameter %s set as %s and parameter %s set as %s." % (risk_model, parameter_name[0], new_value[0], parameter_name[1], new_value[1]))
	if nparam == 3:
		print("Finished job for risk_model %d with parameter %s set as %s, parameter %s set as %s, and parameter %s set as %s." % (risk_model, parameter_name[0], new_value[0], parameter_name[1], new_value[1], parameter_name[2], new_value[2]))
	  
	return iter_data

def parameter_test(nsim, risk_model, parameter_name, value_list):	

	nparam = len(parameter_name)
	
	if nparam == 3:
		parameter_combo = [[value_list[0][x], value_list[1][y], value_list[2][z]] for x in range(len(value_list[0])) for y in range(len(value_list[1])) for z in range(len(value_list[2]))]
		for i in range(len(parameter_combo)):
			iter_data = do_simulations(nsim, risk_model, parameter_name, parameter_combo[i])
			iter_data[str(parameter_name[0])] = pd.Series(np.repeat(parameter_combo[i][0], len(iter_data)))
			iter_data[str(parameter_name[1])] = pd.Series(np.repeat(parameter_combo[i][1], len(iter_data)))
			iter_data[str(parameter_name[2])] = pd.Series(np.repeat(parameter_combo[i][2], len(iter_data)))
			iter_data.to_csv(str(d) + 'rm' + str(risk_model) + '_' + str(parameter_name[0]) + str(parameter_combo[i][0]) + '_' + str(parameter_name[1]) + str(parameter_combo[i][1]) + '_' + str(parameter_name[2]) + str(parameter_combo[i][2]) + '.csv')  
	
	if nparam == 2:	
		parameter_combo = [[value_list[0][x], value_list[1][y]] for x in range(len(value_list[0])) for y in range(len(value_list[1]))]
		for i in range(len(parameter_combo)):
			iter_data = do_simulations(nsim, risk_model, parameter_name, parameter_combo[i])
			iter_data[str(parameter_name[0])] = pd.Series(np.repeat(parameter_combo[i][0], len(iter_data)))
			iter_data[str(parameter_name[1])] = pd.Series(np.repeat(parameter_combo[i][1], len(iter_data)))
			iter_data.to_csv(str(d) + 'rm' + str(risk_model) + '_' + str(parameter_name[0]) + str(parameter_combo[i][0]) + '_' + str(parameter_name[1]) + str(parameter_combo[i][1]) + '.csv')  
	
	if nparam == 1: 
		for i in range(len(value_list[0])):
			results = do_simulations(nsim, risk_model, parameter_name, [value_list[0][i]])
			iter_data = results
			iter_data[str(parameter_name[0])] = pd.Series(np.repeat(value_list[0][i], len(iter_data)))
			iter_data.to_csv(str(d) + 'rm' + str(risk_model) + '_' + str(value_list[0][i]) + str(parameter_name[0]) + '.csv')  

r_list = [-0.5, 0, 0.5, 1, 2, 3, 4]
k_list = [0, 0.155, 0.564] #Figure 3 values  
b_list = [0.0005, 0.5, 0.9]  #Figure 3 values
rr_list = [0, 1, 2, 3, 4] #Figure 3 values
alph_list = [0.6, 0.7, 0.8, 0.9, 1.0]  #Figure 4 values
lam_list = [1, 2.25, 3.5]

#eu
#parameter_test(100, 0, ['r'], [r_list])  #parameter_name and value_list as list

#regret-corrected
#parameter_test(100, 1, ['b','k', 'rr'], [b_list, k_list, rr_list])

#prospect theory
parameter_test(100, 2, ['lam','alph'], [lam_list, alph_list])  

#feather R - Python
#plot results


