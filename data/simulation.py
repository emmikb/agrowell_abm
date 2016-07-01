import os
#os.chdir('/data/emily/WF/SL_ABM')
#delete
# os.chdir('C:\\Users\\Emily Burchfield\\Dropbox\\Vanderbilt\\Dissertation\\Method\\ABM\\project')
import agrowell_abm as abm
import matplotlib.pyplot as plt
import numpy as np
import random

n=20
#random.seed(1818)

###############################################################################
# SIMULATIONS
###############################################################################

fn = "data/payoff.csv"
# sim_dir = "."
# os.chdir(sim_dir)

obj = abm.Iteration(100, fn, risk_parameters = abm.RiskParameters(), 
                    risk_model = 3)  
obj.simulation(n)
#obj.iter_data.to_csv("iter_data.csv")





# PLOT CONSTRUCTION

tank = obj.iter_data['tank']
rf = obj.iter_data['rf']
aw = obj.iter_data['aw']
cd = obj.iter_data['crop']
rev = obj.iter_data['profit']
bm = obj.iter_data['bethma']
itr = obj.iter_data['iter']

def time_prep(var):
    sim = []    
    for i in np.unique(itr):
        season = var[itr == i]
        
        if len(np.unique(season)) > 2:
            sim.append(np.transpose(np.vstack((np.mean(season), np.std(season)))))
    else:
            unique, counts = np.unique(season, return_counts=True)
            if len(unique) == 1 and unique == 0:               
                counts = np.asarray((obj.number_agents, 0))
            if len(unique) == 1 and unique == 1:
                counts = np.asarray((0, obj.number_agents))
            else:
                counts = counts
            per = np.divide(counts.astype(float), obj.number_agents)
            sim.append(per)
    return np.asarray(sim)


time_plots()


################################################################################
##variation across time 
################################################################################
#
def time_plots():
    cd_data = time_prep(cd)
    bm_data = time_prep(bm)
    aw_data = time_prep(aw)
    rev_data = time_prep(rev).reshape((n,2))
    tnk = tank[0::obj.number_agents]
    rfll = rf[0::obj.number_agents]
    
    plt.subplot(231)
    plt.plot(range(n), cd_data[:,0], 'g-', label = "Paddy")   
    plt.plot(range(n), cd_data[:, 1], label = "OFC") 
    plt.ylabel("Percent") 
    plt.ylim([0,1])
    legend = plt.legend(loc = "upper right") 
    plt.title("Crop decisions") 
    plt.xticks(np.arange(0, n, 2))
    
    plt.subplot(232)
    plt.bar(range(n), bm_data[:,1], label = "Bethma")
    plt.xticks(np.arange(0, n, 2))
    plt.ylabel("Percent")  
    plt.title("Bethma") 
    
    plt.subplot(233)
    plt.plot(range(n), rev_data[:,0])  
    plt.fill_between(range(n), rev_data[:,0]+rev_data[:,1], rev_data[:,0]-rev_data[:,1], facecolor='blue', alpha = 0.1)  
    plt.ylabel("Average income") 
    plt.title("Income") 
    plt.xticks(np.arange(0, n, 2))
    
    plt.subplot(234) 
    plt.plot(range(n), aw_data[:, 1], label = "Agrowell ownership") 
    plt.ylim([0,1])
    plt.ylabel("Percent")  
    plt.title("Agrowell") 
    plt.xticks(np.arange(0, n, 2))
    
    plt.subplot(235)
    plt.bar(range(n), tnk)
    plt.xlabel("Iteration")
    plt.ylabel("Tank level")
    plt.xticks(np.arange(0, n, 2))
    plt.yticks((1, 2, 3))
    plt.title("Tank level")
    
    plt.subplot(236)
    plt.bar(range(n), rfll)
    plt.xticks(np.arange(0, n, 2))
    plt.yticks((1, 2, 3))
    plt.xlabel("Iteration")
    plt.ylabel("Rainfall amount")
    plt.title("Rainfall")
    


################################################################################
#variation across water availability 
################################################################################

#matplotlib renderer not working correctly ??  needs work

##categoricals
#data = tank_prep(cd)  
#fig, ax = plt.subplots()
#width = 0.35
#b1 = ax.bar(np.arange(3), data[:,0], width, color='g', yerr=data[:,2])
#x = np.arange(3) + width
#b2 = ax.bar(x, data[:,1], width, color='b', yerr=data[:,3])
#ax.set_xlim(-width, len(x)+width)
#ax.set_ylabel('Counts')
#ax.set_xticks(x + width)
#ax.set_xticklabels(('Low tank', 'Avg. tank', 'High tank'))
#ax.legend((b1[0], b2[0]), ('Paddy', 'OFC'))
#plt.show() 
#
##continuous
#data = tank_prep(rev)  
#fig, ax = plt.subplots()
#width = 0.35
#b1 = ax.bar(np.arange(3), data[0,:], width, color='g', yerr=data[1,:])
#ax.set_ylabel('Counts')
#ax.set_xticklabels(('Low tank', 'Avg. tank', 'High tank'))
#ax.legend((b1[0], b2[0]), ('Paddy', 'OFC'))
#plt.show()

#def tank_prep(var):
#
#    if len(np.unique(var)) > 2:
#        tank_r = np.tile(tank, (var.shape[0],1))
#        mn_tnk =[]
#        sd_tnk = []
#        for i in range(1,4):
#            subset = var[tank_r == i]
#            mn_tnk.append(np.mean(subset, axis=0))
#            sd_tnk.append(np.std(subset, axis=0))
#        return np.vstack((np.asarray(mn_tnk), np.asarray(sd_tnk)))
#    
#    else:
#        mn_tnk = []
#        sd_tnk = []
#        freq = []
#        for i in range(var.shape[1]):
#            unique, counts = np.unique(var[:,i], return_counts=True)  
#            freq.append(counts)
#        freq = np.asarray(freq)
#        
#        for i in range(1,4):
#            idx = np.where(tank == i)[1]
#            subset = freq[idx,:]
#            mn = np.mean(subset, axis=0)
#            sd = np.std(subset, axis=0)
#            mn_tnk.append(mn)
#            sd_tnk.append(sd)
#        return np.hstack((np.asarray(mn_tnk), np.asarray(sd_tnk)))