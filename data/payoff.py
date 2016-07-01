# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 14:51:28 2016

@author: Emily Burchfield
"""
import os
import numpy as np

payoff_dir = "data/"

paddy_price = 100000
ofc_price   = 120000
sd = np.repeat(1000, 72)

#multiplicative effects on price received
agrowell_effect          = 1.2  # 1 is no effect
water_scarce_agrowell    = 0.7  # i.e. reduces profit by 30 percent
water_scarce_no_agrowell = 0.2  # low tank, low rainfall, reduces profit by 80 percent
water_damage_effect      = 0.2  # reduces ofc profit by 0.

bethma_effect            = 0.5  # reduces profit by 50 percent
ofc_bethma_effect        = 0.1  # slashes payoffs, just to ensure no one opts for this option

#additive bethma effects
agrowell_bm_effect       = 0.1  # changes effect on profit by 0.1
high_rf_bm_effect        = 0.1  # changes effect on profit by 0.1

#build dataset
tank_level = np.repeat((1,2,3), 24)
rainfall = np.tile(np.repeat((1,2,3), 8), 3)
crop_choice = np.tile(np.repeat((1,0), 4), 9)
bethma_choice = np.tile((1,0), 36)
agrowell_owner = np.tile((1,0,0,1), 18)
data = np.transpose(np.vstack((tank_level, rainfall, crop_choice, bethma_choice, agrowell_owner)))

#indices
TANK = 0
RF = 1
CROP = 2
BETHMA = 3
AGROWELL = 4
OFC_CROP = 1
PADDY_CROP = 0

#ofc/paddy
ofc = data[data[:,CROP] == OFC_CROP]
paddy = data[data[:,CROP] == PADDY_CROP]
assert ofc.shape[0] == paddy.shape[0]
ofc_payoff = np.zeros(ofc.shape[0])
paddy_payoff = np.zeros(paddy.shape[0])

###############################################################################
#ofc payoffs
###############################################################################

ofc_payoff[(ofc[:,BETHMA] == 1)] = ofc_price * ofc_bethma_effect

########WATER SCARCE CONDITIONS########

#low tank low rf with agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 1) & (ofc[:,RF] == 1) & (ofc[:,AGROWELL] == 1)] = ofc_price * water_scarce_agrowell
#low tank low rf without agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 1) & (ofc[:,RF] == 1) & (ofc[:,AGROWELL] == 0)] = ofc_price * water_scarce_no_agrowell

#low tank avg rf with agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 1) & (ofc[:,RF] == 2) & (ofc[:,AGROWELL] == 1)] = ofc_price * (water_scarce_agrowell + 0.1)
#low tank avg rf without agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 1) & (ofc[:,RF] == 2) & (ofc[:,AGROWELL] == 0)] = ofc_price * (water_scarce_no_agrowell + 0.1)

#low tank high rf with agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 1) & (ofc[:,RF] == 3) & (ofc[:,AGROWELL] == 1)] = ofc_price * (water_scarce_agrowell + 0.2)
#low tank high rf without agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 1) & (ofc[:,RF] == 3) & (ofc[:,AGROWELL] == 0)] = ofc_price * (water_scarce_no_agrowell + 0.2)


########NORMAL CONDITIONS########

#avg tank low rf with agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 2) & (ofc[:,RF] == 1) & (ofc[:,AGROWELL] == 1)] = ofc_price * agrowell_effect
#avg tank low rf without agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 2) & (ofc[:,RF] == 1) & (ofc[:,AGROWELL] == 0)] = ofc_price

#avg tank avg rf with agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 2) & (ofc[:,RF] == 2) & (ofc[:,AGROWELL] == 1)] = ofc_price * agrowell_effect
#avg tank average rf without agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 2) & (ofc[:,RF] == 2) & (ofc[:,AGROWELL] == 0)] = ofc_price

#high tank low rf with agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 3) & (ofc[:,RF] == 1) & (ofc[:,AGROWELL] == 1)] = ofc_price * agrowell_effect
#high tank low rf without agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 3) & (ofc[:,RF] == 1) & (ofc[:,AGROWELL] == 0)] = ofc_price

#high tank avg rf with agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 3) & (ofc[:,RF] == 2) & (ofc[:,AGROWELL] == 1)] = ofc_price * agrowell_effect
#high tank avg rf without agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 3) & (ofc[:,RF] == 2) & (ofc[:,AGROWELL] == 0)] = ofc_price


#########WATER ABUNDANT CONDITIONS########
#high tank high rf with agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 3) & (ofc[:,RF] == 3) & (ofc[:,AGROWELL] == 1)] = ofc_price * (1 - water_damage_effect)
#high tank high rf without agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 3) & (ofc[:,RF] == 3) & (ofc[:,AGROWELL] == 0)] = ofc_price * (1 - water_damage_effect)

#average tank high rf
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 2) & (ofc[:,RF] == 3) & (ofc[:,AGROWELL] == 1)] = ofc_price * (1 - water_damage_effect)
#high tank high rf without agrowell
ofc_payoff[(ofc[:,BETHMA] == 0) & (ofc[:,TANK] == 2) & (ofc[:,RF] == 3) & (ofc[:,AGROWELL] == 0)] = ofc_price * (1 - water_damage_effect)


###############################################################################
#paddy payoffs
###############################################################################

#bethma payoffs half of regular payoffs
paddy_payoff[(paddy[:,BETHMA] == 1)] = paddy_price * bethma_effect
#...unless you have high rainfall...
paddy_payoff[(paddy[:,BETHMA] == 1) & (paddy[:,RF] == 3)] = paddy_price * (bethma_effect + high_rf_bm_effect)
#unlike OFC, agrowells don't help you with paddy cultivation during water scarce periods

#non-bethma payoffs

########WATER SCARCE CONDITIONS########

#low tank low rf
paddy_payoff[(paddy[:,BETHMA] == 0) & (paddy[:,TANK] == 1) & (paddy[:,RF] == 1)] = paddy_price * water_scarce_no_agrowell  #could make this zero

#low tank avg rf
paddy_payoff[(paddy[:,BETHMA] == 0) & (paddy[:,TANK] == 1) & (paddy[:,RF] == 2)] = paddy_price * (water_scarce_no_agrowell + 0.1)

#low tank high rf
paddy_payoff[(paddy[:,BETHMA] == 0) & (paddy[:,TANK] == 1) & (paddy[:,RF] == 3)] = paddy_price * (water_scarce_no_agrowell + 0.2)


########NORMAL CONDITIONS########

#avg tank any rainfall
paddy_payoff[(paddy[:,BETHMA] == 0) & (paddy[:,TANK] == 2)] = paddy_price

#high tank any rainfall
#avg tank avg rf
paddy_payoff[(paddy[:,BETHMA] == 0) & (paddy[:,TANK] == 3)] = paddy_price


#final payoff table
paddy_final = np.hstack((paddy, paddy_payoff.reshape((paddy.shape[0],1))))
ofc_final = np.hstack((ofc, ofc_payoff.reshape((ofc.shape[0],1))))

final_payoff = np.vstack((paddy_final, ofc_final))
final_payoff = np.hstack((final_payoff, sd.reshape((72, 1))))


d = "data/"
np.savetxt(os.path.join(d ,"payoff.csv"), final_payoff, delimiter = ",")
