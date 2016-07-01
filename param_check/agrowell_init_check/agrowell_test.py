# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 15:12:36 2016

@author: Emily Burchfield
"""
import numpy as np
import random
import itertools
import pandas as pd
from scipy.stats import rankdata

#account for market flood potential in payoff table

TANK = 0
RF = 1
CROP = 2
BETHMA = 3
AGROWELL = 4
PROFIT = 5
SD = 6
PADDY_CROP = 0
OFC_CROP = 1
PADDY_BM = 0
PADDY_NBM = 1
OFC_BM = 2
OFC_NBM = 3


class DecisionFunction(object):

    def __init__(self, x1=0.0, p1=0.1, x2=0.25, p2=0.9):
        d = np.log(p1/(1-p1))
        c = np.log(p2/(1-p2))
        self.b = (d-c)/(x1-x2)
        self.a = d - (self.b * x1)

    def p(self, x):
        z = np.exp(self.a + (self.b*x))
        p = z/(1+z)  #p value for given value of x
        return p

    def decide(self, x):
        prob = self.p(x)
        return prob > np.random.uniform(0, 1)

class CropForecast(object):

    def __init__(self, probabilities, risk_profile, income):
        if isinstance(probabilities, (list, tuple)):
            probabilities = np.asarray(probabilities)

        if probabilities.ndim != 1 or income.ndim != 1:
            raise ValueError("Probabilities and income must be 1-dimensional vectors.")

        if probabilities.shape != income.shape:
            raise ValueError("Probabilities and income must have the same length.")
        self.p = probabilities
        self.rp = risk_profile
        self.w = income

    # constant relative  risk aversion utility
    def crra_u(self, r):
        if abs(r - 1) < 1E-3:
            u = np.log(self.w)
        else:
            u = (np.power(self.w, 1.0 - r) - 1) / (1.0 - r)
        return u

    #estimated utility crra
    def eu_crra(self):
        u = self.crra_u(self.rp.r)
        return np.dot(self.p, u)

    #regret aversion
    def regret(self):
        rr = self.rp.rr
        k = self.rp.k
        b = self.rp.b
        u = self.crra_u(rr)
        u_max = max(u)
        delta_u = u_max - u
        g = 1 - np.power(b, delta_u)
        return np.dot(self.p, u - k * g)

    #prospect theory
    def prospect(self, w0):
        alph = self.rp.alph
        lam = self.rp.lam
        gam = self.rp.gam
        delta = self.rp.delta
        delta_w = self.w - w0
        h = np.ones(delta_w.shape)
        h[delta_w < 0] = - lam
        nu = h * np.power(np.absolute(delta_w), alph)
        x = np.repeat(gam, delta_w.shape)
        x[delta_w < 0] = delta
        p_gam = np.power(self.p, x)
        omega = p_gam / np.power(p_gam + np.power(1.0 - self.p, x), 1.0 / x)
        return np.dot(omega, nu) * self.rp.prospect_scale

def calc_profit(payoff, crop):
    po = payoff[payoff[:, CROP] == crop]
    profit = (po[po[:, RF] == 1, PROFIT],
            po[po[:, RF] == 2, PROFIT],
            po[po[:, RF] == 3, PROFIT])
    profit = np.asarray(profit).reshape((3,))
    return profit

def inspect_cf(payoff, risk_profile, tank, agrowell, bethma):
    p_rf = (0.25, 0.50, 0.25)

    po = payoff[(payoff[:, TANK] == 2) & (payoff[:, AGROWELL] == agrowell) & (payoff[:, BETHMA] == 0)]
    baseline = calc_profit(po, PADDY_CROP)[1]

    indices = (payoff[:, TANK] == tank)
    indices = np.logical_and(indices, payoff[:, AGROWELL] == agrowell)
    indices = np.logical_and(indices, payoff[:, BETHMA] == bethma)
    po = payoff[indices]
    crop = PADDY_CROP
    paddy_profit = calc_profit(po, crop)
    paddy_cf = CropForecast(p_rf, risk_profile, paddy_profit)
    crop = OFC_CROP
    ofc_profit = calc_profit(po, crop)
    ofc_cf = CropForecast(p_rf, risk_profile, ofc_profit)
    print ("Paddy income: (%10.0f, %10.0f, %10.0f)" % tuple(paddy_profit))
    print ("OFC income:   (%10.0f, %10.0f, %10.0f)" % tuple(ofc_profit))
    print ("Delta income: (%10.0f, %10.0f, %10.0f)" % tuple(ofc_profit - paddy_profit))
    print ("Paddy EU:     (%10.2f, %10.2f, %10.2f)" % (paddy_cf.eu_crra(),
                           paddy_cf.regret(), paddy_cf.prospect(baseline)))
    print ("OFC EU:       (%10.2f, %10.2f, %10.2f)" % (ofc_cf.eu_crra(),
                           ofc_cf.regret(), ofc_cf.prospect(baseline)))
    print ("Delta EU:     (%+10.2f, %+10.2f, %+10.2f)" % (ofc_cf.eu_crra() - paddy_cf.eu_crra(),
                          ofc_cf.regret() - paddy_cf.regret(),
                          ofc_cf.prospect(baseline) - paddy_cf.prospect(baseline)))

def check_bethma_utility(payoff, risk_profile, tank, agrowell):
    p_rf = (0.25, 0.50, 0.25)

    po = payoff[(payoff[:, TANK] == 2) & (payoff[:, AGROWELL] == agrowell) & (payoff[:, BETHMA] == 0)]
    baseline = calc_profit(po, PADDY_CROP)[1]

    indices = (payoff[:, TANK] == tank)
    indices = np.logical_and(indices, payoff[:, AGROWELL] == agrowell)
    po = payoff[np.logical_and(indices, payoff[:, BETHMA] == 1)]
    crop = PADDY_CROP
    paddy_profit = calc_profit(po, crop)
    paddy_cf = CropForecast(p_rf, risk_profile, paddy_profit)

    po = payoff[np.logical_and(indices, payoff[:, BETHMA] == 0)]
    crop = OFC_CROP
    ofc_profit = calc_profit(po, crop)
    ofc_cf = CropForecast(p_rf, risk_profile, ofc_profit)
    print ("Paddy income: (%10.0f, %10.0f, %10.0f)" % tuple(paddy_profit))
    print ("OFC income:   (%10.0f, %10.0f, %10.0f)" % tuple(ofc_profit))
    print ("Delta income: (%10.0f, %10.0f, %10.0f)" % tuple(ofc_profit - paddy_profit))
    print ("Paddy EU:     (%10.2f, %10.2f, %10.2f)" % (paddy_cf.eu_crra(),
                           paddy_cf.regret(), paddy_cf.prospect(baseline)))
    print ("OFC EU:       (%10.2f, %10.2f, %10.2f)" % (ofc_cf.eu_crra(),
                           ofc_cf.regret(), ofc_cf.prospect(baseline)))
    print ("Delta EU:     (%+10.2f, %+10.2f, %+10.2f)" % (ofc_cf.eu_crra() - paddy_cf.eu_crra(),
                          ofc_cf.regret() - paddy_cf.regret(),
                          ofc_cf.prospect(baseline) - paddy_cf.prospect(baseline)))

class RiskParameters(object):
    def __init__(self, r=1, rr=1, k=0.155, b=0.5, alph=0.88,
                 lam=2.55, gam=0.69, delta=0.61, prospect_scale=5.0E-5):
        self.r = r # CRRA relative risk aversion, nominal value 1
        self.rr = rr # regret-averse relative risk aversion, like r
        self.k = k   # regret aversion coefficient, nominal value 0.155
        self.b = b  # regret coefficient (0 <= b < 1), nominal value 0.5
        self.alph = alph # elasticity of marginal utility, nominal value 0.88
                         # Risk aversion for gain and risk-seeking for loss
        self.lam = lam   # loss aversion (> 1), nominal value 2.25
        self.gam = gam   # probability scaling for gains, nominal value 0.69
        self.delta = delta   # probability scaling for losses, nominal value 0.61
        self.prospect_scale = prospect_scale

class FieldCanal(object):
    def __init__(self, fc_id, farmers=None):
        """
        :param fc_id: ID flagging field canal to which a farmer belongs
        :param farmers: List of farmers who belong to this FC.
       """
        self.fc_id = fc_id
        self.farmers = list()
        if isinstance(farmers, (list, tuple)):
            self.add_farmers(farmers)
        # self.farmer_ids = list()
        self.bethma = None
        self.rainfall = None

    def add_farmer(self, farmer):
        farmer.join_field_canal(self)
        # self.farmer_ids.append(farmer.farmer_id)
        self.farmers.append(farmer)

    def add_farmers(self, new_farmers):
        for f in new_farmers:
            f.join_field_canal(self)
        # self.farmer_ids.extend([f.farmer_id] for f in new_farmers)
        self.farmers.extend(new_farmers)

    def decide_bethma(self, update_farmers=False, tank_level=None):
        if update_farmers:
            for f in self.farmers:
                assert (tank_level is not None and tank_level in (1, 2, 3)), \
                    "FieldCanal.decide_bethma called with update_farmers, but invalid tank_level."
                f.update_utility(tank_level)
        preferences = np.asarray([f.ranked_preferences for f in self.farmers], dtype=np.int)
        mask = np.ma.masked_where(preferences < 4, preferences)
        top_choice = np.divide(np.sum(mask, axis=0), 4)
        top_choice = np.ma.filled(top_choice, 0)

        #if most farmers want to do bethma, do bethma
        yes_votes = top_choice[PADDY_BM] + top_choice[OFC_BM]
        no_votes = top_choice[PADDY_NBM] + top_choice[OFC_NBM]
        self.bethma = (yes_votes > no_votes)

        if (self.fc_id == 1 and tank_level == 1):
            print("Bethma decision: %s. Tank = %d (%d, %d, %d, %d), (%d, %d)" % \
            (self.bethma, tank_level, top_choice[PADDY_BM], top_choice[OFC_NBM],
             top_choice[PADDY_NBM], top_choice[OFC_BM],
             yes_votes,
             no_votes))


    def count_ofc(self):
        return np.sum([f.crop_decision == OFC_CROP for f in self.farmers])

    def mean_ses(self):
        return np.mean([f.ses for f in self.farmers])

    def start_growing_season(self, tank_level, rf_forecast):
        for f in self.farmers:
            f.update_utility(tank_level, rf_forecast)
        self.decide_bethma(False, tank_level)
        for f in self.farmers:
            f.decide_cultivation(tank_level, rf_forecast)

    def finish_growing_season(self, tank_level, rainfall, num_ofc, mean_ses, std_ses):
        for f in self.farmers:
            f.harvest(tank_level, rainfall, num_ofc)
        for f in self.farmers:
            f.invest_in_agrowell(Farmer.SES_mean + 1 * Farmer.SES_std)
        for f in self.farmers:
            f.do_accounting()

class Farmer(object):
    """"
    :param agrowell_duration: indication of how many seasons farmer has owned agrowell
    :param agrowell_cost: cost of installing agrowell
    :param agrowell_fees: fees paid annually by farmer for ten years following installation of agrowell
    :param SES_mean: initial average socio-economic status, computed from survey data (Rs)
    :param SES_std: initial standard deviation of socio-economic status, computed from survey data (Rs)
    """
    SES_mean      = 100000
    SES_std       =  20000
    expense_mean  =  17500
    agrowell_cost =  70000
    agrowell_term =     10
    agrowell_fees = agrowell_cost/agrowell_term

    def __init__(self, farmer_id, payoff, risk_model, risk_parameters,
                 decision_fn, ses, agrowell=False):
        """
        :param farmer_id:  unique ID for each farmer agent
        :param payoff: payoff table
        :param risk_parameters: user input risk parameters or default parameters built with RiskParameter()
        :param risk_model: 0 (expected utility), 1 (regret), 2 (prospect), 3 (random assignment of 0, 1, 2)
        :param agrowell: binary flag indicating if farmer owns an agrowell (1) or not (0)
        :param costs: costs incurred by farmer, deducted seasonally from socio-economic status
        :param crop_decision: binary flag indicating whether the farmer plans paddy (0) or OFC (1)
        :param expected_income: expected income for the next season based on income of the past season
        :param eu: expected utility for one of four scenarios: paddy bethma, paddy non-bethma, ofc bethma, ofc non-bethma
        :param ses: agent socio-economic status (Rs)
        :param income: income earned at the end of each season (Rs)
        :param ranked_preference: matrix indicating prefered cultivation decision of farmers between paddy bethma, paddy non-bethma, ofc bethma, ofc non-bethma
       """
        self.farmer_id = farmer_id
        # self.fc_id = fc_id
        self.fc = None
        self.ses = ses
        self.payoff = payoff
        self.risk_objective = risk_model
        self.risk_parameters = risk_parameters
        self.decision_fn = decision_fn
        self.crop_decision = None
        self.expected_income = None
        self.expected_utility = None
        self.income = None
        self.wealth = None
        self.bethma_vote = None
        self.ranked_preferences = None
        self.agrowell = False
        self.agrowell_duration = 0
        if agrowell:
            self.build_agrowell(Farmer.agrowell_term)

    def join_field_canal(self, fc):
        self.fc = fc

    def build_agrowell(self, duration=0):
        self.agrowell = True
        self.agrowell_duration = duration

    def update_utility(self, tank_level, rf_forecast=np.asarray((0.25, 0.50, 0.25))):
        #paddy bethma
        p_bm = self.payoff[self.payoff[:, BETHMA] == 1]
        p_bm = p_bm[np.logical_and(p_bm[:, TANK] == tank_level, p_bm[:, CROP] == PADDY_CROP)]

        #paddy no bethma
        p_nb = self.payoff[self.payoff[:, BETHMA] == 0]
        p_nb = p_nb[np.logical_and(p_nb[:, TANK] == tank_level, p_nb[:, CROP] == PADDY_CROP)]

        #OFC bethma
        o_bm = self.payoff[self.payoff[:, BETHMA] == 1]
        o_bm = o_bm[np.logical_and(o_bm[:, TANK] == tank_level, o_bm[:, CROP] == OFC_CROP)]

        #OFC no bethma
        o_nb = self.payoff[self.payoff[:, BETHMA] == 0]
        o_nb = o_nb[np.logical_and(o_nb[:, TANK] == tank_level, o_nb[:, CROP] == OFC_CROP)]

        #individual payoff subset
        pb_i  = p_bm[p_bm[:, AGROWELL] == self.agrowell]
        pnb_i = p_nb[p_nb[:, AGROWELL] == self.agrowell]

        ob_i  = o_bm[o_bm[:, AGROWELL] == self.agrowell]
        onb_i = o_nb[o_nb[:, AGROWELL] == self.agrowell]

        #construct forecasts
        paddy_bm = (pb_i[pb_i[:, RF] == 1, PROFIT],
                    pb_i[pb_i[:, RF] == 2, PROFIT],
                    pb_i[pb_i[:, RF] == 3, PROFIT])
        paddy_bm = np.asarray((paddy_bm)).reshape((len(paddy_bm)))
        paddy_bm_fc = CropForecast(rf_forecast, self.risk_parameters, paddy_bm)

        paddy_nb = (pnb_i[pnb_i[:, RF] == 1, PROFIT],
                    pnb_i[pnb_i[:, RF] == 2, PROFIT],
                    pnb_i[pnb_i[:, RF] == 3, PROFIT])
        paddy_nb = np.asarray((paddy_nb)).reshape((len(paddy_nb)))
        paddy_nb_fc = CropForecast(rf_forecast, self.risk_parameters, paddy_nb)

        ofc_bm = (ob_i[ob_i[:, RF] == 1, PROFIT],
                    ob_i[ob_i[:, RF] == 2, PROFIT],
                    ob_i[ob_i[:, RF] == 3, PROFIT])
        ofc_bm = np.asarray((ofc_bm)).reshape((len(ofc_bm)))
        ofc_bm_fc = CropForecast(rf_forecast, self.risk_parameters, ofc_bm)

        ofc_nb = (onb_i[onb_i[:, RF] == 1, PROFIT],
                    ob_i[onb_i[:, RF] == 2, PROFIT],
                    onb_i[onb_i[:, RF] == 3, PROFIT])
        ofc_nb = np.asarray((ofc_nb)).reshape((len(ofc_nb)))
        ofc_nb_fc = CropForecast(rf_forecast, self.risk_parameters, ofc_nb)

        if self.risk_objective == 0:
            eu_paddy_bm = paddy_bm_fc.eu_crra()
            eu_paddy_nb = paddy_nb_fc.eu_crra()
            eu_ofc_bm = ofc_bm_fc.eu_crra()
            eu_ofc_nb = ofc_nb_fc.eu_crra()
        elif self.risk_objective == 1:
            eu_paddy_bm = paddy_bm_fc.regret()
            eu_paddy_nb = paddy_nb_fc.regret()
            eu_ofc_bm = ofc_bm_fc.regret()
            eu_ofc_nb = ofc_nb_fc.regret()
        elif self.risk_objective == 2:
            if self.expected_income is None:
                self.expected_income = paddy_nb[1]
                # if expected income hasn't been initialized,
                # then set it to non-bethma paddy for a normal season
            w0 = self.expected_income # reference = last season's income
            eu_paddy_bm = paddy_bm_fc.prospect(w0)
            eu_paddy_nb = paddy_nb_fc.prospect(w0)
            eu_ofc_bm = ofc_bm_fc.prospect(w0)
            eu_ofc_nb = ofc_nb_fc.prospect(w0)
        else:
            assert False, ("Bad value for risk objective %d for farmer %d on field canal %d" %
                            (self.risk_objective, self.farmer_id, self.fc.fc_id))

        self.expected_utility = (eu_paddy_bm, eu_paddy_nb, eu_ofc_bm, eu_ofc_nb)
        self.ranked_preferences = rankdata(self.expected_utility)

        if False and self.farmer_id == 1:
            print("payoffs: paddy =  %.2f, paddy bethma = %.2f, ofc = %.2f, ofc bethma = %.2f" % \
                (np.dot(paddy_nb, rf_forecast), np.dot(paddy_bm, rf_forecast),
                 np.dot(ofc_nb, rf_forecast), np.dot(ofc_bm, rf_forecast)))
            print("utilities: paddy =  %.2f, paddy bethma = %.2f, ofc = %.2f, ofc bethma = %.2f" % \
                (eu_paddy_nb, eu_paddy_bm, eu_ofc_nb, eu_ofc_bm))


    def decide_cultivation(self, tank_level, rf_forecast=(0.25, 0.50, 0.25)):

        #recompute preferences given information about bethma decision
        po = self.payoff[self.payoff[:, TANK] == tank_level]
        paddy_subset = po[po[:, CROP] == PADDY_CROP]
        ofc_subset = po[po[:, CROP] == OFC_CROP]

        paddy = paddy_subset[np.logical_and(paddy_subset[:, BETHMA] == self.fc.bethma,
                                           paddy_subset[:, AGROWELL] == self.agrowell)]

        ofc = ofc_subset[np.logical_and(ofc_subset[:, BETHMA] == self.fc.bethma,
                                         ofc_subset[:, AGROWELL] == self.agrowell)]

        w_paddy = ( paddy[paddy[:, RF] == 1, PROFIT],
                    paddy[paddy[:, RF] == 2, PROFIT],
                    paddy[paddy[:, RF] == 3, PROFIT] )
        w_paddy = np.asarray((w_paddy)).reshape(len(w_paddy))
        paddy_forecast = CropForecast(rf_forecast, self.risk_parameters, w_paddy)

        w_ofc = (ofc[ofc[:, RF] == 1, PROFIT],
                 ofc[ofc[:, RF] == 2, PROFIT], ofc[ofc[:, RF] == 3, PROFIT])
        w_ofc = np.asarray((w_ofc)).reshape(len(w_ofc))
        ofc_forecast = CropForecast(rf_forecast, self.risk_parameters, w_ofc)

        if self.risk_objective == 0:
            eu_paddy = paddy_forecast.eu_crra()
            eu_ofc   = ofc_forecast.eu_crra()
        elif self.risk_objective == 1:
            eu_paddy = paddy_forecast.regret()
            eu_ofc   = ofc_forecast.regret()
        elif self.risk_objective == 2:
            w0 = self.expected_income
            eu_paddy = paddy_forecast.prospect(w0)
            eu_ofc   = ofc_forecast.prospect(w0)

        #BERNOULLI INVERSE LOGIT HERE
        eu_diff = (eu_ofc - eu_paddy) # /1000  #/1000 to fit logit
        decision = self.decision_fn.decide(eu_diff)
        if False and self.farmer_id == 1:
            print ("risk model = %d, eu_ofc = %.3g, eu_paddy = %.3g, eu_diff = %.3g, p = %.3f, decision = %d" % \
            (self.risk_objective, eu_ofc, eu_paddy, eu_diff, self.decision_fn.p(eu_diff), decision))
        if self.decision_fn.decide(eu_diff):
            self.crop_decision = OFC_CROP
        else:
            self.crop_decision = PADDY_CROP

    def harvest(self, tank_level, rainfall, num_ofc):

        #payoffs for each crop based on rainfall and tank level
        po = self.payoff
        indices = np.logical_and(po[:, RF] == rainfall,
                                 po[:, TANK] == tank_level)
        indices = np.logical_and(indices, po[:, CROP] == self.crop_decision)
        indices = np.logical_and(indices, po[:, AGROWELL] == self.agrowell)
        indices = np.logical_and(indices, po[:, BETHMA] == self.fc.bethma)
        po = po[indices]
        # if (self.farmer_id == 1):
        #    print "Dim po = ",  po.shape
        self.income = po[0, PROFIT] + random.normalvariate(0, po[0, SD])
        # if (self.farmer_id == 1):
        #    print "income = ", self.income

        #OFC market flood
        #if sum(self.crop_decision[self.crop_decision == OFC_CROP]) > self.number_agents/2:
         #   self.income[self.crop_decision[i] == OFC_CROP] = self.income[i] * 0.7

        #expected income for next season as income from this season
        self.expected_income = self.income

    def calc_agrowell_fee(self):
        fees = 0
        if self.agrowell and self.agrowell_duration <= self.agrowell_term:
            fees = self.agrowell_fees
        return fees

    def invest_in_agrowell(self, threshold):
        #adopt agrowell if $$$
        if (not self.agrowell) and self.ses > threshold:  # + np.std(self.SES)):
            self.build_agrowell()
        # Adjust expectations for next season...
        # if self.agrowell_duration == 0:
        #     self.expected_income = self.expected_income - self.agrowell_fees
        # elif self.agrowell_duration == self.agrowell_term:
        #     self.expected_income = self.expected_income + self.agrowell_fees


    def do_accounting(self):
        # add revenue and deduct seasonal expenses from SES
        self.ses = self.ses + self.income * 0.2 - self.calc_agrowell_fee() - random.normalvariate(self.expense_mean, self.SES_std / 10)
        #if adopt agrowell, deduct fee from last season, and add back after 10
        if self.agrowell:
            self.agrowell_duration = self.agrowell_duration + 1


def build_field_canals(num_fcs, num_farmers_per_fc, payoff, risk_model, risk_parameters):
    num_farmers = num_farmers_per_fc * num_fcs

    if isinstance(risk_parameters, RiskParameters):
        risk_parameters = list(itertools.repeat(risk_parameters, num_farmers))
    #if risk objective as 3, randomly assign risk_objective to each agent
    if risk_model == 3:
        risk_objective = np.zeros(num_farmers)
        for i in range(num_farmers):
            risk_objective[i] = weighted_choice(((0, 0.33), (1, 0.33), (2, 0.33)))
    else:
        risk_objective = np.repeat(risk_model, num_farmers)

    fc_id = 0
    farmer_id = 0
    dfn = DecisionFunction()

    fcs = list()
    for i in range(num_fcs):
        fc_id = fc_id + 1
        farmers = list()
        for j in range(num_farmers_per_fc):
            farmer_id = farmer_id + 1
                    #assign initial SES to each farmer
            ses = random.normalvariate(Farmer.SES_mean, Farmer.SES_std)
            agrowell = ses > Farmer.SES_mean + Farmer.SES_std
            f = Farmer(farmer_id, payoff, risk_model=risk_objective[farmer_id - 1],
                       risk_parameters=risk_parameters[farmer_id - 1],
                       decision_fn=dfn,
                       ses=ses, agrowell=agrowell)
            farmers.append(f)
        fc = FieldCanal(fc_id, farmers)
        fcs.append(fc)
    return fcs


class Simulation(object):

    def __init__(self, field_canals):
        self.field_canals = field_canals
        self.number_field_canals = len(field_canals)
        self.number_farmers = np.sum([len(fc.farmers) for fc in field_canals])
        self.farmers = [f for fc in field_canals for f in fc.farmers]
        self.iter_data = None
        self.sim_data = None
        self.tank_level = None
        self.rainfall = None
        self.num_ofc = None
        self.mean_ses = None
        self.std_ses = None
        self.rf_forecast = (0.25, 0.50, 0.25)


    def update_ofc_count(self):
        self.num_ofc = np.sum([f.crop_decision == OFC_CROP for f in self.farmers])

    def update_mean_ses(self):
        self.mean_ses = np.mean([f.ses for f in self.farmers])

    def update_std_ses(self):
        self.std_ses = np.std([f.ses for f in self.farmers])

    def growing_season(self):
        self.tank_level = weighted_choice(((1, 0.25), (2, 0.50), (3, 0.25)))
        self.rf_forecast = (0.25, 0.50, 0.25)
        self.rainfall = weighted_choice(zip(range(1, 4), self.rf_forecast))


        for fc in self.field_canals:
            fc.start_growing_season(tank_level=self.tank_level, rf_forecast=self.rf_forecast)
        self.update_ofc_count()
        self.update_mean_ses()
        self.update_std_ses()
        for fc in self.field_canals:
            fc.finish_growing_season(tank_level=self.tank_level, rainfall=self.rainfall,
                                     num_ofc=self.num_ofc, mean_ses=self.mean_ses,
                                     std_ses=self.std_ses)

    def simulate(self, number_seasons):

        df = pd.DataFrame([])

        for s in range(number_seasons):

            print "Running iteration " + str(s)

            self.growing_season()

            for fc in self.field_canals:
                for f in fc.farmers:
                    agent_data = pd.Series({
                        'iter': s,
                        'rf': self.rainfall,
                        'tank': self.tank_level,
                        'farmer': f.farmer_id,
                        'fc': fc.fc_id,
                        'aw': f.agrowell,
                        'bethma': fc.bethma,
                        'risk_obj': f.risk_objective,
                        'paddy_b': f.expected_utility[0],
                        'paddy_nb': f.expected_utility[1],
                        'ofc_b': f.expected_utility[2],
                        'ofc_nb': f.expected_utility[3],
                        'paddy_b_rank': f.ranked_preferences[0],
                        'paddy_nb_rank': f.ranked_preferences[1],
                        'ofc_b_rank': f.ranked_preferences[2],
                        'ofc_nb_rank': f.ranked_preferences[3],
                        'crop': f.crop_decision,
                        'profit': f.income,
                        'ses': f.ses
                        })
                    df = df.append(agent_data, ignore_index=True)

        self.iter_data = df

        self.sim_data = pd.DataFrame([])

        for fc in self.field_canals:
            for f in fc.farmers:
                agent_data = pd.Series({
                    'aw': f.agrowell,
                    'farmer' : f.farmer_id,
                    'fc': fc.fc_id,
                    'risk_obj': f.risk_objective
                    })
                self.sim_data=self.sim_data.append(agent_data, ignore_index=True)

#functions used in code
def weighted_choice(choices):
    total = sum(w for c, w in choices)
    r = random.uniform(0, total)
    upto = 0
    for c, w in choices:
        if upto + w >= r:
            return c
        upto +=w
    assert False, "Weighted choice function error"


# :param payoff:  path to payoff table
def load_payoff(fn):
    return np.genfromtxt(str(fn), delimiter=",")

