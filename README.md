This repository contains the code for the ABM described in 
E.K. Burchfield and J.M. Gilligan, "Dynamics of Individual and Collective 
Agricultural Adaptation to Water Scarcity," accepted for publication
*Proceedings of the 2016 Winter Simulation Conference,"
T.M.K. Roeder, P.I. Frazier, R. Szechtman, E. Zhou, T. Huschka, and S.E. Chick, eds.

## ABSTRACT

Drought and water scarcity are growing challenges to agriculture around the world. 
Farmers can adapt through both individual and community-based collective actions. 
We draw on extensive field-work conducted with paddy farmers in rural Sri Lanka 
to study several adaptations to water scarcity, including switching to less 
water-intensive crops, farming collectively on shared land, and individually 
turning from surface water to groundwater by digging wells.  We explore how 
variability in climate affects agricultural decision-making at the community and 
individual levels using three types of decision-making, each characterized by an 
objective function:  risk-averse expected utility, regret-adjusted expected 
utility, and prospect theory loss-aversion.   We also assess how the 
introduction of individualized access to irrigation water with wells affects 
community-based drought mitigation practices.  Preliminary results suggest that 
the growth of well-irrigation may produce sudden disruptions to community-based 
adaptations, but that this depends on the mental models farmers use to think 
about risk and make decisions under uncertainty.

## INSTRUCTIONS

The main model is contained in `agrowell_abm.py`. 
Code to generate the payoff table is contained in `payoff.py`
Code to generate the figures from the paper is contained in `read_abm_data.R` 
(Fig. 2--3) and `survey_analysis/log_model_ABM_stan_no_mlm.R`.
The results of the Bayesian regression are contained in `lfit_no_mlm.Rda` 
(Fig. 1)


The command to reproduce the model runs described in the paper is:
```
python parallel_simulation.py
```
