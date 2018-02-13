# statespace-inverts
State space models for inverts

Includes Bayesian versions of the Zonneveld model that can be used to characterize seasonal population dynamics for univoltine invertebrates

Simulation script:
- Assumes detection probability (per individual) is constant across all sites
- Assigns a study-wide mean and standard deviation (among sites) for:
	- total population size (lognormal distribution); sites can vary in their total population size
	- mean emergence date (normal distribution); sites can vary in their mean emergence date
	- sd of emergence date (lognormal distribution); sites can vary in the SD of their emergence date
	- mean daily survival (logit-normal distribution); sites can vary in their daily survival rate
- The simulated true underlying population dynamics include demographic stochasticity in survival and emergence
- The JAGS model does *NOT* currently account for demographic stochasticity in survival or emergence
- Does not estimate *true* population size; does not account for imperfect detection
- Even with imperfect detection, model estimates true variance in population size among sites