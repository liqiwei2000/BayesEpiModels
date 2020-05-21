# BayesGrowthModels
This repository is used for accessing the performance of different population growth models in a Bayesian framework. Before running the code, please install two R packages: Rcpp and RcppArmadillo.

There are two core functions written in C++: growth_model_grm and growth_model_2p. The former can be used to fit the generalized Richards growth model (GRM) and its special cases including Richards growth model, generalized logistic growth model, logistic growth model, and von Bertalanffy growth model. For more details, refer to the paper titled "Analysis of Logistic Growth Models," written by A. Tsoularis and published in Research Letters in the Information and Mathematical Sciences (2001) volume 2, page 23-46. The latter can be used to fit the growth models that have two key parameters: carrying capacity (i.e. final epedemic size) and intrinsic growth rate (transmission rate at the early stage).

For the function growth_model_grm, the inputs are: 1) a non-decreasing count vector (C); 2) scaling of growth (p), if it is -1, then it will estimate p, while a positive value will fix the parameter; 3) curve symmetry (alpha), if it is -1, then it will estimate alpha, while a positive value will fix the parameter; 4) maximum capacity (POP); 5) the number of forecasting days (T_fin); 6) save detailed MCMC or not (store).

For the function growth_model_2p, the inputs are: 1) a non-decreasing count vector (C); 2) maximum capacity (POP); 3) number of forecasting days (T_fin); 4) choice of model: e.g. 3 for Gompeitz growth model; 5) save detailed MCMC or not (store).

For both functions, the main outputs are: 1) C_fit: cumulative case fitted curve; 2) C_upp: cumulative case fitted curve (upper bound); 3) C_lwr: cumulative case fitted curve (lower bound); 4) N_fit: new daily case fitted curve; 5) N_upp: new daily case fitted curve (upper bound); 6) N_lwr: new daily case fitted curve (lower bound).
