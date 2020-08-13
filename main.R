load("texas.Rdata");
# N: Total populuation
# C_obs: Observed cumulative confirmed cases between 03/19/20 and 07/05/20 
# R_obs: Observed cumulative recovery cases between 03/19/20 and 07/05/20 
# D_obs: Observed cumulative death cases between 03/19/20 and 07/05/20 
# C_true: Observed cumulative confirmed cases between 07/06/20 and 07/08/20
# N_true: Observed new daily confirmed cases between 07/06/20 and 07/08/20
# We usually refer C as cumulative Confirmed, R as cumulative Recovery, D as cumulative Death, N as New daily confirmed

# Run generalized Richards growth model
T_pred <- 3; # The number of prediction days
res <- growth_model_grm(C_obs, -1.0, -1.0, N, T_pred, FALSE);
mean(abs(res$C_pred_mean - C_true)/C_true) # MAPE of cumulative confirmed: 0.036, the number may vary due to randomness
mean(abs(res$N_pred_mean - N_true)/N_true) # MAPE of new daily confirmed: 0.501, the number may vary due to randomness

# Run Richards growth model
T_pred <- 3;
res <- growth_model_grm(C_obs, 1.0, -1.0, N, T_pred, FALSE);

# Run generalized logistic growth model
T_pred <- 3;
res <- growth_model_grm(C_obs, -1.0, 1.0, N, T_pred, FALSE);

# Run logistic growth model
T_pred <- 3;
res <- growth_model_grm(C_obs, 1.0, 1.0, N, T_pred, FALSE);

# Run Bertalanffy growth model
T_pred <- 3;
res <- growth_model_grm(C_obs, 2/3, 1/3, N, T_fin, FALSE);

# Run Gompertz growth model
T_pred <- 3;
res <- growth_model_2p(C_obs, N, T_fin, 3, FALSE);
      
# Run ARIMA model
T_pred <- 3;
res <- auto.arima(diff(C_obs));
C_fit <- C_obs[1] + c(0, cumsum(res$fitted));
N_pred_mean <- pmax(0, forecast(res, h = T_pred)$mean);
C_pred_mean <- C_fit[length(C_fit)] + cumsum(N_pred_mean);
mean(abs(C_pred_mean - C_true)/C_true) # MAPE of cumulative confirmed: 0.017, the number may vary due to randomness
mean(abs(N_pred_mean - N_true)/N_true) # MAPE of new daily confirmed: 0.304, the number may vary due to randomness

# Run SIR model with both I and R compartments
T_pred <- 3;
z <- rep(0, length(C_obs)); # z is a group indicator vector to partition the entire period
I_obs <- C_obs - R_obs - D_obs;
res <- sir(I_obs, R_obs + D_obs, N, z, T_pred, FALSE);
mean(abs(res$C_pred_mean - C_true)/C_true) # MAPE of cumulative confirmed: 0.017, the number may vary due to randomness
mean(abs(res$N_pred_mean - N_true)/N_true) # MAPE of new daily confirmed: 0.248, the number may vary due to randomness
mean(res$R0_store) # R0 mean 2.446, the number may vary due to randomness
quantile(res$R0_store, c(0.025, 0.975)) # R0 95$ CI: (2.008, 2.943), the number may vary due to randomness

# Run new SIR model with only cumulative confirmed case input and fixed recovery rate
T_pred <- 3;
z <- rep(0, length(C_obs)); # z is a group indicator vector to partition the entire period
recovery_rate <- 0.1; # Please consider 0.1 and 0.15, if recovery rate = -1, the algorithm will also infer the recovery rate
res <- sir_c(C_obs, N, recovery_rate, z, T_pred, FALSE);
mean(abs(res$C_pred_mean - C_true)/C_true) # MAPE of cumulative confirmed: 0.015, the number may vary due to randomness
mean(abs(res$N_pred_mean - N_true)/N_true) # MAPE of new daily confirmed: 0.267, the number may vary due to randomness
mean(res$R0_store) # R0 mean 1.602, the number may vary due to randomness
quantile(res$R0_store, c(0.025, 0.975)) # R0 95$ CI: (1.409, 1.832), the number may vary due to randomness

# Run new SIR model with only cumulative confirmed case input and unknown recovery rate
T_pred <- 3;
z <- rep(0, length(C_obs)); # z is a group indicator vector to partition the entire period
recovery_rate <- -1; # Please consider 0.1 and 0.15, if recovery rate = -1, the algorithm will also infer the recovery rate
res <- sir_c(C_obs, N, recovery_rate, z, T_pred, FALSE);
mean(abs(res$C_pred_mean - C_true)/C_true) # MAPE of cumulative confirmed: 0.012, the number may vary due to randomness
mean(abs(res$N_pred_mean - N_true)/N_true) # MAPE of new daily confirmed: 0.206, the number may vary due to randomness
mean(res$R0_store) # R0 mean 1.720, the number may vary due to randomness
quantile(res$R0_store, c(0.025, 0.975)) # R0 95$ CI: (1.256, 2.925), the number may vary due to randomness
