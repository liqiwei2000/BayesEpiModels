# Load function
Rcpp::sourceCpp('core_5.cpp')

# Input data
# US cumulative confirmed cases between 01/22 and 05/18, 2020
cumulative_cases <- c(1, 1, 2, 2, 5, 5, 5, 5, 5, 7, 8, 8, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 15, 15, 15, 51, 51, 57, 58, 60, 68, 74, 98, 118, 149, 217, 262, 402, 518, 583, 959, 1281, 1663, 2179, 2727, 3499, 4632, 6421, 7783, 13747, 19273, 25600, 33276, 43843, 53736, 65778, 83836, 101657, 121465, 140909, 161831, 188172, 213242, 243622, 275367, 308650, 336802, 366317, 397121, 428654, 462780, 496535, 526396, 555313, 580619, 607670, 636350, 667592, 699706, 732197, 758809, 784326, 811865, 840351, 869170, 905358, 938154, 965785, 988197, 1012582, 1039909, 1069424, 1103461, 1132539, 1158040, 1180375, 1204351, 1229331, 1257023, 1283929, 1309550, 1329260, 1347881, 1369376, 1390406, 1417774, 1442824, 1467820, 1486757, 1508308);
names(cumulative_cases) <- as.Date("2020-01-22") + 0:(length(cumulative_cases) - 1);
# US estimated population
N <- 328200000;

# Input settings
T_forecast <- 7; # Number of forecasting days
save <- FALSE; # Save detailed MCMC outputs



# ================================================================================
# ================================================================================
# ================================================================================
# Fit the logistic growth model
res <- growth_model_grm(cumulative_cases, 1.0, 1.0, N, T_forecast, save);

# Fit the generalized logistic growth model
res <- growth_model_grm(cumulative_cases, -1.0, 1.0, N, T_forecast, save);

# Fit the Richards growth model
res <- growth_model_grm(cumulative_cases, 1.0, -1.0, N, T_forecast, save);


# Fit the generalized Richards growth model
res <- growth_model_grm(cumulative_cases, -1.0, -1.0, N, T_forecast, save);


# Fit the Bertalanffy growth model
res <- growth_model_grm(cumulative_cases, 2/3, 1/3, N, T_forecast, save);

# Fit the Gompertz growth model
res <- growth_model_2p(cumulative_cases, N, T_forecast, 3, save);
# ================================================================================
# ================================================================================
# ================================================================================


# Plot results
library(ggplot2);
data <- data.frame(date = as.Date(names(cumulative_cases)[1]) + 0:(length(cumulative_cases) + T_forecast - 1),
                   C_obs = c(cumulative_cases, rep(NA, T_forecast)), 
                   C_fit = res$C_fit,
                   C_pred_mean = c(rep(NA, length(cumulative_cases)), res$C_pred_mean),
                   C_pred_upp = c(rep(NA, length(cumulative_cases)), res$C_pred_upp),
                   C_pred_lwr = c(rep(NA, length(cumulative_cases)), res$C_pred_lwr),
                   N_obs = c(0, diff(cumulative_cases), rep(NA, T_forecast)), 
                   N_fit = res$N_fit,
                   N_pred_mean = c(rep(NA, length(cumulative_cases)), res$N_pred_mean),
                   N_pred_upp = c(rep(NA, length(cumulative_cases)), res$N_pred_upp),
                   N_pred_lwr = c(rep(NA, length(cumulative_cases)), res$N_pred_lwr));
ggplot(data) +
  geom_point(mapping = aes(x = date, y = C_obs)) +
  geom_line(mapping = aes(x = date, y = C_fit)) +
  geom_point(mapping = aes(x = date, y = C_pred_mean), pch = 1) +
  geom_errorbar(mapping = aes(x = date, ymin = C_pred_lwr, ymax = C_pred_upp)) +
  labs(x = "Date", y = "Number of cases")



