#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;

static double lfactorial_cpp(int i);
static int max_cpp(int a, int b);
static int min_cpp(int a, int b);
static IntegerVector sort(IntegerVector v);
static double mean_cpp(IntegerVector x);
static double lambda_builder_grm(int K, double gamma, double p, double alpha, int C);
static double lambda_builder(int model, int K, double gamma, int C);

// [[Rcpp::export]]
Rcpp::List growth_model_grm(IntegerVector C, double p, double alpha, int POP, int T_fin, bool store) {
  // Read data information
  int T = C.length();
  bool p_unknown = (p == -1.0);
  bool alpha_unknown = (alpha == -1.0);
  
  // Set algorithm settings
  int iter = 100000;
  int burn = iter*0.5;
  int error = 1; // 1 for NB variance and 0 for Poisson variance
  int K_start = (max(C) + 1)*2;
  double phi_start = 10, gamma_start = 0.5, p_start = 0.5, alpha_start = 1.0;
  double tau_logphi = 1.0, tau_logK = 0.1, tau_loggamma = 0.1, tau_logp = 0.1, tau_logalpha = 0.1;
  
  // Set hyperparameters
  double a_phi = 0.001, b_phi = 0.001; 
  int a_K = max(C) + 1, b_K = POP;
  double a_gamma, b_gamma;
  if (p == 1.0 && alpha == 1.0)
  {
    a_gamma = 1.0;
    b_gamma = 1.0;
  }
  else
  {
    a_gamma = 0.001;
    b_gamma = 0.001;
  }
  double a_p = 1, b_p = 1;
  double a_alpha = 0.001, b_alpha = 0.001;
  
  // Set temporary variables
  int it, t, count = 0, count_2 = 0;
  int K_temp, K_map;
  double phi_temp, gamma_temp, p_temp, alpha_temp, phi_map, gamma_map, p_map, alpha_map;
  double hastings, logposterior = 0, logposterior_map;
  NumericVector lambda(T - 1);
  NumericVector lambda_temp(T - 1);
  IntegerVector flag(T - 1, 0);
  double accept_phi = 0, accept_K = 0, accept_gamma = 0, accept_p = 0, accept_alpha = 0;
  NumericVector logposterior_store(iter - burn, -1.0);
  NumericVector phi_store(iter - burn, -1.0);
  IntegerVector K_store(iter - burn, -1);
  NumericVector gamma_store(iter - burn, -1.0);
  NumericVector alpha_store(iter - burn, -1.0);
  NumericVector p_store(iter - burn, -1.0);
  IntegerMatrix C_predict(200, T_fin);
  IntegerMatrix N_predict(200, T_fin);
  
  // Initialization
  int K = K_start;
  double phi = phi_start, gamma = gamma_start;
  if (p_unknown) {
    p = p_start;
  }
  if (alpha_unknown) {
    alpha = alpha_start;
  }
  phi_map = phi;
  K_map = K;
  gamma_map = gamma;
  p_map = p;
  alpha_map = p;
  IntegerVector DeltaC(T - 1);
  for (t = 0; t < T - 1; t++)
  {
    DeltaC(t) = C(t + 1) - C(t);
    if (DeltaC(t) <= 0)
    {
      flag(t) = 1;
    }
    lambda(t) = lambda_builder_grm(K, gamma, p, alpha, C(t));
    if (flag(t) == 0) 
    {
      if (error == 0)
      {
        logposterior = logposterior + DeltaC(t)*log(lambda(t)) - lambda(t);
      }
      else
      {
        logposterior = logposterior + lgamma(DeltaC(t) + phi) - lgamma(phi) + phi*(log(phi) - log(lambda(t) + phi)) + DeltaC(t)*(log(lambda(t)) - log(lambda(t) + phi));
      }
    }
  }
  if (error == 1) {
    logposterior = logposterior + (a_phi - 1)*log(phi) - b_phi*phi;
  }
  if (alpha == 1.0 && p == 1.0)
  {
    logposterior = logposterior + (a_gamma - 1)*log(gamma) - (b_gamma - 1)*log(1 - gamma);
  }
  else
  {
    logposterior = logposterior + (a_gamma - 1)*log(gamma) - b_gamma*gamma;
  }
  if (p_unknown)
  {
    logposterior = logposterior + (a_p - 1)*log(p) - (b_p - 1)*log(1 - p);
  }
  if (alpha_unknown)
  {
    logposterior = logposterior + (a_alpha - 1)*log(alpha) - b_alpha*alpha;
  }
  if (error == 0)
  {
    phi = -1;
  }
  
  // MCMC
  for (it = 0; it < iter; it++)
  {
    // Update phi
    if (error == 1)
    {
      phi_temp = exp(r_truncnorm(log(phi), tau_logphi, log(1), log(100)));
      hastings = 0;
      for(t = 0; t < T - 1; t++)
      {
        if (flag(t) == 0) 
        {
          hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp + DeltaC(t)) - (phi_temp + DeltaC(t))*log(phi_temp + lambda(t));
          hastings = hastings - (phi*log(phi) - lgamma(phi) + lgamma(phi + DeltaC(t)) - (phi + DeltaC(t))*log(phi + lambda(t)));
        }
      }
      hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
      hastings = hastings - ((a_phi - 1)*log(phi) - b_phi*phi);
      if(hastings >= log(double(rand()%10001)/10000))
      {
        phi = phi_temp;
        logposterior = logposterior + hastings;
        if (it > burn) {
          accept_phi++;
        }
      }
    }
    
    // Update K
    K_temp = exp(r_truncnorm(log(K), tau_logK, log(a_K), log(b_K)));
    for(t = 0; t < T - 1; t++)
    {
      lambda_temp(t) = lambda_builder_grm(K_temp, gamma, p, alpha, C(t));
    }
    hastings = 0;
    for(t = 0; t < T - 1; t++)
    {
      if (flag(t) == 0) 
      {
        if (error == 0)
        {
          hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - lambda_temp(t));
          hastings = hastings - (DeltaC(t)*log(lambda(t)) - lambda(t));
        }
        else
        {
          hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - (phi + DeltaC(t))*log(phi + lambda_temp(t)));
          hastings = hastings - (DeltaC(t)*log(lambda(t)) - (phi + DeltaC(t))*log(phi + lambda(t)));
        }
      }
    }
    if(hastings >= log(double(rand()%10001)/10000))
    {
      K = K_temp;
      logposterior = logposterior + hastings;
      for(t = 0; t < T - 1; t++)
      {
        lambda(t) = lambda_temp(t);
      }
      if (it > burn) {
        accept_K++;
      }
    }
    
    // Update gamma
    if (p == 1.0 && alpha == 1.0) 
    {
      gamma_temp = exp(r_truncnorm(log(gamma), tau_loggamma, log(10e-9), log(1)));
    }
    else
    {
      gamma_temp = exp(rnorm(1, log(gamma), tau_loggamma)(0));
    }
    for(t = 0; t < T - 1; t++)
    {
      lambda_temp(t) = lambda_builder_grm(K, gamma_temp, p, alpha, C(t));
    }
    hastings = 0;
    for(t = 0; t < T - 1; t++)
    {
      if (flag(t) == 0) 
      {
        if (error == 0)
        {
          hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - lambda_temp(t));
          hastings = hastings - (DeltaC(t)*log(lambda(t)) - lambda(t));
        }
        else
        {
          hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - (phi + DeltaC(t))*log(phi + lambda_temp(t)));
          hastings = hastings - (DeltaC(t)*log(lambda(t)) - (phi + DeltaC(t))*log(phi + lambda(t)));
        }
      }
    }
    if (p == 1.0 && alpha == 1.0)
    {
      hastings = hastings + (a_gamma - 1)*log(gamma_temp) - (b_gamma - 1)*log(1 - gamma_temp);
      hastings = hastings - ((a_gamma - 1)*log(gamma) - (b_gamma - 1)*log(1 - gamma));
    }
    else
    {
      hastings = hastings + (a_gamma - 1)*log(gamma_temp) - b_gamma*gamma_temp;
      hastings = hastings - ((a_gamma - 1)*log(gamma) - b_gamma*gamma);
    }
    if(hastings >= log(double(rand()%10001)/10000))
    {
      gamma = gamma_temp;
      logposterior = logposterior + hastings;
      for(t = 0; t < T - 1; t++)
      {
        lambda(t) = lambda_temp(t);
      }
      if (it > burn) {
        accept_gamma++;
      }
    }
    
    // Update p
    if (p_unknown)
    {
      p_temp = exp(r_truncnorm(log(p), tau_logp, log(10e-9), log(1)));
      for(t = 0; t < T - 1; t++)
      {
        lambda_temp(t) = lambda_builder_grm(K, gamma, p_temp, alpha, C(t));
      }
      hastings = 0;
      for(t = 0; t < T - 1; t++)
      {
        if (flag(t) == 0) 
        {
          if (error == 0)
          {
            hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - lambda_temp(t));
            hastings = hastings - (DeltaC(t)*log(lambda(t)) - lambda(t));
          }
          else
          {
            hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - (phi + DeltaC(t))*log(phi + lambda_temp(t)));
            hastings = hastings - (DeltaC(t)*log(lambda(t)) - (phi + DeltaC(t))*log(phi + lambda(t)));
          }    
          
        }
      }
      hastings = hastings + (a_p - 1)*log(p_temp) - (b_p - 1)*log(1 - p_temp);
      hastings = hastings - ((a_p - 1)*log(p) - (b_p - 1)*log(1 - p));
      if(hastings >= log(double(rand()%10001)/10000))
      {
        p = p_temp;
        logposterior = logposterior + hastings;
        for(t = 0; t < T - 1; t++)
        {
          lambda(t) = lambda_temp(t);
        }
        if (it > burn) {
          accept_p++;
        }
      }
    }
    
    // Update alpha
    if (alpha_unknown)
    {
      alpha_temp = exp(rnorm(1, log(alpha), tau_logalpha)(0));
      for(t = 0; t < T - 1; t++)
      {
        lambda_temp(t) = lambda_builder_grm(K, gamma, p, alpha_temp, C(t));
      }
      hastings = 0;
      for(t = 0; t < T - 1; t++)
      {
        if (flag(t) == 0) 
        {
          if (error == 0)
          {
            hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - lambda_temp(t));
            hastings = hastings - (DeltaC(t)*log(lambda(t)) - lambda(t));
          }
          else
          {
            hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - (phi + DeltaC(t))*log(phi + lambda_temp(t)));
            hastings = hastings - (DeltaC(t)*log(lambda(t)) - (phi + DeltaC(t))*log(phi + lambda(t)));
          }
        }
      }
      hastings = hastings + (a_alpha - 1)*log(alpha_temp) - b_alpha*alpha_temp;
      hastings = hastings - ((a_alpha - 1)*log(alpha) - b_alpha*alpha);
      if(hastings >= log(double(rand()%10001)/10000))
      {
        alpha = alpha_temp;
        logposterior = logposterior + hastings;
        for(t = 0; t < T - 1; t++)
        {
          lambda(t) = lambda_temp(t);
        }
        if (it > burn) {
          accept_alpha++;
        }
      }
    }
    
    // Monitor the process
    if(it*100/iter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    if (it >= burn)
    {
      if (it == burn)
      {
        logposterior_map = logposterior;
        phi_map = phi;
        K_map = K;
        gamma_map = gamma;
        p_map = p;
        alpha_map = alpha;
      }
      else
      {
        if(logposterior > logposterior_map)
        {
          logposterior_map = logposterior;
          phi_map = phi;
          K_map = K;
          gamma_map = gamma;
          p_map = p;
          alpha_map = alpha;
        }
      }
      if (store)
      {
        logposterior_store(it - burn) = logposterior;
        phi_store(it - burn) = phi;
        K_store(it - burn) = K;
        gamma_store(it - burn) = gamma;
        p_store(it - burn) = p;
        alpha_store(it - burn) = alpha;
      }
      // Prediction
      if (it % 250 == 0)
      {
        for (t = 0; t < T_fin; t++) 
        {
          if (t == 0)
          {
            if (error == 1)
            {
              N_predict(count_2, t) = rnbinom_mu(1, phi, lambda_builder_grm(K, gamma, p, alpha, C(T - 1)))(0);
              C_predict(count_2, t) = C(T - 1) + N_predict(count_2, t);
            }
            else
            {
              N_predict(count_2, t) = rpois(1, lambda_builder_grm(K, gamma, p, alpha, C(T - 1)))(0);
              C_predict(count_2, t) = C(T - 1) + N_predict(count_2, t);
            }
          }
          else
          {
            if (error == 1)
            {
              N_predict(count_2, t) = rnbinom_mu(1, phi, lambda_builder_grm(K, gamma, p, alpha, C_predict(count_2, t - 1)))(0);
              C_predict(count_2, t) = C_predict(count_2, t - 1) + N_predict(count_2, t);
            }
            else
            {
              N_predict(count_2, t) = rpois(1, lambda_builder_grm(K, gamma, p, alpha, C_predict(count_2, t - 1)))(0);
              C_predict(count_2, t) = C_predict(count_2, t - 1) + N_predict(count_2, t);
            }
          }
        }
        count_2 = count_2 + 1;
      }
    }
  }
  
  // Update C0
  int c, C0 = 1;
  double lambda_0, logprob_max, logprob;
  for (c = 1; c < C(0); c++)
  {
    lambda_0 = lambda_builder_grm(K_map, gamma_map, p_map, alpha_map, c);
    if (c == 1)
    {
      if (error == 0)
      {
        logprob_max = (C(0) - C0)*log(lambda_0) - lambda_0 - lfactorial_cpp(C(0) - C0);
      }
      else
      {
        logprob_max = lgamma(C(0) - C0 + phi_map) - lfactorial_cpp(C(0) - C0) - phi_map*log(lambda_0 + phi_map) + (C(0) - C0)*(log(lambda_0) - log(lambda_0 + phi_map));
      }
    }
    else
    {
      if (error == 0)
      {
        logprob = (C(0) - C0)*log(lambda_0) -lambda_0 - lfactorial_cpp(C(0) - C0);
      }
      if (error == 1)
      {
        logprob = lgamma(C(0) - C0 + phi_map) - lfactorial_cpp(C(0) - C0) - phi_map*log(lambda_0 + phi_map) + (C(0) - C0)*(log(lambda_0) - log(lambda_0 + phi_map));
      }
      if (logprob > logprob_max)
      {
        logprob_max = logprob;
        C0 = c;
      }
    }
  }
  
  // Fit the MAP model
  // int peak; //= K_map*pow(1.0*p_map/(alpha_map + p_map), 1/alpha_map);
  NumericVector C_fit(T + T_fin, 1.0*POP);
  NumericVector N_fit(T + T_fin, 0.0);
  for (t = 0; t < T + T_fin; t++)
  {
    if (t == 0)
    {
      N_fit(t) = lambda_builder_grm(K_map, gamma_map, p_map, alpha_map, C0);
      C_fit(t) = C0 + N_fit(t);
    }
    else
    {
      N_fit(t) = lambda_builder_grm(K_map, gamma_map, p_map, alpha_map, C_fit(t - 1));
      if (C_fit(t - 1) + N_fit(t) >= POP)
      {
        N_fit(t) = POP - C_fit(t - 1);
        break;
      }
      C_fit(t) = C_fit(t - 1) + N_fit(t);
    }
  }
  
  // Prediction
  NumericVector C_mean(T_fin);
  IntegerVector C_upp(T_fin);
  IntegerVector C_lwr(T_fin);
  NumericVector N_mean(T_fin);
  IntegerVector N_upp(T_fin);
  IntegerVector N_lwr(T_fin);
  IntegerVector vtemp;
  for (t = 0; t < T_fin; t++) 
  {
    vtemp = sort(C_predict.column(t));
    C_lwr(t) = vtemp(4);
    C_upp(t) = vtemp(194);
    C_mean(t) = mean_cpp(vtemp);
    vtemp = sort(N_predict.column(t));
    N_lwr(t) = vtemp(4);
    N_upp(t) = vtemp(194);
    N_mean(t) = mean_cpp(vtemp);
  }
  
  // Result wrap-up
  NumericVector accept(5);
  accept(0) = accept_phi/(iter - burn);
  accept(1) = accept_K/(iter - burn);
  accept(2) = accept_gamma/(iter - burn);
  accept(3) = accept_p/(iter - burn);
  accept(4) = accept_alpha/(iter - burn);
  
  if (store)
  {
    return Rcpp::List::create(Rcpp::Named("C_pred") = C_predict, 
                              Rcpp::Named("N_pred") = N_predict, 
                              Rcpp::Named("C_fit") = C_fit, 
                              Rcpp::Named("N_fit") = N_fit, 
                              Rcpp::Named("C0") = C0, 
                              Rcpp::Named("accept") = accept, 
                              Rcpp::Named("iter") = iter, 
                              Rcpp::Named("phi_map") = phi_map, 
                              Rcpp::Named("K_map") = K_map, 
                              Rcpp::Named("gamma_map") = gamma_map, 
                              Rcpp::Named("p_map") = p_map, 
                              Rcpp::Named("alpha_map") = alpha_map, 
                              Rcpp::Named("logposterior_store") = logposterior_store,
                              Rcpp::Named("phi_store") = phi_store, 
                              Rcpp::Named("K_store") = K_store, 
                              Rcpp::Named("gamma_store") = gamma_store, 
                              Rcpp::Named("p_store") = p_store,
                              Rcpp::Named("alpha_store") = alpha_store);
  }
  else
  {
    return Rcpp::List::create(Rcpp::Named("C_pred_mean") = C_mean, 
                              Rcpp::Named("C_pred_upp") = C_upp, 
                              Rcpp::Named("C_pred_lwr") = C_lwr, 
                              Rcpp::Named("N_pred_mean") = N_mean, 
                              Rcpp::Named("N_pred_upp") = N_upp, 
                              Rcpp::Named("N_pred_lwr") = N_lwr, 
                              Rcpp::Named("C_fit") = C_fit, 
                              Rcpp::Named("N_fit") = N_fit, 
                              Rcpp::Named("C0") = C0, 
                              Rcpp::Named("accept") = accept, 
                              Rcpp::Named("iter") = iter, 
                              Rcpp::Named("phi_map") = phi_map, 
                              Rcpp::Named("K_map") = K_map, 
                              Rcpp::Named("gamma_map") = gamma_map, 
                              Rcpp::Named("p_map") = p_map, 
                              Rcpp::Named("alpha_map") = alpha_map);
  }
}

// [[Rcpp::export]]
Rcpp::List growth_model_2p(IntegerVector C, int POP, int T_fin, int model, bool store) {
  // Read data information
  int T = C.length();
  
  // Set algorithm settings
  int iter = 100000;
  int burn = iter*0.5; 
  int error = 1; // 1 for NB variance and 0 for Poisson variance
  int K_start = (max(C) + 1)*2;
  double phi_start = 10, gamma_start = 0.5, p_start = 1, alpha_start = 1.0;
  double tau_logphi = 1.0, tau_logK = 0.1, tau_loggamma = 0.1;
  
  // Set hyperparameters
  double a_phi = 0.001, b_phi = 0.001; 
  int a_K = max(C) + 1, b_K = POP;
  double a_gamma, b_gamma;
  if (model == 0 || model == 1 || model == 3)
  {
    a_gamma = 1;
    b_gamma = 1;
  }
  else
  {
    a_gamma = 0.001;
    b_gamma = 0.001;
  }

  // Set temporary variables
  int it, t, count = 0, count_2 = 0;
  int K_temp, K_map;
  double phi_temp, gamma_temp, phi_map, gamma_map;
  double hastings, logposterior = 0, logposterior_map;
  NumericVector lambda(T - 1);
  NumericVector lambda_temp(T - 1);
  IntegerVector flag(T - 1, 0);
  double accept_phi = 0, accept_K = 0, accept_gamma = 0;
  NumericVector logposterior_store(iter - burn, -1.0);
  NumericVector phi_store(iter - burn, -1.0);
  IntegerVector K_store(iter - burn, -1);
  NumericVector gamma_store(iter - burn, -1.0);
  IntegerMatrix C_predict(200, T_fin);
  IntegerMatrix N_predict(200, T_fin);
  
  // Initialization
  int K = K_start;
  double phi = phi_start, gamma = gamma_start;
  IntegerVector DeltaC(T - 1);
  phi_map = phi;
  K_map = K;
  gamma_map = gamma;
  for (t = 0; t < T - 1; t++)
  {
    DeltaC(t) = C(t + 1) - C(t);
    if (DeltaC(t) <= 0)
    {
      flag(t) = 1;
    }
    lambda(t) = lambda_builder(model, K, gamma, C(t));
    if (flag(t) == 0) 
    {
      if (error == 0)
      {
        logposterior = logposterior + DeltaC(t)*log(lambda(t)) - lambda(t);
      }
      else
      {
        logposterior = logposterior + lgamma(DeltaC(t) + phi) - lgamma(phi) + phi*(log(phi) - log(lambda(t) + phi)) + DeltaC(t)*(log(lambda(t)) - log(lambda(t) + phi));
      }
    }
  }
  if (error == 1) {
    logposterior = logposterior + (a_phi - 1)*log(phi) - b_phi*phi;
  }
  if (model == 0 || model == 1 || model == 3)
  {
    logposterior = logposterior + (a_gamma - 1)*log(gamma) - (b_gamma - 1)*log(1 - gamma);
  }
  else
  {
    logposterior = logposterior + (a_gamma - 1)*log(gamma) - b_gamma*gamma;
  }
  if (error == 0)
  {
    phi = -1;
  }
  
  // MCMC
  for (it = 0; it < iter; it++)
  {
    // Update phi
    if (error == 1)
    {
      phi_temp = exp(r_truncnorm(log(phi), tau_logphi, log(1), log(100)));
      hastings = 0;
      for(t = 0; t < T - 1; t++)
      {
        if (flag(t) == 0) 
        {
          hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp + DeltaC(t)) - (phi_temp + DeltaC(t))*log(phi_temp + lambda(t));
          hastings = hastings - (phi*log(phi) - lgamma(phi) + lgamma(phi + DeltaC(t)) - (phi + DeltaC(t))*log(phi + lambda(t)));
        }
      }
      hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
      hastings = hastings - ((a_phi - 1)*log(phi) - b_phi*phi);
      if(hastings >= log(double(rand()%10001)/10000))
      {
        phi = phi_temp;
        logposterior = logposterior + hastings;
        if (it > burn) {
          accept_phi++;
        }
      }
    }
    
    // Update K
    K_temp = exp(r_truncnorm(log(K), tau_logK, log(a_K), log(b_K)));
    for(t = 0; t < T - 1; t++)
    {
      lambda_temp(t) = lambda_builder(model, K_temp, gamma, C(t));
    }
    hastings = 0;
    for(t = 0; t < T - 1; t++)
    {
      if (flag(t) == 0) 
      {
        if (error == 0)
        {
          hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - lambda_temp(t));
          hastings = hastings - (DeltaC(t)*log(lambda(t)) - lambda(t));
        }
        else
        {
          hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - (phi + DeltaC(t))*log(phi + lambda_temp(t)));
          hastings = hastings - (DeltaC(t)*log(lambda(t)) - (phi + DeltaC(t))*log(phi + lambda(t)));
        }
      }
    }
    if(hastings >= log(double(rand()%10001)/10000))
    {
      K = K_temp;
      logposterior = logposterior + hastings;
      for(t = 0; t < T - 1; t++)
      {
        lambda(t) = lambda_temp(t);
      }
      if (it > burn) {
        accept_K++;
      }
    }
  
    // Update gamma
    if (model == 0 || model == 1 || model == 3)
    {
      gamma_temp = exp(r_truncnorm(log(gamma), tau_loggamma, log(10e-9), log(1)));
    }
    else
    {
      gamma_temp = exp(rnorm(1, log(gamma), tau_loggamma)(0));
    }
    for(t = 0; t < T - 1; t++)
    {
      lambda_temp(t) = lambda_builder(model, K, gamma_temp, C(t));
    }
    hastings = 0;
    for(t = 0; t < T - 1; t++)
    {
      if (flag(t) == 0) 
      {
        if (error == 0)
        {
          hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - lambda_temp(t));
          hastings = hastings - (DeltaC(t)*log(lambda(t)) - lambda(t));
        }
        else
        {
          hastings = hastings + (DeltaC(t)*log(lambda_temp(t)) - (phi + DeltaC(t))*log(phi + lambda_temp(t)));
          hastings = hastings - (DeltaC(t)*log(lambda(t)) - (phi + DeltaC(t))*log(phi + lambda(t)));
        }
      }
    }
    if (model == 0 || model == 1 || model == 3)
    {
      hastings = hastings + (a_gamma - 1)*log(gamma_temp) - (b_gamma - 1)*log(1 - gamma_temp);
      hastings = hastings - ((a_gamma - 1)*log(gamma) - (b_gamma - 1)*log(1 - gamma));
    }
    else
    {
      hastings = hastings + (a_gamma - 1)*log(gamma_temp) - b_gamma*gamma_temp;
      hastings = hastings - ((a_gamma - 1)*log(gamma) - b_gamma*gamma);
    }
    if(hastings >= log(double(rand()%10001)/10000))
    {
      gamma = gamma_temp;
      logposterior = logposterior + hastings;
      for(t = 0; t < T - 1; t++)
      {
        lambda(t) = lambda_temp(t);
      }
      if (it > burn) {
        accept_gamma++;
      }
    }
    
    // Monitor the process
    if(it*100/iter == count)
    {
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    if (it >= burn)
    {
      if (it == burn)
      {
        logposterior_map = logposterior;
        phi_map = phi;
        K_map = K;
        gamma_map = gamma;
      }
      else
      {
        if(logposterior > logposterior_map)
        {
          logposterior_map = logposterior;
          phi_map = phi;
          K_map = K;
          gamma_map = gamma;
        }
      }
      if (store)
      {
        logposterior_store(it - burn) = logposterior;
        phi_store(it - burn) = phi;
        K_store(it - burn) = K;
        gamma_store(it - burn) = gamma;
      }
      //Prediction
      if (it % 250 == 0)
      {
        for (t = 0; t < T_fin; t++) {
          if (t == 0)
          {
            if (error == 1)
            {
              N_predict(count_2, t) = rnbinom_mu(1, phi, lambda_builder(model, K, gamma, C(T - 1)))(0);
              C_predict(count_2, t) = C(T - 1) + N_predict(count_2, t);
            }
            else
            {
              N_predict(count_2, t) = rpois(1, lambda_builder(model, K, gamma, C(T - 1)))(0);
              C_predict(count_2, t) = C(T - 1) + N_predict(count_2, t);
            }
          }
          else
          {
            if (error == 1)
            {
              N_predict(count_2, t) = rnbinom_mu(1, phi, lambda_builder(model, K, gamma, C_predict(count_2, t - 1)))(0);
              C_predict(count_2, t) = C_predict(count_2, t - 1) + N_predict(count_2, t);
            }
            else
            {
              N_predict(count_2, t) = rpois(1, lambda_builder(model, K, gamma, C_predict(count_2, t - 1)))(0);
              C_predict(count_2, t) = C_predict(count_2, t - 1) + N_predict(count_2, t);
            }
          }
        }
        count_2 = count_2 + 1;
      }
    }
  }
  
  // Update C0
  int c, C0 = 1;
  double lambda_0, logprob_max, logprob;
  for (c = 1; c < C(0); c++)
  {
    lambda_0 = lambda_builder(model, K_map, gamma_map, c);
    if (c == 1)
    {
      if (error == 0)
      {
        logprob_max = (C(0) - C0)*log(lambda_0) - lambda_0 - lfactorial_cpp(C(0) - C0);
      }
      else
      {
        logprob_max = lgamma(C(0) - C0 + phi_map) - lfactorial_cpp(C(0) - C0) - phi_map*log(lambda_0 + phi_map) + (C(0) - C0)*(log(lambda_0) - log(lambda_0 + phi_map));
      }
    }
    else
    {
      if (error == 0)
      {
        logprob = (C(0) - C0)*log(lambda_0) -lambda_0 - lfactorial_cpp(C(0) - C0);
      }
      if (error == 1)
      {
        logprob = lgamma(C(0) - C0 + phi_map) - lfactorial_cpp(C(0) - C0) - phi_map*log(lambda_0 + phi_map) + (C(0) - C0)*(log(lambda_0) - log(lambda_0 + phi_map));
      }
      if (logprob > logprob_max)
      {
        logprob_max = logprob;
        C0 = c;
      }
    }
  }
  
  // Fit the MAP model
  NumericVector C_fit(T + T_fin, 1.0*POP);
  NumericVector N_fit(T + T_fin, 0.0);
  for (t = 0; t < T + T_fin; t++)
  {
    if (t == 0)
    {
      N_fit(t) = lambda_builder(model, K_map, gamma_map, C0);
      C_fit(t) = C0 + N_fit(t);
    }
    else
    {
      N_fit(t) = lambda_builder(model, K_map, gamma_map, C_fit(t - 1));
      if (C_fit(t - 1)+ N_fit(t) >= POP)
      {
        N_fit(t) = POP - C_fit(t - 1);
        break;
      }
      C_fit(t) = C_fit(t - 1)+ N_fit(t);
    }
  }
  
  // Prediction
  NumericVector C_mean(T_fin);
  IntegerVector C_upp(T_fin);
  IntegerVector C_lwr(T_fin);
  NumericVector N_mean(T_fin);
  IntegerVector N_upp(T_fin);
  IntegerVector N_lwr(T_fin);
  IntegerVector vtemp;
  for (t = 0; t < T_fin; t++) 
  {
    vtemp = sort(C_predict.column(t));
    C_lwr(t) = vtemp(4);
    C_upp(t) = vtemp(194);
    C_mean(t) = mean_cpp(vtemp);
    vtemp = sort(N_predict.column(t));
    N_lwr(t) = vtemp(4);
    N_upp(t) = vtemp(194);
    N_mean(t) = mean_cpp(vtemp);
  }

  // Result wrap-up
  NumericVector accept(3);
  accept(0) = accept_phi/(iter - burn);
  accept(1) = accept_K/(iter - burn);
  accept(2) = accept_gamma/(iter - burn);
  
  if (store)
  {
    return Rcpp::List::create(Rcpp::Named("C_pred") = C_predict, 
                              Rcpp::Named("C_fit") = C_fit, 
                              Rcpp::Named("N_pred") = N_predict, 
                              Rcpp::Named("N_fit") = N_fit, 
                              Rcpp::Named("C0") = C0, 
                              Rcpp::Named("iter") = iter, 
                              Rcpp::Named("accept") = accept, 
                              Rcpp::Named("phi_map") = phi_map, 
                              Rcpp::Named("K_map") = K_map, 
                              Rcpp::Named("gamma_map") = gamma_map, 
                              Rcpp::Named("logposterior_store") = logposterior_store,
                              Rcpp::Named("phi_store") = phi_store, 
                              Rcpp::Named("K_store") = K_store, 
                              Rcpp::Named("gamma_store") = gamma_store);
  }
  else
  {
    return Rcpp::List::create(Rcpp::Named("C_pred_mean") = C_mean, 
                              Rcpp::Named("C_pred_upp") = C_upp, 
                              Rcpp::Named("C_pred_lwr") = C_lwr, 
                              Rcpp::Named("N_pred_mean") = N_mean, 
                              Rcpp::Named("N_pred_upp") = N_upp, 
                              Rcpp::Named("N_pred_lwr") = N_lwr,
                              Rcpp::Named("C_fit") = C_fit, 
                              Rcpp::Named("N_fit") = N_fit, 
                              Rcpp::Named("accept") = accept, 
                              Rcpp::Named("C0") = C0, 
                              Rcpp::Named("iter") = iter, 
                              Rcpp::Named("phi_map") = phi_map, 
                              Rcpp::Named("K_map") = K_map, 
                              Rcpp::Named("gamma_map") = gamma_map);
  }
}

// [[Rcpp::export]]
double lambda_builder(int model, int K, double gamma, int C) {
  double lambda;
  switch(model) {
  case 0: // Logistic model
    lambda = gamma*C*(1.0 - 1.0*C/K);
    break;
  case 1: // Bertalanffy model
    lambda = gamma*(K - C);
    break;
  case 2: // Bertalanffy model 2
    lambda = gamma*pow(C, 2.0/3.0)*(1 - pow(1.0*C/K, 1.0/3.0));
    break;
  case 3: // Gompeitz model
    lambda = gamma*C*log(1.0*K/C);
    break;
  }
  if (lambda < 0.0) {
    lambda = 0.0;
  }
  return lambda;
}

// [[Rcpp::export]]
double lambda_builder_grm(int K, double gamma, double p, double alpha, int C) {
  
  double lambda = gamma*pow(C, p)*(1.0 - 1.0*pow(1.0*C/K, alpha));
  if (lambda < 0.0) {
    lambda = 0.0;
  }
  return lambda;
}

// [[Rcpp::export]]
double lambda_builder_slr(double R0, double gamma, int NN, int C, int N) {
  double lambda = exp(gamma*(R0*(1.0 - 1.0*C/NN) - 1))*N;
  return lambda;
}

double lfactorial_cpp(int i) {
  double temp = 0;
  for (int j = 0; j < i; j++)
  {
    temp = temp + log(1.0 + j);
  }
  return temp;
}

int max_cpp(int a, int b) {
  if (b > a) 
  {
    return (b);
  } 
  else
  {
    return (a);
  }
}

int min_cpp(int a, int b) {
  if (b <= a) 
  {
    return (b);
  } 
  else
  {
    return (a);
  }
}

IntegerVector sort(IntegerVector v) {
  std::sort(v.begin(), v.end());
  return v;
}

double mean_cpp(IntegerVector x) {
  int n = x.size(); // Size of vector
  double sum = 0; // Sum value
  // For loop, note cpp index shift to 0
  for(int i = 0; i < n; i++){
    // Shorthand for sum = sum + x[i]
    sum = sum + x(i);
  }
  return sum/n; // Obtain and return the Mean
}
