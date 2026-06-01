data {
  // number of sites
  int<lower = 1> n;
  int<lower = 0> ncov_psi;
  matrix[n, ncov_psi] X_psi;

  int<lower = 0> ncov_ord;
  matrix[n, ncov_ord] X_ord;

  int<lower = 0> d;

  // total number of samples
  int<lower = 1> N;
  // total number of marker replicates
  int<lower = 1> N2;
  // total number of technical replicates
  int<lower = 1> N3;

  // samples per site
  int<lower = 0> M[n];
  int<lower = 0> sumM[n];
  // sample-level detection covariates
  int<lower = 1> ncov_theta;
  matrix[N, ncov_theta] X_theta;

  // max number of marker
  int<lower = 1> maxL;

  // marker of samples
  int<lower = 1> L[N];
  int<lower = 0> sumL[N];

  // PCR per marker
  int<lower = 1> K[N2];
  int<lower = 0> sumK[N2];

  // number of species
  int<lower = 1> S;

  // survey level information
  real logy1[N3, S];
  int logy_na[N3, S];
  int y[N3, S];

  // int delta[N, S];

  // priors
  real prior_beta_psi;
  real<lower = 0> prior_beta_psi_sd;
  real prior_beta_theta;
  real<lower = 0> prior_beta_theta_sd;

  real a_sigma0;
  real b_sigma0;
  real a_sigma1;
  real b_sigma1;
  real a_p;
  real b_p;
  real a_q;
  real b_q;
  real a_theta0;
  real b_theta0;

}

parameters {

  matrix[ncov_psi, S] beta_psi;
  matrix[ncov_theta, S] beta_theta;
  matrix[ncov_ord, d] beta_ord;
  vector[S] beta0_psi;

  matrix[n, d] E;
  matrix[d, S] LL;

  real<lower = 0, upper = .4> theta0[S];
  real<lower = 0, upper = 1> p[maxL, S];
  real<lower = 0, upper = 1> q[maxL, S];

}

transformed parameters {

  matrix[n, d] U = X_ord * beta_ord + E;
  matrix[n, S] logit_psi = rep_matrix(beta0_psi', n) + X_psi * beta_psi + U * LL;
  matrix[N, S] logit_theta = X_theta * beta_theta;

  matrix[n, S] log_psi = log_inv_logit(logit_psi);
  matrix[n, S] log1m_psi = log1m_inv_logit(logit_psi);

  real log_theta0[S] = log(theta0);
  real log1m_theta0[S];
  log1m_theta0 = log1m(theta0);

  matrix[N, S] log_theta = log_inv_logit(logit_theta);
  matrix[N, S] log1m_theta = log1m_inv_logit(logit_theta);

}

model {

  real log_p_yz1;
  real log_p_yz0;
  real log_p_ydelta1;
  real log_p_ydelta0;

  for(cov in 1:ncov_ord){
    for(k in 1:d){

      beta_ord[cov, k] ~ normal(prior_beta_psi, prior_beta_psi_sd);
    }
  }

  for(s in 1:S){

    for(cov in 1:ncov_psi){
      beta_psi[cov, s] ~ normal(prior_beta_psi, prior_beta_psi_sd);
    }

    for(cov in 1:ncov_theta){
      beta_theta[cov, s] ~ normal(prior_beta_theta, prior_beta_theta_sd);
    }

    for(l in 1:maxL){
      p[l,s] ~ beta(a_p, b_p);
      q[l,s] ~ beta(a_q, b_q);
    }

    theta0[s] ~ beta(a_theta0, b_theta0);

  }

  for(s in 1:S){
      beta0_psi[s] ~ normal(0, 1);
    }

    // prior on factor loadings
  for(k in 1:d){
    for(s in 1:S){
      LL[k,s] ~ normal(0, 1);
    }
  }

  // prior on factor scores
  for(i in 1:n){
    for(k in 1:d){
      E[i,k] ~ normal(0, 1);
    }
  }

  for(s in 1:S){

    for (i in 1:n) {

      log_p_yz1 = 0;
      log_p_yz0 = 0;

      for(m in 1:M[i]){

        // log probability of y given delta = 1
        log_p_ydelta1 = 0;
        log_p_ydelta0 = 0;

         for(l in 1:L[sumM[i] + m]){

              for(k in 1:K[sumL[sumM[i] + m] + l]){

               log_p_ydelta1 +=
                  binomial_lpmf(
                y[sumK[sumL[sumM[i] + m] + l] + k,s] | 1, p[l,s]);

                  log_p_ydelta0 +=
                  binomial_lpmf(
                y[sumK[sumL[sumM[i] + m] + l] + k,s] | 1, q[l,s]);

              }

         }

          log_p_yz1 += log_sum_exp(
                  log_theta[sumM[i] + m,s] + log_p_ydelta1,
                  log1m_theta[sumM[i] + m,s] + log_p_ydelta0);

          log_p_yz0 += log_sum_exp(
                  log(theta0[s]) + log_p_ydelta1,
                  log(1 - theta0[s]) + log_p_ydelta0);

        }

        target += log_sum_exp(
        log_psi[i, s] + log_p_yz1,
        log1m_psi[i, s] + log_p_yz0
        );

      }

    }


}

generated quantities{

  vector[n] log_lik;

  real log_p_yz1;
  real log_p_yz0;
  real log_p_ydelta1;
  real log_p_ydelta0;

  // compute of loglikelihood

  for (i in 1:n) {

    log_lik[i] = 0;

    for(s in 1:S){

      log_p_yz1 = 0;
      log_p_yz0 = 0;

      for(m in 1:M[i]){

        // log probability of y given delta = 1
        log_p_ydelta1 = 0;
        log_p_ydelta0 = 0;

         for(l in 1:L[sumM[i] + m]){

              for(k in 1:K[sumL[sumM[i] + m] + l]){

               log_p_ydelta1 +=
                  binomial_lpmf(
                y[sumK[sumM[i] + m]+k,s] |
                  1, p[l,s]);

                  log_p_ydelta0 +=
                  binomial_lpmf(
                y[sumK[sumM[i] + m]+k,s] |
                  1, q[l,s]);

              }
         }

          log_p_yz1 += log_sum_exp(
                  log_theta[sumM[i] + m,s] + log_p_ydelta1,
                  log1m_theta[sumM[i] + m,s] + log_p_ydelta0);

          log_p_yz0 += log_sum_exp(
                  log(theta0[s]) + log_p_ydelta1,
                  log(1 - theta0[s]) + log_p_ydelta0);

        }

        log_lik[i] += log_sum_exp(
        log_psi[i, s] + log_p_yz1,
        log1m_psi[i, s] + log_p_yz0
        );

      }

    }


}

