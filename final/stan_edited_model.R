sim.sslniw_gp_rescale = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] R;
  vector[N] y;
  vector[N] mu0;
  matrix[N,k] locs;
  vector[k] prior_lgn_mean;
  vector[k] prior_lgn_var;
  real <lower=0> prior_dof;
  real prior_rescale_mean;
  real prior_rescale_var;
  real <lower=0> diag_add;
}
parameters {
  cov_matrix[k] Q1;
  vector[k] xi;
  real K;
}
transformed parameters {
  corr_matrix[k] Rho; 
  cov_matrix[k] Sigma;
  cov_matrix[N] Sigma_gp;
  vector<lower=0>[k] delta;
  vector<lower=0>[k] delta1;
  real s1;
  real s2;
  real rho;
// Rho is the correlation matrix prior, start with a Q1 ~ IW() and its transformed into
// a correlation matrix with D1*Q1*D1, wehre D1<-diag(delta1), is done with for loops

  for (i in 1:k) delta1[i] <- 1/sqrt(Q1[i,i]);
  for (n in 1:k) {
    for (m in 1:n) {
      Rho[m,n] <- delta1[m] * delta1[n] * Q1[m,n]; 
    }
  }

  for (n in 1:k) {
    for (m in (n+1):k) {
      Rho[m,n] <- Rho[n,m];
    }
  } 

// compute covariance matrix as: Sigma = K * D*Q*D, where D = diag(delta) 
  for (i in 1:k)  delta[i] <- exp( xi[i] );
  //delta[k] <- exp( 1 - sum(xi));
  for (n in 1:k) {
    for (m in 1:n) {
      Sigma[m,n] <- exp(K) * delta[m] * delta[n] * Rho[m,n]/sum(delta^2); 
    }
  }
  for (n in 1:k) {
    for (m in (n+1):k) {
      Sigma[m,n] <- Sigma[n,m];
    }
  }
  s1 <- sqrt(Sigma[1,1]);
  s2 <- sqrt(Sigma[2,2]) ;
  rho <-  Sigma[1,2] /(s1*s2);
  
  for (i in 1:N) {
    for (m in 1:i) {
      Sigma_gp[m,i] <- exp(- (locs[i,]- locs[m,]) * Sigma * (locs[i,]- locs[m,])' /2 ); 
    }
    Sigma_gp[i,i] <- Sigma_gp[i,i] + diag_add;
  }

  for (n in 1:N) {
    for (m in (n+1):N) {
      Sigma_gp[m,n] <- Sigma_gp[n,m];
    }
  } 
}
model {
  Q1 ~ inv_wishart(prior_dof, R);
  for ( i in 1:k) {
      xi[i] ~ normal(prior_lgn_mean[i], sqrt(prior_lgn_var[i]));
  }
  K ~ normal(prior_rescale_mean, sqrt(prior_rescale_var));
  y ~ multi_normal(mu0, Sigma_gp);
}
"



sim.ss_dirichlet_wishart = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] R;
  vector[N] y;
  vector[N] mu0;
  matrix[N,k] locs;
  vector[k] alpha;
  real <lower=0> prior_cor_dof;
  real <lower=0> prior_gamma_a;
  real <lower=0> prior_gamma_b;
  real <lower=0> diag_add;
}
parameters {
  cov_matrix[k] Q1;
  simplex[k] xi;
  real <lower=0> K;
}
transformed parameters {
  corr_matrix[k] Rho; 
  cov_matrix[k] Sigma;
  cov_matrix[N] Sigma_gp;
  vector<lower=0>[k] delta1;
// Rho is the correlation matrix prior, start with a Q1 ~ IW() and its transformed into
// a correlation matrix with D1*Q1*D1, wehre D1<-diag(delta1), is done with for loops

  for (i in 1:k) delta1[i] <- 1/sqrt(Q1[i,i]);
  for (n in 1:k) {
    for (m in 1:n) {
      Rho[m,n] <- delta1[m] * delta1[n] * Q1[m,n]; 
    }
  }

  for (n in 1:k) {
    for (m in (n+1):k) {
      Rho[m,n] <- Rho[n,m];
    }
  } 

// compute covariance matrix as: Sigma = K * D*Q*D, where D = diag(delta) 
  for (n in 1:k) {
    for (m in 1:n) {
      Sigma[m,n] <- K * sqrt(xi[m]) * sqrt(xi[n]) * Rho[m,n]; 
    }
  }
  for (n in 1:k) {
    for (m in (n+1):k) {
      Sigma[m,n] <- Sigma[n,m];
    }
  }
  
  for (i in 1:N) {
    for (m in 1:i) {
      Sigma_gp[m,i] <- exp(- (locs[i,]- locs[m,]) * Sigma * (locs[i,]- locs[m,])' /2 ); 
    }
    Sigma_gp[i,i] <- Sigma_gp[i,i] + diag_add;
  }

  for (n in 1:N) {
    for (m in (n+1):N) {
      Sigma_gp[m,n] <- Sigma_gp[n,m];
    }
  } 
}
model {
  Q1 ~ wishart(prior_cor_dof, R);
  xi ~ dirichlet(alpha);
  K ~ gamma(prior_gamma_a, prior_gamma_b);
  y ~ multi_normal(mu0, Sigma_gp);
}
"

sim.gp_reduce = "
data {
  int <lower=0> N;
  int <lower=0> k;
  int <lower=0> k_reduce;
  vector[N] y;
  vector[N] mu0;
  matrix[N,k] locs;
  real <lower=0> diag_add;
}
parameters {
  vector[k_reduce] log_theta_par;
  vector[k_reduce * k - k_reduce *(k_reduce - 1)/2] theta_mat;
}
transformed parameters {
  matrix[k,k_reduce] W; 
  matrix[k,k] Sigma; 
  matrix[N,k_reduce] Z;
  cov_matrix[N] Sigma_gp;

  
  
  {
    matrix[k,k] Q = diag_matrix(rep_vector(1.0, k)); 

    int ell = 0;
    for (i in 1:k_reduce) {
      int r = ell;

      ell = r + k - i + 1;
      vector[ell - r] v = theta_mat[(r+1):ell];
      vector[ell - r] u = v;
      real sgn = v[1] >= 0 ? 1.0 : -1.0;
      u[1] =  u[1] + sgn * norm2(v);
      u = u/norm2(u);
      
      matrix[ell - r , ell - r] Hhat = -sgn * (diag_matrix(rep_vector(1.0, ell - r)) - 2.0 * u * u');
      matrix[k, k] H = diag_matrix(rep_vector(1.0, k));
      int start_idx = k - (ell - r) + 1;
      H[start_idx:k, start_idx:k] = Hhat;
      Q = H * Q;
    }
    W = Q[,1:k_reduce];
  } 
  Z = locs * W;

  for (i in 1:N) {
    for (m in 1:i) {
      real dist = 0;
      for (zz in 1:k_reduce) {
        dist = dist + square(Z[i,zz] - Z[m,zz]) / square(exp(log_theta_par[zz]));
      }
      Sigma_gp[m,i] = exp(-dist/2); 
    }
    Sigma_gp[i,i] = Sigma_gp[i,i] + diag_add;
  }

  for (n in 1:N) {
    for (m in (n+1):N) {
      Sigma_gp[m,n] = Sigma_gp[n,m];
    }
  }
  
  Sigma = W * diag_matrix(1/((exp(log_theta_par))^2)) * W';
}
model {
  for ( i in 1:k_reduce) {
      log_theta_par[i] ~ normal(0, 1);
  }
  for ( i in 1:(k_reduce * k - k_reduce *(k_reduce - 1)/2)) {
      theta_mat[i] ~ normal(0, 1);
  }
  y ~ multi_normal(mu0, Sigma_gp);
}
"
