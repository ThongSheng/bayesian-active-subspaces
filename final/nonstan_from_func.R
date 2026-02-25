library(mvtnorm)
library(BASS)
library(concordance)
library(activegp)

# --- Simulation Grid Definition ---
grid <- expand.grid(
  'p' = c(2, 10, 20),
  'n' = c(20, 100, 150),
  'seed' = 1:30,
  'func' = c('determ_full', 'determ_2d', 'GP_full', 'GP_2d'),
  'method' = c('bass', 'wycoff', 'mle')
)

# Extract parameters for this specific job array run
array_id <- as.integer(Sys.getenv('THISJOBVALUE'))
p <- grid$p[array_id]
n <- grid$n[array_id]
method <- grid$method[array_id]
set.seed(grid$seed[array_id])

# --- Data Generation  ---
x_obs <- matrix(ncol = p, nrow = n, runif(p*n))

dist_tensor <- array(dim = c(n, n, p))
for (j in 1:p) {
  dist_tensor[,,j] <- matrix(x_obs[,j], nrow = n, ncol = n) - matrix(x_obs[,j], nrow = n, ncol = n, byrow = T)
}
dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
dist_tensor_mat_reduced <- dist_tensor_mat[upper.tri(matrix(nrow =n, ncol = n), diag = T), ]

compute_cov_mat_from_C <- function(C, dist_tensor_mat_reduced, n, diagonal_add = .00001) {
  test_distances <- rowSums((dist_tensor_mat_reduced %*% C) * dist_tensor_mat_reduced)
  dist_mat_C <- matrix(nrow = n, ncol = n)
  dist_mat_C[upper.tri(dist_mat_C, diag = T)] <- test_distances
  cov_mat <- exp(-dist_mat_C/2)
  diag(cov_mat) <- diag(cov_mat) + diagonal_add
  cov_mat
}

if (grid$func[array_id] == 'determ_full') {
  f <- function(x) {sum(x^2)/sqrt(p)}
  C <- matrix(0, nrow = p, ncol = p)
  diag(C) <- 4/(3*p)
  if (p > 1) {
    C[upper.tri(C) | lower.tri(C)] <- 1/p
  }
  y <- apply(x_obs, 1, f) + rnorm(n, sd=sqrt(.00001))
  y <- y - mean(y)
  
} else if (grid$func[array_id] == 'determ_2d') {
  f <- function(x) {sum(x[1:2]^2)/sqrt(p)}
  C <- matrix(0, nrow = p, ncol = p)
  if (p >= 2) {
    C[1,1] <- C[2,2] <- 4/(3*p)
    C[1,2] <- C[2,1] <- 1/p
  }
  y <- apply(x_obs, 1, f) + rnorm(n, sd=sqrt(.00001))
  y <- y - mean(y)
  
} else if (grid$func[array_id] == 'GP_full') {
  W_random <- eigen(crossprod(matrix(rnorm(p * p), nrow = p, ncol = p)))$vectors
  C <- 500/p * (W_random %*% diag(exp((p:1 -p))) %*% t(W_random))
  Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n)
  Cov_mat[lower.tri(Cov_mat)] <- t(Cov_mat)[lower.tri(Cov_mat)]
  y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))
  
} else if (grid$func[array_id] == 'GP_2d') {
  W_random <- eigen(crossprod(matrix(rnorm(p * p), nrow = p, ncol = p)))$vectors
  C <- 500/p * (W_random %*% diag(exp(c(0, -1, rep(-Inf, p-2)))) %*% t(W_random))
  Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n)
  Cov_mat[lower.tri(Cov_mat)] <- t(Cov_mat)[lower.tri(Cov_mat)]
  y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))  
}

# --- MLE functions ---
likelihood_C <- function(C, dist_tensor_mat_reduced, n, diagonal_add, compute_gradient = F, upper_tri_n = NULL, lower_tri_n = NULL, upper_tri_p = NULL, lower_tri_p = NULL) {
  Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n, diagonal_add)
  # The following is the equivalent of mvtnorm::dmvnorm
  # but we do it manually because we only use the upper triangle of Cov_mat
  chol_Cov_mat <- chol(Cov_mat)
  det_Cov_mat <- 2*sum(log(diag(chol_Cov_mat)))
  
  #quadratic form y^T Sigma^-1 y = y^T (L L^T)^-1 y = y^T L^(-T) L^(-1) y 
  backsolved_val <- backsolve(r = chol_Cov_mat, x = y, upper.tri = T, transpose = T)
  cross_Cov_mat <- sum(backsolved_val^2)
  
  # multivariate normal log likelihood:
  -1/2 * ( n * log(2*pi) + det_Cov_mat + cross_Cov_mat)
}

make_C_from_par <- function(par, p) {
  diag_entries_index <- cumsum(1:p)
  diagonal_C <- exp(par[diag_entries_index])
  C_mat_use <- diag(p)
  C_mat_use[upper.tri(C_mat_use)] <- par[-diag_entries_index]
  C_mat_use[lower.tri(C_mat_use)] <- t(C_mat_use)[lower.tri(C_mat_use)]
  diag(sqrt(diagonal_C)) %*% C_mat_use %*% diag(sqrt(diagonal_C)) 
} 

likelihood_function <- function(par, dist_tensor_mat_reduced, n, diagonal_add) {
  C_mat_use <- make_C_from_par(par, ncol(dist_tensor_mat_reduced))
  if (sum(is.na(C_mat_use)) > 0) {
    return(10^6)
  }
  if (min(eigen(C_mat_use)$value) < 10^-6) {
    return(10^6)
  }
  ll_val <- likelihood_C(C_mat_use, dist_tensor_mat_reduced, n, diagonal_add)
  -ll_val
}

# --- Sampling ---
a_time <- Sys.time()
pred_C <- NULL

if (method == 'bass') {
  mod_bass <- BASS::bass(x_obs, y, h2=0.1, verbose=FALSE)
  pred_C <- C_bass(mod_bass)
} else if (method == 'wycoff') {
  prior_scale_matrix <- diag(x = 1, nrow = p)
  Cov_mat_use <- compute_cov_mat_from_C(prior_scale_matrix, dist_tensor_mat_reduced, n)
  Cov_mat_use[lower.tri(Cov_mat_use)] <- t(Cov_mat_use)[lower.tri(Cov_mat_use)]
  
  K_inv <- solve(Cov_mat_use)
  K_inv_y <- K_inv %*% y
  pred_C <- matrix(0, p, p)
  for (j in 1:p) {
    for (k in 1:p) {
      W_active_gp <- activegp::W_kappa_ij(x_obs, sqrt(1/diag(prior_scale_matrix)), j-1, k-1, 1)
      pred_C[j,k] <- prior_scale_matrix[j,k] - sum(K_inv * W_active_gp) + as.vector(t(K_inv_y) %*% W_active_gp %*% K_inv_y)
    }
  }
} else if (method == 'mle') {
  init_matrix <- matrix(nrow = p, ncol = p, 1)
  diag(init_matrix) <- 2/p
  init_par <- log(init_matrix[upper.tri(init_matrix, diag = T)])
  make_C_from_par(init_par, p)
  
  # lower and upper bounds
  upper_mat <-  matrix(nrow = p, ncol = p, exp(.999))
  diag(upper_mat) <- NA
  upper_par <- log(upper_mat[upper.tri(upper_mat, diag = T)])
  lower_par <- - upper_par
  
  likelihood_optim <- optim(init_par, likelihood_function,
                            dist_tensor_mat_reduced = dist_tensor_mat_reduced,
                            n = n, method = 'L-BFGS-B',
                            upper = upper_par,
                            lower = lower_par,
                            diagonal_add = .00001)
  par <- likelihood_optim$par
  pred_C <- make_C_from_par(par, p)
}

b_time <- Sys.time()
time_used <- b_time - a_time
units(time_used) <- 'mins'

# Save results
output_dir <- 'nonstan_results/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
save_file_name <- paste0(output_dir, method, "_", array_id, '.RData')
save(pred_C, C, y, time_used, file = save_file_name)
q(save = 'no')
