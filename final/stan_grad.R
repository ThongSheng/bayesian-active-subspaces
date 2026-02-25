library(mvtnorm)
library(rstan)

# Each prior configuration now only needs:
# - stan_model_name: Path to the compiled STAN model RData file
# - default_params: A list of default prior parameters for this model
# - get_specific_data_params_func: A function to return prior-specific parameters for data_input
# - get_inits_func: A function to return initial values for STAN parameters
# - stan_model_string_var: The name of the model defined in STAN_edited_model.R

prior_configs <- list(
  "conjugate_prior" = list(
    stan_model_name = "stan/conj_prior.RData",
    stan_model_string_var = "sim.conj_prior", 
    default_params = function(p) list(
      R = diag(p)
    ),
    get_specific_data_params_func = function(p, config) {
      list(
        prior_dof = p + 5#,
        #prior_gamma_a = config$default_params$prior_gamma_a,
        #prior_gamma_b = config$default_params$prior_gamma_b
      )
    },
    get_inits_func = function(n_chains, p, config) {
      lapply(1:n_chains, function(x) {
        mat <- matrix(runif(p*p, 0, .5), nrow = p, ncol = p)
        diag(mat) <- 1
        mat <- crossprod(mat)
        list(Sigma = mat)
      })
    }
  ),
  "dirichlet_wishart_grad" = list(
    stan_model_name = "stan/dirichlet_wishart_grad_model.RData",
    stan_model_string_var = "sim.ss_dirichlet_wishart_grad", 
    default_params = list(
      prior_gamma_a = 7,
      prior_gamma_b = .05
    ),
    get_specific_data_params_func = function(p, config) {
      list(
        alpha = rep(1, p),
        prior_cor_dof = p + 5,
        prior_gamma_a = config$default_params$prior_gamma_a,
        prior_gamma_b = config$default_params$prior_gamma_b
      )
    },
    get_inits_func = function(n_chains, p, config) {
      lapply(1:n_chains, function(x) {
        mat <- matrix(runif(p*p, 0, .5), nrow = p, ncol = p)
        diag(mat) <- 1
        mat <- crossprod(mat)
        list(Q1 = mat, xi = rep(1/p, p), K = runif(1, 0.1, 10))
      })
    }
  ),
  "lognormal_inverse_wishart_grad" = list(
    stan_model_name = "stan/lniw_gp_rescale_grad_model.RData",
    stan_model_string_var = "sim.sslniw_gp_rescale_grad",
    default_params = list(
      prior_rescale_mean = 0,
      prior_rescale_var = log(2)
    ),
    get_specific_data_params_func = function(p, config) {
      list(
        prior_lgn_mean = array(rep((0), p)),
        prior_lgn_var = array(rep(log(2), p)),
        prior_dof = p + 5,
        prior_rescale_mean = config$default_params$prior_rescale_mean,
        prior_rescale_var = config$default_params$prior_rescale_var
      )
    },
    get_inits_func = function(n_chains, p, config) {
      lapply(1:n_chains, function(x) {
        q1_mat <- matrix(runif(p*p, 0, .5), nrow = p, ncol = p)
        diag(q1_mat) <- 1
        q1_mat <- crossprod(q1_mat)
        list(Q1 = q1_mat, xi = array(rnorm(p, 0, 0.1)), K = rnorm(1, 0, 0.1))
      })
    }
  ),
  "gp_reduce_grad" = list(
    stan_model_name = "stan/gp_reduce_grad_model.RData",
    stan_model_string_var = "sim.gp_reduce_grad", 
    default_params = list(
      prior_gamma_a = 7,
      prior_gamma_b = .05
    ),
    get_specific_data_params_func = function(p, config) {
      list(
        k_reduce = ifelse(p <= 5, 1, 5),
        #k_reduce = 2,
        alpha = rep(1, p),
        prior_cor_dof = p + 5,
        prior_gamma_a = config$default_params$prior_gamma_a,
        prior_gamma_b = config$default_params$prior_gamma_b
      )
    },
    get_inits_func = function(n_chains, p, config) {
      lapply(1:n_chains, function(x) {
        mat <- matrix(runif(p*p, 0, .5), nrow = p, ncol = p)
        diag(mat) <- 1
        mat <- crossprod(mat)
        list(Q1 = mat, xi = rep(1/p, p), K = runif(1, 0.1, 10))
      })
    }
  )
)

# --- Simulation Grid Definition ---
grid <- expand.grid(
  'p' = c(2, 10, 20),
  'n' = c(20, 100, 150),
  'seed' = 1:3,
  'func' = c('determ_full', 'determ_2d', 'GP_full', 'GP_2d'),
  'prior_choice' = c("conjugate_prior", "dirichlet_wishart_grad",
                     "lognormal_inverse_wishart_grad", 
                     "gp_reduce_grad")
)

# Extract parameters for this specific job array run
array_id <- as.integer(Sys.getenv('THISJOBVALUE'))
p <- grid$p[array_id]
n <- grid$n[array_id]
selected_prior_name <- as.character(grid$prior_choice[array_id])
set.seed(grid$seed[array_id])

# Get the configuration for the selected prior
current_prior_config <- prior_configs[[selected_prior_name]]

# --- Data Generation  ---
x_obs <- matrix(ncol = p, nrow = n, runif(p*n))

dist_tensor <- array(dim = c(n, n, p))
for (j in 1:p) {
  dist_tensor[,,j] <- matrix(x_obs[,j], nrow = n, ncol = n) - matrix(x_obs[,j], nrow = n, ncol = n, byrow = T)
}
dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
dist_tensor_mat_reduced <- dist_tensor_mat[upper.tri(matrix(nrow =n, ncol = n), diag = T), ]

dist_tensor_mat_reduced_crossprod <- lapply(1:n, 
                                            function(i) lapply(1:n, function(j) {
                                              tcrossprod(dist_tensor[i,j,])
                                            }))

compute_grad_cov_mat_from_C <- function(C, dist_tensor_mat, n,
                                        diagonal_add = .00001) {
  test_distances <- matrix(rowSums((dist_tensor_mat %*% C) *
                                     dist_tensor_mat), n, n)
  
  exp_mat <- exp(- test_distances/2)
  NA_mat <- matrix(NA, p, p)
  all_mats <- lapply(1:n, 
                     function(i) {lapply(1:n, function(j) {
                       if (j < i) {
                         return(NA_mat) # for the lower triangle, do not compute
                         # saves computation time
                       } else {
                         mat <- - ((C %*% dist_tensor_mat_reduced_crossprod[[i]][[j]] %*% C - C)*
                                     exp_mat[i,j])
                         if (i == j) {
                           diag(mat) <- diag(mat) + diagonal_add
                         }
                         return(mat)
                       }
                     })
                     })
  # combine all matrices into one pn times pn matrix
  cov_mat <- do.call(rbind, lapply(all_mats, function(x) do.call(cbind, x)))
}

if (grid$func[array_id] == 'determ_full') {
  f <- function(x) {2 * x/sqrt(p)}
  C <- matrix(0, nrow = p, ncol = p)
  diag(C) <- 4/(3*p)
  if (p > 1) {
    C[upper.tri(C) | lower.tri(C)] <- 1/p
  }
  y <- as.vector(apply(x_obs, 1, f))
} else if (grid$func[array_id] == 'determ_2d') {
  f <- function(x) {2 * c(x[1:2], rep(0, p-2))/sqrt(p)}
  C <- matrix(0, nrow = p, ncol = p)
  if (p >= 2) {
    C[1,1] <- C[2,2] <- 4/(3*p)
    C[1,2] <- C[2,1] <- 1/p
  }
  y <- as.vector(apply(x_obs, 1, f))
  
} else if (grid$func[array_id] == 'GP_full') {
  W_random <- eigen(crossprod(matrix(rnorm(p * p), nrow = p, ncol = p)))$vectors
  C <- 500/p * (W_random %*% diag(exp((p:1 -p))) %*% t(W_random))
  Cov_mat <- compute_grad_cov_mat_from_C(C, dist_tensor_mat, n)
  Cov_mat[lower.tri(Cov_mat)] <- t(Cov_mat)[lower.tri(Cov_mat)]
  y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))
} else if (grid$func[array_id] == 'GP_2d') {
  W_random <- eigen(crossprod(matrix(rnorm(p * p), nrow = p, ncol = p)))$vectors
  C <- 500/p * (W_random %*% diag(exp(c(0, -1, rep(-Inf, p-2)))) %*% t(W_random))
  Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n)
  Cov_mat[lower.tri(Cov_mat)] <- t(Cov_mat)[lower.tri(Cov_mat)]
  y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))  
}


# --- STAN Model Loading ---
stan_model_object <- NULL
stan_model_path <- current_prior_config$stan_model_name
stan_code_file_master <- "stan_edited_model_grad.R"
#setwd('~/active-subspace-methods/')
# Ensure the 'stan' directory exists
if (!dir.exists("stan")) {
  dir.create("stan")
}

if (!file.exists(stan_model_path)) {
  if (file.exists(stan_code_file_master)) {
    source(stan_code_file_master)
  } else {
    stop(paste0("Master STAN model code file not found: ", stan_code_file_master))
  }
  
  # Get the specific model string using its variable name
  model_string_var_name <- current_prior_config$stan_model_string_var
  if (!exists(model_string_var_name)) {
    stop(paste0("STAN model string variable '", model_string_var_name, "' not found in ", stan_code_file_master))
  }
  model_string <- get(model_string_var_name)
  
  stan_model_object <- stan_model(model_code = model_string)
  save(stan_model_object, file = stan_model_path)
} else {
  load(stan_model_path)
}


# Construct data_input dynamically
get_data_input <- function(N, k, R, y, mu0, locs, diag_add, specific_params) {
  base_data <- list(
    N = N,
    k = k,
    R = R,
    y = y,
    mu0 = mu0,
    locs = locs,
    diag_add = diag_add
  )
  # Combine base parameters with prior-specific parameters
  c(base_data, specific_params)
}
if (T) {
  mu0_use <- rep(0, p*n)
} else {
  y_mat <- matrix(y, ncol = p, byrow = T)
  y_mean <- colMeans(y_mat)
  mu0_use <- rep(y_mean, times = n)
}

# Generate prior-specific data parameters using the function from prior_configs
specific_data_params <- current_prior_config$get_specific_data_params_func(p, current_prior_config)
data_input <- get_data_input(n, 
                             p, 
                             diag(x = 1, nrow = p), 
                             y, 
                             mu0_use, 
                             x_obs, 
                             .00001, 
                             specific_data_params)

# MCMC initializations
n_chains <- 4
it <- 1500
w <- 500

initial_values <- current_prior_config$get_inits_func(n_chains, p, current_prior_config)
#sapply(initial_values, function(x) min(eigen(x$Q1)$values))


# Run STAN sampling
a_time <- Sys.time()
out_model <- sampling(object = stan_model_object, data = data_input,
                      pars = c('Sigma'), iter = it, chains = n_chains, warmup = w,
                      init = initial_values,
                      cores = n_chains)
b_time <- Sys.time()
time_used <- b_time - a_time
units(time_used) <- 'mins'

extract_vals <- rstan::extract(out_model)
summary_vals <- summary(out_model)

# Save results
output_dir <- 'stan_results/grad/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
save_file_name <- paste0(output_dir, 'grad', '_', selected_prior_name, "_", array_id, '.RData')
save(extract_vals, summary_vals, C, y, time_used, file = save_file_name)
q(save = 'no')
