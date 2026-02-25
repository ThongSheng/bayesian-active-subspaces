load('~/active-subspace-methods/CESM_data.RData')
load('~/active-subspace-methods/CESM_vars.RData')

library(mvtnorm)
library(BASS)
library(concordance)
library(activegp)
n <- length(RESTOM_vals)


compute_cov_mat_from_C <- function(C, dist_tensor_mat_reduced, n, diagonal_add = .00001,
                                   sigma2) {
  test_distances <- rowSums((dist_tensor_mat_reduced %*% C) * dist_tensor_mat_reduced)
  dist_mat_C <- matrix(nrow = n, ncol = n)
  dist_mat_C[upper.tri(dist_mat_C, diag = T)] <- test_distances
  cov_mat <- sigma2 * exp(-dist_mat_C/2)
  diag(cov_mat) <- diag(cov_mat) + diagonal_add
  cov_mat
}

# --- MLE functions ---
likelihood_C <- function(C, dist_tensor_mat_reduced, n, diagonal_add, sigma2, compute_gradient = F, 
                         upper_tri_n = NULL, lower_tri_n = NULL, upper_tri_p = NULL, lower_tri_p = NULL,
                         y_use) {
  Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n, diagonal_add, sigma2)
  # The following is the equivalent of mvtnorm::dmvnorm
  # but we do it manually because we only use the upper triangle of Cov_mat
  chol_Cov_mat <- chol(Cov_mat)
  det_Cov_mat <- 2*sum(log(diag(chol_Cov_mat)))
  
  #quadratic form y^T Sigma^-1 y = y^T (L L^T)^-1 y = y^T L^(-T) L^(-1) y 
  backsolved_val <- backsolve(r = chol_Cov_mat, x = y_use, upper.tri = T, transpose = T)
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

likelihood_function <- function(par, dist_tensor_mat_reduced, n,y_use) {
  C_mat_use <- make_C_from_par(par[1:(length(par)-2)], ncol(dist_tensor_mat_reduced))
  diagonal_add <- exp(par[length(par)-1])
  sigma2 <- exp(par[length(par)])
  print(C_mat_use[1:5, 1:5])
  if (sum(is.na(C_mat_use)) > 0) {
    return(10^6)
  }
  if (min(eigen(C_mat_use)$value) < 10^-6) {
    return(10^6)
  }
  ll_val <- likelihood_C(C_mat_use, dist_tensor_mat_reduced, n, diagonal_add, y_use = y_use,
                         sigma2 = sigma2)
  print(ll_val)
  -ll_val
}
y <- RESTOM_vector
x_obs <- variables_df[, colnames(variables_df) != c('Sample_nmb')]
x_obs <- x_obs[-nrow(x_obs),]
#x_obs <- x_obs[-1,]
p <- ncol(x_obs)

renormalize <- function(x, bounds = c(range(x))) {
  (x - bounds[1])/(bounds[2] - bounds[1])
}
bounds_list <- list('clubb_c1' = c(.4, 3),
                    'clubb_gamma_coef' = c(.25, .35),
                    'zmconv_dmpdz' = c(-.002, -0.0002), 
                    'clubb_C2rt' = c(.2, 2),
                    'micro_mg_autocon_nd_exp' = c(-2, -.8),
                    'micro_mg_dcs' = c(50 * 10^-6, 1000 * 10^-6),
                    'cldfrc_dp1' = c(.05, 0.25),
                    'cldfrc_dp2' = c(100, 1000),
                    'clubb_C6rt' = c(2, 6),
                    'clubb_C6rtb' = c(2, 8),
                    'clubb_C6thl' = c(2, 6),
                    'clubb_C6thlb' = c(2, 8),
                    'clubb_C8' = c(1,5), 
                    'clubb_beta' = c(1.6, 2.5),
                    'clubb_c11' = c(.2, .8),
                    'clubb_c14' = c(.4, 4),
                    'clubb_c_K10' = c(.2, 1.2),
                    'clubb_wpxp_L_thresh' = c(20, 200),
                    'dust_emis_fact' = c(.1, 1),
                    'micro_mg_accre_enhan_fact' = c(.1, 10),
                    'micro_mg_autocon_fact' = c(.005, .2),
                    'micro_mg_autocon_lwp_exp' = c(2.1, 3.3),
                    'micro_mg_berg_eff_factor' = c(.1, 1),
                    'micro_mg_effi_factor' = c(.1, 2),
                    'micro_mg_homog_size' = c(10 * 10^-6, 200 * 10^-6),
                    'micro_mg_iaccr_factor' = c(.2, 1),
                    'micro_mg_max_nicons' = c(1 * 10^5, 10000*10^6),
                    'micro_mg_vtrmi_factor' = c(.2, 5),
                    'microp_aero_npccn_scale' = c(.33, 3),
                    'microp_aero_wsub_min' = c(0, .5),
                    'microp_aero_wsub_scale' = c(.1, 5),
                    'microp_aero_wsubi_min' = c(0, .2),
                    'microp_aero_wsubi_scale' = c(.1, 5),
                    'seasalt_emis_scale' = c(.5, 2.5),
                    'sol_factb_interstitial' = c(.1, 1),
                    'sol_factic_interstitial' = c(.1, 1),
                    'zmconv_c0_lnd' = c(0.002, .1),
                    'zmconv_c0_ocn' = c(.02, .1),
                    'zmconv_capelmt' = c(35, 350),
                    'zmconv_ke' = c(10^-6, 10^-5),
                    'zmconv_ke_lnd' = c(10^-6, 10^-5),
                    'zmconv_momcd' = c(0,1),
                    'zmconv_momcu' = c(0,1),
                    'zmconv_num_cin' = c(1, 5),
                    'zmconv_tiedke_add' = c(0,2))
for (j in 1:p) {
  x_obs[,j] <- renormalize(x_obs[,j], bounds_list[[colnames(x_obs)[j]]]) - .5
}


set.seed(22)
n <- nrow(x_obs)
n_subset <- 1:n
p <- ncol(x_obs)
x_obs <- x_obs#[n_subset, colnames(x_obs) %in% vars_use]
y <- y#[n_subset]
dist_tensor <- array(dim = c(n, n, p))
for (j in 1:p) {
  dist_tensor[,,j] <- matrix(x_obs[,j], nrow = n, ncol = n) -
    matrix(x_obs[,j], nrow = n, ncol = n, byrow = T)
}
dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
dist_tensor_mat_reduced <- dist_tensor_mat[upper.tri(matrix(nrow =n, ncol = n), diag = T), ]

if (F) {
  mod_bass <- BASS::bass(x_obs, y, h1 = 800, h2 = 1,verbose=T)
  pred_C_bass <- coactivity::C_bass(mod_bass)
  # --- Sampling ---
  a_time <- Sys.time()
  
  prior_scale_matrix <- diag(x = 300/p, nrow = p)
  Cov_mat_use <- compute_cov_mat_from_C(prior_scale_matrix, dist_tensor_mat_reduced, n, 
                                        diagonal_add = 4, sigma2 = var(y))
  Cov_mat_use[lower.tri(Cov_mat_use)] <- t(Cov_mat_use)[lower.tri(Cov_mat_use)]
  y0 <- y #- mean(y)
  K_inv <- solve(Cov_mat_use)
  K_inv_y <- K_inv %*% y0
  pred_C_wycoff <- matrix(0, p, p)
  for (j in 1:p) {
    for (k in 1:p) {
      W_active_gp <- activegp::W_kappa_ij(as.matrix(x_obs), sqrt(1/diag(prior_scale_matrix)), j-1, k-1, 1)
      pred_C_wycoff[j,k] <- prior_scale_matrix[j,k] - sum(K_inv * W_active_gp) + as.vector(t(K_inv_y) %*% W_active_gp %*% K_inv_y)
    }
  }
  
  init_matrix <- matrix(nrow = p, ncol = p, 1)
  diag(init_matrix) <- 300/p
  init_par <- log(init_matrix[upper.tri(init_matrix, diag = T)])
  make_C_from_par(init_par, p)
  
  
  init_matrix <- matrix(0, nrow = p, ncol = p)
  diag(init_matrix) <- log(1)
  init_par <- init_matrix[upper.tri(init_matrix, diag = T)]
  make_C_from_par(init_par, p)
  
  # lower and upper bounds
  upper_mat <-  matrix(nrow = p, ncol = p, exp(.999))
  diag(upper_mat) <- 200
  upper_par <- log(upper_mat[upper.tri(upper_mat, diag = T)])
  lower_par <- - upper_par
  parscale_mat <- matrix(nrow = p, ncol = p, 1)
  diag(parscale_mat) <- 1/20
  parscale_par <- c(parscale_mat[upper.tri(parscale_mat, diag = T)], 1)
  
  likelihood_optim <- optim(c(init_par, log(5), log(80)), likelihood_function,
                            dist_tensor_mat_reduced = dist_tensor_mat_reduced,
                            y_use = y0, 
                            n = n, method = 'L-BFGS-B',
                            upper = c(upper_par, log(200)),
                            lower = c(lower_par, log(.00002)))#,
  #control = list('parscale' = parscale_par, 'fnscale' = 1/10^3))
  par <- likelihood_optim$par
  pred_C_likelihood <- make_C_from_par(par[-length(par)], p)
  diagonal_add <- exp(par[length(par)-1])
  sigma2 <- exp(par[length(par)])
  
  
  prior_scale_matrix <- pred_C_likelihood
  Cov_mat_use <- compute_cov_mat_from_C(prior_scale_matrix, dist_tensor_mat_reduced, n,
                                        diagonal_add = diagonal_add, sigma2 = sigma2)
  Cov_mat_use[lower.tri(Cov_mat_use)] <- t(Cov_mat_use)[lower.tri(Cov_mat_use)]
  
  K_inv <- solve(Cov_mat_use)
  K_inv_y <- K_inv %*% y
  pred_C_wycoff2 <- matrix(0, p, p)
  for (j in 1:p) {
    for (k in 1:p) {
      W_active_gp <- activegp::W_kappa_ij(as.matrix(x_obs), sqrt(1/diag(prior_scale_matrix)), j-1, k-1, 1)
      pred_C_wycoff2[j,k] <- prior_scale_matrix[j,k] - sum(K_inv * W_active_gp) + as.vector(t(K_inv_y) %*% W_active_gp %*% K_inv_y)
    }
  }
  
  normalize_matrix <- function(mat) {
    diag(1/sqrt(diag(mat))) %*% mat %*% diag(1/sqrt(diag(mat)))
  }
  
  save(diagonal_add, pred_C_wycoff, pred_C_wycoff2, mod_bass, pred_C_bass,
       pred_C_likelihood, file = 'CESM_nonstan_final.RData')
} else {
  load(file = '~/active-subspace-methods/CESM_nonstan_final.RData')
}
library(fields)

v1_wycoff <- eigen(pred_C_wycoff)$vectors[,1]
v1_wycoff2 <- eigen(pred_C_wycoff2)$vectors[,1]
v1_likelihood <- eigen(pred_C_likelihood)$vectors[,1]
v1_bass <- eigen(pred_C_bass)$vectors[,1]

library(mvtnorm)
library(rstan)

prior_configs <- list(
  "dirichlet_wishart" = list(
    stan_model_name = "~/active-subspace-methods/final/dirichlet_wishart_with_var_model.RData",
    stan_model_string_var = "sim.ss_dirichlet_wishart", 
    default_params = list(
      prior_gamma_a = 7,
      prior_gamma_b = .05
    ),
    get_specific_data_params_func = function(p, config) {
      list(
        alpha = rep(1, p),
        prior_cor_dof = p + 5,
        prior_gamma_a = config$default_params$prior_gamma_a,
        prior_gamma_b = config$default_params$prior_gamma_b,
        prior_sigma2_a = config$default_params$prior_sigma2_a,
        prior_sigma2_b = config$default_params$prior_sigma2_b,
        prior_tau_a = config$default_params$prior_tau_a,
        prior_tau_b = config$default_params$prior_tau_b
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


# Get the configuration for the selected prior
current_prior_config <- prior_configs[["dirichlet_wishart"]]

# --- Data Generation  ---
dist_tensor <- array(dim = c(n, n, p))
for (j in 1:p) {
  dist_tensor[,,j] <- matrix(x_obs[,j], nrow = n, ncol = n) - matrix(x_obs[,j], nrow = n, ncol = n, byrow = T)
}
dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
dist_tensor_mat_reduced <- dist_tensor_mat[upper.tri(matrix(nrow =n, ncol = n), diag = T), ]
# --- STAN Model Loading ---
stan_model_object <- NULL
stan_model_path <- current_prior_config$stan_model_name
stan_code_file_master <- "/home/anyarger/active-subspace-methods/final/stan_edited_model_with_var.R"

# Ensure the 'stan' directory exists
if (!dir.exists("stan")) {
  dir.create("stan")
}

if (F) {
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
get_data_input <- function(N, k, R, y, mu0, locs, specific_params) {
  base_data <- list(
    N = N,
    k = k,
    R = R,
    y = y,
    mu0 = mu0,
    locs = locs
  )
  # Combine base parameters with prior-specific parameters
  c(base_data, specific_params)
}

#diagonal_add <-  0.005

# Generate prior-specific data parameters using the function from prior_configs
#pred_C_bass_xi <- diag(pred_C_bass)^2/sum(diag(pred_C_bass)^2)
init_C <- pred_C_likelihood
init_xi <- sqrt(diag(init_C))/sum(sqrt(diag(init_C)))
init_K <- sum(diag(diag(1/init_xi) %*% init_C %*% diag(1/init_xi)))/p #+ rnorm(n_chains, sd = 20)
init_Q <- p * diag(1/init_xi) %*% 
  (init_C/init_K[1]) %*% 
  diag(1/init_xi)
current_prior_config$default_params$prior_gamma_a
set.seed(22)


init_K <- 2

desired_mean <- init_K * 1.00
desired_var <- (init_K)^(2)
alpha <- desired_mean^2/desired_var
beta <- alpha/desired_mean

plot(seq(0, 20, length.out = 600), dgamma(seq(0, 20, length.out = 600), alpha, beta))

current_prior_config$default_params$prior_gamma_a <- alpha
current_prior_config$default_params$prior_gamma_b <- beta

plot(seq(0, 200, length.out = 500), invgamma::dinvgamma(seq(0, 200, length.out = 500), 
                                                        shape = diagonal_add/8, rate = 8*diagonal_add),
     type = 'l')
prior_mean_sigma <- 25
prior_var_sigma <- 1000
current_prior_config$default_params$prior_sigma2_a <- prior_mean_sigma^2/ prior_var_sigma + 2
current_prior_config$default_params$prior_sigma2_b <- prior_mean_sigma * (prior_mean_sigma^2/ prior_var_sigma + 1)
plot(seq(0, 200, length.out = 500), 
     invgamma::dinvgamma(seq(0, 200, length.out = 500), 
                         shape = current_prior_config$default_params$prior_sigma2_a, 
                         rate =current_prior_config$default_params$prior_sigma2_b ),
     type = 'l')

prior_mean_tau <- 5
prior_var_tau <- 100
current_prior_config$default_params$prior_tau_a <- prior_mean_tau^2/ prior_var_tau + 2
current_prior_config$default_params$prior_tau_b <- prior_mean_tau * 
  (prior_mean_tau^2/ prior_var_tau + 1)
plot(seq(0, 200, length.out = 500), 
     invgamma::dinvgamma(seq(0, 100, length.out = 500), 
                         shape = current_prior_config$default_params$prior_tau_a, 
                         rate =current_prior_config$default_params$prior_tau_b ),
     type = 'l')

specific_data_params <- current_prior_config$get_specific_data_params_func(p, current_prior_config)
specific_data_params$alpha <- rep(1, p)
specific_data_params$prior_cor_dof <- p+ 2

data_input <- get_data_input(n, 
                             p, 
                             diag(x = 1, nrow = p), 
                             y, 
                             mu0 = rep(0, n), 
                             x_obs, 
                             specific_data_params)

# MCMC initializations
n_chains <- 4
it <- 1000
w <- 400
set.seed(26)

initial_values <-  lapply(1:n_chains, function(x) {
  xi_init <- init_xi + abs(rnorm(p, sd = .01))
  xi_init <- xi_init/sum(xi_init)
  K_val <- rnorm(1, mean = init_K, sd = 200)
  while(K_val <= 0) {
    K_val <- rnorm(1, mean = init_K, sd = 200)
  }
  list(Q1 = init_Q, xi = xi_init, K = K_val, 
       sigma2 = .1)
})
sapply(initial_values, function(x) min(eigen(x$Q1)$values))


# Run STAN sampling
a_time <- Sys.time()
out_model <- sampling(object = stan_model_object, data = data_input,
                      pars = c('Sigma', 'sigma2', 'tau'), iter = it, chains = n_chains, warmup = w,
                      init = initial_values,
                      cores = n_chains)
b_time <- Sys.time()
time_used <- b_time - a_time
units(time_used) <- 'mins'

extract_vals <- rstan::extract(out_model, inc_warmup = F, permuted = F)
summary_vals <- summary(out_model)

Sampled_Sigmas <- array(extract_vals[,,1:(p^2)], dim = c((it - w)*4, p, p))
sampled_sigma2 <- as.vector(extract_vals[,,p^2+1])

save(extract_vals, summary_vals, out_model,
     diagonal_add, file = 'CESM_all.RData')
q(save = 'no')
load('CESM_all.RData')

pairs(out_model, pars = c('Sigma[1,1]', 'Sigma[1,2]'))
pairs(out_model, pars = c('Sigma[1,1]', 'Sigma[1,2]', 'Sigma[2,2]'))

it <- 1000
w <- 400
n_chains <- 4
Sampled_Sigmas <- array(extract_vals[,,1:(p^2)], dim = c((it - w)*n_chains, p, p))
sampled_sigma2 <- as.vector(extract_vals[,,p^2+1])
sampled_tau <- as.vector(extract_vals[,,p^2+2])

plot(sampled_sigma2)
plot(sampled_tau)
ttraces <- t(sapply(1:(dim(Sampled_Sigmas)[1]), function(z) diag(Sampled_Sigmas[z,,])/
                      sum(diag(Sampled_Sigmas[z,,]))))
ttraces <- sapply(1:(dim(Sampled_Sigmas)[1]), function(z) 
  sum(diag(Sampled_Sigmas[z,,])))

post_mean <- colMeans(Sampled_Sigmas)
post_mean

post_mean_norm <- diag(1/sqrt(diag(post_mean))) %*% post_mean %*% 
  diag(1/sqrt(diag(post_mean)))
pm_vectors <- eigen(post_mean)$vectors
pm_values <- eigen(post_mean)$values

vectors <- t(sapply(1:(dim(Sampled_Sigmas)[1]), function(x) {
  eigen(Sampled_Sigmas[x,,])$vectors[,1]
}))
values <- t(sapply(1:(dim(Sampled_Sigmas)[1]), function(x) {
  eigen(Sampled_Sigmas[x,,])$values
}))

vectors <- t(sapply(1:nrow(vectors), function(x) {
  if (sum(vectors[x,] * pm_vectors[,1]) > 0) {
    vectors[x,]
  } else {
    -vectors[x,]
  }
}))

prop_var <- t(sapply(1:(dim(Sampled_Sigmas)[1]), function(x) {
  values <- eigen(Sampled_Sigmas[x,,])$values
  cumsum(values)/sum(values)
}))
colnames(prop_var) <- 1:45

library(dplyr)
library(ggplot2)
prop_var_df <- data.frame(prop_var) %>% 
  tidyr::pivot_longer(cols = 1:45, names_to = 'Dimension') %>%
  dplyr::mutate(Dimension =as.integer(substr(Dimension, 2, nchar(Dimension)))) 


df_prop <- rbind(data.frame(Dimension = 1:45, 
                            value = cumsum(eigen(pred_C_bass)$values)/
                              sum(eigen(pred_C_bass)$values), 
                            type = 'BASS'), 
                 data.frame(Dimension = 1:45, 
                            value = cumsum(eigen(pred_C_likelihood)$values)/
                              sum(eigen(pred_C_likelihood)$values), 
                            type = 'MLE'))
ggplot() +
  geom_boxplot(data = prop_var_df, aes(x = Dimension, y = value, group = Dimension)) + 
  geom_point(data = df_prop, 
             aes(x = Dimension, y = value, group = Dimension, color = type,
                 shape = type), size = 1.5) + 
  scale_shape_manual(values = c(15,17)) +
  labs(y= 'Proportion of variance explained', shape = 'Estimate',
       color = 'Estimate') +
  theme_bw() +
  theme(text = element_text(family = 'Arial'))
ggsave('CESM_prop_var_explained.png', height = 4*.8, width = 6.8*.8)


vectors2 <- t(sapply(1:(dim(Sampled_Sigmas)[1]), function(x) {
  eigen(Sampled_Sigmas[x,,])$vectors[,2]
}))
vectors2 <- t(sapply(1:nrow(vectors2), function(x) {
  if (sum(vectors2[x,] * pm_vectors[,2]) > 0) {
    vectors2[x,]
  } else {
    -vectors2[x,]
  }
}))

AS <- t(sapply(1:(dim(Sampled_Sigmas)[1]), function(x) {
  vectors2[x,]^2 * values[x,2] + vectors[x,]^2 * values[x,1]
}))
colnames(AS) <- colnames(x_obs)

ggplot(data = data.frame(AS) %>% tidyr::pivot_longer(cols = colnames(x_obs)) %>%
         mutate(name = factor(name, levels = colnames(x_obs)))) +
  geom_boxplot(aes(x = name, y = value, group = name))+
  theme_bw() + theme(text = element_text(family = 'Arial'), 
                     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + 
  labs(x = 'Input parameter',
       y = 'Posterior samples of activity scores\nbased on two dimensions')
ggsave('CESM_activity_score.png', height = 4*1.4, width = 6.8*1.4)




library(tidyr)
library(dplyr)
df <- left_join(data.frame(pivot_longer(as.data.frame(vectors) * sqrt(values[,1]), 
                                        cols = 1:45, names_to = 'v1'), index = 1:(108000)), 
                data.frame(pivot_longer(as.data.frame(vectors2)* sqrt(values[,2]), 
                                        cols = 1:45, names_to = 'v1'),index = 1:(108000)),
                by = c('index', 'v1'))

projected_points <- as.matrix(x_obs) %*%
  pm_vectors[,1:2] %*% 
  diag(sqrt(eigen(post_mean)$values[1:2]),ncol = 2)
activity_scores <- rowSums(pm_vectors[,1:2]^2 %*% 
                             diag(sqrt(eigen(post_mean)$values[1:2]),ncol = 2))
top_act_scores <- order(-activity_scores)[1:5]
cbind(top_act_scores, colnames(x_obs)[top_act_scores], activity_scores[top_act_scores])

plot(x_obs[['micro_mg_effi_factor']], RESTOM_vector)

set.seed(22)
random_indexes <- sample(1:2400, size = 100)
df_random <- cbind(data.frame(pivot_longer(as.data.frame(vectors[random_indexes,]) * 
                                             sqrt(values[random_indexes,1]),
                                           cols = 1:45, names_to = 'v12', values_to = 'value1')), 
                   data.frame(pivot_longer(as.data.frame(vectors2[random_indexes,])* 
                                             sqrt(values[random_indexes,2]), 
                                           cols = 1:45, names_to = 'v1', values_to = 'value2')),
                   variable = factor(colnames(x_obs), levels = colnames(x_obs)))



ggplot() + 
  geom_point(data = data.frame(projected_points, RESTOM_vector), 
             aes(x = X1, y = X2, color = RESTOM_vector)) +
  coord_equal(xlim = c(-.9, NA)) + 
  # geom_segment(data = df_random %>% filter(variable == 'micro_mg_vtrmi_factor'),
  #              aes(x = value1 * .7, y = value2 *.7),
  #              arrow = grid::arrow(ends = 'first', length = unit(.2, 'cm')),
  #              xend = 0, yend = 0, size = .2,
  #              alpha = .2) + 
  geom_label(data = data.frame(pm_vectors[top_act_scores,1:2]  %*% diag(sqrt(eigen(post_mean)$values[1:2]),ncol = 2), 
                               variable = factor(colnames(x_obs)[top_act_scores], 
                                                 levels = colnames(x_obs)[top_act_scores])), 
             aes(x = X1, y =ifelse(variable == 'clubb_C2rt', X2 + .07, X2), label = variable), 
             nudge_x = -.17, size = 2.5, alpha = .2,nudge_y = .05,
             family = 'Arial') +
  geom_segment(data = data.frame(pm_vectors[top_act_scores,1:2]  %*% diag(sqrt(eigen(post_mean)$values[1:2]),ncol = 2), 
                                 variable = factor(colnames(x_obs)[top_act_scores], 
                                                   levels = colnames(x_obs)[top_act_scores])), 
               aes(x = X1, y = X2, 
                   xend =0, yend = 0, linetype = variable), 
               arrow = grid::arrow(ends = 'first', length = unit(.2, 'cm')), 
               linewidth = .7) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     text = element_text(family = 'Arial'),
                     legend.spacing.y = unit(-0, "cm")) + 
  labs(linetype = 'Variable', color = expression('RESTOM (W/'*m^2*')'),
       x = 'Active subspace direction 1',
       y = 'Active subspace direction 2') + 
  scale_color_viridis_c() + 
  guides(label = guide_legend(order = 0),
         linetype = guide_legend(order = 1),
         color = guide_colorbar(order = 10))
ggsave('CESM_RESTOM_biplot.png', height = 4*.9, width = 6.8*.9)

extract_vals_trace <- rstan::extract(out_model, permuted = F)
Sigma_samples <- extract_vals_trace[,,1:(45)^2]
library(ggplot2)
library(dplyr)
Sigma_df <- data.frame(value = as.vector(Sigma_samples), 
                       chain = rep(1:4, each = dim(Sigma_samples)[1]),
                       Row = rep(rep(1:45, each = 45), each = n_chains*dim(Sigma_samples)[1]),
                       Column = rep(rep(1:45, times = 45), each = n_chains*dim(Sigma_samples)[1]),
                       order = 1:600)

ggplot(data = Sigma_df %>%filter(chain == 1, Row %in% 1:10, 
                                 Column %in% 1:10), aes(x = order, y = value)) + 
  geom_line(linewidth = .2) + 
  facet_grid(Row~Column, labeller = label_both, switch = 'y') +
  labs(x = 'MCMC iteration', y = 'Estimated entry of C') +
  theme_bw() +
  theme(text = element_text(family = 'Arial'))
ggsave('CESM_sampling.png', height = 5.3*1.2, width = 7*1.2)


other_df <- data.frame(value = as.vector(extract_vals_trace[,,((45)^2 + 1):((45)^2 + 2)]), 
                       Chain = rep(1:4, each = dim(Sigma_samples)[1]),
                       Variable = rep(c('sigma2', 'tau2'), each = n_chains*dim(Sigma_samples)[1]),
                       order = 1:600)


ggplot(data = other_df , aes(x = order, y = value)) + 
  geom_line(linewidth = .2) + 
  facet_grid(Variable~Chain, labeller = label_both, switch = 'y', scales = 'free_y') +
  labs(x = 'MCMC iteration', y = 'Sampled value') +
  theme_bw() +
  theme(text = element_text(family = 'Arial'))
ggsave('CESM_sampling_sigma2_tau2.png', height = 5.3 *.8, width = 7*.8)


                  
C_df <- data.frame(C_value = c(as.vector(pred_C_wycoff),
                               as.vector(pred_C_wycoff2),
                               as.vector(pred_C_likelihood),
                               as.vector(pred_C_bass),
                               as.vector(post_mean)),
                   C_value_normalized = c(as.vector(pred_C_wycoff)/ sum(diag(pred_C_wycoff)),
                                          as.vector(pred_C_wycoff2)/ sum(diag(pred_C_wycoff2)),
                                          as.vector(pred_C_likelihood)/ sum(diag(pred_C_likelihood)),
                                          as.vector(pred_C_bass)/ sum(diag(pred_C_bass)),
                                          as.vector(post_mean)/ sum(diag(post_mean))),
                   C_value_normalized2 = c(as.vector(diag(1/sqrt(diag(pred_C_wycoff))) %*% pred_C_wycoff %*% diag(1/sqrt(diag(pred_C_wycoff)))),
                                           as.vector(diag(1/sqrt(diag(pred_C_wycoff2))) %*% pred_C_wycoff2 %*% diag(1/sqrt(diag(pred_C_wycoff2)))),
                                           as.vector(diag(1/sqrt(diag(pred_C_likelihood))) %*% pred_C_likelihood %*% diag(1/sqrt(diag(pred_C_likelihood)))),
                                           as.vector(diag(1/sqrt(diag(pred_C_bass))) %*% pred_C_bass %*% diag(1/sqrt(diag(pred_C_bass)))),
                                           as.vector(diag(1/sqrt(diag(post_mean))) %*% post_mean %*% diag(1/sqrt(diag(post_mean))))),
                   row = factor(as.vector(matrix(colnames(x_obs), p, p)),
                                levels = colnames(x_obs)),
                   column = factor(as.vector(matrix(colnames(x_obs), p, p, byrow = T)), 
                                   rev(colnames(x_obs))),
                   type = rep(c('Wycoff', 'Wycoff2', 'Likelihood', 'BASS', 'Bayes'), 
                              each = p^2))



ggplot(data = C_df  %>% filter(type !='Wycoff') %>%
         mutate(type = ifelse(type == 'Wycoff2', 'GP conditional', type)), aes(x = row, column, fill = C_value_normalized)) + 
  geom_tile(color = 'grey40') +
  facet_wrap(~type) + 
  coord_equal() +
  scale_fill_viridis_c() + 
  theme_bw() + theme(text = element_text(family = 'Arial'),
                     strip.text = element_text(size = 12),
                     legend.title = element_text(size = 12),
                     legend.text = element_text(size = 10),
                     axis.text.x = element_blank(), 
                     axis.title = element_blank(),
                     axis.ticks.x = element_blank(),
                     axis.ticks.y = element_line(linewidth = .2),
                     axis.text.y = element_text(size = 4)) + 
  labs(fill = 'Normalized\nC estimate', x = '', y = 'Column')
ggsave('CESM_C_compare.png', height = 4.5 * 1.3, width = 5.7 * 1.3)


ggplot(data = C_df  %>% filter(type !='Wycoff') %>%
         mutate(type = ifelse(type == 'Wycoff2', 'GP conditional', type)), 
       aes(x = row, column, fill = C_value_normalized2)) + 
  geom_raster() +
  facet_wrap(~type) + 
  coord_equal() +
  scale_fill_gradient2() + 
  theme_bw() + theme(text = element_text(family = 'Arial'),
                     axis.text = element_blank(), 
                     axis.title = element_blank()) + 
  labs(fill = 'Normalized\nC estimate', x = '', y = 'Column')
ggsave('CESM_C_compare.png', height = 4 * 1.3, width = 4.5 * 1.3)


