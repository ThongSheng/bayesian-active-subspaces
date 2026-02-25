load('data/E3SM_RESTOM.RData')
library(mvtnorm)
library(BASS)
library(concordance)
library(activegp)
n <- length(E3SM_RESTOM)


compute_cov_mat_from_C <- function(C, dist_tensor_mat_reduced, n, diagonal_add = .00001) {
  test_distances <- rowSums((dist_tensor_mat_reduced %*% C) * dist_tensor_mat_reduced)
  dist_mat_C <- matrix(nrow = n, ncol = n)
  dist_mat_C[upper.tri(dist_mat_C, diag = T)] <- test_distances
  cov_mat <- exp(-dist_mat_C/2)
  diag(cov_mat) <- diag(cov_mat) + diagonal_add
  cov_mat
}
# --- MLE functions ---
likelihood_C <- function(C, dist_tensor_mat_reduced, n, diagonal_add, compute_gradient = F, 
                         upper_tri_n = NULL, lower_tri_n = NULL, upper_tri_p = NULL, lower_tri_p = NULL,
                         y_use) {
  Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n, diagonal_add)
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
  C_mat_use <- make_C_from_par(par[1:(length(par)-1)], ncol(dist_tensor_mat_reduced))
  diagonal_add <- exp(par[length(par)])
  if (sum(is.na(C_mat_use)) > 0) {
    return(10^6)
  }
  if (min(eigen(C_mat_use)$value) < 10^-6) {
    return(10^6)
  }
  ll_val <- likelihood_C(C_mat_use, dist_tensor_mat_reduced, n, diagonal_add, y_use = y_use)
  -ll_val
}
y <- E3SM_RESTOM
x_obs <- e3sm_parameters
p <- ncol(x_obs)

renormalize <- function(x, bounds = c(range(x))) {
  (x - bounds[1])/(bounds[2] - bounds[1])
}
bounds_list <- list('clubb_c1' = c(1, 5),
                    'clubb_gamma_coef' = c(.1, .5),
                    'zmconv_dmpdz' = c(-.002, -0.0001),
                    'ice_sed_ai' = c(350, 1400),
                    'zmconv_tau' = c(1800, 14400))
for (j in 1:p) {
  x_obs[,j] <- renormalize(x_obs[,j], bounds_list[[colnames(x_obs)[j]]]) - .5
}

dist_tensor <- array(dim = c(n, n, p))
for (j in 1:p) {
  dist_tensor[,,j] <- matrix(x_obs[,j], nrow = n, ncol = n) -
    matrix(x_obs[,j], nrow = n, ncol = n, byrow = T)
}
dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
dist_tensor_mat_reduced <- dist_tensor_mat[upper.tri(matrix(nrow =n, ncol = n), diag = T), ]

mod_bass <- BASS::bass(x_obs, y, verbose=T)
pred_C_bass <- coactivity::C_bass(mod_bass)

# --- Sampling ---
a_time <- Sys.time()
prior_scale_matrix <- diag(x = 300/p, nrow = p)
Cov_mat_use <- compute_cov_mat_from_C(prior_scale_matrix, dist_tensor_mat_reduced, n, diagonal_add = 4)
Cov_mat_use[lower.tri(Cov_mat_use)] <- t(Cov_mat_use)[lower.tri(Cov_mat_use)]
y0 <- y 
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

# lower and upper bounds
upper_mat <-  matrix(nrow = p, ncol = p, exp(.999))
diag(upper_mat) <- NA
upper_par <- log(upper_mat[upper.tri(upper_mat, diag = T)])
lower_par <- - upper_par

likelihood_optim <- optim(c(init_par, log(5)), likelihood_function,
                          dist_tensor_mat_reduced = dist_tensor_mat_reduced,
                          y_use = y0, 
                          n = n, method = 'L-BFGS-B',
                          upper = c(upper_par, log(200)),
                          lower = c(lower_par, log(.00002)))
par <- likelihood_optim$par
pred_C_likelihood <- make_C_from_par(par[-length(par)], p)
diagonal_add <- exp(par[length(par)])


prior_scale_matrix <- pred_C_likelihood
Cov_mat_use <- compute_cov_mat_from_C(prior_scale_matrix, dist_tensor_mat_reduced, n,
                                      diagonal_add = diagonal_add)
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

v1_wycoff <- eigen(pred_C_wycoff)$vectors[,1]
v1_wycoff2 <- eigen(pred_C_wycoff2)$vectors[,1]
v1_likelihood <- eigen(pred_C_likelihood)$vectors[,1]
v1_bass <- eigen(pred_C_bass)$vectors[,1]

#Bayesian
library(mvtnorm)
library(rstan)
prior_configs <- list(
  "dirichlet_wishart" = list(
    stan_model_name = "stan/dirichlet_wishart_with_var_model.RData",
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
stan_model_path <- current_prior_config$stan_model_name
stan_code_file_master <- "stan/stan_edited_model_with_var.R"

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


init_C <- pred_C_likelihood
init_xi <- sqrt(diag(init_C))/sum(sqrt(diag(init_C)))
init_K <- sum(diag(diag(1/init_xi) %*% init_C %*% diag(1/init_xi)))/p #+ rnorm(n_chains, sd = 20)
init_Q <- p * diag(1/init_xi) %*% 
  (init_C/init_K[1]) %*% 
  diag(1/init_xi)
current_prior_config$default_params$prior_gamma_a
set.seed(22)


init_K <- 2

desired_mean <- 2
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

# Generate prior-specific data parameters using the function from prior_configs
specific_data_params <- current_prior_config$get_specific_data_params_func(p, current_prior_config)
data_input <- get_data_input(n, 
                             p, 
                             diag(x = 1, nrow = p), 
                             y, 
                             mu0 = rep(0, n), 
                             x_obs, 
                             diagonal_add, 
                             specific_data_params)

# MCMC initializations
n_chains <- 4
it <- 1500
w <- 500
set.seed(22)

initial_values <-  lapply(1:n_chains, function(x) {
  mat <- matrix(runif(p*p, 0, .01), nrow = p, ncol = p)
  diag(mat) <- 1
  mat <- crossprod(mat)
  list(Q1 = mat, xi = rep(1/p, p), K = runif(1))
})
sapply(initial_values, function(x) min(eigen(x$Q1)$values))


# Run STAN sampling
a_time <- Sys.time()
out_model <- sampling(object = stan_model_object, data = data_input,
                      pars = c('Sigma', 'tau', 'sigma2'), iter = it, chains = n_chains, warmup = w,
                      init = initial_values,
                      cores = n_chains)
b_time <- Sys.time()
time_used <- b_time - a_time
units(time_used) <- 'mins'

extract_vals <- rstan::extract(out_model)
summary_vals <- summary(out_model)

plot(out_model)
#pairs(out_model, vars = c('Sigma[1,1]'))

save(extract_vals, summary_vals, out_model, time_used, file = 'e3sm_samples.RData')

#load('e3sm_samples.RData')
library(ggplot2)
library(tidyr)
library(dplyr)

extract_vals_trace <- rstan::extract(out_model, permuted = F)
Sigma_samples <- extract_vals_trace[,,1:25]

Sigma_df <- data.frame(value = as.vector(Sigma_samples), 
                       chain = rep(1:4, each = dim(Sigma_samples)[1]),
                       Row = rep(rep(1:5, each = 5), each = n_chains*dim(Sigma_samples)[1]),
                       Column = rep(rep(1:5, times = 5), each = n_chains*dim(Sigma_samples)[1]),
                       order = 1:1000)

ggplot(data = Sigma_df %>%filter(chain == 1), aes(x = order, y = value)) + 
  geom_line(linewidth = .2) + 
  facet_grid(Row~Column, labeller = label_both, switch = 'y') +
  labs(x = 'MCMC iteration', y = 'Sampled value of the entry of C') +
  theme_bw() +
  theme(text = element_text(family = 'Arial'))
ggsave('images/E3SM_sampling_C.png', height = 5.3*.8, width = 7*.8)

other_df <- data.frame(value = as.vector(extract_vals_trace[,,26:27]), 
                       Chain = rep(1:4, each = dim(Sigma_samples)[1]),
                       Variable = rep(c('tau2', 'sigma2'), each = n_chains*dim(Sigma_samples)[1]),
                       order = 1:1000)


ggplot(data = other_df , aes(x = order, y = value)) + 
  geom_line(linewidth = .2) + 
  facet_grid(Variable~Chain, labeller = label_both, switch = 'y', scales = 'free_y') +
  labs(x = 'MCMC iteration', y = 'Sampled value') +
  theme_bw() +
  theme(text = element_text(family = 'Arial'))
ggsave('images/E3SM_sampling_sigma2_tau2.png', height = 5.3 *.8, width = 7*.8)


vectors <- t(sapply(1:(dim(extract_vals$Sigma)[1]), function(x) {
  eigen(extract_vals$Sigma[x,,])$vectors[,1]
}))
values <- t(sapply(1:(dim(extract_vals$Sigma)[1]), function(x) {
  eigen(extract_vals$Sigma[x,,])$values
}))

post_mean <- colMeans(extract_vals$Sigma)
pm_vectors <- eigen(post_mean)$vectors
pm_values <- eigen(post_mean)$values

# align first eigenvector with posterior mean first eigenvector
vectors <- t(sapply(1:nrow(vectors), function(x) {
  if (sum(vectors[x,] * pm_vectors[,1]) > 0) {
    vectors[x,]
  } else {
    -vectors[x,]
  }
}))

vectors2 <- t(sapply(1:(dim(extract_vals$Sigma)[1]), function(x) {
  eigen(extract_vals$Sigma[x,,])$vectors[,2]
}))
vectors2 <- t(sapply(1:nrow(vectors2), function(x) {
  if (sum(vectors2[x,] * pm_vectors[,2]) > 0) {
    vectors2[x,]
  } else {
    -vectors2[x,]
  }
}))


AS <- t(sapply(1:(dim(extract_vals$Sigma)[1]), function(x) {
  vectors2[x,]^2 * values[x,2] + vectors[x,]^2 * values[x,1]
}))
colnames(AS) <- colnames(e3sm_parameters)

ggplot(data = data.frame(AS) %>% tidyr::pivot_longer(cols = colnames(e3sm_parameters)) %>%
         mutate(name = factor(name, levels = colnames(e3sm_parameters)))) +
  geom_boxplot(aes(x = name, y = value, group = name), outlier.size = .5)+
  theme_bw() + theme(text = element_text(family = 'Arial')) + 
  labs(x = 'Input parameter',
       y = 'Posterior samples of activity scores\nbased on two dimensions')
ggsave('images/E3SM_activity_score.png', height = 4*.8, width = 6.8*.8)




prop_var <- t(sapply(1:(dim(extract_vals$Sigma)[1]), function(x) {
  cumsum( values[x,])/sum( values[x,])
}))
colnames(prop_var) <- 1:5
prop_var_df <- data.frame(prop_var) %>% 
  tidyr::pivot_longer(cols = 1:5, names_to = 'Dimension') %>%
  dplyr::mutate(Dimension =as.integer(substr(Dimension, 2, nchar(Dimension)))) 


df_prop <- rbind(data.frame(Dimension = 1:5, 
                            value = cumsum(eigen(pred_C_bass)$values)/
                              sum(eigen(pred_C_bass)$values), 
                            type = 'BASS'), 
                 data.frame(Dimension = 1:5, 
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
ggsave('images/E3SM_prop_var_explained.png', height = 4*.8, width = 6.8*.8)

df <- left_join(data.frame(pivot_longer(as.data.frame(vectors) * sqrt(values[,1]), cols = 1:5, names_to = 'v1'), index = 1:(20000)), 
                data.frame(pivot_longer(as.data.frame(vectors2)* sqrt(values[,2]), cols = 1:5, names_to = 'v1'),index = 1:(20000)),
                by = c('index', 'v1'))

sum(v1_bass * pm_vectors[,1])
sum(v1_wycoff2 * pm_vectors[,1])
sum(v1_wycoff * pm_vectors[,1])
sum(v1_likelihood * pm_vectors[,1])

projected_points <- (x_obs) %*% pm_vectors[,1:2] %*% diag(sqrt(eigen(post_mean)$values[1:2]),ncol = 2)

set.seed(22)
random_indexes <- sample(1:4000, size = 100)
df_random <- cbind(data.frame(pivot_longer(as.data.frame(vectors[random_indexes,]) * 
                                             sqrt(values[random_indexes,1]),
                                           cols = 1:5, names_to = 'v12', values_to = 'value1')), 
                   data.frame(pivot_longer(as.data.frame(vectors2[random_indexes,])* 
                                             sqrt(values[random_indexes,2]), 
                                           cols = 1:5, names_to = 'v1', values_to = 'value2')),
                   variable = factor(colnames(x_obs), levels = colnames(x_obs)))

geom_point(data = data.frame(projected_points, E3SM_RESTOM), 
             aes(x = X1, y = X2, color =E3SM_RESTOM)) +
  scale_color_viridis_c() +
  coord_equal(xlim = c(-1, 1.02)) + 
  #geom_point(data = df, aes(x = value.x/10, y = value.y/10, fill = v1), size = .2)+
  # geom_segment(data = df_random %>% filter(variable == 'clubb_gamma_coef'),
  #              aes(x = value1 * .7, y = value2 *.7),
  #              arrow = grid::arrow(ends = 'first', length = unit(.2, 'cm')),
  #              xend = 0, yend = 0, size = .2,
  #              alpha = .2) + 
  geom_label(data = data.frame(pm_vectors[,1:2]  %*% diag(sqrt(eigen(post_mean)$values[1:2]),ncol = 2), 
                               variable = factor(colnames(x_obs), levels = colnames(x_obs))), 
             aes(x = X1, y = ifelse(variable == 'zmconv_dmpdz', X2+ .1, 
                                       X2), label = variable), nudge_x = -.3, size = 3, alpha = .2,
             family = 'Arial') +
  geom_segment(data = data.frame(pm_vectors[,1:2]  %*% diag(sqrt(eigen(post_mean)$values[1:2]),ncol = 2), 
                                 variable = factor(colnames(x_obs), levels = colnames(x_obs))), 
               aes(x = X1, y = X2, xend =0, yend = 0, linetype = variable), 
               arrow = grid::arrow(ends = 'first', length = unit(.2, 'cm')), 
               linewidth = .7) +
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     text = element_text(family = 'Arial'),
                     legend.spacing.y = unit(-0, "cm")) + 
  labs(linetype = 'Variable', color = expression('RESTOM (W/'*m^2*')'),
       x = 'Active subspace direction 1',
       y = 'Active subspace direction 2') +
  guides(label = guide_legend(order = 0),
         linetype = guide_legend(order = 1),
         color = guide_colorbar(order = 10))
ggsave('E3SM_RESTOM_biplot.png', height = 4*.9, width = 6.8*.9)

ggsave('E3SM_RESTOM_biplot.png', height = 4*.9, width = 6.8*.9)

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
                   row = factor(as.vector(matrix(colnames(x_obs), p, p)),
                                levels = colnames(x_obs)),
                   column = factor(as.vector(matrix(colnames(x_obs), p, p, byrow = T)), 
                                   rev(colnames(x_obs))),
                   type = rep(c('Wycoff', 'Wycoff2', 'Likelihood', 'BASS', 'Bayes'), 
                              each = p^2))

ggplot(data = C_df  %>% filter(type !='Wycoff') %>%
         mutate(type = ifelse(type == 'Wycoff2', 'Wycoff', type)), aes(x = row, column, fill = C_value_normalized)) + 
  geom_raster() +
  facet_wrap(~type) + 
  coord_equal() +
  scale_fill_gradient2() + 
  theme_bw() + theme(text = element_text(family = 'Arial'),
                     axis.text.x = element_text(angle = 90, vjust = .5, hjust = .5), axis.title = element_blank()) + 
  labs(fill = 'Normalized\nC estimate', x = '', y = 'Column')
ggsave('E3SM_C_compare.png', height = 4 * 1.3, width = 4.5 * 1.3)

# other priors


specific_data_params$alpha <- c(.4, .3, .3, .05, .05) * 40
specific_data_params$prior_cor_dof <- p+ 50


prior_cor_matrix <- matrix(.5, p, p) 
diag(prior_cor_matrix) <- 1
data_input <- get_data_input(n, 
                             p, 
                             prior_cor_matrix, 
                             y, 
                             mu0 = rep(0, n), 
                             x_obs, 
                             diagonal_add, 
                             specific_data_params)

# MCMC initializations
set.seed(22)

# Run STAN sampling
a_time <- Sys.time()
out_model_p2 <- sampling(object = stan_model_object, data = data_input,
                      pars = c('Sigma', 'tau', 'sigma2'), iter = it, chains = n_chains, warmup = w,
                      init = initial_values,
                      cores = n_chains)
b_time <- Sys.time()
time_used_p2 <- b_time - a_time
units(time_used) <- 'mins'

extract_vals_p2 <- rstan::extract(out_model_p2)
summary_vals_p2 <- summary(out_model_p2)

save(extract_vals_p2, summary_vals_p2, out_model_p2, time_used_p2, file = 'e3sm_samples_prior2.RData')
post_mean_2 <- colMeans(extract_vals_p2$Sigma)


specific_data_params$alpha <- c(1,1,1,1,1)
specific_data_params$prior_cor_dof <- p+ 180

prior_cor_matrix <- matrix(-.2, p, p) 
diag(prior_cor_matrix) <- 1

data_input <- get_data_input(n, 
                             p, 
                             prior_cor_matrix, 
                             y, 
                             mu0 = rep(0, n), 
                             x_obs, 
                             diagonal_add, 
                             specific_data_params)

# MCMC initializations
set.seed(22)

# Run STAN sampling
a_time <- Sys.time()
out_model_p3 <- sampling(object = stan_model_object, data = data_input,
                         pars = c('Sigma', 'tau', 'sigma2'), iter = it, chains = n_chains, warmup = w,
                         init = initial_values,
                         cores = n_chains)
b_time <- Sys.time()
time_used_p3 <- b_time - a_time
units(time_used) <- 'mins'

extract_vals_p3 <- rstan::extract(out_model_p3)
summary_vals_p3 <- summary(out_model_p3)

save(extract_vals_p3, summary_vals_p3, out_model_p3, time_used_p3, file = 'e3sm_samples_prior3.RData')
post_mean_3 <- colMeans(extract_vals_p3$Sigma)



C_df <- data.frame(C_value = c(as.vector(pred_C_wycoff),
                               as.vector(pred_C_wycoff2),
                               as.vector(pred_C_likelihood),
                               as.vector(pred_C_bass),
                               as.vector(post_mean),
                               as.vector(post_mean_2),
                               as.vector(post_mean_3)),
                   C_value_normalized = c(as.vector(pred_C_wycoff)/ sum(diag(pred_C_wycoff)),
                                          as.vector(pred_C_wycoff2)/ sum(diag(pred_C_wycoff2)),
                                          as.vector(pred_C_likelihood)/ sum(diag(pred_C_likelihood)),
                                          as.vector(pred_C_bass)/ sum(diag(pred_C_bass)),
                                          as.vector(post_mean)/ sum(diag(post_mean)),
                                          as.vector(post_mean_2)/ sum(diag(post_mean_2)),
                                          as.vector(post_mean_3)/ sum(diag(post_mean_3))),
                   row = factor(as.vector(matrix(colnames(x_obs), p, p)),
                                levels = colnames(x_obs)),
                   column = factor(as.vector(matrix(colnames(x_obs), p, p, byrow = T)), 
                                   rev(colnames(x_obs))),
                   type = rep(c('Wycoff', 'Wycoff2', 'Likelihood', 'BASS', 'Bayes: Prior 1', 
                                'Bayes: Prior 2',  'Bayes: Prior 3'), 
                              each = p^2))
ggplot(data = C_df  %>% filter(type !='Wycoff') %>%
         mutate(type = ifelse(type == 'Wycoff2', 'Wycoff', type),
                type = factor(type, levels = c('BASS', 'Likelihood', 'Wycoff',
                                               'Bayes: Prior 1', 'Bayes: Prior 2', 'Bayes: Prior 3'))), aes(x = row, column, fill = C_value_normalized)) + 
  geom_raster() +
  facet_wrap(~type) + 
  coord_equal() +
  scale_fill_gradient2() + 
  theme_bw() + theme(text = element_text(family = 'Arial'),
                     axis.text.x = element_text(angle = 90, vjust = .5, hjust = .5, size = 6),
                     axis.text.y = element_text(size = 6), 
                     axis.title = element_blank()) + 
  labs(fill = 'Normalized\nC estimate', x = '', y = 'Column')
ggsave('images/E3SM_C_compare_all.png', height = 3 * 1.2, width = 4.5 * 1.2)
