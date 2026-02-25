library(rstan)
library(ggplot2)

prior_configs <- list(
  "dirichlet_wishart" = list(
    stan_model_name = "final/dirichlet_wishart_model.RData",
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
  "gp_reduce" = list(
    stan_model_name = "final/gp_reduce_model.RData",
    stan_model_string_var = "sim.gp_reduce",
    get_specific_data_params_func = function(p, config) {
      list(
        k_reduce = 2
      )
    },
    get_inits_func = function(n_chains, p, config) {
      lapply(1:n_chains, function(x) {
        mat <- matrix(runif(p*p, 0, .5), nrow = p, ncol = p)
        diag(mat) <- 1
        mat <- crossprod(mat)
      })
    }
  ),
  "lognormal_inverse_wishart" = list(
    stan_model_name = "final/lniw_gp_rescale_model.RData",
    stan_model_string_var = "sim.sslniw_gp_rescale",
    default_params = list(
      prior_rescale_mean = 0,
      prior_rescale_var = log(2)
    ),
    get_specific_data_params_func = function(p, config) {
      list(
        prior_lgn_mean = array(rep(0, p)),
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
  )
)

# --- Simulation Grid Definition ---
grid <- expand.grid(
  'p' = c(2, 10, 20),
  'n' = c(20, 100, 150),
  'seed' = 1:30,
  'func' = c('determ_full', 'determ_2d', 'GP_full', 'GP_2d'),
  'prior_choice' = c("dirichlet_wishart", "gp_reduce", "lognormal_inverse_wishart")
)
array_id <- 544
p <- 4
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


# --- STAN Model Loading ---
stan_model_object <- NULL
stan_model_path <- current_prior_config$stan_model_name
stan_code_file_master <- "final/stan_edited_model.R"

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

# Generate prior-specific data parameters using the function from prior_configs
specific_data_params <- current_prior_config$get_specific_data_params_func(p, current_prior_config)
data_input <- get_data_input(n, 
                             p, 
                             diag(x = 1, nrow = p), 
                             y, 
                             rep(0, n), 
                             x_obs, 
                             .00001, 
                             specific_data_params)

# MCMC initializations
n_chains <- 4
it <- 1500
w <- 500

initial_values <- current_prior_config$get_inits_func(n_chains, p, current_prior_config)
if (selected_prior_name == 'gp_reduce') {
  initial_values <- 'random'
} else {
  print(sapply(initial_values, function(x) min(eigen(x$Q1)$values)))
}

# Run STAN sampling
a_time <- Sys.time()
out_model <- sampling(object = stan_model_object, data = data_input,
                      pars = c('Sigma'), iter = it, chains = n_chains, warmup = w,
                      init = initial_values,
                      cores = n_chains)
b_time <- Sys.time()
time_used <- b_time - a_time
units(time_used) <- 'mins'

extract_vals <- extract(out_model, inc_warmup=F, permuted=F)
summary_vals <- summary(out_model)



# Visuals for 4.3
detach("package:rstan", unload = TRUE)
library(reshape2)
library(dplyr)
library(tidyr)

C
(pred_C <- matrix(colMeans(extract_vals[, , 1:16], dims=2), ncol = 4))

# Heatmaps
melted_C <- melt(C)
melted_pred_C <- melt(pred_C)
colnames(melted_C) <- colnames(melted_pred_C) <- c("Row", "Column", "Value")

p1 <- ggplot(melted_C, aes(x = Column, y = Row, fill = Value)) +
  geom_tile(color = "black") +
  scale_y_reverse() + 
  scale_fill_gradient2(low = "steelblue", high = "tomato", limits = c(-20, 110)) +
  labs(title = "Heatmap of true C") +
  theme_minimal() + 
  theme(text = element_text(family = 'Arial')) +
  coord_equal()

p2 <- ggplot(melted_pred_C, aes(x = Column, y = Row, fill = Value)) +
  geom_tile(color = "black") +
  scale_y_reverse() + 
  scale_fill_gradient2(low = "steelblue", high = "tomato", limits = c(-20, 110)) +
  labs(title = expression(paste("Heatmap of estimated ", hat(C)))) +
  theme_minimal() +
  theme(text = element_text(family = 'Arial')) +
  coord_equal()

(heat <- p1 + p2)
ggsave("images/heatmaps.png", plot = heat, width = 10, height = 3.5, units = "in", dpi = 300)

# Posterior Distribution
plot_data <- melt(extract_vals) %>%
  filter(parameters != "lp__") %>%
  extract(parameters, into = c("Row", "Column"), regex = "\\[(\\d+),(\\d+)\\]", convert = TRUE)

(posterior <- ggplot(plot_data, aes(x = value)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black") +
  geom_vline(data = melted_C, aes(xintercept = Value), linetype = "dashed", color = "red") +
  facet_grid(Row ~ Column, labeller = label_both) + 
  labs(x="", y="") +
  theme_minimal() +
  theme(text = element_text(family = "Arial")))

ggsave("images/posterior.png", plot = posterior, width = 10, height = 5, units = "in", dpi = 300)

# Trace plots
p1 <- ggplot(data.frame(Value = extract_vals[,1,1]), aes(x = 1:1000, y = Value)) +
  geom_line() +
  geom_hline(yintercept = C[1,1], linetype = "dashed", color = "red") +
  labs(title = paste0("C[1,1] (R_hat = ", round(summary_vals$summary[1,10]), ", n_eff = ", round(summary_vals$summary[1,9]), ")"), x = "Iterations", y="") +
  theme_minimal() +
  theme(text = element_text(family = 'Arial'))

p2 <- ggplot(data.frame(Value = extract_vals[,1,2]), aes(x = 1:1000, y = Value)) +
  geom_line() +
  geom_hline(yintercept = C[1,2], linetype = "dashed", color = "red") +
  labs(title = paste0("C[1,2] (R_hat = ", round(summary_vals$summary[2,10]), ", n_eff = ", round(summary_vals$summary[2,9]), ")"), x = "Iterations", y="") +  theme_minimal() +
  theme(text = element_text(family = 'Arial'))

p3 <- ggplot(data.frame(Value = extract_vals[,1,3]), aes(x = 1:1000, y = Value)) +
  geom_line() +
  geom_hline(yintercept = C[1,3], linetype = "dashed", color = "red") +
  labs(title = paste0("C[1,3] (R_hat = ", round(summary_vals$summary[3,10]), ", n_eff = ", round(summary_vals$summary[3,9]), ")"), x = "Iterations", y="") +  theme_minimal() +
  theme(text = element_text(family = 'Arial'))

p4 <- ggplot(data.frame(Value = extract_vals[,1,4]), aes(x = 1:1000, y = Value)) +
  geom_line() +
  geom_hline(yintercept = C[1,4], linetype = "dashed", color = "red") +
  labs(title = paste0("C[1,4] (R_hat = ", round(summary_vals$summary[4,10]), ", n_eff = ", round(summary_vals$summary[4,9]), ")"), x = "Iterations", y="") +  theme_minimal() +
  theme(text = element_text(family = 'Arial'))

(trace <- p1 + p2 + p3 + p4)
ggsave("images/trace.png", plot = trace, width = 10, height = 3.5, units = "in", dpi = 300)

# Cosine similarity
C_eigen <- eigen(C)$vectors[,1]
Sigma_eigen <- eigen(pred_C)$vector[,1]
(cos_sim_C_Sigma <- abs(sum(C_eigen * Sigma_eigen) / (sqrt(sum(C_eigen^2)) * sqrt(sum(Sigma_eigen^2)))))

# Marginal scatterplot of input vs output
p1 <- ggplot(data = data.frame(x_obs, y), aes(x = x_obs[,1], y = y)) +
  geom_point() +
  labs(x = expression(theta[i1]), y = expression(Y[i])) +
  theme_bw() +
  theme(text = element_text(family = 'Arial'))

p2 <- ggplot(data = data.frame(x_obs, y), aes(x = x_obs[,2], y = y)) +
  geom_point() +
  labs(x = expression(theta[i2]), y = expression(Y[i])) +
  theme_bw() +
  theme(text = element_text(family = 'Arial'))

p3 <- ggplot(data = data.frame(x_obs, y), aes(x = x_obs[,3], y = y)) +
  geom_point() +
  labs(x = expression(theta[i3]), y = expression(Y[i])) +
  theme_bw() +
  theme(text = element_text(family = 'Arial'))

p4 <- ggplot(data = data.frame(x_obs, y), aes(x = x_obs[,4], y = y)) +
  geom_point() +
  labs(x = expression(theta[i4]), y = expression(Y[i])) +
  theme_bw() +
  theme(text = element_text(family = 'Arial'))

(scatter <- p1 + p2 + p3 + p4 + plot_layout(ncol=4))
ggsave("images/marg_scat.png", plot = scatter, width = 10, height = 3, units = "in", dpi = 300)
