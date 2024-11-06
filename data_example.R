library(tidyverse)

### Setup data and Stan functions

## Read in and prepare data
d <- read.delim("A5055data.txt")
d <- tibble(id = d$Subject, t = d$day, rna = d$rna)
d$y <- log10(d$rna)
d$y[d$rna < 50] <- log10(50)
d$censored <- d$rna < 50
d <- arrange(d, id)
n <- length(unique(d$id))
all(unique(d$id) == 1:max(d$id)) # check that all numbers appear as id
## Flip y to turn left censoring into right censoring
d <- mutate(d, y = -y)

## Expose Stan functions to R
stan_funs <- cmdstanr::cmdstan_model("functions.stan",
                                     force_recompile = TRUE,
                                     compile_standalone = TRUE)$functions

### Fit our model on the data

do <- filter(d, !censored)
dc <- filter(d, censored)

yo <- do$y
to <- do$t

yc <- dc$y
tc <- dc$t

tc <- c()
Jc <- c()
to <- c()
Jo <- c()
for (i in 1:max(d$id)) {
    d |>
        filter(id == i, !censored) |>
        pull(t) ->
        to_i
    to <- c(to, to_i)
    Jo <- c(Jo, length(to_i))
    d |>
        filter(id == i, censored) |>
        pull(t) ->
        tc_i
    tc <- c(tc, tc_i)
    Jc <- c(Jc, length(tc_i))
}
tc_is <- c(0, cumsum(Jc))
to_is <- c(0, cumsum(Jo))

## params: c(magnitude_mu, length_scale_mu, magnitude_eta, length_scale_eta, sigma)
marg_log_lik <- function(params)
{
    prep_list <- stan_funs$prep_multi_cens_log_lik(yo,
                                                   to, tc,
                                                   to_is, tc_is,
                                                   Jo, Jc,
                                                   params[[1]], params[[2]],
                                                   params[[3]], params[[4]],
                                                   params[[5]])
    prep_list[[1]] + log(mvtnorm::pmvnorm(lower = yc,
                                          mean = prep_list[[2]],
                                          sigma = prep_list[[3]]))
}

neg_marg_log_lik <- function(params) -marg_log_lik(params)

set.seed(1627)
deopt_res <- DEoptim::DEoptim(neg_marg_log_lik,
                              c(0.1, 0.1, 0.1, 0.1, 0.1),
                              c(10, 200, 10, 200, 10))

param_max <- deopt_res$optim$bestmem

tpred <- unique(sort(d$t))
Jpred <- length(tpred)

post_list <-
    stan_funs$prep_multi_cond_post(
                  yo, yc,
                  to, tc,
                  tpred,
                  to_is, tc_is,
                  Jo, Jc,
                  Jpred,
                  param_max[1], param_max[2],
                  param_max[3], param_max[4],
                  param_max[5]
              )

names(post_list) <- c(
    "mean_c_cond_o",
    "cov_c_cond_o",
    "mean_mueta_cond_o",
    "p_factor_mueta",
    "cov_q_mueta"
)

n_draw <- 100

Npred_except_1_group <- (n-1)*Jpred
Npred <- n * Jpred

q_draws <- t(mvtnorm::rmvnorm(n_draw, mean = rep(0, Npred),
                              sigma = post_list$cov_q_mueta))

p_draws <- t(TruncatedNormal::rtmvnorm(n = n_draw,
                                       mu = rep(0, nrow(post_list$cov_c_cond_o)),
                                       sigma = post_list$cov_c_cond_o,
                                       lb = yc - post_list$mean_c_cond_o))

mueta_draws <- post_list$mean_mueta_cond_o +
    post_list$p_factor_mueta %*% p_draws +
    q_draws

mu_draws <- mueta_draws[1:Jpred, ]
eta_draws <- mueta_draws[(Jpred+1):Npred, ]
eta_draws <- array(eta_draws, dim = c(Jpred, n-1, n_draw))
eta_draws_last <- -apply(eta_draws, c(1, 3), sum)
eta_draws_full <- abind::abind(eta_draws, eta_draws_last, along = 2)

## Flip back data

d <- mutate(d, id_char = as.character(id), y = -y)

## Calculate posterior means and intervals (multiplied by -1 to flip around
## again)

post_mean_mu <- -rowMeans(mu_draws)
post_ci_mu <- -apply(mu_draws, 1, \(x) quantile(x, c(0.05, 0.95)))
post_means_eta <- -apply(eta_draws_full, c(1,2), mean)
post_ci_eta <- -apply(eta_draws_full, c(1,2), \(x) quantile(x, c(0.05, 0.95)))

f_draws <- array(NA, dim(eta_draws_full))
for (i in 1:n) {
    f_draws[, i, ] <- mu_draws
}
f_draws <- f_draws + eta_draws_full
post_means_f <- -apply(f_draws, c(1,2), mean)
post_ci_f <- -apply(f_draws, c(1,2), \(x) quantile(x, c(0.05, 0.95)))

## Plot results

censoring_val <- log10(50)
ggplot() +
    geom_hline(aes(yintercept = censoring_val), lty = 1, col = 8) +
    geom_line(aes(x = t, y = y, group = id_char), d, col = 1) +
    geom_point(aes(x = t, y = y, shape = censored, group = id_char), d, col = 1, size = 2) +
    scale_shape_manual(values = c(16, 3),
                       breaks = c(FALSE, TRUE),
                       labels = c("Observed", "Censored")) +
    geom_line(aes(x = tpred, y = post_mean_mu), col = 2) +
    geom_ribbon(aes(x = tpred, ymin = post_ci_mu[1,], ymax = post_ci_mu[2,]),
                alpha = 0.2, fill = 2) +
    scale_linetype_manual(values = c(1, 2)) +
    scale_fill_manual(values = c(1, 8)) +
    labs(x = "Study time [days]",
         y = expression(log[10]*"(HIV-1 RNS)")) +
    theme(legend.position = "bottom") +
    ggtitle(expression("Data and posterior of " * mu))

f_id <- 1
ggplot() +
    geom_hline(aes(yintercept = censoring_val), lty = 1, col = 8) +
    geom_line(aes(x = t, y = y, group = id_char),
              filter(d, id == f_id), col = 1) +
    geom_point(aes(x = t, y = y, shape = censored, group = id_char),
               filter(d, id == f_id), col = 1, size = 3) +
    geom_line(aes(x = tpred, y = post_means_f[,f_id]), col = 2) +
    geom_ribbon(aes(x = tpred,
                    ymin = post_ci_f[1,,f_id],
                    ymax = post_ci_f[2,,f_id]),
                fill = 2, alpha = 0.2) +
    labs(x = "Study time [days]",
         y = expression(log[10]*"(HIV-1 RNS)")) +
    theme(legend.position = "bottom") +
    ggtitle(expression("Data and posterior of " * f[1]))    
