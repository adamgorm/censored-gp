# Exact censoring model: density/Jacobian identity, bounds, and HMC smoke tests.
# Run from the repository root: Rscript research/validate-censoring.R
library(cmdstanr)
output <- tempfile("censored-gp-censor-validation-")
dir.create(output)
file.copy("research/censored-gp.stan", output)
file.copy("research/censor-transform.staninc", output)
model <- cmdstan_model(file.path(output, "censored-gp.stan"),
  stanc_options = list("allow-undefined" = TRUE),
  user_header = normalizePath("research/cpp/gaussian_stan.hpp"),
  quiet = TRUE)
# CmdStanR's standalone wrappers cannot expose externally declared functions.
# Compile just the shared Stan transformation to validate its Jacobian.
transform_file <- file.path(output, "censor-transform.stan")
writeLines(c("functions {", readLines("research/censor-transform.staninc"), "}"), transform_file)
transform_model <- cmdstan_model(transform_file, compile_standalone = TRUE, quiet = TRUE)
transform <- transform_model$functions$censor_transform
set.seed(463)
A <- matrix(rnorm(25),5,5)
K <- tcrossprod(A) + diag(.5,5)
L <- t(chol(K)); yo <- c(-.3,.7)
lower <- c(-Inf,-.2,-.5); upper <- c(.3,Inf,.8)
kind <- c(1L,2L,3L); u <- c(.37,.61,.44)
result <- transform(L,yo,u,lower,upper,kind)
J <- numDeriv::jacobian(function(v) transform(L,yo,v,lower,upper,kind)[[1]],u)
lhs <- mvtnorm::dmvnorm(c(yo,result[[1]]),sigma=K,log=TRUE) + sum(log(abs(diag(J))))
rhs <- mvtnorm::dmvnorm(yo,sigma=K[1:2,1:2],log=TRUE) + result[[2]]
stopifnot(abs(lhs-rhs) < 1e-7)
cat("Censoring Jacobian log-density error:", abs(lhs-rhs), "\n")
for(bounds in list(c(-40,-39.9),c(39.9,40),c(-1e-7,1e-7))) {
  rr <- transform(matrix(1,1,1),numeric(),.4,bounds[1],bounds[2],3L)
  stopifnot(is.finite(rr[[2]]),rr[[1]]>bounds[1],rr[[1]]<bounds[2])
}
# Mixed censoring with observation-specific limits, including no uncensored data.
for(No in c(0L,3L)) for(parameterization in 0:1) {
  dat <- list(N=No+3L,G=2L,No=No,Nl=1L,Nr=1L,Ni=1L,
    parameterization=parameterization, custom_gradient=1L,
    yo=rep(.1,No),left_upper=-.1,right_lower=.2,
    interval_lower=-.3,interval_upper=.4,
    t=runif(No+3L),group=rep(1:2,length.out=No+3L))
  init <- list(alpha_mu=1,rho_mu=.3,alpha_eta=.5,rho_eta=.2,sigma=.3,
    yc_left=if(parameterization==0) -.3 else numeric(),
    yc_right=if(parameterization==0) .4 else numeric(),
    yc_interval=if(parameterization==0) .1 else numeric(),
    u=if(parameterization==1) rep(.5,3) else numeric())
  fit <- model$sample(data=dat,init=list(init),chains=1,seed=617,
    iter_warmup=75,iter_sampling=50,refresh=0,output_dir=output)
  draws <- fit$draws("y_cens", format="matrix")
  stopifnot(all(draws[,1] <= dat$left_upper),
    all(draws[,2] >= dat$right_lower),
    all(draws[,3] >= dat$interval_lower),all(draws[,3] <= dat$interval_upper))
}
cat("PASS: mixed, tail, narrow-interval, and entirely censored cases\n")
