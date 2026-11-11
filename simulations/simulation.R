## Simulation study: smoothness selection for anisotropic tensor-product interactions
##
## Data-generating process: 2-state HMM with normally distributed observations, where
## both off-diagonal transition probabilities depend on two covariates x1 and x2
## through a bivariate parametric surface on the logit scale.
##
## The surfaces are amplitude-modulated waves: each oscillates in one covariate with
## an amplitude that decays in the other. This makes the interaction essential (an
## additive model cannot vary an amplitude) and the anisotropy explicit (the
## oscillating margin needs a small smoothing parameter, the modulating margin a large
## one). The two transition probabilities differ in frequency, modulation, amplitude
## and intercept, so all four smoothing parameters should be distinct. From the
## wiggliness of the marginal shapes, the estimated smoothing parameters should come
## out in the order
##
##   x2 | 2->1 (sin, 1.5 cycles) < x1 | 1->2 (sin, 1 cycle) < x1 | 2->1 (1 - x^3) < x2 | 1->2 (1 - x^2)
##
## A linear modulation was avoided on purpose, as a straight line lies in the null
## space of the second-order difference penalty, such that its smoothing parameter
## would simply run to the upper bound.
##
## Model fitted: the same HMM, with both linear predictors modelled by a tensor-product
## P-spline te(x1, x2), and all four smoothing parameters (one per margin and transition
## probability) selected by qreml(), starting from the mgcv heuristic.
## With k = 10 per margin, the best approximation of the true linear predictors within
## the basis has an RMS error below 0.01 on the logit scale (about 0.04 for k = 8), so
## the approximation error is negligible compared with the estimation error.
##
## The study is repeated for time series of length 2000, 5000 and 10000, to check that
## the estimation error decreases with the sample size.
##
## We use byrow = TRUE throughout, so the first row of beta (and the first column of
## Eta) is the transition 1 -> 2 and the second is 2 -> 1.
##
## Each run is saved as soon as it is fitted, so the simulation can be interrupted and
## resumed. All paths are relative to the project root.


## packages
library(LaMa) # >= 3.1.0
library(scales)
library(parallel)

## color vector for plotting later
color <- LaMaColors(3) # c("orange", "deepskyblue", "seagreen2")

res_dir <- "./simulations/results"
fig_dir <- "./simulations/figs"


# True model --------------------------------------------------------------

N <- 2 # number of states

## state-dependent parameters (well separated -> states are easy to identify)
mu_true <- c(1, 5)      # normal means
sigma_true <- c(1, 1.5) # normal standard deviations

## true linear predictors on the logit scale
# the intercepts imply persistent states (staying probabilities 0.88 and 0.82 where the
# interaction is zero), while switching is more likely than staying only near the wave peaks
# entry (1,2): one full cycle in x1, amplitude decaying quadratically in x2
eta12 <- function(x1, x2) -2.0 + 2.5 * sin(2 * pi * x1) * (1 - x2^2)
# entry (2,1): one and a half cycles in x2, amplitude decaying cubically in x1
eta21 <- function(x1, x2) -1.5 + 2.2 * sin(3 * pi * x2) * (1 - x1^3)

## true transition probability matrices for given covariate values
# handing tpm() a matrix of linear predictors makes it return one t.p.m. per row
make_Gamma <- function(x1, x2) {
  tpm(Eta = cbind(eta12(x1, x2), eta21(x1, x2)), byrow = TRUE)
}

## simulating one data set of length n
simulate_data <- function(n) {
  # covariates: independent and uniform on the unit square
  x1 <- runif(n)
  x2 <- runif(n)

  Gamma <- make_Gamma(x1, x2)
  delta <- stationary(Gamma[, , 1]) # stationary distribution of the first tpm

  # state sequence
  s <- numeric(n)
  s[1] <- sample(1:N, 1, prob = delta)
  for (t in 2:n) s[t] <- sample(1:N, 1, prob = Gamma[s[t - 1], , t])

  # observations
  y <- rnorm(n, mu_true[s], sigma_true[s])

  data.frame(y = y, x1 = x1, x2 = x2, state = s)
}


## prediction grid, used throughout to compare true and estimated surfaces
L <- 51 # grid resolution per margin
x_seq <- seq(0, 1, length = L)
grid <- expand.grid(x1 = x_seq, x2 = x_seq)

Gamma_grid <- make_Gamma(grid$x1, grid$x2)
Pr_true <- cbind(Gamma_grid[1, 2, ], Gamma_grid[2, 1, ]) # one column per entry

## range of the true transition probabilities
round(apply(Pr_true, 2, range), 3)

## visualising the two true surfaces on the probability scale
par(mfrow = c(1, 2))
for (j in 1:2) {
  image(x_seq, x_seq, matrix(Pr_true[, j], L, L),
        col = hcl.colors(100, "YlGnBu", rev = TRUE), zlim = c(0, 1),
        xlab = "x1", ylab = "x2",
        main = paste0("true Pr(", c("1 -> 2", "2 -> 1")[j], ")"))
  contour(x_seq, x_seq, matrix(Pr_true[, j], L, L), add = TRUE, col = "white")
}


# Fitted model ------------------------------------------------------------

## penalised negative log-likelihood function
pnll <- function(par) {
  getAll(par, dat)

  # state-dependent parameters
  REPORT(mu)
  sigma <- exp(logsigma); REPORT(sigma)

  # state process: intercept + tensor-product surface, one row per off-diagonal
  Gamma <- tpm(cbind(beta0, beta_te), Z, byrow = TRUE)

  # state-dependent densities
  allprobs <- statedens(y, dnorm, mu, sigma)

  # forward algorithm, starting in the stationary distribution of the first tpm,
  # + penalty2(), which handles the tensor-product penalties
  -forward("stationary", Gamma, allprobs) +
    penalty2(list(beta_te), S, sp)
}

## fitting one simulated data set
fit_sim <- function(data, k) {
  ## model matrices: anisotropic tensor product of two P-splines
  # te() already contains the two main effects, so no extra s() terms needed
  m <- modmat(~ te(x1, x2, bs = "ps", k = k), data = data)

  ## initial parameters: true values, and zero for the tensor-product coefficients
  # the true interaction terms average to zero over the unit square, so the
  # intercepts of the true linear predictors are the true intercepts in the
  # (centred) tensor-product model
  par <- list(mu = mu_true,                  # state-dependent means
              logsigma = log(sigma_true),    # state-dependent sds
              beta0 = c(-2.0, -1.5),         # state process intercepts
              beta_te = matrix(0, N * (N - 1), m$pardim[["te(x1,x2)"]]))

  ## data and hyperparameter list
  dat <- list(y = data$y, Z = m$Z, S = m$S,
              sp = rep_sp(m, 2)) # initial penalty strengths from mgcv, one per row and margin

  # fitting time excludes building the model matrices
  time <- system.time(
    mod <- try(qreml(pnll, par, dat, random = "beta_te", silent = 2), silent = TRUE)
  )
  if (inherits(mod, "try-error")) return(list(failed = TRUE, error = as.character(mod)))

  # fitted transition probabilities on the prediction grid
  Gamma_hat <- tpm(mod$beta, predict(m, grid), byrow = TRUE)

  list(failed = FALSE,
       Pr_hat = cbind(Gamma_hat[1, 2, ], Gamma_hat[2, 1, ]),
       mu = mod$mu, sigma = mod$sigma,
       sp0 = dat$sp, sp = mod$sp,
       edf0 = mod$edf_hist[, 1], edf = unlist(mod$edf),
       llk = as.numeric(mod$llk), df = attr(logLik(mod), "df"),
       iter = mod$iter, n_fits = mod$n_fits, converged = mod$converged,
       time = time[["elapsed"]])
}


# Running the simulation --------------------------------------------------

nSim <- 100               # number of simulation runs per sample size
ns <- c(2000, 5000, 10000) # lengths of the simulated time series
k <- c(10, 10)            # basis dimension per margin of the tensor product

nCores <- 3 # parallel runs (one fit with n = 10000 needs about 2.2 GB of memory)

## all runs not fitted yet (resume after interruption), largest n first for load balancing
runs_todo <- expand.grid(i = 1:nSim, n = rev(ns))
runs_todo$file <- file.path(res_dir, sprintf("run_n%05d_%03d.rds", runs_todo$n, runs_todo$i))
runs_todo <- runs_todo[!file.exists(runs_todo$file), ]

run_one <- function(r) {
  n <- runs_todo$n[r]; i <- runs_todo$i[r]
  set.seed(n + i) # each run has its own seed, so runs can be fitted in any order
  data <- simulate_data(n)
  res <- c(list(n = n, i = i), fit_sim(data, k))
  saveRDS(res, runs_todo$file[r])

  cat("n =", n, "run", i, "of", nSim, ":",
      if (res$failed) "failed" else paste(round(res$time, 1), "s"), "\n")
  rm(data, res); gc() # RTMB memory piles up across fits otherwise
  NULL
}

# each run in its own forked process, so memory is released after every fit
invisible(mclapply(seq_len(nrow(runs_todo)), run_one,
                   mc.cores = nCores, mc.preschedule = FALSE))


# Results -----------------------------------------------------------------

## reading all runs
runs <- lapply(list.files(res_dir, pattern = "^run_.*\\.rds$", full.names = TRUE), readRDS)
failed <- sapply(runs, `[[`, "failed")
table(n = sapply(runs, `[[`, "n"), failed = failed)
runs <- runs[!failed]

n_run <- sapply(runs, `[[`, "n")
Pr_hat <- simplify2array(lapply(runs, `[[`, "Pr_hat")) # grid x entry x run
sp_hat <- t(sapply(runs, `[[`, "sp"))
edf_hat <- t(sapply(runs, `[[`, "edf"))
edf0 <- t(sapply(runs, `[[`, "edf0"))
mu_hat <- t(sapply(runs, `[[`, "mu"))
sigma_hat <- t(sapply(runs, `[[`, "sigma"))
times <- sapply(runs, `[[`, "time")
iter <- sapply(runs, `[[`, "iter")

lambda_names <- c("x1 | 1->2", "x2 | 1->2", "x1 | 2->1", "x2 | 2->1")
colnames(sp_hat) <- lambda_names


## accuracy of the estimated surfaces, one value per run
# the probability scale is the one we care about, but it compresses differences
# where Pr is near 0 or 1, so we also report the linear predictor (logit) scale
eta_true <- qlogis(Pr_true)
eta_hat <- qlogis(Pr_hat)

rmse_prob <- sqrt(apply((Pr_hat - as.vector(Pr_true))^2, c(2, 3), mean))  # entry x run
rmse_eta <- sqrt(apply((eta_hat - as.vector(eta_true))^2, c(2, 3), mean)) # entry x run

## summary table, one row per sample size
results <- do.call(rbind, lapply(ns, function(n) {
  ind <- n_run == n
  data.frame(n = n,
             "failed fits" = nSim - sum(ind),
             "median time (s)" = median(times[ind]),
             "median iterations" = median(iter[ind]),
             "RMSE prob 1->2" = mean(rmse_prob[1, ind]),
             "RMSE prob 2->1" = mean(rmse_prob[2, ind]),
             "RMSE eta 1->2" = mean(rmse_eta[1, ind]),
             "RMSE eta 2->1" = mean(rmse_eta[2, ind]),
             # median, as the edf of single runs can be broken by a nearly singular Hessian
             "median edf 1->2" = median(edf_hat[ind, 1]),
             "median edf 2->1" = median(edf_hat[ind, 2]),
             check.names = FALSE)
}))
print(results, digits = 3, row.names = FALSE)
write.csv(results, file.path(res_dir, "summary.csv"), row.names = FALSE)

## initial edf implied by the mgcv heuristic
apply(edf0, 2, range)


## state-dependent parameters
par_fit <- NULL
for (n in ns) {
  ind <- n_run == n
  for (j in 1:N) {
    for (nm in c("mu", "sigma")) {
      est <- if (nm == "mu") mu_hat[ind, j] else sigma_hat[ind, j]
      tru <- if (nm == "mu") mu_true[j] else sigma_true[j]
      par_fit <- rbind(par_fit, data.frame(
        n = n, par = paste0(nm, j), true = tru,
        mean = mean(est), bias = mean(est) - tru, sd = sd(est),
        rmse = sqrt(mean((est - tru)^2))))
    }
  }
}
print(par_fit, digits = 3, row.names = FALSE)
write.csv(par_fit, file.path(res_dir, "par_fit.csv"), row.names = FALSE)


## did the anisotropy come out as predicted?
predicted <- c("x2 | 2->1", "x1 | 1->2", "x1 | 2->1", "x2 | 1->2")
for (n in ns) {
  med <- apply(sp_hat[n_run == n, ], 2, median)
  cat("n =", n, "median sp:", round(med, 2), "\n")
  cat("  predicted order matched:", identical(predicted, names(sort(med))), "\n")
}


## accuracy against the sample size
par(mfrow = c(1, 2))
for (j in 1:2) {
  boxplot(rmse_eta[j, ] ~ n_run, xlab = "T", ylab = "RMSE (logit scale)",
          main = paste0("Pr(", c("1 -> 2", "2 -> 1")[j], ")"),
          col = alpha(color[j], 0.5), bty = "n")
}


## slices through the surfaces: all individual fits vs. the truth
# each entry is sliced along its oscillating direction, once where the amplitude
# is at its largest and once where the modulation has damped it down, so that
# the interaction itself is visible
at <- function(v) x_seq[which.min(abs(x_seq - v))] # nearest grid value

# finer grid for drawing only: the truth is evaluated on it directly, the fitted
# slices (stored on the coarser evaluation grid) are interpolated by a cubic spline,
# which is visually exact for these smooth P-spline fits
L_plot <- 150
x_plot <- seq(0, 1, length = L_plot)
smooth_line <- function(y) spline(x_seq, y, xout = x_plot, method = "fmm")$y

pdf(file.path(fig_dir, "sim_slices.pdf"), width = 8, height = 9)
par(mfrow = c(4, 3), mar = c(4, 4.5, 2.5, 0.5))
for (j in 1:2) {
  for (v in c(0, 0.75)) {
    # entry 1 oscillates in x1 and is modulated by x2, entry 2 the other way round
    ind <- if (j == 1) which(grid$x2 == at(v)) else which(grid$x1 == at(v))
    for (n in ns) {
      plot(NA, xlim = c(0, 1), ylim = c(0, 1), bty = "n",
           xlab = if (j == 1) expression(x[1]) else expression(x[2]),
           ylab = if (j == 1) expression(gamma[12]) else expression(gamma[21]),
           main = if (j == 1) bquote(list(T == .(n), x[2] == .(at(v))))
                  else bquote(list(T == .(n), x[1] == .(at(v)))))
      for (r in which(n_run == n)) {
        lines(x_plot, smooth_line(Pr_hat[ind, j, r]), col = alpha(color[j], 0.15))
      }
      Gamma_line <- if (j == 1) make_Gamma(x_plot, at(v)) else make_Gamma(at(v), x_plot)
      lines(x_plot, if (j == 1) Gamma_line[1, 2, ] else Gamma_line[2, 1, ], lwd = 2)
    }
  }
}
dev.off()


## runtime and number of outer iterations
par(mfrow = c(1, 2))
boxplot(times ~ n_run, xlab = "T", ylab = "runtime per fit (s)",
        col = alpha(color[1], 0.5), bty = "n")
boxplot(iter ~ n_run, xlab = "T", ylab = "outer iterations",
        col = alpha(color[1], 0.5), bty = "n")
