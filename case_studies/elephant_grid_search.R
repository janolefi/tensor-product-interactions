## Grid search for the smoothing parameters of the elephant te() model
# Baseline for qREML (R2): each grid point is one penalized fit at fixed smoothing
# parameters, started from the same initial values as qREML (no warm starts), and is
# scored by the same restricted log-likelihood that qREML maximizes.
# The fits are saved one by one, so the grid search can be interrupted and resumed.

## packages
library(LaMa) # >= 3.1.0


## data (as in african_elephant.R)
data = read.csv("./data/elephant_data.csv")
data$tod = 2 * data$tod # 2-hourly time of day in 2, 4, ..., 24
data$step[data$step == 0] = NA # two zero step lengths

# transition probabilities on the grid of unique (tod, doy) combinations
grid = expand.grid(tod = 1:12 * 2, doy = 1:366)
data$idx = (data$doy - 1) * 12 + data$tod / 2 # row of grid for each observation


## penalised likelihood function (as in african_elephant.R)
pnll = function(par) {
  getAll(par, dat)

  # state-dependent parameters
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)

  # state process
  beta_spline = mget(random) # list of spline coefficient matrices
  Gamma = tpm(do.call(cbind, c(list(beta0), beta_spline)), Z) # one tpm per row of grid

  # state-dependent densities
  allprobs = statedens(step, dgamma2, mu, sigma) * statedens(angle, dvm, 0, kappa)

  # forward algorithm + quadratic penalty
  -forward("stationary", Gamma, allprobs, tod = idx) +
    penalty2(beta_spline, S, sp)
}


## model matrices, initial parameters and data for the te() model
k = 12
m = modmat(~ te(tod, doy, bs = "cc", k = k), data = grid,
           knots = list(tod = c(0, 24), doy = c(0, 366)))
random = "beta_te"

par = list(logmu = log(c(0.3, 1.1)),
           logsigma = log(c(0.25, 0.75)),
           logkappa = log(c(0.2, 0.7)),
           beta0 = c(-2, -2),
           beta_te = matrix(0, 2, m$pardim[["te(tod,doy)"]]))

dat = list(step = data$step, angle = data$angle, idx = data$idx,
           Z = m$Z, S = m$S, random = random)


## grid of smoothing parameters
# four smoothing parameters: tod and doy margins for both off-diagonal entries of the t.p.m.
# The same coarse grid is used for all four, equidistant on the log scale. Its range is
# informed by the qREML estimates (0.031, 23.1, 0.019, 8.08), i.e. it is a best case for
# grid search, but the grid points do not contain the estimates themselves.
lsp_grid = seq(log(10^-2.5), log(10^2.5), length.out = 5)
sp_grid = as.matrix(exp(expand.grid(tod.1 = lsp_grid, doy.1 = lsp_grid,
                                    tod.2 = lsp_grid, doy.2 = lsp_grid)))
nrow(sp_grid) # 625 penalized fits


## fitting all grid points
out_dir = "./case_studies/objects/elephant_grid"
dir.create(out_dir, showWarnings = FALSE)

fit_fixed = function(sp) {
  dat$sp = sp
  # all smoothing parameters fixed: a single penalized fit, scored by the restricted log-likelihood
  time = system.time(
    mod <- qreml(pnll, par, dat, random = random,
                 map = list(sp = factor(rep(NA, length(sp)))), silent = 2)
  )
  # conditional AIC and BIC, based on the effective number of parameters (as in Table 1)
  data.frame(t(sp), llk_restricted = mod$llk_restricted, llk = as.numeric(mod$llk),
             df = attr(logLik(mod), "df"), AIC = AIC(mod), BIC = BIC(mod),
             edf.1 = mod$edf[[1]][1], edf.2 = mod$edf[[1]][2], time = time[["elapsed"]])
}

for (i in seq_len(nrow(sp_grid))) {
  file = file.path(out_dir, sprintf("fit_%03d.rds", i))
  if (file.exists(file)) next # resume after interruption
  res = fit_fixed(sp_grid[i, ])
  saveRDS(res, file)
  cat(i, "/", nrow(sp_grid), ": l_r =", round(res$llk_restricted, 2),
      "time =", round(res$time, 1), "s\n")
  rm(res); gc() # RTMB memory piles up across fits otherwise
}


## results
grid_res = do.call(rbind, lapply(list.files(out_dir, full.names = TRUE), readRDS))

# best grid point by restricted log-likelihood, AIC and BIC
best = grid_res[which.max(grid_res$llk_restricted), ]
best
grid_res[which.min(grid_res$AIC), ]
grid_res[which.min(grid_res$BIC), ]

# comparison with qREML (te model refitted in african_elephant.R)
# mod_te = readRDS("./case_studies/objects/elephant_mods_v3.rds")$te
# data.frame(method = c("qREML", "grid search"),
#            fits = c(mod_te$n_fits, nrow(grid_res)),
#            time = c(mod_te$time[["elapsed"]], sum(grid_res$time)),
#            llk_restricted = c(max(mod_te$llk_restricted[mod_te$llk_restricted < 0]), best$llk_restricted),
#            llk = c(as.numeric(mod_te$llk), best$llk))
