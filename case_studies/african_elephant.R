## Ivory coast elephant case study
# The smooth models are saved and loaded by default. Set refit = TRUE to fit them
# again (a few minutes) and to recompute the Monte Carlo quantiles of the stationary distribution.
refit = FALSE
# Set redraw = TRUE to save the paper figures to ./case_studies/figs/.
redraw = FALSE

## packages
library(LaMa) # >= 3.1.0


## data
data = read.csv("./data/elephant_data.csv")
data$tod = 2 * data$tod # 2-hourly time of day in 2, 4, ..., 24
data$step[data$step == 0] = NA # two zero step lengths, which have zero density under the gamma distribution
head(data)

# number of observations
nrow(data)


## color vector for plotting later
color <- LaMaColors(2)



# Homogeneous model -------------------------------------------------------

# likelihood function for homogeneous model
nll_hom = function(par){
  getAll(par, dat)
  ## transforming parameters
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  Gamma = tpm(eta)
  ## state-dependent densities
  allprobs = statedens(step, dgamma2, mu, sigma) * statedens(angle, dvm, 0, kappa)
  ## forward algorithm, starting in the stationary distribution
  -forward("stationary", Gamma, allprobs)
}

# initial parameter list
par = list(logmu = log(c(0.2, 2)),
           logsigma = log(c(0.2, 2)),
           logkappa = log(c(0.2, 1)),
           eta = rep(-2, 2))

# data list
dat = list(step = data$step, angle = data$angle)

# fitting the model
obj_hom = MakeADFun(nll_hom, par, silent = TRUE)
opt_hom = nlminb(obj_hom$par, obj_hom$fn, obj_hom$gr)
mod_hom = lama(obj_hom)
summary(mod_hom)

# decoding states
mod_hom$states = viterbi(mod_hom)

# state-dependent densities weighted by the local state probabilities
plot(mod_hom, type = "statedens")


## plotting state-dependent distributions (paper figure)
mu = mod_hom$mu
sigma = mod_hom$sigma
kappa = mod_hom$kappa
delta = mod_hom$delta

if (redraw) pdf("./case_studies/figs/elephant_marginal.pdf", width = 8, height = 4)

par(mfrow = c(1,2))
hist(data$step, breaks = 100, prob = T, bor = "white", xlim = c(0,5), main = "", xlab = "step length", ylab = "density")
curve(delta[1]*dgamma2(x, mu[1], sigma[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dgamma2(x, mu[2], sigma[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[1]*dgamma2(x, mu[1], sigma[1])+delta[2]*dgamma2(x, mu[2], sigma[2]),
      add = T, lwd = 2, lty = 2, n = 500)
legend("topright", legend = c("encamped", "exploratory", "marginal"), col = c(color[1], color[2], "black"),
       lty = c(1,1,2), bty = "n", lwd = 2)

hist(data$angle, breaks = 20, prob = T, bor = "white", main = "", xlab = "turning angle", ylab = "density")
curve(delta[1]*dvm(x, 0, kappa[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dvm(x, 0, kappa[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[1]*dvm(x, 0, kappa[1])+delta[2]*dvm(x, 0, kappa[2]),
      add = T, lwd = 2, lty = 2, n = 500)

if (redraw) dev.off()



# Smooth models -----------------------------------------------------------

# All smooth models depend on time of day and/ or day of year only, so the
# transition probabilities are computed once for each of the 366 * 12 unique
# combinations, i.e. the model matrices are built on this grid. Indexing them by idx
# makes the state process periodic over one year, such that
# forward("stationary", ..., tod = idx) starts in the periodically stationary distribution.
grid = expand.grid(tod = 1:12 * 2, doy = 1:366)
data$idx = (data$doy - 1) * 12 + data$tod / 2 # row of grid for each observation

## penalised likelihood function for all smooth models
# the spline coefficients are those named in dat$random, one matrix per smooth
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

## fitting a smooth model given the formula for the state process
# names of the coefficient matrices for each smooth term
coef_names = c("s(tod)" = "beta_tod", "s(doy)" = "beta_doy", "ti(tod,doy)" = "beta_ti", "te(tod,doy)" = "beta_te")

fit_smooth = function(formula, silent = 0) {
  m <- modmat(formula, data = grid,
              knots = list(tod = c(0, 24), doy = c(0, 366))) # where to wrap the cyclic bases
  random = unname(coef_names[names(m$S)])

  par = c(list(logmu = log(c(0.3, 1.1)),      # state-dependent step mean
               logsigma = log(c(0.25, 0.75)), # state-dependent step sd
               logkappa = log(c(0.2, 0.7)),   # state-dependent angle concentration
               beta0 = c(-2, -2)),            # state process intercept
          setNames(lapply(m$pardim[names(m$S)], function(p) matrix(0, 2, p)), random)) # spline coefficients

  dat = list(step = data$step, angle = data$angle, idx = data$idx,
             Z = m$Z, S = m$S, random = random,
             sp = rep_sp(m, 2)) # initial penalty strengths from mgcv, one per row and penalty

  # fitting time excludes building the model matrices
  time = system.time(mod <- qreml(pnll, par, dat, random = random, silent = silent))
  mod$time = time
  mod$modmat = m # for predictions later
  mod
}

k = 12

if (refit) {
  ## time of day only
  mod_tday <- fit_smooth(~ s(tod, bs = "cc", k = k))
  print(mod_tday$time)

  ## day of year only
  mod_julian <- fit_smooth(~ s(doy, bs = "cc", k = k))
  print(mod_julian$time)

  ## additive model
  mod_add <- fit_smooth(~ s(tod, bs = "cc", k = k) + s(doy, bs = "cc", k = k))
  print(mod_add$time)

  ## tensor-product interaction model (poorly identified, so not used)
  # mod_ti <- fit_smooth(~ s(tod, bs = "cc", k = k) + s(doy, bs = "cc", k = k) +
  #                        ti(tod, doy, bs = "cc", k = k))
  # mod_ti$time
  #
  # summary(mod_ti)

  ## full tensor-product model, main effects and interaction in one smooth
  mod_te <- fit_smooth(~ te(tod, doy, bs = "cc", k = k))
  print(mod_te$time)

  print(summary(mod_te))

  # outer Hessian for uncertainty in the penalty strengths, evaluated now
  # as outer_hessian() does not work any more after saving and loading the model
  mod_te$H_outer = mod_te$outer_hessian()
  # joint precision of the coefficients and log(sp), for intervals that also account
  # for the uncertainty in sp (first-order correction of Wood, Pya and Saefken, 2016)
  mod_te$Q_joint = mod_te$joint_uncertainty(H_rho = mod_te$H_outer)

  mods = list(hom = mod_hom, tday = mod_tday, julian = mod_julian, add = mod_add, te = mod_te)
  saveRDS(mods, "./case_studies/objects/elephant_mods_v4.rds")
}

## loading saved models
mods = readRDS("./case_studies/objects/elephant_mods_v4.rds")
for (nm in names(mods)) assign(paste0("mod_", nm), mods[[nm]])



# Model comparison --------------------------------------------------------

IC_table = data.frame(llk = sapply(mods, logLik),
                      AIC = sapply(mods, AIC),
                      BIC = sapply(mods, BIC))
round(IC_table, 2)

# differences to the best model
diff_IC_table = IC_table
diff_IC_table$AIC = IC_table$AIC - min(IC_table$AIC)
diff_IC_table$BIC = IC_table$BIC - min(IC_table$BIC)
round(diff_IC_table, 2)



# Model checking ----------------------------------------------------------

# pseudo-residuals, one row per data stream
plot(mod_add)
plot(mod_te)

# state-dependent densities weighted by the local state probabilities
plot(mod_te, type = "statedens")

# uncertainty in the estimated penalty strengths, from the Hessian of the
# negative restricted log-likelihood with respect to log(sp)
rbind(sp = mod_te$sp, se = sqrt(diag(solve(mod_te$H_outer))) * mod_te$sp) # delta method



# Periodically stationary distribution over the year ----------------------

modmat_te = mod_te$modmat

# The periodically stationary distribution at time of day t is computed from the 12
# transition probability matrices at t + 2, t + 4, ..., t + 24, the transitions into
# these times. Each such cycle gives the distribution at t, t + 2, ..., t + 22,
# so 16 offsets t in [0, 2) give 192 time points per day.
offsets = seq(0, 2, length = 17)[-17]
tod_cycle = c(outer(1:12 * 2, offsets, "+")) %% 24 # 12 times per offset
tod_plot = c(outer(0:11 * 2, offsets, "+"))        # times of the distributions
ord = order(tod_plot)
tod_seq = tod_plot[ord]

# periodically stationary distribution of state 2 on one day, for given beta
Delta_day = function(beta, Z_day) {
  Delta2 = sapply(seq_along(offsets), function(o) {
    Gamma = tpm(beta, Z_day[(o - 1) * 12 + 1:12, ])
    stationary_p(Gamma)[, 2]
  })
  c(Delta2)[ord]
}

days = c(1, 61, 121, 182, 243, 304) # days to plot
daynames = paste(c("Jan", "March", "May", "Jul", "September", "November"), 1)
Z_days = lapply(days, function(day) predict(modmat_te, data.frame(tod = tod_cycle, doy = day)))

# point estimate
Delta = sapply(Z_days, Delta_day, beta = mod_te$beta)

# uncertainty: Monte Carlo samples of the coefficients and log(sp) from their joint
# approximate posterior, so that the intervals also account for the uncertainty in sp;
# the log(sp) part of each sample is dropped
if (refit) {
  set.seed(123)
  samples = MCreport(mod_te, nSamples = 2000, Q = mod_te$Q_joint) # or joint = TRUE
  betas = Map(cbind, samples$beta0, samples$beta_te)
  Delta_samples = sapply(betas, function(beta) sapply(Z_days, Delta_day, beta = beta),
                         simplify = "array") # times x days x samples
  Delta_q = apply(Delta_samples, 1:2, quantile, probs = c(0.025, 0.975))
  saveRDS(Delta_q, file = "./case_studies/objects/elephant_delta_quantiles_v4.rds")
}

Delta_q = readRDS("./case_studies/objects/elephant_delta_quantiles_v4.rds")


## plotting the estimated state distribution over the year

if (redraw) pdf("./case_studies/figs/elephant_stationary.pdf", width = 7, height = 5)

par(mfrow = c(2,3))
source("./utils/gen_sun_colors.R")
sun_cycle_colors <- gen_sun_cycle_colors()
for(d in seq_along(days)){
  plot(NA, xlim = c(0, 24), ylim = c(0,1), bty = "n",
       xlab = "Time of day", ylab = "Pr(exploratory)",
       main = daynames[d])
  polygon(x = c(0, 6.5, 6.5, 0), y = c(0, 0, 1, 1), col = "gray95", border = "white")
  polygon(x = c(18.5, 24, 24, 18.5), y = c(0, 0, 1, 1), col = "gray95", border = "white")
  for(t in 0:47){
    polygon(x = c(t/2, (t+1)/2, (t+1)/2, t/2), y = c(-0.04, -0.04, -0.01, -0.01), col = sun_cycle_colors[t+1], border = sun_cycle_colors[t+1])
  }
  polygon(c(tod_seq, rev(tod_seq)),
          c(Delta_q[1,,d], rev(Delta_q[2,,d])), col = "#00000020", border = NA)
  lines(tod_seq, Delta[,d], lwd = 2)
}

if (redraw) dev.off()



# More plots: visualising the fitted surfaces in 3D -----------------------

library(plotly)

tod = seq(0, 24, length = 70)
doy = seq(0, 366, length = 70)
Z_surf = predict(modmat_te, expand.grid(tod = tod, doy = doy))
Gamma_surf = tpm(mod_te$beta, Z_surf)

gamma_12 = matrix(Gamma_surf[1,2,], length(tod), length(doy))
gamma_21 = matrix(Gamma_surf[2,1,], length(tod), length(doy))

# plotly expects z with rows corresponding to y
plot_ly(x = tod, y = doy, z = t(gamma_12), type = "surface") %>%
  layout(title = "Pr(1 -> 2)",
         scene = list(xaxis = list(title = "Time of day"),
                      yaxis = list(title = "Day of year"),
                      zaxis = list(title = "Pr(1 -> 2)")))

plot_ly(x = tod, y = doy, z = t(gamma_21), type = "surface") %>%
  layout(title = "Pr(2 -> 1)",
         scene = list(xaxis = list(title = "Time of day"),
                      yaxis = list(title = "Day of year"),
                      zaxis = list(title = "Pr(2 -> 1)")))
