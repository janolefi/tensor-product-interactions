## Fruit fly case study
# The models are saved and loaded by default. Set refit = TRUE to fit them again
# (a few minutes each).
refit = FALSE
# Set redraw = TRUE to save the paper figures to ./case_studies/figs/.
redraw = FALSE

## packages
library(LaMa) # >= 3.1.0
library(scales)
library(RTMBdist) # for dnbinom2 (mean-size parameterisation)


## loading the data: wild type and modified genotype
data <- list(
  wild = readRDS("./data/fruitflies_genotype0.rds"),
  modified = readRDS("./data/fruitflies_genotype1.rds")
)

## number of observations
sapply(data, nrow)

## ranges
lapply(data, function(d) range(d$activity, na.rm = TRUE))

## unique animal IDs
aniIDs <- lapply(data, function(d) unique(d$aniID))
(nAnimals <- sapply(aniIDs, length))


## color vector for plotting later
color <- LaMaColors(2) # c("orange", "deepskyblue")



# Model fitting -----------------------------------------------------------

## penalised log-likelihood function for the tensor-product model
pnll <- function(par) {
  getAll(par, dat) # extract everything from parameter and data list

  # state-dependent process parameters (mean and dispersion)
  mu <- exp(logmu); REPORT(mu)
  phi <- exp(logphi); REPORT(phi)

  # state process: intercepts (LD, DD) and random intercepts, smooth time of day effects
  # and interactions of time of day and individual for LD and DD, in the order of dat$random
  beta <- cbind(beta0, do.call(cbind, mget(random)))
  Gamma <- tpm(beta, Z)

  # periodically stationary initial distribution for each track, from the 48
  # transition probability matrices into its first 48 time points
  Delta <- AD(matrix(0, length(trackInd), 2))
  for (tr in seq_along(trackInd)) Delta[tr, ] <- stationary_p(Gamma[, , trackInd[tr] + 0:47], t = 2)

  # state-dependent probabilities (negative binomial with mean mu and dispersion phi)
  allprobs <- statedens(activity, dnbinom2, mu, 1 / phi)

  # forward algorithm + penalty
  -forward(Delta, Gamma, allprobs, trackID) +
    penalty2(mget(random), S, sp)
}

## names of the coefficient matrices, e.g. "s(aniID):conditionLD" -> "beta_riLD"
coef_name <- function(term) {
  smooth <- c("s(aniID)" = "ri", "s(tod)" = "time", "ti(aniID,tod)" = "ti")[sub(":condition.*", "", term)]
  paste0("beta_", smooth, sub(".*:condition", "", term))
}

## fitting the model to the data of one genotype
fit_genotype <- function(d, k = 15, silent = 0) {
  nA <- length(unique(d$aniID)) # inserted into the formula, as mgcv cannot see local variables
  m <- modmat(
    eval(bquote(~ condition +
      s(aniID, bs = "re", by = condition) +
      s(tod, bs = "cc", by = condition, k = .(k)) + # same basis dimension for tod as in the interaction
      ti(aniID, tod, bs = c("re", "cc"), by = condition, k = c(.(nA), .(k))))),
    data = d,
    knots = list(tod = c(0, 24))
  )
  random <- coef_name(names(m$S))

  # initial parameters
  par <- c(
    list(logmu = log(c(4, 55)),     # state-dependent mean
         logphi = log(c(10, 0.5)),  # state-dependent dispersion
         beta0 = matrix(-2, 2, 2)), # state process intercepts
    setNames(lapply(m$pardim[names(m$S)], function(p) matrix(0, 2, p)), random) # spline coefficients
  )

  # data and hyperparameter list
  dat <- list(
    activity = d$activity, # observations
    trackID = d$trackID, # track ID
    trackInd = calc_trackInd(d$trackID), # index of the first observation of each track (for initial distributions)
    Z = m$Z, # design matrix
    S = m$S, # penalty matrices for s(ID), s(time) and lists of penalty matrices for ti(ID, time)
    random = random,
    sp = rep_sp(m, 2) # initial penalty strengths from mgcv: one per row and penalty, all 16 estimated freely
  )

  # fitting time excludes building the model matrices
  time <- system.time(
    mod <- qreml(pnll, par, dat,
      random = random,
      silent = silent,
      spHess = TRUE # exploits sparsity pattern in Hessian -> much faster
    )
  )
  mod$time <- time
  mod$modmat <- m # for predictions later
  # outer Hessian for uncertainty in the penalty strengths, evaluated now
  # as outer_hessian() does not work any more after saving and loading the model
  mod$H_outer <- mod$outer_hessian()
  mod
}

if (refit) {
  ## fitting both models with qREML (a few minutes each)
  mod_wild <- fit_genotype(data$wild)
  print(mod_wild$time)
  mod_modified <- fit_genotype(data$modified)
  print(mod_modified$time)

  mods <- list(wild = mod_wild, modified = mod_modified)
  saveRDS(mods, "./case_studies/objects/fruitflies_mods_v4.rds")
}

## loading saved models
mods <- readRDS("./case_studies/objects/fruitflies_mods_v4.rds")

lapply(mods, summary)

# uncertainty in the estimated variance parameters 1 / sp, from the Hessian of the
# negative restricted log-likelihood with respect to log(sp)
lapply(mods, function(m) {
  sp <- m$sp
  se <- sqrt(diag(solve(m$H_outer))) # standard errors of log(sp)
  rbind(var = 1 / sp, se = se / sp) # delta method
})

# information criteria
sapply(mods, AIC)
sapply(mods, BIC)

# pseudo-residuals
plot(mods$wild)
plot(mods$modified)

# state-dependent distributions weighted by the local state probabilities
plot(mods$wild, type = "statedens")
plot(mods$modified, type = "statedens")


## plotting state-dependent distributions (paper figure)
if (redraw) pdf("./case_studies/figs/fruitflies_marginal.pdf", width = 8, height = 4)

par(mfrow = c(1,2))
titles <- c("Wild type", "Modified genotype")
x <- 0:150
for (g in 1:2) {
  states <- viterbi(mods[[g]])
  delta_hat <- prop.table(table(factor(states, levels = 1:2)))

  hist(data[[g]]$activity, breaks = 30, bor = "white",
       xlim = c(0, 150), ylim = c(0, 0.03), prob = TRUE,
       main = titles[g], xlab = "Activity count", ylab = "Frequency")
  pmf <- sapply(1:2, function(j) delta_hat[j] * dnbinom2(x, mods[[g]]$mu[j], 1 / mods[[g]]$phi[j]))
  for (j in 1:2) lines(x, pmf[, j], lwd = 2, col = color[j])
  lines(x, rowSums(pmf), lwd = 2, col = "black", lty = 2)
  if (g == 1) legend("topright", legend = c("Inactive", "Active", "Marginal"), col = c(color[1:2], "black"),
                     lty = c(1,1,2), bty = "n", lwd = 2)
}

if (redraw) dev.off()



# Periodically stationary distribution ------------------------------------

L <- 48 # number of unique times of day
conds <- c("LD", "DD")

# periodically stationary distribution at times of day 1:L / 2:
# Gamma[,,s] is the transition into time s / 2, hence row s + 1 of stationary_p()
Delta_tod <- function(mod, aniID, cond) {
  Z_pred <- predict(mod$modmat, data.frame(tod = 1:L / 2, aniID = aniID, condition = cond))
  stationary_p(tpm(mod$beta, Z_pred))[c(2:L, 1), ]
}

# estimated periodically stationary distribution for each fly (time x state x condition x fly)
Delta <- lapply(1:2, function(g)
  sapply(aniIDs[[g]], function(a) sapply(conds, function(cond) Delta_tod(mods[[g]], a, cond),
                                         simplify = "array"), simplify = "array"))

# for an average fly (all random effects zero)
Delta_mean <- lapply(1:2, function(g)
  sapply(conds, function(cond) Delta_tod(mods[[g]], "newID", cond), simplify = "array"))


## plotting
state <- 2

if (redraw) pdf("./case_studies/figs/fruitflies_stationary.pdf", width = 8, height = 3.5)
par(mfrow = c(1, 4))

for (g in 1:2) {
  for (cond in 1:2) {
    plot(NA,
         xlim = c(0, 24), ylim = c(0, 1), xlab = "Time of day", ylab = "Pr(active)",
         bty = "n", main = conds[cond], xaxt = "n"
    )

    # indicate light
    polygon(x = c(0, 8, 8, 0), y = c(-0.005, -0.005, 0.005, 0.005), col = "black", border = "black")
    polygon(x = c(8, 20, 20, 8), y = c(-0.005, -0.005, 0.005, 0.005),
            col = if (cond == 1) "white" else "#00000098", border = "black")
    polygon(x = c(20, 24, 24, 20), y = c(-0.005, -0.005, 0.005, 0.005), col = "black", border = "black")

    axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
    # individual flies
    for (a in 1:nAnimals[g]) {
      lines(1:L / 2, Delta[[g]][, state, cond, a], col = alpha(color[g], 0.3))
    }
    # average fly
    lines(1:L / 2, Delta_mean[[g]][, state, cond], lwd = 2, col = "#000000")
    points(1:L / 2, Delta_mean[[g]][, state, cond], pch = 16, col = "#000000")
  }
}
if (redraw) dev.off()

