## Muskox case study
# Fitting the smooth models takes hours, hence they are saved and loaded by default.
# Set refit = TRUE to fit them again.
refit = FALSE

## packages

# model fitting
library(LaMa) # >= 3.1.0

# store matrix products as single operations on the AD tape (RTMB >= 2.0 expands them by default):
# with the large dense design matrices here, the expanded tape needs more memory than available
TapeConfig(matmul = "atomic")

# spatial stuff
library(sf)
library(sp)
library(rnaturalearth)


## loading the data
data <- readRDS("./data/muskox_winter.rds")


## color palette
color <- LaMaColors(3) # c("orange", "deepskyblue", "seagreen2")


## EDA and spatial prep

# Extract the island boundaries
# Create an sf object from UTM data
utm_sf = st_as_sf(data[,c("x", "y")], coords = c("x", "y"),
                  crs = paste0("+proj=utm +zone=", 27, " +datum=WGS84"),
                  na.fail = FALSE)
latlon_sf = st_transform(utm_sf, crs = 4326)

# create a bounding box
bbox = st_bbox(latlon_sf)

# download map data for greenland only
world <- ne_countries(scale = "large", returnclass = "sf")

# filter for greenland
greenland = world[world$name == "Greenland",]
greenland_geom = st_geometry(greenland)

# restrict to bounding box
greenland_bbox = st_intersection(greenland_geom, st_as_sfc(bbox))

# Extract boundary points from the land polygon
boundary_coords <- st_coordinates(greenland_bbox)[, 1:2]

# Convert to a data frame
boundary_df <- data.frame(x = boundary_coords[,1], y = boundary_coords[,2])

# Ensure it's in UTM coordinates
boundary_df <- st_transform(st_as_sf(boundary_df, coords = c("x", "y"), crs = 4326), crs = "+proj=utm +zone=27 +datum=WGS84")
boundary_df <- st_coordinates(boundary_df)
boundary_df <- data.frame(x = boundary_df[,1], y = boundary_df[,2])

# split into two islands
bnd_main = boundary_df[25:142, ]
bnd_island = boundary_df[143:203, ]



# Homogeneous model -------------------------------------------------------

## negative log-likelihood function
nll <- function(par){
  getAll(par, dat) # make everything available without $
  ## transforming unconstrained parameters to working
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  rho = plogis(logitrho); REPORT(rho)
  ## state process
  Gamma = tpm(beta) # inverse mlogit
  ## state-dependent densities, resting with turning angles around pi
  allprobs = statedens(step, dgamma2, mu, sigma) *
    statedens(angle, dwrpcauchy, c(pi, 0, 0), rho)
  ## forward algorithm, each track starting in the stationary distribution
  -forward("stationary", Gamma, allprobs, ID)
}

## initial parameters
par = list(
  logmu = log(c(5, 65, 400)), # step mean
  logsigma = log(c(5, 60, 400)), # step sd
  logitrho = qlogis(c(0.25, 0.1, 0.5)), # turning angle correlation
  beta = rep(-2, 6)
)

## data
dat = list(
  step = data$step,
  angle = data$angle,
  ID = data$ID
)

## model fitting
obj <- MakeADFun(nll, par)
opt <- nlminb(obj$par, obj$fn, obj$gr)
mod_hom <- lama(obj)

summary(mod_hom)

## state decoding
states <- viterbi(mod_hom)

## state-dependent densities weighted by local state probabilities
plot(mod_hom, type = "statedens")


## visualising fitted state-dependent distributions (paper figure)
mu <- mod_hom$mu
sigma <- mod_hom$sigma
rho <- mod_hom$rho
delta <- stationary(mod_hom$Gamma)

# pdf("./case_studies/figs/muskox_marginal.pdf", width = 8, height = 4)

par(mfrow = c(1,2))
# step length
hist(data$step, breaks = 200, xlim = c(0, 600), ylim = c(0, 0.01), prob = TRUE, bor = "white",
     main = "", xlab = "Step length (m)", ylab = "Density")
for(j in 1:3){
  curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE,
        col = color[j], lwd = 2, n = 500)
}
curve(delta[1] * dgamma2(x, mu[1], sigma[1]) +
        delta[2] * dgamma2(x, mu[2], sigma[2]) +
        delta[3] * dgamma2(x, mu[3], sigma[3]), add = TRUE,
      lty = 2, lwd = 2, n = 500)
legend("topright", legend = c("Resting", "Foraging", "Travelling", "Marginal"), col = c(color[1:3], "black"),
       lty = c(1,1,1,2), bty = "n", lwd = 2)
# turning angle
br <- seq(-pi, pi, length = 30)
hist(data$angle, breaks = br, prob = TRUE, bor = "white",
     main = "", xlab = "Turning angle (radians)", ylab = "Density")
for(j in 1:3){
  curve(delta[j] * dwrpcauchy(x, c(pi,0,0)[j], rho[j]), add = TRUE,
        col = color[j], lwd = 2, n = 500)
}
curve(delta[1] * dwrpcauchy(x, pi, rho[1]) +
        delta[2] * dwrpcauchy(x, 0, rho[2]) +
        delta[3] * dwrpcauchy(x, 0, rho[3]), add = TRUE,
      lty = 2, lwd = 2, n = 500)

# dev.off()



# Smooth models -----------------------------------------------------------

# Only the transition probability travelling -> foraging (Gamma[3,2,], the 4th
# off-diagonal entry) depends on covariates, the others are constant.

## penalised negative log-likelihood function for all smooth models
# the spline coefficients are those named in dat$random, one vector per smooth
pnll <- function(par){
  getAll(par, dat, warn = FALSE)

  ## transforming parameters
  mu <- exp(logmu); REPORT(mu)
  sigma <- exp(logsigma); REPORT(sigma)
  rho <- plogis(logitrho); REPORT(rho)

  ## state process
  beta_spline <- mget(random) # list of spline coefficient vectors
  beta <- AD(matrix(0, 6, ncol(Z)))
  beta[-4, 1] <- eta # intercepts of the constant entries
  beta[4, ] <- do.call(c, c(list(beta0), beta_spline)) # travelling -> foraging
  Gamma <- tpm(beta, Z)

  ## state-dependent densities
  allprobs <- statedens(step, dgamma2, mu, sigma) * statedens(angle, dwrpcauchy, c(pi, 0, 0), rho)

  ## forward algorithm, each track starting in the stationary distribution
  # of the transition probabilities at its first observation
  -forward("stationary", Gamma, allprobs, ID) +
    penalty2(beta_spline, S, sp)
}

## fitting a smooth model given the formula for travelling -> foraging
# names of the coefficient vectors for each smooth term
coef_names <- c("s(x,y)" = "beta_xy", "s(tday)" = "beta_t", "ti(x,y,tday)" = "beta_xy_t")

fit_smooth <- function(formula, sp = NULL,
                       logmu = log(c(5, 66, 392)),
                       logsigma = log(c(4, 60, 380)),
                       logitrho = qlogis(c(0.23, 0.11, 0.53)),
                       silent = 0, start = NULL) { # start: optional warm start for par
  m <- modmat(formula, data = data, knots = list(tday = c(0, 24)))
  random <- unname(coef_names[names(m$S)])
  if (is.null(sp)) sp <- m$sp0 # initial penalty strengths from mgcv

  par <- c(list(logmu = logmu, logsigma = logsigma, logitrho = logitrho,
                eta = rep(-2, 5),  # intercepts of the constant entries
                beta0 = -2),       # intercept of travelling -> foraging
           setNames(lapply(m$pardim[names(m$S)], numeric), random)) # spline coefficients
  if (!is.null(start)) par <- start

  dat <- list(step = data$step, angle = data$angle, ID = data$ID,
              Z = m$Z, S = m$S, random = random,
              sp = sp) # one initial penalty strength per penalty

  # fitting time excludes building the model matrices
  time <- system.time(mod <- qreml(pnll, par, dat, random = random, silent = silent))
  mod$time <- time
  mod$modmat <- m # for predictions later
  mod
}

if (refit) {
  # initial penalty strengths: the heuristic sp0 from mgcv overshoots by several orders of
  # magnitude here (e.g. 32400 for s(x,y) and the cap 1e6 for s(tday), estimates 0.4 and 73),
  # so the marginal smooths start conservatively at about ten times the estimates of a previous
  # additive fit, and the interaction strongly penalised
  ## spatial smooth only
  mod_space <- fit_smooth(~ s(x, y, k = 70), sp = 5)
  print(mod_space$time)
  saveRDS(mod_space, file = "./case_studies/objects/muskox_mod_space_v4.rds")
  rm(mod_space); gc() # RTMB memory piles up across fits otherwise

  ## time of day only
  mod_time <- fit_smooth(~ s(tday, bs = "cc", k = 10), sp = 500)
  print(mod_time$time)
  # joint precision of the coefficients and log(sp), for intervals that also account for the
  # uncertainty in sp, evaluated now as outer_hessian() does not work any more after loading
  mod_time$H_outer <- mod_time$outer_hessian()
  mod_time$Q_joint <- mod_time$joint_uncertainty(H_rho = mod_time$H_outer)
  saveRDS(mod_time, file = "./case_studies/objects/muskox_mod_time_v4.rds")
  rm(mod_time); gc()

  ## additive space and time of day
  mod_add <- fit_smooth(~ s(x, y, bs = "tp", k = 70) + s(tday, bs = "cc", k = 10),
                        sp = c(5, 500))
  print(mod_add$time)
  saveRDS(mod_add, file = "./case_studies/objects/muskox_mod_add_v4.rds")
  rm(mod_add); gc()

  ## space-time interaction
  mod_spacetime <- fit_smooth(~ s(x, y, bs = "tp", k = 70) + s(tday, bs = "cc", k = 10) +
                                ti(x, y, tday, d = c(2,1), bs = c("tp", "cc"), k = c(30, 8)),
                              sp = c(5, 500, 1e3, 1e3))
  print(mod_spacetime$time)
  saveRDS(mod_spacetime, file = "./case_studies/objects/muskox_mod_spacetime_v4.rds")
  rm(mod_spacetime); gc()
}

mod_space <- readRDS("./case_studies/objects/muskox_mod_space_v4.rds")
mod_time <- readRDS("./case_studies/objects/muskox_mod_time_v4.rds")
mod_add <- readRDS("./case_studies/objects/muskox_mod_add_v4.rds")
mod_spacetime <- readRDS("./case_studies/objects/muskox_mod_spacetime_v4.rds")

summary(mod_space)
summary(mod_spacetime)

## convergence check: space-time model from below, i.e. with a weak penalty on the interaction,
# warm started at the final estimates. Takes hours and is not run by default
if (FALSE) {
  rm(mod_space, mod_time, mod_add); gc() # free memory for the large design matrices
  mod_below <- fit_smooth(~ s(x, y, bs = "tp", k = 70) + s(tday, bs = "cc", k = 10) +
                            ti(x, y, tday, d = c(2,1), bs = c("tp", "cc"), k = c(30, 8)),
                          sp = c(0.15, 70, 1, 0.1), start = mod_spacetime$par)
  print(mod_below$time)
  saveRDS(mod_below, file = "./case_studies/objects/muskox_mod_spacetime_v4_below.rds")
  rbind(above = mod_spacetime$sp, below = mod_below$sp)
  rbind(above = unlist(mod_spacetime$edf), below = unlist(mod_below$edf))
  c(above = mod_spacetime$llk_restricted[mod_spacetime$best_iter],
    below = mod_below$llk_restricted[mod_below$best_iter])
}


# Model comparison --------------------------------------------------------

mods <- list(hom = mod_hom, space = mod_space, time = mod_time,
             add = mod_add, spacetime = mod_spacetime)

IC_table <- data.frame(llk = sapply(mods, logLik),
                       AIC = sapply(mods, AIC),
                       BIC = sapply(mods, BIC))
round(IC_table, 2)

# differences to the best model
diff_IC_table <- IC_table
diff_IC_table$AIC <- IC_table$AIC - min(IC_table$AIC)
diff_IC_table$BIC <- IC_table$BIC - min(IC_table$BIC)
round(diff_IC_table, 2)



# Estimated fields --------------------------------------------------------

## prediction grid, restricted to land within the convex hull of the data
x_seq <- seq(min(data$x), max(data$x), length.out = 150)
y_seq <- seq(min(data$y), max(data$y), length.out = 150)
grid_xy <- expand.grid(x = x_seq, y = y_seq)

hull_idx <- chull(data$x, data$y)  # convex hull indices
hull_coords <- cbind(data$x[hull_idx], data$y[hull_idx])  # hull boundary points

inside_chull <- point.in.polygon(grid_xy$x, grid_xy$y, hull_coords[,1], hull_coords[,2])
inside_main <- point.in.polygon(grid_xy$x, grid_xy$y, bnd_main$x, bnd_main$y)
inside_island <- point.in.polygon(grid_xy$x, grid_xy$y, bnd_island$x, bnd_island$y)
inside <- (inside_main | inside_island) & inside_chull

## Pr(travelling -> foraging) on the grid, NA outside, as a matrix for image()
pr_field <- function(mod, newdata, beta = mod$beta) {
  pr <- tpm(beta, predict(mod$modmat, newdata))[3,2,]
  pr[!inside] <- NA
  matrix(pr, length(x_seq), length(y_seq))
}

## image of a field with a colour gradient legend above
plot_field <- function(Pr, main, zlim = c(0, 1)) {
  image(x_seq, y_seq, Pr, col = hcl.colors(100), zlim = zlim,
        xlim = c(465000, 560000), ylim = c(8230000, 8300000),
        xlab = "UTM easting", ylab = "UTM northing",
        main = main, bty = "n", asp = 1, useRaster = FALSE)

  legend_x <- c(480000, 545000)  # x range for legend
  legend_y <- c(8300000-4000, 8300000-2000)  # above the plot
  rasterImage(as.raster(matrix(hcl.colors(100), nrow = 1)),
              legend_x[1], legend_y[1], legend_x[2], legend_y[2])
  by <- diff(zlim) / 5
  ticks <- round(seq(ceiling(zlim[1] / by) * by, floor(zlim[2] / by) * by, by = by), 2)
  text(legend_x[1] + (ticks - zlim[1]) / diff(zlim) * diff(legend_x), legend_y[2] + 2500,
       labels = ticks, cex = 0.9)
}


## spatial model
# pdf("./case_studies/figs/muskox_space.pdf", width = 7, height = 7)
# png("./case_studies/figs/muskox_space.png", width = 3000, height = 3000, res = 450)
par(mfrow = c(1,1))
plot_field(pr_field(mod_space, grid_xy), main = expression("Pr(travelling" %->% "foraging)"))
# dev.off()


## space-time model
Pr_t <- lapply(1:24, function(t) pr_field(mod_spacetime, cbind(grid_xy, tday = t)))

# animation over the day
par(mfrow = c(1,1))
for(t in 1:24){
  plot_field(Pr_t[[t]], main = paste0("Pr(travelling → foraging) - time of day: ", t, ":00"))
  Sys.sleep(0.25)
}

## contribution of the interaction to Pr(travelling -> foraging): the fit vs. the same fit without ti()
n_ti <- mod_spacetime$modmat$pardim[["ti(x,y,tday)"]]
ti_ind <- ncol(mod_spacetime$beta) - n_ti + seq_len(n_ti) # last coefficients of beta[4, ]
beta_noti <- mod_spacetime$beta
beta_noti[4, ti_ind] <- 0
Pr_ti <- lapply(1:24, function(t) {
  Pr_t[[t]] - pr_field(mod_spacetime, cbind(grid_xy, tday = t), beta = beta_noti)
})

## the two hours most different in that contribution: largest mean absolute difference over the study region
d_t <- outer(1:24, 1:24, Vectorize(function(s, t) mean(abs(Pr_ti[[s]] - Pr_ti[[t]]), na.rm = TRUE)))
hours <- sort(arrayInd(which.max(d_t), dim(d_t))[1, ])
print(hours); print(max(d_t))
Pr_diff <- Pr_ti[[hours[2]]] - Pr_ti[[hours[1]]]
range(Pr_diff, na.rm = TRUE)


## time of day model: Pr(travelling -> foraging) over the day with pointwise intervals
tday_seq <- seq(0, 24, length.out = 200)
Z_t <- predict(mod_time$modmat, data.frame(tday = tday_seq))
gamma_t <- tpm(mod_time$beta, Z_t)[3,2,]

set.seed(123)
samples <- MCreport(mod_time, nSamples = 5000, Q = mod_time$Q_joint)
gamma_t_samples <- sapply(seq_along(samples$beta0), function(s) {
  beta <- mod_time$beta
  beta[-4, 1] <- samples$eta[[s]]
  beta[4, ] <- c(samples$beta0[[s]], samples$beta_t[[s]])
  tpm(beta, Z_t)[3,2,]
})
gamma_t_q <- apply(gamma_t_samples, 1, quantile, probs = c(0.025, 0.975))


## paper figure: time of day model, spatial model and interaction difference. Map and difference
# use the same palette over colour ranges of the same width, so that the same difference in
# colour means the same difference in probability
Pr_space <- pr_field(mod_space, grid_xy)
zlim_width <- 0.5
zlim_map <- c(0, zlim_width)
zlim_diff <- c(-1, 1) * zlim_width / 2

png("./case_studies/figs/muskox_models.png", width = 4600, height = 1800, res = 600)
layout(matrix(1:3, nrow = 1), widths = c(0.75, 1, 1))
par(mar = c(5, 4, 4, 2) + 0.1)
plot(NA, xlim = c(0, 24), ylim = c(0.05, 0.3), bty = "n", xaxt = "n",
     xlab = "time of day", ylab = expression("Pr(travelling" %->% "foraging)"),
     main = "Time of day model")
axis(1, at = seq(0, 24, by = 6))
polygon(c(tday_seq, rev(tday_seq)), c(gamma_t_q[1, ], rev(gamma_t_q[2, ])),
        col = adjustcolor(color[2], alpha.f = 0.3), border = NA)
lines(tday_seq, gamma_t, col = color[2], lwd = 2)
plot_field(Pr_space, main = "Spatial model", zlim = zlim_map)
plot_field(Pr_diff, main = paste0("Interaction: ", hours[2], ":00 - ", hours[1], ":00"),
           zlim = zlim_diff)
dev.off()
