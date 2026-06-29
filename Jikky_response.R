### Jikky response

### Traversa fit
# ============================================================

# Traversa SIDS age curve: LOWESS + candidate smooth baselines

# ============================================================


# Morris/Traversa digitised anchors

anchor_day <- c(
  
  6,13,20,27,34,41,48,55,62,69,76,83,90,97,104,111,118,125,132,139,
  
  146,153,160,167,174,181,188,195,205,220,240,260,280,300,330,360,
  
  400,440,480,520,560,600,640,680,715
  
)


anchor_cnt <- c(
  
  12,16,20,24,31,23,30,28,33,28,37,24,21,24,23,16,15,14,12,18,
  
  21,13,11,9,8,7,6,5,4,3,2.6,2.2,1.9,1.7,1.5,1.4,1.2,1.1,1.0,0.9,
  
  0.9,0.8,0.9,0.9,0.9
  
)


dat <- data.frame(day = anchor_day, count = anchor_cnt)


# LOWESS reference, as in Morris

lo <- lowess(dat$day, dat$count, f = 0.2)

days <- 1:720

lowess_fit <- pmax(approx(lo$x, lo$y, xout = days, rule = 2)$y, 0)


# ------------------------------------------------------------

# Fit helper

# ------------------------------------------------------------


safe_nls <- function(formula, data, start, lower = NULL, upper = NULL) {
  
  tryCatch(
    
    nls(
      
      formula,
      
      data = data,
      
      start = start,
      
      algorithm = "port",
      
      lower = lower,
      
      upper = upper,
      
      control = nls.control(maxiter = 1000, warnOnly = TRUE)
      
    ),
    
    error = function(e) NULL
    
  )
  
}


fit_stats <- function(obs, pred, k) {
  
  rss <- sum((obs - pred)^2)
  
  n <- length(obs)
  
  rmse <- sqrt(mean((obs - pred)^2))
  
  aic <- n * log(rss / n) + 2 * k
  
  bic <- n * log(rss / n) + log(n) * k
  
  c(RSS = rss, RMSE = rmse, AIC = aic, BIC = bic)
  
}


fits <- list()


# ------------------------------------------------------------

# 1. Linear decline

# ------------------------------------------------------------


fit_linear <- lm(count ~ day, data = dat)

pred_linear <- pmax(predict(fit_linear, newdata = data.frame(day = days)), 0)

fits$Linear <- list(
  
  pred = pred_linear,
  
  pred_anchor = pmax(predict(fit_linear, newdata = dat), 0),
  
  k = 2
  
)


# ------------------------------------------------------------

# 2. Exponential decline: y = c + A exp(-k day)

# ------------------------------------------------------------


fit_exp <- safe_nls(
  
  count ~ c0 + A * exp(-k * day),
  
  data = dat,
  
  start = list(c0 = 0.8, A = 35, k = 0.01),
  
  lower = c(c0 = 0, A = 0, k = 0),
  
  upper = c(c0 = 10, A = 100, k = 1)
  
)


if (!is.null(fit_exp)) {
  
  pred <- predict(fit_exp, newdata = data.frame(day = days))
  
  fits$Exponential <- list(
    
    pred = pmax(pred, 0),
    
    pred_anchor = pmax(predict(fit_exp, newdata = dat), 0),
    
    k = 3
    
  )
  
}


# ------------------------------------------------------------

# 3. Weibull-shaped curve with offset:

# y = c + A * dweibull(day, shape, scale)

# ------------------------------------------------------------


fit_weibull <- safe_nls(
  
  count ~ c0 + A * dweibull(day, shape = shape, scale = scale),
  
  data = dat,
  
  start = list(c0 = 0.8, A = 3000, shape = 1.8, scale = 90),
  
  lower = c(c0 = 0, A = 0, shape = 0.2, scale = 10),
  
  upper = c(c0 = 10, A = 100000, shape = 10, scale = 1000)
  
)


if (!is.null(fit_weibull)) {
  
  pred <- predict(fit_weibull, newdata = data.frame(day = days))
  
  fits$Weibull <- list(
    
    pred = pmax(pred, 0),
    
    pred_anchor = pmax(predict(fit_weibull, newdata = dat), 0),
    
    k = 4
    
  )
  
}


# ------------------------------------------------------------

# 4. Gamma-shaped curve with offset:

# y = c + A * dgamma(day, shape, scale)

# ------------------------------------------------------------


fit_gamma <- safe_nls(
  
  count ~ c0 + A * dgamma(day, shape = shape, scale = scale),
  
  data = dat,
  
  start = list(c0 = 0.8, A = 3000, shape = 2.0, scale = 45),
  
  lower = c(c0 = 0, A = 0, shape = 0.2, scale = 1),
  
  upper = c(c0 = 10, A = 100000, shape = 20, scale = 500)
  
)


if (!is.null(fit_gamma)) {
  
  pred <- predict(fit_gamma, newdata = data.frame(day = days))
  
  fits$Gamma <- list(
    
    pred = pmax(pred, 0),
    
    pred_anchor = pmax(predict(fit_gamma, newdata = dat), 0),
    
    k = 4
    
  )
  
}


# ------------------------------------------------------------

# 5. Log-normal-shaped curve with offset:

# y = c + A * dlnorm(day, meanlog, sdlog)

# ------------------------------------------------------------


fit_lognormal <- safe_nls(
  
  count ~ c0 + A * dlnorm(day, meanlog = meanlog, sdlog = sdlog),
  
  data = dat,
  
  start = list(c0 = 0.8, A = 3000, meanlog = log(80), sdlog = 0.7),
  
  lower = c(c0 = 0, A = 0, meanlog = log(5), sdlog = 0.1),
  
  upper = c(c0 = 10, A = 100000, meanlog = log(500), sdlog = 3)
  
)


if (!is.null(fit_lognormal)) {
  
  pred <- predict(fit_lognormal, newdata = data.frame(day = days))
  
  fits$Lognormal <- list(
    
    pred = pmax(pred, 0),
    
    pred_anchor = pmax(predict(fit_lognormal, newdata = dat), 0),
    
    k = 4
    
  )
  
}


# ------------------------------------------------------------

# 6. Natural cubic spline, deliberately modest df

# ------------------------------------------------------------


library(splines)


fit_spline3 <- lm(count ~ ns(day, df = 3), data = dat)

pred_spline3 <- pmax(predict(fit_spline3, newdata = data.frame(day = days)), 0)


fits$Spline_df3 <- list(
  
  pred = pred_spline3,
  
  pred_anchor = pmax(predict(fit_spline3, newdata = dat), 0),
  
  k = 4
  
)


fit_spline5 <- lm(count ~ ns(day, df = 5), data = dat)

pred_spline5 <- pmax(predict(fit_spline5, newdata = data.frame(day = days)), 0)


fits$Spline_df5 <- list(
  
  pred = pred_spline5,
  
  pred_anchor = pmax(predict(fit_spline5, newdata = dat), 0),
  
  k = 6
  
)


# ------------------------------------------------------------

# Fit table

# ------------------------------------------------------------


stats <- do.call(
  
  rbind,
  
  lapply(names(fits), function(nm) {
    
    c(Model = nm, fit_stats(dat$count, fits[[nm]]$pred_anchor, fits[[nm]]$k))
    
  })
  
)


stats <- as.data.frame(stats)

stats[, -1] <- lapply(stats[, -1], as.numeric)

stats <- stats[order(stats$AIC), ]


print(stats, row.names = FALSE)


# ------------------------------------------------------------

# Plot overlay

# ------------------------------------------------------------


plot(
  
  dat$day, dat$count,
  
  pch = 19,
  
  xlim = c(0, 720),
  
  ylim = c(0, max(dat$count) * 1.15),
  
  xlab = "Age at death (days)",
  
  ylab = "Digitised Traversa count units",
  
  main = "Traversa SIDS age distribution: smooth candidate baselines"
  
)


lines(days, lowess_fit, lwd = 4, col = "black")


cols <- c(
  
  Linear = "grey50",
  
  Exponential = "blue",
  
  Weibull = "red",
  
  Gamma = "darkgreen",
  
  Lognormal = "purple",
  
  Spline_df3 = "orange",
  
  Spline_df5 = "brown"
  
)


for (nm in names(fits)) {
  
  lines(days, fits[[nm]]$pred, lwd = 2, col = cols[nm])
  
}


legend(
  
  "topright",
  
  legend = c("Traversa anchors", "Morris LOWESS", names(fits)),
  
  col = c("black", "black", cols[names(fits)]),
  
  pch = c(19, NA, rep(NA, length(fits))),
  
  lty = c(NA, 1, rep(1, length(fits))),
  
  lwd = c(NA, 4, rep(2, length(fits))),
  
  bty = "n",
  
  cex = 0.8
  
)



readline("Press enter for next plot:")


# ------------------------------------------------------------

# Residual plot

# ------------------------------------------------------------


plot(
  
  dat$day, dat$count - approx(days, lowess_fit, xout = dat$day)$y,
  
  pch = 19,
  
  type = "b",
  
  xlim = c(0, 720),
  
  ylim = range(unlist(lapply(fits, function(x) dat$count - x$pred_anchor))),
  
  xlab = "Age at death (days)",
  
  ylab = "Observed - fitted",
  
  main = "Residuals against candidate baselines"
  
)


abline(h = 0, lty = 2, col = "grey50")


for (nm in names(fits)) {
  
  lines(
    
    dat$day,
    
    dat$count - fits[[nm]]$pred_anchor,
    
    type = "b",
    
    pch = 19,
    
    col = cols[nm]
    
  )
  
}


legend(
  
  "topright",
  
  legend = names(fits),
  
  col = cols[names(fits)],
  
  lty = 1,
  
  pch = 19,
  
  bty = "n",
  
  cex = 0.8
  
)
### Constrained Gamma fit
# ============================================================

# Constrained Gamma fit to Morris/Traversa LOWESS

# Forces fit near anchor minima at day ~41 and ~125

# ============================================================


sids_distribution <- function(maxday = 365L, span = 0.2) {
  
  anchor_day <- c(6,13,20,27,34,41,48,55,62,69,76,83,90,97,104,111,118,125,132,139,
                  
                  146,153,160,167,174,181,188,195,205,220,240,260,280,300,330,360,
                  
                  400,440,480,520,560,600,640,680,715)
  
  anchor_cnt <- c(12,16,20,24,31,23,30,28,33,28,37,24,21,24,23,16,15,14,12,18,
                  
                  21,13,11,9,8,7,6,5,4,3,2.6,2.2,1.9,1.7,1.5,1.4,1.2,1.1,1.0,0.9,
                  
                  0.9,0.8,0.9,0.9,0.9)
  
  lo <- lowess(anchor_day, anchor_cnt, f = span)
  
  fit_cnt <- pmax(approx(lo$x, lo$y, xout = 1:maxday, rule = 2)$y, 0)
  
  list(days = 1:maxday, anchor_day = anchor_day,
       
       anchor_cnt = anchor_cnt, fit_cnt = fit_cnt)
  
}


s <- sids_distribution(365, 0.2)


target <- data.frame(
  
  day = s$days,
  
  count = s$fit_cnt
  
)


anchors <- data.frame(
  
  day = s$anchor_day,
  
  count = s$anchor_cnt
  
)


fit_start <- 35

fit_dat <- subset(target, day >= fit_start)


# ------------------------------------------------------------

# Gamma prediction function

# Parameterisation:

# par = c(c0, A, shape, scale)

# ------------------------------------------------------------


gamma_pred_fun <- function(par, day) {
  
  c0 <- par[1]
  
  A <- par[2]
  
  shape <- par[3]
  
  scale <- par[4]
  
  c0 + A * dgamma(day, shape = shape, scale = scale)
  
}


# ------------------------------------------------------------

# Constraint targets

# Use LOWESS target values at day 41 and 125

# ------------------------------------------------------------


anchor_days <- c(41, 125)

anchor_targets <- target$count[match(anchor_days, target$day)]


# ------------------------------------------------------------

# Optional region weights:

# Down-weight suspected excess regions so they do not pull

# the smooth Gamma upward.

# ------------------------------------------------------------


fit_dat$w <- 1

fit_dat$w[fit_dat$day >= 55 & fit_dat$day <= 95] <- 0.25

fit_dat$w[fit_dat$day >= 140 & fit_dat$day <= 175] <- 0.25


# ------------------------------------------------------------

# Penalised objective

# Increase penalty if it does not pass close enough through

# day 41 and day 125.

# ------------------------------------------------------------


penalty <- 1e5


loss_fun <- function(par) {
  
  c0 <- par[1]
  
  A <- par[2]
  
  shape <- par[3]
  
  scale <- par[4]
  
  # enforce sensible parameter bounds manually
  
  if (c0 < 0 || A <= 0 || shape <= 0 || scale <= 0) return(Inf)
  
  if (c0 > 10 || A > 100000 || shape > 20 || scale > 500) return(Inf)
  
  pred <- gamma_pred_fun(par, fit_dat$day)
  
  rss <- sum(fit_dat$w * (fit_dat$count - pred)^2)
  
  pred_anchor <- gamma_pred_fun(par, anchor_days)
  
  constraint_loss <- sum((anchor_targets - pred_anchor)^2)
  
  rss + penalty * constraint_loss
  
}


# ------------------------------------------------------------

# Optimise

# ------------------------------------------------------------


start <- c(
  
  c0 = 1.17,
  
  A = 3373.6,
  
  shape = 2.55,
  
  scale = 37.34
  
)


fit_opt <- optim(
  
  par = start,
  
  fn = loss_fun,
  
  method = "Nelder-Mead",
  
  control = list(maxit = 50000, reltol = 1e-12)
  
)


fit_opt$par

fit_opt$value

fit_opt$convergence


gamma_pred <- pmax(gamma_pred_fun(fit_opt$par, target$day), 0)

resid_gamma <- target$count - gamma_pred


cat("\nAnchor checks:\n")

print(data.frame(
  
  day = anchor_days,
  
  LOWESS = anchor_targets,
  
  Gamma = gamma_pred_fun(fit_opt$par, anchor_days),
  
  residual = anchor_targets - gamma_pred_fun(fit_opt$par, anchor_days)
  
))


cat("\nPositive residual mass / total LOWESS mass:\n")

cat(sum(pmax(resid_gamma, 0)) / sum(target$count), "\n")

gamma_par <- fit_opt$par

print(gamma_par)


# ------------------------------------------------------------

# Plot

# ------------------------------------------------------------


op <- par(mfrow = c(2, 1), mar = c(4.5, 4.5, 3, 1), oma = c(0, 0, 2, 0))


plot(
  
  target$day, target$count,
  
  type = "l",
  
  lwd = 3,
  
  col = "black",
  
  xlim = c(0, 365),
  
  ylim = c(0, max(target$count) * 1.15),
  
  xlab = "Age at death (days)",
  
  ylab = "Morris/Traversa LOWESS count units",
  
  main = "LOWESS target vs constrained Gamma"
  
)


points(
  
  anchors$day[anchors$day <= 365],
  
  anchors$count[anchors$day <= 365],
  
  pch = 16,
  
  cex = 0.55,
  
  col = adjustcolor("black", 0.35)
  
)


lines(target$day, gamma_pred, lwd = 3, col = "blue")


points(anchor_days, anchor_targets, pch = 19, col = "red", cex = 1.3)


abline(v = fit_start, lty = 3, col = "grey50")


legend(
  
  "topright",
  
  legend = c("Morris/Traversa LOWESS", "Original anchors", "Constrained Gamma", "Constraint points"),
  
  col = c("black", adjustcolor("black", 0.35), "blue", "red"),
  
  pch = c(NA, 16, NA, 19),
  
  lty = c(1, NA, 1, NA),
  
  lwd = c(3, NA, 3, NA),
  
  bty = "n",
  
  cex = 0.8
  
)


plot(
  
  target$day, resid_gamma,
  
  type = "l",
  
  lwd = 3,
  
  col = "blue",
  
  xlim = c(35, 365),
  
  ylim = range(resid_gamma[target$day >= 35]) * 1.15,
  
  xlab = "Age at death (days)",
  
  ylab = "LOWESS target - constrained Gamma",
  
  main = "Residual excess over Gamma baseline"
  
)


abline(h = 0, lty = 2, col = "grey50")


polygon(
  
  c(target$day, rev(target$day)),
  
  c(pmax(resid_gamma, 0), rep(0, length(resid_gamma))),
  
  col = adjustcolor("blue", 0.2),
  
  border = NA
  
)


lines(target$day, resid_gamma, lwd = 3, col = "blue")

points(anchor_days, rep(0, length(anchor_days)), pch = 19, col = "red", cex = 1.3)


mtext("Morris/Traversa LOWESS with constrained Gamma natural baseline", outer = TRUE, font = 2)


par(op)



###### Define parameters for derived gamma model ######

#-----------------------------------------------------#

our_gamma_pmf <- function(maxday = 365,
                          
                          c0 = 1.709145,
                          
                          A = 3036.988093,
                          
                          shape = 2.608024,
                          
                          scale = 32.905001) {
  
  
  days <- 1:maxday
  
  
  curve <- pmax(
    
    c0 + A * dgamma(days,
                    
                    shape = shape,
                    
                    scale = scale),
    
    0
    
  )
  
  
  curve / sum(curve)
  
}


###----- For the SIDS model -----###

gamma_pmf <- our_gamma_pmf(365)


### MCSCCS positive control gamma bootstrap function
library(SCCS)


run_one_msccs <- function(N = 5e6,
                          
                          OBS_END = 150L,
                          
                          p_death = 2e-4 * 150 / 365,
                          
                          RR = 1.5,
                          
                          RISK_LEN = 3L,
                          
                          agegrp=NULL,
                          
                          dose1_day = 60L,
                          
                          dose2_day = 120L) {
  
  days <- 1:OBS_END
  
  risk_window <- rep(FALSE, OBS_END)
  
  risk_window[dose1_day:(dose1_day + RISK_LEN - 1)] <- TRUE
  
  risk_window[dose2_day:(dose2_day + RISK_LEN - 1)] <- TRUE
  
  
  ####--- Use derived gamma model as baseline incidence ---####
  
  
  gamma_pmf <- our_gamma_pmf(365)
  
  haz <- gamma_pmf[1:OBS_END]
  
  haz[risk_window] <- haz[risk_window] * RR
  
  ####-----------------------------------------------------####
  
  
  
  pmf_death <- haz / sum(haz)
  
  
  
  n_deaths <- rpois(1, N * p_death)
  
  
  death_day <- sample(days, n_deaths, replace = TRUE, prob = pmf_death)
  
  
  cases <- data.frame(
    
    id = seq_len(n_deaths),
    
    death = death_day
    
  )
  
  
  cases$d1 <- ifelse(dose1_day <= cases$death, dose1_day, NA_integer_)
  
  cases$d2 <- ifelse(dose2_day <= cases$death, dose2_day, NA_integer_)
  
  md <- data.frame(
    
    indiv = cases$id,
    
    astart = 1L,
    
    aend = OBS_END,
    
    aevent = cases$death,
    
    ad1 = cases$d1,
    
    ad2 = cases$d2
    
  )
  
  fit <- tryCatch(
    
    eventdepenexp(
      
      indiv = indiv,
      
      astart = astart,
      
      aend = aend,
      
      aevent = aevent,
      
      adrug = cbind(ad1, ad2),
      
      aedrug = cbind(ad1, ad2) + RISK_LEN,
      
      expogrp = 0,
      
      sameexpopar = TRUE,
      
      agegrp = agegrp,
      
      dataformat = "multi",
      
      data = md
      
    ),
    
    error = function(e) NULL
    
  )
  
  risk_days <- which(risk_window)
  
  
  obs_risk <- sum(cases$death %in% risk_days)
  
  
  exp_frac_gamma <- sum(gamma_pmf[risk_days]) / sum(gamma_pmf[1:OBS_END])
  
  exp_risk_gamma <- exp_frac_gamma * n_deaths
  
  
  oe <- obs_risk / exp_risk_gamma
  
  
  
  if (is.null(fit)) {
    
    return(data.frame(
      
      n_deaths = n_deaths,
      
      obs_risk = obs_risk,
      
      exp_risk = exp_risk_gamma, ### fixed error in exp_risk_flat
      
      OE = oe,
      
      IRR = NA,
      
      lower = NA,
      
      upper = NA,
      
      p = NA
      
    ))
    
  }
  
  irr <- fit$conf.int[1, 1]
  
  lo <- fit$conf.int[1, 3]
  
  hi <- fit$conf.int[1, 4]
  
  se <- (log(hi) - log(lo)) / (2 * qnorm(0.975))
  
  p <- 2 * pnorm(-abs(log(irr) / se))
  
  data.frame(
    
    n_deaths = n_deaths,
    
    obs_risk = obs_risk,
    
    exp_risk = exp_risk_gamma,
    
    OE = oe,
    
    IRR = irr,
    
    lower = lo,
    
    upper = hi,
    
    p = p
    
  )
  
}


run_msccs_sims <- function(K = 200, seed = 20260627, agegrp = NULL, ...) {
  
  
  if (!is.null(seed)) set.seed(seed)
  
  
  agegrp_label <- if (is.null(agegrp)) {
    
    "none"
    
  } else {
    
    paste(agegrp, collapse = ",")
    
  }
  
  
  out <- vector("list", K)
  
  
  for (i in seq_len(K)) {
    
    if (i %% 10 == 0) cat("Simulation", i, "of", K, "\n")
    
    out[[i]] <- run_one_msccs(agegrp = agegrp, ...)
    
    out[[i]]$sim <- i
    
  }
  
  
  results <- do.call(rbind, out)
  
  results <- results[, c("sim", setdiff(names(results), "sim"))]
  
  
  summary <- data.frame(
    
    K = K,
    
    n_age_bins = ifelse(is.null(agegrp), 0, length(agegrp) + 1),
    
    mean_deaths = mean(results$n_deaths),
    
    mean_risk_deaths = mean(results$obs_risk),
    
    mean_OE = mean(results$OE, na.rm = TRUE),
    
    median_OE = median(results$OE, na.rm = TRUE),
    
    mean_IRR = mean(results$IRR, na.rm = TRUE),
    
    median_IRR = median(results$IRR, na.rm = TRUE),
    
    power_p_lt_0.05 = mean(results$p < 0.05, na.rm = TRUE),
    
    power_CI_excludes_1 = mean(results$lower > 1 | results$upper < 1, na.rm = TRUE),
    
    agegrp = agegrp_label
    
  )
  
  
  list(results = results, summary = summary)
  
}






### Simulate age bins
table<-NULL

bins <- list(
  
  none = NULL,
  
  monthly = c(30,60,90,120),
  
  fortnightly = seq(14,140,14),
  
  weekly = seq(7,147,7),
  
  morris = c(60,75,90,105,120,145,170,200,250,300)
  
)


for (i in 1:length(bins)) {
  
  
  sim <- run_msccs_sims(
    
    K = 100,
    
    seed=NULL,
    
    N = 5e6,
    
    OBS_END = 365L,
    
    p_death = 2e-4 * OBS_END /365,
    
    RR = 1.5,
    
    RISK_LEN = 3L,
    
    agegrp = bins[[i]],
    
    dose1_day = 60L,
    
    dose2_day = 120L
    
  )
  
  
  print(sim$summary)
  
  table<-rbind(table,sim$summary)
  
}


print(table)





# Optional plots

hist(sim$results$IRR, breaks = 30, main = "Estimated mSCCS IRR", xlab = "IRR")

abline(v = 1.5, col = "red", lwd = 2)


readline("Press enter for next plot:")


plot(sim$results$OE, sim$results$IRR,
     
     xlab = "Observed / expected",
     
     ylab = "mSCCS IRR",
     
     main = "Crude O/E vs mSCCS IRR",
     
     pch = 19)

abline(0, 1, col = "red", lwd = 2)



#### Rerun using jittered doses like real life.
## ===========================================================================
## SECTION 3: Re-run Jikky's simulation with REALISTIC UK dose-timing
## variability instead of point-mass dose days.
##
## His Section 2 fixes dose1_day = 60 and dose2_day = 120 for EVERY infant, so
## "age" and "days since dose" are the same variable and no SCCS can separate
## them. This section draws each infant's dose-1 age from the published UK
## distribution and spaces dose 2 a realistic interval later, then re-runs his
## exact bin comparison. Everything else (gamma baseline, RR, RISK_LEN, fit,
## summaries) is his, unchanged.
## ===========================================================================

## --- 3.1  UK dose-1 age distribution (per-day pmf over 1..365) --------------
## Two-part model fitted to the MCS Welsh cohort (Cottrell 2022; Hungerford
## 2016): ever-vaccinated ~ Bernoulli(0.986); age|vacc = 6 + Lognormal wk.
## Reproduces the UK bands <8wk 14.2% | 8-12wk 79.6% | >12wk 4.8%.
uk_dose1_pmf <- function(maxday = 365L,
                         p_ever = 0.986, shift = 6,
                         meanlog = 1.1222, sdlog = 0.4039) {
  d   <- 1:maxday
  Fd  <- function(x){ wk <- x/7; ifelse(wk > shift, p_ever*plnorm(wk-shift, meanlog, sdlog), 0) }
  cum <- Fd(d)
  pmf <- c(cum[1], diff(cum))
  pmf / sum(pmf)
}
pmf_v_uk <- uk_dose1_pmf(365)

## quick check that the bands match the UK source
cum_v <- cumsum(pmf_v_uk)
cat(sprintf("UK dose-1 check:  <8wk %.1f%% (~14.2)  |  <12wk %.1f%% (~93.8)  |  peak day %d\n",
            100*cum_v[56], 100*cum_v[84], which.max(pmf_v_uk)))

## --- 3.2  His run_one_msccs, with dose days drawn per-infant ----------------
## CHANGED vs his version: dose1_day / dose2_day are now per-infant VECTORS,
## drawn from pmf_v_uk (dose 1) and dose1 + Normal(gap, gap_sd) (dose 2).
## The risk_window is no longer a fixed pair of age-blocks; each infant's
## risk days are relative to THEIR OWN doses. Everything else is his.
run_one_msccs_uk <- function(N = 5e6,
                             OBS_END = 365L,
                             p_death = 2e-4 * 365 / 365,
                             RR = 1.5,
                             RISK_LEN = 3L,
                             agegrp = NULL,
                             pmf_v = pmf_v_uk,
                             gap = 60, gap_sd = 7) {
  days <- 1:OBS_END
  
  ## ---- per-infant dose schedule (THE fix) ----
  n_inf <- N
  d1 <- sample.int(OBS_END, n_inf, replace = TRUE, prob = pmf_v[1:OBS_END])
  d2 <- d1 + round(rnorm(n_inf, gap, gap_sd)); d2 <- pmax(d2, d1 + 1L)
  
  ## ---- baseline hazard = his gamma; apply RR in each infant's OWN windows ----
  ## We build the death-day pmf per infant by tilting the shared gamma pmf in
  ## that infant's two risk windows. Vectorised over infants via a risk-flag
  ## matrix would be huge, so we sample deaths from the MIXTURE: draw an age
  ## from the gamma baseline, then thin/boost by whether it lands in a window.
  gamma_pmf <- our_gamma_pmf(365)          # his baseline (must be in scope)
  base <- gamma_pmf[1:OBS_END]
  
  ## expected deaths (his Poisson count), then assign each death to an infant
  n_deaths <- rpois(1, N * p_death)
  who      <- sample.int(n_inf, n_deaths, replace = TRUE)      # which infant dies
  dd1 <- d1[who]; dd2 <- d2[who]
  
  ## per-death age pmf: base, x RR on [d1,d1+RL) and [d2,d2+RL); sample one age each
  ## (done per death via inverse-CDF on a per-death tilted vector)
  death_day <- integer(n_deaths)
  ## process in chunks to bound memory (n_deaths is small: ~N*p_death)
  for (i in seq_len(n_deaths)) {
    haz <- base
    w1 <- dd1[i]:min(dd1[i]+RISK_LEN-1L, OBS_END); haz[w1] <- haz[w1]*RR
    w2 <- dd2[i]:min(dd2[i]+RISK_LEN-1L, OBS_END); haz[w2] <- haz[w2]*RR
    death_day[i] <- sample.int(OBS_END, 1L, prob = haz)
  }
  
  cases <- data.frame(id = seq_len(n_deaths), death = death_day,
                      d1 = dd1, d2 = dd2)
  cases$ad1 <- ifelse(cases$d1 <= cases$death, cases$d1, NA_integer_)
  cases$ad2 <- ifelse(cases$d2 <= cases$death, cases$d2, NA_integer_)
  
  md <- data.frame(indiv = cases$id, astart = 1L, aend = OBS_END,
                   aevent = cases$death, ad1 = cases$ad1, ad2 = cases$ad2)
  
  fit <- tryCatch(
    eventdepenexp(indiv = indiv, astart = astart, aend = aend, aevent = aevent,
                  adrug = cbind(ad1, ad2), aedrug = cbind(ad1, ad2) + RISK_LEN,
                  expogrp = 0, sameexpopar = TRUE, agegrp = agegrp,
                  dataformat = "multi", data = md),
    error = function(e) NULL)
  
  ## crude O/E vs the UNtilted gamma baseline, summed over each infant's windows
  in_win <- (cases$death >= cases$d1 & cases$death < cases$d1 + RISK_LEN) |
    (cases$death >= cases$d2 & cases$death < cases$d2 + RISK_LEN)
  obs_risk <- sum(in_win)
  ## expected fraction of person-time in-window, averaged over infants
  exp_frac <- mean(vapply(seq_len(n_inf), function(j) {
    w <- unique(c(d1[j]:min(d1[j]+RISK_LEN-1L,OBS_END),
                  d2[j]:min(d2[j]+RISK_LEN-1L,OBS_END)))
    sum(base[w]) / sum(base)
  }, numeric(1))[sample.int(n_inf, min(n_inf, 5000))])   # subsample for speed
  exp_risk <- exp_frac * n_deaths
  oe <- obs_risk / exp_risk
  
  if (is.null(fit)) {
    return(data.frame(n_deaths=n_deaths, obs_risk=obs_risk, exp_risk=exp_risk,
                      OE=oe, IRR=NA, lower=NA, upper=NA, p=NA))
  }
  irr <- fit$conf.int[1,1]; lo <- fit$conf.int[1,3]; hi <- fit$conf.int[1,4]
  se  <- (log(hi) - log(lo)) / (2*qnorm(0.975))
  p   <- 2*pnorm(-abs(log(irr)/se))
  data.frame(n_deaths=n_deaths, obs_risk=obs_risk, exp_risk=exp_risk,
             OE=oe, IRR=irr, lower=lo, upper=hi, p=p)
}

## --- 3.3  His run_msccs_sims, pointed at the UK version --------------------
run_msccs_sims_uk <- function(K = 100, seed = 20260627, agegrp = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  agegrp_label <- if (is.null(agegrp)) "none" else paste(agegrp, collapse = ",")
  out <- vector("list", K)
  for (i in seq_len(K)) {
    #if (i %% 10 == 0) cat("UK sim", i, "of", K, "\n")
    if (i %% 10 == 0) message("UK sim", i, "of", K)
    out[[i]] <- run_one_msccs_uk(agegrp = agegrp, ...)
    out[[i]]$sim <- i
  }
  results <- do.call(rbind, out)
  data.frame(
    n_age_bins = ifelse(is.null(agegrp), 0, length(agegrp)+1),
    mean_deaths = mean(results$n_deaths),
    mean_OE = mean(results$OE, na.rm=TRUE),
    mean_IRR = mean(results$IRR, na.rm=TRUE),
    median_IRR = median(results$IRR, na.rm=TRUE),
    power_p_lt_0.05 = mean(results$p < 0.05, na.rm=TRUE),
    power_CI_excludes_1 = mean(results$lower > 1 | results$upper < 1, na.rm=TRUE),
    agegrp = agegrp_label)
}

## --- 3.4  Re-run HIS bin comparison, UK timing, RR = 1.5 -------------------
## Same five bin sets he used; run_msccs_sims_uk in place of run_msccs_sims.
## Progress: the "UK sim i of K" lines print every 10 sims (inside the sims
## function); here we add a per-bin-set banner and print the result row as
## soon as each bin set completes.
bins <- list(
  none        = NULL,
  monthly     = c(30,60,90,120),
  fortnightly = seq(14,140,14),
  weekly      = seq(7,147,7),
  morris      = c(60,75,90,105,120,145,170,200,250,300)
)

table_uk <- NULL
n_sets <- length(bins)
for (i in seq_along(bins)) {
  cat(sprintf("\n========== bin set %d of %d: '%s' ==========\n",
              i, n_sets, names(bins)[i]))
  sim <- run_msccs_sims_uk(
    K = 100, seed = NULL,
    N = 5e6, OBS_END = 365L,
    p_death = 2e-4 * 365 / 365,
    RR = 1.5, RISK_LEN = 3L,
    agegrp = bins[[i]],
    pmf_v = pmf_v_uk, gap = 60, gap_sd = 7)
  sim$bins <- names(bins)[i]
  #cat(sprintf(">>> DONE '%s' (%d/%d):  mean_IRR=%.3f  power(p<.05)=%.2f  power(CI)=%.2f\n",
   #           names(bins)[i], i, n_sets,
    #          sim$mean_IRR, sim$power_p_lt_0.05, sim$power_CI_excludes_1))
  message(sprintf(">>> DONE '%s' (%d/%d):  mean_IRR=%.3f  power(p<.05)=%.2f  power(CI)=%.2f",
              names(bins)[i], i, n_sets,
              sim$mean_IRR, sim$power_p_lt_0.05, sim$power_CI_excludes_1))
  print(sim[, c("bins","n_age_bins","mean_IRR","power_p_lt_0.05","power_CI_excludes_1")],
        row.names = FALSE)
  table_uk <- rbind(table_uk, sim)
  #flush.console()
}
#cat("\n=== HIS simulation, RR=1.5, with UK dose-timing variability ===\n")
message("=== HIS simulation, RR=1.5, with UK dose-timing variability ===")

print(table_uk[, c("bins","n_age_bins","mean_IRR","power_p_lt_0.05","power_CI_excludes_1")],
      row.names = FALSE)

## --- 3.5  Side-by-side: his point-mass result vs UK-timing result ----------
## (run his original "Simulate age bins" loop first to get `table`)
if (exists("table")) {
  #cat("\n=== POINT-MASS (his) vs UK-VARIABLE (same bins, same RR=1.5) ===\n")
  message("=== POINT-MASS (his) vs UK-VARIABLE (same bins, same RR=1.5) ===")
  
  cmp <- data.frame(
    bins              = names(bins),
    power_pointmass   = table$power_p_lt_0.05,
    power_uk_variable = table_uk$power_p_lt_0.05)
  print(cmp, row.names = FALSE)
} else {
  #cat("\n[3.5 skipped: run his point-mass 'Simulate age bins' loop first to create `table`]\n")
  message("[3.5 skipped: run his point-mass 'Simulate age bins' loop first to create `table`]")
}

## ===========================================================================
## SECTION 4: NULL SIMULATIONS (RR = 1.0) for his bin sets.
## What he should have run alongside his power table. Adds the two columns he
## dropped: type-I error rate and 95% CI coverage. Run for BOTH his point-mass
## timing AND realistic UK timing, so the contrast is visible.
##
## Requires in scope: our_gamma_pmf(), eventdepenexp (SCCS), and from Section 3:
##   pmf_v_uk, run_one_msccs_uk(). His run_one_msccs() must also be sourced.
## ===========================================================================
library(SCCS)

bins <- list(
  none        = NULL,
  monthly     = c(30,60,90,120),
  fortnightly = seq(14,140,14),
  weekly      = seq(7,147,7),
  morris      = c(60,75,90,105,120,145,170,200,250,300)
)

## --- helper: run K sims at a given RR for one bin set, ONE timing arm -------
## `one_fun` is either his run_one_msccs (point mass) or run_one_msccs_uk.
## Returns type-I/power + coverage + mean/median IRR.
run_null_arm <- function(one_fun, bins_vec, RR, K, dots = list()) {
  irr <- lo <- hi <- p <- nc <- numeric(K)
  for (k in seq_len(K)) {
    if (k %% 10 == 0) message(sprintf("   sim %d of %d", k, K))
    args <- c(list(agegrp = bins_vec, RR = RR), dots)
    r <- do.call(one_fun, args)
    irr[k] <- r$IRR; lo[k] <- r$lower; hi[k] <- r$upper
    p[k]   <- r$p;   nc[k] <- r$n_deaths
  }
  data.frame(
    mean_IRR   = mean(irr, na.rm=TRUE),
    median_IRR = median(irr, na.rm=TRUE),
    reject_rate = mean(p < 0.05, na.rm=TRUE),                 # type-I if RR=1, power if RR>1
    ci_coverage = mean(lo <= RR & hi >= RR, na.rm=TRUE),      # the column he omitted
    mean_deaths = mean(nc, na.rm=TRUE),
    n_valid     = sum(!is.na(irr)))
}

## ===========================================================================
## 4a: His exact setup (point mass at 60/120), now WITH a null arm and the
## coverage column. RR=1.0 gives type-I error; RR=1.5 reproduces his power.
## NOTE: his run_one_msccs() takes dose1_day/dose2_day scalars (default 60/120).
## ===========================================================================
set.seed(20260628)
K4 <- 100
pm_dots <- list(N = 5e6, OBS_END = 365L,
                p_death = 2e-4 * 365 / 365, RISK_LEN = 3L,
                dose1_day = 60L, dose2_day = 120L)

null_pm <- NULL
for (RR in c(1.0, 1.5)) {
  for (bn in names(bins)) {
    message(sprintf(">> point-mass  RR=%.1f  bins=%s", RR, bn))
    s <- run_null_arm(run_one_msccs, bins[[bn]], RR, K4, pm_dots)
    null_pm <- rbind(null_pm, cbind(timing="point_mass", RR=RR, bins=bn, s))
    cat(sprintf("   reject=%.2f  coverage=%.2f  meanIRR=%.3f\n",
                s$reject_rate, s$ci_coverage, s$mean_IRR)); flush.console()
  }
}
write.csv(null_pm, "section4a_null_pointmass.csv", row.names = FALSE)

cat("\n=== POINT-MASS: type-I (RR=1.0) and power (RR=1.5), WITH coverage ===\n")
print(null_pm[, c("timing","RR","bins","mean_IRR","reject_rate","ci_coverage")],
      row.names = FALSE)

## ===========================================================================
## 4b: Same null + power, but realistic UK dose timing (run_one_msccs_uk from
## Section 3). This is the apples-to-apples contrast: does restoring timing
## variability fix the type-I inflation AND the coverage collapse?
## ===========================================================================
set.seed(20260628)
uk_dots <- list(N = 5e6, OBS_END = 365L,
                p_death = 2e-4 * 365 / 365, RISK_LEN = 3L,
                pmf_v = pmf_v_uk, gap = 60, gap_sd = 7)

null_uk <- NULL
for (RR in c(1.0, 1.5)) {
  for (bn in names(bins)) {
    message(sprintf(">> UK-variable  RR=%.1f  bins=%s", RR, bn))
    s <- run_null_arm(run_one_msccs_uk, bins[[bn]], RR, K4, uk_dots)
    null_uk <- rbind(null_uk, cbind(timing="uk_variable", RR=RR, bins=bn, s))
    cat(sprintf("   reject=%.2f  coverage=%.2f  meanIRR=%.3f\n",
                s$reject_rate, s$ci_coverage, s$mean_IRR)); flush.console()
  }
}
write.csv(null_uk, "section4b_null_ukvariable.csv", row.names = FALSE)

cat("\n=== UK-VARIABLE: type-I (RR=1.0) and power (RR=1.5), WITH coverage ===\n")
print(null_uk[, c("timing","RR","bins","mean_IRR","reject_rate","ci_coverage")],
      row.names = FALSE)

## ===========================================================================
## 4c: Assemble the honest version of his table. For each timing arm and bin
## set: null type-I error, null coverage, power, and the bias at RR=1.5.
## A valid method needs type-I near 0.05 AND coverage near 0.95 — THEN power
## is meaningful. High power with type-I ~1 and coverage ~0 is a broken test.
## ===========================================================================
both <- rbind(null_pm, null_uk)

honest <- do.call(rbind, lapply(split(both, list(both$timing, both$bins), drop=TRUE),
                                function(d) {
                                  n  <- d[d$RR==1.0, ]; a <- d[d$RR==1.5, ]
                                  data.frame(timing = n$timing, bins = n$bins,
                                             null_typeI    = n$reject_rate,    # want ~0.05
                                             null_coverage = n$ci_coverage,    # want ~0.95
                                             null_meanIRR  = n$mean_IRR,       # want ~1.00
                                             power_RR1.5   = a$reject_rate,    # only meaningful if null is calibrated
                                             bias_RR1.5    = a$mean_IRR)       # want ~1.50
                                }))
honest <- honest[order(honest$timing, match(honest$bins, names(bins))), ]
rownames(honest) <- NULL
write.csv(honest, "section4c_honest_table.csv", row.names = FALSE)

cat("\n=== THE TABLE HE SHOULD HAVE PUBLISHED ===\n")
cat("(valid = null_typeI ~0.05, null_coverage ~0.95, null_meanIRR ~1.0,\n")
cat(" THEN power is interpretable)\n\n")
print(honest, row.names = FALSE)


## ===========================================================================
## plot_jikky_age_pmf.R   (revised)
## Background-age SIDS-risk pmf with dose-age overlays and age-bin shading.
##
## Figures produced (to a multi-page PDF):
##   Fig 1  background age pmf alone
##   Fig 2  STOCHASTIC overall (no bins): dose-1 & dose-2 marginal pmfs
##          PLUS per-case dose-1->dose-2 timings
##   Fig 3  DETERMINISTIC overall (no bins): dose 1 @ 60d, dose 2 @ 120d
##   Page 4a  DETERMINISTIC + 4 bin sets, 2x2
##   Page 4b  STOCHASTIC   + 4 bin sets, 2x2
## Each 2x2 page carries a two-line outer title:
##   line 1  "Background Age SIDS Risk"
##   line 2  "Deterministic Dose Ages (60d/120d)" | "Stochastic Dose Ages (UK data)"
##
## Base R only. Source after your simulation objects are in the workspace.
## ===========================================================================

## ---------------------------------------------------------------------------
## 0.  INPUTS
## ---------------------------------------------------------------------------
MAXDAY    <- 365L      # age axis upper limit (use 729L for 2 yr)

## dose-2 spacing for the STOCHASTIC model -----------------------------------
## YOUR model: dose 1 ~ UK distribution; dose 2 = dose 1 + Normal(60, 7).
## (Each infant's dose 2 is 60 d after their own dose 1, plus/minus variance.)
GAP2_MEAN <- 60
GAP2_SD   <- 7

## (A) JIKKY'S SMOOTH BACKGROUND AGE PMF -------------------------------------
##     Per-day pmf over 1..MAXDAY (sums to 1).
##     >>> set `dens_jikky <- <his fitted vector>` and delete the fallback.
jikky_age_pmf <- function(maxday = MAXDAY, shape = 2.7, scale = 35) {
  d <- seq_len(maxday); w <- dgamma(d, shape = shape, scale = scale); w / sum(w)
}
if (!exists("dens_jikky")) {
  dens_jikky <- jikky_age_pmf()
  message("NOTE: PLACEHOLDER gamma baseline in use; set `dens_jikky` to Jikky's ",
          "actual fitted pmf before publishing.")
}

## (B) UK DOSE-1 AGE DISTRIBUTION (per-day pmf) ------------------------------
uk_dose1_pmf <- function(maxday = MAXDAY, p_ever = 0.986, shift = 6,
                         meanlog = 1.1222, sdlog = 0.4039) {
  d  <- 1:maxday
  Fd <- function(x){ wk <- x/7; ifelse(wk > shift, p_ever*plnorm(wk-shift, meanlog, sdlog), 0) }
  cum <- Fd(d); pmf <- c(cum[1], diff(cum)); pmf / sum(pmf)
}
if (!exists("pmf_v_uk")) pmf_v_uk <- uk_dose1_pmf()

## (C) DETERMINISTIC DOSE DAYS -----------------------------------------------
DET_DOSES <- c(60, 120)

## (D) BIN SETS (interior breakpoints, verbatim from the results table) ------
bin_sets <- list(
  monthly     = c(30,60,90,120),
  fortnightly = c(14,28,42,56,70,84,98,112,126,140),
  weekly      = c(7,14,21,28,35,42,49,56,63,70,77,84,91,98,105,112,119,126,133,140,147),
  Morris      = c(60,75,90,105,120,145,170,200,250,300)
)

## ---------------------------------------------------------------------------
## 1.  STOCHASTIC dose ages: use your per-infant vectors if present, else draw
##     dose1 ~ UK pmf, dose2 = dose1 + Normal(GAP2_MEAN, GAP2_SD).
## ---------------------------------------------------------------------------
pmf_from_samp <- function(x, maxday = MAXDAY) tabulate(pmin(pmax(round(x),1),maxday), maxday)/length(x)
smooth_pmf <- function(p, k = 5) {
  f <- stats::filter(p, rep(1/k, k), sides = 2); f[is.na(f)] <- p[is.na(f)]; as.numeric(f/sum(f))
}
make_stochastic_doses <- function(maxday = MAXDAY, n = 20000L, pmf_d1 = pmf_v_uk,
                                  gap_mean = GAP2_MEAN, gap_sd = GAP2_SD, seed = 1L) {
  if (exists("dose1_age") && exists("dose2_age")) {     # your real simulated vectors
    d1 <- dose1_age; d2 <- dose2_age
  } else {
    set.seed(seed)
    d1 <- sample(seq_len(maxday), n, replace = TRUE, prob = pmf_d1)
    d2 <- pmin(pmax(round(d1 + rnorm(n, gap_mean, gap_sd)), 1L), maxday)
  }
  list(d1 = d1, d2 = d2,
       pmf1 = smooth_pmf(pmf_from_samp(d1, maxday)),
       pmf2 = smooth_pmf(pmf_from_samp(d2, maxday)))
}
stoch <- make_stochastic_doses()
stoch_overlays <- list(
  list(pmf = pmf_v_uk,   lab = "dose-1 age", col = "#1B7837"),   # green
  list(pmf = stoch$pmf2, lab = "dose-2 age", col = "#D95F02")    # orange
)
stoch_pc <- list(d1 = stoch$d1, d2 = stoch$d2)

## ---------------------------------------------------------------------------
## 2.  HELPERS: bin shading, and per-case dose-timing strip
## ---------------------------------------------------------------------------
add_bins <- function(breaks, xlim = c(1, MAXDAY),
                     col_a = "#F7E2E2", col_b = "#E2EAF7",
                     border = "grey60", lty = 3, lwd = 0.8) {
  br <- breaks[breaks > xlim[1] & breaks < xlim[2]]
  edges <- sort(unique(c(xlim[1], br, xlim[2])))
  usr <- par("usr")
  for (i in seq_len(length(edges) - 1L))
    rect(edges[i], usr[3], edges[i + 1L], usr[4],
         col = if (i %% 2L == 1L) col_a else col_b, border = NA)
  abline(v = br, col = border, lty = lty, lwd = lwd)
}

## per-case strip: a sample of infants, each a faint dose1->dose2 segment with
## colored endpoints, drawn in a band near the bottom (left-axis coords).
add_per_case <- function(d1, d2, n = 110L, col1 = "#1B7837", col2 = "#D95F02", seed = 7L) {
  usr <- par("usr"); ymax <- usr[4]
  set.seed(seed); idx <- sample(length(d1), min(n, length(d1)))
  ys <- runif(length(idx), 0.015 * ymax, 0.11 * ymax)
  segments(d1[idx], ys, d2[idx], ys, col = adjustcolor("grey45", 0.30), lwd = 0.6)
  points(d1[idx], ys, pch = 16, cex = 0.35, col = adjustcolor(col1, 0.75))
  points(d2[idx], ys, pch = 16, cex = 0.35, col = adjustcolor(col2, 0.75))
}

## ---------------------------------------------------------------------------
## 3.  MAIN PLOT FUNCTION
## ---------------------------------------------------------------------------
plot_age_pmf <- function(dens, maxday = MAXDAY,
                         overlays = NULL,        # list(list(pmf=,lab=,col=), ...) on right axis
                         doses = NULL, dose_lab = paste0("dose ", seq_along(doses)),
                         per_case = NULL,        # list(d1=, d2=) -> per-case strip
                         bins = NULL, bin_name = NULL,
                         main = "", bg_col = "#1F6FB2") {
  day <- seq_len(maxday); d <- dens[day]
  plot(day, d, type = "n", xlim = c(1, maxday), ylim = c(0, max(d) * 1.12),
       xlab = "Age (days)", ylab = "Background SIDS-risk pmf (prop./day)",
       main = main, xaxs = "i", yaxs = "i")
  if (!is.null(bins)) add_bins(bins, xlim = c(1, maxday))
  polygon(c(day[1], day, day[maxday]), c(0, d, 0), col = adjustcolor(bg_col, 0.18), border = NA)
  lines(day, d, lwd = 2.2, col = bg_col)
  
  if (!is.null(per_case)) add_per_case(per_case$d1, per_case$d2)
  
  if (!is.null(doses)) {
    abline(v = doses, col = "#B2182B", lwd = 2)
    text(doses, par("usr")[4] * 0.96, labels = dose_lab, col = "#B2182B", pos = 4, cex = 0.75, xpd = NA)
  }
  
  if (!is.null(overlays)) {
    rmax <- max(vapply(overlays, function(o) max(o$pmf[day]), numeric(1)))
    par(new = TRUE)
    plot(day, overlays[[1]]$pmf[day], type = "n", axes = FALSE, xlab = "", ylab = "",
         xlim = c(1, maxday), ylim = c(0, rmax * 1.12), xaxs = "i", yaxs = "i")
    for (o in overlays) lines(day, o$pmf[day], col = o$col, lwd = 2.2)
    axis(4, col = "grey30", cex.axis = 0.8)
    mtext("dose-age pmf (prop./day)", side = 4, line = 2.3, cex = 0.7)
  }
  
  ## legend
  lt <- "background SIDS-risk pmf"; lc <- bg_col; ll <- 1; lw <- 2.2; lp <- NA
  if (!is.null(overlays)) for (o in overlays) { lt <- c(lt,o$lab); lc <- c(lc,o$col); ll <- c(ll,1);  lw <- c(lw,2.2); lp <- c(lp,NA) }
  if (!is.null(per_case)) { lt <- c(lt,"per-case dose1->dose2"); lc <- c(lc,"grey45"); ll <- c(ll,1); lw <- c(lw,0.8); lp <- c(lp,16) }
  if (!is.null(doses))    { lt <- c(lt,"deterministic doses");   lc <- c(lc,"#B2182B"); ll <- c(ll,1); lw <- c(lw,2);   lp <- c(lp,NA) }
  if (!is.null(bins))     { lt <- c(lt,paste0(bin_name," bins")); lc <- c(lc,"grey60"); ll <- c(ll,3); lw <- c(lw,0.8); lp <- c(lp,NA) }
  legend("topright", legend = lt, col = lc, lty = ll, lwd = lw, pch = lp,
         bty = "n", cex = 0.70, inset = 0.01, seg.len = 1.6)
  box()
}

## ---------------------------------------------------------------------------
## 4.  2x2 BIN PAGE (one dosing model, 4 bin sets) with two-line outer title
## ---------------------------------------------------------------------------
make_bin_page <- function(kind = c("deterministic", "stochastic")) {
  kind <- match.arg(kind)
  op <- par(mfrow = c(2, 2), oma = c(0, 0, 5, 0), mar = c(4, 4, 2.4, 4))
  for (nm in names(bin_sets)) {
    if (kind == "deterministic")
      plot_age_pmf(dens_jikky, doses = DET_DOSES, dose_lab = c("dose 1","dose 2"),
                   bins = bin_sets[[nm]], bin_name = nm, main = paste0(nm, " bins"))
    else
      plot_age_pmf(dens_jikky, overlays = stoch_overlays,
                   bins = bin_sets[[nm]], bin_name = nm, main = paste0(nm, " bins"))
  }
  sub <- if (kind == "deterministic") "Deterministic Dose Ages (60d/120d)" else "Stochastic Dose Ages (UK data)"
  mtext("Background Age SIDS Risk", outer = TRUE, line = 2.6, cex = 1.25, font = 2)
  mtext(sub,                        outer = TRUE, line = 0.8, cex = 1.00, font = 2)
  par(op)
}

## ===========================================================================
## 5.  DRIVER
## ===========================================================================
pdf("jikky_age_bins.pdf", width = 11, height = 8.5)   # or dev.new() per figure

## Fig 1: background alone
plot_age_pmf(dens_jikky, main = "Background Age SIDS Risk")

## Fig 2: STOCHASTIC overall (no bins) -- dose1 & dose2 pmfs + per-case timing
plot_age_pmf(dens_jikky, overlays = stoch_overlays, per_case = stoch_pc,
             main = "Stochastic Dose Ages (UK data): dose 1 & dose 2")

## Fig 3: DETERMINISTIC overall (no bins)
plot_age_pmf(dens_jikky, doses = DET_DOSES, dose_lab = c("dose 1 (60d)","dose 2 (120d)"),
             main = "Deterministic Dose Ages (60d / 120d)")

## Page 4a / 4b: each dosing model x 4 bin sets, 2x2
make_bin_page("deterministic")
make_bin_page("stochastic")

dev.off()

## ---------------------------------------------------------------------------
## Notes
## - Overlays (dose-1/dose-2 pmfs) share a RIGHT axis so they keep their true
##   per-day scale; left axis is always the background SIDS-risk pmf.
## - Per-case strip appears only in Fig 2 (the no-bins overall stochastic plot):
##   each faint segment is one infant's dose-1 (green) -> dose-2 (orange) ages.
## - dose2 spacing is GAP2_MEAN/GAP2_SD; or it uses your `dose1_age`/`dose2_age`
##   vectors directly if they exist in the workspace.
## - Swap the placeholder `dens_jikky` for Jikky's actual fitted pmf.
## ===========================================================================