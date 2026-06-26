# ===========================================================================
# Age-distribution models + diagnostics for the multidose SIDS/SCCS simulation.
#   dose_distribution()   -> UK DTP MULTIDOSE schedule (dose 1/2/3) + pmf_v
#                            + gen_schedule(N): the single source of dose timing
#   sids_distribution()   -> SIDS age-at-death null (LOWESS of Traversa) + pmf_s
#   plot_sim_diagnostics()-> 4 diagnostic plots (separate pages) from a dataset
# pmf_v / pmf_s and gen_schedule() feed directly into the simulation, so dose
# generation is defined ONCE and never duplicated.
# ===========================================================================

# ===========================================================================
# 1. DOSE DISTRIBUTION  (UK DTP, three-dose primary schedule)
# ---------------------------------------------------------------------------
# Dose 1: ever-vacc ~ Bernoulli(0.986); age|vacc = 6 + Lognormal(1.1222,0.4039) wk,
#   fitted to MCS Welsh cohort (born 2000-01, N=1,782): <8wk 14.2% | 8-12wk 79.6%
#   | >12wk 4.8% | never 1.4%.  Cottrell et al. medRxiv 2022.04.04.22273336 ;
#   Hungerford et al. Vaccine 2016 (PMC5720480).  UK schedule due 8/12/16 wk.
# Doses 2,3: each GAP days (~2 months) after the previous, +/- GAP_SD (Normal).
# Returns pmf_v (dose-1 per-day density), gen_schedule(N) -> d1/d2/d3 (days),
#   and per-day densities of each dose + all-doses (for plotting).
# ===========================================================================
dose_distribution <- function(GAP = 60, GAP_SD = 7, maxday = 365L, Nmc = 2e6) {
  p_ever <- 0.986; shift <- 6; meanlog <- 1.1222; sdlog <- 0.4039
  days   <- 1:maxday
  F_days <- function(d){ wk <- d/7; ifelse(wk > shift, p_ever*plnorm(wk-shift,meanlog,sdlog), 0) }
  cum_v  <- F_days(days)
  pmf_v  <- c(cum_v[1], diff(cum_v)); pmf_v <- pmf_v/sum(pmf_v)   # dose-1 per-day density
  
  ## single source of dose-schedule generation (reused by the simulation)
  gen_schedule <- function(N) {
    d1 <- sample(days, N, replace = TRUE, prob = pmf_v)
    d2 <- pmax(d1 + round(rnorm(N, GAP, GAP_SD)), d1 + 1L)
    d3 <- pmax(d2 + round(rnorm(N, GAP, GAP_SD)), d2 + 1L)
    data.frame(d1 = d1, d2 = d2, d3 = d3)
  }
  ## per-day densities (dose1/2/3 + all) via MC, for plotting/overlay
  s  <- gen_schedule(Nmc)
  dd <- function(x){ x <- x[x>=1 & x<=maxday]; as.numeric(table(factor(x, levels=days)))/Nmc }
  pd1 <- dd(s$d1); pd2 <- dd(s$d2); pd3 <- dd(s$d3)
  
  list(params = c(p_ever=p_ever, shift=shift, meanlog=meanlog, sdlog=sdlog,
                  GAP=GAP, GAP_SD=GAP_SD),
       days = days, pmf_v = pmf_v, gen_schedule = gen_schedule,
       dose1 = pd1, dose2 = pd2, dose3 = pd3, all_doses = pd1+pd2+pd3,
       cum_v = cum_v)
}

# ===========================================================================
# 2. SIDS DISTRIBUTION  (theoretical null age-at-death)
# ---------------------------------------------------------------------------
# Nonparametric (LOWESS) smooth of the digitized envelope of Traversa et al.,
#   PLoS ONE 2011;6(1):e16363, Fig 3 (604 SUD events). Year-1-normalized.
# Returns pmf_s (per-day density, sums to 1 over 1..365), cum_s, sampler,
#   and the anchors + fitted curve (count units) for the fit-check plot.
# NOTE: You can substitute whatever presumed background incidence by age effect
#       you want -- just be sure to output the list with pmf_s being vector of
#       lenght 365 indicating background incidence for that age in days.
# ===========================================================================
sids_distribution <- function(maxday = 365L, span = 0.2) {
  anchor_day <- c(6,13,20,27,34,41,48,55,62,69,76,83,90,97,104,111,118,125,132,139,
                  146,153,160,167,174,181,188,195,205,220,240,260,280,300,330,360,
                  400,440,480,520,560,600,640,680,715)
  anchor_cnt <- c(12,16,20,24,31,23,30,28,33,28,37,24,21,24,23,16,15,14,12,18,
                  21,13,11,9,8,7,6,5,4,3,2.6,2.2,1.9,1.7,1.5,1.4,1.2,1.1,1.0,0.9,
                  0.9,0.8,0.9,0.9,0.9)
  lo      <- lowess(anchor_day, anchor_cnt, f = span)
  fit_cnt <- pmax(approx(lo$x, lo$y, xout = 1:maxday, rule = 2)$y, 0)   # events/7-day-bin units
  pmf_s   <- fit_cnt / sum(fit_cnt)                                     # per-day density
  list(days = 1:maxday, pmf_s = pmf_s, cum_s = cumsum(pmf_s),
       rage_sids = function(n) sample(1:maxday, n, replace = TRUE, prob = pmf_s),
       anchor_day = anchor_day, anchor_cnt = anchor_cnt, fit_cnt = fit_cnt)
}

# ===========================================================================
# 3. plot_theory_diagnostics() — 4 diagnostic plots (separate pages) for the
# THEORETICAL distributions, to plug in after dose_distribution() /
# sids_distribution() are built.
# ---------------------------------------------------------------------------
# Pages: (1) doses by age           (2) SIDS by age
#        (3) doses & SIDS overlaid  (4) theoretical SIDS RATE by days since dose
# Plot 4 is the EXPECTED rate = expected hazard / person-time, integrated over
# the schedule distribution (which dose is most-recent at each age), with the
# VAX_effect multiplier on days 0..RISK_LEN-1. Under VAX=1 it shows the smooth
# age-confounding decline; VAX>1 adds the 0..RISK_LEN-1 spike on top.
#
# ARGS: dose = dose_distribution(); sids = sids_distribution();
#       VAX (effect to display, default 1 = null), RISK_LEN, since_max,
#       Nmc (schedule draws for plot 4).
# Tip: wrap in pdf("theory.pdf"); plot_theory_diagnostics(dose,sids); dev.off()
# ===========================================================================
plot_theory_diagnostics <- function(dose, sids, VAX = 1.0, RISK_LEN = 3L,
                                    since_max = 60L, Nmc = 2e5) {
  days   <- dose$days; maxday <- length(days)
  pmf_s  <- sids$pmf_s
  all_d  <- dose$all_doses                      # expected doses/child/day (sums ~3)
  
  ## --- plot 4 input: theoretical rate by days-since-most-recent-dose ---
  s  <- dose$gen_schedule(Nmc)
  D  <- list(s$d1, s$d2, s$d3); nextD <- list(s$d2, s$d3, rep(maxday + 1L, Nmc))
  PT <- numeric(since_max + 1L); EV <- numeric(since_max + 1L)
  for (k in 1:3) {
    dk <- D[[k]]; nk <- nextD[[k]]
    for (off in 0:since_max) {
      a     <- dk + off
      valid <- a <= maxday & a < nk             # dose k is the most-recent dose at age a
      mult  <- if (off < RISK_LEN) VAX else 1.0
      PT[off + 1L] <- PT[off + 1L] + sum(valid)
      EV[off + 1L] <- EV[off + 1L] + sum(valid * pmf_s[pmax(pmin(a, maxday), 1L)] * mult)
    }
  }
  since <- 0:since_max; rate <- ifelse(PT > 0, EV / PT, NA)
  
  op <- par(mfrow = c(1,1), mar = c(4.5,4.5,3,1)); on.exit(par(op))
  # (1) doses by age (all doses + components)
  plot(days, 100*all_d, type="l", lwd=2.5, col="#1F6FB2", xlim=c(0,maxday),
       xlab="Age (days)", ylab="% of children dosed on day", main="Doses by age (theoretical)")
  lines(days, 100*dose$dose1, lwd=1.4, lty=2, col="#1F6FB2")
  lines(days, 100*dose$dose2, lwd=1.4, lty=4, col="#2E86C1")
  lines(days, 100*dose$dose3, lwd=1.4, lty=3, col="#5DADE2")
  legend("topright", c("all doses","dose 1","dose 2","dose 3"),
         col=c("#1F6FB2","#1F6FB2","#2E86C1","#5DADE2"),
         lwd=c(2.5,1.4,1.4,1.4), lty=c(1,2,4,3), bty="n")
  # (2) SIDS by age
  plot(days, 100*pmf_s, type="l", lwd=2.5, col="#C0392B", xlim=c(0,maxday),
       xlab="Age (days)", ylab="% of SIDS on day", main="SIDS by age (theoretical)")
  polygon(c(days, rev(days)), c(100*pmf_s, rep(0,maxday)),
          col=adjustcolor("#C0392B",.2), border=NA)
  # (3) overlay (area-normalized; blue=doses, red=SIDS)
  dn <- all_d/sum(all_d)
  plot(days, 100*dn, type="l", lwd=2.2, col="#1F6FB2", xlim=c(0,maxday),
       ylim=c(0, max(100*c(dn, pmf_s))), xlab="Age (days)", ylab="% on each day",
       main="Doses vs SIDS by age (overlaid)")
  lines(days, 100*pmf_s, lwd=2.2, col="#C0392B")
  legend("topright", c("vaccination (all doses)","SIDS deaths"),
         col=c("#1F6FB2","#C0392B"), lwd=2.2, bty="n")
  # (4) theoretical SIDS rate by days since most recent dose
  plot(since, rate, type="l", lwd=2.5, col="#8E44AD",
       xlab="Days since most recent dose", ylab="Expected SIDS rate (per person-day)",
       main=sprintf("SIDS rate by days since dose (theoretical, VAX=%.2g)", VAX))
  abline(v=RISK_LEN-0.5, lty=2, col="grey50")   # end of the 0..RISK_LEN-1 risk window
  
  invisible(list(doses_by_age = data.frame(day=days, all=all_d,
                                           dose1=dose$dose1, dose2=dose$dose2, dose3=dose$dose3),
                 sids_by_age = data.frame(day=days, pmf_s=pmf_s),
                 rate_since_dose = data.frame(days_since=since, person_time=PT,
                                              exp_events=EV, rate=rate)))
}

  # plot_theory_diagnostics_age() — as plot_theory_diagnostics(), but overlays
  # the SCCS age-bin structure on the AGE-axis plots (1-3): vertical lines at each
  # bin boundary, with bins shaded in alternating very-light red / very-light blue.
  # Plot 4 (days-since-dose axis) is NOT age-shaded — age bins don't live on that
  # axis; it keeps the RISK_LEN risk-window boundary instead.
  #
  # NEW ARG: age_bins = interior cut points (e.g. P$age_bins), in days of life.
  #          Edges 0 and maxday are added automatically to close the end bins.
  # ===========================================================================
plot_SIDS_age <- function(dose, sids, age_bins, title, VAX = 1.0,
                                        RISK_LEN = 3L, since_max = 60L, Nmc = 2e5) {
  days   <- dose$days; maxday <- length(days)
  pmf_s  <- sids$pmf_s
  all_d  <- dose$all_doses
  
  ## ---- local helper: shade alternating age bins + draw boundary lines -------
  ## Call AFTER plot(..., type="n") so par("usr") is valid; data goes on top.
  shade_age_bins <- function(col_odd = "#C0392B", col_even = "#1F6FB2",
                             alpha = 0.07) {
    b   <- sort(unique(age_bins[age_bins > 0 & age_bins < maxday]))
    edg <- c(0, b, maxday)                    # region edges
    usr <- par("usr"); yb <- usr[3]; yt <- usr[4]
    for (i in seq_len(length(edg) - 1L)) {
      col <- if (i %% 2L == 1L) adjustcolor(col_odd, alpha) else adjustcolor(col_even, alpha)
      rect(edg[i], yb, edg[i + 1L], yt, col = col, border = NA)
    }
    abline(v = b, col = "grey55", lty = 3, lwd = 0.8)   # boundary at each cut point
    box()                                                # redraw frame over shading
  }
  
  ## --- plot 4 input: theoretical rate by days-since-most-recent-dose ---
  s  <- dose$gen_schedule(Nmc)
  D  <- list(s$d1, s$d2, s$d3); nextD <- list(s$d2, s$d3, rep(maxday + 1L, Nmc))
  PT <- numeric(since_max + 1L); EV <- numeric(since_max + 1L)
  for (k in 1:3) {
    dk <- D[[k]]; nk <- nextD[[k]]
    for (off in 0:since_max) {
      a     <- dk + off
      valid <- a <= maxday & a < nk
      mult  <- if (off < RISK_LEN) VAX else 1.0
      PT[off + 1L] <- PT[off + 1L] + sum(valid)
      EV[off + 1L] <- EV[off + 1L] + sum(valid * pmf_s[pmax(pmin(a, maxday), 1L)] * mult)
    }
  }
  since <- 0:since_max; rate <- ifelse(PT > 0, EV / PT, NA)
  
  op <- par(mfrow = c(1,1), mar = c(4.5,4.5,3,1)); on.exit(par(op))
  png(paste("SIDS_age_", title, "_%d.png", sep = ""), width = 1600, height = 1100, res = 200)  
  # (2) SIDS by age
  plot(days, 100*pmf_s, type="n", xlim=c(0,maxday),
       xlab="Age (days)", ylab="% of SIDS on day", 
       main=paste("SIDS by age (theoretical),",title,"age bands"))
  shade_age_bins()
  polygon(c(days, rev(days)), c(100*pmf_s, rep(0,maxday)),
          col=adjustcolor("#C0392B",.25), border=NA)
  lines(days, 100*pmf_s, lwd=2.5, col="#C0392B")
  box()
  dev.off()
}

# ===========================================================================
# 4. DIAGNOSTICS FROM A SIMULATED DATASET  -> 4 plots, separate pages
# ---------------------------------------------------------------------------
# sim: data.frame with columns death, d1, d2, d3 (NA = dose not received).
# Pages: (1) observed doses by age  (2) observed SIDS by age
#        (3) the two overlaid (area-normalized; blue=doses, red=SIDS)
#        (4) SIDS RATE by days since most recent dose (events / person-time).
# Tip: wrap in pdf("diag.pdf"); plot_sim_diagnostics(sim); dev.off() for pages,
#      or set par(ask=TRUE) for sequential on-screen pages.
# ===========================================================================
plot_sim_diagnostics <- function(sim, maxday = 365L, since_max = 60L) {
  days  <- 1:maxday
  doses <- as.matrix(sim[, c("d1","d2","d3")]); death <- sim$death
  
  ## per-age counts
  dvec <- doses[!is.na(doses) & doses >= 1 & doses <= maxday]
  vax  <- as.numeric(table(factor(dvec, levels = days)))                 # doses by age
  sids <- as.numeric(table(factor(death[death>=1 & death<=maxday], levels = days)))
  
  ## SIDS rate by days since most recent dose (events / person-time)
  PT <- numeric(since_max+1L); EV <- numeric(since_max+1L)
  for (i in seq_len(nrow(sim))) {
    dk <- sort(doses[i, !is.na(doses[i, ])]); if (!length(dk)) next
    d <- death[i]; if (d < dk[1]) next
    for (a in dk[1]:d) { s <- a - dk[max(which(dk <= a))]
    if (s <= since_max) { if (a==d) EV[s+1L] <- EV[s+1L]+1 else PT[s+1L] <- PT[s+1L]+1 } }
  }
  since <- 0:since_max; rate <- ifelse(PT>0, EV/PT, NA)
  
  op <- par(mfrow = c(1,1), mar = c(4.5,4.5,3,1)); on.exit(par(op))
  # (1) doses by age
  plot(days, vax, type="h", col="#1F6FB2", xlim=c(0,maxday),
       xlab="Age (days)", ylab="doses given", main="Observed vaccination by age")
  # (2) SIDS by age
  plot(days, sids, type="h", col="#C0392B", xlim=c(0,maxday),
       xlab="Age (days)", ylab="SIDS deaths", main="Observed SIDS by age")
  # (3) overlay (area-normalized so timing is comparable)
  plot(days, vax/sum(vax)*100, type="l", lwd=2, col="#1F6FB2", xlim=c(0,maxday),
       ylim=c(0, max(vax/sum(vax), sids/sum(sids))*100),
       xlab="Age (days)", ylab="% on each day", main="Doses vs SIDS by age (overlaid)")
  lines(days, sids/sum(sids)*100, lwd=2, col="#C0392B")
  legend("topright", c("doses","SIDS"), col=c("#1F6FB2","#C0392B"), lwd=2, bty="n")
  # (4) SIDS rate by days since most recent dose
  plot(since, rate, type="h", lwd=2, col="#8E44AD",
       xlab="Days since most recent dose", ylab="SIDS rate (events / person-day)",
       main="SIDS rate by days since dose")
  abline(v=2.5, lty=2, col="grey50")   # end of the 0-2 day risk window
  
  invisible(list(vax_by_age = data.frame(day=days, doses=vax),
                 sids_by_age = data.frame(day=days, deaths=sids),
                 rate_since_dose = data.frame(days_since=since, events=EV,
                                              person_time=PT, rate=rate)))
}

# ===========================================================================
# INITIALIZATION  (dose generation defined once, reused everywhere)
# ===========================================================================
dose <- dose_distribution()          # pmf_v + gen_schedule + dose densities
sids <- sids_distribution()          # pmf_s
pmf_v <- dose$pmf_v; pmf_s <- sids$pmf_s
gen_schedule <- dose$gen_schedule    # <- the simulation should call THIS (no duplication)

## verification of properties of vaccination and SIDS distributions for chosen settings
cat(sprintf("VACC: by 8wk %.1f%% (t14.2) | 12wk %.1f%% (t93.8) | dose peaks %d/%d/%d\n",
            100*dose$cum_v[56], 100*dose$cum_v[84],
            which.max(dose$dose1), which.max(dose$dose2), which.max(dose$dose3)))
cat(sprintf("SIDS: median %d d | peak %d d\n",
            which(sids$cum_s>=0.5)[1], which.max(pmf_s)))

## save the per-day distributions
write.csv(data.frame(day=1:365, age_weeks=round((1:365)/7,2),
                     vacc_pct_on_day=round(100*pmf_v,5), vacc_cum_pct=round(100*dose$cum_v,3),
                     sids_pct_on_day=round(100*pmf_s,5), sids_cum_pct=round(100*sids$cum_s,3)),
          "dist_by_day_1_365.csv", row.names = FALSE)

## ---- usage (after the theoretical distributions are generated) ----
 dose <- dose_distribution(); sids <- sids_distribution()
 plot_theory_diagnostics(dose, sids, VAX = 1.0)   # null: smooth age-confounding decline
 plot_theory_diagnostics(dose, sids, VAX = 1.5)   # shows the 0-2 day spike on top

## diagnostics on an individual simulated dataset (uses gen_schedule internally):
 sim <- simulate_multidose(N = 2e6, VAX = 2.0)   # from the simulation script
 plot_sim_diagnostics(sim)                        # 4 pages
 
 
 library(SCCS)
 # ===========================================================================
 # UNIFIED MULTIDOSE TERMINAL-EVENT SCCS — illustration + simulation study.
 # Requires (sourced from the distribution/diagnostics block):
 #   dose_distribution(), sids_distribution(), plot_sim_diagnostics(),
 #   plot_theory_diagnostics().
 # Dose schedule + age distributions come from the theoretical objects, so
 # nothing is duplicated.  Unvaccinated SIDS (death < dose 1) are EXCLUDED from
 # the SCCS; only their count is kept, to report total population SIDS.
 # ===========================================================================
 
 dose <- dose_distribution()      # pmf_v + gen_schedule(N) (d1/d2/d3) + dose densities
 sids <- sids_distribution()      # pmf_s (SIDS age-at-death, null)
 pmf_v <- dose$pmf_v; pmf_s <- sids$pmf_s; gen_schedule <- dose$gen_schedule
 
 
 ## ===========================================================================
 ## CORE FUNCTIONS
 ## ===========================================================================
 ## simulate one population -> ALL SIDS cases (full cohort) + SIDS counts.
 ## Hazard h(a)=p_SIDS*pmf_s(a), x VAX on [d_k, d_k+RISK_LEN-1] for each dose.
 ## Dose coding (the eventdepenexp convention): a dose is NA unless it was
 ## actually received, i.e. d_k <= death AND d_k <= OBS_END. Therefore:
 ##   - pre-dose-1 deaths (death < d1) are KEPT with all three doses = NA
 ##   - any dose scheduled after death is curtailed -> NA
 ## Keeping the pre-dose-1 deaths repopulates the pre-vaccination period with
 ## real events, removing the immortal-time artifact that astart=1 otherwise
 ## creates on a vaccinated-only cohort.
 ## ===========================================================================
 simulate_population <- function(N, VAX, pmf_s, gen_schedule, P) {
   RL <- P$RISK_LEN; maxday <- P$OBS_END; p_SIDS <- P$p_SIDS
   base <- p_SIDS * pmf_s; bcum <- c(0, cumsum(base))
   s <- gen_schedule(N); D <- cbind(s$d1, s$d2, s$d3)
   wmass <- function(dk){ hi <- pmin(dk+RL-1L,maxday); m <- bcum[hi+1L]-bcum[pmax(dk,1L)]
   m[dk<1 | dk>maxday] <- 0; m }
   Wd <- cbind(wmass(D[,1]), wmass(D[,2]), wmass(D[,3])); Wm <- rowSums(Wd)
   Ptot <- p_SIDS + (VAX-1)*Wm
   die <- which(runif(N) < Ptot)
   if (!length(die)) return(list(cases=NULL, n_before=0L, n_after=0L, n_total=0L))
   nC <- length(die); Dc <- D[die,,drop=FALSE]; Wdc <- Wd[die,,drop=FALSE]
   Wmc <- Wm[die]; Ptc <- Ptot[die]
   inwin <- runif(nC) < (VAX*Wmc)/Ptc; death <- integer(nC)
   iw <- which(inwin)
   if (length(iw)) {                                  # in-window: pick dose, then day
     cc <- t(apply(Wdc[iw,,drop=FALSE]/Wmc[iw], 1, cumsum))
     dk <- Dc[cbind(iw, max.col(cc >= runif(length(iw)), ties.method="first"))]
     cand <- outer(dk, 0:(RL-1L), "+"); wb <- matrix(0, nrow(cand), ncol(cand))
     ok <- cand <= maxday; wb[ok] <- base[cand[ok]]
     cb <- t(apply(wb,1,cumsum)); cb <- cb/cb[,RL]
     death[iw] <- dk + (max.col(cb >= runif(length(iw)), ties.method="first") - 1L)
   }
   ow <- which(!inwin)
   if (length(ow)) {                                  # outside-window: base, reject windows
     dr <- sample.int(maxday, length(ow), replace=TRUE, prob=base)
     repeat { b <- rep(FALSE,length(ow))
     for (col in 1:3){ dk <- Dc[ow,col]
     b <- b | (dr>=dk & dr<=pmin(dk+RL-1L,maxday) & dk<=maxday) }
     if (!any(b)) break
     dr[b] <- sample.int(maxday, sum(b), replace=TRUE, prob=base) }
     death[ow] <- dr
   }
   ## counts use ORIGINAL dose dates (before NA-coding)
   n_before <- sum(death <  Dc[,1])                   # SIDS before dose 1 (kept, all doses NA)
   n_after  <- sum(death >= Dc[,1])                   # vaccinated cases (>=1 dose received)
   ## NA-code: keep dose k only if received before death AND within observation
   keep_dose <- function(col) ifelse(Dc[,col] <= death & Dc[,col] <= maxday, Dc[,col], NA_integer_)
   cases <- data.frame(case = seq_len(nC), death = death,
                       d1 = keep_dose(1), d2 = keep_dose(2), d3 = keep_dose(3))
   list(cases = cases, n_before = n_before, n_after = n_after, n_total = nC)
 }
 
 ## format vaccinated cases for eventdepenexp (curtail un-received doses -> NA)
 build_md_data <- function(vax, OBS_END) {
   cut <- function(dk, dd) ifelse(dk <= dd & dk <= OBS_END, dk, NA_integer_)
   data.frame(case = vax$case, aevent = vax$death, sta = 1L, end = OBS_END,
              rv1 = cut(vax$d1, vax$death), rv2 = cut(vax$d2, vax$death),
              rv3 = cut(vax$d3, vax$death))
 }
 
 ## fit modified SCCS; start start followup  day 1; 
 ##   keep pre-vax SIDS to help estimate age effect
 ##   return shared risk-window IRR/CI/p (+ n_cases)
 fit_multidose <- function(vax, P) {
   if (is.null(vax) || nrow(vax) < 30) return(c(IRR=NA,lower=NA,upper=NA,p=NA,n_cases=0))
   md <- build_md_data(vax, P$OBS_END)
   md$sta=P$sta
   tryCatch({
     mod <- eventdepenexp(indiv=case, astart=sta, aend=end, aevent=aevent,
                          adrug=cbind(rv1,rv2,rv3), aedrug=cbind(rv1,rv2,rv3)+P$RISK_LEN,
                          expogrp=0, sameexpopar=TRUE, agegrp=P$age_bins,
                          dataformat="multi", data=md)
     ci <- mod$conf.int; irr <- ci[1,1]; lo <- ci[1,3]; up <- ci[1,4]
     se <- (log(up)-log(lo))/(2*qnorm(0.975))
     c(IRR=irr, lower=lo, upper=up, p=2*pnorm(-abs(log(irr)/se)), n_cases=nrow(md))
   }, error=function(e) c(IRR=NA,lower=NA,upper=NA,p=NA,n_cases=nrow(md)))
 }
 
 ### fit modified SCCS; but only start follow up at first vaccine dose, and 
 ###  eliminate any SIDS before vaccination from modeling.
 fit_multidose_d1 <- function(vax, P) {
   if (is.null(vax) || nrow(vax) < 30) return(c(IRR=NA,lower=NA,upper=NA,p=NA,n_cases=0))
   vax <- vax[!is.na(vax$d1),]
   md <- build_md_data(vax, P$OBS_END)
   md$sta=P$sta
   tryCatch({
     mod <- eventdepenexp(indiv=case, astart=rv1, aend=end, aevent=aevent,
                          adrug=cbind(rv1,rv2,rv3), aedrug=cbind(rv1,rv2,rv3)+P$RISK_LEN,
                          expogrp=0, sameexpopar=TRUE, agegrp=P$age_bins,
                          dataformat="multi", data=md)
     ci <- mod$conf.int; irr <- ci[1,1]; lo <- ci[1,3]; up <- ci[1,4]
     se <- (log(up)-log(lo))/(2*qnorm(0.975))
     c(IRR=irr, lower=lo, upper=up, p=2*pnorm(-abs(log(irr)/se)), n_cases=nrow(md))
   }, error=function(e) c(IRR=NA,lower=NA,upper=NA,p=NA,n_cases=nrow(md)))
 }
 
 
 ## ===========================================================================
 ## (A) PARAMETERS  — all simulation inputs in one place
 ## ===========================================================================
 P <- list(
   N        = 5e6,    # virtual population size per simulated dataset
   p_SIDS   = 2e-4,   # baseline first-year SIDS risk per infant (null hazard scale)
   RISK_LEN = 3L,     # risk window after EACH dose: days 0..RISK_LEN-1 (i.e. 0-2)
   OBS_END  = 365L,   # observation end (days); SCCS astart = 1 (birth)
   age_bins = c(60,75,90,105,120,145,170,200,250,300),  # age strata cut points (days)
   seed     = 20260624,
   sta      = 1L     # starting point of follow up in days.
 )
 
 ## Various settings for age_bins -- important to be fine enough for background
 ##   incidence (pmf_s) to be relatively constant within each age_bin
 ## Main setting is my own choice, none is approximately no age bias adjustmnet
 ## The others are the actual bins used in the six modified SCCS SIDS studies
 ##
 age_bins_main=c(60,75,90,105,120,145,170,200,250,300)
 age_bins_none=c(300,350)
 age_bins_Taiwan=c(31,60,90,120,150,180,210,240,270,300,330)
 age_bins_Italy=c(31,80,100,120,180)
 age_bins_Germany=c(37.4, 59.7, 76.6, 93, 105.9, 127.1, 151, 179.5, 226.6, 270)
 age_bins_UK=c(28,46,59,71,84,101.7,116.2,140,173.2,246.8)
 age_bins_NZ=c(39.6, 52, 62.9, 71, 81.2, 94, 106, 121, 150, 204)
   
 # inter-dose interval (GAP +/- GAP_SD) lives inside dose_distribution()'s gen_schedule.
 VAX_single  <- 1.0                # effect for the single illustrative dataset
 ## ===========================================================================
 ## (B) SINGLE DATASET — diagnostics + one analysis
 ## ===========================================================================
#set.seed(P$seed)
 P$sta=1L
P$age_bins=age_bins_main
one <- simulate_population(P$N, VAX_single, pmf_s, gen_schedule, P)

 cat(sprintf("SINGLE DATASET (VAX=%.2g):  total SIDS = %d  |  before dose-1 = %d  |  after (SCCS) = %d\n",
             VAX_single, one$n_total, one$n_before, one$n_after))
 
 plot_theory_diagnostics(dose, sids, VAX = VAX_single)   # expected curves (4 pages)
 plot_sim_diagnostics(one$cases)                            # observed, this dataset (4 pages)
 
 est1 <- fit_multidose(one$cases, P)
 cat(sprintf("Modified SCCS risk-window IRR = %.3f (%.3f, %.3f)  p = %.4g  [n_cases=%d]\n",
             est1["IRR"], est1["lower"], est1["upper"], est1["p"], est1["n_cases"]))
 est1_d1 <- fit_multidose_d1(one$cases, P)
 cat(sprintf("Modified SCCS sta=d1 risk-window IRR = %.3f (%.3f, %.3f)  p = %.4g  [n_cases=%d]\n",
             est1_d1["IRR"], est1_d1["lower"], est1_d1["upper"], est1_d1["p"], est1_d1["n_cases"]))

 
 ## Set up parameters for SIMULATION STUDY
 VAX_effects <- c(1.0, 1.25, 1.5, 1.75, 2.0) # effects for the simulation study
 K           <- 100                          # number of single data sets per effect
 ## ===========================================================================
 ## (C) SIMULATION STUDY — K datasets per VAX_effect, with SIDS-count summaries
 ## ===========================================================================
 simulate=function(P,pmf_s,VAX_effects=c(1.0,1.25,1.5,1.75,2.0),K=100,file_label="main")
 {
 set.seed(P$seed)
 rows <- vector("list", length(VAX_effects)*K); r <- 0L
 rows_d1 <-rows
 for (ve in VAX_effects) {
   for (k in seq_len(K)) {
     out <- simulate_population(P$N, ve, pmf_s, gen_schedule, P)
     est <- fit_multidose(out$cases, P)
     est_d1 <- fit_multidose_d1(out$cases, P)
     r <- r + 1L
     rows[[r]] <- data.frame(VAX_effect=ve, sim=k,
                             IRR=est["IRR"], lower=est["lower"], upper=est["upper"],
                             p=est["p"], n_cases=est["n_cases"],
                             n_before=out$n_before, n_after=out$n_after,
                             total_SIDS=out$n_total)
     rows_d1[[r]] <- data.frame(VAX_effect=ve, sim=k,
                             IRR=est_d1["IRR"], lower=est_d1["lower"], upper=est_d1["upper"],
                             p=est_d1["p"], n_cases=est_d1["n_cases"],
                             n_before=out$n_before, n_after=out$n_after,
                             total_SIDS=out$n_total)
     cat("Data set",k,"\n")
   }
   cat(sprintf("done VAX_effect = %.2f\n", ve))
 }
 results <- do.call(rbind, rows)
 results_d1 <- do.call(rbind, rows_d1)
 write.csv(results, paste("multidose_sccs_results_",file_label,".csv",sep=""), row.names = FALSE)
 write.csv(results, paste("multidose_sccs_d1_results_",file_label,".csv",sep=""), row.names = FALSE)
 
 ## summary: IRR behavior, error rates, AND population SIDS counts
 summ <- do.call(rbind, lapply(split(results, results$VAX_effect), function(x)
   data.frame(VAX_effect      = x$VAX_effect[1],
              mean_IRR        = mean(x$IRR, na.rm=TRUE),
              median_IRR      = median(x$IRR, na.rm=TRUE),
              pct_p_lt_.05    = mean(x$p < 0.05, na.rm=TRUE),     # type-I (VAX=1) / power
              ci_coverage     = mean(x$lower <= x$VAX_effect & x$upper >= x$VAX_effect, na.rm=TRUE),
              mean_total_SIDS = mean(x$total_SIDS),               # whole-population SIDS
              mean_SIDS_before= mean(x$n_before),                 # pre-dose-1 (excluded)
              mean_SIDS_after = mean(x$n_after))))                # post-dose-1 (in SCCS)
 print(summ, row.names = FALSE)
 write.csv(summ, paste("multidose_sccs_summary_",file_label,".csv",sep=""), row.names = FALSE)
 summ_d1<- do.call(rbind, lapply(split(results_d1, results_d1$VAX_effect), function(x)
   data.frame(VAX_effect      = x$VAX_effect[1],
              mean_IRR        = mean(x$IRR, na.rm=TRUE),
              median_IRR      = median(x$IRR, na.rm=TRUE),
              pct_p_lt_.05    = mean(x$p < 0.05, na.rm=TRUE),     # type-I (VAX=1) / power
              ci_coverage     = mean(x$lower <= x$VAX_effect & x$upper >= x$VAX_effect, na.rm=TRUE),
              mean_total_SIDS = mean(x$total_SIDS),               # whole-population SIDS
              mean_SIDS_before= mean(x$n_before),                 # pre-dose-1 (excluded)
              mean_SIDS_after = mean(x$n_after))))                # post-dose-1 (in SCCS)
 print(summ_d1, row.names = FALSE)
 write.csv(summ_d1, paste("multidose_sccs_d1_summary_",file_label,".csv",sep=""), row.names = FALSE)
 list(modSCCS=summ,modSCCS_d1=summ_d1)
 }
 
### Run simulation with main settings
P$age_bins=age_bins_main
simulate(P=P,VAX_effects=c(1.0,1.25,1.5,1.75,2.0),K=100,pmf_s=pmf_s,file_label="main")
### Run simulation naively not adjusting for age bias
P$age_bins=age_bins_none
simulate(P=P,pmf_s=pmf_s,file_label="no_age")

##. Now consider the age_bins used in the six actual modified SCCS SIDS studies
### Taiwan
###.  677 and 1847 deaths
###   3.45e-4 to 5.47e-4
###. Age bins
###.  Post vax risk 0-1, 2-7, 8-14, 15-30
### Follow up 31d. through 364 days
### age_bins=c(31,60,90,120,150,180,210,240,270,300,330)
###
###.  Italy, hexavalent
####. 604 through 2 years (most yr1)  
###  1.7e-4
###. age_bins=c(31,80,100,120,180)
###. post vax risk 0-1, 0-7, 0-14
###. 31d to 720.  includes pre-vax.
###
###. GESID (German)
###.  
###.  333 cases
###. age_bins=c(37.4, 59.7, 76.6, 93, 105.9, 127.1, 151, 179.5, 226.6, 270)
###. post vax risk 1-3d, 4-7, 8-14, 15-21, 22-28, 29-max
###
###. CESDI (UK)
###. 
###. 303 cases
###.  age_bins=c(28,46,59,71,84,101.7,116.2,140,173.2,246.8)
###. 1-3d, 4-7, 8-14, 15-21, 22-28, 29-max
###  only first and last vaccines recorded, so imputed.
###. imputation underrepresents risk and bias upwards
###
###  NZCD New Zealand
###
###. 393 cases
###. age_bins=c(39.6, 52, 62.9, 71, 81.2, 94, 106, 121, 150, 204)
###. Dosing 1-3, 4-7, 8-14, 15-21, 22-28, 29-max
###.
P$age_bins=age_bins_Taiwan
simulate(P=P,VAX_effects=c(1.0,1.25,1.5,1.75,2.0),K=100,pmf_s=pmf_s,file_label="Taiwan")
P$age_bins=age_bins_Italy
simulate(P=P,VAX_effects=c(1.0,1.25,1.5,1.75,2.0),K=100,pmf_s=pmf_s,file_label="Italy")
P$age_bins=age_bins_Germany
simulate(P=P,VAX_effects=c(1.0,1.25,1.5,1.75,2.0),K=100,pmf_s=pmf_s,file_label="Germany")
P$age_bins=age_bins_UK
simulate(P=P,VAX_effects=c(1.0,1.25,1.5,1.75,2.0),K=100,pmf_s=pmf_s,file_label="UK")
P$age_bins=age_bins_NZ
simulate(P=P,VAX_effects=c(1.0,1.25,1.5,1.75,2.0),K=100,pmf_s=pmf_s,file_label="New_Zealand")
 

### Plot age bins for each study
plot_SIDS_age(dose, sids, age_bins_main,"main",VAX = 1.0)   
plot_SIDS_age(dose, sids, age_bins_Taiwan,"Taiwan",VAX = 1.0)   
plot_SIDS_age(dose, sids, age_bins_Italy,"Italy",VAX = 1.0)   
plot_SIDS_age(dose, sids, age_bins_Germany,"Germany",VAX = 1.0)   
plot_SIDS_age(dose, sids, age_bins_UK,"United Kingdom",VAX = 1.0)   
plot_SIDS_age(dose, sids, age_bins_NZ,"New Zealand",VAX = 1.0)   

