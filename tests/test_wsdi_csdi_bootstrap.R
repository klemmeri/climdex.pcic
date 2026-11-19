library(climdex.pcic)
library(RUnit)

## WSDI and CSDI use the Zhang et al. (2005) in-base bootstrap thresholds for
## years inside the base period, as the TX90p-type indices do.
climdex.pcic.test.wsdi.csdi.bootstrap <- function() {
  tmax.dates <- as.PCICt(do.call(paste, ec.1018935.tmax[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  tmin.dates <- as.PCICt(do.call(paste, ec.1018935.tmin[,c("year", "jday")]), format="%Y %j", cal="gregorian")
  ci <- climdexInput.raw(ec.1018935.tmax$MAX_TEMP, ec.1018935.tmin$MIN_TEMP, NULL,
                         tmax.dates, tmin.dates, NULL, base.range=c(1981, 1990))

  for(index in c("wsdi", "csdi")) {
    fun <- match.fun(paste("climdex", index, sep="."))
    for(span in c(FALSE, TRUE)) {
      old <- fun(ci, spells.can.span.years=span, bootstrap=FALSE)
      new <- fun(ci, spells.can.span.years=span)
      years <- as.numeric(names(new))
      in.base <- years >= 1981 & years <= 1990

      ## Years outside the base period are unchanged.
      checkEquals(old[!in.base], new[!in.base], paste(index, "changed outside the base period"))
      ## Years inside the base period use the bootstrap thresholds.
      checkTrue(any(abs(old[in.base] - new[in.base]) > 0, na.rm=TRUE), paste(index, "unchanged inside the base period"))
      checkTrue(all(new[in.base] >= 0, na.rm=TRUE))
    }
  }

  ## Check one base year (1986; not masked for missing data) directly against the definition: the mean, over the
  ## (n - 1) bootstrap threshold sets for that year, of the number of days in
  ## spells of at least 6 days.
  q90 <- ci@quantiles$tmax$inbase$q90
  yr <- as.numeric(format(ci@dates, "%Y"))
  idx <- which(yr == 1986)
  expected <- mean(sapply(1:dim(q90)[3], function(k) {
    sum(select.blocks.gt.length(ci@data$tmax[idx] > q90[ci@jdays[idx], 1986 - 1981 + 1, k], 5))
  }))
  checkEqualsNumeric(expected, climdex.wsdi(ci)["1986"])
}

## Synthetic data with warm and cold spells that cross 31 December between base
## years, so that spells.can.span.years matters for the bootstrap.
make.spanning.ci <- function(base.range=c(1981, 1990)) {
  dates <- seq(as.PCICt("1979-01-01", cal="gregorian"), as.PCICt("1994-12-31", cal="gregorian"), by=86400)
  doy <- as.numeric(format(dates, "%j"))
  set.seed(42)
  clim <- 15 - 10 * cos(2 * pi * (doy - 15) / 365)
  tx <- clim + rnorm(length(dates), sd=3)
  ## Warm spells of 10 days (27 Dec - 5 Jan) and cold spells around the turn of
  ## other years, set close to the 90th (10th) percentile so that whether a day
  ## exceeds the threshold depends on the bootstrap sample.
  for(y in c(1982, 1984, 1986, 1988)) {
    i <- which(dates >= as.PCICt(paste0(y, "-12-27"), cal="gregorian") & dates <= as.PCICt(paste0(y + 1, "-01-05"), cal="gregorian"))
    tx[i] <- clim[i] + 3.85
  }
  tn <- tx - 10
  for(y in c(1983, 1985, 1987)) {
    i <- which(dates >= as.PCICt(paste0(y, "-12-28"), cal="gregorian") & dates <= as.PCICt(paste0(y + 1, "-01-06"), cal="gregorian"))
    tn[i] <- clim[i] - 10 - 3.85
  }
  climdexInput.raw(tx, tn, NULL, dates, dates, NULL, base.range=base.range)
}

## Reviewer's point: the result must not depend on the order of the bootstrap
## samples, including for spells that cross a boundary between two base years.
climdex.pcic.test.wsdi.csdi.bootstrap.sample.order <- function() {
  ci <- make.spanning.ci()
  for(v in c("tmax", "tmin")) {
    q <- ci@quantiles[[v]]$inbase[[if(v == "tmax") "q90" else "q10"]]
    op <- if(v == "tmax") ">" else "<"
    out <- ci@quantiles[[v]]$outbase[[if(v == "tmax") "q90" else "q10"]]
    set.seed(7)
    qp <- q
    for(y in 1:dim(q)[2]) qp[, y, ] <- q[, y, sample(dim(q)[3])]
    for(span in c(FALSE, TRUE)) {
      a <- threshold.exceedance.duration.index(ci@data[[v]], ci@date.factors$annual, ci@jdays, out, op, spells.can.span.years=span,
                                               max.missing.days=ci@max.missing.days['annual'], dates=ci@dates, base.thresholds=q, base.range=ci@base.range)
      b <- threshold.exceedance.duration.index(ci@data[[v]], ci@date.factors$annual, ci@jdays, out, op, spells.can.span.years=span,
                                               max.missing.days=ci@max.missing.days['annual'], dates=ci@dates, base.thresholds=qp, base.range=ci@base.range)
      checkEqualsNumeric(a, b, paste(v, "depends on the order of the bootstrap samples; span =", span))
    }
  }

  ## Engineered thresholds: in sample 1 only, the last 5 days of 1984 and the
  ## first 5 days of 1985 are exceeded, so a 10-day spell exists only if the
  ## same threshold set is used on both sides of 31 December. Moving 1985's
  ## sample 1 to another position must not change the result.
  dat <- rep(10, length(ci@dates))
  yr <- as.numeric(format(ci@dates, "%Y"))
  q <- array(20, dim(ci@quantiles$tmax$inbase$q90))
  q[361:365, 1984 - 1981 + 1, 1] <- 0
  q[1:5, 1984 - 1981 + 1, 1] <- 0
  q[361:365, 1985 - 1981 + 1, 1] <- 0
  q[1:5, 1985 - 1981 + 1, 1] <- 0
  qp <- q
  qp[, 1985 - 1981 + 1, ] <- q[, 1985 - 1981 + 1, c(2:dim(q)[3], 1)]
  run <- function(qq) threshold.exceedance.duration.index(dat, ci@date.factors$annual, ci@jdays, rep(20, 365), ">", spells.can.span.years=TRUE,
                                                          max.missing.days=ci@max.missing.days['annual'], dates=ci@dates, base.thresholds=qq, base.range=ci@base.range)
  checkEqualsNumeric(run(q)[c("1984", "1985")], run(qp)[c("1984", "1985")])
  ## Each year's own sample 1 gives a 10-day spell over its year boundaries:
  ## 5 days in the year itself at each end, averaged over the 9 samples.
  checkEqualsNumeric(c(10, 10) / 9, run(q)[c("1984", "1985")])
}

## A spell crossing from one base year into the next, with spells.can.span.years
## = TRUE, is evaluated with one threshold set per (year, replacement) pair.
## Check one year directly against that definition.
climdex.pcic.test.wsdi.bootstrap.spell.across.base.years <- function() {
  ci <- make.spanning.ci()
  q <- ci@quantiles$tmax$inbase$q90
  yr <- as.numeric(format(ci@dates, "%Y"))
  y <- 1984
  win <- which(yr >= y - 1 & yr <= y + 1)
  expected <- mean(sapply(1:dim(q)[3], function(k) {
    s <- select.blocks.gt.length(ci@data$tmax[win] > q[ci@jdays[win], y - 1981 + 1, k], 5)
    sum(s[yr[win] == y])
  }))
  checkEqualsNumeric(expected, climdex.wsdi(ci, spells.can.span.years=TRUE)[as.character(y)])
  ## The spell over 31 Dec 1984 / 1 Jan 1985 adds days to both years.
  checkTrue(climdex.wsdi(ci, spells.can.span.years=TRUE)["1985"] >= climdex.wsdi(ci, spells.can.span.years=FALSE)["1985"])
}

## A baseline too short for any bootstrap sample (one year) inside a longer
## record must not fail; the out-of-base thresholds are used.
climdex.pcic.test.wsdi.csdi.bootstrap.one.year.base <- function() {
  ci <- make.spanning.ci(base.range=c(1985, 1985))
  for(index in c("wsdi", "csdi")) {
    fun <- match.fun(paste("climdex", index, sep="."))
    for(span in c(FALSE, TRUE))
      checkEquals(fun(ci, spells.can.span.years=span, bootstrap=FALSE), fun(ci, spells.can.span.years=span))
  }
}

## Bootstrap samples with missing thresholds are not used for that period, and
## a period with no usable sample is NA. Mismatched dimensions are an error.
climdex.pcic.test.wsdi.bootstrap.missing.thresholds <- function() {
  ci <- make.spanning.ci()
  q <- ci@quantiles$tmax$inbase$q90
  out <- ci@quantiles$tmax$outbase$q90
  run <- function(qq) threshold.exceedance.duration.index(ci@data$tmax, ci@date.factors$annual, ci@jdays, out, ">", spells.can.span.years=FALSE,
                                                          max.missing.days=ci@max.missing.days['annual'], dates=ci@dates, base.thresholds=qq, base.range=ci@base.range)
  yr <- as.numeric(format(ci@dates, "%Y"))
  idx <- which(yr == 1986)
  per.sample <- sapply(1:dim(q)[3], function(k) sum(select.blocks.gt.length(ci@data$tmax[idx] > q[ci@jdays[idx], 6, k], 5)))
  q1 <- q
  q1[360:365, 6, 1] <- NA
  checkEqualsNumeric(mean(per.sample[-1]), run(q1)["1986"])
  q2 <- q
  q2[100, 6, ] <- NA
  checkTrue(is.na(run(q2)["1986"]))
  checkException(run(q[, 1:5, ]), silent=TRUE)
}
