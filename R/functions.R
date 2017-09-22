#' Sample size for same target design.
#'
#' \code{staggerN} returns the sample size required to acheive a specified power.
#'
#' This function calculates the sample size per dose group required to acheive a specified power.
#' The calculations are based on a same-target design; that is, the doses used for the control and
#' intervention groups will produce the same target lethality probabilities. It is assumed a probit
#' model will be used to analyze the data.  Total sample size is $2*.Npergrp*.ngrps$.
#'
#' @param .alpha Numeric scalar or vector.  The Type I error rate (default is 0.05).
#' @param .power Numeric scalar or vector.  The target power of the test. The default is 0.05.
#' @param .slope Numeric scalar or vector.  The slope(s) of the $log_10$ doses.  A probit models
#' is assumed.
#' @param .rho Numeric scalar or vector.  The dose reduction factor.
#' @param .ngrps Numeric scalar or vector.  The number of dosing groups.  This value is ignored if .probs
#' is not NULL.
#' @param .minp Numeric scalar.  The minimum lethality probability.
#' 
staggerN <- function(.alpha = 0.05, .power = 0.8, .slope = 2, .rho = 1.1, .ngrps = 5, .minp = 0.05,
                     .probs = NULL) {
    if (is.null(.probs))
        probs <- seq(.minp, 1 - .minp, length.out = .ngrps) 
    else {
        probs <- .probs
        .ngrps <- length(probs)
        .minp <- NA
    }

    logit <- probs * (1 - probs)
    density.sq <- dnorm(qnorm(probs))^2
    sum.wt <- sum(density.sq/logit)
    df <- 2 * .ngrps - 3

    x <- expand.grid(.alpha, .power, .slope, .rho)
    .alpha <- x[, 1]
    .power <- x[, 2]
    .slope <- x[, 3]
    .rho <- x[, 4]

    t.alpha <- qt(1 - .alpha, df)
    t.power <- qt(.power, df)
    log.rho <- log10(.rho)
    num <- 2 * (t.alpha + t.power)^2
    denom <- sum.wt * (.slope * log.rho)^2
    n.per.grp <- ceiling(num/denom)
    n.total <- 2 * .ngrps * n.per.grp
    list(results = cbind(.alpha, .power, .slope, .rho, .ngrps, n.per.grp, n.total), probs = probs)
}

stagger.Power <- function(.Npergrp = 1033, .alpha = 0.05, .slope = 2, .rho = 1.1, .ngrps = 5, .minp = 0.05,
    .probs = NULL) {
    if (is.null(.probs))
        probs <- seq(.minp, 1 - .minp, length.out = .ngrps) else {
        probs <- .probs
        .ngrps <- length(probs)
        .minp <- NA
    }

    logit <- probs * (1 - probs)
    density.sq <- dnorm(qnorm(probs))^2
    sum.wt <- sum(density.sq/logit)
    df <- 2 * .ngrps - 3

    x <- expand.grid(.Npergrp, .alpha, .slope, .rho)
    .Npergrp <- x[, 1]
    .alpha <- x[, 2]
    .slope <- x[, 3]
    .rho <- x[, 4]

    t.alpha <- qt(1 - .alpha, df)
    log.rho <- log10(.rho)
    num <- .Npergrp * sum.wt * (.slope * log.rho)^2
    t.power <- sqrt(num/2) - t.alpha
    power <- pt(t.power, df)
    list(results = cbind(.Npergrp, .alpha, .slope, .rho, .ngrps, t.power, df, power), probs = probs)
}

#'  I modified the *Stagger.Power* function to accept doses instead of target lethality probabilities.
#'  This function assumes a *Probit* model.
#'
#'  **Function Inputs**:
#'
#'  *  $.Npergrp$ is the number of animals per dose group,
#'  *  $.alpha$ is the Type 1 error rate (one-sided),
#'  *  $.slope$ is the slope of the probit model,
#'  *  $.rho$ is the expected dose reduction factor,
#'  *  $.ed50$ is the expected dose resulting in 50% lethality in control animals, and
#'  *  $.dose0$ is the target doses for the controls.
#'
#'  **Function Outputs**:
#'
#'  In addition to the user-supplied parameters listed above, the function outputs the following information:
#'
#'  *  $t.power$ is the $t$-statistic corresponding to the power estimate,
#'  *  $df$ is the degrees-of-freedom, and
#'  *  $power$ is the power estimate.
#'
#'  Note that the total number of animals required is $2*.Npergrp*.ngrps$.
#'
#'  Here is the code:
#'

dose.ST.Power <- function(.Npergrp = 5, .alpha = 0.05, .slope = 25, .rho = 1.2, .ed50 = 7.5, .dose0 = 6:10) {
    .ngrps <- length(.dose0)
    df <- 2 * .ngrps - 3

    x <- expand.grid(.Npergrp, .alpha, .slope, .rho, .ed50)
    .Npergrp <- x[, 1]
    .alpha <- x[, 2]
    .slope <- x[, 3]
    .rho <- x[, 4]
    .ed50 <- x[, 5]

    y <- cbind(.slope, .ed50)
    sum.wt <- apply(y, 1, sumW.fn, .dose0)

    t.alpha <- qt(1 - .alpha, df)
    log.rho <- log10(.rho)
    num <- .Npergrp * sum.wt * (.slope * log.rho)^2
    t.power <- sqrt(num/2) - t.alpha
    power <- pt(t.power, df)
    cbind(.slope, .rho, .ed50, .Npergrp, .ngrps, .alpha, t.power, df, power)
}

stagger.Rho <- function(.Npergrp = 1033, .alpha = 0.05, .power = 0.8, .slope = 2, .ngrps = 5, .minp = 0.05,
    .probs = NULL) {

    if (is.null(.probs))
        probs <- seq(.minp, 1 - .minp, length.out = .ngrps) else {
        probs <- .probs
        .ngrps <- length(probs)
        .minp <- NA
    }

    logit <- probs * (1 - probs)
    density.sq <- dnorm(qnorm(probs))^2
    sum.wt <- sum(density.sq/logit)
    df <- 2 * .ngrps - 3

    x <- expand.grid(.Npergrp, .alpha, .power, .slope)
    .Npergrp <- x[, 1]
    .alpha <- x[, 2]
    .power <- x[, 3]
    .slope <- x[, 4]

    t.alpha <- qt(1 - .alpha, df)
    t.power <- qt(.power, df)
    num <- 2 * (t.alpha + t.power)^2
    denom <- .Npergrp * sum.wt * (.slope)^2
    log.rho <- sqrt(num/denom)
    rho <- 10^(log.rho)
    list(results = cbind(.Npergrp, .alpha, .power, .slope, .ngrps, rho), probs = probs)
}

stagger.Doses <- function(.slope = 2, .rho = 1.1, .ed50 = 10, .ngrps = 5, .minp = 0.05, .probs = NULL) {
    if (is.null(.probs))
        probs <- seq(.minp, 1 - .minp, length.out = .ngrps) else {
        probs <- .probs
        .ngrps <- length(probs)
        .minp <- NA
    }
    C.log.dose <- qnorm(probs)/.slope + log10(.ed50)
    C.dose <- 10^C.log.dose
    T.dose <- .rho * C.dose
    cbind(probs = probs, control = C.dose, treat = T.dose)
}

