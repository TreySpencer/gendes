
Power <- function(.dose0 = 6:10, .delta = 1, .rho = 1.2, .ed50 = 7.5, .slope = 25, .Npergrp = 5, .alpha = 0.05) {
    probs <- response.targets(slope = .slope, rho = .rho, ed50 = .ed50, dose0 = .dose0, delta = .delta)
    p0 <- probs$p0
    p1 <- probs$p1
    .dose1 <- probs$.dose1
    # Computing the degrees of freedom
    df <- 2 * nrow(probs) - 3
    # Computing critical t under H0
    t.alpha <- qt(1 - .alpha, df)
    # Computes the probits for the given target responses
    z0 <- qnorm(p0)
    z1 <- qnorm(p1)
    # Computes probit weights
    logit0 <- p0 * (1 - p0)
    logit1 <- p1 * (1 - p1)
    phi.z0 <- dnorm(z0)
    phi.z1 <- dnorm(z1)
    probit.wt0 <- phi.z0^2/logit0
    probit.wt1 <- phi.z1^2/logit1
    sum.wt0 <- sum(probit.wt0)
    sum.wt1 <- sum(probit.wt1)
    #---- Computing xbar & ybar & log(dose) sq err for User-specified CTRL
    wz0 <- probit.wt0 * z0
    wlogx0 <- probit.wt0 * log10(.dose0)
    xbar0 <- sum(wlogx0)/sum(probit.wt0)
    ybar0 <- sum(wz0)/sum(probit.wt0)
    logx0.sqerr <- (log10(.dose0) - xbar0)^2
    #---- Computing xbar & ybar & log(dose) sq err for User-specified TRT
    wz1 <- probit.wt1 * z1
    wlogx1 <- probit.wt1 * log10(.dose1)
    xbar1 <- sum(wlogx1)/sum(probit.wt1)
    ybar1 <- sum(wz1)/sum(probit.wt1)
    logx1.sqerr <- (log10(.dose1) - xbar1)^2
    #---- From Eqn 8 in Kodell et al. (2010), computing SUM( w_i*n_i*(x_i - x.bar)^2 ) for CTRL & TRT
    sumsq.logx0 <- sum(.Npergrp * probit.wt0 * logx0.sqerr)
    sumsq.logx1 <- sum(.Npergrp * probit.wt1 * logx1.sqerr)
    #---- Variance of theta-hat - Eqn 8 in Kodell et al. (2010)
    # var.rhohat <- (tmp.slope^-2)*( (1/sum( probit.wt0*.Npergrp )) + (1/sum( probit.wt1*.Npergrp )) + (((ybar1
    # - ybar0)/tmp.slope)^2)*(1/(sumsq.logx0 + sumsq.logx1)) )
    var.rhohat <- (.slope^-2) * ((1/sum(probit.wt0 * .Npergrp)) + (1/sum(probit.wt1 * .Npergrp)) + (((ybar1 - 
        ybar0)/.slope)^2) * (1/(sumsq.logx0 + sumsq.logx1)))
    #---- t-statistic evaluated at alternative hypothesized value of theta
    t.beta <- log10(.rho)/sqrt(var.rhohat)
    .power <- pt(t.beta - t.alpha, df, lower.tail = TRUE)
    .delta <- ifelse(is.null(.delta), NA, .delta)
    list(cbind(.delta, .rho, .ed50, .slope, .Npergrp, .alpha, .power), probs)
}
