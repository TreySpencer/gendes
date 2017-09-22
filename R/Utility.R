
sumW.fn <- function(x, dose) {
    z.dose <- x[1] * log10(dose/x[2])
    probs <- pnorm(z.dose)
    logit <- probs * (1 - probs)
    density.sq <- dnorm(z.dose)^2
    sum(density.sq/logit)
}


response.targets <- function(slope = 25, rho = 1.2, ed50 = 7.5, dose0 = 6:10, delta = 1) {
    p0 <- pnorm(((log10(dose0) - log10(ed50)) * slope))
    if (!is.null(delta)) {
        dose1 <- dose0 + delta
        p1 <- pnorm(((log10(dose1) - log10(ed50 * rho)) * slope))
    } else {
        p1 <- p0
        dose1 <- rho * dose0
    }
    data.frame(.dose0 = dose0, p0 = p0, .dose1 = dose1, p1 = p1)
}
