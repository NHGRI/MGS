#########################################################
# Tingfen Yan and Narisu Narisu
# Script to perform binomial, beta binomial tests
# 09/30/2023
#########################################################

library(VGAM)
library(tools)

options(scipen = 999)

# File containing list of ASEReadCounter output files
in.files        <- commandArgs(T)[1]
# Outfile
out.files       <- commandArgs(T)[2]
# Minimum number of total reads per locus
min.total       <- as.numeric(commandArgs(T)[3])

# binom.test()/Cassie way - we will use this way
betabinom.test_2 <- function(x,n,p,rho, alternative = c("two.sided", "less", "greater"), conf.level = 0.95) {
    DNAME <- deparse(substitute(x))
    xr <- round(x)
    if (any(is.na(x) | (x < 0)) || max(abs(x - xr)) > 1e-07) 
        stop("'x' must be nonnegative and integer")
    x <- xr
    if (length(x) == 2L) {
        n <- sum(x)
        x <- x[1L]
    }
    else if (length(x) == 1L) {
        nr <- round(n)
        if ((length(n) > 1L) || is.na(n) || (n < 1) || abs(n - nr) > 1e-07 || (x > nr)) 
            stop("'n' must be a positive integer >= 'x'")
        DNAME <- paste(DNAME, "and", deparse(substitute(n)))
        n <- nr
    }
    else stop("incorrect length of 'x'")
    if (!missing(p) && (length(p) > 1L || is.na(p) || p < 0 || p > 1)) 
        stop("'p' must be a single number between 0 and 1")
    alternative <- match.arg(alternative)
    if (!((length(conf.level) == 1L) && is.finite(conf.level) && 
        (conf.level > 0) && (conf.level < 1))) 
        stop("'conf.level' must be a single number between 0 and 1")
    PVAL <- switch(alternative, less = pbetabinom(x, size=n, p=p, rho=rho), 
        greater = { 
            if (pbetabinom(x, size=n, p=p, rho=rho) > 1 ) {
                0 
            } else {           
                1 - pbetabinom(x, size=n, p=p, rho=rho) 
            }
        }, 
        two.sided = {
        if (p == 0) 
            (x == 0) 
        else if (p == 1) 
            (x == n) 
        else {
            relErr <- 1 + 1e-07
            d <- dbetabinom(x, size = n, p = p, rho = rho)
            m <- n * p
            if (x == m) 1 else if (x < m) {
                i <- seq.int(from = ceiling(m), to = n)
                y <- sum(dbetabinom(i, size = n, p = p, rho = rho) <= d * relErr)
                if (pbetabinom(n - y, size = n, p = p, rho = rho) > 1 ) {
                pbetabinom(x, size = n, p = p, rho = rho)
                } else {
                pbetabinom(x, size = n, p = p, rho = rho) + (1 - pbetabinom(n - y, size = n, p = p, rho= rho))}
            } else {
                i <- seq.int(from = 0, to = floor(m))
                y <- sum(dbetabinom(i, size = n, p = p, rho = rho) <= d * relErr) 
                if (pbetabinom(x - 1, size = n, p = p, rho = rho) > 1 ) {
                  pbetabinom(y - 1, size = n, p = p, rho = rho)
                } else {
                pbetabinom(y - 1, size = n, p = p, rho = rho) + (1 - pbetabinom(x - 1, size = n, p = p, rho = rho))} 
            }
        }   
    })
    p.L <- function(x, alpha) {
        if (x == 0) 
            0
        else 
            qbeta(alpha, x, n - x + 1)
    }
    p.U <- function(x, alpha) {
        if (x == n) 
            1
        else 
            qbeta(1 - alpha, x + 1, n - x)
    }
    CINT <- switch(alternative, less = c(0, p.U(x, 1 - conf.level)), 
        greater = c(p.L(x, 1 - conf.level), 1), two.sided = {
            alpha <- (1 - conf.level)/2
            c(p.L(x, alpha), p.U(x, alpha))
        })
    attr(CINT, "conf.level") <- conf.level
    ESTIMATE <- x/n
    names(x) <- "number of successes"
    names(n) <- "number of trials"
    names(ESTIMATE) <- names(p) <- "probability of success"
    structure(list(statistic = x, parameter = n, p.value = PVAL, 
        conf.int = CINT, estimate = ESTIMATE, null.value = p, 
        alternative = alternative, method = "Exact betabinomial test", 
        data.name = DNAME), class = "htest")
}
betabinom.test_flag <- function(x,n,p,rho, alternative = c("two.sided", "less", "greater"), conf.level = 0.95) {
    DNAME <- deparse(substitute(x))
    xr <- round(x)
    if (any(is.na(x) | (x < 0)) || max(abs(x - xr)) > 1e-07)
        stop("'x' must be nonnegative and integer")
    x <- xr
    if (length(x) == 2L) {
        n <- sum(x)
        x <- x[1L]
    }
    else if (length(x) == 1L) {
        nr <- round(n)
        if ((length(n) > 1L) || is.na(n) || (n < 1) || abs(n - nr) > 1e-07 || (x > nr))
            stop("'n' must be a positive integer >= 'x'")
        DNAME <- paste(DNAME, "and", deparse(substitute(n)))
        n <- nr
    }
    else 
        stop("incorrect length of 'x'")
    if (!missing(p) && (length(p) > 1L || is.na(p) || p < 0 || p > 1))
        stop("'p' must be a single number between 0 and 1")
    alternative <- match.arg(alternative)
    if (!((length(conf.level) == 1L) && is.finite(conf.level) &&
        (conf.level > 0) && (conf.level < 1)))
        stop("'conf.level' must be a single number between 0 and 1")
     PVAL <- switch(alternative, less = pbetabinom(x, size=n, p=p, rho=rho), greater = 1 - pbetabinom(x, size=n, p=p, rho=rho),
        two.sided = {
        if (p == 0) 
            (x == 0) 
        else if (p == 1) 
            (x == n) 
        else {
            relErr <- 1 + 1e-07
            d <- dbetabinom(x, size = n, p = p, rho = rho)
            m <- n * p
            if (x == m) 1 else if (x < m) {
                i <- seq.int(from = ceiling(m), to = n)
                y <- sum(dbetabinom(i, size = n, p = p, rho = rho) <= d * relErr)
                pbetabinom(x, size = n, p = p, rho = rho) + (1 - pbetabinom(n - y, size = n, p = p, rho = rho))
            } else {
                i <- seq.int(from = 0, to = floor(m))
                y <- sum(dbetabinom(i, size = n, p = p, rho = rho) <= d * relErr)
                pbetabinom(y - 1, size = n, p = p, rho = rho) + (1 - pbetabinom(x - 1, size = n, p = p, rho = rho))
            }
        }
    })
    p.L <- function(x, alpha) {
        if (x == 0)
            0
        else 
            qbeta(alpha, x, n - x + 1)
    }
    p.U <- function(x, alpha) {
        if (x == n) 
            1
        else 
            qbeta(1 - alpha, x + 1, n - x)
    }
    CINT <- switch(alternative, less = c(0, p.U(x, 1 - conf.level)),
        greater = c(p.L(x, 1 - conf.level), 1), two.sided = {
            alpha <- (1 - conf.level)/2
            c(p.L(x, alpha), p.U(x, alpha))
        })
    attr(CINT, "conf.level") <- conf.level
    ESTIMATE <- x/n
    names(x) <- "number of successes"
    names(n) <- "number of trials"
    names(ESTIMATE) <- names(p) <- "probability of success"
    structure(list(statistic = x, parameter = n, p.value = PVAL,
        conf.int = CINT, estimate = ESTIMATE, null.value = p,
        alternative = alternative, method = "Exact betabinomial test",
        data.name = DNAME), class = "htest")
} 
data <- read.table(in.files,header=T,comment.char="")

# Modify only require totalCount >= min.total, not individual alleles
data.subset <- data[data$totalCount >= min.total ,]

num.draws <- data.subset[,8]
successes <- data.subset[,6]

# Estimate binomial probability of success (ref allele ratio)
# using maximum likelihood
fit.b <- vglm(cbind(successes,num.draws-successes)~1,binomialff)
binom.rr.est <- signif(logitlink(attr(fit.b,"coef")[1],inverse=T),digits=4)
binom.prob <- round(signif(binom.rr.est,digits=4),4)

fit.bb <- vglm(cbind(successes,num.draws-successes)~1,betabinomial)
beta.binom.rr.est <- signif(logitlink(attr(fit.bb,"coef")[1],inverse=T),digits=4)
beta.binom.rho.est <- signif(logitlink(attr(fit.bb,"coef")[2],inverse=T),digits=4)
beta.binom.prob <- round(signif(beta.binom.rr.est,digits=4),4)
beta.binom.rho <- round(signif(beta.binom.rho.est,digits=4),4)

binom.p.vals <- rep(0,length(successes))
beta.binom_2.p.vals <- rep(0,length(successes))
beta.binom_2.p.vals.flag <- rep(0,length(successes))
RR.ref <- rep(0,length(successes))
RR.alt <- rep(0,length(successes))

for (j in 1:length(successes)) {
    binom.p.vals[j] <- binom.test(successes[j],num.draws[j],p=binom.prob,alternative="two.sided")$p.value
    # binom.test()/Cassie way 
    beta.binom_2.p.vals[j] <- betabinom.test_2(successes[j],num.draws[j],p=beta.binom.prob,rho=beta.binom.rho,alternative="two.sided")$p.value
    # flag those extreme cases that result in a negative p-value
    beta.binom_2.p.vals.flag[j] <- betabinom.test_flag(successes[j],num.draws[j],p=beta.binom.prob,rho=beta.binom.rho,alternative="two.sided")$p.value
    
    ref.ratio <- successes[j]/num.draws[j]
    alt.ratio <- (num.draws[j] - successes[j])/num.draws[j]
    # effect size(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6)
    RR.ref[j] <- 0.5-ref.ratio
}

binom.q.vals <- p.adjust(binom.p.vals,method="BH")
beta.binom_2.q.vals <- p.adjust(beta.binom_2.p.vals,method="BH")
# Make the first three columns of out file in bed format.
# This will make intersections easy
data.subset[,1] = gsub("chr", "", data.subset[,1])
data.out <- data.frame("#chr"=data.subset[,1],check.names=F)
data.out["start"] <- data.subset[,2]-1
data.out["stop"] <- data.subset[,2]
data.out["variant_ID"] <- data.subset[,3]
data.out["ref_base"] <- data.subset[,4]
data.out["alt_base"] <- data.subset[,5]
data.out["ref_count"] <- data.subset[,6]
data.out["alt_count"] <- data.subset[,7]
data.out["total_count"] <- data.subset[,8]
# betabinomial test
data.out["beta_binom_p_val"] <- beta.binom_2.p.vals
data.out["beta_binom_q_val"] <- beta.binom_2.q.vals
data.out["beta.binom.flag"] <- beta.binom_2.p.vals.flag
data.out["ref.effectSize"] <- RR.ref
# beta-binomial test: binom.test()/Cassie way 
# binomial test: binom.test()/Cassie way when rho == 0
data.out["beta.binom.prob"] <- beta.binom.prob
data.out["beta.binom.rho"] <- beta.binom.rho
# binomial test
data.out["binom_p_val"] <- binom.p.vals
data.out["binom_q_val"] <- binom.q.vals
write.table(data.out, file=out.files, sep="\t",quote=F,row.names=F,col.names=T)
