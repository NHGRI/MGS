##################################################
# Tingfen Yan and Narisu Narisu
# Script to perform 1) betabinomial test and meta analysis as
# well as 2) generate pseudo bulk allelic counts across samples
# and then perform betabinomial test and meta analysis.
# Method 1) is used for the paper.
# 03/10/2024
##################################################

in.file <- commandArgs(T)[1]
out.file <- commandArgs(T)[2]

.libPaths("/cluster/ifs/projects/collins/users/yant/anaconda3/envs/seurat4/lib/R/library")
library(metap)
library(Rmpfr)
library(VGAM)
library(tools)
set.seed(123)

# Beta-binomial test: modified binom.test() after extensive
# discussion with Cassie Robertson. This is the function we
# used for final draft.
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
    else 
        stop("incorrect length of 'x'")
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
            if (p == 0) (x == 0) else if (p == 1) (x == n) else {
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
                    pbetabinom(y - 1, size = n, p = p, rho = rho) + (1 - pbetabinom(x - 1, size = n, p = p, rho = rho))
                } 
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

# Record those extreme cases with negative p value for the upper tail 
# This value should be the same as betabinomial p value, but
# not forcing very out of boundary p values (very rare cases)
# in 0 to 1 range. So this p value could be slightly off
# from the actual betabinom p value.
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
     PVAL <- switch(alternative, less = pbetabinom(x, size=n, p=p, rho=rho), 
        greater = 1 - pbetabinom(x, size=n, p=p, rho=rho),
        two.sided = {
            if (p == 0) (x == 0) else if (p == 1) (x == n) else {
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
# Read in the input file
in.data <- read.table(in.file,sep="\t",header=F, stringsAsFactors = F)

# A function to avoid 0 pvalue 
pvalue_one2two <- function(p_vector, log_space = FALSE) {
    filt_na <- is.na(p_vector)
    if (log_space) {
        p_vector_twotailed <- p_vector
        p_vector[filt_na] <- log(0.5)
        p_vector_twotailed <- log(2) + p_vector_twotailed
        # NOTE: Because these are not small p-values, we will not loose
        # precision when we move out of log space. This will throw
        # an Inf though in the case where we have a pvalue of 1.
        filt <- p_vector > log(0.5)
        tmp <- Rmpfr::mpfr(p_vector[filt], 64*15)
        p_vector_twotailed[filt] <- as.numeric(log(2) + log(1 - exp(tmp)))
    } else {
        p_vector_twotailed <- p_vector
        # p_vector[filt_na] <- 0.5 # NAs are only an issue in log space
        p_vector_twotailed <- 2 * p_vector_twotailed
        filt <- p_vector > 0.5
        p_vector_twotailed[filt] <- 2 * (1 - p_vector[filt])
    }
    p_vector_twotailed[filt_na] <- NA
    return(p_vector_twotailed)
}

stouffer.z.vec <- c()
stouffer.logP.metap <- c()
stouffer.logP.twoTail_sumz <- c()

# Assign empty vector for calculating pseudocount across sample for each SNP
pseudo.ref.count <- c()
pseudo.alt.count <- c()
pseudo.total.count <- c()

for (i in 1:nrow(in.data)) {
    r=unlist(strsplit(in.data[i,7],":"))
    r=unlist(strsplit(r,";"))
    endn=length(r)
    ref.count.cols <- seq(1,endn,5)
    alt.count.cols <- seq(2,endn,5)
    total.count.cols <- seq(3,endn,5)
    p.val.cols <- seq(4,endn,5)
  
    ref.count.mat <- as.numerbetabinom.test_2ic(r[ref.count.cols])
    alt.count.mat <- as.numeric(r[alt.count.cols])
    total.count.mat <- as.numeric(r[total.count.cols])
    p.val.mat <- as.numeric(r[p.val.cols])

    # Calculating pseudocount across sample for each SNP
    pseudo.ref.count[i] <- sum(ref.count.mat)
    pseudo.alt.count[i] <- sum(alt.count.mat)
    pseudo.total.count[i] <- sum(total.count.mat)
    
    # Convert two-tailed p-values to one-tailed p-values.
    # meta analysis
    one.tail.p.mat <- matrix(0,nrow=1,ncol=length(p.val.mat))
    # This needs to be done for each individual p-value (for every sample and every locus).
    for (j in 1:length(p.val.mat)) {
        ref.count <- as.numeric(ref.count.mat[j])
        alt.count <- as.numeric(alt.count.mat[j])
        # precalculated beta binomial p value
        p.val <- as.numeric(p.val.mat[j])
    
        if (ref.count < alt.count) {
            one.tail.p.mat[1,j] <- p.val/2
        }
        # I think the this should be ">=" instead of ">"
        # because beta-binomial distribution can have 0 value.
        # so 5 successes out of 10 draws has 0-5 on one end (length 6)
        # and 6-11 on the other (length 5).
        # I had some issues with numbers very slightly above 1,
        # so had to add in the min statement.
    
        if (ref.count >= alt.count) {
            one.tail.p.mat[1,j] <- min(1-(p.val/2),1)
        }
    }
    # Fix issue when we have 0 or 1 pvalues that would otherwise be dropped from
    # the meta-analysis
    one.tail.p.mat[one.tail.p.mat == 1.0] <- 1.0 - (1e-7)
    one.tail.p.mat[one.tail.p.mat == 0.0] <- .Machine$double.xmin

    ##################
    # Using sumz
    ##################
    # endn > 5 means there are at least two replicates in the file
    # endn <=5 means there is only one data, no need to do meta analysis, simply take the pvalue
    if (endn > 5) {
        z.matS <- sumz(p=one.tail.p.mat[1,], weights=total.count.mat, log.p = TRUE, log.input = FALSE)
        sumz=z.matS$z[[1]][1]
        stouffer.z.vec[i] <- sumz
        logp=z.matS$p[[1]][1]
        stouffer.logP.metap[i] <- logp
    } else {
        stouffer.logP.metap[i] <- log(one.tail.p.mat)
    }
}

### Adjust Stouffer p-values for multiple testing using the BH procedure
#########################################################################
stouffer.logP.twoTail_sumz = pvalue_one2two(stouffer.logP.metap, log_space = TRUE)

stouffer.q.twoTail.sumz <- p.adjust(exp(stouffer.logP.twoTail_sumz), method = "BH")
# Organize the output data in a dataframe
out.data0 <- in.data[,1:6]
out.data0["meta_betabinomial_p"] <- exp(stouffer.logP.twoTail_sumz)
out.data0["meta_betabinomial_q"] <- stouffer.q.twoTail.sumz

# Calculating p value for binom.test
out.data0["pseudo.ref.count"] <- pseudo.ref.count
out.data0["pseudo.total.count"] <- pseudo.total.count

out.data <- out.data0
successes <- out.data$pseudo.ref.count
num.draws <- out.data$pseudo.total.count

fit.bb <- vglm(cbind(successes,num.draws-successes)~1,betabinomial)
beta.binom.rr.est <- signif(logitlink(attr(fit.bb,"coef")[1],inverse=T),digits=4)
beta.binom.rho.est <- signif(logitlink(attr(fit.bb,"coef")[2],inverse=T),digits=4)
beta.binom.prob <- round(signif(beta.binom.rr.est,digits=4),4)
beta.binom.prob
beta.binom.rho <- round(signif(beta.binom.rho.est,digits=4),4)
beta.binom.rho

beta.binom.p.vals.flag<-rep(0,length(successes))
beta.binom.p.vals <- rep(0,length(successes))
RR.ref <- rep(0,length(successes))
# for reference allele and overall test
for (j in 1:length(successes)) {
    beta.binom.p.vals[j] <- betabinom.test_2(successes[j],num.draws[j],p=beta.binom.prob,rho=beta.binom.rho,alternative="two.sided")$p.value
    beta.binom.p.vals.flag[j] <- betabinom.test_flag(successes[j],num.draws[j],p=beta.binom.prob,rho=beta.binom.rho,alternative="two.sided")$p.value
    ref.ratio <- successes[j]/num.draws[j]
    alt.ratio <- (num.draws[j] - successes[j])/num.draws[j]
    # Effect size(https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0762-6)
    RR.ref[j] <- 0.5-ref.ratio
}
# binom.q.vals <- p.adjust(binom.p.vals,method="BH")
beta.binom.q.vals <- p.adjust(beta.binom.p.vals,method="BH")
out.data["pseudo_betabinomial_p"] <- beta.binom.p.vals
out.data["pseudo_betabinomial_q"] <- beta.binom.q.vals
# Very very rare cases, pseudo_betabinomial_p_flag could be out
# of range of 0-1. For these cases, we forced the
# beta.binom.p.vals in the pseudo betabionomial p values to
# be in this range. When pseudo_betabinomial_p and
# pseudo_betabinomial_p_flag are not equal, the test might
# be less reliable. 
out.data["pseudo_betabinomial_p_flag"] <- beta.binom.p.vals.flag
out.data["ref.effect"] <- RR.ref

out.data1 <- cbind(out.data,in.data[in.data$V4 %in% out.data$V4,"V7"])
out.names <- c("#chr","start","stop","var_id","ref_base","alt_base",
    "meta_betabinomial_p","meta_betabinomial_q","pseudo.ref.count",
    "pseudo.total.count","pseudo_betabinomial_p","pseudo_betabinomial_q","pseudo_betabinomial_p_flag","ref.effect","variantInf")

names(out.data1) <- out.names

write.table(out.data1,out.file,quote=F,sep="\t",row.names=F)
