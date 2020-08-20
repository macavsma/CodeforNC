# specify data:
data <- cbind.data.frame("study" = c("Veterans", "Truven", "PearlDiver",
                                     "Medicare","Clinformatics"),
                         "HR"    = c(0.665, 0.614, 0.738, 0.828, 0.727),
                         "lower" = c(0.625, 0.524, 0.600, 0.646, 0.572),
                         "upper" = c(0.708, 0.718, 0.908, 1.062, 0.924),
                         "logHR" = NA_real_,
                         "SE"    = NA_real_)

# derive logarithmic HRs and their standard errors:
data$logHR <- log(data$HR)
data$SE    <- (log(data$upper) - log(data$lower)) / (2*qnorm(0.975))


#################################
# run frequentist meta-analysis;
# load "metafor" R package:
require("metafor")
# perform analysis:
fma <- rma.uni(yi=data$logHR, sei=data$SE, slab=data$study)
# show results (summary):
fma

# show heterogeneity estimate and CI:
confint(fma)

# show effect estimates on log-HR and HR scales:
esti <- c("estimate"=fma$b, "lower"=fma$ci.lb, "upper"=fma$ci.ub)
print(rbind("log-HR"=esti, "HR"=exp(esti)))
# show heterogeneity estimate:
fma$tau2
sqrt(fma$tau2)


##############################
# run Bayesian meta-analyses;
# load "bayesmeta" R package:
require("bayesmeta")
# perform analyses:
bm1 <- bayesmeta(y=data$logHR, sigma=data$SE, labels=data$study,
                 mu.prior=c(0,2),
                 tau.prior=function(t){dhalfcauchy(t,scale=0.280)})
bm2 <- bayesmeta(y=data$logHR, sigma=data$SE, labels=data$study,
                 mu.prior=c(0,2),
                 tau.prior=function(t){dhalfcauchy(t,scale=0.587)})

# show results:
bm1
bm2

# show heterogeneity posteriors:
par(mfrow=c(2,1))
 plot(bm1, which=4, prior=TRUE, taulim=1)
 plot(bm2, which=4, prior=TRUE, taulim=1)
par(mfrow=c(1,1))

# show effect estimates;
# log-HRs:
rbind(bm1$summary[,"mu"], bm2$summary[,"mu"])
# HRs:
round(exp(rbind(bm1$summary[c(2,5:6),"mu"], bm2$summary[c(2,5:6),"mu"])),3)

# show heterogeneity estimates;
# tau:
rbind(bm1$summary[,"tau"], bm2$summary[,"tau"])
# tau-squared:
rbind(bm1$summary[c(2,5:6),"tau"], bm2$summary[c(2,5:6),"tau"])^2


# illustrate log-HRs:
forestplot(bm1)

# illustrate HRs:
forestplot(bm1, expo=TRUE, xlog=TRUE, digits=1)

# illustrate HRs, omit shrinkage intervals:
forestplot(bm1, expo=TRUE, xlog=TRUE, shrinkage=FALSE)

# illustrate HRs, adjust x-axis size etc.:
forestplot(bm1, expo=TRUE, xlog=TRUE,
           xticks=c(0.5,0.7,1.0),
           xlab="hazard ratio (HR)",
           txt_gp = fpTxtGp(ticks = gpar(cex=1), xlab = gpar(cex=1)))


# compute posterior tail probabilities
# (probability of non-beneficial effect):
(1 - bm1$pposterior(mu=0))
(1 - bm2$pposterior(mu=0))

# show shrinkage estimates;
# log-HRs:
bm1$theta
# HRs:
exp(bm1$theta[c(4,7:8),])
