

#TIMESTAMP SESSION
cat(paste("Session run on:",date(),"\n\n"))
#REMOVE ALL OBJECTS IN R ENVIRONMENT
rm(list=ls())
#SET CURRENT WORKING DIRECTORY
file.dir <- "ENTER DIR"
setwd(file.dir)
#SET DEFAULT GRAPHICAL PARAMETERS
default.par <- par(no.readonly=TRUE)
Sys.setenv(TZ="America/Los_Angeles")	#need to add this to overcome issue with local (EDT) timezone and save re-specifying each time when working with date-time objects

options(max.print=999999)

# adjustments to CJS apparent survival estimates to give actual survival estimates


rescale_RSS <- function(tmp){ # Reproductive Success and Survival?
  
  Rmax <- 0.9 # upper reference point for Rmax
  meanr <- mean(tmp) # mean RSS from draws after adjustment and scaling
  sdr <- sd(tmp) # standard deviation in RSS from draws
  CVr <- sdr/meanr # CV in RSS from draws
  medianr <- median(tmp) # median RSS from draws
  sigr <- sqrt(log(CVr^2+1)) # calculated value for sigma of lognormal density function for RSS
  lnmedr <- log(medianr) # natural logarithm of median RSS
  lnRmax <- log(Rmax) # natural logarithm of upper reference point for RSS
  deltalnR <- lnRmax-lnmedr # Normal deviate in ln RSS where P(RSS > 0.9)=0.0001
  Z.99.99th <- qnorm(0.9999,mean=0,sd=1) #99.99th percentile of Normal(0,1) density function
  Sdnew <- deltalnR/Z.99.99th # adjusted lognormal sigma in RSS
  Corr <- Sdnew/sigr # correction factor to be applied to ln(RSS) deviates
  
  lnRSS <- log(tmp)
  DevlnRSS <- lnRSS-lnmedr
  adjDevlnRSS <- DevlnRSS*Corr
  newlnRSS <- lnmedr+adjDevlnRSS
  adjRSS <- exp(newlnRSS)
  return(adjRSS)
}


n.sim <- 10000

#tag loss
rss.loss.min <- 0.97
rss.loss.max <- 0.99
rss.loss <- runif(n.sim,rss.loss.min,rss.loss.max)

sas.loss.min <- 0.85
sas.loss.max <- 0.98
sas.loss <- runif(n.sim,sas.loss.min,sas.loss.max)

#tag-induced mortality
rss.mort.min <- 0.85
rss.mort.max <- 0.95
rss.mort <- runif(n.sim,rss.mort.min,rss.mort.max)

sas.mort.min <- 0.67
sas.mort.max <- 0.96
sas.mort <- runif(n.sim,sas.mort.min,sas.mort.max)

#hatchery effect
hatch.min <- 0.5
hatch.max <- 0.95
hatch <- runif(n.sim,hatch.min,hatch.max)


par(mfcol=c(2,3))
hist(rss.loss)
hist(rss.mort)
hist(sas.loss)
hist(sas.mort)
hist(hatch)
par(default.par)


#calculate rss adjustment factor for a hatchery release
rss.h.prod <- 1/(rss.loss*rss.mort*hatch)
hist(rss.h.prod)
rss.h.lnprod <- log(rss.h.prod)
min.rss.h <- min(rss.h.prod)
max.rss.h <- max(rss.h.prod)
mean.rss.h <- mean(rss.h.prod)
mean.ln.rss.h <- mean(rss.h.lnprod)
sd.rss.h <- sd(rss.h.prod)
sd.ln.rss.h <- sd(rss.h.lnprod)
cv.rss.h <- sd.rss.h / mean.rss.h
med.rss.h <- exp(mean.ln.rss.h)

#calculate sas adjustment factor for a hatchery release
sas.h.prod <- 1/(sas.loss*sas.mort*hatch)
hist(sas.h.prod)
sas.h.lnprod <- log(sas.h.prod)
min.sas.h <- min(sas.h.prod)
max.sas.h <- max(sas.h.prod)
mean.sas.h <- mean(sas.h.prod)
mean.ln.sas.h <- mean(sas.h.lnprod)
sd.sas.h <- sd(sas.h.prod)
sd.ln.sas.h <- sd(sas.h.lnprod)
cv.sas.h <- sd.sas.h/mean.sas.h
med.sas.h <- exp(mean.ln.sas.h)

#png("./RSS_SAS_adjustment_factors_hatchery_release.png", height=16, width=16, units="cm", res=300, pointsize=10, antialias="cleartype")
par(mfcol=c(2,1),mai=c(0.7,0.7,0.2,0.1))
hist(rss.h.prod,xlab="RSS adjustment factor (hatchery release)",main=NA)
hist(sas.h.prod,xlab="SAS adjustment factor (hatchery release)",main=NA)
par(default.par)
#dev.off()

#calculate rss adjustment factor for a natural-origin release
rss.n.prod <- 1/(rss.loss*rss.mort)
hist(rss.n.prod)
rss.n.lnprod <- log(rss.n.prod)
min.rss.n <- min(rss.n.prod)
max.rss.n <- max(rss.n.prod)
mean.rss.n <- mean(rss.n.prod)
mean.ln.rss.n <- mean(rss.n.lnprod)
sd.rss.n <- sd(rss.n.prod)
sd.ln.rss.n <- sd(rss.n.lnprod)
cv.rss.n <- sd.rss.n / mean.rss.n
med.rss.n <- exp(mean.ln.rss.n)

#calculate sas adjustment factor for a natural-origin release
sas.n.prod <- 1/(sas.loss*sas.mort*hatch)
hist(sas.n.prod)
sas.n.lnprod <- log(sas.n.prod)
min.sas.n <- min(sas.n.prod)
max.sas.n <- max(sas.n.prod)
mean.sas.n <- mean(sas.n.prod)
mean.ln.sas.n <- mean(sas.n.lnprod)
sd.sas.n <- sd(sas.n.prod)
sd.ln.sas.n <- sd(sas.n.lnprod)
cv.sas.n <- sd.sas.n/mean.sas.n
med.sas.n <- exp(mean.ln.sas.n)


#calculate adjustment factor for hatchery vs. natural-origin (>1 increases survival of hatchery fish)
hatch.prod <- 1/(hatch)
hist(hatch.prod)
hatch.lnprod <- log(hatch.prod)
min.hatch <- min(hatch.prod)
max.hatch <- max(hatch.prod)
mean.hatch <- mean(hatch.prod)
mean.ln.hatch <- mean(hatch.lnprod)
sd.hatch <- sd(hatch.prod)
sd.ln.hatch <- sd(hatch.lnprod)
cv.hatch <- sd.hatch/mean.hatch
med.hatch <- exp(mean.ln.hatch)

#create data frame of adjustment factors for CJS results from hatchery release
adj.h <- data.frame(rss.adj=rss.h.prod,pSUJ=rep(1,n.sim),sas.adj=sas.h.prod,pWFF=rep(1,n.sim))

#read in posterior samples from CJS model
# CH_CSM_2013_BCLTR - hatchery release
CJS_posterior <- read.csv("./North Santiam/CJS_RWM_s1sbeta(2.44,3.67)_s2dbeta(3.19,277.76)/CH_CSM_2013_BCLTR_posterior_samples.csv")
colnames(CJS_posterior) <- c("RSS","pSUJ","SAS","pWFF")

ind.10k <- sample(1:nrow(CJS_posterior),n.sim,replace=FALSE)
CJS_posterior_10k <- CJS_posterior[ind.10k,]

CJS_posterior_adj <- CJS_posterior_10k * adj.h


#need to scale adjusted RSS for each life history group in each basin
#use Zabel 2015 (p.5.11) differences between RSS for each age to scale estimate 
#obtained from release that were majority subyearling migrants
#North Santiam
rss_zabel_f_NSANT <- 0.15 #fry-spring sub
rss_zabel_s_NSANT <- 0.425 #sub to fall sub
rss_zabel_y_NSANT <- 0.6 #yrlg to 2nd spring yrlg
rss_f_scalar_NSANT <- rss_zabel_f_NSANT / rss_zabel_s_NSANT
rss_s_scalar_NSANT <- rss_zabel_s_NSANT / rss_zabel_s_NSANT
rss_y_scalar_NSANT <- rss_zabel_y_NSANT / rss_zabel_s_NSANT

#use Zabel 2015 (p.7.11) differences between RSS for each age to scale estimate 
#obtained from subyearling release that were majority subyearling migrants
#South Santiam
rss_zabel_f_SSANT <- 0.135 #fry-spring sub
rss_zabel_s_SSANT <- 0.375 #sub to fall sub
rss_zabel_y_SSANT <- 0.5 #yrlg to 2nd spring yrlg
rss_f_scalar_SSANT <- rss_zabel_f_SSANT / rss_zabel_s_SSANT
rss_s_scalar_SSANT <- rss_zabel_s_SSANT / rss_zabel_s_SSANT
rss_y_scalar_SSANT <- rss_zabel_y_SSANT / rss_zabel_s_SSANT

RSS_F_NSANT <- CJS_posterior_adj$RSS * rss_f_scalar_NSANT
RSS_S_NSANT <- CJS_posterior_adj$RSS * rss_s_scalar_NSANT
RSS_Y_NSANT <- CJS_posterior_adj$RSS * rss_y_scalar_NSANT
RSS_F_SSANT <- CJS_posterior_adj$RSS * rss_f_scalar_SSANT
RSS_S_SSANT <- CJS_posterior_adj$RSS * rss_s_scalar_SSANT
RSS_Y_SSANT <- CJS_posterior_adj$RSS * rss_y_scalar_SSANT




#rescale RSS only if max draw is >=1
print(summary(RSS_F_NSANT))
if(max(RSS_F_NSANT)>=1){
  RSS_F_NSANT <- rescale_RSS(RSS_F_NSANT)
}
print(summary(RSS_F_NSANT))

print(summary(RSS_S_NSANT))
if(max(RSS_S_NSANT)>=1){
  RSS_S_NSANT <- rescale_RSS(RSS_S_NSANT)
}
print(summary(RSS_S_NSANT))

print(summary(RSS_Y_NSANT))
if(max(RSS_Y_NSANT)>=1){
  RSS_Y_NSANT <- rescale_RSS(RSS_Y_NSANT)
}
print(summary(RSS_Y_NSANT))

print(summary(RSS_F_SSANT))
if(max(RSS_F_SSANT)>=1){
  RSS_F_SSANT <- rescale_RSS(RSS_F_SSANT)
}
print(summary(RSS_F_SSANT))

print(summary(RSS_S_SSANT))
if(max(RSS_S_SSANT)>=1){
  RSS_S_SSANT <- rescale_RSS(RSS_S_SSANT)
}
print(summary(RSS_S_SSANT))

print(summary(RSS_Y_SSANT))
if(max(RSS_Y_SSANT)>=1){
  RSS_Y_SSANT <- rescale_RSS(RSS_Y_SSANT)
}
print(summary(RSS_Y_SSANT))


#crude metohod of ensuring RSS <1
#CJS_posterior_adj$RSS <- ifelse(CJS_posterior_adj$RSS<1,CJS_posterior_adj$RSS,1)
CJS_posterior_adj$RSS_F_NSANT <- RSS_F_NSANT
CJS_posterior_adj$RSS_S_NSANT <- RSS_S_NSANT
CJS_posterior_adj$RSS_Y_NSANT <- RSS_Y_NSANT
CJS_posterior_adj$RSS_F_SSANT <- RSS_F_SSANT
CJS_posterior_adj$RSS_S_SSANT <- RSS_S_SSANT
CJS_posterior_adj$RSS_Y_SSANT <- RSS_Y_SSANT


apply(CJS_posterior,2,median)
apply(CJS_posterior_10k,2,median)


apply(adj.h,2,median)
apply(CJS_posterior_adj,2,median)
apply(CJS_posterior_adj,2,max)



write.csv(CJS_posterior_adj,"ENTER OUTPUT DIR",
          quote=FALSE,row.names=FALSE)


# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = abs(cex.cor * r))
}

# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19, col = "black", cex = 0.5)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="tan", ...)
}


png("./CH_CSM_2013_BCLTR_unadjusted_survivals.png", height=16, width=16, units="cm", res=300, pointsize=10, antialias="cleartype")
pairs(CJS_posterior_10k,
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      diag.panel=panel.hist,
      gap=0, main="CH_CSM_2013_BCLTR: unadjusted survivals")
dev.off()

png("./CH_CSM_2013_BCLTR_adjusted_survivals.png", height=16, width=16, units="cm", res=300, pointsize=10, antialias="cleartype")
pairs(CJS_posterior_adj,
      lower.panel = panel.cor,
      upper.panel = upper.panel,
      diag.panel=panel.hist,
      gap=0, main="CH_CSM_2013_BCLTR: adjusted survivals")
dev.off()



