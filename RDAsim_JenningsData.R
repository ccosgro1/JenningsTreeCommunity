setwd("~OneDrive/Documents/School/Research/Data/Trees")
library(vegan)
library(boot)
library(metafor)
logaxis = function(minlog,maxlog,side){
  pow <- seq(minlog,maxlog,by=1)
  ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
  axis(side, 10^pow, las=1, cex.axis=0.9)
  axis(side, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
}

logaxis.b = function(minlog,maxlog,side){
  pow <- seq(minlog,maxlog,by=1)
  ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
  axis(side, 10^pow,labels=NA, las=1)
  axis(side, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
}

treedata=read.csv("C:/Users/Colleen/Dropbox/Jennings 2016/New Adult Basal Area.csv")
soildata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Grad School/SpatialStats/soil data May 08 for tree plots.csv")
names(treedata)
dim(treedata)
names(soildata)
dim(soildata)

treedata[is.na(treedata)] <- 0

# The tree data was collected in three rings for each plot
# We are going to exclude stream plots (no soil data) and then sum the three rings
treedata=treedata[treedata$ecosys!="S",]
tree.ringsum=treedata[1:80,]
#This [1:80,] refers to all plots within the first ring.

for(i in 1:80) {
  tree.ringsum[i,10:43]=treedata[i,10:43]+treedata[i+83,10:43]+treedata[i+166,10:43]}
rownames(tree.ringsum)=tree.ringsum[,1]
# Let's look at a rank abbundance plot for fun
tree.colsums=apply(tree.ringsum[,10:43],2,sum)
plot(rank(tree.colsums),tree.colsums)
print(tree.colsums)

# Perform Hellinger transformation on tree species data so that RDA is based on Hellinger distance
tree.rowsums=apply(tree.ringsum[,10:43],1,sum)
tree.hel=tree.ringsum
tree.hel[,10:43]=tree.hel[,10:43]/tree.rowsums
tree.hel[,10:43]=sqrt(tree.hel[,10:43])


#### Forward selection using regionalized soil variables
# First we have to merge the datasets
tree.hel.soil=merge(tree.hel,soildata)
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
colnames(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil = rda(tree.hel.soil[,10:43], tree.hel.soil[,52:76])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared

# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,10:43], tree.hel.soil[,46:79], adjR2thresh=global.thresh)
colnames(tree.hel.soil)
hel.rda.for
names(tree.hel.soil)
hel.rda.selsoil = rda(tree.hel.soil[,10:42] ~ as.matrix(tree.hel.soil[,c(49,57,61,51,77,56)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
plot(hel.rda.selsoil)


###3. Calculating standard error of R2 based on solution found on Stats Exchange

#Function to calculate standard error of Rsqr
SER2 <- function(model){
  R2 <- (summary(model))$r.squared
  k <- summary(model)$df[1]
  n <- summary(model)$df[1]+summary(model)$df[2]
  se <- sqrt(4*R2*(1-R2)^2*(n-k-1)^2/((n^2-1)*(n+3)))
  R2.ci <- round(data.frame(R2-(abs(qt(0.05/2, n-2))*se),
                            R2,
                            R2+(abs(qt(0.05/2, n-2))*se),
                            se),
                 digits=3)
  colnames(R2.ci) <- c("Lower","R2","Upper","SE")
  return(R2.ci)
}

SER2(npp.lbbp)

#Bootstrapped R2 and confidence interval from regression
boot.R2 <- boot(biof.LBBP,function(data,indices)
  summary(lm(NPP~log(TOTMe),data[indices,]))$r.squared, R=10000)
boot.R2
boot.R2$t0
quantile(boot.R2$t,c(0.025,0.975))

#Function to calculate standard error of Rsqr
SER2.rda <- function(model){
  RDA.R2 <- RsquareAdj(model)$r.squared
  amod <- anova(model)
  n <- amod[1]$Df[1]+amod[1]$Df[2]
  k <- amod[1]$Df[1]
  se <- sqrt((4*RDA.R2*((1-RDA.R2)^2)*((n-k-1)^2))/((n^2-1)*(n+3)))
  RDA.R2.ci <- round(data.frame(RDA.R2-(abs(qt(0.05/2, n-2))*se),
                                RDA.R2,
                                RDA.R2+(abs(qt(0.05/2, n-2))*se),
                                se),
                     digits=3)
  colnames(RDA.R2.ci) <- c("Lower","R2","Upper","SE")
  return(RDA.R2.ci)
}

SER2.rda(rda.hbbp)

#Bootstrapped R2 and confidence interval from RDA
boot.R2 <- boot(bugs.LBBP,function(data,indices)
  RsquareAdj(rda(bugs.LBBP[,9:51]~log(TOTMe),data[indices,]))$r.squared, 
  R=10000)
boot.R2
boot.R2$t0
quantile(boot.R2$t,c(0.025,0.975))
sd(boot.R2$t)


###4. Testing R2 effect size using metafor


#Models for FUN (lm) and STR (RDA)
FUNmod <- list(npp.lbbp,npp.lbst,npp.hbbp,npp.hbst)
STRmod <- list(rda.lbbp,rda.lbst,rda.hbbp,rda.hbst)
#Datasets for FUN and STR
FUNdat <- list(biof.LBBP,biof.LBST,biof.HBBP,biof.HBST) 
STRdat <- list(bugs.LBBP,bugs.LBST,bugs.HBBP,bugs.HBST)

#Matrix with FUN R2, FUN SE formula, FUN SE bootstrap, STR R2, STR SE formula, STR SE bootstrap
r2.metafor <- matrix(NA, nrow=4, ncol=9, 
                     dimnames=list(c("Low-binding Big Pup","Low-binding Salmon Trout",
                                     "High binding Big Pup","High-binding Salmon Trout"),
                                   c("FUN R2","FUN SE form","FUN SE boot",
                                     "STR R2","STR SE form","STR SE boot","R2 diff",
                                     "SE form pooled","SE boot pooled")))

r2.metafor[1:4,1] <- t(sapply(FUNmod, function(x) RsquareAdj(x)$r.squared)) #FUN R2
r2.metafor[1:4,2] <- t(sapply(FUNmod, function(x) SER2(x)$SE)) #FUN SE formula
r2.metafor[1:4,3] <- t(sapply(FUNdat, function(x) sd(boot(x,function(data,indices)
  summary(lm(NPP~log(TOTMe),data[indices,]))$r.squared, R=10000)$t))) #FUN SE bootstrap
r2.metafor[1:4,4] <- t(sapply(STRmod, function(x) RsquareAdj(x)$r.squared)) #STR R2
r2.metafor[1:4,5] <- t(sapply(STRmod, function(x) SER2.rda(x)$SE)) #STR SE formula
r2.metafor[1:4,6] <- t(sapply(STRdat, function(x) sd(boot(x,function(data,indices)
  RsquareAdj(rda(x[,9:51]~log(TOTMe),data[indices,]))$r.squared, R=10000)$t))) #STR SE bootstrap
r2.metafor[1:4,7] <- r2.metafor[,1]-r2.metafor[,4] #R2 diff as FUN-STR (i.e., if negative, FUN smaller than STR)
r2.metafor[1:4,8] <- sqrt(r2.metafor[,2]+r2.metafor[,5]) #Pooled SE formula as sqrt(SE1+SE2)
r2.metafor[1:4,9] <- sqrt(r2.metafor[,3]+r2.metafor[,6]) #Pooled SE bootstrap

r2.metafor
