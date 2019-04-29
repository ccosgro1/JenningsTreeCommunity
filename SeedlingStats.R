##Cosgrove 20 March 2019
#Seedling Data Exploration

library(data.table)
library(gtools)
library(vegan)
library(packfor)
library(MASS)
##Load the Data
juvies=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Research/Data/Trees/JuvenileMatrix.csv")
seeds=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Research/Data/Trees/SeedlingMatrix.csv")
soildata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Grad School/SpatialStats/soil data May 08 for tree plots.csv")
oldtreedata=read.csv("C:/Users/Colleen/Dropbox/Jennings 2016/tree basal areas.csv")
treedata=read.csv("C:/Users/Colleen/Dropbox/Jennings 2016/New Adult Basal Area.csv")


#Define rownames
dim(juvies)
rownames(juvies)=juvies[,1]
rownames(seeds)=seeds[,1]
rownames(seeds)
for(i in 1:83) {
  tree.ringsum[i,10:42]=oldtreedata[i,10:42]+oldtreedata[i+83,10:42]+oldtreedata[i+166,10:42]  }
rownames(tree.ringsum)=tree.ringsum[,1]
tree.ringsum

#Set NAs to 0 in seedling/juvenile datasets
juvies[is.na(juvies)]<-0
seeds[is.na(seeds)]<-0
treedata[is.na(treedata)]<-0

#to use for later analyses, after things get messed up -- subsetted,etc..
juvies1=juvies
seeds1=seeds

###RDAs####
#####JUVENILE#####
#Hellinger Transform data
colnames(juvies1)
tree.rowsums=apply(juvies1[5:35],1,sum)
juvies1[,5:35]=juvies1[,5:35]/tree.rowsums
juvies1[,5:35]=sqrt(juvies1[,5:35])

#### Forward selection using regionalized soil variables
# First we have to merge the datasets
#Combine ecosystem type and soil variables into the seedling and juvenile datasets
tree.hel.soil=merge(juvies1,soildata)
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
colnames(tree.hel.soil)
dim(tree.hel.soil)
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], tree.hel.soil[,47:71])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,47:71], adjR2thresh=global.thresh)
hel.rda.for
names(tree.hel.soil)
hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(49,50,58,57,67,55)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
plot(hel.rda.selsoil)

#Core Plots
tree.hel.soil=merge(juvies1,soildata)
tree.hel.soil=tree.hel.soil[tree.hel.soil$Transect=="Core",]
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
colnames(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], tree.hel.soil[,c(49,50,58,57,67,55)])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared

# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,c(49,50,58,57,67,55)], adjR2thresh=global.thresh)
colnames(tree.hel.soil)
hel.rda.for
names(tree.hel.soil)
hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(49,50)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
plot(hel.rda.selsoil)


#Core + edge Plots
tree.hel.soil=merge(juvies1,soildata)
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=tree.hel.soil[tree.hel.soil$Transect=="Core"|tree.hel.soil$Transect=="Edge",]
tree.hel.soil=na.exclude(tree.hel.soil)
colnames(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], tree.hel.soil[,47:71])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared

# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,47:71], adjR2thresh=global.thresh)
colnames(tree.hel.soil)
hel.rda.for
names(tree.hel.soil)
hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(49,50,58,57)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
plot(hel.rda.selsoil)


####SEEDLING######
#Hellinger Transform data
colnames(seeds1)
tree.rowsums=apply(seeds1[5:38],1,sum)
seeds1[,5:38]=seeds1[,5:38]/tree.rowsums
seeds1[,5:38]=sqrt(seeds1[,5:38])

#### Forward selection using regionalized soil variables
# First we have to merge the datasets
#Combine ecosystem type and soil variables into the seedling and juvenile datasets
tree.hel.soil=merge(seeds1,soildata)
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
colnames(tree.hel.soil)
hel.rda.fulsoil = rda(tree.hel.soil[,5:38], tree.hel.soil[,49:73])
anova(hel.rda.fulsoil) #NOT SIGNIFICANT SO CANNOT RUN RDA ON THIS DATA
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
#hel.rda.for=forward.sel(tree.hel.soil[,39:72], tree.hel.soil[,11:35], adjR2thresh=global.thresh)
#hel.rda.for
#names(tree.hel.soil)
#hel.rda.selsoil = rda(tree.hel.soil[,39:72] ~ as.matrix(tree.hel.soil[,c(18,13,31)]))
#anova(hel.rda.selsoil)
#RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#Core Plots
tree.hel.soil=merge(seeds1,soildata)
tree.hel.soil=tree.hel.soil[tree.hel.soil$Transect=="Core",]
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
colnames(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil_long = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(49:56)]) #variables_long. significant: p=0.006
anova(hel.rda.fulsoil_long)
#hel.rda.fulsoil_nug = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(57:65)]) #variables_nug. not significant.
#anova(hel.rda.fulsoil_nug)
#hel.rda.fulsoil_short = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(66:73)]) #variables_long. significant: p=0.006
#anova(hel.rda.fulsoil_short)
RsquareAdj(hel.rda.fulsoil_long)
global.thresh=RsquareAdj(hel.rda.fulsoil_long)$adj.r.squared

# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,5:38], tree.hel.soil[,49:56], adjR2thresh=global.thresh)
colnames(tree.hel.soil)
hel.rda.for
hel.rda.selsoil = rda(tree.hel.soil[,5:38] ~ as.matrix(tree.hel.soil[,c(56,51)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
plot(hel.rda.selsoil)

#Core + edge Plots
scolsums=apply(seeds1[,5:38],2,sum)
tree.hel.soil=merge(seeds1,soildata)
tree.hel.soil=tree.hel.soil[tree.hel.soil$Transect=="Core"|tree.hel.soil$Transect=="Edge",]
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
colnames(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(47,49:53,55,57:62,64,66:70)])
anova(hel.rda.fulsoil) #not significant
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared

# Here is the forward selection command and then constructing a model with selected variables.
#hel.rda.for=forward.sel(tree.hel.soil[,5:38], tree.hel.soil[,11:35], adjR2thresh=global.thresh)
#colnames(tree.hel.soil)
#hel.rda.for
#hel.rda.selsoil = rda(tree.hel.soil[,5:38] ~ as.matrix(tree.hel.soil[,c(18,13,31)]))
#anova(hel.rda.selsoil)
#RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

