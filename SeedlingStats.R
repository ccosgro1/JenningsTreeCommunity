##Cosgrove 20 March 2019
#Seedling Data Exploration

library(data.table)
library(gtools)
library(vegan)
library(packfor)
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

##ANOVAs: Juvenile abundance by ecosystem then by core/edge BY PLOT
dim(juvies)
colnames(juvies)
juvie.rowsums=apply(juvies[,5:35],1,sum)
jecoabund=aov(juvies$Grand.Total~juvies$TerrType) #Significant. p=0.0183
summary(jecoabund)
plot(jecoabund)
TukeyHSD(jecoabund)

jtransabund=aov(juvies$Grand.Total~juvies$Transect)
summary(jtransabund) #not significant. p=0.6

##ANOVAs: Seedling abundance by ecosystem then by core/edge BY PLOT
dim(seeds)
colnames(seeds)
seed.rowsums=apply(seeds[,5:38], 1, sum)
secoabund=aov(seeds$Grand.Total~seeds$TerrType) #Not significant. p=0.4
summary(secoabund)

stransabund=aov(seeds$Grand.Total~seeds$Transect) #Also, not significant. p=0.3
summary(stransabund)

##Rank abudance graphs
jcolsums=apply(juvies[,5:35],2,sum)
plot(rank(jcolsums),jcolsums, main="Juvenile Rank Abundance")
barplot(jcolsums)
scolsums=apply(seeds[,5:38],2,sum)
barplot(scolsums)
plot(rank(scolsums),scolsums, main="Seedling Rank Abundance") #This distribution is crazy!! 

#Correlation between seedlings and juveniles
cor(juvies.tree[,5:25], seeds.tree[,5:28], method="pearson")


####Species Regressions####
#Juvenile Regression
#FIRST: change names in juvies so that merging with adult dataset there aren't duplicate columns.
dim(juvies)
colnames(juvies)
addseeds=cbind(juvies[,c(1,5:35)])
#add _seed to all of the juvenile species columns
colnames(addseeds)=paste(colnames(addseeds), "seed", sep = "_")
colnames(addseeds)[1]="Plot"
colnames(addseeds)
#put everything back together and then add the adult dataset
juvies=merge(juvies[,c(1:4)], addseeds)
juvies.tree=merge(juvies,tree.ringsum, by="Plot")
head(juvies.tree)

##Linear Regression of seedling abundance by basal area proportion
#Remove columns that don't have both juveniles and adult trees.
colnames(juvies.tree)
juvies.tree=juvies.tree[,c(1:4,6:12,15:18,20:24,26:30,31:45,48:50,52:60,62:63,66:70,73:76)]
colnames(juvies.tree)
tcolsums=apply(juvies.tree[,39:63],2,sum)
#Remove columns==0 (found using tcolsums)
juvies.tree=juvies.tree[,c(1:9,11:17,19:25,27,29,30:43,45:51,53:60,62)]

#proportion basal area
colnames(juvies.tree)
tcolsums=apply(juvies.tree[,35:55],2,sum)
rownames(juvies.tree)=juvies.tree[,1]
trowsums=apply(juvies.tree[,35:55],1,sum)

juvab=decostand(juvies.tree[,35:55],method = "total" ,MARGIN = 1)
colnames(juvies.tree)
tcolsums
trowsums
#Global Comparison without proportion of juveniles.
for(i in 5:25){
  propBA=(juvies.tree[,i+30]/tcolsums[i-4])
  regBA=lm(juvies.tree[,i]~propBA)
  summary(regBA)
  graph=plot(juvies.tree[,i]~propBA, xlab="Proportion BA", ylab="Juvenile Abundance", main=colnames(juvies.tree)[i+30])
  abline(regBA)
}

#Local comparison without proportion of juveniles
for(i in 5:25){
  propBA=(juvies.tree[,i+30]/tcolsums[i-4])
  regBA=lm(juvies.tree[,i]~juvab[,i-4])
  summary(regBA)
  graph=plot(juvies.tree[,i]~juvab[,i-4], xlab="Proportion BA", ylab="Juvenile Abundance", main=colnames(juvies.tree)[i+30])
  abline(regBA)
}

#Global Comparison WITH proportion juveniles
jcolsums=apply(juvies.tree[,5:25],2,sum)
for(i in 5:25){
  propjuveniles=juvies.tree[,i]/jcolsums[i-4]
  propBA=(juvies.tree[,i+30]/tcolsums[i-4])
  regBA=lm(propjuveniles~propBA)
  summary(regBA)
  graph=plot(propjuveniles~propBA, xlab="Proportion BA", ylab="Juvenile Proportion", main=colnames(juvies.tree)[i+30])
  abline(regBA)
}

#Local Comparison WITH juvenile proportion
for(i in 5:25){
  propjuveniles=juvies.tree[,i]/jcolsums[i-4]
  regBA=lm(propjuveniles~juvab[,i-4])
  summary(regBA)
  graph=plot(propjuveniles~juvab[,i-4], xlab="Proportion BA", ylab="Juvenile Abundance", main=colnames(juvies.tree)[i+30])
  abline(regBA)
}


#Seedling Regression
dim(seeds)
colnames(seeds)
addseeds=cbind(seeds[,c(1,5:38)])
#add _seed to all of the juvenile species columns
colnames(addseeds)=paste(colnames(addseeds), "seed", sep = "_")
colnames(addseeds)[1]="Plot"
colnames(addseeds)
#put everything back together and then add the adult dataset
seeds=merge(seeds[,c(1:4)], addseeds)
seeds.tree=merge(seeds,tree.ringsum, by="Plot")
head(seeds.tree)

##Linear Regression of seedling abundance by basal area proportion
#Remove columns that don't have both juveniles and adult trees.
colnames(seeds.tree)
seeds.tree=seeds.tree[,c(1:4,6:11,13,15,18:21,23:26,28:32,34:48,50:53,55:62,64:65,69:73,75:76,78:80)]
colnames(seeds.tree)
tcolsums=apply(seeds.tree[,39:64],2,sum)
#Remove columns==0 (found using tcolsums)
seeds.tree=seeds.tree[,c(1:10,12:28,30:44,46:62,64)]

#proportion basal area
colnames(seeds.tree)
tcolsums=apply(seeds.tree[,37:60],2,sum)

colnames(seeds.tree)
tcolsums
#Global Comparison with no proportion seedlings
for(i in 5:28){
  propBA=(seeds.tree[,i+32]/tcolsums[i-4])
  regBA=lm(seeds.tree[,i]~propBA)
  summary(regBA)
  graph=plot(seeds.tree[,i]~propBA, xlab="Proportion BA", ylab="Seedling Abundance", main=colnames(seeds.tree)[i+32])
  abline(regBA)
}

seedab=decostand(seeds.tree[,37:60],method = "total" ,MARGIN = 1)
#Local comparison with no proportion seedlings
for(i in 5:28){
  regBA=lm(seeds.tree[,i]~seedab[,i-4])
  summary(regBA)
  graph=plot(seeds.tree[,i]~seedab[,i-4], xlab="Proportion BA", ylab="Seedling Abundance", main=colnames(seeds.tree)[i+32])
  abline(regBA)
}

#Global comparison WITH proportion seedlings
colnames(seeds.tree)
scolsums=apply(seeds.tree[5:28],2,sum)
for(i in 5:28){
  propseedling=seeds.tree[,i]/scolsums[i-4]
  propBA=(seeds.tree[,i+32]/tcolsums[i-4])
  regBA=lm(propseedling~propBA)
  summary(regBA)
  graph=plot(propseedling~propBA, xlab="Proportion BA", ylab="Seedling Proportion", main=colnames(seeds.tree)[i+32])
  abline(regBA)
}

#Local comparison WITH proportion seedlings
for(i in 5:28){
  propseedling=seeds.tree[,i]/scolsums[i-4]
  regBA=lm(propseedling~seedab[,i-4])
  summary(regBA)
  graph=plot(propseedling~seedab[,i-4], xlab="Proportion BA", ylab="Seedling Proportion", main=colnames(seeds.tree)[i+32])
  abline(regBA)
}

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
