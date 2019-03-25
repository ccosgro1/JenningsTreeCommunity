##Cosgrove 20 March 2019
#Seedling Data Exploration

library(data.table)
library(gtools)
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


#Set NAs to 0 in seedling/juvenile datasets
juvies[is.na(juvies)]<-0
seeds[is.na(seeds)]<-0
treedata[is.na(treedata)]<-0

#Combine ecosystem type and soil variables into the seedling and juvenile datasets
#juvies=merge(juvies,soildata)
#seeds=merge(seeds,soildata)


#Stats to do:
#species by species way seedling to the adult and juvenile communities
  #regression of seedling abundance/height(?) to adult basal area proportion by plot
  #core v edge
#RDAs

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
plot(rank(scolsums),scolsums, main="Seedling Rank Abundance") #This distribution is crazy!! 


##Species Regressions
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
juvies.tree=merge(juvies,treedata, by="Plot")
juvies.tree

#Linear Regression of seedling abundance by basal area proportion
#proportion basal area
colnames(juvies.tree)
tcolsums=apply(juvies.tree[,44:78],2,sum)
propAS=(juvies.tree$ACESAC/tcolsums["ACESAC"])
beginfor=lm(juvies.tree$ACESAC_seed~propAS)
plot(beginfor)
summary(beginfor)
graph=plot(juvies.tree$ACESAC_seed~propAS, ylab="Juvenile Abundance", xlab="Proportion BA")
abline(beginfor)

colnames(juvies.tree)
for(i in 1:n){
  propBA=(juvies.tree[,6]/tcolsums["ACERUB"])
  regBA=lm(juvies.tree[,6]~propBA)
  plot(regBA)
  summary(regBA)
  graph=plot(juvies.tree[,6]~propBA, xlab="Proportion BA", ylab="Juvenile Abundance")
  coef(graph)
  abline(regBA)
}

