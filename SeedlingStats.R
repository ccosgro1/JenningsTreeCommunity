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

seeds[is.na(seeds)]<-0

#sum rings
tree.ringsum=oldtreedata[1:83,]
for(i in 1:83) {
  tree.ringsum[i,10:42]=oldtreedata[i,10:42]+oldtreedata[i+83,10:42]+oldtreedata[i+166,10:42]  }
rownames(tree.ringsum)=tree.ringsum[,1]
tree.ringsum

treedata=treedata[c(1:266,268:290),]
newtree.ringsum=treedata[1:96,]
colnames(treedata)
for(i in 1:96) {
  newtree.ringsum[i,10:44]=treedata[i,10:44]+treedata[i+96,10:44]+treedata[i+192,10:44]  }
rownames(newtree.ringsum)=newtree.ringsum$Plot
newtree.ringsum[is.na(newtree.ringsum)]<-0
newtree.ringsum

#add the name _juv to name for differentation from adult columns
addseeds=cbind(juvies[,c(1,5:35)])
#add _seed to all of the juvenile species columns
colnames(addseeds)=paste(colnames(addseeds), "juv", sep = "_")
colnames(addseeds)[1]="Plot"
colnames(addseeds)

#put everything back together and then add the adult dataset
juvies=merge(juvies[,c(1:4)], addseeds)
juvies[is.na(juvies)]<-0
juvies.tree=merge(juvies,soildata, by="Plot")
rownames(juvies.tree)=juvies.tree[,1]
head(juvies.tree)
juvies.tree=merge(juvies.tree,newtree.ringsum,by="Plot")


#to use for later analyses, after things get messed up -- subsetted,etc..
#juvies1=juvies
#seeds1=seeds

###RDAs####
#####JUVENILE#####
#Hellinger Transform data
juvies1=juvies.tree
colnames(juvies1)

tree.rowsums=apply(juvies1[5:35],1,sum)
juvies1[,5:35]=juvies1[,5:35]/tree.rowsums
juvies1[,5:35]=sqrt(juvies1[,5:35])

###RDA using ecosystem
#juvies1=juvies
colnames(juvies1)
#juvies1[is.na(juvies1)]<-0
#tree.rowsums=apply(juvies1[5:35],1,sum)
#juvies1[,5:35]=juvies1[,5:35]/tree.rowsums
#juvies1[,5:35]=sqrt(juvies1[,5:35])
#ecosyseff=rda(juvies1[,5:35],juvies1$Ecosystem)
#anova(ecosyseff)
#RsquareAdj(ecosyseff)

#juvies1=juvies1[juvies1$Transect!="Trans",]
#juvies1=juvies1[juvies1$Transect=="Core",]


#### Forward selection using regionalized SOIL variables
tree.hel.soil=juvies1
tree.hel.soil=na.exclude(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
colnames(tree.hel.soil)
dim(tree.hel.soil)
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], tree.hel.soil[,45:69])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,45:69], adjR2thresh=global.thresh)
hel.rda.for
names(tree.hel.soil)
hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(47,48,56,55,65,53)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
plot(hel.rda.selsoil)

#### Forward selection using TREE BA as variables
colnames(tree.hel.soil)
juv.pca=prcomp(tree.hel.soil[,c(78,79,84,85,97,100,103)])
juv.pca.scores=scores(juv.pca)[,1:2]
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], juv.pca.scores)
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], tree.hel.soil[,c(78,79,84,85,97,100,103)])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
#hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,45:69], adjR2thresh=global.thresh)
#hel.rda.for
#names(tree.hel.soil)
#hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(47,48,56,55,65,53)]))
#anova(hel.rda.selsoil)
#RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,76:77],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,5:35], tree.hel.all.pcnm[,114:134])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared
# now run step 2, the selection part
hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,5:35], tree.hel.all.pcnm[,114:134], adjR2thresh=global.thresh)
hel.rda.all.forpcnm
# Let's store a model using the selected variables only
names(tree.hel.all.pcnm)
hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,5:35] ~ as.matrix(tree.hel.all.pcnm[,c(114,116,122,117,133,128,125,115,118)]))
anova(hel.rda.all.selpcnm)
RsquareAdj(hel.rda.all.selpcnm)

# Now that we have these parsimonious models, NOW WE CAN PERFORM VARIANCE PARTITIONING!  Hurray!
varpart(tree.hel.soil[,5:35], tree.hel.soil[,c(47,48,56,55,65,53)],tree.hel.all.pcnm[,c(114,116,122,117,133,128,125,115,118)])


#Core Plots
tree.hel.soil=juvies1
tree.hel.soil=tree.hel.soil[tree.hel.soil$Transect.x=="Core",]
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
dim(tree.hel.soil)
colnames(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], tree.hel.soil[,c(47,48,56,55,65,53)])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared

# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,c(47,48,56,55,65,53)], adjR2thresh=global.thresh)
colnames(tree.hel.soil)
hel.rda.for
names(tree.hel.soil)
hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(47,48)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using TREE BA as variables
colnames(tree.hel.soil)
juv.pca=prcomp(tree.hel.soil[,c(78,79,84,85,97,100,103)])
juv.pca.scores=scores(juv.pca)[,1:2]
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], juv.pca.scores)
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], tree.hel.soil[,c(78,79,84,85,97,100,103)])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
#hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,45:69], adjR2thresh=global.thresh)
#hel.rda.for
#names(tree.hel.soil)
#hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(47,48,56,55,65,53)]))
#anova(hel.rda.selsoil)
#RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,76:77],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,5:35], tree.hel.all.pcnm[,114:123])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared
# now run step 2, the selection part
hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,5:35], tree.hel.all.pcnm[,114:123], adjR2thresh=global.thresh)
hel.rda.all.forpcnm
# Let's store a model using the selected variables only
names(tree.hel.all.pcnm)
hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,5:35] ~ as.matrix(tree.hel.all.pcnm[,c(114,117,115)]))
anova(hel.rda.all.selpcnm)
RsquareAdj(hel.rda.all.selpcnm)

## Variance partitioning!!
varpart(tree.hel.soil[,5:35], tree.hel.soil[,c(47,48)],tree.hel.all.pcnm[,c(114,117,115)])


#Core + edge Plots
tree.hel.soil=juvies1
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=tree.hel.soil[tree.hel.soil$Transect.x=="Core"|tree.hel.soil$Transect.x=="Edge",]
tree.hel.soil=na.exclude(tree.hel.soil)
dim(tree.hel.soil)
colnames(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], tree.hel.soil[,45:69])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared

# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,45:69], adjR2thresh=global.thresh)
colnames(tree.hel.soil)
hel.rda.for
names(tree.hel.soil)
hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(47,48,56,55)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using TREE BA as variables
colnames(tree.hel.soil)
juv.pca=prcomp(tree.hel.soil[,c(78,79,84,85,97,100,103)])
juv.pca.scores=scores(juv.pca)[,1:2]
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], juv.pca.scores)
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], tree.hel.soil[,c(78,79,84,85,97,100,103)])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
#hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,45:69], adjR2thresh=global.thresh)
#hel.rda.for
#names(tree.hel.soil)
#hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(47,48,56,55,65,53)]))
#anova(hel.rda.selsoil)
#RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,76:77],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,5:35], tree.hel.all.pcnm[,114:130])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared
# now run step 2, the selection part
#hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,5:35], tree.hel.all.pcnm[,114:123], adjR2thresh=global.thresh)
#hel.rda.all.forpcnm
# Let's store a model using the selected variables only
#names(tree.hel.all.pcnm)
#hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,5:35] ~ as.matrix(tree.hel.all.pcnm[,c(114,117,115)]))
#anova(hel.rda.all.selpcnm)
#RsquareAdj(hel.rda.all.selpcnm)

## Variance partitioning!!
#varpart(tree.hel.soil[,5:35], tree.hel.soil[,c(47,48,56,55,65,53)])


########new tree data 2016; bonus code########
###### Some Bonus tree code
# Here is a quick analysis to test for the edge effect on tree communities when Ecosystem type is the explanatory factor, just like we did before with soil data
# First the core plots
tree.hel.core=tree.hel[tree.hel$EP=="Core",]
rownames(tree.hel.core)
tree.hel.core.rda.ecosys = rda(tree.hel.core[,10:43] ~ tree.hel.core[,3])
anova(tree.hel.core.rda.ecosys)
core.arsq=RsquareAdj(tree.hel.core.rda.ecosys)
# Then the core+edge plots
tree.hel.edge=tree.hel[tree.hel$EP=="Core"|tree.hel$EP=="Edge",]
names(tree.hel.edge)
tree.hel.edge.rda.ecosys = rda(tree.hel.edge[,10:43] ~ tree.hel.edge[,3])
anova(tree.hel.edge.rda.ecosys)
edge.arsq=RsquareAdj(tree.hel.edge.rda.ecosys)
# Then the permmutation test of the edge effect
edge.eff=core.arsq$adj.r.squared/edge.arsq$adj.r.squared
# The permutation test for significance of edge effect
nperm=999
res.vec=matrix(nrow=nperm,ncol=1,0)
for(i in 1:nperm){
  dat.ec.rand=tree.hel.edge
  dat.ec.rand[dat.ec.rand$ecosys=="B",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="B",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="B",]))
  dat.ec.rand[dat.ec.rand$ecosys=="R",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="R",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="R",]))
  dat.ec.rand[dat.ec.rand$ecosys=="U",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="U",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="U",]))
  dat.ec.rand[dat.ec.rand$ecosys=="D",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="D",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="D",]))
  dat.core.rand=dat.ec.rand[dat.ec.rand$EP=="Core",]
  datcore.rda.ecosys = rda(dat.core.rand[,10:43] ~ dat.core.rand[,3])
  core.arsq.rand=RsquareAdj(datcore.rda.ecosys)
  res.vec[i]=core.arsq.rand$adj.r.squared/edge.arsq$adj.r.squared
}
tiles=rank(rbind(res.vec,edge.eff))
1-tiles[nperm+1]/(nperm+1)
# How well is each species fitted by these models?  Find out below.
goodness(tree.hel.core.rda.ecosys)
goodness(tree.hel.edge.rda.ecosys)

##Same idea as above, but Jaccard distance analysis by PCoA/distance-based RDA
jac.core.dist=vegdist(tree.hel.core[,10:43],method="jaccard",binary=TRUE)
jac.core.pcoa=cmdscale(jac.core.dist,eig=TRUE,k=41,add=TRUE)
jac.core.rda.ecosys = rda(jac.core.pcoa$points ~ tree.hel.core[,3])
anova(jac.core.rda.ecosys)
jac.core.arsq=RsquareAdj(jac.core.rda.ecosys)
jac.edge.dist=vegdist(tree.hel.edge[,10:43],method="jaccard",binary=TRUE)
jac.edge.pcoa=cmdscale(jac.edge.dist,eig=TRUE,k=70,add=TRUE)
jac.edge.pcoa.points=cbind(jac.edge.pcoa$points,tree.hel.edge)
jac.edge.rda.ecosys = rda(jac.edge.pcoa.points[,1:70] ~ jac.edge.pcoa.points[,73])
anova(jac.edge.rda.ecosys)
jac.edge.arsq=RsquareAdj(jac.edge.rda.ecosys)
jac.edge.eff=jac.core.arsq$adj.r.squared/jac.edge.arsq$adj.r.squared
nperm=999
res.vec=matrix(nrow=nperm,ncol=1,0)
for(i in 1:nperm){
  dat.ec.rand=jac.edge.pcoa.points
  dat.ec.rand[dat.ec.rand$ecosys=="B",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="B",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="B",]))
  dat.ec.rand[dat.ec.rand$ecosys=="R",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="R",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="R",]))
  dat.ec.rand[dat.ec.rand$ecosys=="U",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="U",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="U",]))
  dat.ec.rand[dat.ec.rand$ecosys=="D",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="D",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="D",]))
  dat.core.rand=dat.ec.rand[dat.ec.rand$EP=="Core",]
  datcore.rda.ecosys = rda(dat.core.rand[,1:70] ~ dat.core.rand[,73])
  core.arsq.rand=RsquareAdj(datcore.rda.ecosys)
  res.vec[i]=core.arsq.rand$adj.r.squared/jac.edge.arsq$adj.r.squared
}
tiles=rank(rbind(res.vec,jac.edge.eff))
1-tiles[nperm+1]/(nperm+1)







####SEEDLING######

dim(seeds)
colnames(seeds)
addseeds=cbind(seeds[,c(1,5:38)])
#add _seed to all of the juvenile species columns
colnames(addseeds)=paste(colnames(addseeds), "seed", sep = "_")
colnames(addseeds)[1]="Plot"
colnames(addseeds)

###RDA using ecosystem
colnames(seeds)
#seeds[is.na(seeds)]<-0
#tree.rowsums=apply(seeds[5:38],1,sum)
#seeds[,5:38]=seeds[,5:38]/tree.rowsums
#seeds[,5:38]=sqrt(seeds[,5:38])
#ecosyseff=rda(seeds[,5:38],seeds$Ecosystem)
#anova(ecosyseff)
#RsquareAdj(ecosyseff)

#seeds=seeds[seeds$Transect!="Trans",]
#seeds=seeds[seeds$Transect=="Core",]

#put everything back together and then add the adult dataset
seeds=merge(seeds[,c(1:4)], addseeds)
seeds[is.na(seeds)]<-0
dim(seeds)
seeds=merge(seeds,soildata)
dim(seeds)
seeds.tree=merge(seeds,newtree.ringsum, by="Plot")
head(seeds.tree)
#Hellinger Transform data
colnames(seeds.tree)
tree.rowsums=apply(seeds.tree[5:38],1,sum)
seeds.tree[,5:38]=seeds.tree[,5:38]/tree.rowsums
seeds.tree[,5:38]=sqrt(seeds.tree[,5:38])

#### Forward selection using regionalized soil variables
# First we have to merge the datasets
#Combine ecosystem type and soil variables into the seedling and juvenile datasets
tree.hel.soil=seeds.tree
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
dim(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
colnames(tree.hel.soil)
hel.rda.fulsoil_long = rda(tree.hel.soil[,5:38], tree.hel.soil[,48:55]) #variables_long. significant: p=0.006
anova(hel.rda.fulsoil_long)
hel.rda.fulsoil_nug=rda(tree.hel.soil[,5:38], tree.hel.soil[,c(56:64)]) #variables_nug. not significant.
anova(hel.rda.fulsoil_nug)
hel.rda.fulsoil_short = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(65:72)]) #variables_long. significant: p=0.006
anova(hel.rda.fulsoil_short)
RsquareAdj(hel.rda.fulsoil_long)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,39:72], tree.hel.soil[,48:55], adjR2thresh=global.thresh)
hel.rda.for
names(tree.hel.soil)
hel.rda.selsoil = rda(tree.hel.soil[,39:72] ~ as.matrix(tree.hel.soil[,c(51)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using TREE BA as variables
colnames(tree.hel.soil)
juv.pca=prcomp(tree.hel.soil[,c(81,82,87,91,92,100,103,107)])
juv.pca.scores=scores(juv.pca)[,1:2]
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], juv.pca.scores)
hel.rda.fulsoil = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(81,82,87,91,92,100,103,107)])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
#hel.rda.for=forward.sel(tree.hel.soil[,5:38], tree.hel.soil[,c(81,82,87,91,92,100,103,107)], adjR2thresh=global.thresh)
#hel.rda.for
#names(tree.hel.soil)
#hel.rda.selsoil = rda(tree.hel.soil[,5:38] ~ as.matrix(tree.hel.soil[,c(91,81)]))
#anova(hel.rda.selsoil)
#RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,79:80],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,5:38], tree.hel.all.pcnm[,117:135])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared
# now run step 2, the selection part
#hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,5:38], tree.hel.all.pcnm[,114:123], adjR2thresh=global.thresh)
#hel.rda.all.forpcnm
# Let's store a model using the selected variables only
#names(tree.hel.all.pcnm)
#hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,5:38] ~ as.matrix(tree.hel.all.pcnm[,c(114,117,115)]))
#anova(hel.rda.all.selpcnm)
#RsquareAdj(hel.rda.all.selpcnm)

## Variance partitioning!!
#varpart(tree.hel.soil[,5:35], tree.hel.soil[,c(47,48,56,55,65,53)],tree.hel.soil[,c(91,81)])

#Core Plots
tree.hel.soil=tree.hel.soil[tree.hel.soil$Transect.x=="Core",]
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
dim(tree.hel.soil)
colnames(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil_long = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(48:55)]) #variables_long. significant: p=0.006
anova(hel.rda.fulsoil_long)
#hel.rda.fulsoil_nug=rda(tree.hel.soil[,5:38], tree.hel.soil[,c(56:64)]) #variables_nug. not significant.
#anova(hel.rda.fulsoil_nug)
hel.rda.fulsoil_short = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(65:72)]) #variables_long. significant: p=0.006
anova(hel.rda.fulsoil_short)
RsquareAdj(hel.rda.fulsoil_long)
global.thresh=RsquareAdj(hel.rda.fulsoil_long)$adj.r.squared

# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,5:38], tree.hel.soil[,c(48:55,65:72)], adjR2thresh=global.thresh)
colnames(tree.hel.soil)
hel.rda.for
hel.rda.selsoil = rda(tree.hel.soil[,5:38] ~ as.matrix(tree.hel.soil[,c(49)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using TREE BA as variables
colnames(tree.hel.soil)
juv.pca=prcomp(tree.hel.soil[,c(81,82,87,91,92,100,103,107)])
juv.pca.scores=scores(juv.pca)[,1:2]
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], juv.pca.scores)
hel.rda.fulsoil = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(81,82,87,91,92,100,103,107)])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
#hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,45:69], adjR2thresh=global.thresh)
#hel.rda.for
#names(tree.hel.soil)
#hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(47,48,56,55,65,53)]))
#anova(hel.rda.selsoil)
#RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,79:80],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,5:38], tree.hel.all.pcnm[,117:125])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared
# now run step 2, the selection part
hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,5:38], tree.hel.all.pcnm[,117:125], adjR2thresh=global.thresh)
hel.rda.all.forpcnm
# Let's store a model using the selected variables only
names(tree.hel.all.pcnm)
hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,5:38] ~ as.matrix(tree.hel.all.pcnm[,c(120)]))
anova(hel.rda.all.selpcnm)
RsquareAdj(hel.rda.all.selpcnm)

## Variance partitioning!!
varpart(tree.hel.soil[,5:38], tree.hel.soil[,c(47,48,56,55,65,53)],tree.hel.all.pcnm[,c(120)])

#Core + edge Plots
tree.hel.soil=seeds.tree
tree.hel.soil=tree.hel.soil[tree.hel.soil$Transect.x=="Core"|tree.hel.soil$Transect.x=="Edge",]
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
colnames(tree.hel.soil)
dim(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil_long = rda(tree.hel.soil[,5:38], tree.hel.soil[,48:55]) #variables_long. significant: p=0.006
anova(hel.rda.fulsoil_long)
hel.rda.fulsoil_nug=rda(tree.hel.soil[,5:38], tree.hel.soil[,c(56:64)]) #variables_nug. not significant.
anova(hel.rda.fulsoil_nug)
hel.rda.fulsoil_short = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(65:72)]) #variables_long. significant: p=0.006
anova(hel.rda.fulsoil_short)

RsquareAdj(hel.rda.fulsoil_long)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared

# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,5:38], tree.hel.soil[,48:55], adjR2thresh=global.thresh)
colnames(tree.hel.soil)
hel.rda.for
hel.rda.selsoil = rda(tree.hel.soil[,5:38] ~ as.matrix(tree.hel.soil[,c(55,48)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using TREE BA as variables
colnames(tree.hel.soil)
juv.pca=prcomp(tree.hel.soil[,c(81,82,87,91,92,100,103,107)])
juv.pca.scores=scores(juv.pca)[,1:2]
hel.rda.fulsoil = rda(tree.hel.soil[,5:35], juv.pca.scores)
hel.rda.fulsoil = rda(tree.hel.soil[,5:38], tree.hel.soil[,c(81,82,87,91,92,100,103,107)])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
#hel.rda.for=forward.sel(tree.hel.soil[,5:35], tree.hel.soil[,45:69], adjR2thresh=global.thresh)
#hel.rda.for
#names(tree.hel.soil)
#hel.rda.selsoil = rda(tree.hel.soil[,5:35] ~ as.matrix(tree.hel.soil[,c(47,48,56,55,65,53)]))
#anova(hel.rda.selsoil)
#RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,79:80],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,5:38], tree.hel.all.pcnm[,117:132])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared
# now run step 2, the selection part
#hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,5:38], tree.hel.all.pcnm[,117:125], adjR2thresh=global.thresh)
#hel.rda.all.forpcnm
# Let's store a model using the selected variables only
#names(tree.hel.all.pcnm)
#hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,5:38] ~ as.matrix(tree.hel.all.pcnm[,c(120)]))
#anova(hel.rda.all.selpcnm)
#RsquareAdj(hel.rda.all.selpcnm)

## Variance partitioning!!
#varpart(tree.hel.soil[,5:35], tree.hel.soil[,c(47,48,56,55,65,53)])
