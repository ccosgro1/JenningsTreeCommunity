# Citation for this R code
# Blackwood, C.B., K.A. Smemo, M.W. Kershner, L.M. Feinstein, O.J. Valverde-Barrantes. Decay of ecosystem differences and decoupling of tree community-soil environment relationships at ecotones. Ecological Monographs 83:403-417.
#This is to demonstrate forward selection on the Jennings tree data using soil variables and PCNM vectors as predictors
#attempt to push to github

library(vegan)
library(packfor)
# Although the code below does not, you should use the better PCNM library if you can, to select significantly structured PCNMs before anything else with them
#library(PCNM)

#### OLD DATA 2008 ######
treedata=read.csv("C:/Users/Colleen/Dropbox/Jennings 2016/tree basal areas.csv")
soildata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Grad School/SpatialStats/soil data May 08 for tree plots.csv")
names(treedata)
dim(treedata)
names(soildata)
dim(soildata)
# The tree data was collected in three rings for each plot
# We are going to exclude stream plots (no soil data) and then sum the three rings
treedata=treedata[treedata$ecosys!="S",]
tree.ringsum=treedata[1:83,]
for(i in 1:83) {
	tree.ringsum[i,10:42]=treedata[i,10:42]+treedata[i+83,10:42]+treedata[i+166,10:42]  }
rownames(tree.ringsum)=tree.ringsum[,1]
# Let's look at a rank abbundance plot for fun
tree.colsums=apply(tree.ringsum[,10:42],2,sum)
plot(rank(tree.colsums),tree.colsums, xlab="tree species")

# Perform Hellinger transformation on tree species data so that RDA is based on Hellinger distance
tree.rowsums=apply(tree.ringsum[,10:42],1,sum)
tree.hel=tree.ringsum
tree.hel[,10:42]=tree.hel[,10:42]/tree.rowsums
tree.hel[,10:42]=sqrt(tree.hel[,10:42])

#### Forward selection using regionalized soil variables
# First we have to merge the datasets
tree.hel.soil=merge(tree.hel,soildata)
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil = rda(tree.hel.soil[,10:42], tree.hel.soil[,52:76])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,10:42], tree.hel.soil[,52:76], adjR2thresh=global.thresh)
hel.rda.for
names(tree.hel.soil)
hel.rda.selsoil = rda(tree.hel.soil[,10:42] ~ as.matrix(tree.hel.soil[,c(52,54,56,57,58)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
plot(hel.rda.selsoil)
# Cool analysis, but what a terrible plot.  Much better to get the scores and plot them in a real graphing program, or do some fancy graphing in R if that is your thing.

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,8:9],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
# If you can get the PCNM library to run, you could do the line below to select only significantly structured PCNMs
#tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs[,which(all.pcnm1$Moran_I[,2]<0.05)])
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,77:120])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared
# now run step 2, the selection part
hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,77:120], adjR2thresh=global.thresh)
hel.rda.all.forpcnm
# Let's store a model using the selected variables only
names(tree.hel.all.pcnm)
hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,10:42] ~ as.matrix(tree.hel.all.pcnm[,c(77:81,83:86)]))
anova(hel.rda.all.selpcnm)
RsquareAdj(hel.rda.all.selpcnm)

# Now that we have these parsimonious models, NOW WE CAN PERFORM VARIANCE PARTITIONING!  Hurray!
varpart(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(52,54,56,57,58)])
# The variance explained by soil variables and spatial structure clearly overlap. 
# It also appears that there is more independent spatial structure in the trees than independent structuring by soil.
# Are the soil variables significant after accounting for spatial structure, and vice versa?
# Let's see.
hel.rda.all.partpcnm=rda(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(52,54,56,57,58)])
anova(hel.rda.all.partpcnm)
# Yes.
hel.rda.all.partsoil=rda(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,c(52,54,56,57,58)], tree.hel.all.pcnm[,c(77:81,83:86)])
anova(hel.rda.all.partsoil)
# Also yes.  But only by p-value; that is really not much variance explained anymore.

### Can you try to run this on core only plots or core+edge plots and see how the importance of soil versus spatial structure shifts?


######Old Data 2008; Some Bonus tree code####
# Here is a quick analysis to test for the edge effect on tree communities when Ecosystem type is the explanatory factor, just like we did before with soil data
# First the core plots
tree.hel.core=tree.hel[tree.hel$EP=="Core",]
tree.hel.core.rda.ecosys = rda(tree.hel.core[,10:42] ~ tree.hel.core[,3])
anova(tree.hel.core.rda.ecosys)
core.arsq=RsquareAdj(tree.hel.core.rda.ecosys)
# Then the core+edge plots
tree.hel.edge=tree.hel[tree.hel$EP=="Core"|tree.hel$EP=="Edge",]
tree.hel.edge.rda.ecosys = rda(tree.hel.edge[,10:42] ~ tree.hel.edge[,3])
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
	datcore.rda.ecosys = rda(dat.core.rand[,10:42] ~ dat.core.rand[,3])
	core.arsq.rand=RsquareAdj(datcore.rda.ecosys)
	res.vec[i]=core.arsq.rand$adj.r.squared/edge.arsq$adj.r.squared
}
tiles=rank(rbind(res.vec,edge.eff))
1-tiles[nperm+1]/(nperm+1)
# How well is each species fitted by these models?  Find out below.
goodness(tree.hel.core.rda.ecosys)
goodness(tree.hel.edge.rda.ecosys)

##Same idea as above, but Jaccard distance analysis by PCoA/distance-based RDA
jac.core.dist=vegdist(tree.hel.core[,10:42],method="jaccard",binary=TRUE)
jac.core.pcoa=cmdscale(jac.core.dist,eig=TRUE,k=41,add=TRUE)
jac.core.rda.ecosys = rda(jac.core.pcoa$points ~ tree.hel.core[,3])
anova(jac.core.rda.ecosys)
jac.core.arsq=RsquareAdj(jac.core.rda.ecosys)
jac.edge.dist=vegdist(tree.hel.edge[,10:42],method="jaccard",binary=TRUE)
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

goodness(jac.core.rda.ecosys)
goodness(jac.edge.rda.ecosys)



######## NEW DATA 2016 ######
treedata=read.csv("C:/Users/Colleen/Dropbox/Jennings 2016/New Adult Basal Area.csv")
soildata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Grad School/SpatialStats/soil data May 08 for tree plots.csv")
names(treedata)
dim(treedata)
names(soildata)
dim(soildata)
treedata$ecosys

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
# Cool analysis, but what a terrible plot.  Much better to get the scores and plot them in a real graphing program, or do some fancy graphing in R if that is your thing.

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,8:9],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
# If you can get the PCNM library to run, you could do the line below to select only significantly structured PCNMs
#tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs[,which(all.pcnm1$Moran_I[,2]<0.05)])
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,80:122])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared
# now run step 2, the selection part
hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,80:121], adjR2thresh=global.thresh)
hel.rda.all.forpcnm
# Let's store a model using the selected variables only
names(tree.hel.all.pcnm)
hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,10:43] ~ as.matrix(tree.hel.all.pcnm[,c(80:82,85:88)]))
anova(hel.rda.all.selpcnm)
RsquareAdj(hel.rda.all.selpcnm)

# Now that we have these parsimonious models, NOW WE CAN PERFORM VARIANCE PARTITIONING!  Hurray!
colnames(tree.hel.all.pcnm)
varpart(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(86,80,88,82,81,90)])
# The variance explained by soil variables and spatial structure clearly overlap. 
# It also appears that there is more independent spatial structure in the trees than independent structuring by soil.
# Are the soil variables significant after accounting for spatial structure, and vice versa?
# Let's see.
hel.rda.all.partpcnm=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(86,80,88,82,81,90)])
anova(hel.rda.all.partpcnm)
# Yes.
hel.rda.all.partsoil=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(86,80,88,82,81,90)], tree.hel.all.pcnm[,c(77:81,83:86)])
anova(hel.rda.all.partsoil)
# Also yes.  But only by p-value; that is really not much variance explained anymore.


colnames(tree.hel.all.pcnm)
varpart(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(86,80,88,82,81,90)])


names(tree.hel.all.pcnm)
plot(hel.rda.all.partpcnm)

### Can you try to run this on core only plots or core+edge plots and see how the importance of soil versus spatial structure shifts?

#Core Plots
tree.hel.soil=merge(tree.hel,soildata)
tree.hel.soil=tree.hel.soil[tree.hel.soil$EP=="Core",]
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
hel.rda.selsoil = rda(tree.hel.soil[,10:42] ~ as.matrix(tree.hel.soil[,c(55,78,61,77,73)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
plot(hel.rda.selsoil)
# Cool analysis, but what a terrible plot.  Much better to get the scores and plot them in a real graphing program, or do some fancy graphing in R if that is your thing.

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,8:9],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
# If you can get the PCNM library to run, you could do the line below to select only significantly structured PCNMs
#tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs[,which(all.pcnm1$Moran_I[,2]<0.05)])
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,80:122])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared
# now run step 2, the selection part
hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,80:121], adjR2thresh=global.thresh)
hel.rda.all.forpcnm
# Let's store a model using the selected variables only
names(tree.hel.all.pcnm)
hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,10:43] ~ as.matrix(tree.hel.all.pcnm[,c(80:82,85:88)]))
anova(hel.rda.all.selpcnm)
RsquareAdj(hel.rda.all.selpcnm)

# Now that we have these parsimonious models, NOW WE CAN PERFORM VARIANCE PARTITIONING!  Hurray!
colnames(tree.hel.all.pcnm)
varpart(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(86,80,88,82,81,90)])
# The variance explained by soil variables and spatial structure clearly overlap. 
# It also appears that there is more independent spatial structure in the trees than independent structuring by soil.
# Are the soil variables significant after accounting for spatial structure, and vice versa?
# Let's see.
hel.rda.all.partpcnm=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(86,80,88,82,81,90)])
anova(hel.rda.all.partpcnm)
# Yes.
hel.rda.all.partsoil=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(86,80,88,82,81,90)], tree.hel.all.pcnm[,c(77:81,83:86)])
anova(hel.rda.all.partsoil)
# Also yes.  But only by p-value; that is really not much variance explained anymore.


#Core+Edge
tree.hel.soil=merge(tree.hel,soildata)
tree.hel.soil=tree.hel.soil[tree.hel.soil$EP=="Core"|tree.hel.soil$EP=="Edge",]
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
hel.rda.selsoil = rda(tree.hel.soil[,10:42] ~ as.matrix(tree.hel.soil[,c(49,57,61,78)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
plot(hel.rda.selsoil)
# Cool analysis, but what a terrible plot.  Much better to get the scores and plot them in a real graphing program, or do some fancy graphing in R if that is your thing.

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,8:9],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
# If you can get the PCNM library to run, you could do the line below to select only significantly structured PCNMs
#tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs[,which(all.pcnm1$Moran_I[,2]<0.05)])
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,80:118])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared
# now run step 2, the selection part
hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,80:118], adjR2thresh=global.thresh)
hel.rda.all.forpcnm
# Let's store a model using the selected variables only
names(tree.hel.all.pcnm)
hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,10:43] ~ as.matrix(tree.hel.all.pcnm[,c(80:81,83,85:89)]))
anova(hel.rda.all.selpcnm)
RsquareAdj(hel.rda.all.selpcnm)

# Now that we have these parsimonious models, NOW WE CAN PERFORM VARIANCE PARTITIONING!  Hurray!
colnames(tree.hel.all.pcnm)
varpart(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(86,80,88,82,81,90)])
# The variance explained by soil variables and spatial structure clearly overlap. 
# It also appears that there is more independent spatial structure in the trees than independent structuring by soil.
# Are the soil variables significant after accounting for spatial structure, and vice versa?
# Let's see.
hel.rda.all.partpcnm=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(86,80,88,82,81,90)])
anova(hel.rda.all.partpcnm)
# Yes.
hel.rda.all.partsoil=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(86,80,88,82,81,90)], tree.hel.all.pcnm[,c(77:81,83:86)])
anova(hel.rda.all.partsoil)
# Also yes.  But only by p-value; that is really not much variance explained anymore.


###This is repeat code from up above - and is not actually helpful for any reason whatsoever right now.
# Now that we have these parsimonious models, NOW WE CAN PERFORM VARIANCE PARTITIONING!  Hurray!
varpart(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(52,54,56,57,58)])
# The variance explained by soil variables and spatial structure clearly overlap. 
# It also appears that there is more independent spatial structure in the trees than independent structuring by soil.
# Are the soil variables significant after accounting for spatial structure, and vice versa?
# Let's see.
hel.rda.all.partpcnm=rda(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(52,54,56,57,58)])
anova(hel.rda.all.partpcnm)
# Yes.
hel.rda.all.partsoil=rda(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,c(52,54,56,57,58)], tree.hel.all.pcnm[,c(77:81,83:86)])
anova(hel.rda.all.partsoil)
# Also yes.  But only by p-value; that is really not much variance explained anymore.


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
	datcore.rda.ecosys = rda(dat.core.rand[,10:42] ~ dat.core.rand[,3])
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






###Species level decoupling###
#Mortality

treedata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Research/Data/Trees/Species Mortality Matrix.csv")
soildata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Grad School/SpatialStats/soil data May 08 for tree plots.csv")
names(treedata)
rownames(treedata)=treedata[,1]
dim(treedata)
names(soildata)
rownames(soildata)=soildata[,1]
dim(soildata)
# The tree data was collected in three rings for each plot
# We are going to exclude stream plots (no soil data) and then sum the three rings - there's no need to sum the rings here.
treedata=treedata[treedata$Ecosystem!="S",]
dim(treedata)
treedata[is.na(treedata)] <- 0
treedata

tree.ringsum=treedata

##for(i in 1:83) {
##	tree.ringsum[i,10:42]=treedata[i,10:42]+treedata[i+83,10:42]+treedata[i+166,10:42]  }

#
rownames(tree.ringsum)=tree.ringsum[,1]
# Let's look at a rank abbundance plot for fun

tree.colsums=apply(tree.ringsum[,10:42],2,sum)
plot(rank(tree.colsums),tree.colsums)

# Perform Hellinger transformation on tree species data so that RDA is based on Hellinger distance
tree.rowsums=apply(tree.ringsum[,10:42],1,sum)
tree.hel=tree.ringsum
tree.hel[,10:42]=tree.hel[,10:42]/tree.rowsums
tree.hel[,10:42]=sqrt(tree.hel[,10:42])

#### Forward selection using regionalized soil variables
# First we have to merge the datasets
tree.hel.soil=merge(tree.hel,soildata)
rownames(tree.hel.soil)=tree.hel.soil[,1]
colnames(tree.hel.soil)
dim(tree.hel.soil)
is.na(tree.hel.soil)
tree.hel.soil=na.exclude(tree.hel.soil)

# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
#hel.rda.fulsoil = rda(tree.hel.soil[,10:42], tree.hel.soil[,43:76])
#anova(hel.rda.fulsoil)
#RsquareAdj(hel.rda.fulsoil)
#global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
#hel.rda.for=forward.sel(tree.hel.soil[,10:42], tree.hel.soil[,52:76], adjR2thresh=global.thresh)
#hel.rda.for
#names(tree.hel.soil)
#hel.rda.selsoil = rda(tree.hel.soil[,10:42] ~ as.matrix(tree.hel.soil[,c(55,54,68)]))
#anova(hel.rda.selsoil)
#RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)
# Cool analysis, but what a terrible plot.  Much better to get the scores and plot them in a real graphing program, or do some fancy graphing in R if that is your thing.

#### Forward selection using PCNM vectors
## First we perform PCNM to get the vectors
#all.geog.dist=vegdist(tree.hel.soil[,8:9],method="euclid")
#all.pcnm1=pcnm(all.geog.dist)
#all.pcnmvecs=all.pcnm1$vectors
## If you can get the PCNM library to run, you could do the line below to select only significantly structured PCNMs
#tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs[,which(all.pcnm1$Moran_I[,2]<0.05)])
#tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
#names(tree.hel.all.pcnm)
## Now run step 1, global analysis
#hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,10:44], tree.hel.all.pcnm[,80:117])
#anova(hel.rda.all.fulpcnm)
#RsquareAdj(hel.rda.all.fulpcnm)
#global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared

###Not signficant, so we don't run the following code:
##now run step 2, the selection part
#hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,10:44], tree.hel.all.pcnm[,80:117], adjR2thresh=global.thresh)
#hel.rda.all.forpcnm
## Let's store a model using the selected variables only
#names(tree.hel.all.pcnm)
#hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,10:44] ~ as.matrix(tree.hel.all.pcnm[,c(82,99)]))
#anova(hel.rda.all.selpcnm)
#RsquareAdj(hel.rda.all.selpcnm)

## Now that we have these parsimonious models, NOW WE CAN PERFORM VARIANCE PARTITIONING!  Hurray!
#varpart(tree.hel.all.pcnm[,10:44], tree.hel.all.pcnm[,c(55,54,68)], tree.hel.all.pcnm[,c(52,54,56,57,58)])
## The variance explained by soil variables and spatial structure clearly overlap. 
## It also appears that there is more independent spatial structure in the trees than independent structuring by soil.
## Are the soil variables significant after accounting for spatial structure, and vice versa?
#hel.rda.all.partpcnm=rda(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(52,54,56,57,58)])
#anova(hel.rda.all.partpcnm)
## 
#hel.rda.all.partsoil=rda(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,c(52,54,56,57,58)], tree.hel.all.pcnm[,c(77:81,83:86)])
#anova(hel.rda.all.partsoil)
## 


###### Some Bonus tree code
# Here is a quick analysis to test for the edge effect on tree communities when Ecosystem type is the explanatory factor, just like we did before with soil data
# First the core plots
is.na(tree.hel)
tree.hel[is.na(tree.hel)] <- 0
tree.hel.core=tree.hel[tree.hel$EP=="Core",]
tree.hel.core.rda.ecosys = rda(tree.hel.core[,10:42] ~ tree.hel.core[,3])
anova(tree.hel.core.rda.ecosys)
core.arsq=RsquareAdj(tree.hel.core.rda.ecosys)
# Then the core+edge plots
tree.hel.edge=tree.hel[tree.hel$EP=="Core"|tree.hel$EP=="Edge",]
tree.hel.edge.rda.ecosys = rda(tree.hel.edge[,10:42] ~ tree.hel.edge[,3])
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
	datcore.rda.ecosys = rda(dat.core.rand[,10:42] ~ dat.core.rand[,3])
	core.arsq.rand=RsquareAdj(datcore.rda.ecosys)
	res.vec[i]=core.arsq.rand$adj.r.squared/edge.arsq$adj.r.squared
}
tiles=rank(rbind(res.vec,edge.eff))
1-tiles[nperm+1]/(nperm+1)
# How well is each species fitted by these models?  Find out below.
goodness(tree.hel.core.rda.ecosys)
goodness(tree.hel.edge.rda.ecosys)

##Same idea as above, but Jaccard distance analysis by PCoA/distance-based RDA
jac.core.dist=vegdist(tree.hel.core[,10:42],method="jaccard",binary=TRUE)
jac.core.pcoa=cmdscale(jac.core.dist,eig=TRUE,k=41,add=TRUE)
jac.core.rda.ecosys = rda(jac.core.pcoa$points ~ tree.hel.core[,3])
anova(jac.core.rda.ecosys)
jac.core.arsq=RsquareAdj(jac.core.rda.ecosys)
jac.edge.dist=vegdist(tree.hel.edge[,10:42],method="jaccard",binary=TRUE)
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

goodness(jac.core.rda.ecosys)
goodness(jac.edge.rda.ecosys)

##Recruitment
treedata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Research/Data/Trees/General transition graph files/Species Recruitment Matrix.csv")
soildata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Grad School/SpatialStats/soil data May 08 for tree plots.csv")
names(treedata)
dim(treedata)
names(soildata)
dim(soildata)
# The tree data was collected in three rings for each plot
# We are going to exclude stream plots (no soil data) and then sum the three rings
treedata=treedata[treedata$Ecosystem!="S",]

treedata[is.na(treedata)] <- 0

dim(treedata)
tree.ringsum=treedata

##for(i in 1:83) {
##	tree.ringsum[i,10:42]=treedata[i,10:42]+treedata[i+83,10:42]+treedata[i+166,10:42]  }

#
rownames(tree.ringsum)=tree.ringsum[,1]
# Let's look at a rank abbundance plot for fun
tree.colsums=apply(tree.ringsum[,10:35],2,sum)
plot(rank(tree.colsums),tree.colsums)

# Perform Hellinger transformation on tree species data so that RDA is based on Hellinger distance
tree.rowsums=apply(tree.ringsum[,10:35],1,sum)
tree.hel=tree.ringsum
tree.hel[,10:35]=tree.hel[,10:35]/tree.rowsums
tree.hel[,10:35]=sqrt(tree.hel[,10:35])


#### Forward selection using regionalized soil variables
# First we have to merge the datasets
rownames(soildata)=soildata$Plot
rownames(tree.hel)=treedata$Plot
tree.hel.soil=merge(tree.hel,soildata)
rownames(tree.hel.soil)=tree.hel.soil$Plot
is.na(tree.hel.soil)
tree.hel.soil=na.exclude(tree.hel.soil)
dim(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil = rda(tree.hel.soil[,10:35], tree.hel.soil[,38:71])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
##Not significant. So we do not run the following code:
## Here is the forward selection command and then constructing a model with selected variables.
#hel.rda.for=forward.sel(tree.hel.soil[,10:42], tree.hel.soil[,52:76], adjR2thresh=global.thresh)
#hel.rda.for
#names(tree.hel.soil)
#hel.rda.selsoil = rda(tree.hel.soil[,10:42] ~ as.matrix(tree.hel.soil[,c(55,54,68)]))
#anova(hel.rda.selsoil)
#RsquareAdj(hel.rda.selsoil)
#plot(hel.rda.selsoil)
## Cool analysis, but what a terrible plot.  Much better to get the scores and plot them in a real graphing program, or do some fancy graphing in R if that is your thing.

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,8:9],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
# If you can get the PCNM library to run, you could do the line below to select only significantly structured PCNMs
#tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs[,which(all.pcnm1$Moran_I[,2]<0.05)])
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,10:35], tree.hel.all.pcnm[,72:107])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared

###Not signficant, so we don't run the following code:
##now run step 2, the selection part
#hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,10:44], tree.hel.all.pcnm[,80:117], adjR2thresh=global.thresh)
#hel.rda.all.forpcnm
## Let's store a model using the selected variables only
#names(tree.hel.all.pcnm)
#hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,10:44] ~ as.matrix(tree.hel.all.pcnm[,c(82,99)]))
#anova(hel.rda.all.selpcnm)
#RsquareAdj(hel.rda.all.selpcnm)

# Now that we have these parsimonious models, NOW WE CAN PERFORM VARIANCE PARTITIONING!  Hurray!
varpart(tree.hel.all.pcnm[,10:44], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(52,54,56,57,58)])
# The variance explained by soil variables and spatial structure clearly overlap. 
# It also appears that there is more independent spatial structure in the trees than independent structuring by soil.
# Are the soil variables significant after accounting for spatial structure, and vice versa?
# Let's see.
hel.rda.all.partpcnm=rda(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,c(77:81,83:86)], tree.hel.all.pcnm[,c(52,54,56,57,58)])
anova(hel.rda.all.partpcnm)
# Yes.
hel.rda.all.partsoil=rda(tree.hel.all.pcnm[,10:42], tree.hel.all.pcnm[,c(52,54,56,57,58)], tree.hel.all.pcnm[,c(77:81,83:86)])
anova(hel.rda.all.partsoil)
# Also yes.  But only by p-value; that is really not much variance explained anymore.


###### Some Bonus tree code
# Here is a quick analysis to test for the edge effect on tree communities when Ecosystem type is the explanatory factor, just like we did before with soil data
# First the core plots
tree.hel.core=tree.hel[tree.hel$EP=="Core",]
tree.hel.core.rda.ecosys = rda(tree.hel.core[,10:35] ~ tree.hel.core[,3])
anova(tree.hel.core.rda.ecosys)
core.arsq=RsquareAdj(tree.hel.core.rda.ecosys)
# Then the core+edge plots
tree.hel.edge=tree.hel[tree.hel$EP=="Core"|tree.hel$EP=="Edge",]
tree.hel.edge.rda.ecosys = rda(tree.hel.edge[,10:35] ~ tree.hel.edge[,3])
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
	datcore.rda.ecosys = rda(dat.core.rand[,10:35] ~ dat.core.rand[,3])
	core.arsq.rand=RsquareAdj(datcore.rda.ecosys)
	res.vec[i]=core.arsq.rand$adj.r.squared/edge.arsq$adj.r.squared
}
tiles=rank(rbind(res.vec,edge.eff))
1-tiles[nperm+1]/(nperm+1)
# How well is each species fitted by these models?  Find out below.
goodness(tree.hel.core.rda.ecosys)
goodness(tree.hel.edge.rda.ecosys)

##Same idea as above, but Jaccard distance analysis by PCoA/distance-based RDA
jac.core.dist=vegdist(tree.hel.core[,10:35],method="jaccard",binary=TRUE)
jac.core.pcoa=cmdscale(jac.core.dist,eig=TRUE,k=34,add=TRUE)
jac.core.rda.ecosys = rda(jac.core.pcoa$points ~ tree.hel.core[,3])
anova(jac.core.rda.ecosys)
jac.core.arsq=RsquareAdj(jac.core.rda.ecosys)
jac.edge.dist=vegdist(tree.hel.edge[,10:42],method="jaccard",binary=TRUE)
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

goodness(jac.core.rda.ecosys)
goodness(jac.edge.rda.ecosys)




###Species Growth
treedata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Research/Data/Trees/Species Growth Matrix.csv")
soildata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Grad School/SpatialStats/soil data May 08 for tree plots.csv")
names(treedata)
dim(treedata)
names(soildata)
dim(soildata)
# The tree data was collected in three rings for each plot
# We are going to exclude stream plots (no soil data) and then sum the three rings - there's no need to sum the rings here.
treedata=treedata[treedata$ecosys!="S",]
dim(treedata)
treedata[is.na(treedata)] <- 0
tree.ringsum=treedata

##for(i in 1:83) {
##	tree.ringsum[i,10:42]=treedata[i,10:42]+treedata[i+83,10:42]+treedata[i+166,10:42]  }


#
rownames(tree.ringsum)=tree.ringsum[,1]
# Let's look at a rank abbundance plot for fun
tree.colsums=apply(tree.ringsum[,10:35],2,sum)
plot(rank(tree.colsums),tree.colsums)

# Perform Hellinger transformation on tree species data so that RDA is based on Hellinger distance
tree.rowsums=apply(tree.ringsum[,10:43],1,sum)
tree.hel=tree.ringsum
tree.hel[,10:43]=tree.hel[,10:43]/tree.rowsums
tree.hel[,10:43]=sqrt(tree.hel[,10:43])

#### Forward selection using regionalized soil variables
# First we have to merge the datasets
dim(tree.hel)
dim(soildata)
rownames(soildata)=soildata$Plot
tree.hel.soil=merge(tree.hel,soildata)
dim(tree.hel.soil)
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
names(tree.hel.soil)
dim(tree.hel.soil)
rownames(tree.hel.soil)
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil = rda(tree.hel.soil[,10:43], tree.hel.soil[,45:78])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared
# Here is the forward selection command and then constructing a model with selected variables.
hel.rda.for=forward.sel(tree.hel.soil[,10:43], tree.hel.soil[,45:78], adjR2thresh=global.thresh)
hel.rda.for
names(tree.hel.soil)
hel.rda.selsoil = rda(tree.hel.soil[,10:43] ~ as.matrix(tree.hel.soil[,c(58,59,73,56)]))
anova(hel.rda.selsoil)
RsquareAdj(hel.rda.selsoil)
plot(hel.rda.selsoil)
# Cool analysis, but what a terrible plot.  Much better to get the scores and plot them in a real graphing program, or do some fancy graphing in R if that is your thing.

#### Forward selection using PCNM vectors
# First we perform PCNM to get the vectors
all.geog.dist=vegdist(tree.hel.soil[,8:9],method="euclid")
all.pcnm1=pcnm(all.geog.dist)
all.pcnmvecs=all.pcnm1$vectors
# If you can get the PCNM library to run, you could do the line below to select only significantly structured PCNMs
#tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs[,which(all.pcnm1$Moran_I[,2]<0.05)])
tree.hel.all.pcnm=cbind(tree.hel.soil,all.pcnmvecs)
names(tree.hel.all.pcnm)
# Now run step 1, global analysis
hel.rda.all.fulpcnm=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,79:112])
anova(hel.rda.all.fulpcnm)
RsquareAdj(hel.rda.all.fulpcnm)
global.thresh=RsquareAdj(hel.rda.all.fulpcnm)$adj.r.squared

##now run step 2, the selection part
hel.rda.all.forpcnm=forward.sel(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,79:112], adjR2thresh=global.thresh)
hel.rda.all.forpcnm
## Let's store a model using the selected variables only
names(tree.hel.all.pcnm)
hel.rda.all.selpcnm = rda(tree.hel.all.pcnm[,10:44] ~ as.matrix(tree.hel.all.pcnm[,c(80,82,85,87,86)]))
anova(hel.rda.all.selpcnm)
RsquareAdj(hel.rda.all.selpcnm)

# Now that we have these parsimonious models, NOW WE CAN PERFORM VARIANCE PARTITIONING!  Hurray!
varpart(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(80,82,85,87,86)], tree.hel.all.pcnm[,c(58,59,73,56)])
# The variance explained by soil variables and spatial structure clearly overlap. 
# It also appears that there is more independent spatial structure in the trees than independent structuring by soil.
# Are the soil variables significant after accounting for spatial structure, and vice versa?
# Let's see.
hel.rda.all.partpcnm=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(80,82,85,87,86)], tree.hel.all.pcnm[,c(58,59,73,56)])
anova(hel.rda.all.partpcnm)
# Yes.
hel.rda.all.partsoil=rda(tree.hel.all.pcnm[,10:43], tree.hel.all.pcnm[,c(58,59,73,56)], tree.hel.all.pcnm[,c(80,82,85,87,86)])
anova(hel.rda.all.partsoil)
# Also yes.  But only by p-value; that is really not much variance explained anymore.

for(i in 1:nperm){
	dat.ec.rand=tree.hel.edge
	dat.ec.rand[dat.ec.rand$ecosys=="B",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="B",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="B",]))
	dat.ec.rand[dat.ec.rand$ecosys=="R",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="R",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="R",]))
	dat.ec.rand[dat.ec.rand$ecosys=="U",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="U",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="U",]))
	dat.ec.rand[dat.ec.rand$ecosys=="D",5]=sample(dat.ec.rand[dat.ec.rand$ecosys=="D",5],size=nrow(dat.ec.rand[dat.ec.rand$ecosys=="D",]))
	dat.core.rand=dat.ec.rand[dat.ec.rand$EP=="Core",]
	datcore.rda.ecosys = rda(dat.core.rand[,10:35] ~ dat.core.rand[,3])
	core.arsq.rand=RsquareAdj(datcore.rda.ecosys)
	res.vec[i]=core.arsq.rand$adj.r.squared/edge.arsq$adj.r.squared
}


###### Some Bonus tree code
# Here is a quick analysis to test for the edge effect on tree communities when Ecosystem type is the explanatory factor, just like we did before with soil data
# First the core plots
tree.hel.core=tree.hel[tree.hel$EP=="Core",]
ncol(tree.hel.core)
tree.hel.core=na.omit(tree.hel.core)
tree.hel.core.rda.ecosys = rda(tree.hel.core[,10:43] ~ tree.hel.core[,3])
anova(tree.hel.core.rda.ecosys)
core.arsq=RsquareAdj(tree.hel.core.rda.ecosys)
# Then the core+edge plots
tree.hel.edge=tree.hel[tree.hel$EP=="Core"|tree.hel$EP=="Edge",]
tree.hel.edge=na.omit(tree.hel.edge)
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
	datcore.rda.ecosys = rda(dat.core.rand[,10:44] ~ dat.core.rand[,3])
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
jac.core.pcoa=cmdscale(jac.core.dist,eig=TRUE,k=42,add=TRUE)
jac.core.rda.ecosys = rda(jac.core.pcoa$points ~ tree.hel.core[,3])
anova(jac.core.rda.ecosys)
jac.core.arsq=RsquareAdj(jac.core.rda.ecosys)
jac.edge.dist=vegdist(tree.hel.edge[,10:42],method="jaccard",binary=TRUE)
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

goodness(jac.core.rda.ecosys)
goodness(jac.edge.rda.ecosys)



####Rank specialists by ecosystem####
treedata=read.csv("C:/Users/Colleen/Dropbox/Jennings 2016/tree basal areas.csv")
soildata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Grad School/SpatialStats/soil data May 08 for tree plots.csv")
treedata[is.na(treedata)] <- 0
treedata[,c(1,3)]

tree.ringsum=treedata[1:80,]
#This [1:80,] refers to all plots within the first ring.
names(treedata)
dim(tree.ringsum)
for(i in 1:80) {
    tree.ringsum[i,10:42]=treedata[i,10:42]+treedata[i+83,10:42]+treedata[i+166,10:42]}
rownames(tree.ringsum)=tree.ringsum[,1]
# Let's look at a rank abbundance plot for fun
tree.colsums=apply(tree.ringsum[,10:42],2,sum)
plot(rank(tree.colsums),tree.colsums)
print(tree.colsums)


# Perform Hellinger transformation on tree species data so that RDA is based on Hellinger distance
tree.rowsums=apply(tree.ringsum[,10:42],1,sum)
tree.hel=tree.ringsum
tree.hel[,10:42]=tree.hel[,10:42]/tree.rowsums
tree.hel[,10:42]=sqrt(tree.hel[,10:42])

tree.ringsum=treedata[treedata$ecosys=="R",]
treedata[,c(1,3)]
tree.ringsum[,c(1,3)]
tree.colsums=apply(tree.ringsum[,10:42],2,sum)
plot(rank(tree.colsums),tree.colsums)
print(tree.colsums)
ncol(tree.ringsum)

tree.hel.soil=merge(tree.ringsum,soildata)
dim(tree.hel.soil)
rownames(tree.hel.soil)=tree.hel.soil[,1]
tree.hel.soil=na.exclude(tree.hel.soil)
names(tree.hel.soil)
dim(tree.hel.soil)
rownames(tree.hel.soil)
colnames(tree.hel.soil)
tree.hel.soil$Plot
# Here is the global analysis. Can proceed if significant.  Use adj. R square as additional stopping criterion.
hel.rda.fulsoil = rda(tree.hel.soil[,10:42], tree.hel.soil[,52:76])
anova(hel.rda.fulsoil)
RsquareAdj(hel.rda.fulsoil)
global.thresh=RsquareAdj(hel.rda.fulsoil)$adj.r.squared

apply(tree.hel.soil[,10:42],2,sum)

