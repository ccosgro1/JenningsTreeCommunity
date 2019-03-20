##Cosgrove 20 March 2019
#Seedling Data Exploration

##Load the Data
juvies=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Research/Data/Trees/JuvenileMatrix.csv")
seeds=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Research/Data/Trees/SeedlingMatrix.csv")
soildata=read.csv("C:/Users/Colleen/OneDrive/Documents/School/Grad School/SpatialStats/soil data May 08 for tree plots.csv")
treedata=read.csv("C:/Users/Colleen/Dropbox/Jennings 2016/tree basal areas.csv")

#Define rownames
dim(juvies)
rownames(juvies)=juvies[,1]
rownames(seeds)=seeds[,1]
rownames(seeds)

#Combine ecosystem type and soil variables into the seedling and juvenile datasets
#juvies=merge(juvies,soildata)
#seeds=merge(seeds,soildata)

#Some prelim exploratory stats
colsums(juvies)
juveco=aov()
