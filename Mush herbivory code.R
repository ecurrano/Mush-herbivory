#MUSH HERBIVORY PAPER
#Purpose: Recreate all analyses in Currano & Jacobs manuscript
#R Studio version: 1.0.153
#Script Author: Ellen D. Currano

#######################LOAD REQUIRED PACKAGES###################################
library(tidyverse)
library(vegan)
library(MASS)
library(agricolae)

#######################LOAD REQUIRED DATA FILES#################################
mask      <- read.csv("mask.csv", header=T) #file defines each DT: whether it is herbivory, specialized, gall, mine
ffg_def   <- read.csv("ffgdef.csv", header=T)$x #file that matches DTs to FFGs: Oviposition=0, Galling=3, Hole Feeding=4, Mining =5, Margin Feeding =6,  Piercing and Sucking=7, Skeletonization=1, Surface feeding =2, 
Mush_raw <- read.csv("Mush_full dataset.csv", header=T) #raw Mush herbivory data: every leaf analysed, the host ID, and presence/absence of every DT
raw.hosts <- read.csv("Mush host data.csv", header=T, row.names = 1) #plant host data for all Mush taxa, including potential driving factors and damage metrics 
strat <- read.csv("Mush Strat Level Data.csv", header=T, row.names=1)

################CREATE ADDITIONAL MATRICES FOR ANALYSIS#########################
#The code below 
  #1. subsets the data by stratigraphic level and creates a full raw matrix for each (called Strat_A, etc)
  #2. creates a matrix with the percent of leaves with each FFG at each stratigraphic level (called ffg_comp)
ffg_comp <- matrix(ncol = 7, nrow=6) #create a matrix with the relative abundance of each ffg at each stratigraphic level
colnames(ffg_comp) <- c("Skeletonization", "Surface feeding", "Gall", "Hole Feeding", "Mining", "Margin Feeding", "Piercing")
rownames(ffg_comp) <- c ("StratA","StratB","StratC","StratD","StratE","StratF")

strat <- unique(unlist(Mush_raw$Stratigraphic.Level))
strat <- sort(strat)
for (i in 1:length(strat)){
  data_2 <- subset(Mush_raw, Stratigraphic.Level==strat[i])
  data_2 <- data_2[,6:231]
  data_3 <- data_2[,2:226]
  ffg_data  <- matrix (ncol=max(ffg_def), nrow=nrow(data_3), data=0)
  for (k in 1:max(ffg_def)) {
    ffg_data[,k]<-as.numeric(as.logical(rowSums(data_3[,ffg_def==k])))
    ffg_comp[i,] <- colSums(ffg_data)/nrow(data_3)*100 }
  assign(paste0("Strat_",strat[i]), data_2)
}

#####################CONDUCT ANALYTICAL RAREFACTION#############################
#The analytic rarefaction solution for herbivory was provided by Torsten Wappler and was first described and published in Gunkel & Wappler 2015.
#calculates total, specialized, mine, and gall richness on bulk floras and individual species
#The summary .csv file output includes total, specialized, mine, and gall damage frequency

dataname <- "Strat_A" #Name appended to output files. Choose yourself
cutoff <- 20 #cutoff for resampling data from single species

dataall   <- Strat_A #choose matrix on which do analytical rarefaction will be performed
datatable <- dataall[,2:ncol(dataall)]
species   <- as.factor(dataall[,1])

ffg_data  <- matrix (ncol=max(ffg_def), nrow=nrow(datatable), data=0)
for (i in 1:max(ffg_def)) ffg_data[,i]<-as.numeric(as.logical(rowSums(datatable[,ffg_def==i])))

#Counter Functions

counta <- function (dat, classa) {
  v <- dat[,classa]
  count <- sum(v)
  count}

countb <- function (dat, classa, classb) {
  va <- dat[,classa]
  vb <- dat[,classb]
  va <- 1-va
  vb <- 1-vb
  v  <- va*vb
  v <- 1-v
  count <- sum(v)
  count}

rarefy <- function(data, dname) {
  
  if (!(is.vector(data[,colSums(data)>0]))) data <-data[,colSums(data)>0]
  
  N <- nrow(data) # find the number of individuals
  m <- ncol(data) # find the number of classes
  
  k <- array(0, c(m))
  l <- array(0, c(m,m))
  
  rarefaction <- matrix(data=NA, nrow=N, ncol=3) #This is your final data matrix 
  barefaction <- matrix(data=NA, nrow=N, ncol=3) #Intermediate results
  
  colnames(rarefaction) <- c("E(X)", "Var(X)", "SD")
  
  #Main Calculations
  
  lister <- array (0, c(N,m)) #Start: 1st Moment
  for (i in 1:m) k[i] <- counta (data, i)
  for (j in 1:m) lister[1,j] <- ((N-k[j])/N)
  for (s in 2:N) for (j in 1:m) if (N-s-k[j]>0) lister[s,j] <- lister[s-1,j]*(N-s-k[j])/(N-s) else lister [s,j]=0
  for (s in 1:N) barefaction [s,1] <- sum(lister[s,])
  for (s in 1:N) rarefaction[s,1] <- m-barefaction [s,1] #End: 1st Moment
  
  listerb <- array (0, c(N,m,m)) #Start: 2nd Moment
  for (i in 2:m) for (j in 1:(i-1)) l[i,j] <- countb (data, i, j)
  for (i in 2:m) for (j in 1:(i-1)) listerb[1,i,j] <- ((N-l[i,j])/N)
  for (s in 2:N) for (i in 2:m) for (j in 1:(i-1)) if (N-s-l[i,j]>0) listerb[s,i,j] <- listerb[s-1,i,j]*(N-s-l[i,j])/(N-s) else listerb[s,i,j] <- 0
  for (s in 1:N) barefaction [s,2] <- sum(listerb[s,,])
  for (s in 1:N) rarefaction[s,2] <- (barefaction [s,1] +2*barefaction[s,2]-barefaction[s,1]**2) #End: 2nd Moment
  for (s in 1:N) rarefaction[s,3] <- rarefaction[s,2]**(1/2)
  
  write.table(rarefaction, paste(dname,".csv"), sep=",")
}


rarefy (datatable[,mask$All],paste(dataname,"_All"))
rarefy (datatable[,mask$Spec],paste(dataname,"_Spec"))
rarefy (datatable[,mask$Gall],paste(dataname,"_Gall"))
rarefy (datatable[,mask$Mine],paste(dataname,"_Mine"))
rarefy (ffg_data, paste(dataname,"_FFG"))

specnames<-levels(species)
species<-as.numeric(species)

for (i in 1:max(species)) if (sum(species==i)>=cutoff & sum(datatable[(species==i), mask$All])>0) rarefy (datatable[(species==i),mask$All],paste(dataname, "_", specnames[i], "_All"))
for (i in 1:max(species)) if (sum(species==i)>=cutoff & sum(datatable[(species==i), mask$Spec])>0) rarefy (datatable[(species==i),mask$Spec],paste(dataname, "_", specnames[i], "_Spec"))
for (i in 1:max(species)) if (sum(species==i)>=cutoff & sum(datatable[(species==i), mask$Gall])>0) rarefy (datatable[(species==i),mask$Gall],paste(dataname, "_", specnames[i], "_Gall"))
for (i in 1:max(species)) if (sum(species==i)>=cutoff & sum(datatable[(species==i), mask$Mine])>0) rarefy (datatable[(species==i),mask$Mine],paste(dataname, "_", specnames[i], "_Mine"))
for (i in 1:max(species)) if (sum(species==i)>=cutoff & sum(ffg_data[species==i,])>0) rarefy (ffg_data[(species==i),], paste(dataname, "_", specnames[i],"_FFG"))

specdata<-matrix(nrow=length(species), ncol=max(species), data=0)
for (i in 1:length(species)) specdata[i,(species[i])]<-1
rarefy (specdata, paste(dataname, "_Species"))

overview<-matrix (nrow=length(specnames)+1, ncol=11)
colnames(overview)<-c("Species", "#leaves", "%DMG", "%Spec", "%Gall", "%Mine", "DTs", "SpecDTs", "GDTs", "MDTs", "#FFGs")
overview[1:length(specnames),1]<-specnames
overview[1:length(specnames),2]<-colSums(specdata)

for (i in 1:length(specnames)) overview[i,3]  <-mean(rowSums(datatable[species==i, mask$All])>0)
for (i in 1:length(specnames)) overview[i,4]  <-mean(rowSums(datatable[species==i, mask$Spec])>0)
for (i in 1:length(specnames)) overview[i,5]  <-mean(rowSums(datatable[species==i, mask$Gall])>0)
for (i in 1:length(specnames)) overview[i,6]  <-mean(rowSums(datatable[species==i, mask$Mine])>0)

for (i in 1:length(specnames)) overview[i,7]  <-sum(colSums(datatable[species==i, mask$All])>0)
for (i in 1:length(specnames)) overview[i,8]  <-sum(colSums(datatable[species==i, mask$Spec])>0)
for (i in 1:length(specnames)) overview[i,9]  <-sum(colSums(datatable[species==i, mask$Gall])>0)
for (i in 1:length(specnames)) overview[i,10] <-sum(colSums(datatable[species==i, mask$Mine])>0)

for (i in 1:length(specnames)) overview[i,11] <-sum(colSums(as.matrix(ffg_data[species==i, ]))>0)

overview[length(specnames)+1,1]<-"Total"
overview[length(specnames)+1,2]<-sum(as.numeric(overview[1:length(specnames),2]))
overview[length(specnames)+1,3]<-mean(rowSums(datatable[,mask$All])>0)
overview[length(specnames)+1,4]<-mean(rowSums(datatable[,mask$Spec])>0)
overview[length(specnames)+1,5]<-mean(rowSums(datatable[,mask$Gall])>0)
overview[length(specnames)+1,6]<-mean(rowSums(datatable[,mask$Mine])>0)
overview[length(specnames)+1,7]<-sum(colSums(datatable[,mask$All])>0)
overview[length(specnames)+1,8]<-sum(colSums(datatable[,mask$Spec])>0)
overview[length(specnames)+1,9]<-sum(colSums(datatable[,mask$Gall])>0)
overview[length(specnames)+1,10]<-sum(colSums(datatable[,mask$Mine])>0)
overview[length(specnames)+1,11]<-sum(colSums(ffg_data)>0)

write.table(overview, paste0("Output_",dataname,"_summary.csv"), sep=",")

#################PLOTS & LINEAR MODELS, PLANT HOST DATA#########################
hosts.20 <- subset(raw.hosts, raw.hosts$No_Leaves_Damage_Census>20) #reduce the dataset to only taxa with at least 20 leaves in the damage census 
hosts.20$Trichomes <- factor(hosts.20$Trichomes) #Get rid of character states that do not exist in the subsetted matrix.

#Nonmetric multidimensional scaling ordination Mush plant hosts by FFGs (Figure 3C)
MtaxaFFG <- hosts.20[,32:38]/100
MtaxaFFG.asin <- asin(sqrt(MtaxaFFG))
MtaxaFFG.asin.nms <- metaMDS(MtaxaFFG.asin, distance="bray")
s1 <- scores(MtaxaFFG.asin.nms, display=c("sites"), choices=1)
s2 <- scores(MtaxaFFG.asin.nms, display=c("sites"), choices=2)
v1 <- scores(MtaxaFFG.asin.nms, display=c("species"), choices=1)
v2 <- scores(MtaxaFFG.asin.nms, display=c("species"), choices=2)

plot(MtaxaFFG.asin.nms, type="n")
points(s1,s2, col="black", pch=19, cex=1)
identify(s1,s2, labels=rownames(MtaxaFFG), cex=1)
points(v1,v2, col="black", pch=3, cex=1)
identify(v1,v2, labels=colnames(MtaxaFFG), cex=1,col="red")
title("Mush Species by Functional Feeding Group")

#Plots of damage metrics on individual hosts, looking at leaf mass per area and legumes (Figure 4)
legsplit <- split(hosts.20, hosts.20$Legume)
trisplit <- split(hosts.20, hosts.20$Trichomes)
par(mfrow=c(2,2))
plot(hosts.20$Leaf_Mass_per_Area, hosts.20$Damage_Frequency, type="n", xlab="Leaf mass per area (g/mm2)", ylab="% Leaves Damaged")
points(legsplit$Y$Leaf_Mass_per_Area, legsplit$Y$Damage_Frequency, col="black", pch=16, cex=2.5)
points(legsplit$N$Leaf_Mass_per_Area, legsplit$N$Damage_Frequency, col="black", pch=1, cex=2.5)
points(trisplit$Y$Leaf_Mass_per_Area, trisplit$Y$Damage_Frequency, col="gray", pch=3, cex=1.5)
plot(hosts.20$Leaf_Mass_per_Area, hosts.20$Spec_Damage_Freq, type="n", xlab="Leaf mass per area (g/mm2)", ylab="% Leaves with Spec. Dam.")
points(legsplit$Y$Leaf_Mass_per_Area, legsplit$Y$Spec_Damage_Freq, pch=16, cex=2.5)
points(legsplit$N$Leaf_Mass_per_Area, legsplit$N$Spec_Damage_Freq, pch=1, cex=2.5)
points(trisplit$Y$Leaf_Mass_per_Area, trisplit$Y$Spec_Damage_Freq, col="gray", pch=3, cex=1.5)
plot(hosts.20$Leaf_Mass_per_Area, hosts.20$Total_DTs_at20, type="n", xlab="Leaf mass per area (g/mm2)", ylab="DTs at 20 leaves", ylim = c(0,10))
points(legsplit$Y$Leaf_Mass_per_Area, legsplit$Y$Total_DTs_at20, pch=16, cex=2.5)
points(legsplit$N$Leaf_Mass_per_Area, legsplit$N$Total_DTs_at20, pch=1, cex=2.5)
points(trisplit$Y$Leaf_Mass_per_Area, trisplit$Y$Total_DTs_at20, pch=3, col="gray", cex=1.5)
plot(hosts.20$Leaf_Mass_per_Area, hosts.20$Spec_DTs_at20, type="n", xlab="Leaf mass per area (g/mm2)", ylab="Spec. DTs at 20 leaves")
points(legsplit$Y$Leaf_Mass_per_Area, legsplit$Y$Spec_DTs_at20, pch=16, cex=2.5)
points(legsplit$N$Leaf_Mass_per_Area, legsplit$N$Spec_DTs_at20, pch=1, cex=2.5)
points(trisplit$Y$Leaf_Mass_per_Area, trisplit$Y$Spec_DTs_at20, pch=3, col="gray", cex=1.5)

###Linear models, testing the influence of leaf mass per area, legume status, presence/absence of trichomes, and proportional abundance of the host plant taxon in the flora
###to view model results, use the summary function
###note that the code below is for single drivers. linear models can also be run (and were) considering multiple factors.

#Damage frequency models
m1.legume <- lm(hosts.20$Damage_Frequency~hosts.20$Legume)
m1.trichome <- lm(hosts.20$Damage_Frequency~hosts.20$Trichomes)
m1.LMA <- lm(hosts.20$Damage_Frequency~hosts.20$Leaf_Mass_per_Area)
m1.abund <- lm(hosts.20$Damage_Frequency~hosts.20$Percent_of_flora)

#Specialized damage frequency models
m2.legume <- lm(hosts.20$Spec_Damage_Freq~hosts.20$Legume)
m2.trichome <- lm(hosts.20$Spec_Damage_Freq~hosts.20$Trichomes)
m2.LMA <- lm(hosts.20$Spec_Damage_Freq~hosts.20$Leaf_Mass_per_Area)
m2.abund <- lm(hosts.20$Spec_Damage_Freq~hosts.20$Percent_of_flora)

#Damage diversity models
m3.legume <- lm(hosts.20$Total_DTs_at20~hosts.20$Legume)
m3.trichome <- lm(hosts.20$Total_DTs_at20~hosts.20$Trichomes)
m3.LMA <- lm(hosts.20$Total_DTs_at20~hosts.20$Leaf_Mass_per_Area)
m3.abund <- lm(hosts.20$Total_DTs_at20~hosts.20$Percent_of_flora)

#Specialized damage diversity models
m4.legume <- lm(hosts.20$Spec_DTs_at20~hosts.20$Legume)
m4.trichome <- lm(hosts.20$Spec_DTs_at20~hosts.20$Trichomes)
m4.LMA <- lm(hosts.20$Spec_DTs_at20~hosts.20$Leaf_Mass_per_Area)
m4.abund <- lm(hosts.20$Spec_DTs_at20~hosts.20$Percent_of_flora)

#Boxplots & t-tests of herbivory on Legumes vs. Non-legumes
par(mfrow=c(2,2))
bp.leg1 <- boxplot(hosts.20$Damage_Frequency~hosts.20$Legume, xlab="Legume", ylab = "% Damage", )
t.test(Damage_Frequency~Legume, data = hosts.20)
bp.leg2 <- boxplot(hosts.20$Spec_Damage_Freq~hosts.20$Legume, xlab="Legume", ylab = "% Spec Damage", ylim=c(0,15))
t.test(Spec_Damage_Freq~Legume, data = hosts.20)
bp.leg3 <- boxplot(hosts.20$Total_DTs_at20~hosts.20$Legume, xlab="Legume", ylab = "Total DTs", ylim=c(3,9))
t.test(Total_DTs_20~Legume, data = hosts.20)
bp.leg4 <- boxplot(hosts.20$Spec_DTs_at20~hosts.20$Legume, xlab="Legume", ylab = "Specialized DTs")
t.test(Spec_DTs_20~Legume, data = hosts.20)

############LINEAR MODELS, BULK FLORAS AT EACH STRATIGRAPHIC LEVEL##############
###Linear models, testing the influence of legume abundance, abundance of plants with trichomes, and proportional abundance of the host plant taxon in the flora
###to view model results, use the summary function
###note that the code below is for single drivers. linear models can also be run considering multiple factors.
###Results summarized in Table 2

#Predictor = Proportional abundance of plant hosts with trichomes
m1.strat.tri <- lm(strat$Dam_Freq~strat$Perc_Trichomes)
m2.strat.tri <- lm(strat$Spec_Dam_Freq~strat$Perc_Trichomes)
m3.strat.tri <- lm(strat$DT_at_80~strat$Perc_Trichomes)
m4.strat.tri <- lm(strat$Spec_DTs_80~strat$Perc_Trichomes)

#Predictor = Proportional abundance of plant hosts with N-fixing symbionts
m1.strat.N <- lm(strat$Dam_Freq~strat$Perc_N.fixers)
m2.strat.N <- lm(strat$Spec_Dam_Freq~strat$Perc_N.fixers)
m3.strat.N <- lm(strat$DT_at_80~strat$Perc_N.fixers)
m4.strat.N <- lm(strat$Spec_DTs_80~strat$Perc_N.fixers)

#Predictor = Mean Annual Precipitation
m1.strat.MAP <- lm(strat$Dam_Freq~strat$Mean_Annual_Precipitation)
m2.strat.MAP <- lm(strat$Spec_Dam_Freq~strat$Mean_Annual_Precipitation)
m3.strat.MAP <- lm(strat$DT_at_80~strat$Mean_Annual_Precipitation)
m4.strat.MAP <- lm(strat$Spec_DTs_80~strat$Mean_Annual_Precipitation)

#Predictor = Floral Evenness (Pielou's J)
m1.strat.ev <- lm(strat$Dam_Freq~strat$Evenness)
m2.strat.ev <- lm(strat$Spec_Dam_Freq~strat$Evenness)
m3.strat.ev <- lm(strat$DT_at_80~strat$Evenness)
m4.strat.ev <- lm(strat$Spec_DTs_80~strat$Evenness)

#Predictor = Floral Diversity at 80 leaves
m1.strat.pl <- lm(strat$Dam_Freq~strat$Plant_diversity_80)
m2.strat.pl <- lm(strat$Spec_Dam_Freq~strat$Plant_diversity_80)
m3.strat.pl <- lm(strat$DT_at_80~strat$Plant_diversity_80)
m4.strat.pl <- lm(strat$Spec_DTs_80~strat$Plant_diversity_80)
