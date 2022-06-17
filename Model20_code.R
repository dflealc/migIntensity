### MST & Migration Intensity
### code written by Diego F. Leal, edited by Nicolas L. Harder
### Last update: 5/13/2020
### Code and outputs are in UMass server: C:\Users\dleal\Documents\Net_Weights

## clear all
rm(list=ls())

##using R siena
library(RSiena)



#setwd("~/Dropbox/Leal_Harder/Abstract for new paper")
setwd("~/JEMS Replication folder")

#load("Abstract workspace.RData")    #Model 2
load("Abstract workspace .33 - cutoff 2 - binarized.RData") #Model 2


setwd("~/Output")

##################################
#Create the node level attributes#
##################################

#Create nodeset
nodeset <- as.vector(node_attributes$Master)

#GDP node Var for time change#
node_attributes$aici_e90s  <-((GDPnode$`1991`)+(GDPnode$`1995`)) / 2 # aici in early 1990s
node_attributes$aici_l90s  <-((GDPnode$`1996`)+(GDPnode$`2000`)) / 2 # aici in late 1990s
node_attributes$aici_e00s  <-((GDPnode$`2001`)+(GDPnode$`2005`)) / 2 # aici in early 2000s
node_attributes$aici_l00s  <-((GDPnode$`2006`)+(GDPnode$`2010`)) / 2 # aici in late 2000s 
node_attributes$aici_e10s  <-((GDPnode$`2011`)+(GDPnode$`2015`)) / 2 # aici in early 2010s 

#Pop node Var for time change#
node_attributes$pop_e90s <-((Popfullnode$`1991`)+(Popfullnode$`1995`)) / 2 # Population in early 1990s
node_attributes$pop_l90s  <-((Popfullnode$`1996`)+(Popfullnode$`2000`)) / 2 # Population in late 1990s
node_attributes$pop_e00s  <-((Popfullnode$`2001`)+(Popfullnode$`2005`)) / 2 # Population in early 2000s
node_attributes$pop_l00s  <-((Popfullnode$`2006`)+(Popfullnode$`2010`)) / 2 # Population in late 2000s 
node_attributes$pop_e10s  <-((Popfullnode$`2011`)+(Popfullnode$`2015`)) / 2 # Population in early 2010s

node_attributes$dep_e90s <- ((OldYoungfull$`1990`)+(OldYoungfull$`1995`)) / 2 # Old/Young dependency in early 1990s
node_attributes$dep_l90s <- ((OldYoungfull$`1995`)+(OldYoungfull$`2000`)) / 2 # Old/Young dependency in late 1990s
node_attributes$dep_e00s <- ((OldYoungfull$`2000`)+(OldYoungfull$`2005`)) / 2 # Old/Young dependency in early 2000s
node_attributes$dep_l00s <- ((OldYoungfull$`2005`)+(OldYoungfull$`2010`)) / 2 # Old/Young dependency in late 2000s
node_attributes$dep_e10s <- ((OldYoungfull$`2010`)+(OldYoungfull$`2015`)) / 2 # Old/Young dependency in early 2010s

#GDP CoVar for time change#
GDP_e90s  <-((absGDPlist[[1]])+(absGDPlist[[2]])) / 2 # aici in early 1990s
GDP_l90s  <-((absGDPlist[[3]])+(absGDPlist[[4]])) / 2 # aici in late 1990s
GDP_e00s  <-((absGDPlist[[5]])+(absGDPlist[[6]])) / 2 # aici in early 2000s
GDP_l00s  <-((absGDPlist[[7]])+(absGDPlist[[8]])) / 2 # aici in late 2000s 
GDP_e10s  <-((absGDPlist[[9]])+(absGDPlist[[10]])) / 2 # aici in early 2010s 

#Pop CoVar for time change#
POP_e90s  <-((absPOPlist[[1]])+(absPOPlist[[2]])) / 2 # Population in early 1990s
POP_l90s  <-((absPOPlist[[3]])+(absPOPlist[[4]])) / 2 # Population in late 1990s
POP_e00s  <-((absPOPlist[[5]])+(absPOPlist[[6]])) / 2 # Population in early 2000s
POP_l00s  <-((absPOPlist[[7]])+(absPOPlist[[7]])) / 2 # Population in late 2000s 
POP_e10s  <-((absPOPlist[[9]])+(absPOPlist[[10]])) / 2 # Population in early 2010s 

#Dependency CoVar for time change#
DEP_e90s <- ((absDEPlist[[1]])+(absDEPlist[[2]])) / 2 # Old/Young dependency in early 1990s
DEP_l90s <- ((absDEPlist[[2]])+(absDEPlist[[3]])) / 2 # Old/Young dependency in late 1990s
DEP_e00s <- ((absDEPlist[[3]])+(absDEPlist[[4]])) / 2 # Old/Young dependency in early 2000s
DEP_l00s <- ((absDEPlist[[4]])+(absDEPlist[[5]])) / 2 # Old/Young dependency in late 2000s
DEP_e10s <- ((absDEPlist[[5]])+(absDEPlist[[6]])) / 2 # Old/Young dependency in early 2010s


#### SIENA MODEL #####

## Create a generic RSIENA model object, call it interaction model
interactionsModel <- sienaModelCreate(fn=simstats0c,           #use the Robbins-Monro estimation algorithm, do not change. For technical details see: http://www.stats.ox.ac.uk/~snijders/siena/Siena_algorithms.pdf
                                      projname="migration",    #give a name to the project. The output file from the simulation will be called projname.out (e.g., migration.out)
                                      nsub=3,                  #number of subphases in the estimation algorith
                                      n3=3000,                  #number of simulated data sets, determines the precision of the standard errors
                                      seed=5)                  #random seed generator

##labels of the node Attributes to be used in the simulation
vtxDat <- node_attributes[,c("Region_Num", 
                             "Coups",
                             "pop_e90s",
                             #"pop_l90s",
                             #"pop_e00s",
                             #"pop_l00s",
                             #"pop_e10s", 
                             "aici_e90s",
                             #"aici_l90s",
                             #"aici_e00s",
                             #"aici_l00s",
                             #"aici_e10s",
                             "dep_e90s", 
                             #"dep_l90s",
                             #"dep_e00s",
                             #"dep_l00s",
                             #"dep_e10s",
                             #"coup_ave",
                             "Access"
                             )]

## Create three arrays to represent different flows, one for North-North flows, one for South-South flows, and one for North-South flows.
SmallArray<-array(cbind(SmallMatrices[[1]],SmallMatrices[[2]]),
                  dim=c(dim(SmallMatrices[[1]]),2))

MediumArray<-array(cbind(MediumMatrices[[1]],MediumMatrices[[2]]),
                  dim=c(dim(MediumMatrices[[1]]),2))

LargeArray<-array(cbind(LargeMatrices[[1]],LargeMatrices[[2]]),
                   dim=c(dim(LargeMatrices[[1]]),2))


## Create RSiena network objects based on the two arrways above
interactionNetS <- sienaDependent(SmallArray)
interactionNetM <- sienaDependent(MediumArray)
interactionNetL <- sienaDependent(LargeArray)

# Define as time-constant actor covars
region     <- coCovar(as.vector(vtxDat[,"Region_Num"]))
#Coup       <- coCovar(as.vector(vtxDat[,"coup_ave"]))

#is.vector(vtxDat$Access)
access     <- coCovar(as.vector(vtxDat[,"Access"]))

POPNODE     <- coCovar(as.vector(vtxDat[,"pop_e90s"]))

yy<-vtxDat[,"aici_e90s"]
colnames(yy)<-"aici_e90s"

GDPNODE     <- coCovar(as.vector(yy[,"aici_e90s"]))

#Define time-constant Dyad covars
comlang <- coDyadCovar(CEPII.comlang)
contig  <- coDyadCovar(CEPII.contig)

#Define time-varying actor covariates
DEPNODE     <- coCovar(as.vector(vtxDat[,"dep_e90s"]))


zz<-as.matrix(IGOshared[[3]])

igo     <- coDyadCovar(zz)


#crate a siena object with all three actor arrays and all covars and codyadvars
friendshipDat <- sienaDataCreate(interactionNetS,interactionNetM,
                                 interactionNetL,
                                 region,DEPNODE,GDPNODE,POPNODE,comlang,contig,igo)

#Decalring migration flows matrices as disjoint (i.e., mutually exclusive)
sienaDataConstraint(friendshipDat, interactionNetS,interactionNetM, 
                    type="disjoint", value = TRUE)
sienaDataConstraint(friendshipDat, interactionNetM,interactionNetS,
                    type="disjoint", value = TRUE)
sienaDataConstraint(friendshipDat, interactionNetL,interactionNetS,
                    type="disjoint", value = TRUE)
sienaDataConstraint(friendshipDat, interactionNetS,interactionNetL,
                    type="disjoint", value = TRUE)
sienaDataConstraint(friendshipDat, interactionNetL,interactionNetM,
                    type="disjoint", value = TRUE)
sienaDataConstraint(friendshipDat, interactionNetM,interactionNetL,
                    type="disjoint", value = TRUE)

##Create effects object to specify the model
friendshipEff <- getEffects(friendshipDat)

##list of all effects
effNames <- friendshipEff$functionName
effSNames <- friendshipEff$shortName

#setwd("~/Dropbox/Leal_Harder/Model Setup")
write.csv(effNames,file="effectsDescriptionmigrationflows.csv")
write.csv(effSNames,file="effectsShortNamemigrationflows.csv")

#Selection of effects to be used in RSiena Model#


inEffs <- c( 
  #dyadic effects
  "interactionNetS: Number of reciprocated ties",
  "interactionNetM: Number of reciprocated ties",
  "interactionNetL: Number of reciprocated ties",
  #tradic effcts
  "interactionNetS: Number of transitive triplets",
  "interactionNetM: Number of transitive triplets",
  "interactionNetL: Number of transitive triplets",
  "interactionNetS: 3-cycles",   
  "interactionNetM: 3-cycles",
  "interactionNetL: 3-cycles",
  "interactionNetS: Number of transitive recipr. triplets",
  "interactionNetM: Number of transitive recipr. triplets",
  "interactionNetL: Number of transitive recipr. triplets",
  #monadic (degree-related) effects
  "interactionNetS: Sum of indegrees x sqrt(indegree)", 
  "interactionNetM: Sum of indegrees x sqrt(indegree)",
  "interactionNetL: Sum of indegrees x sqrt(indegree)",
  "interactionNetS: Sum of indegrees x sqrt(outdegree)",
  "interactionNetM: Sum of indegrees x sqrt(outdegree)", 
  "interactionNetL: Sum of indegrees x sqrt(outdegree)",
  "interactionNetS: Sum of outdegrees^(1.5)",
  "interactionNetM: Sum of outdegrees^(1.5)",
  "interactionNetL: Sum of outdegrees^(1.5)",
  #Exogenous effects: MST linkages
  "interactionNetS: Sum of ties x comlang",
  "interactionNetM: Sum of ties x comlang",
  "interactionNetL: Sum of ties x comlang",
  "interactionNetS: Sum of ties x contig",
  "interactionNetM: Sum of ties x contig",
  "interactionNetL: Sum of ties x contig",
  "interactionNetS: Same values on region",
  "interactionNetM: Same values on region",
  "interactionNetL: Same values on region",
  "interactionNetS: Sum of ties x igo",
  "interactionNetM: Sum of ties x igo",
  "interactionNetL: Sum of ties x igo",
  "interactionNetS: Sum outdegrees x GDPNODE",
  "interactionNetM: Sum outdegrees x GDPNODE",
  "interactionNetL: Sum outdegrees x GDPNODE",
  "interactionNetS: Sum indegrees x GDPNODE",
  "interactionNetM: Sum indegrees x GDPNODE",
  "interactionNetL: Sum indegrees x GDPNODE",
  #exogenous controls
  "interactionNetS: Sum outdegrees x DEPNODE",
  "interactionNetM: Sum outdegrees x DEPNODE",
  "interactionNetL: Sum outdegrees x DEPNODE",
  "interactionNetS: Sum indegrees x DEPNODE",
  "interactionNetM: Sum indegrees x DEPNODE",
  "interactionNetL: Sum indegrees x DEPNODE",
  "interactionNetS: Sum outdegrees x POPNODE",
  "interactionNetM: Sum outdegrees x POPNODE",
  "interactionNetL: Sum outdegrees x POPNODE",
  "interactionNetS: Sum indegrees x POPNODE",
  "interactionNetM: Sum indegrees x POPNODE",
  "interactionNetL: Sum indegrees x POPNODE"
)

effs <- match(inEffs,effNames)

friendshipEff$include[effs] <- T

friendshipEff

est_siena <- siena07(interactionsModel,
                     data=friendshipDat,
                     effects=friendshipEff,
                     batch=F,
                     returnDeps=T,
                     useCluster=T,
                     nbrNodes=35) #35

setwd("~/Net_Weights/Model20")

save(est_siena, file=paste("Model20_iteration",1,"_conv_",round(est_siena$tconv.max,3),".RData",sep=""))

#this loop will save every iteration of the model as a separatee (numerated) workspace

tm <- est_siena$tconv.max
print(tm[1])

for (i in 2:30)  # repeat this loop 49 times,for a total of 50 iterations
{
  if (tm[1] > 0.17) # if overall max convergence ration is > 0.1, estimate the model one more time
  {
    est_siena <- siena07(interactionsModel,
                         data=friendshipDat,
                         effects=friendshipEff,
                         batch=F,
                         returnDeps=T,
                         useCluster=T,
                         prevAns = est_siena,   #run the model based on the most recent estimates
                         nbrNodes=35
    )
  }
  save(est_siena, file=paste("Model20_iteration",i,"_conv_",round(est_siena$tconv.max,3),".RData",sep="")) # save an image with current estimates
  tm <- est_siena$tconv.max              ## retrieve new version of overall max convergence ratio
  print(i)
  print(tm[1])
  if (tm[1] < 0.17) {break}                 ## if the overall max convergence ratio is < 0.1,the model has converged. Then stop
}


#sienaTimeTest(est_siena)   ## test for time heterogeneity
est_siena         ## summarize the model


dev.off()

#GOF indegree S

pdf(paste("Model20_","conv_",round(est_siena$tconv.max,3),"_IndegreeDistribution","interactionNetS",".pdf",sep=""),width=15,height = 15)
gofa.idegS <- sienaGOF(est_siena, IndegreeDistribution, verbose=TRUE, join=TRUE,varName="interactionNetS")
summary(gofa.idegS)
plot(gofa.idegS)
dev.off()

#GOF indegree M

pdf(paste("Model20_","conv_",round(est_siena$tconv.max,3),"_IndegreeDistribution","interactionNetM",".pdf",sep=""),width=15,height = 15)
gofa.idegM <- sienaGOF(est_siena, IndegreeDistribution, verbose=TRUE, join=TRUE,varName="interactionNetM")
summary(gofa.idegM)
plot(gofa.idegM)
dev.off()

#GOF indegree L

pdf(paste("Model20_","conv_",round(est_siena$tconv.max,3),"_IndegreeDistribution","interactionNetL",".pdf",sep=""),width=15,height = 15)
gofa.idegL <- sienaGOF(est_siena, IndegreeDistribution, verbose=TRUE, join=TRUE,varName="interactionNetL")
summary(gofa.idegL)
plot(gofa.idegL)
dev.off()

#GOF outdegree S

pdf(paste("Model20_","conv_",round(est_siena$tconv.max,3),"_OutdegreeDistribution","interactionNetS",".pdf",sep=""),width=15,height = 15)
gofa.odegS <- sienaGOF(est_siena, OutdegreeDistribution, verbose=TRUE, join=TRUE,varName="interactionNetS")
summary(gofa.odegS)
plot(gofa.odegS)
dev.off()

#GOF outdegree M

pdf(paste("Model20_","conv_",round(est_siena$tconv.max,3),"_OutdegreeDistribution","interactionNetM",".pdf",sep=""),width=15,height = 15)
gofa.odegM <- sienaGOF(est_siena, OutdegreeDistribution, verbose=TRUE, join=TRUE,varName="interactionNetM")
summary(gofa.odegM)
plot(gofa.odegM)
dev.off()

#GOF outdegree L

pdf(paste("Model20_","conv_",round(est_siena$tconv.max,3),"_OutdegreeDistribution","interactionNetL",".pdf",sep=""),width=15,height = 15)
gofa.odegL <- sienaGOF(est_siena, OutdegreeDistribution, verbose=TRUE, join=TRUE,varName="interactionNetL")
summary(gofa.odegL)
plot(gofa.odegL)
dev.off()

save.image(paste("final_Model20_",round(est_siena$tconv.max,3),"_.RData",sep=""))
