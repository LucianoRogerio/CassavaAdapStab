# PREPARE RAW DATA FOR STAGE 1 ANALYSIS -----------------------------------------------------------------------
dir("Analises Doutorado/Desregressao")
setwd("Analises Doutorado/Desregressao")
phen<-read.table("Dados_Finais_2020_15_locais.txt", head=T)
str(phen)
table(phen$LOCATION,phen$YEAR)
table(phen$TRIAL,phen$YEAR)
table(phen$CLONE,phen$YEAR)
table(phen$CLONE,phen$TRIAL)

# Histograms for each trait
hist(phen$FRY)
hist(phen$DMC)
hist(phen$StY)
hist(phen$PPD)
hist(phen$BroLS)
hist(phen$BliLS)
hist(phen$WLS)
hist(phen$ANT)
hist(phen$CBB)
hist(phen$Fus.Peel)
hist(phen$Fus.Pulp)
hist(phen$Phy.Peel)
hist(phen$Phy.Pulp)
hist(phen$Scy.Peel)
hist(phen$Scy.Pulp)

###IF NECESSARY TRANSFORM DATA
library(forecast)
lambda = BoxCox.lambda(phen$FRY)
phen[,"tFRY"]<-BoxCox(phen$FRY, lambda)
hist(phen$tFRY)

lambda = BoxCox.lambda(phen$StY)
phen[,"tStY"]<-BoxCox(phen$StY, lambda)
hist(phen$tStY)

lambda = BoxCox.lambda(phen$PPD)
phen[,"tPPD"]<-BoxCox(phen$PPD, lambda)
hist(phen$tPPD)

lambda = BoxCox.lambda(phen$WLS)
phen[,"tWLS"]<-BoxCox(phen$WLS, lambda)
hist(phen$tWLS)

lambda = BoxCox.lambda(phen$ANT)
phen[,"tANT"]<-BoxCox(phen$ANT, lambda)
hist(phen$tANT)

# Trait list (this be frequently useful to iterate an analysis over many traits)
traits<-c("PTRaiz", "PPAerea", "Aplanta", "DMCont")

##########
# Create design variables that explicitly nest factors
# LOC.YEAR = YEAR nested in LOCATION
phen[,"LOC.YEAR"] <- phen$LOCATION
# LOC.YEAR.REP = REP nested within LOC.YEAR
phen[,"LOC.YEAR.REP"]<-paste(phen$LOC.YEAR,phen$REP,sep=":")

# STAGE 1 ANALYSIS: Mixed model using lmer(), no markers (CLONE is i.i.d.), make de-regressed BLUPs --------------------------------------------------
# When the data are very large, it is much faster to do a two-stage analysis
# In the first stage, experimental design factors will be accounted for and we will reduce the
# data down to a single observation per CLONE / genotype
# No markers will be used at this stage
# In order to avoid bias from fitting CLONE as a random effect in both stage 1 and stage 2
# We de-regress the BLUPs for CLONE and then weight the data input into the genomic prediction (stage 2) using information
# from the first stage that accounts for imbalance (e.g. checks with hundreds of observations and entries with very few)

# Function to de-regress BLUPs from the first stage
# This Model: Y ~ (1|CLONE) + LOC.YEAR + (1|LOC.YEAR.REP)
# This is a function we will run on the output of a model fit in the lmer() function [lme4 R package]
# It follows work by Garrick et al. 2009, which has been widely adopted
# The function has to be changed every time you fit a different model / dataset
deregress<-function(model, trait){
  BLUP <- ranef(model, condVar=TRUE)$CLONE
  PEV <- c(attr(BLUP, "postVar"))
  CLONE.var <- c(VarCorr(model)$CLONE) # Extract the variance component for CLONE
  ResidVar = (attr(VarCorr(model),"sc"))^2 # Extract the residual variance component
  LOC.YEAR.REP = c(VarCorr(model)$LOC.YEAR.REP) # Extract the variance component for REP.RANGE
  # You will need a line like the one above for every random effect (not for fixed effects)
  out <- BLUP/(1-(PEV/CLONE.var)) # This is the actual de-regress part (the BLUP for CLONE is divided by (1 - PEV/CLONE.var))
  # PEV is a vector of error variances associated with each individual BLUP...
  # it tells you about how confident you should be in the estimate of an individual CLONE's BLUP value.
  # For now you can just use the code, but if you get a chance, think about what the above formula does.
  r2 <- 1-(PEV/CLONE.var) # Reliability: a confidence value for a BLUP (0 to 1 scale)
  H2 = CLONE.var/(CLONE.var + LOC.YEAR.REP + ResidVar) # An estimate of the broad-sense heritability, must change this formula when you change the model analysis
  wt = (1-H2)/((0.1 + (1-r2)/r2)*H2) # Weights for each de-regressed BLUP
  # There is a paper the determined this crazy formula, Garrick et al. 2009. I wouldn't pay much attn. to it.
  # These weights will be used in the second-step (e.g. cross-validation) to account for what we've done in this step
  # The weights will be fit as error variances associated with each residual value
  VarComps<-as.data.frame(VarCorr(model))
  return(list(Trait=trait, drgBLUP=out, BLUP=BLUP, weights=wt, varcomps=VarComps,H2=H2,Reliability=r2,model=model))
}
save(deregress,file="deregress.Rdata");

# Now to fit the model
# I like to take advantage of multi-core processers wherever possible
# Below I use a parallel loop (foreach) to run multiple traits at once
# You can modify this code to use a simple non-parallel code if you prefer
library(foreach)
library(doParallel)
proctime <- proc.time(); # This is just to store the start time
registerDoParallel(cores=4) # Specify the number of cores (my lab computer has 8; I will use 6 of them)
DRG <- foreach(a=traits, i=icount(), .inorder=TRUE) %dopar% {
  require(lme4)
  model<-lmer(data=phen,formula = get(traits[i]) ~ (1|CLONE) + LOC.YEAR + (1|LOC.YEAR.REP)) # CLONE and REP are random LOC.YEAR is fixed
  drg<-deregress(model,traits[i]) }
save(DRG,file="DRGblups.Rdata");
proc.time() - proctime # Print the end time minus start time = time to complete analysis (about 2.5 mins using 6 cores and Revolutions R)
str(DRG)
drgBLUP.PTRaiz <- as.matrix(DRG[[1]]$drgBLUP)
drgBLUP.PPAerea <- as.matrix(DRG[[2]]$drgBLUP)
drgBLUP.APlanta <- as.matrix(DRG[[3]]$drgBLUP)
drgBLUP.DMCCont <- as.matrix(DRG[[4]]$drgBLUP)

write.table (drgBLUP.PTRaiz, "Fen?tipos Desregressados - PTRaiz.txt", row.names = T, col.names = F, quote = F)
write.table (drgBLUP.PPAerea, "Fen?tipos Desregressados - PPAerea.txt", row.names = T, col.names = F, quote = F)
write.table (drgBLUP.APlanta, "Fen?tipos Desregressados - APlanta.txt", row.names = T, col.names = F, quote = F)
write.table (drgBLUP.DMCCont, "Fen?tipos Desregressados - DMCont.txt", row.names = T, col.names = F, quote = F)






# PROCESS SNP DATA, MAKE KINSHIP MATRIX FOR STAGE 2 ANALYSIS --------------------------------------------------------------
library(tidyverse)
load("Dados GBS Atualizados.RData")
dim(GBSdata)
dnames<-read.table("Clones_GBS_MHPRVG.txt", head=T)

## Select GBS clones
GBSdata2 <- GBSdata[rownames(GBSdata)%in%dnames[,2],]
##Replace rownames
rownames(GBSdata2)
dnames2 <- dnames[dnames[,2]%in%rownames(GBSdata),]
dnames3 <- dnames2[match(rownames(GBSdata2),dnames2[,2]),]

rownames(GBSdata2) <- dnames3[,3]

grep("Blank",rownames(GBSdata),invert=t)
dim(GBSdata)
table(rownames(GBSdata2) %in% phen$CLONE) # All genotyped clones are phenotyped

# Calculate allele frequencies at each locus
freq <- colMeans(GBSdata2, na.rm=T)/2
# Convert raw allele frequencies to MINOR allele frequencies (MAF)
maf <- freq
maf[which(maf > 0.5)] <- 1 - maf[which(maf > 0.5)]
hist(maf,xlab="MAF Distribution",ylab="Frequency")
# Subset to only SNP with MAF>1%
length(which(maf > 1/sqrt(nrow(GBSdata2)))) # 17650... of these SNPs pass
GS.sel.z <- GBSdata2[,which(maf > 1/sqrt(nrow(GBSdata2)))]
# Save SNP file
save(snps,file="snps_MAFfiltered.Rdata")

# Construct kinship matrix with A.mat function (rrBLUP package)
library(rrBLUP)
K<-A.mat(GS.sel.z-1) # A.mat expects SNPs to be coded -1, 0, 1 but these are 0, 1, 2 so we subtract 1 for the correct result
save(K,file="Kinship.Rdata")

# PCA on Kinship using prcomp, scale = FALSE, center = FALSE ------------------------------
# This is a useful step for exploring the population structure in the data
pca<-prcomp(K, scale=FALSE, center=FALSE)
pcs<-data.frame(pca$x)
pcs<-pcs[,1:100]
varexplained<-summary(pca)$importance[,1:30]
varexplained
save(pcs,varexplained,file="PCA_OnKinship.Rdata")

# Plot: PCA -------------
load("PCA_OnKinship.Rdata")
varexplained[,1:4]
library(ggplot2) # Awesome R package for plotting... just my opinion.
ggplot(pcs,aes(x=PC1,y=PC2))+geom_point(size=1.5) + theme_bw() + labs(x="PC1 (46%)",y="PC2 (15%)",title="PCA on Kinship Matrix")

# LOAD AND COMBINE RESULTS FROM STAGE 1 ANALYSIS ---------------------------------------------------------------------------
load("DRGblups.Rdata")
class(DRG) # since I used "foreach" to run each stage 1 analysis in parallel, each trait is in a separate element of a list
# We need to process the object DRG into a data.frame or matrix for further analysis
load("Kinship.Rdata")
drgphen<-data.frame(CLONE = rownames(K),stringsAsFactors=F)
for(i in 1:length(DRG)){
  drg<-data.frame(CLONE = rownames(DRG[[i]]$drgBLUP),stringsAsFactors=F)
  drg[,DRG[[i]]$Trait]<-DRG[[i]]$drgBLUP
  drg[,paste(DRG[[i]]$Trait,"egv",sep=".")]<-DRG[[i]]$BLUP
  drg[,paste(DRG[[i]]$Trait,"wt",sep=".")]<-DRG[[i]]$weights
  drgphen<-merge(drgphen,drg,by="CLONE",all.x=T)
}
rownames(drgphen)<-drgphen$CLONE

# The dataframe made by the loop above has 3 columns for each trait:
# DMC = the de-regressed BLUP
# DMC.egv = the BLUP (not de-regressed); egv stands for "Estimated Genetic Value"
# DMC.wt = the weights

# STAGE 2 ANALYSIS: OBTAIN GENOMIC ESTIMATED BREEDING VALUES (G-BLUP method) ------------------------------------------------
# Will generate predicted breeding values for all individuals with and without phenotypes
# Note how we incorporate the weights, it is the same as fitting 1/sqrt(weight) in the error/residual variance-covariance matrix (R)
# Here I use the R package EMMREML
# I used a regular loop over the 3 traits here, but you could also try using the parallel for loop (foreach)
#traits<-c("DMC","FRY","StY", "PPD", "BroLS", "BliLS", "WLS", "ANT", "CBB")
#traits<-c("FRY", "DMC", "StY")#, "PPD", "BroLS", "BliLS", "WLS", "ANT", "CBB", "FusPeel", "FusPulp", "PhyPeel", "PhyPulp", "ScyPeel", "ScyPulp")
ModelResults<-list()
for(i in 1:length(traits)){
  trait<-traits[i]
  require(EMMREML)
  data1<-drgphen[which(drgphen$CLONE %in% rownames(K)),] # Only genotyped clones can be analyzed
  data1<-data1[!is.na(data1[,trait]),] # Must remove rows where there are missing values for the current trait
  data1$CLONE <- factor(as.character(data1$CLONE),levels=rownames(K))
  Z = as.matrix(model.matrix(~data1$CLONE-1)) # The incidence matrix for random effects; only CLONE here
  X = model.matrix(~1,data=data1) # The incidence matrix for fixed effects
  y = data1[,trait] # The vector of phenotypes
  weight<-sqrt(data1[,paste(trait,"wt",sep=".")])
  Z = Z/weight
  X = X/weight
  y = y/weight
  funout <- emmreml(y=y, X=X, Z=Z, K=K, varuhat=T,PEVuhat=T) # The actual prediction
  funout[["Trait"]]<-trait # Add an element to the output that stores the trait that was analyzed
  ModelResults[[trait]]<-funout # Add the model results to the list "ModelResults"
}
save(ModelResults,file="ModelResults.Rdata")

# Combine results into a dataframe with GEBVs -------------------------------------------------------------------------------------------
load("ModelResults.Rdata")
load("Kinship.Rdata")
gebvs<-data.frame(CLONE=rownames(K))
for(i in 1:length(ModelResults)){
  thistrait<-ModelResults[[i]]
  uhat<-thistrait$uhat; colnames(uhat)[1]<-thistrait$Trait
  gebvs<-merge(gebvs,uhat,by.x="CLONE",by.y="row.names",all=T);
}
write.csv(gebvs,file="GEBVS.csv",row.names=F)
hist(gebvs$DMC)
hist(gebvs$FRY)
hist(gebvs$StY)
hist(gebvs$PPD)
hist(gebvs$BroLS)
hist(gebvs$BliLS)
hist(gebvs$WLS)
hist(gebvs$ANT)
hist(gebvs$CBB)

# FOLD CROSS-VALIDATION -------------------------
# Below are two functions below designed to conduct cross-validation with EMMREML
# The first is the main part, it chooses the folds and sets up a dataset for input to the second function
# The second function will do the prediction requested by the first one
# The first / main function collects the output from the second one, calculates accuracy
# and puts everything together into a useful output format
# The code is pretty complicated only because I set it up so you could fit models with multiple random effects if you want to.
# And also because it checks the dataset for missing and other things to make it as general as possible
# It works only with de-regressed BLUPs arranged in a dataset like the one in the above section (LOAD AND COMBINE RESULTS FROM STAGE 1 ANALYSIS)

# Fold cross-validation function for emmreml multikernel 3.0 (V3: saves loglik + BLUPs for each rep) --------------------------------------------------------------------
FoldCrossValidation.V3.emmreml <- function(dataset, trait, genoID, Klist, nFolds, nRepeats){
  data<-dataset[which(dataset[ ,genoID] %in% rownames(Klist[[1]])),] # Only genotyped indivs in phenotype dataset
  data<-data[!is.na(data[,trait]),] # Remove anyone with a missing value for the trait of interest
  rownames(data)<-data$CLONE
  # Only phenotyped individuals in genotype dataset (loop below)
  for(kernels in 1:length(Klist)){
    K<-Klist[[kernels]]
    K<-K[rownames(K) %in% data[ ,genoID], rownames(K) %in% data[ ,genoID]]
    Klist[[kernels]]<-K }
  nInd <- dim(data)[1] # How many individuals?
  # Create empty dataframes where results will be stored across reps
  accuracies<-data.frame()
  varcomps<-data.frame()
  blups<-data.frame()
  # The rep loop
  for (rep in 1:nRepeats){
    print(paste("Rep ",rep," of ",nRepeats,sep=""))
    folds <- sample(rep(1:nFolds, length.out=nInd)) # Individuals are assigned fold membership
    # Empty dataframes to hold results from each fold within reps
    BLUPSthisRep<-data.frame()
    varcomps.thisrep<-data.frame()
    for (fold in 1:nFolds){
      print(paste("Fold ",fold," of ",nFolds,sep=""))
      # index of individuals in the current fold (validation set)
      indInFold <- which(folds == fold)
      # index of individuals NOT in the current fold (training set)
      indNotInFold <- which(folds != fold)
      # List of clone names in current Validation set
      ValidSet<-data[indInFold,genoID]
      # List of clone names in current Training set
      TrainSet<-data[indNotInFold,genoID]
      # Prediction occurs using EMMREML wrapper function (see below)
      out<-emmremlMultiKernel.CrossVal.mod(data,train=TrainSet,test=ValidSet,trait,Klist=Klist)
      # Calculate AIC
      out[["AIC"]]<-2*length(Klist)-2*out$loglik
      # If there is more than one kernel (kinship matrix)
      if(length(Klist)>1){
        colnames(out$uhatmat)<-paste("K",1:length(Klist),sep="")
        # Store the predictions (BLUPs) for the Validation set this current fold
        BLUPSthisFold<-data.frame(CLONE=rownames(out$uhatmat[ValidSet,]))
        BLUPSthisFold<-cbind(BLUPSthisFold,out$uhatmat[ValidSet,])
        # Combine the BLUPs for this fold's validation set with those from previous folds
        BLUPSthisRep<-rbind(BLUPSthisRep,BLUPSthisFold)
        names(out$weights)<-colnames(out$uhatmat)
        # Also store the variance components estimated this fold
        varcomps.thisfold<-data.frame(Trait=trait,Rep=rep,Fold=fold,Vu=out$Vu, Ve=out$Ve, loglik=out$loglik, AIC=out$AIC)
        varcomps.thisfold<-cbind(varcomps.thisfold,t(out$weights)) }
      # Else if there was just one kinsihp matrix
      else {
        BLUPSthisFold<-data.frame(CLONE=names(out$uhatmat[ValidSet,]))
        BLUPSthisFold<-cbind(BLUPSthisFold,out$uhatmat[ValidSet,])
        colnames(BLUPSthisFold)<-c(genoID,"K1")
        BLUPSthisRep<-rbind(BLUPSthisRep,BLUPSthisFold)
        varcomps.thisfold<-data.frame(Trait=trait,Rep=rep,Fold=fold,Vu=out$Vu, Ve=out$Ve, loglik=out$loglik, AIC=out$AIC)
      }
      varcomps.thisrep<-rbind(varcomps.thisrep,varcomps.thisfold)
    }
    BLUPSthisRep[,"Rep"]<-rep
    # Calc accuracy after predicting each fold to complete the rep
    # Code to calculate accuracy if more than one kernel
    if(length(Klist)>1){
      BLUPSthisRep[,"TotGen"]<-rowSums(BLUPSthisRep[,c(paste("K",1:length(Klist),sep=""))],na.rm=T)
      BLUPSthisRep<-merge(BLUPSthisRep,data[,c("CLONE",paste(trait,".egv",sep=""))],by="CLONE")
      blup.names<-c(paste("K",1:length(Klist),sep=""),"TotGen")
      accuracies.thisrep<-vector()
      for(tacos in 1:length(blup.names)){
        acc = cor(BLUPSthisRep[,blup.names[tacos]],BLUPSthisRep[,paste(trait,".egv",sep="")], use="complete.obs")
        accuracies.thisrep<-c(accuracies.thisrep,acc) }
      names(accuracies.thisrep)<-paste("Acc.",blup.names,sep="")
      AccuracyThisRep<-data.frame(Trait=trait,Rep=rep)
      AccuracyThisRep<-cbind(AccuracyThisRep,t(accuracies.thisrep))
      accuracies<-rbind(accuracies,AccuracyThisRep)
      varcomps<-rbind(varcomps,varcomps.thisrep)
      # Code to calculate accuracy if only one kernel
    } else {
      # The dataframe BLUPSthisRep stores the predictions from this replication
      # That were made on each individual clone when it was in the validation set
      # Merge is with the BLUP from step 1 (*.ebv in drgphenos)
      BLUPSthisRep<-merge(BLUPSthisRep,data[,c("CLONE",paste(trait,".egv",sep=""))],by="CLONE")
      blup.names<-c(paste("K",1:length(Klist),sep=""))
      accuracies.thisrep<-vector()
      for(tacos in 1:length(blup.names)){
        # Calculate accuracy for this rep
        # correlation(SPROUT.gebv,SPROUT.ebv)
        acc = cor(BLUPSthisRep[,blup.names[tacos]],BLUPSthisRep[,paste(trait,".egv",sep="")], use="complete.obs")
        accuracies.thisrep<-c(accuracies.thisrep,acc) }
      names(accuracies.thisrep)<-paste("Acc.",blup.names,sep="")
      AccuracyThisRep<-data.frame(Trait=trait,Rep=rep)
      AccuracyThisRep<-cbind(AccuracyThisRep,t(accuracies.thisrep))
      accuracies<-rbind(accuracies,AccuracyThisRep)
      varcomps<-rbind(varcomps,varcomps.thisrep)
    }
    blups<-rbind(blups,BLUPSthisRep)
  }
  return(list(accuracies=accuracies,varcomps=varcomps, blups=blups))
}
# Wrapper for multikernel emmreml 3.0 (V5: Used in Cross-validation function) ------------------------------
emmremlMultiKernel.CrossVal.mod<-function(data,train,test,trait,Klist){
  require(EMMREML)
  # Only data from phenotyped individuals can be considered
  data1 <- data[which(data$CLONE %in% union(train,test)),c("CLONE",trait,paste(trait,"egv",sep="."),paste(trait,"wt",sep="."))]
  # No missing value is allowed for the trait
  data1 <- data1[!is.na(data1[,trait]),]
  rownames(data1) <- data1$CLONE
  # Training and test lists are supplied by the cross-validation function(above)
  # Check that they don't contain unphenotyped individuals and remove them
  train = train[train %in% data1$CLONE]
  test = test[test %in% data1$CLONE]
  # For each kernel, make sure only phenotyped inividuals are allowed
  for(kernels in 1:length(Klist)){ Klist[[kernels]]<-Klist[[kernels]][union(train,test),union(train,test)] }
  datatrain <- data1[train,]
  datatest <- data1[test,]
  # Even though datatrain only contains the training set, we set the levels of
  # the factor CLONE equal to the rownames of the kinship matrix
  # There should be more individuals in the kinship matrix than in the
  # phenotypes because kinship matrix contains vaidation and training sets
  datatrain$CLONE <- factor(as.character(datatrain$CLONE),levels=rownames(Klist[[1]]))
  # We do this because it ensures the model.matrix() function below creates the
  # correct design/incidence matrices for prediction unphenotyped individuals
  Z = as.matrix(model.matrix(~datatrain$CLONE-1))
  y = datatrain[,trait]
  X = model.matrix(~1,data=datatrain)
  # The next four lines are where the weights from Step 1 come in
  # We divided the Z, X and y by the sqrt(weights)
  weight<-sqrt(data1[train,paste(trait,"wt",sep=".")])
  Z = Z/weight
  X = X/weight
  y = y/weight
  Zlist = list()
  for(kernels in 1:length(Klist)){ Zlist[[paste("Z",kernels,sep="")]] <- Z }
  # The function decides whether to use the single or multikernel version of EMMREML
  # based on how many kinship matrices are in Klist
  if(length(Klist) == 1){
    funout <- emmreml(y=y, X=X, Z=Zlist[[1]], K=Klist[[1]]) } else {
      funout <- emmremlMultiKernel(y=y, X=X, Z=Zlist, K=Klist, varbetahat=F,varuhat=F, PEVuhat=F, test=F) }
  uhatmat<-matrix(funout$uhat, ncol=length(Klist))
  rownames(uhatmat)<-rownames(Klist[[1]])
  totvals<-rowSums(uhatmat,na.rm=T)
  names(totvals)<-rownames(Klist[[1]])
  funout[["Train"]]<-train
  funout[["Test"]]<-test
  funout[["Trait"]]<-trait
  funout[["uhatmat"]]<-uhatmat
  funout[["TotVals"]]<-totvals
  funout[["Dataset"]]<-data1
  return(funout)
}

# Load and prepare data for cross-validation --------------------------------------
load("DRGblups.Rdata")
class(DRG) # since I used "foreach" to run each stage 1 analysis in parallel, each trait is in a separate element of a list
# We need to process the object DRG into a data.frame or matrix for further analysis
load("Kinship.Rdata")
drgphen<-data.frame(CLONE = rownames(K),stringsAsFactors=F)
for(i in 1:length(DRG)){
  drg<-data.frame(CLONE = rownames(DRG[[i]]$drgBLUP),stringsAsFactors=F)
  drg[,DRG[[i]]$Trait]<-DRG[[i]]$drgBLUP
  drg[,paste(DRG[[i]]$Trait,"egv",sep=".")]<-DRG[[i]]$BLUP
  drg[,paste(DRG[[i]]$Trait,"wt",sep=".")]<-DRG[[i]]$weights
  drgphen<-merge(drgphen,drg,by="CLONE",all.x=T)
}
rownames(drgphen)<-drgphen$CLONE

# Cross-validate
CrossVal<-list()
for(l in 1:length(traits)){
  require(EMMREML); crossval<-FoldCrossValidation.V3.emmreml(drgphen,traits[l],"CLONE",list(K),5,3)
  CrossVal[[l]]<-crossval
}
save(CrossVal,file="CrossVal.Rdata")

# Combine and plot results
setwd("~/Google Drive/EMBRAPA/TAMU_GS_Workshop/")
load("CrossVal.Rdata")
library(plyr)
accuracies<-ldply(CrossVal,function(x){ out<-x$accuracies; return(out) })

library(ggplot2)
ggplot(accuracies,aes(x=Trait,y=Acc.K1,fill=Trait)) + geom_boxplot() + theme_bw()

