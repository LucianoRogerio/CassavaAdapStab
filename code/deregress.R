
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
