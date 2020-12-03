library(leaps)
library(glmnet)
library(pls)
library(knitr)

###########regsubsets
LBP <- read.table("LBPnumerical.dat.txt",header=TRUE)

varnames <- names(LBP[2:103])
response <- "y"

LBP<-LBP[sample(nrow(LBP)),]
indexes <- cut(seq(1,nrow(LBP)),breaks=10,labels=FALSE)
fold=list()
for(i in 1:10){
  fold[[i]] <- which(indexes==i,arr.ind=TRUE)
  
}

nmodels <- 50
P=102
response <- 1
cvdata <- LBP
names(cvdata)[response] <- "y"
full.model=as.formula(paste("y","~", paste(varnames,collapse=" + ")))
n=773
yhat <- besti_sub <- numeric(0)
fmbesti_sub <- mbesti_sub <- bestmodel <- modforcvi <- list()
yhati <- sqrlossi <- matrix(NA,nrow=nmodels,ncol=n)
foldnumber=10
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  yhati <- matrix(NA,nrow=nmodels,ncol=n)
  modforcv <- list()
  for (j in (1:foldnumber)[-i]){
    cat("Inner fold ",j,"\n")
    bestsubi <- regsubsets(full.model,data=cvdata[-c(fold[[i]],fold[[j]]),],
                           nvmax=nmodels,method="forward") # leaving folds i and j out
    sbesti_sub <- summary(bestsubi)
    for (k in 1:nmodels){
      modforcv[[k]] <- as.formula(paste("y~", paste(varnames[sbesti_sub$which[k,2:(P+1)]],
                                                    collapse=" + "))) # extract best model
      fmi <- lm(modforcv[[k]], data=cvdata[-c(fold[[i]],fold[[j]]),]) # fit it
      yhati[k,fold[[j]]] <- predict(fmi,cvdata[fold[[j]],])
      # predict fold j for finding best k
    } # end for k (models)
  } # end for j (inner loop) outer loop still running
  for (k in 1:nmodels){
    sqrlossi[k,i] <- sqrt(mean((yhati[k,]-cvdata$y)^2,na.rm=TRUE))
  }
  besti_sub[i] <- which.min(sqrlossi[,i]) # Best model chosen without fold i
  bestmodel[[i]] <- regsubsets(full.model,data=cvdata[-fold[[i]],],
                               nvmax=besti_sub[i],method="forward") # run forward selection on data without fold i
  sbesti_sub <- summary(bestmodel[[i]])
  modforcvi[[i]] <- as.formula(paste("y~", paste(varnames[sbesti_sub$which[besti_sub[[i]],2:(P+1)]],
                                                 collapse=" + ")))
  # Extract best model as found in inner loop
  fmbesti_sub[[i]] <- lm(modforcvi[[i]], data=cvdata[-fold[[i]],]) # Fit this
  yhat[fold[[i]]] <- predict(fmbesti_sub[[i]],cvdata[fold[[i]],])
  # Predict fold i data with best model selected without fold i.
}
sqrlossbest <- sqrt(mean((yhat-cvdata$y)^2))
besti_sub
#########LASSO
# Use existing folds for 10-fold for outer loop, to make result
# comparable with forward:
foldidorig <- rep(NA,n)
for(i in 1:foldnumber){
  foldidorig[fold[[i]]] <- i # cv.glmnet allows fold indicating vector as input.
}

# The outer loop needs to be the same, but not the inner loop!
response <- 1
yhat <- bestp <- bestlambda<-besti_LASSO <- bestcvm <- numeric(0)
lmodel <- list()
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  lmodel[[i]] <- glmnet(x=as.matrix(cvdata[varnames])[-fold[[i]],],y=cvdata$y[-fold[[i]]])
  # Need this in oder to unify lambdas for all fitted models.
  lambdai <- lmodel[[i]]$lambda
  foldidx <- foldidorig # Preparation of removing fold i
  foldidx[foldidx>=i] <- foldidx[foldidx>=i]-1 # Update fold numbers if larger than i
  foldidi <- foldidx[-fold[[i]]] # remove fold i from foldid vector
  foldilasso <- cv.glmnet(x=as.matrix(cvdata[varnames])[-fold[[i]],],
                          y=cvdata$y[-fold[[i]]],foldid=foldidi,lambda=lambdai)
  besti_LASSO[i] <- which.min(foldilasso$cvm)
  bestp[i] <- foldilasso$nzero[besti_LASSO[i]] # Number of nonzero variables
  bestlambda[i] <- foldilasso$lambda[besti_LASSO[i]] # lambda of best model
  bestcvm[i] <- foldilasso$cvm[besti_LASSO[i]] # cross validation loss of best model
  yhat[fold[[i]]] <- predict(foldilasso$glmnet.fit,
                             as.matrix(cvdata[varnames][fold[[i]],]),s=bestlambda[i])
}

sqrlosslasso <- sqrt(mean((yhat-cvdata$y)^2))
besti_LASSO
##############PCR
response <- 1
yhat<- besti_PCR<-numeric(0)
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  pcrLBP <- pcr(full.model, ncomp=50,data=cvdata[-fold[[i]],],scale=TRUE,validation="CV")
  besti_PCR[i] <- which.min(as.data.frame(RMSEP(pcrLBP)$val)[1,])-1

  yhat[fold[[i]]] <-predict(pcrLBP,cvdata[fold[[i]],] , ncomp = besti_PCR[i])
}

sqrlossPCR <- sqrt(mean((yhat-cvdata$y)^2))
besti_PCR

#########mean
yhat=numeric(0)
for (i in 1:foldnumber){
  model=lm(y~1,data=cvdata[-fold[[i]],],)
  yhat[fold[[i]]] <-predict(model,cvdata[fold[[i]],] , ncomp = besti_PCR[i])
  
}
sqrlossmean<- sqrt(mean((yhat-cvdata$y)^2))

########## Exercise 4 ############################################################
cvdata=read.table("covidweeklygrowths.txt",header = T)
cvdata <- cvdata[,5:170]
varnames <- names(cvdata)[-1]
response <- names(cvdata)[1]

cvdata<-cvdata[sample(nrow(cvdata)),]
indexes <- cut(seq(1,nrow(cvdata)),breaks=10,labels=FALSE)
fold=list()
for(i in 1:10){
  fold[[i]] <- which(indexes==i,arr.ind=TRUE)
  
}

nmodels <- 50
P=165
response <- 1
names(cvdata)[response] <- "y"
cvdata$y=log(cvdata$y)
full.model=as.formula(paste("y","~", paste(varnames,collapse=" + ")))
n=175
yhat <- besti_sub <- numeric(0)
fmbesti_sub <- mbesti_sub <- bestmodel <- modforcvi <- list()
yhati <- sqrlossi <- matrix(NA,nrow=nmodels,ncol=n)
foldnumber=10
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  yhati <- matrix(NA,nrow=nmodels,ncol=n)
  modforcv <- list()
  for (j in (1:foldnumber)[-i]){
    cat("Inner fold ",j,"\n")
    bestsubi <- regsubsets(full.model,data=cvdata[-c(fold[[i]],fold[[j]]),],
                           nvmax=nmodels,method="forward") # leaving folds i and j out
    sbesti_sub <- summary(bestsubi)
    for (k in 1:nmodels){
      modforcv[[k]] <- as.formula(paste("y~", paste(varnames[sbesti_sub$which[k,2:(P+1)]],
                                                    collapse=" + "))) # extract best model
      fmi <- lm(modforcv[[k]], data=cvdata[-c(fold[[i]],fold[[j]]),]) # fit it
      yhati[k,fold[[j]]] <- predict(fmi,cvdata[fold[[j]],])
      # predict fold j for finding best k
    } # end for k (models)
  } # end for j (inner loop) outer loop still running
  for (k in 1:nmodels)
    sqrlossi[k,i] <- sqrt(mean((yhati[k,]-cvdata$y)^2,na.rm=TRUE))
  besti_sub[i] <- which.min(sqrlossi[,i]) # Best model chosen without fold i
  bestmodel[[i]] <- regsubsets(full.model,data=cvdata[-fold[[i]],],
                               nvmax=besti_sub[i],method="forward") # run forward selection on data without fold i
  sbesti_sub <- summary(bestmodel[[i]])
  modforcvi[[i]] <- as.formula(paste("y~", paste(varnames[sbesti_sub$which[besti_sub[[i]],2:(P+1)]],
                                                 collapse=" + ")))
  # Extract best model as found in inner loop
  fmbesti_sub[[i]] <- lm(modforcvi[[i]], data=cvdata[-fold[[i]],]) # Fit this
  yhat[fold[[i]]] <- predict(fmbesti_sub[[i]],cvdata[fold[[i]],])
  # Predict fold i data with best model selected without fold i.
}
sqrlossbest4 <- sqrt(mean((yhat-cvdata$y)^2))
besti_sub
#########LASSO
# Use existing folds for 10-fold for outer loop, to make result
# comparable with forward:
foldidorig <- rep(NA,n)
for(i in 1:foldnumber){
  foldidorig[fold[[i]]] <- i # cv.glmnet allows fold indicating vector as input.
}

# The outer loop needs to be the same, but not the inner loop!
response <- 1
yhat <- bestp <- bestlambda<-besti_LASSO <- bestcvm <- numeric(0)
lmodel <- list()
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  lmodel[[i]] <- glmnet(x=as.matrix(cvdata[varnames])[-fold[[i]],],y=cvdata$y[-fold[[i]]])
  # Need this in oder to unify lambdas for all fitted models.
  lambdai <- lmodel[[i]]$lambda
  foldidx <- foldidorig # Preparation of removing fold i
  foldidx[foldidx>=i] <- foldidx[foldidx>=i]-1 # Update fold numbers if larger than i
  foldidi <- foldidx[-fold[[i]]] # remove fold i from foldid vector
  foldilasso <- cv.glmnet(x=as.matrix(cvdata[varnames])[-fold[[i]],],
                          y=cvdata$y[-fold[[i]]],foldid=foldidi,lambda=lambdai)
  besti_LASSO[i] <- which.min(foldilasso$cvm)
  bestp[i] <- foldilasso$nzero[besti_LASSO[i]] # Number of nonzero variables
  bestlambda[i] <- foldilasso$lambda[besti_LASSO[i]] # lambda of best model
  bestcvm[i] <- foldilasso$cvm[besti_LASSO[i]] # cross validation loss of best model
  yhat[fold[[i]]] <- predict(foldilasso$glmnet.fit,
                             as.matrix(cvdata[varnames][fold[[i]],]),s=bestlambda[i])
}

sqrlosslasso4 <- sqrt(mean((yhat-cvdata$y)^2))
besti_LASSO
##############PCR
response <- 1
yhat<- besti_PCR<-numeric(0)
for (i in 1:foldnumber){
  cat("Outer fold ",i,"\n")
  pcrLBP <- pcr(full.model, ncomp=50,data=cvdata[-fold[[i]],],scale=TRUE,validation="CV")
  besti_PCR[i] <- which.min(as.data.frame(RMSEP(pcrLBP)$val)[1,])-1
  
  yhat[fold[[i]]] <-predict(pcrLBP,cvdata[fold[[i]],] , ncomp = besti_PCR[i])
}

sqrlossPCR4 <- sqrt(mean((yhat-cvdata$y)^2))
besti_PCR

#########mean
yhat=numeric(0)
for (i in 1:foldnumber){
  model=lm(y~1,data=cvdata[-fold[[i]],],)
  yhat[fold[[i]]] <-predict(model,cvdata[fold[[i]],] , ncomp = besti_PCR[i])
  
}
sqrlossmean4<- sqrt(mean((yhat-cvdata$y)^2))
