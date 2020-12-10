# Homework 8
Weekly assignment about dimension reduction for linear regression

## Exercise 1

On the Virtuale page of the course (under Data Sets) you find the dataset
LBPnumerical.dat. The 773 observations are lower back pain patients. The Variable y is the Roland-Morris
score, a measurement of the severity of back pain after two weeks of treatment on a scale
between 0 and 100. It is of interest to predict this from the 102 x-variables, which all
contain information about the patients before their treatment was started. This dataset
is a subset (with some observations and variables removed because of too many missing
values, and some imputed values for further missing values, so that no missing values are
left in the data) of the dataset documented on

http://ifcs.boku.ac.at/repository/challenge2/

so you can go there if you like to find out more about the variables, which is not strictly
necessary for this exercise.
The main interest here is to construct a good prediction model for the Roland-Morris
score. Find such a model from a comparison of approaches that should involve at least the
lasso, principal component regression, and a comparison by double cross-validation. Also
demonstrate that (or whether) your model is better than a model that predicts y from the
mean only. You can try out more if you want.


## Exercise 2 

For the Covid-data, try to find a model that is substantially better than prediction by mean only using the methods that you have learnt up to now, using a logtransformed y-variable (note that for comparison you also need to assess the mean-only
model using the transformed y). If you think it is useful, you can also transform x-variables.




# Exercise 1

```r
library(leaps)
library(glmnet)
library(pls)
library(knitr)
```



```r
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
```


|                 | regsubsets|  LASSO|    PCR|   mean|
|:----------------|----------:|------:|------:|------:|
|prediction error |     22.103| 21.662| 21.567| 25.526|

```

```



|           |   |   |   |   |   |   |   |   |   |   |
|:----------|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|regsubsets | 23|  8| 14|  6| 29| 27| 24| 14| 20| 17|
|LASSO      |  7|  9| 25| 27| 23| 29| 24|  7| 18| 18|
|PCR        |  8| 13| 15| 29| 33| 31|  9|  8| 18|  8|

As we can see the best model is the one selected by the PCR method in terms of prediction error. A model selected by PCR is also clearly better the predictig y using only the mean. The problem is that we get very unstable results in terms of how many components we should introduce in the model.






# Exercise 2 

```r
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
```



|                 | regsubsets| LASSO|   PCR| mean|
|:----------------|----------:|-----:|-----:|----:|
|prediction error |      1.609| 1.645| 1.598| 1.68|

```
## Warning in kable_markdown(x = structure(c("regsubsets", "LASSO", "PCR", : The
## table should have a header (column names)
```



|           |   |   |   |   |   |   |   |   |   |   |
|:----------|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
|regsubsets |  1|  1|  1|  1| 62| 15| 16|  4|  8|  5|
|LASSO      |  3|  1|  1| 16| 25| 35| 32| 13|  4|  5|
|PCR        |  1|  1|  1| 51| 19| 27|  5|  5|  6|  5|

Again the suggested selection method is PCR, that gives better results in terms of prediction error. The most frequent suggested number of components is 1, but again we get pretty unstable results, expecially because of the fourth, fifth and sixth iteration.
