######################################

##### RS Project 2020       ######      DFA by area using FIM data
##### Rec fishing PAPER     ######
##### FWC FIM data

######################################



#Load files ######################
# here we will load the MRIP angler effort data for each time series (combo of state, fishing mode, and fishing area)


dat = readRDS("Data outputs/MRIP4DFA.rds")
head(dat)



###################
## Plot          ##
###################

par(mfrow=c(1,1))
plot(dat$Year, dat$`ALABAMA Federal Charter`, type="l", col="red", lwd=2, ylim=c(0,1000000))
lines(dat$Year, dat$`ALABAMA Federal Private`, type="l", col="orange", lwd=2)
lines(dat$Year, dat$`ALABAMA Inland Charter`, type="l", col="green", lwd=2)
lines(dat$Year, dat$`ALABAMA Inland Private`, type="l", col="blue", lwd=2)
lines(dat$Year, dat$`ALABAMA Inland Shore`, type="l", col="purple", lwd=2)
lines(dat$Year, dat$`ALABAMA State Charter`, type="l", col="black", lwd=2)
lines(dat$Year, dat$`ALABAMA State Private`, type="l", col="gray", lwd=2)




###########################
## Moving on to the DFA  ##
###########################

# Libraries

library(MARSS)
library(xtable)
library(bbmle)
library(beepr)


##Get rid of the year column in the dataset
dat1 <- dat[,2:ncol(dat)]

# transpose data so that each aggregate group is a row, and each time step is a column
dat.T <- t(dat1)
N.ts <- nrow(dat.T); N.ts
TT <- ncol(dat.T); TT

# Z score the data
dat.z <- zscore(dat.T)
trial <- zscore(dat.z)

checks <- as.vector(trial[6,])
mean(checks, na.rm=TRUE)
var(checks, na.rm=TRUE)

## Plot the z scored data
regions = rownames(dat.z)
years = seq(1981,2021,1)
par(mfcol=c(4,8), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in regions){
  plot(dat.z[i,]~years,xlab="",ylab="Angler effort", bty="L", pch=16, col="blue", type="b")
  axis(1,1:nrow(dat))
  title(i)
}



###################
## DFA Code      ##  
## Model Fitting ##
###################

## Set up the controls
cntl.list = list(minit=200, maxit=100000, allow.degen=FALSE, conv.test.slope.tol = 0.1, trace=1, safe=TRUE)


## Note you control the number of common trends you want to fit in the model statement by controling the m = part.
## m = 1 (one common trend), m = 2 (two common trends), etc
## Options for the observation error matrix (R =) are 'diagonal and equal', 'diagonal and unequal', 'equalvarcov', and 'unstructured'


#We will fit the models in the following order. For now we do not have covariates, we will look at those later.


# 1. No covariates
#     1.1 One trend
#           1.1.1 diagonal and equal
#           1.1.2 diagonal and unequal
#           1.1.3 equalvarcov
#           1.1.4 unconstrained
#     1.2 Two trends
#           1.2.1 diagonal and equal
#           1.2.2 diagonal and unequal
#           1.2.3 equalvarcov
#           1.2.4 unconstrained
#     1.3 Three trends
#           1.3.1 diagonal and equal
#           1.3.2 diagonal and unequal
#           1.3.3 equalvarcov
#           1.3.4 unconstrained
# 2. Plus covariate x
#     2.1 One trend
#           2.1.1 diagonal and equal
#           2.1.2 diagonal and unequal
#           2.1.3 equalvarcov
#           2.1.4 unconstrained
#     2.2 Two trends
#           2.2.1 diagonal and equal
#           2.2.2 diagonal and unequal
#           2.2.3 equalvarcov
#           2.2.4 unconstrained
#     2.3 Three trends
#           2.3.1 diagonal and equal
#           2.3.2 diagonal and unequal
#           2.3.3 equalvarcov
#           2.3.4 unconstrained
#           ...






###
##
#    1. NO COVARIATES
##
###


###
### 1.1. One Trend
###



###  1.1.1

dfa.1.1.1 = MARSS(dat.z, model=list(m=1, R="diagonal and equal"), z.score=TRUE, form="dfa", control=cntl.list, silent = 2); beep()  

#Converged!
summary(dfa.1.1.1)
MARSSparamCIs(dfa.1.1.1)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.1.1, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.1.1$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.6198877 -  mean fit values should be close to 0, less than 0.6. 

##Plot of common trend
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.1.1$states, type = 'l' )





###  1.1.2

dfa.1.1.2 = MARSS(dat.z, model=list(m=1, R="diagonal and unequal"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#Converged!
summary(dfa.1.1.2)
MARSSparamCIs(dfa.1.1.2)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.1.2, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.1.2$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.6256117

##Plot of common trend
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.1.2$states, type = 'l' )





###  1.1.3

dfa.1.1.3 = MARSS(dat.z, model=list(m=1, R="equalvarcov"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  
#Converged!
summary(dfa.1.1.3)
MARSSparamCIs(dfa.1.1.3)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.1.3, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.1.3$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.6220803 

##Plot of common trend
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.1.3$states, type = 'l' )





###  1.1.4

dfa.1.1.4 = MARSS(dat.z, model=list(m=1, R="unconstrained"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#didn't work well
summary(dfa.1.1.4)
MARSSparamCIs(dfa.1.1.4)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.1.4, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.1.4$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.7405027

##Plot of common trend
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.1.4$states, type = 'l' )






###
### 1.2. Two Trends
###



###  1.2.1

dfa.1.2.1 = MARSS(dat.z, model=list(m=2, R="diagonal and equal"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#Converged!
summary(dfa.1.2.1)
MARSSparamCIs(dfa.1.2.1)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.2.1, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.2.1$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.4941556

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.2.1$states[1,], type = 'l', ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.2.1$states[2,], type = 'l', col = 'blue' )





###  1.2.2

dfa.1.2.2 = MARSS(dat.z, model=list(m=2, R="diagonal and unequal"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#abstol only
summary(dfa.1.2.2)
MARSSparamCIs(dfa.1.2.2)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.2.2, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.2.2$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.5169711

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.2.2$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.2.2$states[2,], type = 'l', col = 'blue' )





###  1.2.3

dfa.1.2.3 = MARSS(dat.z, model=list(m=2, R="equalvarcov"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  
#Converged!
summary(dfa.1.2.3)
MARSSparamCIs(dfa.1.2.3)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.2.3, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.2.3$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.4934001

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.2.3$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.2.3$states[2,], type = 'l', col = 'blue' )





###  1.2.4

dfa.1.2.4 = MARSS(dat.z, model=list(m=2, R="unconstrained"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#didn't really work
summary(dfa.1.2.4)
MARSSparamCIs(dfa.1.2.4)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.2.4, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.2.4$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.6758041 

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.2.4$states[1,], type = 'l' , ylim=c(-8,4))
lines(1:nrow(dat), dfa.1.2.4$states[2,], type = 'l', col = 'blue' )








###
### 1.3. Three Trends
###



###  1.3.1

dfa.1.3.1 = MARSS(dat.z, model=list(m=3, R="diagonal and equal"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#converged!
summary(dfa.1.3.1)
MARSSparamCIs(dfa.1.3.1)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.3.1, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.3.1$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.4104119

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.3.1$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.3.1$states[2,], type = 'l', col = 'blue' )
lines(1:nrow(dat), dfa.1.3.1$states[3,], type = 'l', col = 'green' )





###  1.3.2

dfa.1.3.2 = MARSS(dat.z, model=list(m=3, R="diagonal and unequal"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#converged!
summary(dfa.1.3.2)
MARSSparamCIs(dfa.1.3.2)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.3.2, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.3.2$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.4530755

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.3.2$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.3.2$states[2,], type = 'l', col = 'blue' )
lines(1:nrow(dat), dfa.1.3.2$states[3,], type = 'l', col = 'green' )





###  1.3.3

dfa.1.3.3 = MARSS(dat.z, model=list(m=3, R="equalvarcov"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  
#Converged!
summary(dfa.1.3.3)
MARSSparamCIs(dfa.1.3.3)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.3.3, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.3.3$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.4162422

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.3.3$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.3.3$states[2,], type = 'l', col = 'blue' )
lines(1:nrow(dat), dfa.1.3.3$states[3,], type = 'l', col = 'green' )





###  1.3.4

dfa.1.3.4 = MARSS(dat.z, model=list(m=3, R="unconstrained"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#unstable
summary(dfa.1.3.4)
MARSSparamCIs(dfa.1.3.4)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.3.4, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.3.4$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.4299 

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.3.4$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.3.4$states[2,], type = 'l', col = 'blue' )
lines(1:nrow(dat), dfa.1.3.4$states[3,], type = 'l', col = 'green' )






###
### 1.4. Four Trends
###




###  1.4.1

dfa.1.4.1 = MARSS(dat.z, model=list(m=4, R="diagonal and equal"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#converged!
summary(dfa.1.4.1)
MARSSparamCIs(dfa.1.4.1)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.4.1, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.4.1$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.3103533

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.4.1$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.4.1$states[2,], type = 'l', col = 'blue' )
lines(1:nrow(dat), dfa.1.4.1$states[3,], type = 'l', col = 'green' )
lines(1:nrow(dat), dfa.1.4.1$states[4,], type = 'l', col = 'red' )



###  1.4.2

dfa.1.4.2 = MARSS(dat.z, model=list(m=4, R="diagonal and unequal"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#abstol only
summary(dfa.1.4.2)
MARSSparamCIs(dfa.1.4.2)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.4.2, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.4.2$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.3823471

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.4.2$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.4.2$states[2,], type = 'l', col = 'blue' )
lines(1:nrow(dat), dfa.1.4.2$states[3,], type = 'l', col = 'green' )
lines(1:nrow(dat), dfa.1.4.2$states[4,], type = 'l', col = 'red' )



###  1.4.3

dfa.1.4.3 = MARSS(dat.z, model=list(m=4, R="equalvarcov"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  
#abstol only
summary(dfa.1.4.3)
MARSSparamCIs(dfa.1.4.3)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.4.3, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.4.3$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.3121077

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.4.3$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.4.3$states[2,], type = 'l', col = 'blue' )
lines(1:nrow(dat), dfa.1.4.3$states[3,], type = 'l', col = 'green' )
lines(1:nrow(dat), dfa.1.4.3$states[4,], type = 'l', col = 'red' )




###  1.4.4

dfa.1.4.4 = MARSS(dat.z, model=list(m=4, R="unconstrained"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#didn't work
summary(dfa.1.4.4)
MARSSparamCIs(dfa.1.4.4)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.4.4, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.4.4$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.4558999

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.4.4$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.4.4$states[2,], type = 'l', col = 'blue' )
lines(1:nrow(dat), dfa.1.4.4$states[3,], type = 'l', col = 'green' )
lines(1:nrow(dat), dfa.1.4.4$states[4,], type = 'l', col = 'red' )






###
### 1.5. Five Trends
###




###  1.5.1

dfa.1.5.1 = MARSS(dat.z, model=list(m=5, R="diagonal and equal"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#converged!
summary(dfa.1.5.1)
MARSSparamCIs(dfa.1.5.1)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.5.1, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.5.1$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.2393809

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.5.1$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.5.1$states[2,], type = 'l', col = 'blue' )
lines(1:nrow(dat), dfa.1.5.1$states[3,], type = 'l', col = 'green' )
lines(1:nrow(dat), dfa.1.5.1$states[4,], type = 'l', col = 'red' )


###  1.5.2

dfa.1.5.2 = MARSS(dat.z, model=list(m=5, R="diagonal and unequal"), z.score=TRUE, form="dfa", control=cntl.list, silent=2); beep()  

#abstol only
summary(dfa.1.5.2)
MARSSparamCIs(dfa.1.5.2)

##This gives predicted values for each aggregate group
par.mat = coef(dfa.1.5.2, type = "matrix") 
fit.mod = par.mat$Z  %*%  dfa.1.5.2$states  

##Fit ratios
sumResids = rowSums((dat.z - fit.mod)^2, na.rm=TRUE)
sumObserved = rowSums(dat.z^2, na.rm=TRUE)
FitRatio = sumResids / sumObserved
mean(FitRatio) #0.381521

##Plot of common trends
par(mfrow = c(1,1))
plot(1:nrow(dat), dfa.1.5.2$states[1,], type = 'l' , ylim=c(-5,4))
lines(1:nrow(dat), dfa.1.5.2$states[2,], type = 'l', col = 'blue' )
lines(1:nrow(dat), dfa.1.5.2$states[3,], type = 'l', col = 'green' )
lines(1:nrow(dat), dfa.1.5.2$states[4,], type = 'l', col = 'red' )


###
##
# Save the models so we don't have to re-run
##
###

saveRDS(dfa.1.1.1, file = "DFA models/dfa.1.1.1") #good
saveRDS(dfa.1.1.2, file = "DFA models/dfa.1.1.2") #good
saveRDS(dfa.1.1.3, file = "DFA models/dfa.1.1.3") #good
saveRDS(dfa.1.1.4, file = "DFA models/dfa.1.1.4") #no convergence
saveRDS(dfa.1.2.1, file = "DFA models/dfa.1.2.1") #good
saveRDS(dfa.1.2.2, file = "DFA models/dfa.1.2.2") #abstol only
saveRDS(dfa.1.2.3, file = "DFA models/dfa.1.2.3") #good
saveRDS(dfa.1.2.4, file = "DFA models/dfa.1.2.4") #no convergence
saveRDS(dfa.1.3.1, file = "DFA models/dfa.1.3.1") #good
saveRDS(dfa.1.3.2, file = "DFA models/dfa.1.3.2") #good
saveRDS(dfa.1.3.3, file = "DFA models/dfa.1.3.3") #good
saveRDS(dfa.1.3.4, file = "DFA models/dfa.1.3.4") #no convergence
saveRDS(dfa.1.4.1, file = "DFA models/dfa.1.4.1") #abstol only
saveRDS(dfa.1.4.2, file = "DFA models/dfa.1.4.2") #abstol only
saveRDS(dfa.1.4.3, file = "DFA models/dfa.1.4.3") #abstol only
saveRDS(dfa.1.4.4, file = "DFA models/dfa.1.4.4") #no convergence
saveRDS(dfa.1.5.1, file = "DFA models/dfa.1.5.1") #good
saveRDS(dfa.1.5.2, file = "DFA models/dfa.1.5.2") #abstol only



###
##
# Load the saved models
##
###


dfa.1.1.1 = readRDS("DFA models/dfa.1.1.1")
dfa.1.1.2 = readRDS("DFA models/dfa.1.1.2")
dfa.1.1.3 = readRDS("DFA models/dfa.1.1.3")
#dfa.1.1.4 = readRDS("DFA models/dfa.1.1.4")
dfa.1.2.1 = readRDS("DFA models/dfa.1.2.1")
dfa.1.2.2 = readRDS("DFA models/dfa.1.2.2")
dfa.1.2.3 = readRDS("DFA models/dfa.1.2.3")
#dfa.1.2.4 = readRDS("DFA models/dfa.1.2.4")
dfa.1.3.1 = readRDS("DFA models/dfa.1.3.1")
dfa.1.3.2 = readRDS("DFA models/dfa.1.3.2")
dfa.1.3.3 = readRDS("DFA models/dfa.1.3.3")
#dfa.1.3.4 = readRDS("DFA models/dfa.1.3.4")
dfa.1.4.1 = readRDS("DFA models/dfa.1.4.1")
dfa.1.4.2 = readRDS("DFA models/dfa.1.4.2")
dfa.1.4.3 = readRDS("DFA models/dfa.1.4.3")
#dfa.1.4.4 = readRDS("DFA models/dfa.1.4.4")
dfa.1.5.1 = readRDS("DFA models/dfa.1.5.1")
dfa.1.5.2 = readRDS("DFA models/dfa.1.5.2")




###
##
#  BEST MODEL
##
###


#We will use AICc to determine which is the best model

#AIC.1 = c(dfa.1.1.1$AICc, dfa.1.1.2$AICc, dfa.1.1.3$AICc, dfa.1.1.4$AICc, dfa.1.2.1$AICc, dfa.1.2.2$AICc, dfa.1.2.3$AICc, dfa.1.2.4$AICc, dfa.1.3.1$AICc, dfa.1.3.2$AICc, dfa.1.3.3$AICc, dfa.1.3.4$AICc, dfa.1.4.1$AICc, dfa.1.4.2$AICc, dfa.1.4.3$AICc, dfa.1.4.4$AICc, dfa.1.5.1$AICc, dfa.1.5.2$AICc)

AIC.1 = c(dfa.1.1.1$AICc, dfa.1.1.2$AICc, dfa.1.1.3$AICc, dfa.1.2.1$AICc, dfa.1.2.2$AICc, dfa.1.2.3$AICc, dfa.1.3.1$AICc, dfa.1.3.2$AICc, dfa.1.3.3$AICc, dfa.1.4.1$AICc, dfa.1.4.2$AICc, dfa.1.4.3$AICc, dfa.1.5.1$AICc, dfa.1.5.2$AICc)

delAIC <- AIC.1 - min(AIC.1)
relLik <- exp(-0.5 * delAIC)
aicweight <- relLik/sum(relLik)
#And this leads to our model weights table:

aic.table <- data.frame(AICc = AIC.1, delAIC = delAIC, relLik = relLik, weight = aicweight)
rownames(aic.table) <- c("1 CT, d&e", "1 CT d&u", "1 CT, ev", "2 CT, d&e", "2 CT d&u", "2 CT, ev", "3 CT, d&e", "3 CT d&u", "3 CT, ev", "4 CT, d&e", "4 CT d&u", "4 CT, ev","5 CT, d&e","5 CT, d&u")
#Here the table is printed using round() to limit the number of digits shown.

round(aic.table, digits = 3)

#              AICc  delAIC relLik weight
# 1 CT, d&e 424.309   0.000  1.000  0.594
# 1 CT d&u  437.562  13.253  0.001  0.001
# 1 CT, ev  425.087   0.778  0.678  0.402
# 1 CT, unc 480.403  56.094  0.000  0.000
# 2 CT, d&e 435.650  11.341  0.003  0.002
# 2 CT d&u  447.933  23.624  0.000  0.000
# 2 CT, ev  436.534  12.225  0.002  0.001
# 2 CT, unc 499.834  75.525  0.000  0.000
# 3 CT, d&e 447.236  22.927  0.000  0.000
# 3 CT d&u  460.924  36.615  0.000  0.000
# 3 CT, ev  449.831  25.522  0.000  0.000
# 3 CT, unc 519.703  95.394  0.000  0.000
# 4 CT, d&e 458.276  33.967  0.000  0.000
# 4 CT d&u  473.032  48.723  0.000  0.000
# 4 CT, ev  460.973  36.664  0.000  0.000
# 4 CT, unc 537.601 113.292  0.000  0.000


#So the best model is 1 common trend d&e, (1.1.1) BUT the fit ratio was really bad (above 0.6) so the next best model with a good fit ratio is 1.3.1




###
##
#  WORKING ON THE BEST MODEL
##
###





###############################
## NO VARIMAX ROTATION (YET) ##
###############################

##Plotting the common trends###

bestdfa = dfa.1.3.1

#Common trends
c.trends <- bestdfa$states

##Confidence intervals
ct.UCI <- bestdfa$states + bestdfa$states.se *qnorm(1-0.05/2)
ct.LCI <- bestdfa$states - bestdfa$states.se *qnorm(1-0.05/2)


par.mat = coef(bestdfa, type = "matrix") 
fit.mod = par.mat$Z  %*%  bestdfa$states 

par(mfrow=c(3,1))

#CT1
plot(1:nrow(dat), c.trends[1,], type = 'l' , ylim = c(-6,4), main = 'Common Trend 1', xaxt = 'n', xlab = '', 
     ylab = 'meanCPUE', font.lab=2, yaxt = 'n', lwd = 4)
axis(1, at = c(0,5,10,15,20,25), labels = c('0','5', '10', '15', '20', '25'), font = 2, cex.axis = 0.8)
axis(2, at = c(-8, -6, -4, -2, 0, 2, 4), labels = c('-8', '-6', '-4', '-2', '0', '2', '4'), las = 2, font = 2)
abline(h=0)
lines(1:nrow(dat), ct.UCI[1,], type = 'l', lwd = 4, col = 'gray', lty = 3 )
lines(1:nrow(dat), ct.LCI[1,], type = 'l', lwd = 4, col = 'gray', lty = 3 )

#CT2
plot(1:nrow(dat), c.trends[2,], type = 'l' , ylim = c(-6,4), main = 'Common Trend 2', xaxt = 'n', xlab = '', 
     ylab = 'meanCPUE', font.lab=2, yaxt = 'n', lwd = 4)
axis(1, at = c(0,5,10,15,20,25), labels = c('0','5', '10', '15', '20', '25'), font = 2, cex.axis = 0.8)
axis(2, at = c(-8, -6, -4, -2, 0, 2, 4), labels = c('-8', '-6', '-4', '-2', '0', '2', '4'), las = 2, font = 2)
abline(h=0)
lines(1:nrow(dat), ct.UCI[2,], type = 'l', lwd = 4, col = 'gray', lty = 3 )
lines(1:nrow(dat), ct.LCI[2,], type = 'l', lwd = 4, col = 'gray', lty = 3 )

#CT3
plot(1:nrow(dat), c.trends[3,], type = 'l' , ylim = c(-8,4), main = 'Common Trend 3', xaxt = 'n', xlab = '', 
     ylab = 'Carbon', font.lab=2, yaxt = 'n', lwd = 4)
axis(1, at = c(0,5,10,15,20,25), labels = c('0','5', '10', '15', '20', '25'), font = 2, cex.axis = 0.8)
axis(2, at = c(-8, -6, -4, -2, 0, 2, 4), labels = c('-8', '-6', '-4', '-2', '0', '2', '4'), las = 2, font = 2)
abline(h=0)
lines(1:nrow(dat), ct.UCI[3,], type = 'l', lwd = 4, col = 'gray', lty = 3 )
lines(1:nrow(dat), ct.LCI[3,], type = 'l', lwd = 4, col = 'gray', lty = 3 )

#CT4
plot(1:nrow(dat), c.trends[4,], type = 'l' , ylim = c(-8,4), main = 'Common Trend 3', xaxt = 'n', xlab = '', 
     ylab = 'Carbon', font.lab=2, yaxt = 'n', lwd = 4)
axis(1, at = c(0,5,10,15,20,25), labels = c('0','5', '10', '15', '20', '25'), font = 2, cex.axis = 0.8)
axis(2, at = c(-8, -6, -4, -2, 0, 2, 4), labels = c('-8', '-6', '-4', '-2', '0', '2', '4'), las = 2, font = 2)
abline(h=0)
lines(1:nrow(dat), ct.UCI[4,], type = 'l', lwd = 4, col = 'gray', lty = 3 )
lines(1:nrow(dat), ct.LCI[4,], type = 'l', lwd = 4, col = 'gray', lty = 3 )



## Plotting the fits to the aggregate groups ##

fit = bestdfa
regions = rownames(dat.z)
par(mfrow=c(3,3), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:nrow(dat)){
  plot(dat.z[i,],xlab="",ylab="",xaxt="n",yaxt="n", bty="L",  pch=16, col="black")
  axis(1, at = c(0,5,10,15,20,25,30), labels = c('0','5', '10', '15', '20', '25', '30'), font = 2, cex.axis = 0.8)
  axis(2, at=c(-1, 0, 1, 2, 3), labels = c('-1', '0', '1', '2', '3'), las = 2, font = 2)
  mtext(side = 2, line = 2, text = 'meanCPUE', cex = 0.7, font = 2)
  par.mat=coef(fit,type="matrix")
  lines(as.vector(par.mat$Z[i,,drop=FALSE]%*%fit$states+par.mat$A[i,]) , lwd=2)
  title(regions[i])
}
## If have covariates, add '+ matrix(par.mat$D, nrow = nrow(nm.dat.z))[i,] %*% rbind(zSST,zCOM))' in line 307

1:length(regions)

## Ploting the factor loadings  ##

change <- par.mat$Z

change.1 <- change[seq(1,length(regions)),]

t.z <- t(par.mat$Z)


regions.1 <- 1:length(regions)
minZ = 0.00
m=4
ylims = c(-1.1*max(abs(change.1)), 1.1*max(abs(change.1)))
par(mfrow=c(ceiling(m/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
par(mfrow = c(2,1))

##Common trend 1
plot(c(1:ncol(dat1)), as.vector(change.1),
     type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,ncol(dat1)+1), yaxt = 'n')
axis(side = 2, at=c(-1, 0, 1, 2), labels = c('-1', '0', '1', '2'), font = 2, las = 2)
abline(h=0, lwd=1.5, col="gray")
abline(h=0.2, lty = 2, lwd = 2, col = 'gray')
abline(h=-0.2, lty = 2, lwd = 2, col = 'gray')
for(j in 1:ncol(dat1)) {
  
  if(change.1[j] > minZ) {text(j, -0.05, regions[j], srt=90, adj=1, cex=1, font = 2)}
  if(change.1[j] < -minZ) {text(j, 0.05, regions[j], srt=90, adj=0, cex=1 , font = 2)}
  if(change.1[j] == minZ) {text(j, -0.05, regions[j], srt=90, adj=1, cex=1 , font = 2)}
  
} # end j loop
mtext(paste("Factor loadings on trend 1"),side=3,line=.5, font = 2)

##Common trend 2  
plot(c(1:ncol(dat1)), as.vector(change.1[,2]),
     type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,ncol(dat1)+1), yaxt = 'n')
axis(side = 2, at=c(-1, 0, 1, 2), labels = c('-1', '0', '1', '2'), font = 2, las = 2)
abline(h=0, lwd=1.5, col="gray")
abline(h=0.2, lty = 2, lwd = 2, col = 'gray')
abline(h=-0.2, lty = 2, lwd = 2, col = 'gray')
for(j in 1:ncol(dat1)) {
  
  if(change.1[j,2] > minZ) {text(j, -0.05, regions[j], srt=90, adj=1, cex=1, font = 2)}
  if(change.1[j,2] < -minZ) {text(j, 0.05, regions[j], srt=90, adj=0, cex=1 , font = 2)}
  if(change.1[j,2] == minZ) {text(j, -0.05, regions[j], srt=90, adj=1, cex=1 , font = 2)}
  
} # end j loop
mtext(paste("Factor loadings on trend 2"),side=3,line=.5, font = 2)


##Common trend 3  
plot(c(1:ncol(dat1)), as.vector(change.1[,3]),
     type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,ncol(dat1)+1), yaxt = 'n')
axis(side = 2, at=c(-0.8, -0.4, 0, 0.4, 0.8), labels = c('-0.8', '-0.4', '0', '0.4', '0.8'), font = 2, las = 2)
abline(h=0, lwd=1.5, col="gray")
abline(h=0.2, lty = 2, lwd = 2, col = 'gray')
abline(h=-0.2, lty = 2, lwd = 2, col = 'gray')
for(j in 1:19) {
  
  if(change.1[j,3] > minZ) {text(j, -0.05, fishes[j], srt=90, adj=1, cex=1, font = 2)}
  if(change.1[j,3] < -minZ) {text(j, 0.05, fishes[j], srt=90, adj=0, cex=1 , font = 2)}
  if(change.1[j,3] == minZ) {text(j, -0.05, fishes[j], srt=90, adj=1, cex=1 , font = 2)}
  
} # end j loop
mtext(paste("Factor loadings on trend 3"),side=3,line=.5, font = 2)



############################
### VARIMAX ROTATION    ####
############################

# get the inverse of the rotation matrix
Z.est = coef(bestdfa, type="matrix")$Z
H.inv = 1
if(ncol(Z.est)>1) H.inv = varimax(coef(bestdfa, type="matrix")$Z)$rotmat

# rotate factor loadings (from equation 10.10 in vingette)
Z.rot = Z.est %*% H.inv   
# rotate trends (from equation 10.10 in vingette)
trends.rot = solve(H.inv) %*% bestdfa$states

##Plot rotated common trends

##First, Confidence intervals
ctr.UCI <- trends.rot + bestdfa$states.se *qnorm(1-0.05/2)
ctr.LCI <- trends.rot - bestdfa$states.se *qnorm(1-0.05/2)

##Checks here - gives me confidence in the CIs
trends.rot-c.trends #Same as two below
ctr.UCI-ct.UCI # Same as either side
ctr.LCI-ct.LCI # Same as two above


##Rotated common trends

par(mfrow=c(3,1))
#CT1
plot(1:nrow(dat), trends.rot[1,], type = 'l' , ylim = c(-6,4), main = 'Common Trend 1', xaxt = 'n', xlab = '', 
     ylab = 'Carbon', font.lab=2, yaxt = 'n', lwd = 4, cex.lab = 1.3,  cex.main = 1.4)
axis(1, at = c(0,5,10,15,20,25), labels = c('0','5', '10', '15', '20', '25'), font = 2, cex.axis = 0.8)
axis(2, at = c(-8, -6, -4, -2, 0, 2, 4), labels = c('-8', '-6', '-4', '-2', '0', '2', '4'), las = 2, font = 2)
abline(h=0)
lines(1:nrow(dat), ctr.UCI[1,], type = 'l', lwd = 4, col = 'gray', lty = 3 )
lines(1:nrow(dat), ctr.LCI[1,], type = 'l', lwd = 4, col = 'gray', lty = 3 )

#CT2
plot(1:nrow(dat), trends.rot[2,], type = 'l' , ylim = c(-6,4), main = 'Common Trend 2', xaxt = 'n', xlab = '', 
     ylab = 'Carbon', font.lab=2, yaxt = 'n', lwd = 4, cex.lab = 1.3,  cex.main = 1.4)
axis(1, at = c(0,5,10,15,20,25), labels = c('0','5', '10', '15', '20', '25'), font = 2, cex.axis = 0.8)
axis(2, at = c(-8, -6, -4, -2, 0, 2, 4), labels = c('-8', '-6', '-4', '-2', '0', '2', '4'), las = 2, font = 2)
abline(h=0)
lines(1:nrow(dat), ctr.UCI[2,], type = 'l', lwd = 4, col = 'gray', lty = 3 )
lines(1:nrow(dat), ctr.LCI[2,], type = 'l', lwd = 4, col = 'gray', lty = 3 )

#CT3
plot(1:nrow(dat), trends.rot[3,], type = 'l' , ylim = c(-6,4), main = 'Common Trend 3', xaxt = 'n', xlab = '', 
     ylab = 'Carbon', font.lab=2, yaxt = 'n', lwd = 4, cex.lab = 1.3,  cex.main = 1.4)
axis(1, at = c(0,5,10,15,20,25), labels = c('0','5', '10', '15', '20', '25'), font = 2, cex.axis = 0.8)
axis(2, at = c(-8, -6, -4, -2, 0, 2, 4), labels = c('-8', '-6', '-4', '-2', '0', '2', '4'), las = 2, font = 2)
abline(h=0)
lines(1:nrow(dat), ctr.UCI[3,], type = 'l', lwd = 4, col = 'gray', lty = 3 )
lines(1:nrow(dat), ctr.LCI[3,], type = 'l', lwd = 4, col = 'gray', lty = 3 )


## Plot the rotated factor loadings

##  Z.rot is the new matrix of factor loadings

change <- Z.rot

change.1 <- change[seq(1,length(fishes)),]

t.z <- t(Z.rot)

minZ = 0.00
m=2
ylims = c(-1.1*max(abs(change.1)), 1.1*max(abs(change.1)))
par(mfrow=c(ceiling(m/2),2), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
par(mfrow = c(3,1))

##Common trend 1
plot(c(1:ncol(dat1)), as.vector(change.1[,1]),
     type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,ncol(dat1)+1), yaxt = 'n')
axis(side = 2, at=c(-0.8, -0.4, 0, 0.4, 0.8), labels = c('-0.8', '-0.4', '0', '0.4', '0.8'), font = 2, las = 2)
abline(h=0, lwd=1.5, col="gray")
abline(h=0.2, lty = 2, lwd = 2, col = 'gray')
abline(h=-0.2, lty = 2, lwd = 2, col = 'gray')
for(j in 1:ncol(dat1)) {
  
  if(change.1[j,1] > minZ) {text(j, -0.05, fishes[j], srt=90, adj=1, cex=1, font = 2)}
  if(change.1[j,1] < -minZ) {text(j, 0.05, fishes[j], srt=90, adj=0, cex=1 , font = 2)}
  if(change.1[j,1] == minZ) {text(j, -0.05, fishes[j], srt=90, adj=1, cex=1 , font = 2)}
  
} # end j loop
mtext(paste("Factor loadings on trend 1"),side=3,line=.5, font = 2)

##Common trend 2  
plot(c(1:ncol(dat1)), as.vector(change.1[,2]),
     type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,ncol(dat1)+1), yaxt = 'n')
axis(side = 2, at=c(-0.8, -0.4, 0, 0.4, 0.8), labels = c('-0.8', '-0.4', '0', '0.4', '0.8'), font = 2, las = 2)
abline(h=0, lwd=1.5, col="gray")
abline(h=0.2, lty = 2, lwd = 2, col = 'gray')
abline(h=-0.2, lty = 2, lwd = 2, col = 'gray')
for(j in 1:ncol(dat1)) {
  
  if(change.1[j,2] > minZ) {text(j, -0.05, fishes[j], srt=90, adj=1, cex=1, font = 2)}
  if(change.1[j,2] < -minZ) {text(j, 0.05, fishes[j], srt=90, adj=0, cex=1 , font = 2)}
  if(change.1[j,2] == minZ) {text(j, -0.05, fishes[j], srt=90, adj=1, cex=1 , font = 2)}
  
} # end j loop
mtext(paste("Factor loadings on trend 2"),side=3,line=.5, font = 2)


##Common trend 3  
plot(c(1:ncol(dat1)), as.vector(change.1[,3]),
     type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,ncol(dat1)+1), yaxt = 'n')
axis(side = 2, at=c(-0.8, -0.4, 0, 0.4, 0.8), labels = c('-0.8', '-0.4', '0', '0.4', '0.8'), font = 2, las = 2)
abline(h=0, lwd=1.5, col="gray")
abline(h=0.2, lty = 2, lwd = 2, col = 'gray')
abline(h=-0.2, lty = 2, lwd = 2, col = 'gray')
for(j in 1:ncol(dat1)) {
  
  if(change.1[j,3] > minZ) {text(j, -0.05, fishes[j], srt=90, adj=1, cex=1, font = 2)}
  if(change.1[j,3] < -minZ) {text(j, 0.05, fishes[j], srt=90, adj=0, cex=1 , font = 2)}
  if(change.1[j,3] == minZ) {text(j, -0.05, fishes[j], srt=90, adj=1, cex=1 , font = 2)}
  
} # end j loop
mtext(paste("Factor loadings on trend 3"),side=3,line=.5, font = 2)



##  Getting data needed to plot species aggregates trends, and then plot 

## Get DFA fits


getDFAfits <- function(MLEobj, alpha=0.05, covariates=NULL) {
  fits <- list()
  Ey <- MARSShatyt(MLEobj) # for var() calcs
  ZZ <- coef(MLEobj, type="matrix")$Z # estimated Z
  nn <- nrow(ZZ) # number of obs ts
  mm <- ncol(ZZ) # number of factors/states
  TT <- ncol(Ey$ytT)  # number of time steps
  ## check for covars
  if(!is.null(covariates)) {
    DD <- coef(MLEobj, type="matrix")$D
    cov_eff <- DD %*% covariates
  } else {
    cov_eff <- matrix(0, nn, TT)
  }
  ## model expectation
  fits$ex <- ZZ %*% MLEobj$states + cov_eff
  ## Var in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for(tt in 1:TT) {
    RZVZ <- coef(MLEobj, type="matrix")$R - ZZ%*%VtT[,,tt]%*%t(ZZ)
    SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop=FALSE] %*% t(MLEobj$states[,tt,drop=FALSE])
    VV <- cbind(VV,diag(RZVZ + SS%*%t(ZZ) + ZZ%*%t(SS)))
  }
  SE <- sqrt(abs(VV))
  ## upper & lower (1-alpha)% CI
  fits$up <- qnorm(1-alpha/2)*SE + fits$ex
  fits$lo <- qnorm(alpha/2)*SE + fits$ex
  return(fits)
}


##Gives the model fits for each time series of aggregate abundance with 95% CIs

fit.b = getDFAfits(bestdfa)


fit.b$ex[,1]



fit = bestdfa
par(mfrow=c(4,5), mar=c(3,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:length(fishes)){
  mn <- fit.b$ex[i,]
  plot(dat.z[i,],xlab="",ylab="",bty="L", xaxt="n", yaxt = 'n', ylim=c(-3,3), pch=16, col="blue")
  axis(1, at = c(0,5,10,15,20,25), labels = c('0','5', '10', '15', '20', '25'), font = 2, cex.axis = 0.8)
  axis(2, at = c(-8, -6, -4, -2, 0, 2, 4), labels = c('-8', '-6', '-4', '-2', '0', '2', '4'), las = 2, font = 2)
  mtext(side = 2, line = 2, text = 'Carbon', cex = 0.7, font = 2)
  par.mat=coef(fit,type="matrix")
  lines(as.vector(mn) , lwd=2)
  title(fishes[i])
}





###Assuming SEs of the estimates don't change for states or loadings with rotation

Z.rot
CIS = MARSSparamCIs(bestdfa)


#                   ML.Est Std.Err    low.CI   up.CI
# Z.11             0.35504  0.0643  0.229112  0.4810  sig pos
# Z.21             0.29617  0.0797  0.139980  0.4523  sig pos
# Z.31             0.33672  0.0791  0.181739  0.4917  sig pos
# Z.41             0.34487  0.0682  0.211292  0.4785  sig pos
# Z.51             0.32271  0.0692  0.186982  0.4584  sig pos
# Z.61             0.34273  0.0671  0.211165  0.4743  sig pos
# Z.71             0.35214  0.1327  0.092012  0.6123  sig pos
# Z.81             0.34754  0.1036  0.144556  0.5505  sig pos
# Z.91             0.37192  0.0982  0.179442  0.5644  sig pos
# Z.101            0.26301  0.0967  0.073538  0.4525  sig pos
# Z.111            0.24441  0.0949  0.058372  0.4305  sig pos
# Z.121            0.12237  0.1348 -0.141841  0.3866  not sig
# Z.131           -0.01571  0.1016 -0.214932  0.1835  not sig
# Z.141            0.10518  0.1208 -0.131614  0.3420  not sig
# Z.151            0.37233  0.0878  0.200329  0.5443  sig pos
# Z.161            0.10241  0.1192 -0.131235  0.3361  not sig
# Z.171           -0.00859  0.1384 -0.279777  0.2626  not sig
# Z.181            0.36130  0.1416  0.083675  0.6389  sig pos
# Z.191           -0.11982  0.1243 -0.363437  0.1238  not sig

# Z.22            -0.07978  0.0951 -0.266229  0.1067  not sig
# Z.32             0.15163  0.0822 -0.009528  0.3128  not sig
# Z.42             0.11925  0.0552  0.011067  0.2274  sig pos
# Z.52            -0.13904  0.0514 -0.239793 -0.0383  sig neg
# Z.62             0.10266  0.0592 -0.013344  0.2187  not sig
# Z.72             0.38528  0.1206  0.148829  0.6217  sig pos
# Z.82             0.28417  0.0910  0.105884  0.4625  sig pos
# Z.92             0.26053  0.0787  0.106180  0.4149  sig pos
# Z.102           -0.28567  0.0737 -0.430138 -0.1412  sig neg
# Z.112            0.21129  0.1225 -0.028871  0.4515  not sig
# Z.122           -0.41145  0.1290 -0.664295 -0.1586  sig neg
# Z.132            0.23005  0.1534 -0.070670  0.5308  not sig
# Z.142            0.33260  0.1499  0.038729  0.6265  sig pos
# Z.152            0.20360  0.0747  0.057226  0.3500  sig pos
# Z.162           -0.36558  0.1149 -0.590681 -0.1405  sig neg
# Z.172           -0.43798  0.1372 -0.706828 -0.1691  sig neg
# Z.182            0.36978  0.1338  0.107443  0.6321  sig pos
# Z.192           -0.24915  0.1840 -0.609750  0.1114  not sig



ct1.urot.LCI = as.vector(CIS$par.lowCI$Z[(1:ncol(dat1))])
ct1.urot.UCI = as.vector(CIS$par.upCI$Z[(1:ncol(dat1))])

ct2.urot.LCI = as.vector(CIS$par.lowCI$Z[(20:37)])
ct2.urot.UCI = as.vector(CIS$par.upCI$Z[(20:37)])



ct1.CIdiff <- ct1.urot.UCI - ct1.urot.LCI; ct1.CIdiff
ct2.CIdiff <- ct2.urot.UCI - ct2.urot.LCI; ct2.CIdiff


rot.Z1 <- Z.rot[,1]
rot.Z2 <- Z.rot[,2]


se.ct1 <- as.vector(CIS$par.se$Z[(1:ncol(dat1))])
se.ct2 <- as.vector(c(0,CIS$par.se$Z[(20:37)]))


LCI.rot.Z1 <- rot.Z1 - se.ct1*qnorm(1-0.05/2)
UCI.rot.Z1 <- rot.Z1 + se.ct1*qnorm(1-0.05/2)

ct1.CIdiff.rot <-  UCI.rot.Z1-LCI.rot.Z1; ct1.CIdiff.rot

sum(ct1.CIdiff- ct1.CIdiff.rot)


LCI.rot.Z2 <- rot.Z2 - se.ct2*qnorm(1-0.05/2)
UCI.rot.Z2 <- rot.Z2 + se.ct2*qnorm(1-0.05/2)

ct2.CIdiff.rot <-  UCI.rot.Z2-LCI.rot.Z2; ct2.CIdiff.rot

sum(ct2.CIdiff-ct2.CIdiff.rot[2:ncol(dat1)])


rotloads1 = round(as.data.frame(cbind(LCI.rot.Z1,rot.Z1,UCI.rot.Z1)),digits=3)
rotloads2 = round(as.data.frame(cbind(LCI.rot.Z2,rot.Z2,UCI.rot.Z2)),digits=3)

cbind(fishes,rotloads1)
cbind(fishes,rotloads2)


##ML rotated

# Common Trend 1
#    fishes LCI.rot.Z1 rot.Z1 UCI.rot.Z1
# 1     JB2      0.208  0.334      0.460  sig pos
# 2     JB3      0.095  0.251      0.407  sig pos
# 3     JB4      0.213  0.368      0.523  sig pos
# 4     JB5      0.231  0.365      0.498  sig pos
# 5     JB6      0.120  0.256      0.392  sig pos
# 6     JB7      0.226  0.357      0.489  sig pos
# 7   STUR1      0.202  0.463      0.723  sig pos
# 8   STUR2      0.221  0.424      0.627  sig pos
# 9   STUR3      0.246  0.439      0.631  sig pos
# 10  STUR4     -0.040  0.150      0.339  not sig
# 11   UK25      0.116  0.302      0.488  sig pos
# 12   UK22     -0.290 -0.025      0.239  not sig
# 13  FTL17     -0.135  0.064      0.263  not sig
# 14  FTL26     -0.024  0.212      0.449  not sig
# 15    AL1      0.247  0.419      0.591  sig pos
# 16  AL468     -0.262 -0.029      0.205  not sig
# 17    AL9     -0.429 -0.158      0.114  not sig
# 18  AL181      0.188  0.466      0.743  sig pos
# 19   AL95     -0.441 -0.198      0.046  not sig


# Common Trend 2
#    fishes LCI.rot.Z2 rot.Z2 UCI.rot.Z2
# 1     JB2     -0.121 -0.121     -0.121  sig neg
# 2     JB3     -0.363 -0.176      0.010  not sig
# 3     JB4     -0.134  0.028      0.189  not sig
# 4     JB5     -0.114 -0.006      0.103  not sig
# 5     JB6     -0.342 -0.241     -0.140  sig neg
# 6     JB7     -0.137 -0.021      0.095  not sig
# 7   STUR1      0.005  0.242      0.478  sig pos
# 8   STUR2     -0.030  0.148      0.327  not sig
# 9   STUR3     -0.036  0.118      0.272  not sig
# 10  STUR4     -0.503 -0.358     -0.214  sig neg
# 11   UK25     -0.125  0.115      0.355  not sig
# 12   UK22     -0.681 -0.429     -0.176  sig neg
# 13  FTL17     -0.079  0.222      0.522  not sig
# 14  FTL26     -0.017  0.277      0.571  not sig
# 15    AL1     -0.082  0.064      0.211  not sig
# 16  AL468     -0.604 -0.379     -0.153  sig neg
# 17    AL9     -0.678 -0.409     -0.140  sig neg
# 18  AL181     -0.038  0.224      0.487  not sig
# 19   AL95     -0.554 -0.193      0.167  not sig






