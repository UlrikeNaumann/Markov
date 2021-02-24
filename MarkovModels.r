
# functions for estimating transition probs by brute force optimisation of likelihood
# likelihood, p-values are based on
# M. Ataharul Islam, Rafiqul Islan Chowdhury, Shahariar Huda, Markov Models With Covariate Dependence For Repeated Measures, 2009, Nova Science Publishers, Inc. 
# M. Ataharul Islam, Rafiqul Islan Chowdhury, Analysis of Repeated Measures Data, 2017, Springer, pp. 51 - 66

# running time increases with the number of time points, but roughly linear.
# running time icnreases with number of covariates in X! - more then linear!


# before each run, make sure to load the correct version of the function.
# versions: general: all states estimated
# IQVIA: 4: only estimates transitions from first 3 states, and not trans_3_4
# functions_row: estimate transitions from each state separetly

# library(seqHMM)

# simulate some data, 2 states
y = round(runif(100))+1
for (i in 1:6)
  y = cbind( y, round(runif(300))+1)
head(y)
# 7 time points

# simulate some data, 4 states
y = round(runif(100)*3)+1
for (i in 1:6)
y = cbind( y, round(runif(300)*3)+1)
# 7 time points

 head(y)
 for (i in 1:6)
print(table(y[,i], y[,(i+1)]) )
 
 
x = round(runif(100))+1
head(x);
y = round(runif(100))+1
nhist =c(t(table(x,y))); (nhist)

cov = round(runif(300)*2)+1; table(cov)
X = cbind(rep(1,times = 300), as.numeric(cov == 2),as.numeric(cov == 3))
head(X) # X1 is Intercept, X2 is dummy for group 2, X3 is dummy for group 3

dim(X); dim(y);

Y = cbind(x,y); head(Y)

Z = apply(Y, MARGIN=1, function(x) paste0(x[1], x[2]) )
head(Z)

# test like2states
y = round(runif(100))+1
for (i in 1:6)
  y = cbind( y, round(runif(100))+1)
head(y)


a = matrix(1:12,4,3)
a %*% 1:3
# to do: depending on value of const, and X --> set some warning message if param does not have correct length

# test wrapper function - 1 run
like2states(param=rep(0,2), n=nrows(y),K=2,nhist=y, X=NULL, const=FALSE) # because our starting values
# are the same, and we just add up everything!
like2states(param=rep(0,2), n=nrows(y),K=2,nhist=y, X=NULL, const=TRUE)

# with covariates!, starting parameter needs to have correct dimension!
like2states(param=c(0, -5,-5,0,-5,-5), n=nrows(y),K=2,nhist=y, X=X, const=FALSE) # because our starting values
# are the same, and we just add up everything!
like2states(param=c(0, -5,-5,0,-5,-5), n=nrows(y),K=2,nhist=y, X=X, const=TRUE)

# fitmarkov function is not really needed if const = TRUE
# choo
opt1 <- fitmarkov(init=rep(0,2) , nhist=y,n=nrow(nhist),K=max(y),X=NULL, const=FALSE )

opt2 <- fitmarkov(init=c(0, -5,-5,0,-5,-5) , nhist=y,n=nrow(nhist),K=max(y),X=X, const=FALSE )
cisopt2 = getCI(opt1)
transCI(CIS=getCI(opt1),X=NULL, states=2, categorical = TRUE)
param.unpack1(param=opt1$estimate)[1,] # all works

dim(X)
# test optimisation of the wrapper function
# to do. put nlm optimisation into the getCI
#const - no covariate, don't forget to set const=TRUE, otherwise p needs to have different dimension!
opt1 <- nlm(like2states, p=rep(0,2), n=100,K=2,nhist=y,X=NULL,const=TRUE, hessian=TRUE)
getCI(opt1)
transCI(CIS=getCI(opt1),X=NULL, states=2, categorical = TRUE)
param.unpack1(param=opt1$estimate)

#opt2 <- nlm(liktwostateparest, p=c(0, -5,-5,0,-5,-5), n=100,K=2,nhist=y,X=X, hessian=TRUE)
opt2 <- nlm(like2states, p=c(0, -5,-5,0,-5,-5), n=100,K=2,nhist=y,X=X,const=TRUE, hessian=TRUE)
# Warnmeldung:In nlm(like2states, p = c(0, -5, -5, 0, -5, -5), n = 100, K = 2,  : NA/Inf durch größte positive Zahl ersetzt
getCI(opt2)
opt2a <- nlm(like2states, p=c(0, -5,-5,0,-5,-5), n=100,K=2,nhist=y[,1:2],X=X,const=TRUE, hessian=TRUE)

head(nhist); max(nhist)
param[1:12]
out = liktfourstateparest(param, n=nrow(nhist),K=4,nhist,X=NULL) # works

nhistT =c(t(table(nhist[,1],nhist[,2])))
head(ynew);tail(ynew)
out = liktwostateest(param=param[1:12], n=sum(nhistT),K=4,nhist=nhistT)

head(y)

# like2states works also for 4 states if we have  X=NULL, const=TRUE
out = like4states(param=rep((0),12), n=nrow(y),K=4,nhist=y, X=NULL, const=TRUE)
opt1 <- nlm(like4states, p=param[1:12], n=nrow(y),K=4,nhist=y,X=NULL,const=TRUE, hessian=TRUE)

opt1 <- nlm(like4states, p=rep((0),12), n=nrow(y),K=4,nhist=y,X=NULL,const=TRUE, hessian=TRUE)
opt1$estimate
par = param.unpackX(param=opt1$estimate, X=NULL, states=4,n=nrow(y)) ; # if X = NULL, we need to provide n
cis =getCI(opt1)
trCI=transCI(CIS=cis,X=NULL, states=4, categorical = TRUE) 
# some bug occurs
# Fehler in Xbeta1[, k] <- X %*% (param[(k - 1) * (l/p) + 1:(l/p)]) : 
# Anzahl der zu ersetzenden Elemente ist kein Vielfaches der Ersetzungslänge

out = like4states(param=rep(0,36), n=nrow(y),K=4,nhist=y, X=X, const=TRUE)
opt1 <- nlm(like4states, p=rep(0,36), n=nrow(y),K=4,nhist=y,X=X,const=TRUE, hessian=TRUE)
opt1$estimate
par = param.unpackX(param=opt1$estimate, X=X, states=4) ; head(par)
cis =getCI(opt1)
trCI=transCI(CIS=cis,X=X, states=4, categorical = TRUE) # this works!


#the time dept version for nlm needs param as matrix
param = matrix(0, ncol = ncol(y)-1, nrow = 12)
# like2states works also for 4 states if we have  X=NULL, const=TRUE
out = like4states(param=rep((0),12), n=nrow(y),K=4,nhist=y, X=NULL, const=FALSE)
out = like4states(param=param, n=nrow(y),K=4,nhist=y, X=NULL, const=FALSE) # both work

# opt1 <- nlm(like4states, p=param, n=nrow(y),K=4,nhist=y,X=NULL,const=FALSE, hessian=TRUE)
# for the time dependent version we need to use fitmarkov
# opt1 <- fitmarkov(init=rep(0,2) , nhist=y,n=nrow(nhist),K=max(y),X=NULL, const=FALSE )

# opt2 <- fitmarkov(init=c(0, -5,-5,0,-5,-5) , nhist=y,n=nrow(nhist),K=max(y),X=X, const=FALSE )


opt1x <- like4states(param =rep(0,12) , nhist=y[,1:2],n=nrow(y),K=max(y),X=NULL, const=TRUE )
opt1x <- nlm(like4states, p=rep(0,12), n=nrow(y),K=4,nhist=y[,1:2],X=NULL,const=TRUE, hessian=TRUE)
cis =getCI(opt1x)
trCI=transCI(CIS=cis,X=NULL, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf


head(y)
opt1 <- fitmarkov(init=rep(0,12) , nhist=y,n=nrow(y),K=max(y),X=NULL, const=FALSE )
opt1$estimate
# we need a param.unpack function for time dependent
par = param.unpackX(param=opt1$estimate, X=NULL, states=4,n=nrow(y))[1,] ; # if X = NULL, we need to provide n
# here, each time point is just stacked next to each other
# better advice: do not directly use param.unpackX for time dependent
cis =getCI(opt1)
# check cis
# we need a version for time dependent!
trCI=transCI(CIS=cis,X=NULL, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf
is.vector(1:5)
head(y); head(X)
opt2 = NULL;
opt2 <- fitmarkov(init=rep(0,36) , nhist=y,n=nrow(y),K=max(y),X=X, const=FALSE )
opt2$estimate
#par = param.unpackX(param=opt2$estimate, X=X, states=4,n=nrow(y))[1,] ; # if X = NULL, we need to provide n

par = transEst(opt2$estimate,X=X, states=4, categorical = TRUE)
# better advice: do not directly use param.unpackX for time dependent
# better to create a separate function!!
# that is not correct!
cis =getCI(opt2)
trCI=transCI(CIS=cis,X=X, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf
names(trCI)

# add a convenience function that transforms these either into transition matrix form or into
# big matrix # with variable trt, transition estimate, lower CI, upper CI


expit(logit(c(0,0,0)))

# does it work without adjustment for our type of data, where 4 -- 4
head(y)
# adjust y: once state = 4, it will stay at 4
yadj = y;
for(t in 1:6){  # t = 1
  wh4 <- which(y[,t] == 4)
  for(tj in (t+1):7)
    yadj[wh4,tj] <- 4
}
yadj[1:30,]
# subjects can't dropout at the baseline
# delete lines with not allowed things, eg, start with 4 not allowed
yadj <- yadj[- which(yadj[,1] == 4),]
dim(yadj)

for (i in 1:6)
  print(table(yadj[,i], yadj[,(i+1)]) )

opt1x <- like4states(param =rep(0,12) , nhist=yadj[,1:2],n=nrow(yadj),K=max(yadj),X=NULL, const=TRUE )
opt1x <- nlm(like4states, p=rep(0,12), n=nrow(yadj),K=4,nhist=yadj[,1:2],X=NULL,const=TRUE, hessian=TRUE)

opt1x$estimate # estimation runs
# but for not-occuring transitions, estimate stays at the initial value
cis =getCI(opt1x)
transEst(opt1x$estimate, X=NULL, states=4, categorical = TRUE)
# make an adjustion: if opt1x$estimate == 0 (and opt1x$gradient == 0 / diag(opt1x$hessian) == 0 ) then
# then set cis$lowerCI and cis$upperCI to NA
par = param.unpackX(param=opt1x$estimate, X=NULL, states=4, n=nrow(yadj)) ; head(par)

trCI=transCI(CIS=cis,X=NULL, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf

# can I fix it by specifying the starting values?
par = c(rep(0,9), logit(c(0,0,0)) )
opt1x <- like4states(param =par , nhist=yadj[,1:2],n=nrow(yadj),K=max(yadj),X=NULL, const=TRUE, fixed=4 )

opt1x <- nlm(like4states, p=par, n=nrow(yadj),K=4,nhist=yadj[,1:2],X=NULL,const=TRUE, fixed=4, hessian=TRUE)
# since we don't estimate some parts of par, maybe we need to leave it out?
# no, does not work like this
K = 4
par4 = param.unpackX(param=par, X=NULL, states=4, n=nrow(yadj)); # that does work
par = par4
head(parK)

head(parK)# parK part also works!
dim(parK)
param = rep(0,9)
par[10:12] = -Inf
param[]


######################### test whether we get the same result for testing 2 and 7 time points
param
# change a little for test
param = c(0.25, 0.3, 0.3, 0.15, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, 0.25, 0.25)
# param = matrix(0.25, ncol= 6, nrow = 3 )
# param = rep(0, times = 12)

head(y)
opt1 <- fitmarkov(init=param , nhist=y,n=nrow(y),K=max(y),X=NULL, const=FALSE )
opt1$estimate[,1] # Achtung: for each = 0, each time point gets one column
dim(opt1$estimate)
# we need a param.unpack function for time dependent
par = param.unpackX(param=opt1$estimate, X=NULL, states=4,n=nrow(y))[1,] ; # if X = NULL, we need to provide n
# here, each time point is just stacked next to each other
# better advice: do not directly use param.unpackX for time dependent
cis =getCI(opt1)
# check cis
# we need a version for time dependent!
trCI=transCI(CIS=cis,X=NULL, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf


opt2 <- fitmarkov(init=param , nhist=y[,1:3],n=nrow(y),K=max(y),X=NULL, const=FALSE )
opt2$estimate # estimates are the same
# we need a param.unpack function for time dependent
par2 = param.unpackX(param=opt2$estimate, X=NULL, states=4,n=nrow(y))[1,] ; # if X = NULL, we need to provide n
# here, each time point is just stacked next to each other
# better advice: do not directly use param.unpackX for time dependent
cis2 =getCI(opt2)
# check cis
# we need a version for time dependent!
trCI2=transCI(CIS=cis2,X=NULL, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf

trCI[[1]]$lowerCI
trCI[[2]]$upperCI # confirmed, these are the same




param = rep( opt1$estimate[,1], times = 3)
is.vector(1:5)
head(y); head(X); nrow(X)
opt2 = NULL;
opt2 <- fitmarkov(init=param , nhist=y,n=nrow(y),K=max(y),X=X, const=FALSE )
opt2$estimate
#par = param.unpackX(param=opt2$estimate, X=X, states=4,n=nrow(y))[1,] ; # if X = NULL, we need to provide n

par2 = transEst(opt2$estimate,X=X, states=4, categorical = TRUE)
# better advice: do not directly use param.unpackX for time dependent
# better to create a separate function!!
# that is not correct!
cis2 =getCI(opt2)
trCI2=transCI(CIS=cis2,X=X, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf
names(trCI)
trCI2[[1]]$lowerCI
trCI2[[2]]$upperCI

# now try only the first 3 time points
opt4 <- fitmarkov(init=param , nhist=y[,1:3],n=nrow(y),K=max(y),X=X, const=FALSE )
opt4$estimate # check whether first 2 time points the same result; yes
#par = param.unpackX(param=opt2$estimate, X=X, states=4,n=nrow(y))[1,] ; # if X = NULL, we need to provide n

par4 = transEst(opt4$estimate,X=X, states=4, categorical = TRUE)
# better advice: do not directly use param.unpackX for time dependent
# better to create a separate function!!
# that is not correct!
cis4 =getCI(opt4)
trCI4=transCI(CIS=cis4,X=X, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf
names(trCI)

trCI2[[1]]$lowerCI ; trCI4[[1]]$lowerCI
trCI2[[2]]$upperCI; trCI4[[2]]$upperCI
# check whether these are the same results: yes

############# when using constant time, does it matter if I provide an nhist matrix
# or whether I just append all time points. First try without X then with X (appended=)
param = rep(0, times = 12)
opt2 <- fitmarkov(init=param , nhist=y[,1:3],n=nrow(y),K=max(y),X=NULL, const=TRUE )
opt2$estimate # estimates are the same
# we need a param.unpack function for time dependent
par2 = param.unpackX(param=opt2$estimate, X=NULL, states=4,n=nrow(y))[1,] ; # if X = NULL, we need to provide n
# here, each time point is just stacked next to each other
# better advice: do not directly use param.unpackX for time dependent
cis2 =getCI(opt2)
# check cis
# we need a version for time dependent!
trCI2=transCI(CIS=cis2,X=NULL, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf

ynew = rbind( y[,1:2], y[,2:3])
dim(ynew)
opt2new <- fitmarkov(init=param , nhist=ynew,n=nrow(ynew),K=max(ynew),X=NULL, const=TRUE )
opt2new$estimate # estimates are the same
# we need a param.unpack function for time dependent
par2new = param.unpackX(param=opt2new$estimate, X=NULL, states=4,n=nrow(y))[1,] ; # if X = NULL, we need to provide n
# here, each time point is just stacked next to each other
# better advice: do not directly use param.unpackX for time dependent
cis2new =getCI(opt2new)
# check cis
# we need a version for time dependent!
trCI2new=transCI(CIS=cis2new,X=NULL, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf
# estimates are the same!!!

########### with X #############
param = rep(0, times = 36)
opt2 <- fitmarkov(init=param , nhist=y[,1:3],n=nrow(y),K=max(y),X=X, const=TRUE )
opt2$estimate # estimates are the same
# we need a param.unpack function for time dependent
par2 = param.unpackX(param=opt2$estimate, X=X, states=4,n=nrow(y))[1,] ; # if X = NULL, we need to provide n
# here, each time point is just stacked next to each other
# better advice: do not directly use param.unpackX for time dependent
cis2 =getCI(opt2)
# check cis
# we need a version for time dependent!
trCI2=transCI(CIS=cis2,X=X, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf

ynew = rbind( y[,1:2], y[,2:3])
Xnew = rbind(X, X)
dim(ynew)
opt2new <- fitmarkov(init=param , nhist=ynew,n=nrow(ynew),K=max(ynew),X=Xnew, const=TRUE )
opt2new$estimate # estimates are the same
# we need a param.unpack function for time dependent
par2new = param.unpackX(param=opt2new$estimate, X=Xnew, states=4,n=nrow(y))[1,] ; # if X = NULL, we need to provide n
# here, each time point is just stacked next to each other
# better advice: do not directly use param.unpackX for time dependent
cis2new =getCI(opt2new)
# check cis
# we need a version for time dependent!
trCI2new=transCI(CIS=cis2new,X=Xnew, states=4, categorical = TRUE) # gives a mistake due to Inf, -Inf
# estimates are the same!!!

opt2$estimate
opt2new$estimate # there are small differences at the 6th decimal only
trCI2   # there are small differences at the 6th decimal only
trCI2new
# FAzit: yes, it works