
# test whether estimating three groups separately and together gives different results
# answer: yes
# simulate some data, 4 states
y = round(runif(100)*3)+1
for (i in 1:6)
  y = cbind( y, round(runif(300)*3)+1)
# 7 time points

head(y)
for (i in 1:6)
  print(table(y[,i], y[,(i+1)]) )

cov = round(runif(300)*2)+1; table(cov)
X = cbind(rep(1, times=300), as.numeric(cov == 2),as.numeric(cov == 3))
head(X) # X1 is Intercept, X2 is dummy for group 2, X3 is dummy for group 3

dim(X); dim(y);

# test it for time constant function

# in case of probs test for 2 state version

param= rep(0, times = 12)
out = liktfourstateparest(param, n=nrow(y),K=4,nhist=y,X=NULL) # works
# head(out)
# out = like4states(param=rep((0),12), n=nrow(y),K=4,nhist=y, X=NULL, const=TRUE)
opt1 <- nlm(like4states, p=param[1:12], n=nrow(y),K=4,nhist=y,X=NULL,const=TRUE, hessian=TRUE)

length(whX1)
whX1 = which(cov == 1)
whX2 = which(cov == 2)
whX3 = which(cov == 3)

opt1 = nlm(like4states, p=param[1:12], n=length(whX1),K=4,nhist=y[whX1 ,],X=NULL,const=TRUE, hessian=TRUE)
opt1$estimate
par = param.unpackX(param=opt1$estimate, X=NULL, states=4,n=nrow(y)) ; # if X = NULL, we need to provide n
par[1,]
cis1 =getCI(opt1)
trCI1=transCI(CIS=cis1,X=NULL, states=4, categorical = TRUE) 
opt2 = nlm(like4states, p=param[1:12], n=length(whX2),K=4,nhist=y[whX2 ,],X=NULL,const=TRUE, hessian=TRUE)
opt3 = nlm(like4states, p=param[1:12], n=length(whX3),K=4,nhist=y[whX3 ,],X=NULL,const=TRUE, hessian=TRUE)


out = like4states(param=rep(0,36), n=nrow(y),K=4,nhist=y, X=X, const=TRUE)
optX <- nlm(like4states, p=rep(0,36), n=nrow(y),K=4,nhist=y,X=X,const=TRUE, hessian=TRUE)
optX$estimate
parX = param.unpackX(param=optX$estimate, X=X, states=4) ; head(parX)
getparX(parX=parX, X=X) # get estimates for each group

cis =getCI(optX)
trCI=transCI(CIS=cis,X=X, states=4, categorical = TRUE) # this works!
trCI$lowerCI$`100`
trCI$upperCI$`100`
head(X)
# result: estimating for each group separately gives wider CI's, as expected!


# test looking across time- non -stationary models
# looking at just one group, but time dependent, vs. modelling time as a covariate
head(y)
dim(y)


dim(y)
# first look at no groups - time dependent
out = like4states(param=matrix(0,ncol=6,nrow=12), n=nrow(y),K=4,nhist=y, X=matrix(1,ncol=1,nrow=300), const=FALSE)
out =fitmarkov(init=matrix(0,ncol=6,nrow=12) , nhist=y,n=nrow(y),K=max(y),
                      X=NULL, const=FALSE )
length(out)
names(out)
# opt1 <- nlm(like4states, p=matrix(0,ncol=6,nrow=12), n=nrow(y),K=4,nhist=y,X=matrix(1,ncol=1,nrow=300),const=FALSE, hessian=TRUE)
out$estimate
par = param.unpackX(param=out$estimate, X=NULL, states=4,n=nrow(y)) ; # if X = NULL, we need to provide n
head(par)
cis =getCI(out)
trCI=transCI(CIS=cis,X=NULL, states=4, categorical = TRUE) 
# gives me estimates and CIs for each time point

# version 2 using time dependent covariates
ytime = y[,1:2]
for (t in 2:6) ytime = rbind(ytime, y[,t:(t+1)])
dim(ytime)
covtime = rep(1:6, each = 300); length(covtime)
Xtime = cbind(rep(1, times=1800), as.numeric(covtime == 2),as.numeric(covtime == 3),
              as.numeric(covtime == 4), as.numeric(covtime == 5),as.numeric(covtime == 6) )  ;
head(Xtime) # X1 is Intercept, X2 is dummy for group 2, X3 is dummy for group 3
tail(Xtime) 
outt = like4states(param=rep(0,times=6*12), n=nrow(ytime),K=4,nhist=ytime, X=Xtime, const=TRUE)

outt =fitmarkov(init=rep(0,times=6*12) , nhist=ytime,n=nrow(ytime),K=max(ytime),
               X=Xtime, const=TRUE )
# takes a long time to fit
outt$estimate
part = param.unpackX(param=outt$estimate, X=NULL, states=4,n=nrow(y)) ; # if X = NULL, we need to provide n
cist =getCI(outt)
trCIt=transCI(CIS=cist,X=Xtime, states=4, categorical = TRUE) 

head(part)

# no, different results!

?nlm
