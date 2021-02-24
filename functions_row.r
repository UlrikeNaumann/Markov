
# adaption of the markov model with covariates
# new idea: do the estimation for each time point and for each starting transition separately
# that allows us to feed into the code just the number of transitions that are possible

# function to obtain only the subset of subjects (ts and cov) where at time t, start value is start
# cov needs to be a matrix
# time: starting time point
# start: starting state
# data: time series dataset
gettransrow = function(time, start, data, cov){

  datatst <- data[data[,time]==start ,]
  covst   <-  cov[data[,time]==start ,]
    
  datat = datatst[ , time:(time+1)]
  
  return(list(out = datat, cov = covst) )
}

# for time constant estimation; just append each time point into a 2 column matrix for 
# the time series,
# and for the cov: copy each cov (t-1) times, and append onto each other

ttransconst <- function(data, cov, stop = 7){
  
  stp <- stop - 1;
  covnew = cov;
  datnew <- data[, 1:2];
  for (t in 2:stp){
    covnew <- rbind(covnew, cov);
    datnew <- rbind(datnew, data[, t:(t+1)] );
    
  }
  return( list(out =datnew, cov = covnew ) );
}
# head(y); dim(y)
# ttrandat =ttransconst(data=y, cov=X, stop = 7); head(ttrandat$out); dim(ttrandat$out)
# now we can separately for each startin state, fit the model
# obtain data and cov with starting state 1
# hist1 = gettransrow(time=1, start=1, data=ttrandat$out, cov=ttrandat$cov); dim(hist1$out); dim(hist1$cov)
table(Z[hist1$cov[,2]== 0 & hist1$cov[,3]== 0])
table(Z[hist1$cov[,2]== 1 & hist1$cov[,3]==0])
table(Z[hist1$cov[,2]== 0 & hist1$cov[,3]==1])
# hist2 = gettransrow(time=1, start=2, data=ttrandat$out, cov=ttrandat$cov); dim(hist2$out)
# hist3 = gettransrow(time=1, start=3, data=ttrandat$out, cov=ttrandat$cov); dim(hist3$out)
# hist4 = gettransrow(time=1, start=4, data=ttrandat$out, cov=ttrandat$cov); dim(hist4$out)
param = c( c(12, 19, 20)/sum(c(12, 19, 20, 12)),  c(13, 40, 56)/sum(c(13, 40, 56, 27 )) ,c(8, 30, 26)/sum(c(8, 30, 26,  8  )) )
out =likXstateparest(param, n=length(which(y[,1]==1)),1:4 ,K=4, start=1, nhist=y[which(y[,1]==1),1:2],X=NULL)
out =likXstateparest(rep(0, times = 9), n=nrow(hist1$cov),1:4 ,K=4, start=1, nhist=hist1$out,X=hist1$cov)
out =likXstateparest(param, n=nrow(hist1$cov),1:4 ,K=4, start=1, nhist=hist1$out,X=hist1$cov)

opt1 <- nlm(likXstateparest, p=param, n=nrow(hist1$cov),whichst=1:4,K=4,nhist=hist1$out,X=hist1$cov, hessian=TRUE)
opt1
# now we need to adapt the fct for getting CI and transfCI
# could test by caclulating by hand the prob for each value of X, and using correct value as init
 #par = param.unpackKX(param=rep(0, times = 9), X=hist1$cov, states=4); head(par)
# table(Z) # each one is freq often
# Fehler in f(x, ...) : some problem with the probability limits - probably bad starting values
# ie some estimated transitions < 0 or >1
# normally this should give us the time constant transition probs

# logit function
logit=function(p){
  
  if(length(p)== 1){
    logitp = log(p/(1-p))
    return(logitp)
  } else {
    
    # if length(p) > 1: we assume that the last p is missing, eg the first 3 of 4 are provided
    pK <- 1 - sum(p);
    sum_exp_bj <- 1/pK;
    beta <- log(p * sum_exp_bj);
    return(beta)
  }
}

# for multinomial prob with K classes, supply beta for the first 3 (not the one which is set to 0)
expit = function(beta){
  
  if(length(beta) == 1)
    return( exp(beta)/(1+exp(beta)) )
  
  # we assume that the value betaK = 0 is not included!
  if(length(beta)>1){
    return( exp(beta)/(1+sum(exp(beta))) )
  }
}




# this version will not work
like4states <- function(param, n=nrows(nhist),K=4,nhist, X=NULL, const=FALSE){
  
  dims <- dim(nhist);
  T <- dims[2];
  Tm1 = T-1;
  likelihood <- c();
  
 if(!const & !is.null(X) ) { #time dependent, implementation with covariates
    
    for(i in 1:Tm1){
      # nhistT <- nhist[ ,i:(i+1)];
      nhistT =  nhist[,i:(i+1)];
      if(is.null(dim(param))) paramT <- param else {
        paramT <- param[,i]   }
      
      likT <-  liktfourstateparest(param=paramT, n=dims[1],K=K,X=X,nhist=nhistT);
      likelihood <- c(likelihood, likT);
    }
    #  likelihood <- sum(likelihood) ;  # check whether this works for time dept.
    
  } 
  
  return(likelihood);
  #  return(sum(likelihood));
}




# get log likelihood for >2 state model, incl design matrix with covariates X
# here we only look at one starting state
# whichst: specifies which end states we look at, ie don't include any which are fixed at 0
# for categorical covariates!
# param for >2 states has length = ncol(X)*2
# param = (b00,b01, b0p, b10, b11, b1p)
# nhist is in this case (general form) a 2 dim matrix with first and 2nd time point for each subject

# out =likXstateparest(param, n=length(which(y[,1]==1)),1:4 ,K=4, start=1, nhist=y[which(y[,1]==1),1:2],X=NULL)
### start = 3; whichst = c(1,3,4);
likXstateparest <- function(param, n=nrow(nhist),whichst ,K=length(whichst), start=NULL, nhist,X=NULL){
  
  if(is.null(X)) X <- matrix(1, ncol =1, nrow =n) ;
  if(is.null(start)) start = unique(nhist[,1]);
  if (length(start ) != 1) stop("error: length of start not correct");
  
  # create a short form for the 2 state model history a la 11,12,21,22
  Z = apply(nhist, MARGIN=1, function(x) paste0(x[1], x[2]) )
  # assuming states are a and 2
  # n1 =  sum(Z ==  "11") +  sum(Z ==  "12") ;
  #   n2 =  sum(Z ==  "21") +  sum(Z ==  "2                                                                                                                                                       2") ;
  # assuming we have an intercept and dummy variables in X
  p <- ncol(X); # therefore we assume that we have p groups
  
  # here we get back a matrix with entries for each subject head(par); dim(par)
  par = param.unpackKX(param=param, X=X, states=K);
  
  # to make life easier in the next step, create the full prob matrix (as vector)
  # here we assume that param.unpackX returns a matrix with (K-1) columns!
  parK <- matrix(0, ncol = K,nrow = n);
  for(i in 1:1){ # could be simplified, here: only calculate for one starting state
    parK[ ,   (1:(K-1)) ] <-   par[,  (1:(K-1)) ]
    parK[ ,   K ] <- 1 - rowSums( par[,   (1:(K-1)) ] );
  }
  if(any(parK< 0 | parK >1) ) stop("some problem with the probability limits")
  # if we have only categorical covariates,we could do this more efficient, but this version is more general
  
  Lj = 0;
  # unhist <- unique(Z);
  # unhist <-paste0(rep(1:K, each=K-1), 1:(K-1)) 
  # unhist <-paste0(rep(1:K, each=K), 1:(K))
  unhist <-paste0(rep(start, K), whichst ) # this has the length of whichst
  
  
  nj <- length(unhist); 
  # L0 = 0; L1 = 0; # initialise
  # problem: if nj 
  
  for(j in 1:nj){
    for(i in 1:nrow(X)){
      Lj <- Lj + as.numeric(Z[i]== unhist[j]) *log(parK[i,j]) # +  as.numeric(Z[i]== unhist[j])*log(1-par[i,j]) ;# slightly different from formula in book
    }
  }
  
  if(Lj < -100000) Lj = -100000
  return(-Lj);
}



# works with X with covariates states=K, X = NULL, not allowed 
# (for creating X, we would need to know n)
# look only at one starting state
param.unpackKX <- function(param,X, states,n ){
  # length of param is assumed to be p*(K-1) with K = 4
  # if(is.null(X))
  lk = length(param)
  if(is.null(X)) X <- matrix(1, ncol =1, nrow =n) ; # needed below
  p <- ncol(X);
  # if (lk = 2) states = 2; if(lk = 12) states = 4;
  l = length(param); # only looking at one starting state
  
  
  # for the first state dim(Xbeta1) dim(output)
  # if X has just one column, then it is simply each value of p times X
  if(is.null(X) | p == 1){
    k = 1;
   #  paramx <- param[(k-1)*(l/p) +  1:(l/p)  ];
    paramx <- param[ 1:lk  ];
    
    # each of the below are vectors for each subject
    p1_ = exp(paramx)/(1+ sum(exp(paramx)) );
    output <- p1_ ;
   # for(j in 2:states){ 
   #    paramx <- param[( (j-1)*l+1):(j*l)];
   #    pj_ = exp(paramx)/(1+ sum(exp(paramx)) );
   #    output <- c( output,pj_ )
   #  }
    outputx = matrix(rep(output, times=nrow(X)), nrow =nrow(X),ncol= lk, byrow =TRUE)
    return(outputx)
  }
  
  
  # I need to work through this again, with an example to get the dimensions right
  # for the first state
  param1 <- param[1:(p*(states-1))  ]; # param for the first row of the trans.matrix= (p * (states-1))
  # now loop through each element - for each element, we have p params
  # Xbeta1 <- matrix(0, nrow=nrow(X), ncol=(states-1)) ; # each element s1i of the transition matrix row gets a column
  
  # instead of this loop, we could do a matrix multiplication.
  # arrange each of the (states-1) times p elements of paramn1 into a p x (states-1) matrix
  parS <- matrix(param1, ncol = (states-1), nrow =p  ,byrow=T);
  Xbeta1 = X%*%parS;
  

  # for p = 1 the dimension here is incorrect
  # each of the below are vectors for each subject
  p1_ = exp(Xbeta1)/(1+ rowSums(exp(Xbeta1)) ); # I think these should be rowSums, as we have one row for each subject
  output <- p1_ ; # this is the first state result
  
  # for(j in 2:states){ # loop through each remaining state
  #  paramj <- param[(j-1)*(p*(states-1)) +  1:(p*(states-1))  ]; 
  #  parS <- matrix(paramj, ncol = (states-1), nrow =p  );
  #  Xbetaj = X%*%parS;
    # each of the below are vectors for each subject
  #  pj_ = exp(Xbetaj)/(1+ rowSums(exp(Xbetaj)) );
  #  output <- cbind( output,pj_ ) 
  # }
  
  return( output);  #  as output we get the matrix with K*K-1 columns
}



# for X not NULL, and time dependent: init =matrix(0,ncol=6,nrow=2)
# optimisation is done for each time change separately for time dependent
# output: diaghessian is a matrix of diagonals for each time point
fitmarkov <- function(init , nhist,n=nrow(nhist),K=max(nhist),
                      X=NULL, const=FALSE ){
  # opt3 <- nlm(like2states, p=par, n=100,K=2,nhist=y,X=NULL,const=FALSE, hessian=TRUE)
  #instead: loop through each time point, possible here: as we have markov model of order 1
  # h <- matrix(0,ncol=2*(ncol(nhist)-1),nrow=2*(ncol(nhist)-1));
  Tm1 <- ncol(nhist)-1
  
  # some dimension checks
  if(!const){
    if( is.null(nrow( init)) ) 
      init = matrix(rep(init, times = Tm1), nrow =length(init),ncol=Tm1)
    
    lpp <- nrow(init) ;
    
  } else {
    lpp <- length(init) ;
  }
  if(is.null(X)){
    if(lpp != K*(K-1)) stop("dimension of init is not correct")
  } else {
    if(lpp != ncol(X)*K*(K-1)) stop("dimension of init is not correct")
  }
  
  if(const){
    if(K == 2) {  opt2 <- nlm(like2states, p=init, n=n,K=K,nhist=nhist,X=X,const=TRUE, hessian=TRUE)
    } else if (K > 2) {
      
      opt2 <- nlm(like4states, p=init, n=n,K=K,nhist=nhist,X=X,const=TRUE, hessian=TRUE)
    }
    
    return(opt2)
  } # else
  
  h <- array(0,dim = c(lpp,lpp , Tm1  ))
  
  # class(h) <- c("matrix", "array" )
  opt3 <- list(minimum= c(),estimate= matrix(0,nrow=lpp, ncol=Tm1),
               gradient = matrix(0,nrow=lpp, ncol=Tm1), pseudohessian = h,
               diaghessian = matrix(0,nrow=lpp, ncol=Tm1),code =c(), iterations = c());
  
  #names(opt3) <- c("minimum", "estimate", "gradient", "hessian", "code", "iterations");
  
  for(t in 1:Tm1){
    
    if(K == 2){
      opt3a <- nlm(like2states, p=init[,t], n=n,K=K,nhist=nhist[,t:(t+1)],X=X,const=TRUE, hessian=TRUE)
    } else {
      opt3a <- nlm(like4states, p=init[,t], n=n,K=K,nhist=nhist[,t:(t+1)],X=X,const=TRUE, hessian=TRUE)
    }
    # ie fit the small model - could also just use the sub-model
    opt3$minimum <- c(opt3$minimum, opt3a$minimum);  
    opt3$estimate[,t] <-  opt3a$estimate;
    opt3$gradient[,t]  <- opt3a$gradient;  
    opt3$code <- c(opt3$code, opt3a$code);
    opt3$iterations <- c(opt3$iterations, opt3a$iterations); 
    # opt3$pseudohessian[(t*2-2)+ ( 1:2) ,(t*2-2)+( 1:2)] <- opt3a$hessian
    opt3$pseudohessian[, , t] <- opt3a$hessian; 
    opt3$diaghessian[, t] <- diag(opt3a$hessian); 
  }
  return(opt3)
}
