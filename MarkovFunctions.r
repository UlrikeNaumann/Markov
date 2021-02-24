
# next steps: test: 4 state Markov with covariates!
# next: include in programs: some transitions are fixed to a pre-set value. e.g 0 or 1
# for 4 states: if I use such an example, will the resulting estimates automatically be correct?
# or adjust afterwards?


# functions for estimating transition probs by brute force optimisation of likelihood
# likelihood, p-values are based on
# M. Ataharul Islam, Rafiqul Islan Chowdhury, Shahariar Huda, Markov Models With Covariate Dependence For Repeated Measures, 2009, Nova Science Publishers, Inc. 
# M. Ataharul Islam, Rafiqul Islan Chowdhury, Analysis of Repeated Measures Data, 2017, Springer, pp. 51 - 66


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


# works t = 1

# pp <- c(0, -5,-5,0,-5,-5) # assuming X with 3 columns, and 2 states
# par = matrix(rep(pp, times = 6), nrow =6,ncol=6) # assuming 7 time points

# for X not NULL, and time dependent
# fitmarkov <- function(init , n=100,K=2,nhist, X, const=FALSE, ){
# lpp <- nrow(init) ;
#  h <- array(0,dim = c(lpp,lpp , ncol(nhist)-1  ))
# opt4 <- list(minimum= c(),estimate= matrix(0,nrow=lpp, ncol=ncol(nhist)-1),
#              gradient = matrix(0,nrow=lpp, ncol=ncol(nhist)-1), pseudohessian = h,
#              diaghessian = matrix(0,nrow=lpp, ncol=ncol(nhist)-1), code =c(), iterations = c());
# for(t in 1:(ncol(nhist)-1)){
#   opt4a <- nlm(like2states, p=init[,t] , n=100,K=2,nhist=nhist[,t:(t+1)],X=X,const=TRUE, hessian=TRUE)
#   opt4$minimum <- c(opt4$minimum, opt4a$minimum);  
#   opt4$estimate[,t] <-  opt4a$estimate;
#   opt4$gradient[,t]  <- opt4a$gradient;  
#   opt4$code <- c(opt4$code, opt4a$code);
#   opt4$iterations <- c(opt4$iterations, opt4a$iterations); 
#   opt4$pseudohessian[, , t] <- opt4a$hessian; 
#   opt4$diaghessian[, t] <- diag(opt4a$hessian); 
# }
# return(opt4)
# }


# opt4
# use diaghessian: include only the diagonals, per time point
# to do: apply the same for X = NULL, ie diaghessian and pseudohessian


# estimate one state matrix

#tmat = matrix(c(a11,a12,a13,1 - a11-a12-a13,   a21,a22,a23,1 - a21-a22-a23,  
 # dim(1:5)              a31,a32,a33,1 - a31-a32-a33,   0,0,0,1   )  )

# like2states: wrapper function for liktwostateest/liktwostateparest, for several time points, no covariates
# this function calculates the likelihood
# nhist: 2 states are assumed to be named 1 and 2!
# nhist being the entire time series, subjects in rows, time columns
# param here being a matrix, with time in columns, the parameters across rows, or a 
#       vector, where parameters identical across all time points
# if const = TRUE, the time constant transition probs will be calculated, only 1 likel. value returned
like2states <- function(param, n=nrows(nhist),K=2,nhist, X=NULL, const=FALSE){
  
  dims <- dim(nhist);
  T <- dims[2];
  Tm1 = T-1;
  likelihood <- c();
  
  # fuer const == FALSE
  if(is.null(X) & !const){ # implementation without covariates
  
  for(i in 1:Tm1){
    # nhistT <- nhist[ ,i:(i+1)];
    nhistT =c(t(table(nhist[,i],nhist[,i+1])))
    if(is.null(dim(param))) paramT <- param else {
    paramT <- param[,i]   }
    
  likT <-  liktwostateest(param=paramT, n=dims[1],K=K,nhist=nhistT);
  likelihood <- c(likelihood, likT);
  }
 # likelihood <- sum(likelihood) ;  # check whether this works for time dept.
  } else if(!const & !is.null(X) ) { #time dependent, implementation with covariates
    
    for(i in 1:Tm1){
      # nhistT <- nhist[ ,i:(i+1)];
      nhistT =  nhist[,i:(i+1)];
      if(is.null(dim(param))) paramT <- param else {
        paramT <- param[,i]   }
      
      likT <-  liktwostateparest(param=paramT, n=dims[1],K=K,X=X,nhist=nhistT);
      likelihood <- c(likelihood, likT);
    }
   #  likelihood <- sum(likelihood) ;  # check whether this works for time dept.
    
  } else if(is.null(X) & const){ # constant, no covariates, param needs to be a vector
    # the sum across state and the sum across time points can be exchanged
    likelihood <- 0;
    for(i in 1:Tm1){
      # nhistT <- nhist[ ,i:(i+1)];
      nhistT =c(t(table(nhist[,i],nhist[,i+1])))
        likT <-  liktwostateest(param=param, n=dims[1],K=K,nhist=nhistT);
      likelihood <- likelihood + likT;
    }
  }  else if(!is.null(X) & const){ # constant, with covariates, param needs to be a vector
    # the sum across state and the sum across time points can be exchanged
    likelihood <- 0;
    for(i in 1:Tm1){  # i = 1
      # nhistT <- nhist[ ,i:(i+1)];# head(nhistT)
      nhistT =  nhist[,i:(i+1)];

      likT <-  liktwostateparest(param=param, n=dims[1],K=K,X=X,nhist=nhistT);
      likelihood <- likelihood + likT;
    }
  }
    
  return(likelihood);
  #  return(sum(likelihood));
}
# this is a vector of - loglikelihoods. best would be to optimise the sum of them?
# to do: test. optimise each time point separetely, vs. optimise the sum
# to do: 




like4states <- function(param, n=nrows(nhist),K=4,nhist, X=NULL, const=FALSE){
  
  dims <- dim(nhist);
  T <- dims[2];
  Tm1 = T-1;
  likelihood <- c();
  
  # fuer const == FALSE
  if(is.null(X) & !const){ # implementation without covariates
    
    for(i in 1:Tm1){
      nhistT <- nhist[ ,i:(i+1)];
      # nhistT =c(t(table(nhist[,i],nhist[,i+1])))
      if(is.null(dim(param))) paramT <- param else {
        paramT <- param[,i]   }

      likT <-  
      liktfourstateparest(param=paramT, n=dims[1],K=K,nhist=nhistT,X=NULL);
      likelihood <- c(likelihood, likT);
    }
    # likelihood <- sum(likelihood) ;  # check whether this works for time dept.
  } else if(!const & !is.null(X) ) { #time dependent, implementation with covariates
    
    for(i in 1:Tm1){
      # nhistT <- nhist[ ,i:(i+1)];
      nhistT =  nhist[,i:(i+1)];
      if(is.null(dim(param))) paramT <- param else {
        paramT <- param[,i]   }
      
      likT <-  liktfourstateparest(param=paramT, n=dims[1],K=K,X=X,nhist=nhistT);
      likelihood <- c(likelihood, likT);
    }
    #  likelihood <- sum(likelihood) ;  # check whether this works for time dept.
    
  } else if(is.null(X) & const){ # constant, no covariates, param needs to be a vector
    # the sum across state and the sum across time points can be exchanged # i = 1
    likelihood <- 0;
    for(i in 1:Tm1){
      nhistT <- nhist[ ,i:(i+1)];
      # nhistT =c(t(table(nhist[,i],nhist[,i+1])))
      likT <-  liktfourstateparest(param=param, n=dims[1],K=K,nhist=nhistT,X=NULL);
      likelihood <- likelihood + likT;
    }
  }  else if(!is.null(X) & const){ # constant, with covariates, param needs to be a vector
    # the sum across state and the sum across time points can be exchanged
    likelihood <- 0;
    for(i in 1:Tm1){  # i = 1
      # nhistT <- nhist[ ,i:(i+1)];# head(nhistT)
      nhistT =  nhist[,i:(i+1)];
      
      likT <-  liktfourstateparest(param=param, n=dims[1],K=K,X=X,nhist=nhistT);
      likelihood <- likelihood + likT;
    }
  }
  
  return(likelihood);
  #  return(sum(likelihood));
}



# liktwostateest: brute force algorithm without covariates
# get log likelihood for 2 state model
# for two states 1,2, nhist = c(n11,n12,n21,n22), ie nhist is a vector!
# n = sum(nhist), K = 2
# param = starting values on logit scale
liktwostateest <- function(param, n=sum(nhist),K=2,nhist){
  
  par = param.unpack1(param);
  # par <- c( 1- param[1], 1-param[2]);
  L0 <- nhist[1]*log(par[1]) +  (nhist[2])*log(1-par[1]) ;# slightly different from formula in book
  L1 <- nhist[3]*log(par[2]) +  (nhist[4])*log(1-par[2])  ;
return(-L0-L1);
}


# liktwostateparest: get log likelihood for 2 state model, incl design matrix with covariates X
# for categorical covariates!
# param for 2 states has length = ncol(X)*2
# param = (b00,b01, b0p, b10, b11, b1p)
# nhist is in this case (general form) a 2 dim matrix with first and 2nd time point for each subject
liktwostateparest <- function(param, n=nrow(nhist),K=2,nhist,X){
  
  # create a short form for the 2 state model history a la 11,12,21,22
  Z = apply(nhist, MARGIN=1, function(x) paste0(x, collapse="") ) # paste0(x[1], x[2])
  # assuming states are a and 2
 # n1 =  sum(Z ==  "11") +  sum(Z ==  "12") ;
 #   n2 =  sum(Z ==  "21") +  sum(Z ==  "2                                                                                                                                                       2") ;
  # assuming we have an intercept and dummy variables in X
  p <- ncol(X); # therefore we assume that we have p groups
  
  par = param.unpack2(param, X=X); # is a matrix with entry of delta11 and delta 21 for each subject
  # par <- c( 1- param[1], 1-param[2]);
  
  # if we have only categorical covariates,we could do this more efficient, but this version is more general
  L0 = 0; L1 = 0; # initialise
  for(i in 1:nrow(X)){
  L0 <- L0 + as.numeric(Z[i]== "11") *log(par[i,1]) +  as.numeric(Z[i]== "12")*log(1-par[i,1]) ;# slightly different from formula in book
  }
  
  for(i in 1:nrow(X)){
    L1 <- L1 + as.numeric(Z[i]== "21") *log(par[i,2]) +  as.numeric(Z[i]== "22")*log(1-par[i,2]) ;# slightly different from formula in book
  }
  # L1 <- nhist[3]*log(par[2]) +  (nhist[4])*log(1-par[2])  ;
  return(-L0-L1);
}


# get log likelihood for >2 state model, incl design matrix with covariates X
# for categorical covariates!
# param for >2 states has length = ncol(X)*2
# param = (b00,b01, b0p, b10, b11, b1p)
# nhist is in this case (general form) a 2 dim matrix with first and 2nd time point for each subject
liktfourstateparest <- function(param, n=nrow(nhist),K=4,nhist,X=NULL){
  
  if(is.null(X)) X <- matrix(1, ncol =1, nrow =n) ;
  # create a short form for the 2 state model history a la 11,12,21,22
  Z = apply(nhist, MARGIN=1, function(x) paste0(x[1], x[2]) )
  # assuming states are a and 2
  # n1 =  sum(Z ==  "11") +  sum(Z ==  "12") ;
  #   n2 =  sum(Z ==  "21") +  sum(Z ==  "2                                                                                                                                                       2") ;
  # assuming we have an intercept and dummy variables in X
  p <- ncol(X); # therefore we assume that we have p groups
  
  # here we get back a matrix with entries for each subject head(par); dim(par)
  if(K== 4)  par = param.unpackX(param=param, X=X, states=K);
    #par = param.unpack4(param=param, X=X); # is a matrix eith entry of delta11 and delta 21 for each subject
  # par <- c( 1- param[1], 1-param[2]); 
  if(K != 4)  par = param.unpackX(param=param, X=X, states=K);
  
  # to make life easier in the next step, create the full prob matrix (as vector)
  # here we assume that param.unpackX returns a matrix with K*(K-1) columns!
  parK <- matrix(0, ncol = K*K,nrow = n);
  for(i in 1:K){
    parK[ , (i-1)*(K) +  (1:(K-1)) ] <-   par[, (i-1)*(K-1) +  (1:(K-1)) ]
    parK[ , (i-1)*(K) +  K ] <- 1 - rowSums( par[, (i-1)*(K-1) +  (1:(K-1)) ] );
  }
  if(any(parK< 0 | parK >1) ) stop("some problem with the probability limits")
  # if we have only categorical covariates,we could do this more efficient, but this version is more general
  
  Lj = 0;
  # unhist <- unique(Z);
  # unhist <-paste0(rep(1:K, each=K-1), 1:(K-1))
  unhist <-paste0(rep(1:K, each=K), 1:(K))
  
  nj <- length(unhist); # so, we don't use the transitions that don't occur
  # L0 = 0; L1 = 0; # initialise
  for(j in 1:nj){
    for(i in 1:nrow(X)){
      Lj <- Lj + as.numeric(Z[i]== unhist[j]) *log(parK[i,j]) # +  as.numeric(Z[i]== unhist[j])*log(1-par[i,j]) ;# slightly different from formula in book
    }
  }
  
  return(-Lj);
}


# param.unpack1: for 2 state model without covariates / 
# transformation function for making sure that param values are probabilities betw. 0 an 1
# can also be used for larger number of states, where param has a different length
param.unpack1 <- function(param){
  p = exp(param)/(1+exp(param))
  return(p)
}

# param.unpack2: for 2 state model with covariates, general version
# returns a matrix with probs p11 and p21 as columns, rows for each subject
# (if X is just the intercept Xbeta0 and Xbeta1 are matrix versions of param in param.unpack1 with an entry for each individual)
# as final output, we can get transition probs for each parameter combinations, by selecting the index for individuals
# with a combination of parameters
param.unpack2 <- function(param,X){
  l = length(param)/2
  Xbeta0 = X%*%(param[1:l]);
  Xbeta1 = X%*%(param[(l+1):length(param)]);
  
  p01 = exp(Xbeta0)/(1+exp(Xbeta0))
  p11 = exp(Xbeta1)/(1+exp(Xbeta1))
  return( cbind(p01,p11 ));
}

#!!!: need to check the matrix operations!
# param.unpack4: for 4 state model with covariates, general version
#returns a matrix with vectors probs p1.,  p2., p3.,p4. in combinded columns, rows for each subject
#(if X is just the intercept Xbeta0 and Xbeta1 are matrix versions of param in param.unpack1 with an entry for each individual)
# as final output, we can get transition probs for each parameter combinations, by selecting the index for individuals
# with a combination of parameters
# param should include the full parameters, incl all transitions (even though one per row of the trans.matrix follows from the other row elements) 


# this function could be generalised to a larger number of transition states, 
# by looping through states and using lists as interim parameters

# param needs to have the following form 
# beta_s1_s2_. is the beta for the transition from s1 to s2, for X1,X2,X3
# beta_s1_s2_. = c(beta_s1_s2_x1,  beta_s1_s2_x2,beta_s1_s2_x3)  
# beta_s1_. = c(beta_s1_s1_. , beta_s1_s2_. , ... , beta_s1_smax_.) [max = number of states]
# param = c(beta_s1_., beta_s2_., ..., beta_smax_.)



# works with X with covariates states=K, X = NULL, not allowed 
# (for creating X, we would need to know n)
param.unpackX <- function(param,X, states,n ){
  # eg, for 4 states, length of param is assumed to be K*(K-1) with K = 4
  # if(is.null(X))
  lk = length(param)
  if(is.null(X)) X <- matrix(1, ncol =1, nrow =n) ; # needed below
  p <- ncol(X);
  # if (lk = 2) states = 2; if(lk = 12) states = 4;
 #  l = length(param)/states
  l = length(param)/(states) ; # i think that is the correct version,
  # this is the length for transition from one specific state
  
  
  # for the first state dim(Xbeta1) dim(output)
  # if X has just one column, then it is simply each value of p times X
  if(is.null(X) | p == 1){
    k = 1;
    paramx <- param[(k-1)*(l/p) +  1:(l/p)  ];
    p1_ = exp(paramx)/(1+ sum(exp(paramx)) );
    output <- p1_ ;
    for(j in 2:states){ 
      paramx <- param[( (j-1)*l+1):(j*l)]
      # Xbetaj <- matrix(0, nrow=nrow(X), ncol=l/p)
      # for(k in 1:(l/p) )    paramxt = paramx[(k-1)*(l/p) +  1:(l/p)  ] ;
      #  Xbetaj[,k] = X%*%(paramx[(k-1)*(l/p) +  1:(l/p)  ]   );
      # each of the below are vectors for each subject
      pj_ = exp(paramx)/(1+ sum(exp(paramx)) );
      
      output <- c( output,pj_ )
    }
    outputx = matrix(rep(output, times=nrow(X)), nrow =nrow(X),ncol= lk, byrow =TRUE)
    return(outputx)
  }
  

  # I need to work through this again, with an example to get the dimensions right
  # for the first state
  param1 <- param[1:(p*(states-1))  ]; # param for the first row of the trans.matrix= (p * (states-1))
  # now loop through each element - for each element, we have p params
  Xbeta1 <- matrix(0, nrow=nrow(X), ncol=(states-1)) ; # each element s1i of the transition matrix row gets a column
  
  # instead of this loop, we could do a matrix multiplication.
  # arrange each of the (states-1) times p elements of paramn1 into a p x (states-1) matrix
  parS <- matrix(param1, ncol = (states-1), nrow =p ,byrow=T );
  Xbeta1 = X%*%parS;
  
#  for(k in 1:(states-1) ){  
  #     paramx <- param1[(k-1)*(p) +  1:p  ];
  #    Xbeta1[,k] = X%*%(paramx );
  #   }
  # for p = 1 the dimension here is incorrect
   # each of the below are vectors for each subject
   p1_ = exp(Xbeta1)/(1+ rowSums(exp(Xbeta1)) ); # I think these should be rowSums, as we have one row for each subject
   output <- p1_ ; # this is the first state result
   
   for(j in 2:states){ # loop through each remaining state
     paramj <- param[(j-1)*(p*(states-1)) +  1:(p*(states-1))  ]; 
     
     parS <- matrix(paramj, ncol = (states-1), nrow =p  );
     Xbetaj = X%*%parS;
     
     #  Xbetaj <- matrix(0, nrow=nrow(X), ncol=(states-1))
     #   for(k in 1:(states-1) ){
     #     Xbetaj[,k] = X%*%(paramj[(k-1)*(states-1) +  1:p  ]   );
     #   } 
       # each of the below are vectors for each subject
     pj_ = exp(Xbetaj)/(1+ rowSums(exp(Xbetaj)) );
     
     output <- cbind( output,pj_ ) 
   }
   
  return( output);  #  as output we get the matrix with K*K-1 columns
}

#?rowSums

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




# getCI: get se, p values, and CI's (for beta)
# input: result of nlm function, incl $estimate and $hessian
# idea: put the nlm optimisation into this function
# this function works only for hessian full matrix
# otherwise use diaghessian
# if opt has time dependent estimates, all outputs are matrices, else vectors
getCI <- function(opt){
  if(is.null(opt$hessian)){ 
    se = 1/sqrt(opt$diaghessian) # se is a matrix with a diag for each time point
  } else se = 1/sqrt(diag(opt$hessian)) # standard errors
W = opt$estimate /se # Wald test statistic
pvalues = pchisq((W^2), df=1, lower.tail=FALSE) # get p values # check whether using abs is correct
# here I assume that the Wald test statistics is chi2(df=1) distributed

lowerCI = opt$estimate - 1.96*se # CI simply uses a t-distr, normality assumption on beta
higherCI = opt$estimate + 1.96*se # CI simply uses a t-distr, normality assumption on beta

return(list(se=se, p=pvalues, lowerCI=lowerCI,upperCI =higherCI ))
}



# alternative to param.unpack for > 2 states, but no X, only 1 transition!
# or average transition
transexpit <- function(param, states=4){
  out <- c();
  for(i in 1:states){
    par <- param[ (i-1)*(states-1) +  (1:(states-1))  ]
    out <- c(out, expit(par) );
  }
  return(out);
}



# transCI: transform CI estimate for beta into CI estimates for transition matrices
# CIS is assumed to be a list of se, p, lowerCI, upperCI
# categorical: whether covariates X are all categorical. if not (X ne NULL),
# returns lowerCI and upperCI with a row for each subject
# currently only implemented for 2 and 4 states transition matrices
transCI <- function(CIS,X, states, categorical = TRUE){
  
  if(!is.vector(CIS$se)){ # if each element is a matrix, each column = one time point
    output <- list();
    for (i in 1:ncol(CIS$se)){
      cit <- list();         cit$se <- CIS$se[,i];cit$p <- CIS$p[,i];
      cit$lowerCI <- CIS$lowerCI[,i];cit$upperCI <- CIS$upperCI[,i];
      output[[i]] <- transCI(cit,X, states=states, categorical = categorical)
    }
    return(output) # a list element for each time point
  }
  
 
  if( is.null(X) & states ==2) { # tested only for states == 2 
    
 #   if(is.vector(CIS$se)){
  return (list(lowerCI = param.unpack1(CIS$lowerCI),
               upperCI = param.unpack1(CIS$upperCI) ) )
    #  } else {
    #     output <- list();
    #     for (i in 1:ncol(CIS$se)){
    #       output[[i]]$lowerCI = param.unpack1(CIS$lowerCI[,i])
    #       output[[i]]$upperCI = param.unpack1(CIS$upperCI[,i])
    #     }
    #  return(output) # a list element for each time point
    # }
    # upperCI[is.na(upperCI)] <- 1; # from Inf/Inf
  } else if(is.null(X) & states > 2) {
    #  if(is.vector(CIS$se)){
    low = transexpit(param=CIS$lowerCI, states=states) 
    high = transexpit(param=CIS$upperCI, states=states) 
   # high[is.na(high)] <- 1;
    # sometimes what is the higher and lower CI is not correct
    # exception for NAs
    corr <- low < high;
    newlow <- newhigh <- rep(0, times = length(low));
    newlow[is.na(low)] <- low[is.na(low)]; newhigh[is.na(high)] <- high[is.na(high)];
    corr1 = corr; corr1[is.na(corr)] <- FALSE; corr2 = !corr; corr2[is.na(corr)] <- FALSE;
    newlow[corr1] <- low[corr1] ; newlow[corr2] <- high[corr2] ; 
    newhigh[corr1] <- high[corr1] ; newhigh[corr2] <- low[corr2] ;
    if(any(!corr)) warning("CI needed to be sorted by lower and higher limit")
    if(any(is.na(corr))) warning("NA values of CI could not be back-transformed")
    
    return(list(lowerCI = newlow,
                upperCI = newhigh)) ;
    
    # } else {
    #     output <- list();
    #   for (i in 1:ncol(CIS$se)){
    #      cit <- list();         cit$se <- CIS$se[,i];cit$p <- CIS$p[,i];
    #      cit$lowerCI <- CIS$lowerCI[,i];cit$upperCI <- CIS$upperCI[,i];
    #      output[[i]] <- transCI(cit,X, states=states, categorical = categorical)
    #    }
    #     return(output) # a list element for each time point
    #   }
  }
  
  if(!is.null(X)) {
    
    if(states == 2 ){
    lowerCIs = param.unpack2(CIS$lowerCI,X=X);
    upperCIs = param.unpack2(CIS$upperCI,X=X) ;
    } else if(states >2 ){
      lowerCIs = param.unpackX(CIS$lowerCI,X=X, states=states);
      upperCIs = param.unpackX(CIS$upperCI,X=X, states=states) ;
    }
    
    if(!categorical){
      return(list(lowerCI = lowerCIs, upperCI= upperCIs))
      # for continuous covariates, return CI's for each subject
    }
    
    # now find a representative of X for each group category
    # for this purpose, concatenate X
    Xconc <- apply(X,  MARGIN=1, function(x) paste0(x,collapse=""))
    unX <- unique(Xconc);
    lun <- length(unX)
    CISlower <- CISupper  <- list();
    for( i in 1:lun){
      CISlower[[i]] <- lowerCIs[which(Xconc == unX[i])[1], ]
      CISupper[[i]] <- upperCIs[which(Xconc == unX[i])[1], ]
      # it is possible that we need to resort, as above
      if(any(CISlower[[i]] > CISupper[[i]])) warning("CI needed to be sorted by lower and higher limit")
      
    }
    names(CISlower) <- names(CISupper) <-  unX;
    return(list(lowerCI = CISlower, upperCI= CISupper) )    
    # this is a list of lists, lowerCI contains lists for each different vector of X values
  }
}



# provide the output of param.unpackX (matrix with estimated transiitons row for each subject), and X
# return estimate for each group member, indicated by design matrix X
getparX <- function(parX, X){
Xconc <- apply(X,  MARGIN=1, function(x) paste0(x,collapse=""))
unX <- unique(Xconc);
lun <- length(unX)
estimate <-  list();
for( i in 1:lun){
  estimate[[i]] <- parX[which(Xconc == unX[i])[1], ]
}
names(estimate) <-  unX;
return(estimate  )    
}


# transEst: transform est estimate for beta into estimates for transition matrices
# est is assumed to be a  ?list of se, p, lowerCI, upperCI
# categorical: whether covariates X are all categorical. if not (X ne NULL),
# returns Estimate with a row for each subject
# currently only implemented for 2 and 4 states transition matrices
transEst <- function(est,X, states, categorical = TRUE){
  
  if(!is.vector(est)){ # if each element is a matrix, each column = one time point
    output <- list();
    for (i in 1:ncol(est)){
      output[[i]] <- transEst(est[,i],X, states=states, categorical = categorical)
    }
    return(output) # a list element for each time point
  }
  
  if( is.null(X) & states ==2) { # tested only for states == 2 
     return (param.unpack1(est) );
  } else if(is.null(X) & states > 2) {
    return(transexpit(est, states=states) ) ;
  }
  if(!is.null(X)) {
    
    if(states == 2 ){
      Estimate = param.unpack2(est,X=X);
    } else if(states >2 ){
      Estimate = param.unpackX(est,X=X, states=states);
    }
    
    if(!categorical){
      return(Estimate)
      # for continuous covariates, return CI's for each subject
    }
    
    # now find a representative of X for each group category
    # for this purpose, concatenate X
    Xconc <- apply(X,  MARGIN=1, function(x) paste0(x,collapse=""))
    unX <- unique(Xconc);
    lun <- length(unX)
    output <-  list();
    for( i in 1:lun){
      output[[i]] <- Estimate[which(Xconc == unX[i])[1], ]

    }
    names(output) <- unX;
    return(output )    
    # this is a list of lists, lowerCI contains lists for each different vector of X values
  }
}



