
# for bayesian - use other code .. see big PC

# likelihood closed CMR model, constant p, as in Pledger
# but for male and female, version 2
# y=c(nj,k) with

dat = simulatedat(n=1000)

# ------------------------------------------------------------------------------
# simulate capture history for male, female individuals, prop delta, with misclassification
# ------------------------------------------------------------------------------

# better: use version big PC
simulatedat <- function(n=1000,k=3, pm=0.3, pf=0.25,delta=0.6, cm=c(0.8, 0.1,0.1),
                        cf=c(0.2,0.7,0.1)){
  
  # simulate not capture as binomial
  nm= (delta)*n
  nf= (1-delta)*n
  nm0 = rbinom(n=nm,size=k,prob=1-pm)
  nf0 = rbinom(n=nf,size=k,prob=1-pf)
  # take out those with k not observed
  n0m = sum(nm0 != k)
  n0f = sum(nf0 != k)
  
  nm0x = nm0[nm0 != k]
  nf0x = nf0[nf0 != k]
  dat = matrix(0, ncol = 4, nrow = length(nm0x)+length(nf0x))
  dat[ 1:length(nm0x), 4] <- nm0x
  dat[(length(nm0x)+1):(length(nm0x)+length(nf0x)), 4] <- nf0x
  
  for(i in 1:length(nm0x)){
    dat[ i, 1:3] <- rmultinom(1,k-nm0x[i],cm)
     }
  for(i in 1:length(nf0x)){
    dat[ length(nm0x)+i, 1:3] <- rmultinom(1,k-nf0x[i],cf)
  }
  
  return(list(dat=dat, n0m = n0m, n0f=n0f))  
}

# check
dim(dat$dat)
all(rowSums(dat$dat) == 3)
# now mix up the data. 
datnew = dat$dat[order(rnorm(nrow(dat$dat))) , ]
dim(datnew)


# ------------------------------------------------------------------------------
# unified loglike version for CMR with gender uncertainty
# ------------------------------------------------------------------------------
# n - number of observed animals
# nj vector of length k with n1,...,nk, number of animals observed j times
# xj vector 1:k
# k number of sampling occassions
# x vector of (initial) values . length 8, log(n0), logit(pm) ,logit(pf), logit(delta_m)
#  cmm, cmf,  cfm, cff
# x[5] = logit(cmm), x[6] = logit(cmf/ (1-cmm)),x[7] = logit(cfm), x[8] = logit(cff/ (1-cfm))
# gender = c(n_m, n_f, n_u, n_not) matrix with columns not out of K=3 occs, of female, male, unknown, not caught
# gender needs to sum to k, for each , n nrows

# does not optimise correctly
# if plim not wanted, set it to 0. the min values for p_m and p_f
likgenderunv2 <- function(x, gender, plim = 0.05 ){
  
  # n0_m =  exp(x[1]) #change this
  #  n0_f =  exp(x[2]) # i suppose: n0_f = n0_m + gamma 
 # n0 = n0_m + n0_f
  # N_m = n0_m + n_m 
  n0 = exp(x[1])
   n = nrow(gender)
   #  n = sum(nj)
   N = n + n0;
   delta_m = 1/(1+ exp(- x[4] ))
   N_m = delta_m * N
   N_f = N -  N_m
   

  K = sum(gender[1,])
  if(! all(rowSums(gender) == K) ) stop(paste0("All observed rows need to sum to " , K))
  # c(cfm, cff, cfu)
  

  # N must be larger then n! and larger then 0
  # p must be between 0 and 1
  p_m = 1/(1+ exp(- x[2] ))  #  male detection prob
  if(p_m < plim) p_m <- plim
 #  p_f = 1/(1+ exp(- x[3]-x[4] ))  # x[4] is tau
  p_f = 1/(1+ exp(-x[3] ))
  if(p_f < plim) p_f <- plim
  
  # for this we need delta_notobs_m
  # consider using logs here to avoid underflow
 delta_notobs_m = (delta_m * (1-p_m)^K )/(delta_m * (1-p_m)^K + (1-delta_m) * (1-p_f)^K)
 # delta_notobs_m = exp( log(delta_m * (1-p_m)^K ) - log(delta_m * (p_m)^K + (1-delta_m) * (1-p_f)^K) )
 # alternative: 
 #  d_obs = (delta *(1-(1-pdetmale)^timep))/ ((delta *(1-(1-pdetmale)^timep)) + ((1-delta) *(1-(1-pdetfemale)^timep) ))
  
  n0_m = delta_notobs_m * n0
  n0_f = n0 - n0_m
  # we do not use n_m and n_f
  
  
 #  cmm, cmf,  cfm, cff
  cmm = 1/(1+ exp(- x[5] )) 
 #  x[6] = logit(cmf/ (1-cmm))
  cmf <- ( 1-cmm)/(1+ exp(- x[6] ))
#   if(1 - cmf - cmm)
  cmu = 1 - cmf - cmm 
  # we need to ensure that cmu >= 0
  # convention: eta[3] = 0. formally it is p_ij = exp(eta_ij)/sum_k(exp(eta_ik)) - we use a different parameterisation
  # cmu =  (1- p_m) + cmm+ cmf # this was incorrect
  cm  = c(cmm, cmf, cmu)
  # cfu =  (1- p_f) + cfm+ cff
  if(sum(cm) != 1){ # actually we ensured it with the transformation used
    warning("sum(cm)!=1; adjusted")
    # cmu = 1 - cmm- cmf
    tmp = (cmm + cmf + cmu);
    cmm = cmm/tmp; cmf = cmf/tmp;cmu = cmu/tmp;
    cm  = c(cmm, cmf, cmu)
  }  
 
  # this formula is incorrect
  cfm = 1/(1+ exp(- x[7] ))  
  cff <- ( 1-cfm)/(1+ exp(- x[8] ))
  cfu = 1 - cfm - cff 
  cf  = c(cfm, cff, cfu)
  if(sum(cf) != 1){
    warning("sum(cf)!=1; adjusted")
   # cfu = 1 - cfm- cff
    tmp = (cfm + cff + cfu );
    cfm = cfm/tmp; cff = cff/tmp;cfu = cfu/tmp;
    cf  = c(cfm, cff, cfu)
  }
  # the transformation used for cmm cmf needs to ensure that cmu can't be smaller then 0
  # the transformation used for cfm cff needs to ensure that cmu can't be smaller then 0
  
  # gender = c(n_m, n_f, n_u, n_not)

 # N_m = delta_m * N  # sum(n_m) + n0_m ;
 #  N_f =  N - N_m # sum(n_f) + n0_f ;
  #  k = length(nj)

  
  # check these 
  # if(p_m == 1 & (n_m == 0 & n0_m>0))   return(lik <- 10000)
  if(p_m == 1 & (n0_m> 0))   return(lik <- 10000)
 #  if(p_m == 0 & (n_m) >0)    return(lik <- 10000)
  if(p_f == 1 & (n0_f > 0))   return(lik <- 10000)
  # if(p_f == 0 & (n_f) >0)    return(lik <- 10000)
  if(p_f == 0 & p_m == 0 & n >0)    return(lik <- 10000)
 # if( any(rowsums(gender) != K ) ) stop("gender matrix misspecified")
  # rowSums(matrix(1:10, 2,5))
  
  # the sumnj part is independent of p, n0, delta etc,
  # ie a constant, that we don't really need for optimisation
  sumnj <- 0
  #  for (i in 1:k) sumnj <- sumnj + lfactorial(nj[i])
  part1 = lfactorial(N)-lfactorial(n0_m) - lfactorial(n0_f)-lfactorial(n) # - sumnj

  # part for L0
   # for non - observed animals
  #  lfactorial(K) -  lfactorial(K)
 #  part2 = n0_m * log((delta_m *(1-p_m)^(K)) ) + 
 #         n0_f * log((1-delta_m) *(1-p_f)^(K))
   part2 =  log(n0_m *(delta_m *(1-p_m)^(K))  + 
     n0_f*((1-delta_m) *(1-p_f)^(K))  )
  
  logliki = rep(0, times =n)
  # for observed animals
 
  p11 = p_m*cmm
  p12 = p_m*cmf
  p13 = p_m*cmu
  p21 = p_f*cfm
  p22 = p_f*cff
  p23 = p_f*cfu

  for(i in 1:n){
    l1 = delta_m * (p11^gender[i,1]) * (p12^gender[i,2])* (p13^gender[i,3])*(1-p_m)^(gender[i,4])
    l2 = (1-delta_m) * (p21^gender[i,1]) * (p22^gender[i,2])* (p23^gender[i,3])*(1-p_f)^(gender[i,4])
    
    lfac = lfactorial(K)-lfactorial(gender[i,1])-lfactorial(gender[i,2])-lfactorial(gender[i,3])-lfactorial(gender[i,4])
    
    logliki[i] = lfac +  log(l1+l2 ) # there is some danger of underflow here
    # as we only log after sum up l1 andl2
    
  }

  
  loglik = part1 +sum(logliki) + part2
  
  return(loglik)
}





# n - number of observed animals
# nj vector of length k with n1,...,nk, number of animals observed j times
# xj vector 1:k
# k number of sampling occassions
# x vector of (initial) values  log(n0_m),log(n0f), logit(pm) ,logit(pf), logit(delta_m)
# logit(delta_m), cmm, cmf,  cfm, cff
# gender = c(n_m, n_f, n_u, n_not) # not out of 3 occs, of female, male, unknown, not caught
# gender needs to sum to k

# this likelihood potentially still has errors!
likgenderunc <- function(x, y=c(nj,k), gender ){
  
  n0_m =  exp(x[1])
  n0_f =  exp(x[2]) # i suppose: n0_f = n0_m + gamma 
  n0 = n0_m + n0_f
  
  c(cfm, cff, cfu)
  
  delta_m = 1/(1+ exp(- x[5] ))
  # N must be larger then n! and larger then 0
  # p must be between 0 and 1
  p_m = 1/(1+ exp(- x[3] ))  #  male detection prob
  p_f = 1/(1+ exp(- x[3]-x[4] ))  # x[4] is tau
  
  #  cmu =  (1- p_m) + cmm+ cmf # this is not correct!
  cmu =  1-   cmm- cmf
   cm  = c(cmm, cmf, cmu)
   cfu =  1- - cfm- cff 
   cf  = c(cfm, cff, cfu)
   
   # the transformation used for cmm cmf needs to ensure that cmu can't be smaller then 0
   # the transformation used for cfm cff needs to ensure that cmu can't be smaller then 0
   
  # gender = c(n_m, n_f, n_u, n_not)
   k =y[  length(y)]
   nj = y[1:k]
  n = sum(nj)
  N = n + n0;
  N_m = delta_m * N  # sum(n_m) + n0_m ;
  N_f =  N - N_m # sum(n_f) + n0_f ;
#  k = length(nj)
  
  # check these - what is  n_m, n_f ?
#?  if(p_m == 1 & !(n_m[1] == 0 & n_m[2]==0))   return(lik <- 10000)
  #  if(p_m == 0 & sum(n_m) >0)    return(lik <- 10000)
  #  if(p_f == 1 & !(n_f[1] == 0 & n_f[2]==0))   return(lik <- 10000)
  #  if(p_f == 0 & sum(n_f) >0)    return(lik <- 10000)
  if( any(rowsums(gender) != k ) ) stop("gender matrix misspecified")
   # rowSums(matrix(1:10, 2,5))
  
  # the sumnj part is independent of p, n0, delta etc,
  # ie a constant, that we don't really need for optimisation
  sumnj <- 0
 #  for (i in 1:k) sumnj <- sumnj + lfactorial(nj[i])
  part1 = lfactorial(N)-(lfactorial(n0)) # - sumnj
  
  
  
  # part for L0
  # proportion of male in observed group for 1 occassion
#  delta_obs_m = (delta_m * p_m)/((delta_m * p_m) + ((1-delta_m)* p_f ) )  
 # delta_obs_f = 1 - delta_obs_m
   # the conditonal prob of G = male, given animal never obs
  a =   (delta_m * (1- p_m)^k)
   P_G1_y0 = (a)/(a + (1- delta_m)*((1- p_f)^k) )  
   P_G0_y0  = 1- P_G1_y0 
   suma = delta_m * ((1 - p_m)*P_G1_y0 )^k
   sumb = (1 - delta_m)* ((1 - p_f)*P_G0_y0 )^k
  part2 = (n0)* log(suma + sumb) # + log(w.n0)
  
  # determine how many different structures for ever observed 
  
  
  # L part for all the data avail
  gender_text = apply(gender,  margin = 1,paste0, collapse = "")
  # how many different 
  unique_gen <- unique(gender_text);
  which_gu <- rep(0, times = length(gender_text))
  unique_index <- rep(0, times = length(unique_gen))
  # calc P_G1_Y
  P_G1_Y = rep(0, times = length(unique_gen))
  part3 <- rep(0, times = length(unique_gen))
  
  for (i == 1:(length(unique_gen)) ) {
    which_gu[gender_text ==  unique_gen[i]   ] <- i
    unique_index[i] <- sum(gender_text ==  unique_gen[i])
    
    genderx <- gender[ (which_gu == 1)[1], ]
     # use delta_obs_m or delta_m
    A =  delta_m * (cmm^genderx[1]) * (cmf^genderx[2]) * (cmu^genderx[3])*( (1 - p_m)^genderx[4])
    B =  (1-delta_m) * (cfm^genderx[1]) * (cff^genderx[2]) * (cfu^genderx[3])*( (1 - p_f)^genderx[4])
    P_G1_Y[i] = A / (A+ B)
    
     A2 = delta_m *( (p_m*P_G1_Y[i])^sum(gender[1:3])  )* ((1-p_m)*P_G1_Y[i]) ^(gender[4])
     B2 = (1-delta_m) *( (p_f*(1-P_G1_Y[i]))^sum(gender[1:3])  )* ((1-p_f)*(1-P_G1_Y[i]))^(gender[4])

     # instead of summing over all n, use that we have counted how often each different combo appears
    part3[i]  <- unique_index[i] * log(A2 + B2)
  }  # delta_obs_m
 
  
  
  loglik = part1 + part2 + sum(part3)
  
  return(-loglik)
}


# ------------------------------------------------------------------------------
# likelihood closed CMR model, constant p, as in Pledger
# ------------------------------------------------------------------------------
# y=c(n,nj,xj,k) with
# n - number of observed animals
# nj vector of length k with n1,...,nk, number of animals observed j times
# xj vector 1:k
# k number of sampling occassions
# x vector of (initial) values log(n0), logit(p)
likelihoodclosedcmr <- function(x, y=c(n,nj,xj,k), weights =1){ # x = N,p
  # actually y could be shortened to nj as
  # n = sum(nj), k=length(nj) xj = 1:k, 
  
  # might be easier to estimate the number of non-observed
  n0 =  exp(x[1])
  # N must be larger then n! and larger then 0
  # p must be between 0 and 1
  p = 1/(1+ exp(- x[2] ))  #   s <- 1/(1+exp(-logits))
  n=y[1]
  N = n + n0;
  
  k=y[length(y)]
  nj = y[2:(1+k)]
  xj = y[(2+k):(1+2*k)]
  
  if(any(weights != 1)){
    # calculation of average weights over the respective k groups
    if(length(weights)== k){ w.k = weights ; } else if(length(weights)== 1) {
      w.k = rep(weights, times = k);} else 
    if(length(weights)== n){ # this works only for k = 3
      w.k = c(sum(weights[1:nj[1]])/nj[1] , sum(weights[(nj[1]+1):(nj[1]+nj[2])])/nj[2], 
              sum(weights[(nj[1]+nj[2]+1):n])/nj[3] );
    } else stop("weights should have length n")
    w.n0 = sum(w.k * nj )/n 
    # assume the average for the observed equals that of the unobserved
    } else {
     w.k = rep(1, times =k); w.n0 = 1;
      
  }
  
  if(p == 1 & !(nj[1] == 0 & nj[2]==0))   return(lik <- 10000)
  if(p == 0 & n>0)    return(lik <- 10000)
  J = length(nj)
  sumnj <- 0
  for (i in 1:k) sumnj <- sumnj + lfactorial(nj[i])
  part1 = lfactorial(N)-(lfactorial(n0)) - sumnj
  part2 = (n0)*k* log(1-p) + log(w.n0)
  
  # the 2 sums will give trouble for p = 0 or p = 1 
  sum1 = 0
  for (i in 1:J) sum1 = sum1 + nj[i]*xj[i]*log(p)
  
  value = 0
  for (i in 1:J) value <- 
    value+nj[i]*( lfactorial(k)-lfactorial(xj[i])-lfactorial(k-xj[i]))
  sum2 = 0
  for (i in 1:J) sum2 = sum2 + nj[i]*(k-xj[i])*log(1-p)
  
  part3 = sum1 + sum2 + sum(log(w.k) )
  
  lik = -(part1 + part2 + part3 + value)
  
  if (lik < 0)  {
    lik <- 10000
  }
  
  # this is - log likelihood
  return(lik)
}



# ------------------------------------------------------------------------------
# for 2 groups with weights. This version is not yet tested
# likelihood closed CMR model, constant p, as in Pledger
# ------------------------------------------------------------------------------
# g = 2 ( groups, not yet implemented for more then 2)
# y as in the 1 group case
# # n - number of observed animals
# # nj vector of length k with n1,...,nk, number of animals observed j times
# # xj vector 1:k
# k number of sampling occassions
# x vector of (initial) values log(n0_m), log(n0_f), logit(pmale), tau 
# with (logit(pfemale) = logit(pmale) + tau)

 mixedclosedCMR <- function(x, y=c(n,nj,xj,k),  weights, g = 2){
  # for now we assume g = 2 groups
# for more then 2 groups, the weights would need to be a matrix, that would
  # sum to 1 over all rows (=groups), columns 1,..., k for times captured
  n0_m =  exp(x[1])
  n0_f =  exp(x[2]) # i suppose: n0_f = n0_m + gamma 
  n0 = n0_m + n0_f
  # N must be larger then n! and larger then 0
  # p must be between 0 and 1
  p_m = 1/(1+ exp(- x[3] ))  #  male detection prob
  p_f = 1/(1+ exp(- x[3]-x[4] ))  # x[4] is tau
  
  # get - loglik for group 1
  likgroup1 = likelihoodclosedcmr(x=c(x[1], x[3] ) ,  y=c(n,nj,xj,k),weights = weights) 
  likgroup2 = likelihoodclosedcmr(x=c(x[2], log(p_f/(1-p_f) ) ) ,  y=c(n,nj,xj,k),weights = 1-weights) 
  
  return(likgroup1 + likgroup2 )
}
# p_f = 1/(1+ exp(- init[3]-init[4] ))
# log(p_f/(1 - p_f))



# ------------------------------------------------------------------------------
# this model works when group membership is known, at the moment only 
# for 2 groups. can adapt for more groups
# ------------------------------------------------------------------------------
likelihoodclosedcmr.g <- function(x, y, g=2, weights =1){ # x = N,p
  # actually y could be shortened to nj as
  # n = sum(nj), k=length(nj) xj = 1:k, 
  
  n0_m =  exp(x[1])
  n0_f =  exp(x[2]) # i suppose: n0_f = n0_m + gamma 
  n0 = n0_m + n0_f
  # N must be larger then n! and larger then 0
  # p must be between 0 and 1
  p_m = 1/(1+ exp(- x[3] ))  #  male detection prob
  p_f = 1/(1+ exp(- x[3]-x[4] ))  # x[4] is tau

  n = sum(y)
  
  k=length(y)/g
  xj = 1:k # the same for each group
  
  # n_g = (1:k)+seq(0,g*k, by=k)
  n_m=y[1:k]; n_f = y[(k+1):(2*k)] # this could be done into a matrix
  # not yet implemented for more then 2 groups
  N = n + n0;
  N_m = sum(n_m) + n0_m ;
  N_f = sum(n_f) + n0_f ;
  
  if(p_m == 1 & !(n_m[1] == 0 & n_m[2]==0))   return(lik <- 10000)
  if(p_m == 0 & sum(n_m) >0)    return(lik <- 10000)
  if(p_f == 1 & !(n_f[1] == 0 & n_f[2]==0))   return(lik <- 10000)
  if(p_f == 0 & sum(n_f) >0)    return(lik <- 10000)
  
  
  # easier: calculate lik for male, and lik for female 
  # and return sum
  # p=init[c(1,3)] ,y= c(sum(y[1:3]),y[1:3],1:k,k)
  # in order to get the weights, sort expected  prob male by number of captures / indiv.
  # , weights = c(0.8, 0.7,0.7)
  likmale = likelihoodclosedcmr(x=x[c(1,3)], y=c(sum(n_m),n_m,1:k,k) )
  likfemale = likelihoodclosedcmr(x=c(x[2], log(p_f/(1-p_f)) ), y=c(sum(n_f),n_f,1:k,k))
  liksum = likmale+likfemale
  if(g > 2) liksum = liksum  # + likelihoodclosedcmr(x=c(x[2], log(p_f/(1-p_f)) ), y=c(sum(n_f),n_f,1:k,k))
  # further adapt code for > 2 groups
 
   return(liksum)
  # ok, using this shortcut works!
  # the code below probably has some mistake in it
  
  J = length(n_m)
  sumnj <- 0
  for (i in 1:k) sumnj <- sumnj + lfactorial(n_m[i])+ lfactorial(n_f[i])
  part1 = lfactorial(N)-(lfactorial(n0_f)) -(lfactorial(n0_m)) - sumnj
  # ignore this part
  part1 = lfactorial(N)-(lfactorial(n0)) - sumnj
  
   part2 = (n0_m)*k* log(1-p_m)+(n0_f)*k* log(1-p_f)
  
  sum1 = 0
  for (i in 1:k) sum1 = sum1 + n_m[i]*xj[i]*log(p_m) + n_f[i]*xj[i]*log(p_f)
 # print(sum1)
  value = 0
  for (i in 1:k) {value <- 
    value+n_m[i]*( lfactorial(k)-lfactorial(xj[i])-lfactorial(k-xj[i])) +
          n_f[i]*( lfactorial(k)-lfactorial(xj[i])-lfactorial(k-xj[i])) 
  }
  sum2 = 0
  for (i in 1:k) {sum2 = sum2 + ( n_m[i]*(k-xj[i])*log(1-p_m) )+
   ( n_f[i]*(k-xj[i])*log(1-p_f))
  }
 # print(sum2)
  part3 = sum1 + sum2
  # print(c(part1, part2 , part3 , value))
  # print(c( (n0_m), (n0_f), (p_m), p_f ) )
  lik = -(part1 + part2 + part3 + value)
  # print(lik)
#  if (lik < 0)  {     lik <- 10000   }
  
  # this is - log likelihood
  return(lik)
}



# ------------------------------------------------------------------------------
# this model works when group membership is known, at the moment only 
# for 2 groups. can adapt for more groups
# ------------------------------------------------------------------------------
# here: n_m, n_f vorgegeben. allow instead to specify delta_m = N_m/N

# y should be a matrix with 5 columns. first 4 columns. no of observations as
# male, female, unknown, not observed, out of k, sum(y[1:4] = k)
# 5th columns. no of observations with this characteristics

likelihoodclosedcmr.unc <- function(x, y, g=2, weights =1){ # x = N,p
  # actually y could be shortened to nj as
  # n = sum(nj), k=length(nj) xj = 1:k, 
  
  n0_m =  exp(x[1])
  n0_f =  exp(x[2]) # i suppose: n0_f = n0_m + gamma 
  n0 = n0_m + n0_f
  # N must be larger then n! and larger then 0
  # p must be between 0 and 1
  p_m = 1/(1+ exp(- x[3] ))  #  male detection prob
  p_f = 1/(1+ exp(- x[3]-x[4] ))  # x[4] is tau
  
 #  n = sum(y)
  nj_m = y[,1] # no of male out of k observations
  nj_f = y[,2] # no of female out of k observations
  nj_u = y[,3] # no of unknown out of k observations
  nj_n = y[,4] # no of not observed out of k observations
  
  l = nrow(y) # number of different combinations in the dataset
  
  njs = y[,5] # how often this was observed - how many individuals
  n = sum(y[5])
  
  k=length(y)/g  # number of occassions
  xj = 1:k # the same for each group
  
  # n_g = (1:k)+seq(0,g*k, by=k)
 #  n_m=y[1:k]; n_f = y[(k+1):(2*k)] # this could be done into a matrix
  # not yet implemented for more then 2 groups
  N = n + n0;
  N_m = sum(n_m) + n0_m ; # actually, we don't know that
  N_f = sum(n_f) + n0_f ;
  
  delta_m = N_m/(N_m + N_f) # actually we don't know that
  
  if(p_m == 1 & !(n_m[1] == 0 & n_m[2]==0))   return(lik <- 10000)
  if(p_m == 0 & sum(n_m) >0)    return(lik <- 10000)
  if(p_f == 1 & !(n_f[1] == 0 & n_f[2]==0))   return(lik <- 10000)
  if(p_f == 0 & sum(n_f) >0)    return(lik <- 10000)
  
  
  # easier: calculate lik for male, and lik for female 
  # and return sum
  # p=init[c(1,3)] ,y= c(sum(y[1:3]),y[1:3],1:k,k)
  # in order to get the weights, sort expected  prob male by number of captures / indiv.
  # , weights = c(0.8, 0.7,0.7)
 # likmale = likelihoodclosedcmr(x=x[c(1,3)], y=c(sum(n_m),n_m,1:k,k) )
 #  likfemale = likelihoodclosedcmr(x=c(x[2], log(p_f/(1-p_f)) ), y=c(sum(n_f),n_f,1:k,k))
 #  liksum = likmale+likfemale
#  if(g > 2) liksum = liksum  # + likelihoodclosedcmr(x=c(x[2], log(p_f/(1-p_f)) ), y=c(sum(n_f),n_f,1:k,k))
  # further adapt code for > 2 groups
  
  # return(liksum)
  # ok, using this shortcut works!
  # the code below probably has some mistake in it
  
  J = length(n_m)
#  sumnj <- 0
#  for (i in 1:k) sumnj <- sumnj + lfactorial(n_m[i])+ lfactorial(n_f[i])
 #  part1 = lfactorial(N)-(lfactorial(n0_f)) -(lfactorial(n0_m)) - sumnj
  # ignore this part
  part1 = lfactorial(N)-(lfactorial(n0)) # - sumnj
  
  # never observed animals
  # conditional prob of gender for never observed animals
  
  pgender <-   (delta_m * (1 -p_m )^k)/ (delta_m * (1 -p_m )^k  + (1-delta_m) * (1 -p_f )^k)
  # part of the likelihood for never observed animals
  part2_m <- delta_m * ( (1- p_m) *pgender)^k 
  part2_f <- (1-delta_m) * ((1- p_f) *pgender)^k 
  # look again at this part. when we know n0_m and n0_f, then we can use the following

   part2 = (n0_m)*k* log(1-p_m)+(n0_f)*k* log(1-p_f)
  
   # for each combination of (n_{(i),m} +  n_{(i),f} +  n_{(i),u} +  n_{(i),0}
   # where $(n_{(i),m} +  n_{(i),f} +  n_{(i),u} +  n_{(i),0}) = K$.
   # max 20 comb - 1 (0,0,0) = 19 combinations max
   
   
   delta_obs_m <- (delta_m * p_m )/(delta_m * p_m + (1-delta_m) * p_f)
   delta_obs_f <- 1 - delta_obs_m
    likl <- rep(0, times = l)
   for(j in 1:l){
     A <-    delta_obs_m * ( cmm^nj_m  )* ( cmf^nj_f  )* ( cmu^nj_u  )* ( (1-p_m )^nj_n  )
     B <-   (1- delta_obs_m )* ( cfm^nj_m  )* ( cff^nj_f  )* ( cfu^nj_u  )* ( (1-p_f )^nj_n  )
     P_G.Y1 <- A / (A + B) # conditional prob of being male, given the observed combination
     P_G.Y0 <- 1 - P_G.Y1 # conditional prob of being male, given the observed combination
     likl.j_m <-   delta_m *( ((p_m*P_G.Y1)^(nj_m+nj_f+nj_u ) )*((1-p_m*P_G.Y1)^nj_n)   ) 
     likl.j_f <-   (1-delta_m )*( ((p_f*P_G.Y0)^(nj_m+nj_f+nj_u ) )*((1-p_f*P_G.Y0)^nj_n)   ) 
     
     likl[j] <- njs * log(likl.j_m + likl.j_f)
   }     
    
   part3 <- sum(c(likl) )   
#   sum1 = 0
#  for (i in 1:k) sum1 = sum1 + n_m[i]*xj[i]*log(p_m) + n_f[i]*xj[i]*log(p_f)

   #  value = 0
  # for (i in 1:k) {value <- 
  #  value+n_m[i]*( lfactorial(k)-lfactorial(xj[i])-lfactorial(k-xj[i])) +
  #  n_f[i]*( lfactorial(k)-lfactorial(xj[i])-lfactorial(k-xj[i])) 
  #}
#   sum2 = 0
#  for (i in 1:k) {sum2 = sum2 + ( n_m[i]*(k-xj[i])*log(1-p_m) )+
#    ( n_f[i]*(k-xj[i])*log(1-p_f))
#  }
 # part3 = sum1 + sum2
   lik = -(part1 + part2 + part3 )
  
  # this is - log likelihood
  return(lik)
}


# for t = 6
y = c(125, 3,0 , 76,4,0) # c(nj_male, nj_female)
init = c(log(5), log(5), log(0.1/(1-0.1)),log(0.08/(1-0.08)) -log(0.1/(1-0.1)))
# c(log(n0_m), log(n0_f), logit(pmale), tau )
likelihoodclosedcmr.g(x=init, y= y, g=2)# works
opt <- nlm(f=likelihoodclosedcmr.g, p=init ,y=y, iterlim=10000,
           hessian=TRUE, print.level=2, steptol = 1e-10)
?nlm
# now it works
opt$estimate
n0m = exp(opt$estimate[1]);n0m ; sum(c(125, 3,0))+ n0m # exp(gamma)
n0f = exp(opt$estimate[2]);n0f ; sum(c(76,4,0))+ n0m # exp(gamma)
pm = 1/(1+ exp(- opt$estimate[3] ));pm  #  male detection prob
pf = 1/(1+ exp(- opt$estimate[3]-opt$estimate[4] ));pf  # x[4] is tau

# ok, so the separate versions work well
opt2 <- nlm(f=likelihoodclosedcmr, p=init[c(1,3)] ,y= c(sum(y[1:3]),y[1:3],1:k,k), iterlim=10000,
           hessian=TRUE, print.level=2, steptol = 1e-10)
n0m = exp(opt2$estimate[1]);n0m ; sum(c(125, 3,0))+ n0m # exp(gamma)
pm = 1/(1+ exp(- opt2$estimate[2] ));pm  #  male detection prob
exp(-3.733812) # this should be p
likmale = likelihoodclosedcmr(x=c(log(n0m), log(1/(1-pm))), y=c(sum(n_m),n_m,1:k,k))


opt3 <- nlm(f=likelihoodclosedcmr, p=init[c(2,4)] ,y= c(sum(y[4:6]),y[4:6],1:k,k), iterlim=10000,
            hessian=TRUE, print.level=2, steptol = 1e-10)
n0f = exp(opt3$estimate[1]);n0f ; sum(y[4:6])+ n0m # exp(gamma)
pf = 1/(1+ exp(- opt3$estimate[2] ));pf  #  male detection prob

likmale = likelihoodclosedcmr(x=x[1:2], y=c(sum(n_m),n_m,1:k,k))
likfemale = likelihoodclosedcmr(x=x[3:4], y=c(sum(n_f),n_f,1:k,k))

opt <- optim(par=init,fn= likelihoodclosedcmr.g,y=y,  hessian =TRUE, method="BFGS") # starting with optim
n0m = exp(opt$par[1]);n0m ; sum(c(125, 3,0))+ n0m # exp(gamma)
n0f = exp(opt$par[2]);n0f ; sum(c(76,4,0))+ n0m # exp(gamma)
pm = 1/(1+ exp(- opt$par[3] ));pm  #  male detection prob
pf = 1/(1+ exp(- opt$par[3]-opt$par[4] )); pf  # x[4] is tau
# same

# for estimation ignoring gender
nj = c(125, 3,0) + c(76,4,0)
nj =c(76,4,0)
# following wikipedia: try out this code for estimation of the close model
# example: t = 11
nj=c(247 + 406 +1 , 35 + 22+ 1 + 1 + 4, 2 ); k = 3
xj = 1:k 
y = nj; sum(nj)
exposure =  factorial(k)/(factorial(xj)*factorial(k-xj))# the offset term
designmatx = cbind(rep(1, times=k), xj)

glm(y ~ offset(log(exposure)) + designmatx, family=poisson(link=log) )

model = glm(y ~ offset(log(exposure)) + 1+ xj, family=poisson(link=log) )
names(model)
model$coefficients
n0 = exp(model$coefficients[1]);n0 ; sum(nj)+ n0 # exp(gamma)
obsocc(pstar.beta(model$coefficients[2])) # estimate for p

# excellent - this is the result from using the Rcapture package

# ------------------------------------------------------------------------------
# In Poisson regression this is handled as an offset, where the exposure 
# variable enters on the right-hand side of the equation, but with a parameter estimate (for log(exposure)) constrained to 1. 
# ------------------------------------------------------------------------------


# theta = c(gamma, beta), with n0 = exp(gamma), ie gamma=log(n0) and pstar = fctx(beta)
# nj are the observed freqs of animals with j = 1:k observations, k = no of samples
likelihoodclosedpois <- function(theta, nj ){
  n = sum(nj)
  k=length(nj)
  xj = 1:k 
  offs =  factorial(k)/(factorial(xj)*factorial(k-xj))
  # offset = -log( factorial(k)/(factorial(xj)*factorial(k-xj)))
  
  #now: log(E(Y|x)) = log(y) = log(exposure) + theta' x
  # --> log( y/offs) = theta' x
   y = nj/offs
    designmatx = cbind(rep(1, times=k), xj)
    # calculate the loglikelihood
   #  temp = y*theta*x - exp(theta * x)
    # temp = nj*theta*designmatx - exp(theta * designmatx)
    theta.x = rowSums(  t(apply(X=designmatx,1 , FUN=function(x)(theta * x) ) ) )
    temp = y*theta.x 
    temp = temp - exp(theta.x )
   # ?apply
  lik = sum(temp)
  lik = -lik 
 # if (lik < 0)  {     lik <- 10000   } # do not use this here!
  # this is - log likelihood
  # However, the negative log-likelihood, ??? ??? ( ?? ??? X , Y), is a convex function,
  # and so standard convex optimization techniques such as gradient descent can be applied to find the optimal value of ??. 
  return(lik)
  # it turns out that -lik is a neg value! ok
  
}



# now it works!
beta =beta.pstar(pstar=everobs(0.1), k=3)
startval = c(log(1000), beta )
likelihoodclosedpois(theta=startval ,nj=nj) 
opt <- nlm(f=likelihoodclosedpois, p=startval ,nj=nj, iterlim=1000,
           hessian=TRUE, print.level=2, steptol = 1e-7)
?nlm
# gives an error for unknown reason
opt$estimate
n0 = exp(opt$estimate[1]);n0 ; sum(nj)+ n0 # exp(gamma)
obsocc(pstar.beta(opt$estimate[2]))


opt <- optim(par=startval,fn= likelihoodclosedpois,nj=nj,  hessian =TRUE, method="BFGS") # starting with optim
opt$par
# backtransform
?optim
n0 = exp(opt$par[1]);n0 ; sum(nj)+ n0 # exp(gamma)
obsocc(pstar.beta(opt$par[2]))
opt <- optim(par=opt$par,fn= likelihoodclosedpois,nj=nj,  hessian =TRUE, method="BFGS") # starting with optim
# 720*0.443 ; 2983*0.088
# calculate prob ever obs in k = 3 occas from p observed once and vice versa
 

# ------------------------------------------------------------------------------
# calculate lambda from theta and x
# ------------------------------------------------------------------------------
lambdafct <- function(theta, x=NULL,k=3){
  
  if(is.null(x)) x = cbind(rep(1, times=k),1:k)
  
  theta.x = rowSums(  t(apply(X=designmatx,1 , FUN=function(x)(theta * x) ) ) )
  
  lambda =    exp(theta.x)
  return(lambda)
}

# ------------------------------------------------------------------------------
# calculate theta  from lambda and x
# ------------------------------------------------------------------------------

thetafct <- function(lambda, x=NULL,k=3){
  
  if(is.null(x)) x = cbind(rep(1, times=k),1:k)
  k = dim(x)[1]
  
  loglambda = log(lambda) # = thetax
  sumj = sum(1:k)
  loglambda = k*gamma + sumj*beta
  
  # we don't know for either of these two values whether these are large enough
  exact = 0.1
  mat <- as.matrix(expand.grid(seq(0,10,exact), seq(-10, 10, exact)))
  dim(mat)   ; head(mat)   
  brutemethod = mat[,1] * mat[,2]
 
  # a = seq(0,10,exact) # log n0
  # b = seq(-10, 10, exact) # beta 
  # brutemethod = matrix(0, ncol = length(b), nrow = length(a))
  # for(i in 1:length(a)) brutemethod[i,] <- a[i] * b
  minimize = (brutemethod - loglambda)^2; # dim(minimize)
  wh = which( minimize == min(minimize))
  update = mat[wh, , drop = FALSE] # mat[13864,]
  # exp(update[,1]); xt= pstar.beta(update[,2], k=3); obsocc(xt[1],k=3);obsocc(xt[2],k=3) 
  # gamma = log(n0). problem. this is what we want to estimate
  # if we assume n0 >=1, gamma >= 0 up to 9.21 = log(10000)
  # beta.pstar(everobs(0.00001)); beta.pstar(everobs(0.99999)) . # -10: 10

  # no, we need to use loglambda_j = gamma + beta*j, for j = 1:k
}
# beta= beta.pstar(everobs(0.1)); gamma = log(2000);gamma;beta

# ------------------------------------------------------------------------------
# calculate pstar from beta
# ------------------------------------------------------------------------------
pstar.beta = function(beta, k=3){
  temp = 1/(1 + exp(beta))
  pstar = 1 - temp^k
  return(pstar)
}


# ------------------------------------------------------------------------------
# calculate beta from pstar
# ------------------------------------------------------------------------------
beta.pstar <- function(pstar, k=3){
  # pstar = 1 - (1 + exp(beta))^(-k)
  require(pracma)
  temp = nthroot(x=(1/(1-pstar)), n =k)
  
  beta <- log(temp - 1)
  return(beta)
} 
# beta.pstar(0.875) # eg for p = 0.5 -> pstar = 0.875 -> beta = 0
# pstar.beta(0)

# ------------------------------------------------------------------------------
# calculate prob ever obs from constant p and k
# ------------------------------------------------------------------------------
everobs = function(p, k=3){
  # this function works now for any positive k
 #  pstar = p^3 - 3*p^2 + 3*p
  temp = 1:k
    fct = function(p, m,k) {
      factorial(k)/(factorial(m)* factorial(k-m))*p^m * (1-p)^(k-m)
    }  
    pstar= sum( fct(p,temp,k) )
  return(pstar)
}
# pever = everobs(p=0.5)

# ------------------------------------------------------------------------------
# calculate prob at one occassion with constant p and k, from prob ever obs 
# ------------------------------------------------------------------------------
obsocc <- function(pever, exact = 0.001, k =3){
  # for 3 occassions
  # use a brute force method, up to 3 decimal places exact
  p.grid = seq(0,1,exact)
  
#   pstar = p.grid^3 - 3*p.grid^2 + 3*p.grid
  # for general k
  temp = 1:k
  # function for prob to be observed m out of k times with prob p
  fct = function(p, m,k) { 
    factorial(k)/(factorial(m)* factorial(k-m))*p^m * (1-p)^(k-m)
  } 
  p.x <- matrix(0, ncol = length(p.grid), nrow = k)
  for(i in 1:k)  p.x[i,] = fct(p.grid,temp[i],k)
  pstar= colSums( p.x ) # sum over all values from 1:k to obtain prob ever observed
  
  # now find which value of p is closest to pever
   tmp = (pever - pstar)
   # which tmp is closest to 0, some values will be neg, some pos
   index = which(tmp^2 == min(tmp^2))
  # index = which(pever == pstar ) # this only works if there is an exact match
  return( p.grid[index])
}
# obsocc(0.875)


# ------------------------------------------------------------------------------
# bayesian version closed model CMR. 
# MCMC using Gibbs sampling
# ------------------------------------------------------------------------------
# p = 1/(1+ exp( -x[2] ))
# in practise the function is used within an optimiser such as nlm, optim
# to obtain MLE estimates of p and n0, thus of N
# y=c(n,nj,xj,k) with
# n - number of observed animals
# nj vector of length k with n1,...,nk, number of animals observed j times
# xj vector 1:k
# k number of sampling occassions
# x vector of (initial) values log(n0), logit(p)
# prior for N can be either eg N(2000, sd = 1500) = normal
# or gamma: Gamma(4,0.002)

GibbsclosedCMR <- function(x, y=c(n,nj,xj,k), iter = 1000, nprior = "normal"){

  iteration = 1;
  
  n0 =  exp(x[1]) # initial value
  logn0 =x[1]
  
  # N must be larger then n! and larger then 0
  # p must be between 0 and 1
  logitp = x[2]
  p = 1/(1+ exp(- x[2] ))  #  initial value
  n=y[1]
  N = n + n0;
  
  k=y[length(y)]
  nj = y[2:(1+k)]
  xj = y[(2+k):(1+2*k)]
  
  while(iteration <= iter){
    
    # define hyperpriors
    mu_a  <- rnorm(1, mean = 0, sd = 1000)
    tau_a <- rgamma(1, 0.001, 0.001)

    # prior:
    logit_p_prior <- rnorm(mean = mu_a, sd = 1/tau_a)
    p_prior <- (1/(1+ exp(-logit_p_prior) ))
    
    # prior for N can be either eg N(2000, sd = 1500) = normal
    # or gamma: Gamma(4,0.002)
    if(nprior == "normal") {
      N_prior <- rnorm(1, mean =2000, sd = 1500)
    }  else N_prior <- rgamma(1, 4, 0.002)
        
##################################################################
# Gibbs-Sampler für p
##################################################################
    
    logn0.iter <- logn0[iteration]
    # some Gibbs sampler algorithm
    
    # likelihood * prior prob (p) * prior prob(N)
    newlogitp =   
    logitp <- c(logitp, newlogitp)
     
    logitp.iter <- logitp[iteration+ 1]

##################################################################
# Gibbs-Sampler für n0 (from n0 get N)
##################################################################
  
    # some algorithm
    newlogn0 = 
    logn0 <- c(logn0, newlogn0)
    
    iteration = iteration + 1
  }
  
  output = list(samplep = samplep, samplen0 = samplen0)
  return(output)
}  



# bayesian version. looking for a prior
var(c(500,5000)); var(500:(2*2700))

testprior = rnorm(1000,2000,sqrt(var(500:(2*2700)) ) )
testprior = rgamma(1000,4, 0.002)
b = 1/testprior
testprior = rgamma(1000,0.5,1)
hist(sort(testprior) )
hist(b)




plot(samplep, type = "l")
plot( samplen0, type = "l")
plot(1:length(samplep), samplep)
plot(1:length(samplen0), samplen0)

##################################################################
# this is a function to apply thinning, and get rid of the burnin
##################################################################
postcorrect.mcmc <- function(samplep, samplen0, burnin = 1000, thinning=10){
  
  samplepnew = samplep[-(1:burnin)]
  samplen0new = samplen0[-(1:burnin)]
  l = length(samplepnew)
  samplepnew  <- samplepnew[ seq(1,l, by=thinning)]
  l = length(samplen0new)
  samplen0new <- samplen0new[ seq(1,l, by=thinning)]
  
  output = list(samplep = samplepnew, samplen0 = samplen0new)
  return(output)
  
}


# check autocorrelation
par(mfrow = c(1,2))
acf(samplep)
acf(samplen0)



