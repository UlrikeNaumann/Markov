
# to do: test all of these functions!

# function to transform matrix Zg (capture history, with info on group observation for each indiv.) 
# into matrix Numg (matrix with no. of captures (1st col: no of captures), subsequ cols: no of observed group memberships among all recaptures)

ZgintoNumg <- function(Gs = NULL, Zg,g=3){
  # ie the first column of Num.g is for the number of occassions unequal to 0 # careful as this could cause issues  
  # g is the number of groups, ie not including "0"
  G = unique(c(Zg)); # if Gs is not provided, we just use the groups which were actually found 
  G = G[G!= "0"];
  
        
  # default assumption: Gs = c("m", "f", "u")
  if(is.null(Gs)) Gs = G;
  if(g != length(Gs)) stop("number of groups provided does not match Zg / Gs")
  Num.g = matrix(nrow = nrow(Zg), ncol= g+1) # + 1 for 0
  
  colnames(Num.g) = c("captures", Gs);
  Num.g[,1] = rowSums(Zg !="0") # ie the first column of Num.g is for the number of occassions unequal to 0 # careful as this could cause issues
  for(i in 1:g){
    
    Num.g[,1+i] = rowSums(Z.g == Gs[i]);
  }
  
  # wenn g == 3 & length(G) == 4 & paste(sort(G), collapse="") == "0fmu" # however i change this to 3 above by excluding "0" in G
 # dann nehme struktur an mit 3 gruppen male, female, unknown. 
  
  # kann aber sein, dass wir diese Struktur wollen, und eins der Gruppen gar nicht vorkommt. was dann?
  }


# ------------------------------------------------------------------------------
# Name: num.capt.gender
# Objective: To find the number and frequency of unique histories given multiple periods of (robus) capture history matrices with several groups
# Inputs: Z.g - group (gender) encounter histories
#         Num.g - alternative to Z.g. matrix with no of captures, no of captures per group for each group
#         T - number of years
#         K - no of secondary time points
#         g - no of groups. Default: 3
#         Gs - definition (names) of the groups
# Outputs: H - number of unique (gender) histories in matrix
#          nh - frequency of each unique history
#          Z.g - (if povided) copy of the input group (gender) encounter histories
#          Num.g - (if povided) copy of the (alternative) matrix with no of captures, no of captures per group for each group
#          unique.num - matrix containing encounter histories associated with each unique (gender) capture history 
#          Gs - definition (names) of the groups

num.capt.gender = function(Z.g=NULL, Num.g=NULL, T, K,g, Gs=NULL){
# Z.g = gender capt hist
# we are only interested in no. of captures, and no of male, female, unknown
# alternative: instead of Z.g, provide Num.g (matrix with 4 cols, captures, male, female, unknown)
## transform Num.g into a vector of character strings- of length 4  

# if Z.g is provided, and not Num.g, then first transform Z.g into Num.g
# the other way round is not possible
  if(!is.null(Z.g) &is.null(Num.g)) Num.g = ZgintoNumg(Gs = Gs, Zg=Zg,g=g);
  
char = apply(Num.g, 1, paste, collapse ="")
unique.num = unique(Num.g);
rownames(unique.num) <- NULL
unique.char = apply(unique.num, 1, paste, collapse = "");
H = length(unique.num[,1]);
nh = rep(0, H);
for(h in 1:H){
  test = unique.char[h] == char;
  nh[h] = sum(test);
}
return(list('H'= H, 'nh'=nh, 'Z.g'=Z.g, 'Num.g'=Num.g, 'unique.num'=unique.num, 'Gs' = Gs))
}

# ------------------------------------------------------------------------------
# Function to find number and frequency of unique histories for multiple periods
# ------------------------------------------------------------------------------
# Name: num.histories.multi
# Objective: To find the number and frequency of unique histories given multiple periods of capture history matrices
# Inputs: X - list of capture history matrices
#         Z - encounter histories
#         T - number of years
# Outputs: H - number of unique histories in matrix
#          nh - frequency of each unique history
#          unique.x - matrix containing the unique histories (full histories)
#          X.t - list containing capture matrices for each period
#          unique.z - matrix containing encounter histories associated with each unique capture history

num.histories.multi <- function(X, Z, T, K)  {
  
  # combine matrices
  X.full <- as.matrix(X[[1]])
  #if(K==1)  X.full <- as.matrix(X.full)
  for (t in 2:T)  {
    X.full <- cbind(X.full, as.matrix(X[[t]]) )
  }
  # remove all zero capture histories
  if(any(rowSums(X.full) == 0)){ # changes
    capture.0 <- which(rowSums(X.full) == 0)
    X.full <- X.full[-capture.0,]
    Z.full <- Z[-capture.0,]
  } else Z.full <- Z
  # find unique capture histories
  unique.x <- unique(X.full)  
  # make the unique histories into character strings for comparison
  unique.char <- apply(unique.x, 1, paste, collapse='')
  
  # set constants and storage variables
  H <- length(unique.x[,1]) 
  nh <- rep(0, H)  
  unique.z <- matrix(0, nrow=H, ncol=T)
  
  # change capture history matrix into vector of character strings
  x.char <- apply(X.full, 1, paste, collapse='')
  
  # loop over unique capture histories, counting how many times they are seen
  for (h in 1:H)  {
    
    # logical test if capture histories match the unique history being considered
    test <- unique.char[h] == x.char
    
    # sum up logical vector and store
    nh[h] <- sum(test)
    
    # store encounter history
    unique.z[h,] <- Z.full[which(test==1)[1],]
  }
  
  # split capture histories back into separate periods
  X.t <-  vector('list')
  
  unique.hist <- unique.x
  for (t in 1:T)  {
    X.t[[t]] <- as.matrix(unique.hist[,1:K])
    unique.hist <- as.matrix(unique.hist[,-(1:K)] )
  }
  
  # return unique histories, number and frequency
  return(list('H'=H, 'nh'=nh, 'unique.x'=unique.x, 'X.t'=X.t, 'unique.z'=unique.z))
}

# ------------------------------------------------------------------------------
# Function to find number and frequency of unique histories for multiple periods
# ------------------------------------------------------------------------------

# Name: num.histories.multi.gender
# Objective: To find the number and frequency of unique histories given multiple periods of capture history matrices
# Inputs: 
#         Z - encounter histories
#         T - number of years
#         K - vector with no. of secondary periods
# Outputs: H - number of unique histories in matrix
#          nh - frequency of each unique history
#          n.0 - frequency of 0 capture history
#          unique.x - matrix containing the unique histories (full histories)
#          unique.z.g - matrix containing encounter histories associated with each unique capture history

# assume. Z.g has "m", "f", "u" instead of 1 and 0 # Z.g is character matrix
# Z as normal - 0- 1 capture history
# don't input X, don't output X.t
# X.full = Z.g
# use Z to find all Zero capture histories
# leave out unique.z, instead use unique.z.g

num.histories.multi.gender <- function(Z, Z.g,T, K)  {
  
  # combine matrices
  # X.full <- as.matrix(X[[1]])
  #if(K==1)  X.full <- as.matrix(X.full)
  # for (t in 2:T)  {
  #   X.full <- cbind(X.full, as.matrix(X[[t]]) )
  # }
  X.full <- Z.g
  # remove all zero capture histories
  if(any(rowSums(Z) == 0)){ # changes
    capture.0 <- which(rowSums(Z) == 0) # which(rowSums(X.full) == 0)
    X.full <- X.full[-capture.0,]
    Z.full <- Z[-capture.0,]
  } else Z.full <- Z
  n.0 = length(capture.0)
  
  # find unique capture histories
  unique.x <- unique(X.full)  
  # make the unique histories into character strings for comparison
  unique.char <- apply(unique.x, 1, paste, collapse='')
  
  # set constants and storage variables
  H <- length(unique.x[,1])  # no of different histories
  nh <- rep(0, H)            # how often is each history seen
  unique.z.g <- matrix(0, nrow=H, ncol=T)
  
  # change capture history matrix into vector of character strings
  x.char <- apply(X.full, 1, paste, collapse='')
  
  # loop over unique capture histories, counting how many times they are seen
  for (h in 1:H)  {
    # X.full is used as Z.g, ie no of differing capture hisotires over Z.g is looked at!
    
    # logical test if capture histories match the unique history being considered
    test <- unique.char[h] == x.char
    
    # sum up logical vector and store
    nh[h] <- sum(test) # no of times seen for eah of the unique capture histories
    
    # store encounter history
    unique.z.g[h,] <- Z.full[which(test==1)[1],]
  }
  
  if(length(K) == 1) K = rep(K.times = T);
  
  # split capture histories back into separate periods
 #  X.t <-  vector('list')
  
  unique.hist <- unique.x
  for (t in 1:T)  {
   #  X.t[[t]] <- as.matrix(unique.hist[,1:K[t] ])
    unique.hist <- as.matrix(unique.hist[,-(1:K[t] )] )
  }
  
  # return unique histories, number and frequency
  return(list('H'=H, 'nh'=nh, 'n.0' = n.0,'unique.x'=unique.x, 'Z'=Z, 'unique.z.g'=unique.z.g))
}


