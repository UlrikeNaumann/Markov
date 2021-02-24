

require(chron, lib = libloc);


# now create the capture history!
# input variables: id number: vector, survey number: vector
# location: not yet used
# Output: surveyhist ,uniqueid , counthist

caphist <- function(id , surveyno, vt=NULL, location = NULL){
# at the moment we ignore the location
  
# first: how many animals
uniqueid <-   unique(id)
ni <- length(uniqueid)
#nsurveys <- length(unique(surveyno )) # max(surveyno)
nsurveys <-max(surveyno)
idnew <- rep(0 , times = length(id))
# to do: exclude vlaues of id which are NA or 0 

# give each idno a running number instead of the complicated id number
for (i in 1:ni){
  idnew[id == uniqueid[i]] <- i
}

encounterhist <- matrix(0,ncol = nsurveys , nrow = ni)
rownames(encounterhist ) <- 1:ni  # no of animals in total
colnames(encounterhist ) <- 1:nsurveys  # no of occasions
# now go through the data to input where any animal was detected

# now create encounter history
for(i in 1:length(idnew)){
  # the first value is the running no. of the id, the second is the survey at which the animal was encountered
  # if encountered twice in one survey, the value (already 1) will stay at 1
  
	encounterhist[idnew[i] , surveyno[i] ] <- 1 ; # dim(encounterhist) # tail(surveyno)

}


# could put in some check here, whether sum of all animals encountered at each survey is correct
colSums(encounterhist) # should all be = 4 in this case

rownames(encounterhist) <- as.character(uniqueid) # the row names should be again the original id numbers
# now I have them as rownames, I could again leave out the uniqueid in the output

# possibility: use the original dates of the survery in the col names, however these are multiple dates

if(!is.null(vt)){
  counthist <- matrix(0,ncol = length(vt) , nrow = ni)
  rownames(counthist) <- as.character(uniqueid) # the row names should be again the original id numbers
  colnames(counthist ) <- 1:length(vt);
  vtcount <- 1;
  for(j in 1:length(vt)){
    counthist[, j] <- rowSums( encounterhist[, vtcount:(vtcount + vt[j]-1) ]);
    vtcount <- vtcount + vt[j];
  }
  
  return(list(surveyhist = encounterhist,uniqueid = uniqueid, counthist= counthist))
} else return(
  list(surveyhist = encounterhist,uniqueid = uniqueid) )

}


# extension to capthist function
# now create the capture history!
# input variables: id number: vector, survey number: vector
# gender: vector of gender: can be either numeric or string:
# accepts: coding. "male adult" = 1 "female adult"=2  "unknown adult"=3 or 1,2,3
# ! all character/numeric, which are not "male adult" = 1 "female adult" = 2 are set to "unknown adult"=3
# location: not yet used
# Output: surveyhist ,uniqueid , genderhist , 
#         if(!is.null(vt)): counthist, primarygender

caphist.gender <- function(id , surveyno, gender=NULL, vt=NULL, location = NULL){
  # at the moment we ignore the location
  if(is.null(gender)) return(caphist(id , surveyno) )
  
  # make gender numeric, if it isn't. 
  # ! all character, which are not "male adult" = 1 "female adult" = 2 are set to "unknown adult"=3
 if(is.factor(gender)) gender <- as.character(gender)
 # if(!is.na(tryCatch(as.numeric(gender[1]))){ # improve on try
 #   gender.num <- as.numeric(gender);
 #   gender.num[gender.num!= 1 & gender.num!= 2] <- 3 # set all to unknown which are not 1 or 2
 # } 
 
   if(is.character(gender)) {
    gender.num <- rep(3, times = length(id))
    gender.num[gender == "male adult" ]<- 1
    gender.num[gender == "female adult" ]<- 2
   } 
 
 # head(gender.num); head(gender)
    
  # first: how many animals
  uniqueid <-   unique(id)
  ni <- length(uniqueid)
  #nsurveys <- length(unique(surveyno )) # max(surveyno)
  nsurveys <-max(surveyno)
  idnew <- rep(0 , times = length(id))
  # to do: exclude vlaues of id which are NA or 0 
  
  # give each idno a running number instead of the complicated id number
  for (i in 1:ni){
    idnew[id == uniqueid[i]] <- i
  }
  
  encounterhist <- matrix(0,ncol = nsurveys , nrow = ni)
  # rownames(encounterhist ) <- 1:ni  # no of animals in total
  rownames(encounterhist) <- as.character(uniqueid) # the row names should be again the original id numbers
  # now I have them as rownames, I could again leave out the uniqueid in the output
  colnames(encounterhist ) <- 1:nsurveys  # no of occasions
  
  genderhist <- matrix(0,ncol = nsurveys , nrow = ni)
  # rownames(encounterhist ) <- 1:ni  # no of animals in total
  rownames(genderhist) <- as.character(uniqueid) # the row names should be again the original id numbers
  # now I have them as rownames, I could again leave out the uniqueid in the output
  colnames(genderhist ) <- 1:nsurveys  # no of occasions
  
  # now go through the data to input where any animal was detected
  # possibly this was originally intended to be the location
  # withinhist <- matrix(0,ncol = nsurveys , nrow = ni)
  # rownames(withinhist) <- as.character(uniqueid) # the row names should be again the original id numbers
  # colnames(withinhist ) <- 1:nsurveys  # no of occasions
  
  # now create encounter history; idnew is numeric version of the id number vector
  for(i in 1:length(idnew)){
    # the first value is the running no. of the id, the second is the survey at which the animal was encountered
    # if encountered twice in one survey, the value (already 1) will stay at 1
    
    encounterhist[idnew[i] , surveyno[i] ] <- 1 ; # dim(encounterhist) # tail(surveyno)
    genderhist[idnew[i] , surveyno[i] ]    <- gender.num[i]
    
    # for withinhist, I count the number of encounters within one survey
    # withinhist[idnew[i] , surveyno[i] ] <- withinhist[idnew[i] , surveyno[i] ] + 1 ;
  }
  # to do: withinhist does not yet make any sense! - change into number of encounters per survey for 11 survey years
  
  if(!is.null(vt)){
    counthist <- matrix(0,ncol = length(vt) , nrow = ni)
    rownames(counthist) <- as.character(uniqueid) # the row names should be again the original id numbers
    colnames(counthist ) <- 1:length(vt);
    
    primarygender <- matrix(0,ncol = length(vt) , nrow = ni)
    rownames(primarygender) <- as.character(uniqueid) # the row names should be again the original id numbers
    colnames(primarygender ) <- 1:length(vt);
    
    vtcount <- 1;
    for(j in 1:length(vt)){
      counthist[, j] <- rowSums( encounterhist[, vtcount:(vtcount + vt[j]-1) ]);
      vtcount <- vtcount + vt[j];
    }
    # not yet correct!
    vtstart <- 1; for(i in 2:length(vt)) vtstart <- c(vtstart, vtstart[i-1] + vt[i-1])
    vtend <- vtstart + vt - 1;
    # go through each id number # vt = vtnew
    for(i in 1:ni){
      wh <- which( counthist[i, ]>0 )
      primarygender[i, wh] <- 3   # default all observed genders to unknown in primary survey hist
      vtid <- vt[j]
      for(j in (wh)){
        gtemp <- genderhist[i, vtstart[j]:vtend[j]] # sec. gender hist
        gtemp <- gtemp[gtemp != 0]
          # if( counthist[i,j] == 1 )  primarygender[i, j] <- unique(gtemp)
          if(length(unique(gtemp)) == 1) primarygender[i, j] <- unique(gtemp);
        # else if there are several different genders observed in this primary history, we have 
        # unknown gender, which is already the default which we used 
      }
    }  
    # for each id number, check which primary occassion > 0
    # if number of detections in the occassion = 1 --> copy gender into primarygender
    # if number of detections in the occassion > 1 --> are they all the same gender? --> copy gender into primarygender
    # else: set primarygender to 3 (unknown)
    # genderhist[, vtcount:(vtcount + vt[j]-1) ]
    
    return(list(surveyhist = encounterhist,uniqueid = uniqueid, 
                genderhist = genderhist,counthist= counthist, primarygender=primarygender))
  } else  return( list(surveyhist = encounterhist,uniqueid = uniqueid, genderhist = genderhist) )

  # could put in some check here, whether sum of all animals encountered at each survey is correct
  # colSums(encounterhist) # should all be = 4 in this case? why
  
  # possibility: use the original dates of the survery in the col names, however these are multiple dates
  
}


# encounterhist == surveyhist


# at is a vector with places to cut the capture history (eg 2 or 3 parts)
# eg for cutting a history with 9 time point into 3 parts, at =c(4,7)
# used for secondary history
# surveyhist = capture history
# at = points where the cut is. 
caphistcut <- function(surveyhist , at = round(ncol(surveyhist)/2))  {
  
  ll <- length(at);
  newhist <- matrix(0, ncol= length(at)+1, nrow=nrow(surveyhist));
  start <- 1
  end <- at[1]-1
  for(i in 1:(ll+1)){
    survpart <- surveyhist[,start:end, drop =FALSE]
    newhist[,i] <-    as.numeric(apply(survpart==1, any, MARGIN=1) );
    # specify start and end for next round
    start <- at[i];
    if (i ==ll) end <- ncol(surveyhist) else end <- at[i+1]-1;
  }
  return(newhist)
  
}

# reduce the genderhist to a smaller robust history version
# accepted codes are : 0 = not observed. 1 = adult male, 2 = adult female, 3 = adult unknown
caphistcut.gender <- function(genderhist , at = round(ncol(genderhist)/2))  {
  
  ll <- length(at);
  newhist <- matrix(0, ncol= length(at)+1, nrow=nrow(genderhist));
  start <- 1
  end <- at[1]-1
  for(i in 1:(ll+1)){
    survpart <- genderhist[,start:end, drop =FALSE]
    newhist[,i] <- 3* as.numeric(apply(survpart>0, any, MARGIN=1) );
    
    for(j in 1:nrow(survpart)){   # for each individual
      wh <- which( survpart[j, ]>0 )
      un <- unique( survpart[j, wh] )
      if(length(un)== 1) newhist[j,i] <- un
    }
    
#    newhist <-   3*( as.numeric(apply(survpart>0, which, MARGIN=1) ) ); # default to unknown gender
    # specify start and end for next round

    start <- at[i];
    if (i ==ll) end <- ncol(genderhist) else end <- at[i+1]-1;
  }
  return(newhist)
  
}
# to do: test this



# create small capture hist versions for lenghty robust models
# surveyhist = capture history
# at = vector of points where the cut is (same for each primary history)
# atall = points where the cut is (a list with 1 vector per primary history), 0 to fill up
# cut = a list of (atall) at vectors or atall lists with cut points (for different number of cuts)
caphist.min <- function(surveyhist, vt, at=NULL, atall=NULL, cut=NULL ){
  
  idnames <- rownames(surveyhist);
  if(is.null(at) & is.null(atall) & is.null(cut)) stop("Error. Provide cut values!") 
  
  l = 0; # running number old
  lnw = 0;  # running number new
  noprim <- length(vt);
  
  if(!is.null(at)) {
    ll <- length(at);
    vtnew <- rep( ll+1,times = noprim)
    newhist <- matrix(0, ncol= noprim*(ll+1), nrow=nrow(surveyhist));
  } else {
    if(length(atall)!= noprim) stop( paste("Error. atall needs to have the length", noprim) )
    # length of the list needs to be = noprim
    ll <- 0; # ll is defined differently in this case
    vtnew <- rep( NA,times = noprim)
    for(i in 1:noprim) { ll <- ll + length(atall[[i]])+1;
    vtnew[i]  <- length(atall)+1;  # print(sum(vtnew) == ll) # check
    }
    newhist <- matrix(0, ncol= ll, nrow=nrow(surveyhist));
    rownames(newhist) <- idnames;
  }
  sec.histlist <- list();  
  
  for(i in 1:noprim){
    # atx <- c(4,7);
    # if(i == 6) atx <- c(3,6);
      # sec.hist ist the secondary survey history for primary survey i

      sec.histlist[[i]] <- surveyhist[,(l+1):(l+vt[i])];
    if(is.null(at)) atx <- atall[[i]] else atx = at;
    
    # now create one combined capture history out of the 11 small ones
    # (ll+1) is the number of new parts in the secondary capt hist
    newhist[,(lnw+1):(lnw+vtnew[i])] <- caphistcut(sec.histlist[[i]] , at = atx);
    
     l <- l + vt[i] # running number for original capture history
     lnw <- lnw + vtnew[i]  # running number for new capture history
  }
  
  return(newhist);
  
  if(!is.null(cut)) {
    stop("not yet implemented")
    # check whether the following code would work - calls caphist.min from caphist.min
    noversions <- length(cut);
    output <- list()
    for(i in 1:noversions){
      if(is.list(cut[[i]])) {output[[i]] <- caphist.min(surveyhist=surveyhist, vt=vt,  atall=cut[[i]]) ; 
      } else {output[[i]] <- caphist.min(surveyhist=surveyhist, vt=vt, at=cut[[i]]) ;}
    }
    return(output) # is a list of shortened new capture histories
  }
  idnames
}  


# caphistcut(surveyhist=a, at=3)



# now create an M-array out of this data
# or what is done for capture recapture

# for SECR: instead of 1, record the location id of the respective detection location
# instead of 0, record NA
# Each entry (e.g. A9) records the detector at which a known animal (ID) was observed on the given # occasion (sample time). '.' indicates no detection. Each detector has known x-y coordinates.
# For 'proximity detectors' (eg cameras) multiple entries are possible on each occasion. 
# this only works if we have that information! eg, in which area was the animal
# Each detector has known x-y coordinates.
######################################################


# for SECR: instead of 1, record the location id of the respective detection location
# instead of 0, record NA
# surveyno needs to be a matrix with one column being a surveyid, the other the actual time
# Each entry location is a matrix which records the latitude/longitude at which a known animal (ID) was observed on the given 
# occasion (sample time). NA (alternative: '.') indicates no detection. 
# output includes either time code. or time distances between all survey occasions in days
# since this is planned to be used in robust CMR
# ie, in this function, we do not allow within survey multiple encounters --> make check on this

caphist.spat <- function(id , surveyno, location = NULL){
  require(chron);
  
  if(ncol(surveyno)<2) stop("error: surveyno needs to be matrix, with id and date");
  if(ncol(location)<2) stop("error: location needs to be matrix, with latitude and longitude");
  mydateTel <- as.Date(surveyno[,2], format = "%d-%m-%Y")
  
  timedist <- diff(surveyno[,2]  ) # is this function correct? calc the time dist vector
  # test <- c("01-02-2010" , "04-01-2011", "07-02-2012");test2 <- as.Date(test, format = "%d-%m-%Y"); diff(test2);
  
  
  # first: how many animals
  uniqueid <-   unique(id)
  ni <- length(uniqueid)
  nsurveys <- length(unique(surveyno[,1] ))
  surveymat <- as.data.frame(matrix(ncol=3, nrow= ncol(surveyno)) );
  surveymat[,1:2] <- surveyno;
  surveymat[1,3] <- 0; surveymat[2:ncol(surveyno),3] <- timedist;
  
  idnew <- rep(0 , times = length(id))
  # to do: exclude vlaues of id which are NA or 0 
  
  # give each idno a running number instead of the complicated id number
  for (i in 1:ni){
    idnew[id == uniqueid[i]] <- i
  }
  
  # nsurveys= 14; ni = 10;
  encounterhist <- matrix(NA,ncol = 2*nsurveys , nrow = ni) # for each survey occasion, 2 values
  rownames(encounterhist ) <- 1:ni  # no of animals in total
  colnames(encounterhist ) <- paste0(c("latitude","longitude") , rep(1:nsurveys, each=2))  # no of occasions

  # now create encounter history
  # in this function, we do not allow within survey multiple encounters
  for(i in 1:length(idnew)){
    tmp <- surveyno[i,1]
    # the first value is the running no. of the id, the second is the survey at which the animal was encountered
    # if encountered twice in one survey, the value (already 1) will stay at 1
    encounterhist[idnew[i] , c(tmp-1, tmp) ] <- location[i,] ;
  }
  
  rownames(encounterhist) <- as.character(uniqueid) # the row names should be again the original id numbers
  # now I have them as rownames, I could again leave out the uniqueid in the output
  
  list(surveyhist = encounterhist,uniqueid = uniqueid, surveymat = surveymat)
  # it is probably better to return also the matrix of 01... histories
}

# matrix(1:10, nrow=5,ncol=2 , byrow=TRUE)

# ideally: input: an output of capthist.spat, including info on which animal was caught 
# multiple times, output info must include time differences between survey occasions
# input: a matrix consisting of ID_animal, time (day-month-year), location
disttime <- function(surveyhist, surveymat, nsurveys = ncol(surveyhist)/2){
  
  
  capthist <- !is.na(as.data.frame(surveyhist))[,(1:nsurveys)*2];
  tmp = rowSums(!is.na(as.data.frame(surveyhist)));
  keep <- which( tmp/2 > 1); # ie having been caught at least 2 times
  l = length(keep);
  recaps = tmp/2 -1; # total counts of recaptures matrix  
  recaps[recaps<0] <- 0   # we are only interested in the positive values
  # idea: build a matrix with id, dist time, dist space
  # build one matrix with distances to first capture, and distance to the prev capture
  dist1 <- matrix(NA, ncol = 5, nrow = sum(recaps) ) ;
  recatches <- unique(recaps[recaps!=0]);
  
  rownames(dist1) <-  as.character(rep( rownames(capthist), times = recaps));
  dist1[,1] <- rep( rownames(capthist), times = recaps); # bei recaps = 0, faellt das Tier raus
  
  j = 1;
  for(i in 1:l){
    survhisti <- surveyhist[keep[l], ] ;# surveyhist for this specific animal
    # best resort them to a (recaps[i]+1) x 2 matrix
    survhistform <- matrix(survhisti, ncol = 2, nrow = recaps[i]+1, byrow=TRUE)

    capthisti <- capthist[keep[l],];
    times = surveyid[capthisti,]; # which were the survey occasions where this animal was caught
    k = j + recaps[i];
    dist1[j:(k-1),2] <- times[-1,2] - times[1,2] # time differences to the first occasion
    dist1[j:(k-1),3] <- sqrt ( rowSums( matrix( (survhistform[-1,] - survhistform[1,])^2,ncol = 2, nrow = recaps[i] ) ) );
    dist1[j:(k-1),4] <- diff(times[,2]); # does the diff function work on dates?
    dist1[j:(k-1),5] <- sqrt (rowSums(diff(survhistform)^2) ); # does the diff function work on dates?
    
    }
  
  # build one matrix with distances only for recaptures within the same survey
  # rowSums(matrix(1:2,ncol=2))
  # build one matrix only with distance between surveys, to the first recapture within each survey
  return(dist1) # differences to first capture, and to last previous
  
}


# make a summary of sightings, no of sightings per survey period
# type is a vector of either  " " or "sighting", surveyno is vector of surveynos 
sightings <- function(type, surveyno, biogenderage = NULL){
  type <- as.character(type);
  
  survsight = surveyno[type== "sighting"];
  if(!is.null(biogenderage)){
    genderage <- biogenderage[type== "sighting"];
    return(table(survsight, genderage))
  }
  
  return(table(survsight))
  
}
