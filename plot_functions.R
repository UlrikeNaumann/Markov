
# update 24.04.2019: plotCameraCMR: option: wh = 1:max(5,nanimals) (to avoid error when nanimals < 5)

require(ggplot2, lib = libloc)
require(reshape2, lib = libloc)
require(plotly, lib = libloc)

# Working with Spatio-Temporal Data Classes, Pebesma (2012).
# Spatio-Temporal Statistics with R, S. 73
require("sp", lib = libloc)
require("spacetime", lib = libloc)
# require(crosstalk, lib = libloc)

# Multiple time-indexed spatial maps can be plotted from one long-format table
# using the functions facet_grid or facet_wrap in ggplot2 with time as a grouping
# variable.

# Hovmöller plot

# Hovmoller, E. 1949. The trough and ridge diagram. Tellus 1, 62-66. 

Hovmoeller <- function(time, latitude, longitude, x = NULL){
  # provide either time, latitude and longitude or a dataframe x containing these 3 columns
  par(mfrow=c(3,1));
  if(!is.null(x)){
    latitude <- x$latitude; longitude <- x$longitude; time <- x$time;
  }
  dist <- sqrt(diff(latitude)^2 + diff(longitude)^2); # need to generalise
  # prob: this only works if the animal was observed on two consecutive occacsions
  plot(time, latitude, type = "b");
  plot(time, longitude, type = "b");  
  plot(time[2:length(time)], dist, type = "b", main = "time vs. move"); 
  
}

# to plot both cameras and animal centres into one plot, with all camera radiuses and a few animal ranges
plotcameraCMR <- function(activitycentre, cameraCMR, nanimals = nrow(activitycentre), wh = 1:max(5,nanimals),  size=500, radius=30, 
                          cameraRadius = 25, main = "simulated species data"){
  plot(1, type="n", xlim=c(0, size), ylim=c(0, size ), ylab="length in m", xlab="width in m", 
       main = main)
  points(cameraCMR[,2:3] , pch=20, cex = 0.5, col="grey")
  points(activitycentre, col= "red", pch=20, cex = 0.2)
  
  points(activitycentre[wh,], col= "blue", pch=20, cex = 0.2)
  for (i in wh) {# plot the home ranges around the activity centres for each animal
    symbols(x=activitycentre[i,1], y=activitycentre[i,2], circles= radius, inches = FALSE, add = TRUE, fg = "blue")
  }
  for (i in 1:dim(cameraCMR)[1]) # plot the home ranges around the activity centres for each animal
    symbols(x=cameraCMR[i,2], y=cameraCMR[i,3], circles= cameraRadius, inches = FALSE, add = TRUE, fg = "grey")
  
}

# plotcameraCMR(activitycentre=activitycentre, cameraCMR =cameraCMR) 


# lag-0 (top) and lag-1 (bottom) empirical spatial covariance plots

# plot(mysecrdata) # this default function produces slightly strange plot

# interactive plots with time




# get sum of all captures for each proximity detector


# get sum of all animals caught at each prox detector


# get sum of all captures for each survey = time point. 
# here: no of captures = no of animals captured, since we simplify in the simulation one survey as one Event 
# (in practice, it could be over several days and thus could generate more then one capture per animal)


# get sum of all initial / first captures for each survey = time point.


# get sum of all recaptures for each survey = time point.



# model selection using, later add: model selection method = "AIC", "AICc",  "score", "lik"
# model selection for JS or CJS models from marked package,  class crm
# either compare two models or provide a list of models
cmr.selection <- function(crm1= NULL, crm2= NULL, crmlist=NULL, criterion = "AIC" ){
  select = 0;
  if (criterion == "AIC") use <- 5 else stop("method not yet defined");
  
  if (!is.null(crm1) & !is.null(crm2)){ 
    diff <- crm2$results$AIC - crm1$results$AIC;
   if (diff > 0) { bestaic <- crm1$results$AIC;
     select <- 1 } 
    else { select <- 2;
    bestaic <- crm2$results$AIC;
    }
  } else if(!is.null(crmlist) & is.list(crmlist)){
    listl <- length(crmlist)
    aiclist <- vector("list", listl)
    for(i in 1:listl) aiclist[[i]] <- crmlist[[i]]$results$AIC
    aicnum <- unlist(aiclist)
    select <- order(aicnum)[1]
    bestaic <- sort(aicnum)[1];
    diff <- sort(aicnum)[2] - bestaic;
  }
    # for option "AICc" (for small sample sizes), one would need k and sample size
    # aicc <- aic + (2*k^2 + 2*k)/(n - k - 1) # with n sample size,k param numbers
  # return(select) 
  c(select, bestaic, criterion, diff)
  # select is the number of the model with the smallest AIC
  # could also return the AIC itself, and the diff to the second best for AIC 
}  

# model selection for secr.fit models from secr package,  class secr
# either compare two models or provide a list of models with secrlist
# for model selection method = "AIC", "AICc" (small sample sizes) "score", "lik"
secr.selection <- function(secr1= NULL, secr2= NULL, secrlist=NULL, criterion = "AIC" ){
  select = 0;
  if (criterion == "AIC") use <- 5 else if(criterion == "AICc") use <- 6 else stop("criterion not yet defined");
  
  if (!is.null(secr1) & !is.null(secr2)){ 
      # diff <- summary(secr2)$AICtable[use] - summary(secr1)$AICtable[use];
      diff <- AIC(secr2, criterion =criterion)[use] - AIC(secr1, criterion =criterion)[use];
      # if (summary(secr1)$AICtable[5] < summary(secr2)$AICtable[5] )
      if (diff > 0){
        select <- 1 ;
        bestaic <- AIC(secr1, criterion =criterion)[use]
      } else { 
        select <- 2;
        bestaic <- AIC(secr2, criterion =criterion)[use]
      }
  } else if(!is.null(secrlist) & is.list(secrlist)){
    listl <- length(secrlist);
    aiclist <- vector("list", listl);
    for(i in 1:listl) aiclist[[i]] <- AIC(secr2, criterion =criterion)[use]
    # summary(secrlist[[i]])$AICtable[use];
    # alternatively check out usage of secrlist()
    aicnum <- unlist(aiclist);
    select <- order(aicnum)[1];
    bestaic <- sort(aicnum)[1]
    diff <- sort(aicnum)[2] - bestaic;
  }
  
  # return(select) 
  c(select, bestaic, criterion, diff)
  # select is the number of the model with the smallest AIC
  # could also return the crm model, possible 
}
  
# extractAIC(cjsmodel1$results)
