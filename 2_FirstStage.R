#####################################################################
#
#                 MCC HetPoll
#             Part 1: First stage, city level
#
#####################################################################

library(dlnm)
library(splines)

setwd("C:/Users/masselpl/Documents/Recherche/2020 - LSHTM/Projects/MCC-HetPoll")
load("Data/0_Data.RData")

#-------------------------------------
# Set up empty objects to save results
#-------------------------------------

coefall <- matrix(NA, nrow(cities), 1, dimnames = list(cities$city))
redall <- list()

vcovall <- vector("list", nrow(cities))
conv <- rep(NA, nrow(cities))

names(vcovall) <- names(conv) <- cities$city



#-------------------------------------
# Model parameters
#-------------------------------------

maxlagp <- 1
arglagp <- list(fun = "strata") # Equivalent to MA
argvarp <- list(fun = "lin")  # Nonlin in NEJM paper

cen <- 0

timedf <- 7

#-------------------------------------
# Loop on cities
#-------------------------------------

for(i in seq(length(dlist))) {
  cat(i,"")
  citydat <- dlist[[i]]
  
  # Construct crossbasis for temperature confounding
  cbt <- crossbasis(citydat$tmean, lag = 3, 
    arglag = list(fun = "strata"),
    argvar = list(fun = "bs", 
      knots = quantile(citydat$tmean, c(.1, .75, .9), na.rm = T))
  )
  
  # Construct crossbasis for PM2.5
  cbp <- crossbasis(citydat$pm25, lag = maxlagp, 
    argvar = argvarp, arglag = arglagp) 
  
  # Estimate the model
  model <- glm(death ~ cbp + cbt + dow + 
      ns(date, df = round(timedf * length(date) / 365.25)), 
    data = citydat, family = quasipoisson, na.action = "na.exclude")
  
  # Succesful estimation?
  conv[i] <- model$converged
  
  # Store results for 2nd stage meta-analysis
  redall[[i]] <- crosspred(cbp, model, cen = cen, at = 10) # Overall
  coefall[i,] <- redall[[i]]$allfit
  vcovall[[i]] <- redall[[i]]$allse^2
}

# Save results
save.image("Data/2_FirstStageResults.RData")

capture.output(print("Convergence failed"),
  cities[!conv, c("cityname","countryname")],
  file = "Results/2_ConvergenceFailure.txt"
)

