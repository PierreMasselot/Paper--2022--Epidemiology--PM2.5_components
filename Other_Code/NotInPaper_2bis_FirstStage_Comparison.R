#####################################################################
#
#                 MCC HetPoll
#             Part 1: First stage, city level
#
#####################################################################

library(dlnm)
library(splines)
library(mgcv)
library(viridis)

setwd("C:/Users/masselpl/Documents/Recherche/2020 - LSHTM/Projects/MCC-HetPoll")
load("Data/0_Data.RData")

#-------------------------------------
# Set up empty objects to save results
#-------------------------------------

qaicres <- gcvres <- RRres <- matrix(NA, nrow(cities), 7)
rownames(qaicres) <- rownames(gcvres) <- rownames(RRres) <- 
  cities$city
colnames(qaicres) <- colnames(gcvres) <- colnames(RRres) <- 
  c("NEJM", "Unconstrained", 
  "Nonlinear", "Both", "NoHum", "Temp21", "Before2003")


#-------------------------------------
# Loop on cities
#-------------------------------------

for(i in seq(length(dlist))) {
  cat(i,"")
  citydat <- dlist[[i]]
  
  if (is.null(citydat$rhum) || (mean(is.na(citydat$rhum)) > .95)){
    form <- death ~ cbp + cbt + dow +
      ns(date, df = round(7 * length(date) / 365.25))
  } else {
    form <- death ~ cbp + cbt + dow + ns(rhum, df = 3) +
      ns(date, df = round(7 * length(date) / 365.25))
  }
  
  # Exactly as in the NEJM paper
  cbt <- crossbasis(citydat$tmean, lag = 3, arglag = list(fun = "strata"),
    argvar = list(fun="ns", df = 6)
  )
  cbp <- crossbasis(citydat$pm25, lag = 1, 
    argvar = list(fun = "lin"), arglag = list(fun = "strata")) 
  model <- glm(form, 
    data = citydat, family = quasipoisson, na.action = "na.exclude")
  qaicres[i, 1] <- QAIC(model)
  gcvres[i, 1] <- mean((residuals(model)/
    (1 - summary(model)$df[3]/length(model$y)))^2, na.rm = T)
  cp <- crosspred(cbp, model, cen = 0, at = 10)
  RRres[i, 1] <- cp$allfit
  
  # Unconstrained lag for PM2.5
  cbp <- crossbasis(citydat$pm25, lag = 1, 
    argvar = list(fun = "lin"), arglag = list(fun = "integer")) 
  model <- glm(form, 
    data = citydat, family = quasipoisson, na.action = "na.exclude")
  qaicres[i, 2] <- QAIC(model)
  gcvres[i, 2] <- mean((residuals(model)/
    (1 - summary(model)$df[3]/length(model$y)))^2, na.rm = T)
  cp <- crosspred(cbp, model, cen = 0, at = 10)
  RRres[i, 2] <- cp$allfit
  
  # nonlinear relationship for pm25
  cbp <- crossbasis(citydat$pm25, lag = 1, 
    argvar = list(fun = "bs", knots = quantile(citydat$pm25, .25, .75, na.rm = T)), 
    arglag = list(fun = "strata")) 
  model <- glm(form, 
    data = citydat, family = quasipoisson, na.action = "na.exclude")
  qaicres[i, 3] <- QAIC(model)
  gcvres[i, 3] <- mean((residuals(model)/
    (1 - summary(model)$df[3]/length(model$y)))^2, na.rm = T)
  cp <- crosspred(cbp, model, cen = 0, at = 10)
  RRres[i, 3] <- cp$allfit
  
  # Both unconstrained and nonlinear
  cbp <- crossbasis(citydat$pm25, lag = 1, 
    argvar = list(fun = "bs", knots = quantile(citydat$pm25, .25, .75, na.rm = T)), 
    arglag = list(fun = "integer")) 
  model <- glm(form, 
    data = citydat, family = quasipoisson, na.action = "na.exclude")
  qaicres[i, 4] <- QAIC(model)
  gcvres[i, 4] <- mean((residuals(model)/
    (1 - summary(model)$df[3]/length(model$y)))^2, na.rm = T)
  cp <- crosspred(cbp, model, cen = 0, at = 10)
  RRres[i, 4] <- cp$allfit
    
  # Without rhum control
  cbt <- crossbasis(citydat$tmean, lag = 3, arglag = list(fun = "strata"),
    argvar = list(fun="ns", df = 6)
  )
  cbp <- crossbasis(citydat$pm25, lag = 1, 
    argvar = list(fun = "lin"), arglag = list(fun = "strata")) 
  model <- glm(death ~ cbp + cbt + dow +
      ns(date, df = round(7 * length(date) / 365.25)), 
    data = citydat, family = quasipoisson, na.action = "na.exclude")
  qaicres[i, 5] <- QAIC(model)
  gcvres[i, 5] <- mean((residuals(model)/
    (1 - summary(model)$df[3]/length(model$y)))^2, na.rm = T)
  cp <- crosspred(cbp, model, cen = 0, at = 10)
  RRres[i, 5] <- cp$allfit
  
  # With more temperature control
  cbt <- crossbasis(citydat$tmean, lag = 21, 
    arglag = list(fun = "ns", knots = logknots(21, nk = 3)),
    argvar = list(fun = "bs", knots = quantile(citydat$tmean, c(.1, .75, .9), na.rm = T))
  )
  cbp <- crossbasis(citydat$pm25, lag = 1, 
    argvar = list(fun = "lin"), arglag = list(fun = "strata")) 
  model <- glm(form, 
    data = citydat, family = quasipoisson, na.action = "na.exclude")
  qaicres[i, 6] <- QAIC(model)
  gcvres[i, 6] <- mean((residuals(model)/
    (1 - summary(model)$df[3]/length(model$y)))^2, na.rm = T)
  cp <- crosspred(cbp, model, cen = 0, at = 10)
  RRres[i, 6] <- cp$allfit
  
  # With more data (i.e. prior to 2003)
  cbt <- crossbasis(dlist_orig[[i]]$tmean, lag = 3, arglag = list(fun = "strata"),
    argvar = list(fun="ns", df = 6)
  )
  cbp <- crossbasis(dlist_orig[[i]]$pm25, lag = 1, 
    argvar = list(fun = "lin"), arglag = list(fun = "strata")) 
  model <- glm(form, data = dlist_orig[[i]], family = quasipoisson, 
    na.action = "na.exclude")
  qaicres[i, 7] <- QAIC(model)
  gcvres[i, 7] <- mean((residuals(model)/
    (1 - summary(model)$df[3]/length(model$y)))^2, na.rm = T)
  cp <- crosspred(cbp, model, cen = 0, at = 10)
  RRres[i, 7] <- cp$allfit
}

ncountries <- length(unique(cities$country))
pal <- viridis(nrow(countries))
pchs <- rep_len(15:18, nrow(countries))

ncomp <- ncol(qaicres) - 1

# Plot QAIC
x11(width = 15, height = 15)
par(mfrow = n2mfrow(ncomp + 1, asp = 1.5))
for (i in 1:ncomp){
  plot(qaicres[,1], qaicres[,i + 1]/qaicres[,1], 
    pch = pchs[as.numeric(droplevels(cities$country))], 
    xlab = "NEJM", ylab = sprintf("%s / NEJM", colnames(qaicres)[i + 1]),
    col = pal[as.numeric(droplevels(cities$country))])
  abline(h = 1)
}
plot(0,0, col = "white", axes = F, xlab = "", ylab = "")
legend("topright", legend = unique(cities$countryname), col = pal, pch = pchs,
  ncol = 3, bty = "n", cex = .8, xpd = T)
dev.print(png, filename = "Results/2bis_QAIC.png", 
  units = "in", res = 100)

# Plot GCV 
# (warning: theoretically better for squared-error loss, 
#   i.e. data approximately gaussian, which is usually the case here when the
#   mortality count is high enough)
x11(width = 15, height = 15)
par(mfrow = n2mfrow(ncomp + 1, asp = 1.5))
for (i in 1:ncomp){
  plot(gcvres[,1], gcvres[,i + 1]/gcvres[,1],  
    pch = pchs[as.numeric(droplevels(cities$country))],
    xlab = "NEJM", ylab = sprintf("%s / NEJM", colnames(qaicres)[i + 1]),
    col = pal[as.numeric(droplevels(cities$country))])
  abline(h = 1)
}
plot(0,0, col = "white", axes = F, xlab = "", ylab = "")
legend("topright", legend = unique(cities$countryname), col = pal, pch = pchs,
  ncol = 3, bty = "n", cex = .8, xpd = T)
dev.print(png, filename = "Results/2bis_GCV.png", 
  units = "in", res = 100)
  
# Plot RR at 10 ug/m3 
x11(width = 15, height = 15)
par(mfrow = n2mfrow(ncomp + 1, asp = 1.5))
for (i in 1:ncomp){
  plot(RRres[,1], RRres[,i + 1],  
    pch = pchs[as.numeric(droplevels(cities$country))],
    xlab = "NEJM", ylab = colnames(qaicres)[i + 1],
    col = pal[as.numeric(droplevels(cities$country))])
  abline(a = 0, b = 1)
}
plot(0,0, col = "white", axes = F, xlab = "", ylab = "")
legend("topright", legend = unique(cities$countryname), col = pal, pch = pchs,
  ncol = 3, bty = "n", cex = .8, xpd = T)
dev.print(png, filename = "Results/2bis_RR10.png", 
  units = "in", res = 100)