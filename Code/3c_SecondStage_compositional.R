#####################################################################
#
#                 MCC HetPoll
#             Part 2: Second stage meta analysis
#                using a compositional data approach
#
#####################################################################

library(mixmeta)
library(compositions)
library(zCompositions)

load("Data/2_FirstStageResults.RData")

#-------------------------------------
#   Preparing constituent data
#-------------------------------------
spec_inds <- grep("PM25", colnames(dlist_spec[[1]]))
spec_names <- sapply(strsplit(colnames(dlist_spec[[1]])[spec_inds], "_"),
  "[", 2)
spec_pal <- c(2, 6, 4, 1, 3, 5, 7)

# We take the mean composition for the city
mean_comp <- sapply(dlist_spec, function(x){
  # Zero value imputation as in Martín-Fernández et al. (2003)
  imp <- multRepl(x[,spec_inds], label = 0, dl = rep(1e-5, 7))
  # Transformation to compositional object
  xc <- acomp(imp)
  # Compositional mean (pass through the irl transformation)
  mean(xc)
})
mean_comp <- t(mean_comp)
p <- ncol(mean_comp)

#-------------------------------------
#   Naive meta-regression
#-------------------------------------

# We do it one at a time since the design matrix is singular (sum to one)
res_naive <- matrix(NA, p, 4, 
  dimnames = list(spec_names, c("coef", "lo", "up", "sig")))
for (j in seq_len(p)){
  metamod <- mixmeta(coefall ~ mean_comp[,j], vcovall, 
    random = ~ 1|country/city,
    data = cities, method = "reml", subset = conv)
  print(metamod)
  
  # We store only the first coef
  res_naive[j, 1] <- coef(metamod)[2]
  se_naive <- sqrt(diag(vcov(metamod))[2])
  res_naive[j, 2] <- res_naive[j, 1] - 1.96 * se_naive
  res_naive[j, 3] <- res_naive[j, 1] + 1.96 * se_naive
  res_naive[j, 4] <- res_naive[j, 2] > 0 | res_naive[j, 3] < 0
}

# Forestplot
bord <- order(res_naive[,1], decreasing = T)
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(res_naive[bord,1], -seq_len(p), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(res_naive[,2:3]))
segments(res_naive[bord,2], -seq_len(p), res_naive[bord, 3], -seq_len(p), 
  col = "darkgrey", lwd = 2)
points(res_naive[bord,1], -seq_len(p), 
   col = spec_pal[bord], cex = ifelse(res_naive[bord, 4], 1.5, 1), 
   pch = ifelse(res_naive[bord, 4], 15, 16))
axis(1)
axis(2, at = -seq_len(p), labels = spec_names[bord], 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)
abline(v = 0)

dev.print(png, filename = "Results/3c_forestplot_naive.png", 
  units = "in", res = 100)

#-------------------------------------
#   Meta-regression with logratios
#-------------------------------------

# Loop on the constituents to have the isolated contribution of each
# See Hron et al. (2012)
res_comp <- matrix(NA, p, 4, 
  dimnames = list(spec_names, c("coef", "lo", "up", "sig")))
for (j in seq_len(p)){
  design_mat <- ilr(mean_comp[,c(j, setdiff(1:p, j))])
  metamod <- mixmeta(coefall ~ design_mat, vcovall, 
    random = ~ 1|country/city,
    data = cities, method = "reml", subset = conv)
  print(metamod)
  
  # We store only the first coef
  res_comp[j, 1] <- coef(metamod)[2]
  se_spec <- sqrt(diag(vcov(metamod))[2])
  res_comp[j, 2] <- res_comp[j, 1] - 1.96 * se_spec
  res_comp[j, 3] <- res_comp[j, 1] + 1.96 * se_spec
  res_comp[j, 4] <- res_comp[j, 2] > 0 | res_comp[j, 3] < 0
}

# Forestplot
bord <- order(res_comp[,1], decreasing = T)
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(res_comp[bord,1], -seq_len(p), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(res_comp[,2:3]))
segments(res_comp[bord,2], -seq_len(p), res_comp[bord, 3], -seq_len(p), 
  col = "darkgrey", lwd = 2)
points(res_comp[bord,1], -seq_len(p), 
   col = spec_pal[bord], cex = ifelse(res_comp[bord, 4], 1.5, 1), 
   pch = ifelse(res_comp[bord, 4], 15, 16))
axis(1)
axis(2, at = -seq_len(p), labels = spec_names[bord], 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)
abline(v = 0)

dev.print(png, filename = "Results/3c_forestplot_logratio.png", 
  units = "in", res = 100)
  
#-------------------------------------
#  Secondary inorganic aerosol vs other constituents
#-------------------------------------
inor_inds <- sapply(c("SO4", "NH4", "NIT"), grep, colnames(dlist_spec[[1]]))

# Summing secondary inorganic aerosols
mean_inor <- sapply(dlist_spec, function(x){
  xs <- cbind(rowSums(x[,inor_inds]), x[,setdiff(spec_inds, inor_inds)])
  # Zero value imputation as in Martín-Fernández et al. (2003)
  imp <- multRepl(xs, label = 0, dl = rep(1e-5, 5))
  # Transformation to compositional object
  xc <- acomp(imp)
  # Compositional mean (pass through the irl transformation)
  mean(xc)
})
mean_inor <- t(mean_inor)
colnames(mean_inor)[1] <- "INOR"

# Effect modification estimation
inormod <- mixmeta(coefall ~ ilr(mean_inor)[,1], vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)
coef_inor <- coef(inormod)[2]
se_inor <- sqrt(diag(vcov(inormod))[2])
lo_inor <- coef_inor - 1.96 * se_inor
up_inor <- coef_inor + 1.96 * se_inor
sig_inor <- lo_inor > 0 | up_inor < 0

#-------------------------------------
#  Secondary inorganic aerosol as a subcomposition
#    /!\ perhaps too-much zeros
#-------------------------------------
p_inorsub <- length(inor_inds)

mean_inorsub <- sapply(dlist_spec, function(x){
  xs <- x[,inor_inds]
  # Zero value imputation as in Martín-Fernández et al. (2003)
  xs <- xs[apply(xs, 1, function(l) any(l != 0)),]
  if (any(xs == 0)){
    xs <- multRepl(xs, label = 0, dl = rep(1e-5, p_inorsub))
  }
  # Transformation to compositional object
  xc <- acomp(xs)
  # Compositional mean (pass through the irl transformation)
  mean(xc)
})
mean_inorsub <- t(mean_inorsub)

# Effect modifications of each secondary inorganic aerosol vs others
res_inorsub <- matrix(NA, p_inorsub, 4, 
  dimnames = list(sprintf("%svsINOR",c("SO4", "NH4", "NIT")), 
    c("coef", "lo", "up", "sig")))
for (j in seq_len(p_inorsub)){
  design_mat <- ilr(mean_inorsub[,c(j, setdiff(1:p_inorsub, j))])
  metamod <- mixmeta(coefall ~ design_mat, vcovall, 
    random = ~ 1|country/city,
    data = cities, method = "reml", subset = conv)
  print(metamod)
  
  # We store only the first coef
  res_inorsub[j, 1] <- coef(metamod)[2]
  se_spec <- sqrt(diag(vcov(metamod))[2])
  res_inorsub[j, 2] <- res_inorsub[j, 1] - 1.96 * se_spec
  res_inorsub[j, 3] <- res_inorsub[j, 1] + 1.96 * se_spec
  res_inorsub[j, 4] <- res_inorsub[j, 2] > 0 | res_inorsub[j, 3] < 0
}

# Addition of the results for INOR vs other constituents
res_inorsub <- rbind(res_inorsub, 
  INORvsALL = c(coef_inor, lo_inor, up_inor, sig_inor))
  
# Forestplot to summarize these results
p_inor <- nrow(res_inorsub)
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(res_inorsub[,1], -seq_len(p_inor), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(res_inorsub[,2:3]))
segments(res_inorsub[,2], -seq_len(p_inor), res_inorsub[, 3], -seq_len(p_inor), 
  col = "darkgrey", lwd = 2)
points(res_inorsub[,1], -seq_len(p_inor), 
   col = c(spec_pal[1:3], 1), cex = ifelse(res_inorsub[, 4], 1.5, 1), 
   pch = ifelse(res_inorsub[, 4], 15, 16))
axis(1)
axis(2, at = -seq_len(p_inor), labels = rownames(res_inorsub), 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)
abline(v = 0)

dev.print(png, filename = "Results/3c_forestplot_secondaryInorganic.png", 
  units = "in", res = 100)
  