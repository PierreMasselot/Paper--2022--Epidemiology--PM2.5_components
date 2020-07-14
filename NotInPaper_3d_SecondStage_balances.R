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
library(robCompositions)

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
colnames(mean_comp) <- spec_names

#-------------------------------------
#   Scores of MCC indicators
#-------------------------------------

### Which indicators
# Density has too much NAs
# For the UCD indicators, we take both the 2000 and 2014/2015 ones
indic_names <- c("oldpopprop", "GDP", "Poverty", "avgtmean", 
  "totalrange", "E_GR_AV00", "E_GR_AV14", "B00", "B15")
indic_form <- sprintf("~ %s", paste(indic_names, collapse = " + "))

# PCA
pca_indic <- princomp(as.formula(indic_form), data = mcc.indicators,
  na.action = na.exclude)
# screeplot(pca_indic)
# The first component accounts for 99% of the variance
indicator <- pca_indic$scores[,1]

#-------------------------------------
#   Naive meta-regression
#-------------------------------------

# We do it one at a time since the design matrix is singular (sum to one)
#res_naive <- matrix(NA, p, 4, 
#  dimnames = list(spec_names, c("coef", "lo", "up", "sig")))
#for (j in seq_len(p)){
#  metamod <- mixmeta(coefall ~ mean_comp[,j], vcovall, 
#    random = ~ 1|country/city,
#    data = cities, method = "reml", subset = conv)
#  print(metamod)
#  
#  # We store only the first coef
#  res_naive[j, 1] <- coef(metamod)[2]
#  se_naive <- sqrt(diag(vcov(metamod))[2])
#  res_naive[j, 2] <- res_naive[j, 1] - 1.96 * se_naive
#  res_naive[j, 3] <- res_naive[j, 1] + 1.96 * se_naive
#  res_naive[j, 4] <- res_naive[j, 2] > 0 | res_naive[j, 3] < 0
#}
#
## Forestplot
#bord <- order(res_naive[,1], decreasing = T)
#x11()
#par(mar = c(5, 10, 4, 2) + .1)  
#plot(res_naive[bord,1], -seq_len(p), 
#  xlab = "Meta-regression coefficient", ylab = "",
#  axes = F, xlim = range(res_naive[,2:3]))
#segments(res_naive[bord,2], -seq_len(p), res_naive[bord, 3], -seq_len(p), 
#  col = "darkgrey", lwd = 2)
#points(res_naive[bord,1], -seq_len(p), 
#   col = spec_pal[bord], cex = ifelse(res_naive[bord, 4], 1.5, 1), 
#   pch = ifelse(res_naive[bord, 4], 15, 16))
#axis(1)
#axis(2, at = -seq_len(p), labels = spec_names[bord], 
#  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)
#abline(v = 0)
#
#dev.print(png, filename = "Results/3d_forestplot_naive.png", 
#  units = "in", res = 100)

#-------------------------------------
#   Meta-regression with logratios
#-------------------------------------

# Loop on the constituents to have the isolated contribution of each
# See Hron et al. (2012)
res_comp <- matrix(NA, p, 4, 
  dimnames = list(spec_names, c("coef", "lo", "up", "sig")))
for (j in seq_len(p)){
  design_mat <- data.matrix(pivotCoord(mean_comp, pivotvar = j))
  metamod <- mixmeta(coefall ~ design_mat + indicator, vcovall, 
    random = ~ 1|country/city,
    data = cities, method = "reml", subset = conv)
  
  # We store only the first coef
  res_comp[j, 1] <- coef(metamod)[2]
  se_spec <- sqrt(diag(vcov(metamod))[2])
  res_comp[j, 2] <- res_comp[j, 1] - 1.96 * se_spec
  res_comp[j, 3] <- res_comp[j, 1] + 1.96 * se_spec
  res_comp[j, 4] <- res_comp[j, 2] > 0 | res_comp[j, 3] < 0
}

# Forestplot
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(res_comp[,1], -seq_len(p), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(res_comp[,2:3]))
abline(v = 0)
segments(res_comp[,2], -seq_len(p), res_comp[, 3], -seq_len(p), 
  col = "darkgrey", lwd = 2)
points(res_comp[,1], -seq_len(p), 
   col = spec_pal, cex = ifelse(res_comp[, 4], 1.5, 1), 
   pch = ifelse(res_comp[, 4], 15, 16))
axis(1)
axis(2, at = -seq_len(p), labels = spec_names, 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)

dev.print(png, filename = "Results/3d_forestplot_all.png", 
  units = "in", res = 100)
  

#-------------------------------------
#  "Aggregated" composition
#-------------------------------------
# As in Hvidtfeldt et al. (2019)
# Secondary inorganic pollutants, Organic components and Sea Salt

agg_inds <- list(INOR = 1:3, CARBON = 4:5, SS = 6, DUST = 7)
pa <- length(agg_inds)
agg_pal <- c("purple", "darkgreen", 5, 7)

# Create the aggregation
agg_comp <- t(sapply(dlist_spec, function(x){
    x <- x[,spec_inds]
    # Aggregating
    xs <- sapply(agg_inds, function(inds) rowSums(x[,inds, drop = F]))
    # Zero value imputation as in Martín-Fernández et al. (2003)
    if (any(xs == 0)) xs <- multRepl(xs, label = 0, dl = rep(1e-5, pa))
    # Compositional mean (pass through the irl transformation)
    mean(acomp(xs))
}))

# Loop on groups
res_agg <- matrix(NA, pa, 4, 
  dimnames = list(names(agg_inds), c("coef", "lo", "up", "sig")))
for (j in seq_len(pa)){
  design_mat <- data.matrix(pivotCoord(agg_comp, pivotvar = j))  
  metamod <- mixmeta(coefall ~ design_mat + indicator, vcovall, 
    random = ~ 1|country/city,
    data = cities, method = "reml", subset = conv)
  
  # We store only the first coef
  res_agg[j, 1] <- coef(metamod)[2]
  se_spec <- sqrt(diag(vcov(metamod))[2])
  res_agg[j, 2] <- res_agg[j, 1] - 1.96 * se_spec
  res_agg[j, 3] <- res_agg[j, 1] + 1.96 * se_spec
  res_agg[j, 4] <- res_agg[j, 2] > 0 | res_agg[j, 3] < 0
}

# Forest plot
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(res_agg[,1], -seq_len(pa), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(res_agg[,2:3]))
abline(v = 0)
segments(res_agg[,2], -seq_len(pa), res_agg[, 3], -seq_len(pa), 
  col = "darkgrey", lwd = 2)
points(res_agg[,1], -seq_len(pa), cex = ifelse(res_agg[, 4], 1.5, 1), 
   pch = ifelse(res_agg[, 4], 15, 16), col = agg_pal)
axis(1)
axis(2, at = -seq_len(pa), labels = names(agg_inds), 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)

dev.print(png, filename = "Results/3d_forestplot_aggregated.png", 
  units = "in", res = 100)

#-------------------------------------
#  Balances
#-------------------------------------

# Loop on groups
res_bal <- matrix(NA, pa, 4, 
  dimnames = list(names(agg_inds), c("coef", "lo", "up", "sig")))
for (j in seq_len(pa)){
  pj <- length(agg_inds[[j]])
  bal_mean <- sapply(dlist_spec, function(x){
    imp <- multRepl(x[,spec_inds], label = 0, dl = rep(1e-5, 7))
    # balance
    num <- apply(imp[,agg_inds[[j]], drop = F], 1, function(z) exp(mean(log(z))))
    denom <- apply(imp[,-agg_inds[[j]], drop = F], 1, function(z) exp(mean(log(z))))
    xb <- log(num / denom) * sqrt((pj * (p - pj)) / p)
    
    # Compositional mean (pass through the irl transformation)
    mean(xb)
  })
  
  metamod <- mixmeta(coefall ~ bal_mean + indicator, vcovall, 
    random = ~ 1|country/city,
    data = cities, method = "reml", subset = conv)
  
  # We store only the first coef
  res_bal[j, 1] <- coef(metamod)[2]
  se_spec <- sqrt(diag(vcov(metamod))[2])
  res_bal[j, 2] <- res_bal[j, 1] - 1.96 * se_spec
  res_bal[j, 3] <- res_bal[j, 1] + 1.96 * se_spec
  res_bal[j, 4] <- res_bal[j, 2] > 0 | res_bal[j, 3] < 0
}

# Forest plot
x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(res_bal[,1], -seq_len(pa), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(res_bal[,2:3]))
abline(v = 0)
segments(res_bal[,2], -seq_len(pa), res_bal[, 3], -seq_len(pa), 
  col = "darkgrey", lwd = 2)
points(res_bal[,1], -seq_len(pa), cex = ifelse(res_bal[, 4], 1.5, 1), 
   pch = ifelse(res_bal[, 4], 15, 16), col = agg_pal)
axis(1)
axis(2, at = -seq_len(pa), labels = names(agg_inds), 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)

dev.print(png, filename = "Results/3d_forestplot_balances.png", 
  units = "in", res = 100)


#-------------------------------------
#  Fossil-fuel/traffic related subcomposition
#------------------------------------- 
# Sulfate originate from sulfur that mainly comes from fossil fuel combustion (power plants)
# Nitrate is secondary from NOx emission mainly traffic related and gasoline power
# BC is originate from diesel vehicles and traffic more generally

traff_inds <- c(1, 3, 4)
ptr <- length(traff_inds)

traff_bal <- sapply(dlist_spec, function(x){
  imp <- multRepl(x[,spec_inds], label = 0, dl = rep(1e-5, 7))
  # balance
  num <- apply(imp[,traff_inds, drop = F], 1, function(z) exp(mean(log(z))))
  denom <- apply(imp[,-traff_inds, drop = F], 1, function(z) exp(mean(log(z))))
  xb <- log(num / denom) * sqrt((ptr * (p - ptr)) / p)
  
  # Compositional mean (pass through the irl transformation)
  mean(xb)
})

# Regression
traffic_reg <- mixmeta(coefall ~ traff_bal + indicator, vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)
summary(traffic_reg) # Not significant