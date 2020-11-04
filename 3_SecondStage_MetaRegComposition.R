#####################################################################
#
#                 MCC HetPoll
#             Part 2: Second stage meta analysis
#              experiment with mixtures
#
#####################################################################

# Here I simply consider the approach of Aitchison & Bacon-Shone (1984)
# It consists in transforming the constituents with ALR and using them in a regression model.

library(mixmeta)
library(compositions)
library(zCompositions)
library(fields)
library(maps)
library(mapdata)


load("Data/2_FirstStageResults.RData")

#-------------------------------------
#   Preparing constituent objects
#-------------------------------------

spec_inds <- grep("PM25", colnames(dlist_spec[[1]]))
spec_names <- c("SO4", "NH4", "NO3", "BC", "OC", "SS", "DUST")
spec_labs <- c(expression(SO[4]^{"2-"}), expression(NH[4]^{"+"}), expression(NO[3]^{"-"}), 
  "BC", "OC", "SS", "DUST")
spec_pal <- c(2, 6, 4, 1, 3, 5, 7)

# Mean per city
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

# Transforming as ALR
alr_comp <- alr(mean_comp)

#-------------------------------------
#   PC Scores of MCC indicators
#-------------------------------------

indic_form <- sprintf("~ %s", paste(indic_names, collapse = " + "))

# PCA
pca_indic <- prcomp(as.formula(indic_form), data = mcc.indicators,
  na.action = na.exclude, scale. = T)
# screeplot(pca_indic)
# The first two components account for 58% of the variance
indicators <- pca_indic$x[,1:2]

#-------------------------------------
#  Main meta-regression model
#-------------------------------------

# Model with ALR as meta-predictors
metamod <- mixmeta(coefall ~ alr_comp + indicators, vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)
  
summary(metamod) # To obtain Q and I2 statistics

#-------------------------------------
# BLUPS
#-------------------------------------

# Compute RRs
BLUPall <- exp(blup(metamod))

# Extremes
cities[which.max(BLUPall),]
cities[which.min(BLUPall),]

# Number of RR > 1
sum(BLUPall > 1)

#-------------------------------------
# Figure 1: Map with BLUPS
#-------------------------------------

# Compute average PM2.5
mean_pm <- sapply(dlist, function(d) mean(d$pm25, na.rm = T))

# Create colorscale based on RR
cutoff <- seq(0.975, 1.025, by = 0.005)
labels <- paste0(paste0(cutoff[-length(cutoff)], "-", cutoff[-1]))
citycat <- cut(BLUPall, cutoff, labels = labels, include.lowest = T)
pal <- tim.colors(length(labels))

# Create point size based on mean PM2.5
ptsiz_rng <- c(.5, 2.5)
mean_scale <- (mean_pm - min(mean_pm)) / (diff(range(mean_pm)))
pt_size <- mean_scale * diff(ptsiz_rng) + ptsiz_rng[1]
pm_scale <- pretty(mean_pm)
size_scale <- (pm_scale - min(mean_pm)) / (diff(range(mean_pm))) * 
  diff(ptsiz_rng) + ptsiz_rng[1]

#---- Draw map
x11(width = 10)
map("worldHires", mar=c(0,0,0,0), col = grey(0.95),
    myborder = 0, fill = T, border = grey(0.5), lwd = 0.3)
# Add points with size = RR and colorscale for mean PM2.5
points(cities$long, cities$lat, pch = 21, cex = pt_size, bg = pal[citycat])
# Scale and legend
map.scale(-5, -50, ratio = F, cex = 0.7, relwidth = 0.1)
rect(par("usr")[1], 45, -129, -70, border = NA, col = "white")
lg <- legend(-175, 47, labels, pt.cex = 1.2, bg = "white",
             pch = 21, pt.bg = pal, box.col = "white", cex = 0.7, inset = 0.02,
             title = "Predicted RR", title.adj = 0
)
legend(lg$rect$left, lg$rect$top - lg$rect$h, pm_scale, inset = 0.02, 
       pch = 21, pt.cex = size_scale, pt.bg = "grey",box.col = "white", cex = 0.7,   
       #  title = expression(paste("Mean annual ", PM[2.5], "(", mu, "g/", m^3, ")")),
       title = "", bg = "white",
       title.adj = 1, y.intersp = 1.5, xjust = 0
)
text(lg$rect$left, lg$rect$top - lg$rect$h, cex = .7, adj = c(0, 1.3),
     expression(paste("Mean ", PM[2.5], " concentration (", mu, "g/", m^3, ")")))

dev.print(png, filename = "Results/Figure1.png", units = "in", res = 200)
dev.print(pdf, file = "Results/Figure1.pdf")

#-------------------------------------
#  Prediction of the RR for each components separately
#-------------------------------------

# We create a sequence for the component
cseq <- seq(.01, .99, by = .01)

# Overall mean of composition
ov_mean <- mean(acomp(mean_comp))

# Prepare prediction data.frame
newdat <- list(indicators = matrix(0, length(cseq), 2, 
  dimnames = list(NULL, sprintf("PC%i", 1:2))))

# Prepare objects to store predictions, confidence limits and whether the proportion value is observed for the component
preds <- plo <- pup <- is_obs <- matrix(NA, length(cseq), p)
for (j in seq_len(p)){
  # Create a compositional grid. The component of interests varies from 0 to 1 (excluded) and the other are taken as the overall mean adjusted for closure while keeping the subcomposition constant
  xj <- matrix(NA, length(cseq), p, dimnames = list(NULL, spec_names))
  xj[,-j] <- sapply(ov_mean[-j] / sum(ov_mean[-j]), "*", 1 - cseq)
  xj[,j] <- cseq
  
  # Transform by ALR and add to the newdat data.frame
  newdat$alr_comp <- alr(xj)
  
  # Make predictions
  pred_obj <- predict(metamod, newdata = newdat, ci = T)
  
  # Store predictions and prediction intervals only in the range of observed values
  obs_rng <- range(mean_comp[,j])
  is_obs[,j] <- cseq >= obs_rng[1] & cseq <= obs_rng[2]
  
  preds[,j] <- exp(pred_obj[,1])
  plo[,j] <- exp(pred_obj[,2])
  pup[,j] <- exp(pred_obj[,3])
}
preds_obs <- preds
preds_obs[!is_obs] <- NA

#-------------------------------------
#  Figure 4: all predictions on different panels
#-------------------------------------

# Panel matrix
pm <- matrix(rep(c(1:6, 0, 7, 0), each = 2), nrow = 3, byrow = T)

# Plot
x11(height = 10, width = 10)
par(mar = c(5, 4, 3, 1), cex.main = 1.5, cex.lab = 1.2)
layout(pm)
for (j in seq_len(p)){
  # Initiate an empty plot with labels
  plot(0, 0, col = NA, xlim = 100 * range(mean_comp), 
    ylim = c(min(plo, na.rm = T), max(pup, na.rm = T)),
    ylab = "RR", xlab = "Proportion (%)",
    main = bquote(bold(.(spec_labs[j][[1]]))))
  # Draw prediction interval
  polygon(100 * c(cseq[!is_obs[,j]], rev(cseq[!is_obs[,j]])), 
    c(plo[!is_obs[,j],j], rev(pup[!is_obs[,j],j])), 
    col = adjustcolor(spec_pal[j], .3), border = NA, density = 20)
  polygon(100 * c(cseq[is_obs[,j]], rev(cseq[is_obs[,j]])), 
    c(plo[is_obs[,j],j], rev(pup[is_obs[,j],j])), 
    col = adjustcolor(spec_pal[j], .2), border = NA)
  # Add prediction
  lines(100 * cseq, preds[,j], lwd = 1, col = spec_pal[j])
  lines(100 * cseq[is_obs[,j]], preds[is_obs[,j],j], lwd = 3, col = spec_pal[j])
  # Add lines indicating observed range
  abline(v = 100 * range(cseq[is_obs[,j]]), lty = 2)
  abline(h = 1)
}

dev.print(png, filename = "Results/Figure4.png", 
  units = "in", res = 100)
dev.print(pdf, file = "Results/Figure4.pdf")

#-------------------------------------
#  Comparison with nested models
#-------------------------------------

# Apply model with only indicator PC
metaindic <- mixmeta(coefall ~ indicators, vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)

# Apply null model: no meta-predictor
metanull <- mixmeta(coefall, vcovall, random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)

allmodels <- list(full = metamod, indic = metaindic, null = metanull)

#--- create comparison table ---
compar_tab <- data.frame(
  qstat = sapply(allmodels, function(x) summary(x)$qstat$Q),
  i2stat = sapply(allmodels, function(x) summary(x)$i2stat)
)

fwald <- function(full, null) {
  ind <- !names(coef(full)) %in% names(coef(null))
  coef <- coef(full)[ind]
  vcov <- vcov(full)[ind,ind]
  waldstat <- coef %*% solve(vcov) %*% coef
  df <- length(coef)
  pval <- 1 - pchisq(waldstat, df)
  return(list(waldstat = waldstat, pvalue = pval))
}

full_wald <- fwald(metamod, metaindic)
indic_wald <- fwald(metaindic, metanull)

# Add to the table
compar_tab$Wald_stat <- c(full_wald$waldstat, indic_wald$waldstat, NA)
compar_tab$Wald_pvalue <- c(full_wald$pvalue, indic_wald$pvalue, NA)

#---- Export Table 2
write.table(compar_tab, file = "Results/Table2.csv", quote = F,
  row.names = F, sep = ";")
  
#-------------------------------------
#  Save results
#-------------------------------------

save.image("Data/3_MetaResults.RData")
