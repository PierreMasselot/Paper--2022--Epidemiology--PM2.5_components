#####################################################################
#
#                 MCC HetPoll
#             Part 2: Second stage meta analysis
#              experiment with mixtures
#
#####################################################################

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

# Components and labels
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
pca_indic <- prcomp(as.formula(indic_form), data = cities,
  na.action = na.exclude, scale. = T)
# screeplot(pca_indic)
# The first two components account for 58% of the variance
indicators <- pca_indic$x[,1:2]

#-------------------------------------
#  Main meta-regression model
#-------------------------------------

# Model with ALR as meta-predictors
metamod <- mixmeta(coefall ~ alr_comp + indicators, vcovall, 
  random = ~ 1|country/city, data = cities, method = "reml", subset = conv)

# Q and I2 statistics  
summary(metamod) 

#-------------------------------------
# BLUPS
#-------------------------------------

# Compute RRs
BLUPall <- rep(NA, nrow(coefall))
BLUPall[conv] <- exp(blup(metamod))

# Extremes
cities[which.max(BLUPall),]
cities[which.min(BLUPall),]

# Number of RR > 1
sum(BLUPall > 1, na.rm = T)

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
lg <- legend(-175, 40, labels, pt.cex = 1.2, bg = "white",
  pch = 21, pt.bg = pal, box.col = "white", cex = 0.7, inset = 0.02,
  title = expression(paste("Predicted RR (10 ", mu, "g/", m^3, ")")), 
  title.adj = 0
)
legend(lg$rect$left, lg$rect$top - lg$rect$h, pm_scale, inset = 0.02, 
       pch = 21, pt.cex = size_scale, pt.bg = "grey",box.col = "white", cex = 0.7,   
       title = expression(paste("Mean ", PM[2.5], " concentration (", mu, "g/", m^3, ")")),
       bg = "white",
       title.adj = 0, y.intersp = 1.5, xjust = 0, ncol = length(pm_scale)
)
# text(lg$rect$left, lg$rect$top - lg$rect$h, cex = .7, adj = c(0, 1.3),
#      expression(paste("Mean ", PM[2.5], " concentration (", mu, "g/", m^3, ")")))

dev.print(png, filename = "Results/Figure1.png", units = "in", res = 200)
dev.print(pdf, file = "Results/Figure1.pdf")

#-------------------------------------
#  Get meta-coefficients
#-------------------------------------

# Coefficient retrieving
coef_est <- rep(NA, p)
coef_est[-p] <- coef(metamod)[2:p]
coef_est[p] <- -sum(coef_est, na.rm = T) # beta7 = - sum(beta1, ..., beta6)

# Standard errors
se_est <- rep(NA, p)
se_est[-p] <- sqrt(diag(vcov(metamod))[2:p])
se_est[p] <- sqrt(sum(vcov(metamod)[2:p,2:p])) # var(beta7) = sum_j sum_k cov(betaj, betak) (variance of sum of random variables)

# Confidence intervals
lo_est <- coef_est - 1.96 * se_est
up_est <- coef_est + 1.96 * se_est
sig_est <- lo_est > 0 | up_est < 0

# We mulitply everything by ln(2) to interpret it as the impact of doubling the relative proportion of the components
coef_est <- log(2) * coef_est
lo_est <- log(2) * lo_est
up_est <- log(2) * up_est

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
  
  # Exponentiate predictions to put on RR scale
  preds[,j] <- exp(pred_obj[,1])
  plo[,j] <- exp(pred_obj[,2])
  pup[,j] <- exp(pred_obj[,3])
}

# Copy preds and keep only observed proportions
preds_obs <- preds
preds_obs[!is_obs] <- NA

#-------------------------------------
#  Figure 3: Meta-regression resuls
#-------------------------------------

# Panel matrix
pm <- rbind(c(9, rep(10, 3)), 
  cbind(1, matrix(c(2:7, 0, 8, 0), nrow = 3, byrow = T)))

# Initialize plot
x11(height = 10, width = 15)
par(mar = c(5, 5, 3, 1), cex.main = 1.5, cex.lab = 1.2)
layout(pm, width = c(.4, .2, .2, .2), height = c(.1, .3, .3, .3))

#----- Panel A: meta-coefficients

plot(exp(coef_est), -seq_len(7), 
  xlab = "Relative Excess Risk", ylab = "", cex.lab = 1.3,
  axes = F, xlim = range(c(exp(lo_est), exp(up_est))))
abline(v = 1, lty = 2)
segments(exp(lo_est), -seq_len(7), exp(up_est), -seq_len(7), lwd = 2,
  col = grey(.2))
points(exp(coef_est), -seq_len(7), cex = 3, 
  pch = ifelse(sig_est, 15, 16), col = spec_pal)
axis(1, cex.axis = 1.1)
axis(2, at = -seq_len(7), labels = spec_labs, 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0, cex.axis = 1.5)
box()

#----- Panel B: predicted values
par(mar = c(5, 4, 3, 1))
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

#----- Add letters for plots
par(mar = rep(0, 4), cex.main = 1.5, cex.lab = 1.2)
plot.new()
text(par("usr")[1], par("usr")[3], "A", cex = 3, adj = c(0, 0), xpd = T)

plot.new()
text(par("usr")[1], par("usr")[3], "B", cex = 3, adj = c(0, 0), xpd = T)

#----- Save

dev.print(png, filename = "Results/Figure3.png", 
  units = "in", res = 100)
dev.print(pdf, file = "Results/Figure3.pdf")

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
