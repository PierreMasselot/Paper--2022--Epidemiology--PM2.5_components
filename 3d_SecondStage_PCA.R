#####################################################################
#
#                 MCC HetPoll
#             Part 2: Second stage meta analysis
#               with compositional PCA
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
#         Apply PCA on components
#------------------------------------- 

pcares <- princomp(acomp(mean_comp))

#-------------------------------------
#     Use PCs as meta-predictors
#------------------------------------- 

pca_reg <- mixmeta(coefall ~ pcares$scores[,1:2] + indicator, vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)
coef_pca <- coef(pca_reg)[2:3]
se_pca <- sqrt(diag(vcov(pca_reg))[2:3])
lo_pca <- coef_pca - 1.96 * se_pca
up_pca <- coef_pca + 1.96 * se_pca
sig_pca <- lo_pca > 0 | up_pca < 0

x11()
par(mar = c(5, 10, 4, 2) + .1)  
plot(coef_pca, -seq_len(2), 
  xlab = "Meta-regression coefficient", ylab = "",
  axes = F, xlim = range(c(lo_pca, up_pca)))
segments(lo_pca, -seq_len(2), up_pca, -seq_len(2), 
  col = "darkgrey", lwd = 2)
points(coef_pca, -seq_len(2), cex = ifelse(sig_pca, 1.5, 1), 
   pch = ifelse(sig_pca, 15, 16))
axis(1)
axis(2, at = -seq_len(2), labels = sprintf("Score %i", 1:2), 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)
abline(v = 0)

dev.print(png, filename = "Results/3d_forestplot_pca.png", 
  units = "in", res = 100)


#-------------------------------------
#   Derive similar composition
#------------------------------------- 
vs_var <- list(c("BC", "OC", "NH4"), c("SS", "DUST"))

scores_comp <- t(sapply(dlist_spec, function(x){
  # Zero value imputation as in Martín-Fernández et al. (2003)
  imp <- multRepl(x[,spec_inds], label = 0, dl = rep(1e-5, 7))
  colnames(imp) <- spec_names
  imp <- cbind(rowSums(imp[,vs_var[[1]]]), rowSums(imp[,vs_var[[2]]])) 
  # Transformation to compositional object
  xc <- acomp(imp)
  # Compositional mean (pass through the irl transformation)
  mean(xc)
}))
colnames(scores_comp) <- c("PC", "NATURAL")

# Apply meta-reg
newpred_reg <- mixmeta(coefall ~ alr(scores_comp) + indicator, vcovall, 
  random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)
summary(newpred_reg) # Still significant this way
