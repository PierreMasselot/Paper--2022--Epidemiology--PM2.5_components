#####################################################################
#
#                 MCC HetPoll
#             Part 1: Summary statistics fo compositional data
#
#####################################################################

library(viridis)
library(compositions)
library(zCompositions)
library(corrplot)

load("Data/0_Data.RData")

countries_pal <- viridis(nrow(countries))
countries_pch <- rep_len(15:18, nrow(countries))

#------------------------------------------
#   SPEC data as compositional object
#------------------------------------------

# Bind all SPEC data in a single data.frame
tot_spec <- do.call(rbind, dlist_spec)

#Useful objects for later
spec_inds <- grep("PM25", colnames(tot_spec))
spec_names <- sapply(strsplit(colnames(tot_spec)[spec_inds], "_"),
  "[", 2)
spec_pal <- c(2, 6, 4, 1, 3, 5, 7)

# Zeros imputation (/!\ think about it more thoroughly)
imp_spec <- multRepl(tot_spec[,spec_inds], label = 0, dl = rep(1e-5, 7))

# Create compositional data object
comp_spec <- acomp(imp_spec)
colnames(comp_spec) <- spec_names

#------------------------------------------
#   Mean composition per year and country
#------------------------------------------

count_split <- split(as.data.frame(imp_spec), 
  list(country = tot_spec$country, year = tot_spec$year))
agg_spec <- t(sapply(count_split, function(x) mean(acomp(x))))
agg_spec_rows <- Reduce(rbind,strsplit(rownames(agg_spec), "\\."))

x11(width = 15, height = 10)
par(mfrow = n2mfrow(nrow(countries) + 1, asp = 1.5), 
  mar = c(3, 2, 3, 1))
for (i in seq_len(nrow(countries))){
  country_inds <- agg_spec_rows[,1] == countries$country[i]
  country_dat <- agg_spec[country_inds,]
  country_dat <- merge(cbind(as.numeric(agg_spec_rows[country_inds, 2]), country_dat), 
    2003:2017, by = 1, all.y = T)
  
  bp <- barplot(t(data.matrix(country_dat[,-1])), col = spec_pal, 
    border = NA, names.arg = country_dat[,1], 
    main = countries$countryname[i])
}
plot.new()
legend("topleft", spec_names, fill = spec_pal, bty = "n", ncol = 2,
  cex = .8, xpd = T, border = NA)

dev.print(png, filename = "Results/1b_CompCountries.png", 
  units = "in", res = 100)
  
#------------------------------------------
#       Total variation matrix
#------------------------------------------

# On city average
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
var_mat <- variation(acomp(mean_comp)) 

x11()
corrplot.mixed(var_mat, is.corr = F, tl.col = spec_pal, tl.cex = 1.5)

dev.print(png, filename = "Results/1b_VariationMatrix.png", 
  units = "in", res = 100)
  
#------------------------------------------
#       PCA
#------------------------------------------

# PCA on city average
pcares <- princomp(acomp(mean_comp))


# Put everything on the same scale
scores <- scale(pcares$scores[,1:2])
loads <- scale(pcares$loadings[,1:2])

x11()
par(mar = c(5, 4, 4, 7) + .1)
plot(scores, xlim = range(c(scores[,1], loads[,1] * 1.5)),
  ylim = range(c(scores[,2], loads[,2] * 1.2)),
  col = countries_pal[as.numeric(droplevels(cities$country))],
  pch = countries_pch[as.numeric(droplevels(cities$country))]
)
arrows(0, 0, loads[,1], loads[,2], lwd = 2, length = .1)
text(loads * 1.2, labels = spec_names, cex = 1.5)
legend(par("usr")[2], par("usr")[4], legend = countries$countryname, 
  col = countries_pal, pch = countries_pch, bty = "n", xpd = T)

dev.print(png, filename = "Results/1b_PCAbiplot.png", 
  units = "in", res = 100)
  

