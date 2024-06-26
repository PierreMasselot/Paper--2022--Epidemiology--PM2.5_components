#####################################################################
#
#                 MCC HetPoll
#             Part 1: Summary statistics fo compositional data
#
#####################################################################

library(compositions)
library(zCompositions)
library(corrplot)

load("Data/0_Data.RData")

#------------------------------------------
#   SPEC data as compositional object
#------------------------------------------

# Bind all SPEC data in a single data.frame
tot_spec <- do.call(rbind, dlist_spec)

#Useful objects for later
spec_inds <- grep("PM25", colnames(tot_spec))
spec_names <- c("SO4", "NH4", "NO3", "BC", "OC", "SS", "DUST")
spec_labs <- c(expression(SO[4]^{"2-"}), expression(NH[4]^{"+"}), 
  expression(NO[3]^{"-"}), "BC", "OC", "Sea salt", "Dust")
spec_pal <- c(2, 6, 4, 1, 3, 5, 7)

# Zeros imputation
imp_spec <- multRepl(tot_spec[,spec_inds], label = 0, dl = rep(1e-5, 7))

# Create compositional data object
comp_spec <- acomp(imp_spec)
colnames(comp_spec) <- spec_names

#------------------------------------------
#  Proportion of zeros of each component
#------------------------------------------

prop_zero <- apply(tot_spec[,spec_inds], 2, function(x) mean(x == 0))

#------------------------------------------
#   Mean composition per year and country
#------------------------------------------

# Split by country
count_split <- split(as.data.frame(imp_spec), 
  list(country = tot_spec$country, year = tot_spec$year))

# Average by country and year
agg_spec <- t(sapply(count_split, function(x) mean(acomp(x))))
agg_spec_rows <- Reduce(rbind,strsplit(rownames(agg_spec), "\\."))

#----- Plot the average
x11(width = 15, height = 10)
par(mfrow = n2mfrow(nrow(countries) + 1, asp = 1.5), 
  mar = c(3, 2, 3, 1))
for (i in seq_len(nrow(countries))){
  # Select country
  country_inds <- agg_spec_rows[,1] == countries$country[i]
  country_dat <- agg_spec[country_inds,]
  country_dat <- merge(cbind(as.numeric(agg_spec_rows[country_inds, 2]), country_dat), 
    2003:2017, by = 1, all.y = T)
  
  # Barplot
  bp <- barplot(t(data.matrix(country_dat[,-1])), col = spec_pal, 
    border = NA, names.arg = country_dat[,1], 
    main = countries$countryname[i])
}
plot.new()
legend("topleft", spec_labs, fill = spec_pal, bty = "n", ncol = 3,
  cex = 1.8, xpd = NA, border = NA)

dev.print(pdf, file = "Results/Figure2.pdf")

  
#------------------------------------------
#       Total variation matrix
#------------------------------------------

# On city average
mean_comp <- sapply(dlist_spec, function(x){
  # Zero value imputation as in Mart�n-Fern�ndez et al. (2003)
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
spec_labs[6] <- expression("SS")
colnames(var_mat) <- rownames(var_mat) <- paste(":", spec_labs) 


# Color palette
colpal <- colorRampPalette(c('#FFFFFF', '#D1E5F0', '#92C5DE',
  '#4393C3', '#2166AC', '#053061'))(20)

x11()
corrplot.mixed(var_mat, is.corr = F, tl.col = spec_pal, tl.cex = 1.5,
  lower.col = colpal, upper.col = colpal)

dev.print(png, file = "Results/eFigure7.png", res = 1000, units = "in")
  