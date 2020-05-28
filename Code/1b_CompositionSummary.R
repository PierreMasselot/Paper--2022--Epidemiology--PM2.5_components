#####################################################################
#
#                 MCC HetPoll
#             Part 1: Summary statistics fo compositional data
#
#####################################################################

library(viridis)
library(compositions)
library(zCompositions)

load("Data/0_Data.RData")

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

agg_spec <- aggregate(comp_spec, 
  by = list(country = tot_spec$country, year = tot_spec$year),
  mean.acomp)

x11(width = 15, height = 10)
par(mfrow = n2mfrow(nrow(countries) + 1, asp = 1.5), 
  mar = c(3, 2, 3, 1))
for (i in seq_len(nrow(countries))){
  country_dat <- agg_spec[agg_spec$country == countries$country[i], -1]
  country_dat <- merge(country_dat, 2003:2017, by = 1, all.y = T)
  
  bp <- barplot(t(data.matrix(country_dat[,-1])), col = spec_pal, 
    border = NA, names.arg = country_dat[,1], 
    main = countries$countryname[i])
}
plot.new()
legend("topleft", spec_names, fill = spec_pal, bty = "n", ncol = 2,
  cex = .8, xpd = T, border = NA)

dev.print(png, filename = "Results/1b_CompCountries.png", 
  units = "in", res = 100)