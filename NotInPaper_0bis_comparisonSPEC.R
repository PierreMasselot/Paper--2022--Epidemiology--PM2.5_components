#####################################################################
#
#                 MCC HetPoll
#             Comparison of SPEC datasets
#
#####################################################################

library(colorspace)

setwd("C:/Users/masselpl/Documents/Recherche/2020 - LSHTM/Projects/MCC-HetPoll")

#------------------------------------------
# Load datasets
#------------------------------------------

datapath <- "C:/Users/masselpl/Documents/Recherche/2020 - LSHTM/Data/MCC"

# File list
fexp <- "PM2\\.5_SPEC"
flist <- list.files(sprintf("%s/MCC_PM_SPEC_New", datapath), 
  pattern = fexp)

# Prepare data list
dlist_buffer <- dlist_point <-  vector("list", length(flist) / 2)
names(dlist_buffer) <- names(dlist_point) <- as.character(2003:2017)

# Load buffer data
ind_buffer <- grep("buffer", flist)
for (i in seq_along(ind_buffer)){
  dlist_buffer[[i]] <- read.csv(sprintf("%s/MCC_PM_SPEC_New/%s", 
    datapath, flist[ind_buffer[i]]))
  dlist_buffer[[i]]$year <- names(dlist_buffer)[i]  
}
ddf_buffer <- do.call(rbind, dlist_buffer)

# Load point data
ind_point <- setdiff(seq_along(flist),ind_buffer)
for (i in seq_along(ind_point)){
  dlist_point[[i]] <- read.csv(sprintf("%s/MCC_PM_SPEC_New/%s", 
    datapath, flist[ind_point[i]]))
  dlist_point[[i]]$year <- names(dlist_point)[i]  
}
ddf_point <- do.call(rbind, dlist_point)

#------------------------------------------
# Dataset comparison
#------------------------------------------

spec_col <- grep("PM25", colnames(ddf_point))
spec_names <- substring(colnames(ddf_point)[spec_col], 6)

spec_pal <- 1:length(spec_col)
country_pal <- rainbow_hcl(length(unique(ddf_point$countryname)))

# Number of zeros for each costituent
spec_zero_point <- apply(ddf_point[,spec_col], 2, function(x){
  mean(x == 0)
})
spec_zero_buffer <- apply(ddf_buffer[,spec_col], 2, function(x){
  mean(x == 0)
})

x11()
barplot(rbind(spec_zero_point, spec_zero_buffer), 
  col = rep(spec_pal, each = 2), 
  names.arg = spec_names, ylab = "# Zero values",
  beside = T, density = rep(c(-1,20), length(spec_col)))
legend("topleft", c("Point", "Buffer"), density = c(-1, 20),
  bty = "n")
dev.print(png, filename = "Results/0bis_NbZero.png", 
  units = "in", res = 100)
  
# Mean annual concentration
spec_mean_point <- apply(ddf_point[,spec_col], 2, function(x){
  mean(x[x > 0])})
spec_mean_buffer <- apply(ddf_buffer[,spec_col], 2, function(x){
  mean(x[x > 0])})

x11()
barplot(spec_mean_point - spec_mean_buffer, col = spec_pal, 
  names.arg = spec_names, ylab = "Mean (point - buffer)")
dev.print(png, filename = "Results/0bis_MeanConc.png", 
  units = "in", res = 100)
  
# Number of zeros per country
count_zero_point <- aggregate(ddf_point[,spec_col], 
  by = list(country = ddf_point$countryname), 
  function(x) sum(x == 0))
count_zero_buffer <- aggregate(ddf_buffer[,spec_col], 
  by = list(country = ddf_buffer$countryname), 
  function(x) sum(x == 0))

x11(width = 10)
matplot(count_zero_point[,-1] - count_zero_buffer[,-1], col = spec_pal,
  ylab = "# Zero values (point - buffer)", pch = 16, xaxt = "n")
axis(1, at = 1:nrow(count_zero_point), 
  labels = unique(ddf_point$countryname), las = 3)
legend("topleft", spec_names, col = spec_pal, pch = 16,
  bty = "n", ncol = 2)
dev.print(png, filename = "Results/0bis_NbZeroCountry.png", 
  units = "in", res = 100)

# Mean per country
count_mean_point <- aggregate(ddf_point[,spec_col], 
  by = list(country = ddf_point$countryname), 
  function(x) mean(x[x > 0]))
count_mean_buffer <- aggregate(ddf_buffer[,spec_col], 
  by = list(country = ddf_buffer$countryname), 
  function(x) mean(x[x > 0]))

x11(width = 10)
matplot(count_mean_point[,-1] - count_mean_buffer[,-1], col = spec_pal,
  ylab = "Mean concentration (point - buffer)", pch = 16, xaxt = "n")
axis(1, at = 1:nrow(count_mean_point), 
  labels = unique(ddf_point$countryname), las = 3)
legend("bottomleft", spec_names, col = spec_pal, pch = 16,
  bty = "n", ncol = 2)
dev.print(png, filename = "Results/0bis_MeanCountry.png", 
  units = "in", res = 100)