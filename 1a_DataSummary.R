#####################################################################
#
#                 MCC HetPoll
#             Part 1: Summary statistics
#
#####################################################################

library(fields)
library(maps)
library(mapdata)

load("Data/0_Data.RData")

#------------------------------------------
#  Summary tables
#------------------------------------------

summary_tab <- data.frame(country = unique(cities$countryname))

# Number of cities per country
summary_tab$ncities <- table(droplevels(cities$country))
sum(summary_tab$ncities)

# Data period per country
cityperiod <- lapply(dlist, function(x) range(x$year))
countryperiod <- tapply(cityperiod, droplevels(cities$country), 
  function(x) range(unlist(x)))
mean(sapply(countryperiod, diff) + 1) # Average length
summary_tab$period <- sapply(countryperiod, paste, collapse = " - ")

# Total deaths per country
citydeath <- sapply(dlist, function(x) sum(x$death, na.rm = T))
summary_tab$death <- tapply(citydeath, droplevels(cities$country), sum)
sum(summary_tab$death)

# Average PM2.5 concentration per country
citypm <- lapply(dlist, "[", , "tpm25")
summary_tab$pm <- tapply(citypm, droplevels(cities$country), function(x){
  allpm <- na.omit(unlist(x))
  sprintf("%2.1f (%2.1f - %2.1f)", mean(allpm), min(allpm), max(allpm))
})

#--- Export table
write.table(summary_tab, file = "Paper_Figures/Table1.csv", quote = F,
  row.names = F, sep = ";")

#------------------------------------------
#  Map of average PM2.5 concentration
#------------------------------------------

# Average PM2.5 per city
mean_pm <- sapply(citypm, mean, na.rm = T)

# Create colorscale
cutoff <- c(0, 5, 10, 20, 40, 60, 80)
labels <- paste0(paste0(cutoff[-length(cutoff)], "-", cutoff[-1]))
citycat <- cut(mean_pm, cutoff, labels = labels, include.lowest = T)
pal <- tim.colors(length(labels))

# Draw map
x11(width = 10)
map("worldHires", mar=c(0,0,0,0), col = grey(0.95),
  myborder = 0, fill = T, border = grey(0.5), lwd = 0.3)
symbols(cities$long, cities$lat, circles = rep(1.5, nrow(cities)), 
  inches = F, bg = pal[citycat], add = T)
map.scale(-160, -55, ratio = F, cex = 0.8, relwidth = 0.1)
legend(-175, 10, labels, pt.cex = 1.2,
  pch = 21, pt.bg = pal, bty = "n", cex = 0.7, inset = 0.02,
  title = expression(paste("Mean annual ", PM[2.5], "\n concentration (", mu, "g/", m^3, ")")))

dev.print(png, filename = "Results/1a_AvePM_map.png", units = "in", res = 200)
