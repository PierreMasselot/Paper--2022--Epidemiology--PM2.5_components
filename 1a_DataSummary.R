#####################################################################
#
#                 MCC HetPoll
#             Part 1: Summary statistics
#
#####################################################################

load("Data/0_Data.RData")

#------------------------------------------
#  Country level summary table
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
citypm <- lapply(dlist, "[", , "pm25")
summary_tab$pm <- tapply(citypm, droplevels(cities$country), function(x){
  allpm <- na.omit(unlist(x))
  sprintf("%2.1f (%2.1f - %2.1f)", mean(allpm), quantile(allpm, .1), 
    quantile(allpm, .9))
})

#--- Export Table 1
write.table(summary_tab, file = "Results/Table1.csv", quote = F,
  row.names = F, sep = ";")

#------------------------------------------
#  Mean PM2.5 by city for the map
#------------------------------------------

mean_pm <- sapply(dlist, function(d) mean(d$pm25, na.rm = T))

## Used in the map at the end of script 3_SecondStage_MetaRegComposition.R

save(mean_pm, summary_tab, file = "Data/1_Summary.RData")