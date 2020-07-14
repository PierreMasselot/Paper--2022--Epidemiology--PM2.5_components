#####################################################################
#
#                 MCC HetPoll
#             Part 1: Summary statistics
#
#####################################################################

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
  