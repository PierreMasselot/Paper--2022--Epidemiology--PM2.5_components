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

# Initialize summary table with number of cities
summary_tab <- aggregate(city ~ countryname, cities, length)

#----- Add summaries from dlist
dlist_sum <- tapply(dlist, droplevels(cities$countryname), function(d){
  # Data period per country
  period <- range(unlist(lapply(d, "[[", "year")))
  periodstr <- paste(period, collapse = " - ")
  
  # Total deaths per country
  totdeaths <- sum(unlist(lapply(d, "[[", "death")), na.rm = T)
  
  # Average PM2.5 concentration per country
  allpm <- na.omit(unlist(lapply(d, "[[", "pm25")))
  pmstr <- sprintf("%2.1f (%2.1f - %2.1f)", mean(allpm), quantile(allpm, .1), 
    quantile(allpm, .9))
  
  # Return
  c(period = periodstr, death = totdeaths, pm = pmstr)
})

# Add to summary table
summary_tab <- cbind(summary_tab, 
  do.call(rbind, dlist_sum)[match(summary_tab$countryname, names(dlist_sum)),])

#--- Export Table 1
write.table(summary_tab, file = "Results/Table1.csv", quote = F,
  row.names = F, sep = ",")

