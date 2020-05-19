#####################################################################
#
#                 MCC HetPoll
#             Part 2: Second stage meta analysis
#                 No meta-predictor
#
#####################################################################

library(mixmeta)

setwd("C:/Users/masselpl/Documents/Recherche/2020 - LSHTM/Projects/MCC-HetPoll")
load("Data/2_FirstStageResults.RData")

#-------------------------------------
#   Overall association
#-------------------------------------
overall <- mixmeta(coefall ~ 1, S = vcovall, random = ~ 1|city,
  data = cities, method = "reml", subset = conv)
  
cityblup <- blup(overall, vcov = T)

#-------------------------------------
#   Overall with 2-level grouping
#-------------------------------------
overall_country <- mixmeta(coefall ~ 1, vcovall, random = ~ 1|country/city,
  data = cities, method = "reml", subset = conv)
cityblup_country <- blup(overall_country, vcov = T, level = 1, pi = T)
rownames(cityblup_country) <- cities$countryname[conv]
country_blup <- exp(unique(cityblup_country))


# Forest plot ordered from highest blup
country_blup <- country_blup[order(country_blup[,1], decreasing = TRUE),]

par(mar = c(5, 10, 4, 2) + .1)  
plot(country_blup[,1], -(1:nrow(country_blup)), xlab = "RR", ylab = "",
  axes = F, xlim = range(country_blup[,1:3]))
segments(country_blup[,2], -(1:nrow(country_blup)), country_blup[,3],
  -(1:nrow(country_blup)), col = "darkgrey")
points(country_blup[,1], -(1:nrow(country_blup)), pch = 16,
   col = ifelse(country_blup[,2] > 1, "red", "black"))
axis(1)
axis(2, at = -(1:nrow(country_blup)), labels = rownames(country_blup), 
  las = 1, hadj = 1, lwd.ticks = 0, lwd = 0)
abline(v = 1, col = "darkgrey")
abline(v = exp(coef(overall_country)), lwd = 2, lty = 2)
mtext(text = "Overall", at = exp(coef(overall_country)))
dev.print(png, filename = "Results/3a_forestplot_country.png", 
  units = "in", res = 100)