#------------------------------------------
#       Ternary scatterplots
#------------------------------------------
overall_mean <- mean(acomp(mean_comp))

# SIA
x11()
TernaryPlot(atip = spec_labs[1], btip = spec_labs[2], ctip = spec_labs[3], 
            lab.col = spec_pal[1:3], tip.col = spec_pal[1:3], 
            grid.col = "white", grid.minor.col = "white")
TernaryPoints(mean_comp[,1:3])
TernaryPoints(overall_mean[1:3], pch = 4, lwd = 3, cex = 1.5, col = "red")

dev.print(png, filename = "Results/scatter_SO4_NH4_NO3.png", 
          units = "in", res = 200)
dev.print(pdf, file = "Results/scatter_SO4_NH4_NO3.pdf")

# Carbonaceous and nitrate
x11()
TernaryPlot(atip = spec_labs[3], btip = spec_labs[4], ctip = spec_labs[5], 
            lab.col = spec_pal[3:5], tip.col = spec_pal[3:5], 
            grid.col = "white", grid.minor.col = "white")
TernaryPoints(mean_comp[,3:5])
TernaryPoints(overall_mean[3:5], pch = 4, lwd = 3, cex = 1.5, col = "red")

dev.print(png, filename = "Results/scatter_NO3_BC_OC.png", 
          units = "in", res = 200)
dev.print(pdf, file = "Results/scatter_NO3_BC_OC.pdf")

# Naturals 
x11()
TernaryPlot(atip = spec_labs[5], btip = spec_labs[6], ctip = spec_labs[7], 
            lab.col = spec_pal[5:7], tip.col = spec_pal[5:7], 
            grid.col = "white", grid.minor.col = "white")
TernaryPoints(mean_comp[,5:7])
TernaryPoints(overall_mean[5:7], pch = 4, lwd = 3, cex = 1.5, col = "red")

dev.print(png, filename = "Results/scatter_OC_SS_DUST.png", 
          units = "in", res = 200)
dev.print(pdf, file = "Results/scatter_OC_SS_DUST.pdf")

#------------------------------------------
#  Correlation between prop and total PM mass
#------------------------------------------

total_pm <- sapply(dlist_spec, function(d)sum(d$SPEC_total))

cor(log(mean_comp), total_pm)