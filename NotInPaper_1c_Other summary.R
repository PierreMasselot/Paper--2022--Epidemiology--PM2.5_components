#------------------------------------------
#  Time series length compared with whole data
#------------------------------------------

countries_pal <- viridis(nrow(countries))
countries_pch <- rep_len(15:18, nrow(countries))

x11()
plot(sapply(dlist_orig, nrow), sapply(dlist, nrow), 
  col = countries_pal[as.numeric(droplevels(cities$country))],
  pch = countries_pch[as.numeric(droplevels(cities$country))], 
  xlab = "Time series length", 
  ylab = "Common years with SPEC"
)
abline(a = 0, b = 1)
legend("topleft", legend = countries$countryname, 
  col = countries_pal, pch = 15:18,
  ncol = 2, bty = "n")
dev.print(png, filename = "Results/1a_TSlength.png", 
  units = "in", res = 100)

#------------------------------------------
# Some descriptive statistics of SPEC
#------------------------------------------

tot_spec <- do.call(rbind, dlist_spec)

spec_inds <- grep("PM25", colnames(tot_spec))
spec_names <- sapply(strsplit(colnames(tot_spec)[spec_inds], "_"),
  "[", 2)
spec_pal <- c(2, 6, 4, 1, 3, 5, 7)

# Total variation of constituents
x11()
boxplot(tot_spec[,spec_inds], border = spec_pal, lwd = 2, names = spec_names, 
  col = "white")
dev.print(png, filename = "Results/1a_SpecBoxplot.png", 
  units = "in", res = 100)

# Per country
agg_spec <- aggregate(tot_spec[,spec_inds], 
  by = list(country = tot_spec$country, year = tot_spec$year),
  mean)
agg_pm25 <- aggregate(mon_total ~ country + year, data = tot_spec, mean)

x11(width = 15, height = 10)
par(mfrow = n2mfrow(nrow(countries) + 1, asp = 1.5), 
  mar = c(3, 2, 3, 1))
for (i in seq_len(nrow(countries))){
  country_dat <- agg_spec[agg_spec$country == countries$country[i], -1]
  country_dat <- merge(country_dat, 2003:2017, by = 1, all.y = T)
  country_pm25 <- agg_pm25[agg_pm25$country == countries$country[i], -1]
  country_pm25 <- merge(country_pm25, 2003:2017, by = 1, all.y = T)
  
  bp <- barplot(t(data.matrix(country_dat[,-1])), col = spec_pal, 
    border = NA, names.arg = country_dat[,1], 
    main = countries$countryname[i], 
    ylim = c(0, max(country_pm25[,2], na.rm = T)))
  lines(bp, country_pm25$mon_total, type = "b", lwd = 2, 
    col = "black", xpd = T)
}
plot.new()
legend("topleft", spec_names, fill = spec_pal, bty = "n", ncol = 2,
  cex = .8, xpd = T, border = NA)
dev.print(png, filename = "Results/1a_SpecCountries.png", 
  units = "in", res = 100)

# Comparison between mean observed PM2.5 concentration and species
spec_sum <- rowSums(tot_spec[,spec_inds])

x11()
par(mar = c(5, 4, 4, 7) + .1)
plot(as.numeric(tot_spec$year), (tot_spec$mon_total - spec_sum) / tot_spec$mon_total, 
  col = countries_pal[as.numeric(as.factor(tot_spec$country))],
  pch = countries_pch[as.numeric(as.factor(tot_spec$country))], 
  ylab = "Monitor - SPEC sum (Rel. Difference)", 
  xlab = "Year")
abline(h = 0)
legend(par("usr")[2], par("usr")[4], countries$countryname, 
  pch = countries_pch, col = countries_pal, bty = "n", xpd = T)
dev.print(png, filename = "Results/1a_PM25difference.png", 
  units = "in", res = 100)

# Correlation between each constituent and monitored PM2.5
x11()
par(mfrow = n2mfrow(length(spec_names) + 2))
for (i in seq_along(spec_names)){
  plot(tot_spec$mon_total, tot_spec[,spec_inds[i]], 
    col = countries_pal[as.numeric(droplevels(cities$country))],
    pch = countries_pch[as.numeric(droplevels(cities$country))],
    xlab = "Mean annual PM2.5", ylab = spec_names[i]
  )
  form <- as.formula(sprintf("%s ~ mon_total", 
    colnames(tot_spec)[spec_inds[i]]))
  reg <- lm(form, data = tot_spec)
  abline(reg, lwd = 2, lty = 2)
  text(par("usr")[1], par("usr")[4], sprintf("R2 = %1.2f", 
    summary(reg)$r.squared), adj = c(-0.1,1.1))
}
plot(tot_spec$mon_total, rowSums(tot_spec[,spec_inds]), 
  col = countries_pal[as.numeric(droplevels(cities$country))],
  pch = countries_pch[as.numeric(droplevels(cities$country))],
  xlab = "Mean annual PM2.5", ylab = "Constituent sum"
)
reg <- lm(rowSums(tot_spec[,spec_inds]) ~ tot_spec$mon_total, 
  data = tot_spec)
abline(a = 0, b = 1)
abline(reg, lwd = 2, lty = 2)
text(par("usr")[1], par("usr")[4], sprintf("R2 = %1.2f", 
  summary(reg)$r.squared), adj = c(-0.1,1.1))
plot.new()
legend(par("usr")[1], par("usr")[4], countries$countryname, 
  pch = countries_pch, col = countries_pal, bty = "n", xpd = T,
  ncol = 2)
dev.print(png, filename = "Results/1a_SPECcorrelations.png", 
  units = "in", res = 100)
  
  
# Correlation between each constituent and sum of constituents
x11()
par(mfrow = n2mfrow(length(spec_names) + 1))
for (i in seq_along(spec_names)){
  plot(tot_spec$SPEC_total, tot_spec[,spec_inds[i]], 
    col = countries_pal[as.numeric(droplevels(cities$country))],
    pch = countries_pch[as.numeric(droplevels(cities$country))],
    xlab = "Sum of constituents", ylab = spec_names[i]
  )
  form <- as.formula(sprintf("%s ~ SPEC_total", 
    colnames(tot_spec)[spec_inds[i]]))
  reg <- lm(form, data = tot_spec)
  abline(reg, lwd = 2, lty = 2)
  text(par("usr")[1], par("usr")[4], sprintf("R2 = %1.2f", 
    summary(reg)$r.squared), adj = c(-0.1,1.1))
}
plot.new()
legend(par("usr")[1], par("usr")[4], countries$countryname, 
  pch = countries_pch, col = countries_pal, bty = "n", xpd = T,
  ncol = 2)
dev.print(png, filename = "Results/1a_SPECcorrelationsSUM.png", 
  units = "in", res = 100)
  
# Correlation matrix
cormat <- cor(tot_spec[,spec_inds])
rownames(cormat) <- colnames(cormat) <- spec_names

x11()
corrplot(cormat, tl.pos = "d", tl.col = "black", tl.cex = 1.5)
dev.print(png, filename = "Results/1a_SPECcorrMat.png", 
  units = "in", res = 100)
  
  
