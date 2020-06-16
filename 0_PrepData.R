#####################################################################
#
#                 MCC HetPoll
#             Part 0: Preparing data
#
#####################################################################

library(tsModel)
library(colorspace)

file_head <- "0"

#------------------------------------------
# Load data and keep cities with records
#------------------------------------------

# Load pollution data
datapath <- "C:/Users/masselpl/Documents/Recherche/2020 - LSHTM/Data/MCC"

load(sprintf("%s/MCCdata_Pollution_20200110.RData", datapath))
citiespoll <- read.csv(sprintf("%s/Metadatamccfinal20200110.csv", datapath))

# Keep only cities that have an air pollution record
subcities <- as.character(citiespoll[,2])

dlist <- dlist[cities$city %in% subcities]
cities <- cities[cities$city %in% subcities,]

subcoutries <- as.character(unique(cities$country))
countries <- countries[countries$country %in% subcoutries,]

# Then keep only cities with a PM2.5 record
citiespm25 <- citiespoll$pm25 == TRUE

dlist <- dlist[citiespm25]
cities <- cities[citiespm25,]

subcountries <- as.character(unique(cities$country))
countries <- countries[countries$country %in% subcountries,]

# Keep only some codes for each country (ask Francesco)
subcountry <- sort(c("aus8809","can8615","chi9615","chl0414", "est9715","fnl9414",
  "ger9315", "grc0110", "jap1115","mex9814", "por8018", "rom9416","sa9713", "spa9014",
  "sui9513","swe9010","twn9414","uk9016Poll","usa7306"))

dlist <- dlist[cities$country %in% subcountry]
cities <- cities[cities$country %in% subcountry,]
countries <- countries[countries$country %in% subcountry,]

# Load PM constituent data
fexp <- "PM2\\.5_SPEC_10km_buffer_[0-9]{4}"
#fexp <- "PM2\\.5_SPEC_[0-9]{4}"
flist <- list.files(sprintf("%s/MCC_PM_SPEC_New", datapath), 
  pattern = fexp)
dlist_spec <- vector("list", length(flist))
names(dlist_spec) <- as.character(2003:2017)
for (i in seq_along(flist)){
  dlist_spec[[i]] <- read.csv(sprintf("%s/MCC_PM_SPEC_New/%s", 
    datapath, flist[i]))
  dlist_spec[[i]]$year <- names(dlist_spec)[i]
  rownames(dlist_spec[[i]]) <- NULL  
}
ddf_spec <- do.call(rbind, dlist_spec)
dlist_spec <- split(ddf_spec, f = ddf_spec$city)

# Load indicator data
load(sprintf("%s/mcc_indicators_20200504.RData", datapath)) # Socio-economic
load(sprintf("%s/MCC_indicators_UCD_20200504.RData", datapath)) # Environment
mcc.indicators <- merge(mcc.indicators, final.ucd.mcc, 
  by.x = "city", by.y = "citiescode")

# Select common cities
commoncities <- intersect(names(dlist_spec), cities$city)
commoncities <- Reduce(intersect, 
  list(names(dlist_spec), cities$city, mcc.indicators$city))

dlist <- dlist[names(dlist) %in% commoncities]
dlist_spec <- dlist_spec[names(dlist_spec) %in% commoncities]
dlist_spec <- dlist_spec[names(dlist)]
cities <- cities[cities$city %in% commoncities,]
countries <- countries[countries$country %in% cities$country,]
mcc.indicators <- mcc.indicators[mcc.indicators$city %in% commoncities,]

capture.output(print("Cities with pollution record"), table(cities$country),
  file = "Results/0_cities.txt")

#------------------------------------------
#      Compute relevant variables
#------------------------------------------
pmL <- 1
trim <- .05

for (i in seq_len(nrow(cities))){
  # Define the mortality variable (all-cause or non-external depending on the
  #   country)
  if (is.null(dlist[[i]]$all)){
    dlist[[i]]$death <- dlist[[i]]$nonext
  } else {
    dlist[[i]]$death <- dlist[[i]]$all
  }
  
  # Define PM moving average
  dlist[[i]]$mapm <- runMean(dlist[[i]]$pm25, 0:pmL)
  
  # Define trimmed pm2.5
  quants <- quantile(dlist[[i]]$pm25, c(trim, 1 - trim), na.rm = T)
  dlist[[i]]$tpm25 <- dlist[[i]]$pm25
  dlist[[i]]$tpm25[dlist[[i]]$tpm25 < quants[1] | 
    dlist[[i]]$tpm25 > quants[2]] <- NA
    
  # PM2.5 annual mean in the PM constituents dataset
  annmean <- aggregate(dlist[[i]]$pm25, 
    by = list(year = dlist[[i]]$year), mean, na.rm = T)
  dlist_spec[[i]]$mon_total <- annmean$x[
    match(dlist_spec[[i]]$year, annmean$year)]
  
  # Sum of constituents
  spec_inds <- grep("PM25", colnames(dlist_spec[[i]]))
  dlist_spec[[i]]$SPEC_total <- rowSums(dlist_spec[[i]][,spec_inds])
}

#------------------------------------------
#     Filter cities on data quality
#------------------------------------------
dlist_orig <- dlist
addyears <- 1999:2002 # Additional years for 1st stage to have longer records (especially for USA for which data only covers until 2006)
for(i in seq(nrow(cities))) {
  # Set outliers to NAs in order to exclude them from the analysis
  dlist[[i]]$death[dlist[[i]]$outlierm == 1] <- NA
  dlist[[i]]$tmean[dlist[[i]]$outliert == 1] <- NA
  
  # Keep only periods with available PM2.5 data
  dr <- range(which(!is.na(dlist[[i]]$pm25)))
  dlist[[i]] <- dlist_orig[[i]] <- dlist[[i]][dr[1]:dr[2],]
  
  # Keep only common years between dlist and dlist_spec
  # N.B. Here I keep incomplete years as well 
#  commonyears <- intersect(unique(dlist[[i]]$year), 
#    dlist_spec[[i]]$year)
#  dlist[[i]] <- dlist[[i]][dlist[[i]]$year %in% 
#    union(addyears, commonyears),]
#  dlist_spec[[i]] <- dlist_spec[[i]][
#    dlist_spec[[i]]$year %in% commonyears,]
  # We consider that pollution does not vary in time for recent years.
  # Thus we keep all years of species.
  # However we keep pm2.5 data only post-1999
  dlist[[i]] <- dlist[[i]][dlist[[i]]$year >= 1999,]
  
  # Removing lines containing only zeros in constituents
  spec_inds <- grep("PM25", colnames(dlist_spec[[i]]))
  dlist_spec[[i]] <- dlist_spec[[i]][
    apply(dlist_spec[[i]][,spec_inds] != 0, 1, any),]
}

# Remove badly linked values from MCC indicators
mcc.indicators$TotPop[mcc.indicators$TotPop.source != 1] <- NA
mcc.indicators$Density[mcc.indicators$Density.source != 1] <- NA
mcc.indicators$UrbArea[mcc.indicators$UrbArea.source != 1] <- NA
mcc.indicators$GreenArea[mcc.indicators$GreenArea.source != 1] <- NA
mcc.indicators$PopCore[mcc.indicators$PopCore.source != 1] <- NA
mcc.indicators$Sprawl[mcc.indicators$Sprawl.source != 1] <- NA


#---- Exclusion according to a number of criteria
exclusion <- list()
# Cities that have less than ny years of consecutive data
ny <- 1
exclusion$notEnoughYears <- sapply(dlist, nrow) < (ny * 365)

# Cities with a proportion > propm missing death data
propm <- 0.5
exclusion$NAmortality <- sapply(dlist, function(d) mean(is.na(d$death))) > propm

# Cities with a proportion > propp missing PM2.5 data
propp <- 0.95
exclusion$NApoll <- sapply(dlist, function(d) mean(is.na(d$mapm))) > propp

# Cities with a proportion > propt missing PM2.5 data
propt <- 0.95
exclusion$NAtemp<- sapply(dlist, function(d) mean(is.na(d$tmean))) > propt

# Cities with no overlap between pollution and SPEC datasets
exclusion$noOverlap <- sapply(dlist_spec, nrow) == 0

#---- Check excluded cities and perform exclusion
capture.output(print("Excluded cities because of missing values"),
  lapply(exclusion, function(e) cities[e,c("cityname", "countryname")]),
  file = "Results/0_cities.txt", append = TRUE
)

# Finally proceed to exclusion
exclusion <- as.data.frame(exclusion)
excl <- apply(exclusion, 1, any)

dlist <- dlist[!excl]
dlist_spec <- dlist_spec[!excl]
dlist_orig <- dlist_orig[!excl]
cities <- cities[!excl,]
countries <- countries[countries$country %in% cities$country,]
mcc.indicators <- mcc.indicators[!excl,]

# Final count
citycount <- table(cities$country)
capture.output(print("Final count"), citycount[citycount > 0], 
  sprintf("Total: %i", sum(citycount)),
  file = "Results/0_cities.txt", append = TRUE
)

save.image("Data/0_Data.RData")