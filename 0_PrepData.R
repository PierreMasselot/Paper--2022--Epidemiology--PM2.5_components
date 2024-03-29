#####################################################################
#
#                 MCC HetPoll
#             Part 0: Preparing data
#
#####################################################################

library(tsModel)

# Workspace directory must be set at the root of the repo

# Data are not shared. Should be stored in a subfolder "Data"

#------------------------------------------
# Load data and keep cities with records
#------------------------------------------

#----- Load MCC pollution data
load("Data/MCCdata_Pollution_20200504.RData")
citiespoll <- read.csv("Data/Metadatamccfinal20200504.csv")

# Keep only cities that have a PM2.5 record
pmcities <- subset(citiespoll, pm25, citycode, drop = T)

dlist <- dlist[pmcities]
cities <- cities[cities$city %in% pmcities,]

subcountries <- as.character(unique(cities$country))
countries <- countries[countries$country %in% subcountries,]

# Select country codes
subcountry <- sort(c("aus8809","can8615","chi9615","chl0414", "est9715",
  "fnl9414", "ger9315", "grc0110", "jap1115","mex9814", "nor6918", "per0814", 
  "por8018", "rom9416","sa9713", "spa9014", "sui9513","swe9010","twn9414",
  "uk9016Poll", "usa7306"))

cities <- cities[cities$country %in% subcountry,]
dlist <- dlist[cities$city]
countries <- countries[countries$country %in% subcountry,]

#----- Load PM constituent data

# Select files
fexp <- "PM2\\.5_SPEC_10km_buffer_[0-9]{4}"
flist <- list.files("Data/MCC_PM_SPEC_New", pattern = fexp)

# Loop on files to read and store them
dlist_spec <- vector("list", length(flist))
names(dlist_spec) <- as.character(2003:2017)
for (i in seq_along(flist)){
  dlist_spec[[i]] <- read.csv(sprintf("Data/MCC_PM_SPEC_New/%s", flist[i]))
  dlist_spec[[i]]$year <- names(dlist_spec)[i]
  rownames(dlist_spec[[i]]) <- NULL  
}

# Bind them in a big data.frame
ddf_spec <- do.call(rbind, dlist_spec)

# And finally split by city
dlist_spec <- split(ddf_spec, f = ddf_spec$city)

#----- Load indicator data

# Socio-economic indicators (OECD)
load("Data/mcc_indicators_20200504.RData")

# Environmental indicators (Urban Centre Database)
load("Data/MCC_indicators_UCD_20200504.RData")

# Merge them together
mcc.indicators <- merge(mcc.indicators, final.ucd.mcc, 
  by.x = "city", by.y = "citiescode")

# Select indicators of interest
indic_names <- c("oldpopprop", "GDP", "avgtmean", 
  "totalrange", "E_GR_AV00", "E_GR_AV14", "B00", "B15")
mcc.indicators <- mcc.indicators[,c("city", indic_names)]

# Merge with city data.frame
cities <- merge(cities, mcc.indicators, all = F)

#----- Linkage of all datasets
  
# Select common cities
commoncities <- intersect(names(dlist_spec), cities$city)

# Keep only common cities and reorder
cities <- cities[cities$city %in% commoncities,]
cities <- cities[with(cities, order(country, city)),]

# Select cities in other objects
dlist <- dlist[cities$city]
dlist_spec <- dlist_spec[cities$city]
countries <- countries[countries$country %in% cities$country,]

#------------------------------------------
#      Compute relevant variables
#------------------------------------------

# Parameters
pmL <- 1
trim <- .05

#---- Loop on cities to compute variables
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

#----- Curate time series
for(i in seq(nrow(cities))) {
  # Set outliers to NAs in order to exclude them from the analysis
  dlist[[i]]$death[dlist[[i]]$outlierm == 1] <- NA
  dlist[[i]]$tmean[dlist[[i]]$outliert == 1] <- NA
  
  # Keep only periods with available PM2.5 data
  dr <- range(which(!is.na(dlist[[i]]$pm25)))
  dlist[[i]] <- dlist_orig[[i]] <- dlist[[i]][dr[1]:dr[2],]
  
  # Keep 1999 and after althougb spec starts in 2003 to have longer time series
  dlist[[i]] <- dlist[[i]][dlist[[i]]$year >= 1999,]
  
  # Removing lines containing only zeros in constituents
  spec_inds <- grep("PM25", colnames(dlist_spec[[i]]))
  dlist_spec[[i]] <- dlist_spec[[i]][
    apply(dlist_spec[[i]][,spec_inds] != 0, 1, any),]
}

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

# Remove incomplete data
exclusion$missingMCCindicators <- !complete.cases(cities[,indic_names])

# Finally proceed to exclusion
exclusion <- as.data.frame(exclusion)
excl <- apply(exclusion, 1, any)

dlist <- dlist[!excl]
dlist_spec <- dlist_spec[!excl]
dlist_orig <- dlist_orig[!excl]
cities <- cities[!excl,]
countries <- countries[countries$country %in% cities$country,]

save.image("Data/0_Data.RData")