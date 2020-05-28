#####################################################################
#
#                 MCC HetPoll
#             MASTER: One file to rule them all
#
#####################################################################

setwd("C:/Users/masselpl/Documents/Recherche/2020 - LSHTM/Projects/MCC-HetPoll")

#--------------------------------
# Main analysis: current state
#--------------------------------

# Data preparation
source("Code/0_PrepData.R")

# Summary statistics
source("Code/1a_DataSummary.R") # Data in absolute value
source("Code/1b_CompositionSummary.R") # Focuses on composition

# Firts stage: city-level models
source("Code/2_FirstStage.R")

# Second stage: meta-regression
source("Code/3c_SecondStage_compositional.R")