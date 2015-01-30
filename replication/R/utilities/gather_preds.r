##
##		Create data frame of all model predictions (ch and cr) and ensemble
##		Andreas Beger
##		17 July 2014
##
##		Relies on estimated model and ensemble files that are optionally saved
##		when running the replication script runme.R
##


# Directories
setwd("~/Desktop/replication")

library(reshape2)
library(spduration)
library(EBMAforecastbeta)

##
##		Load estimated models and ensemble (saved during runme.R replication run)	
##

load("data/irc_data_mod.rda")
load("data/model_estimates.rda")
load("data/ensemble_data.rda")
load("data/ensemble.rda")

ls()

# Dates:
# train = data used to estimate models
# calib + test = test data for theme model estimates
# calib = calibration period for ensemble
# test  = test period for ensemble
train.start <- as.Date("2001-03-01")
train.end   <- as.Date("2009-12-31")
calib.start <- as.Date("2010-01-01")
calib.end   <- as.Date("2012-04-30")
test.start  <- calib.end+1
test.end    <- as.Date("2014-03-31")


##
##		Functions to make model prediction data (copy theme_models.R from git/R)
##

source("R/utilities/ensemble_forecast.R")

##
##		What do we want? Data frame from data start to end of forecast period,
##		with columns for model ch, cr, and ensemble predictions + forecasts
##
##		Assemble by pieces, then sort the way we want
##

##
##		First and second pieces, predictions for training and test data.
##		This covers 2001 to 3/2014.
##

# Parameters for predData
n.models <- 7
model.names <- paste0("m", 1:7)

# For reference
long.names <- c("Leader char.", "Public disc.", "Classic", "Protest", 
	"Contagion", "Internal conflict", "Financial")

# Wrapper for predData
predsOnly <- function(data, stat) {
	res <- predData(data, stat=stat)
	res <- res[, 3:10]
	prefix <- switch(stat, "conditional hazard"="ch.", "conditional risk"="risk.")
	colnames(res) <- paste0(prefix, colnames(res))
	res
}

# Get risk and CH for training and test data for all models and ensemble
train.CH   <- predsOnly(train, "conditional hazard")
train.risk <- predsOnly(train, "conditional risk") 
test.CH    <- predsOnly(test, "conditional hazard")
test.risk  <- predsOnly(test, "conditional risk")

# Combine into wide frames for each data partition
train.preds <- cbind(train[, c("country", "date", "irr.t")], train.CH, train.risk)
test.preds  <- cbind(test[, c("country", "date", "irr.t")], test.CH, test.risk)

##
##		Third piece, forecasts through 9/2014
##

# Make frame for future period, assuming 6-months forecasts
n.ahead <- 6

# Slice of last time period to base forecasts on
data <- test[format(test$date, "%Y-%m")==format(test.end, "%Y-%m"), ]

# List of forecasts from each thematic model
inputs.raw <- lapply(1:n.models, function(x) 
	modelForecast(get(paste0("model", x)), n.ahead=n.ahead, pred.data=data))
names(inputs.raw) <- model.names[1:n.models]

# Reshape data so that we have records id'd by country and date, with columns
# for ch from each model
foo <- melt(inputs.raw, varnames=c("country", "date"), value.name="ch")
foo$L1 <- paste0("ch.", foo$L1)
fut.CH <- dcast(foo, country + date ~ L1, value.var="ch")

# Now the same for conditional risk:
#
# List of forecasts from each model
inputs.raw <- lapply(1:n.models, function(x) 
	modelForecast(get(paste0("model", x)), n.ahead=n.ahead, pred.data=data, 
		stat="conditional risk"))
names(inputs.raw) <- model.names[1:n.models]

# Reshape data so that we have records id'd by country and date, with columns
# for ch from each model
foo <- melt(inputs.raw, varnames=c("country", "date"), value.name="risk")
foo$L1 <- paste0("risk.", foo$L1)
fut.risk <- dcast(foo, country + date ~ L1, value.var="risk")

# We have two pieces; We need to add the ensemble numbers 
fut.ens <- ensembleForecast(n.ahead, n.models=7)
fut.ens <- melt(fut.ens, varnames=c("country", "date"))

# I don't think it makes sense to include an ensemble version of the risk,
# so it is set to NA below. 

# It's not in the same order as the otehr two pieces. Sort everything.
# After first converting dates to date class
fut.CH   <- fut.CH[order(fut.CH$country, fut.CH$date), ]
fut.risk <- fut.risk[order(fut.risk$country, fut.risk$date), ]
fut.ens  <- fut.ens[order(fut.ens$country, fut.ens$date), ] 

# Change value name here
colnames(fut.ens) <- c("country", "date", "ch.Ensemble")

# Combine into one df for forecasts from all and fix dates
fut.preds <- cbind(fut.ens[, 1:2], irr.t=NA, fut.ens[, 3, drop=FALSE], 
	fut.CH[, 3:9], risk.Ensemble=NA, fut.risk[, 3:9])

# Fix dates
fut.preds$date <- as.Date(paste0(fut.preds$date, "-01"))

# Check that the results sets have the same columns:
dim(train.preds)
dim(test.preds)
dim(fut.preds)

# Combine the results sets
all.preds <- rbind(train.preds, test.preds, fut.preds)
all.preds <- all.preds[order(all.preds$country, all.preds$date), ]

# Save it
save(all.preds, file="data/all_preds.rda")














