##
##		Functions for EBMA-based forecasts
##		12 May 2014
##		Andreas Beger
##

modelForecast <- function(model, n.ahead, pred.data, stat="conditional hazard") {
	# Produce n.aheah months forecasts using model
	# model   - spdur object
	# n.ahead - how many months to forecast ahead
	# date    - which data to use for forecasting?
	data  <- pred.data
	pred  <- forecast(model, pred.data=data, n.ahead=n.ahead, 
		stat=stat)
	rownames(pred) <- data$country
	months <- seq.Date(test.end+1, by="month", length.out=n.ahead)
	months <- format(months, "%Y-%m")
	colnames(pred) <- months
	pred
}

.makeAdj <- function(x, exp){
	# Bias adjustment for raw input probabilities
    .adjPred <- qlogis(x)
    .negative <- .adjPred<0
	.pos <- .adjPred>1
	.adjPred <- ((1+abs(.adjPred))^(1/exp))-1
	.miss <- is.na(.adjPred)
	.negative[.miss] <- FALSE
	.adjPred[.negative] <- .adjPred[.negative]*(-1)
	#.adjPred[.pos] <- NA
	.adjPred[.miss] <- NA
	.adjPred
}

affineTransform <- function(x, a, b) {
	# Affine transform of y = a + x*b
	y <- a + x * b
	y
}

weighInputs <- function(raw.inputs, ensemble, adjust=FALSE) {
	# Create weighted ensemble model predictions from matrix of raw
	# input probabilities. 
		# Use model parameters or not?
	if (adjust==FALSE) {
		inputs <- raw.inputs
	} else if (adjust==TRUE) {
		useModelParams <- TRUE

		# Adjust model predictions
		exp <- ensemble@exp
		inputs.adj <- lapply(inputs.raw, .makeAdj, exp)
		
		# Model parameter transformaton
		inputs.model <- NULL
		for (i in 1:n.models) {
			model.params <- ensemble@modelParams[, i, ]
			inputs.model[[i]] <- sapply(inputs.adj[[i]], affineTransform, 
			model.params[1], model.params[2])
		}
		inputs <- inputs.model
	}
	
	# Multiply by model weights
	weights <- ensemble@modelWeights
	
	# 2019-04-05: replace block below with these two lines
	input_preds <- do.call(cbind, inputs)
  ebma_preds <- (input_preds %*% weights)[, 1]
  
	# # Initialize results matrix
	# pred <- matrix(NA, nrow=nrow(inputs[[1]]), ncol=ncol(inputs[[1]]))
	# for (m in 1:ncol(inputs[[1]])) {
	# 	# Matrix of predictions from inputs for one month
	# 	month <- lapply(1:length(inputs), function(i) inputs[[i]][, m])
	# 	month <- do.call(cbind, month)
	# 	pred[, m] <- month %*% weights
	# }

	if (adjust==TRUE) {
		# Convert to probabilities
		ebma_preds <- plogis(ebma_preds)
	}

	return(ebma_preds)
}

ensembleForecast <- function(n.ahead=6, n.models=7, adjust=FALSE) {
	# Produce n.ahead months forecasts using the ensemble model
	# Implicity input: test, test.end
	# adjust - use EBMA parameters or raw transform using model weights only?
	data <- test[format(test$date, "%Y-%m")==format(test.end, "%Y-%m"), ]
	inputs.raw <- lapply(1:n.models, function(x) 
		modelForecast(get(paste0("model", x)), n.ahead=n.ahead, pred.data=data))
	names(inputs.raw) <- model.names[1:n.models]
	
	# 2019-04-05: 
	# inputs.raw originally was/is a list(7) where each element is a 164x6 matrix
	# with country rows and h/yearmonth columns
	# make this a list of vectors to match the other changes in weighInputs
	# cram back to this matrix list structure after pred so the rest doesn't break
	cc <- rownames(inputs.raw[[1]])
	ym <- colnames(inputs.raw[[1]])
	inputs.raw <- lapply(inputs.raw, as.vector)

	pred <- weighInputs(inputs.raw, ensemble, adjust=adjust)
	
	# 2019-04-05: 
	# cram back into matrix
	pred <- matrix(pred, nrow = length(cc), byrow = FALSE)

	# Pretty formatting
	rownames(pred) <- cc
	colnames(pred) <- ym

	pred
}

# function to determine optimal cutpoint for pred, obs
optCut <- function(obs, pred, statistic="f", plot=FALSE) {
	# Determine cutpoint for 0:1 pred that maximizes statistic
	
	require(ROCR)

	# Data for performance calculations
	pred.rocr <- prediction(pred, obs)

	# Calculate statistic
	perf.rocr <- performance(pred.rocr, statistic, "cutoff")
	opt.cut.no <- which.max(perf.rocr@y.values[[1]])
	opt.cut <- perf.rocr@x.values[[1]][opt.cut.no]

	# Plot if needed
	if (plot) { plot(perf.rocr) }

	return(opt.cut)		
}

# function to create data frame of id, pred, obs for train/calib/test
predData <- function(data.subset, input.names=model.names, stat="conditional hazard") {
	# Implicit input: models[1:7], ensemble
	# List of input model predictions
	inputs.raw <- lapply(1:n.models, function(x) 
		predict(get(paste0("model", x)), newdata = data.subset, stat=stat))
	names(inputs.raw) <- input.names

	# Calculate ensemble prediction
	pred <- weighInputs(inputs.raw, ensemble, adjust=FALSE)

	# Combine input and ensemble predictions with id
	inputs.df <- do.call(cbind, inputs.raw)
	colnames(inputs.df) <- names(inputs.raw)
	res <- data.frame(
		id=data.subset[, "id"], 
		obs=data.subset[, "failure"],
		Ensemble=pred,
		inputs.df,
		stringsAsFactors=FALSE)
	return(res)
}

predsStat <- function(pred.data, stat, cutpoints) {
	# Takes pred.data from predData as input
	# Loops over models in columns 3 on to calculate stat, using cutpoints
	# for each model
	res <- vector("numeric", length=(ncol(pred.data)-2))
	for (i in 1:(ncol(pred.data)-2)) {
		res[i] <- classificationStat(pred.data[, "obs"], pred.data[, (i+2)],
			stat, cutpoints[i])
		}
	return(res)
}

classificationStat <- function(obs, pred, stat, cutpoint=NULL) {
	require(ROCR)

	if (stat=="auc") {
		pred <- prediction(pred, obs)
		auc <- performance(pred, "auc")@y.values[[1]]
		return(auc)
	} else {
		
		bin.pred <- as.numeric(pred >= cutpoint)
	
		# Classifier categories
		n  <- length(obs)
		tp <- sum(bin.pred==1 & obs==1)
		fp <- sum(bin.pred==1 & obs==0)
		tn <- sum(bin.pred==0 & obs==0)
		fn <- sum(bin.pred==0 & obs==1)

		# Accuracy
		acc <- (tp + tn)/n

		# Recall and precision
		rec <- tp / (tp + fn)
		prec <- tp / (tp + fp)

		# Balanced accuracy
		spec <- tn / (fp + tn)
		ba <- 0.5*rec + 0.5*spec

		return(get(stat))
	}
}

annualizePredData <- function(pred.data) {
	require(lubridate)

	preds <- pred.data
	
	# Create country-year id
	ids <- strsplit(preds$id, " ")
	ids <- data.frame(do.call(rbind, ids), stringsAsFactors=FALSE)
	ids[, 1] <- as.Date(ids[, 1])
	preds$cyear<- paste(year(ids[, 1]), ids[, 2])

	# Create country-year data frame and merge observed events
	preds.cy <- data.frame(cyear=unique(preds$cyear))
	failure.cy <- by(preds$obs, preds$cyear, max)
	failure.cy <- data.frame(cyear=names(failure.cy), obs=as.vector(failure.cy))
	preds.cy <- join(preds.cy, failure.cy)

	# Identify non-id columsn in country-month data
	want <- ! names(preds) %in% c("id", "obs", "cyear")
	preds.mod <- preds[, want]

	# For each want column, aggregate probabilities
	for (i in 1:ncol(preds.mod)) {
		pred.cy <- by(preds.mod[, i], preds$cyear, function(x) 1 - prod(1 -x))
		pred.cy <- data.frame(cyear=names(pred.cy), pred=as.vector(pred.cy))
		colnames(pred.cy) <- c("cyear", colnames(preds.mod)[i])
		preds.cy <- join(preds.cy, pred.cy)
	}

	return(preds.cy)
}


summaryFitTable <- function(pred.data, model.names, cutpoints) {
	# Create table to model fit stats for prediction data frame
	cuts <- cutpoints

	res.table <- data.frame(
		Model=model.names,
		AUC=predsStat(pred.data, "auc", cuts),
		Cutpoint=cuts,
		Accuracy=predsStat(pred.data, "acc", cuts),
		Recall=predsStat(pred.data, "rec", cuts),
		Precision=predsStat(pred.data, "prec", cuts)
	)

	return(res.table)
}

# function for false positives that fall within x-time window
rollapply.panel <- function(data, width, FUN, x, x.id, t.id, ...) {
	require(zoo)
	data$orig.order <- 1:nrow(data)
	data <- data[order(data[, x.id], data[, t.id]), ]
	res <- by(data[, x], data[, x.id], 
		function(q) rollapply(q, width=width, FUN=FUN, partial=TRUE, 
			...)
		)
	data$res <- unlist(res)
	data <- data[order(data$orig.order), ]
	return(data$res)
}

# How many false positives fall within k/2 months of an event
fpFuzzyMatch <- function(data, pred, cutpoint, k) {
	# What fraction of false positives is withink k time units of an event
	# needs data with ccode, date, failure
	# pred is 0:1 predictions of same length as data rows
	data$bin.pred <- as.numeric(pred > cutpoint)
	data$fuzzy <- rollapply.panel(data, width=k, max, "failure", "ccode", "date")
	cat("Classification table\n")
	print(table(data[, "failure"], data[, "bin.pred"]))

	data$fp <- as.numeric(data$bin.pred==1 & data[, "failure"]==0)
	cat(paste0("\nHow many FP are within ", ceiling(k/2), " months of failure\n"))
	print(table(with(data[data$fp==1, ], fp==fuzzy)))
}