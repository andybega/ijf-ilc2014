##
##		Functions to help with theme models
##		12 May 2014
##		Andreas Beger
##

# Brier score
brier <- function(obs, pred) {
	# Calculates Brier score of model predictions "pred" given observed values
	# in "obs".
	# Input:  two vectors of same length
	# Output: scalar containing Brier score
  	res <- mean((pred-obs)^2)
  	res
}

# Model summary stats
summaryFit <- function(model, data) {
	# Calculate Brier score, AUC, percent correct, and percent reduction in 
	# error relative to an all-0 model.
	# Inputs: a spdur object ("model") and the data in which to create 
	#	      predictions ("data"). Observed values are assumed to be "failure"
	#		  in "data".
	# Output: Named 1x4 matrix with summary fit statistics

	ch <- predict(model, data=data, stat="conditional hazard")
	model.auc <- auc(data[, "failure"], ch)

	# Find optimum cutpoint
	target <- with(attr(model.auc, "roc"), sensitivities + specificities)
	optim.t <- attr(model.auc, "roc")$thresholds[target==max(target)]

	# Binary predictions
	yhat <- as.numeric(ch>optim.t)
	#yhat <- as.numeric(ch>0.5)

	# Classifier categories
	n  <- nrow(data)
	tp <- sum(yhat==1 & data[, "failure"]==1)
	fp <- sum(yhat==1 & data[, "failure"]==0)
	tn <- sum(yhat==0 & data[, "failure"]==0)
	fn <- sum(yhat==0 & data[, "failure"]==1)

	pC   <- (tp + tn)/n  # aka Accuracy
	pC.null <- sum(data[, "failure"]==0)/nrow(data)
	pre  <- (pC-pC.null)/pC.null

	# Recall and precision
	recall <- tp / (tp + fn)
	prec   <- tp / (tp + fp)

	# Balanced accuracy
	specificity <- tn / (fp + tn)
	ba <- 0.5*recall + 0.5*specificity

	res <- matrix(c(optim.t, brier(data[, "failure"], ch), model.auc, pC, 
		pre, recall, prec, ba), nrow=1)
	colnames(res) <- c("Cut point", "Brier", "AUC", "Accuracy", "pre", "Recall", 
		"Precision", "Bal. Accuracy")
	row.names(res) <- c(substitute(data))
	res
}

modelFitStats <- function(model) {
	# Calculate training and test period fit for a theme model
	res <- rbind(
		summaryFit(model, train),
		summaryFit(model, test)
		)
	res
}

table.spdur <- function(model) {
	# Format spdur estimates summary for pretty printing with xtable
	dur   <- summary(model)$duration[, c(1, 4)]
	risk  <- summary(model)$split[, c(1, 4)]
	alpha <- summary(model)$alpha[, c(1, 4)]
	res <- data.frame(rbind(dur, alpha, risk))
	res
}
