#' Decompose deviance
#'
#' \code{devDecomp} will decompose between and within group deviance for a 
#' variable. It returns a data frame with the group and variable as well
#' as the group means, and within and between deviance.

devDecomp <- function(group, var) {
	# Decompose within and between group deviation
	df <- data.frame(group=group, var=var)
	df$group.mean <- ave(df$var, df$group)
	df$within <- with(df, var - group.mean)
	df$between <- with(df, group.mean - mean(var))
	df
}

#' Decompose variance
#'
#' \code{varDecomp} will decompose the between and within group sum of squared
#' errors (or deviations) for a variable.
#'

varDecomp <- function(group, var) {
	# Calculate total sum of squares given df with var, within and between dev.
	df <- devDecomp(group, var)
	ss.total   <- sum(df$within^2) + sum(df$between^2)
	ss.within  <- sum(df$within^2)
	ss.between <- sum(df$between^2)
	res <- c(total=ss.total, within=ss.within, between=ss.between)
	res
}