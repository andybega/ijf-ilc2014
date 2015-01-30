#
#   Replication for 
#	"Irregular Leadership Changes in 2014: Forecasts using ensemble, split-
#    population duration models", International Journal of Forecasting
#
#	30 January 2015
#	Andreas Beger (andreas.beger@duke.edu)
#

# The "replication" folder, alongside with this script, contains 2 R packages
# that are needed but not available on CRAN, several functions, and several
# data files. The data files include the main data, "irc.data", as well 
# as intermediate steps in the results like saved theme model estimates, 
# the ensemble model data object and ensemble model estimates themselves.
# 
# Results may differ slightly by system, but we have included intermediate
# estimates/data needed to exactly replicate the results reported.
#


##------------------------------------------------------------------------------
##
##		1. Set up directories and install/load needed packages.
##	
##		   You may have to install missing packages from CRAN manually.
##
##------------------------------------------------------------------------------

rm(list=ls())

##
##	CHANGE: path
##

# Set to path containing replication folder
setwd("~/Desktop/replication")

# Install packages in R/packages if not already done
if (!is.element("EBMAforecastbeta", installed.packages()[, 1])) {
	install.packages("R/packages/EBMAforecastbeta_0.44.zip", repos=NULL, 
		type="source")
}
if (!is.element("spduration", installed.packages()[, 1])) {
	install.packages("R/packages/spduration_0.12.zip", repos=NULL, 
		type="source")
}


# Load required libraries
library(cshapes)
library(EBMAforecastbeta)
library(ggplot2)
library(lubridate)
library(maptools)
library(plyr)
library(pROC)
library(reshape2)
library(RColorBrewer)
library(ROCR)
library(sbgcop)
library(scales)
library(sp)
library(spduration)
library(xtable)
library(zoo)

# Source helper functions
source("R/utilities/theme_models.R")
source("R/utilities/ensemble_forecast.R")
source("R/utilities/worldMap.R")


##------------------------------------------------------------------------------
##
##		Figure 1: Map of ILCs between 2001 and 2014.
##
##------------------------------------------------------------------------------

# load data
load(file="data/irc-data-v3.rda")

# Aggregate by country
irc.by.country <- subset(irc.data, select=c('ccode','irr.t'))
irc.by.country <- irc.by.country[complete.cases(irc.by.country),]
irc.by.country <- with(irc.by.country, aggregate(irr.t, by=list(ccode), FUN=sum))
colnames(irc.by.country)<- c('ccode', 'irr.t')

# Plot
dpi <- 300
png("graphics/map_ILC.png", width=3*dpi, height=1.26*dpi, pointsize=20)
data <- irc.by.country
id <- "ccode"
x  <- "irr.t"
world <- cshp(date=as.Date("2012-01-01"))
world@data <- data.frame(world@data, data[match(world@data[, 'GWCODE'], data[, id]), ])

# Set fill colors
colorpal <- rev(brewer.pal(max(data[, x]), 'Reds'))
colors <- ifelse(is.na(world@data[, x])==T, '#B0B0B0', colorpal[match(world@data[, x], sort(unique(world@data[, x]), decreasing=T))])
    
# Plot map
par(mar=c(1, 1, 1, 1))
plot(world, col='gray30', border='gray30', lwd=1)
plot(world, col=colors, border=F, add=T)
    
# Legend
legend.text <- c('No data', rev(unlist(dimnames(table(world@data[, x])))))
legend(x=-170, y=15, legend=legend.text, fill=c('#B0B0B0', colorpal), 
       bty='n')
dev.off()


##
##		Positive rates in data
##

ilc <- irc.data[, c("id", "irr.t")]

# Create country-year id
ids <- strsplit(ilc$id, " ")
ids <- data.frame(do.call(rbind, ids), stringsAsFactors=FALSE)
ids[, 1] <- as.Date(ids[, 1])
ilc$cyear<- paste(year(ids[, 1]), ids[, 2])

# Create country-year data frame and merge observed events
ilc.cy <- by(ilc$irr.t, ilc$cyear, max)
ilc.cy <- data.frame(cyear=names(ilc.cy), ilc=as.vector(ilc.cy))

# Tabulate positives in annual data
table(ilc.cy$ilc)
sum(ilc.cy$ilc)/nrow(ilc.cy)*1000  # positive rate /1000

# Tabulare positives in monthly data
table(irc.data$irr.t)
sum(irc.data$irr.t)/nrow(irc.data)*1000  # positive rate /1000

##------------------------------------------------------------------------------
##
##		Figure 1: Map of ILCs between 2001 and 2014.
##
##		Note: this was not in the IJF version we sent out in September '14
##		AB added it after Cassy pointed out the ILC number reported, 46, was 
##		wrong on 11/5/2014.
##
##------------------------------------------------------------------------------

# List all cases of IRC
cond <- irc.data$irr.t==1
casesAll <- irc.data[cond, c("country", "date", "goelname", "irr.exit", 
	"irr.entry", "mths.in.power", "duration")]

casesIrr <- na.omit(casesAll)
casesIrr$date <- format(casesIrr$date, "%Y-%m")
casesIrr <- casesIrr[order(casesIrr$date, casesIrr$country), ]
colnames(casesIrr) <- c("Country", "Date", "Leader", "Irr. Exit", "Irr. Entry",
	"Mths in power", "TTF")
print(xtable(casesIrr, digits=0), include.rownames=FALSE)


##------------------------------------------------------------------------------
##
##    	Set up data for modeling:
##		   - Impute missing values for the monthly variables from PITF
##		   - Exclude some cases missing due to lack of coverage in CRISP
##			 (Serbia and South Sudan, as well as smaller islands)
##		   - Set dates defining training, calibration, and test periods
##		   - Create training, test data sets for theme models
##
##------------------------------------------------------------------------------


# Impute missing
Y <- irc.data[, c(1,226:235)]
leidos.imputed <- sbgcop.mcmc(Y,seed=123)
plot(leidos.imputed)
new.Y<-data.frame(leidos.imputed$Y.pmean)
colnames(new.Y)<-c("ccode.i","prsincon.i","prsexcon.i","prsmilit.i",
                   "prsethtn.i","eiacoilt.i","ifsinrsv.i","ifs__cpi.i",
                   "iloaunrb.i","dpidtleg.i","dpidtexe.i")

new.Y$dpidtexe.i <- round (new.Y$dpidtexe.i)
new.Y$prsethtn.i <- round (2*new.Y$prsethtn.i)/2
new.Y$dpidtleg.i <- round (new.Y$dpidtleg.i)
new.Y$prsmilit.i <- round (2*new.Y$prsmilit.i)/2
new.Y$prsexcon.i <- round (2*new.Y$prsexcon.i)/2
new.Y$prsincon.i <- round (2*new.Y$prsincon.i)/2     

new.Y <- new.Y[, -1]
irc.data <- cbind(irc.data, new.Y)

# Subset some countries missing due to CRISP
#irc.data <- irc.data[!is.na(irc.data$insurgency.l1), ]
irc.data <- irc.data[!is.na(irc.data$exclpop.l1), ]
irc.data <- irc.data[!is.na(irc.data$W.gower.econ.opp_resistance.l1), ]
irc.data <- irc.data[!is.na(irc.data$EXREC.l1),]


# irc.data  is from 2001 to 2013
# base.data is from 1955 to 2013

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

# Tabulate failures
table(irc.data$failure[with(irc.data, date>=train.start & date<=train.end)])
table(irc.data$failure[with(irc.data, date>=calib.start & date<=calib.end)])
table(irc.data$failure[with(irc.data, date>=test.start & date<=test.end)])

# Subset training and test data
# We have to manually set the end.spell indicator
train <- irc.data[irc.data$date>=train.start & irc.data$date<calib.start, ]
train$end.spell[train$date==max(train$date)] <- 1

test <- irc.data[irc.data$date>=calib.start & irc.data$date<=test.end, ]
test$end.spell[test$date==max(test$date)] <- 1

# Optionally save the processed data
save(train, test, irc.data, file="data/irc_data_mod.rda")
#load("data/irc_data_mod.rda")

##------------------------------------------------------------------------------
##
##		Figure 2: Variance decomposition, all variables
##
##------------------------------------------------------------------------------


# Subset some countries missing due to CRISP
#irc.data <- irc.data[!is.na(irc.data$insurgency.l1), ]
var.data <- irc.data[!is.na(irc.data$DEMOC.l1), ]

var.data$i.conf.GOVtGOV.l1 <- with(var.data, i.verb.conf.GOVtGOV.l1 + 
	i.matl.conf.GOVtGOV.l1)
var.data$i.coop.GOVtGOV.l1 <- with(var.data, i.verb.coop.GOVtGOV.l1 + 
	i.verb.coop.GOVtGOV.l1)

library(grid)
library(plyr)

source("R/utilities/varDecomp.R")

var.names <- names(var.data)[c(14, 19:225)]

ss <- lapply(var.names, function(x) varDecomp(group=var.data$ccode, 
	var=var.data[, x]))
ss <- do.call("rbind", ss)
ss <- data.frame(variable=var.names, ss)

ss$within.p <- ss$within/ss$total
ss$between.p <- ss$between/ss$total

ss$type <- ""
ss$type[c(1:42, 206:208)] <- "Structural"
ss$type[c(43:67, 178:205)] <- "Behavioral"
ss$type[68:111] <- "Spatial weights, Centroid or knn"
ss$type[112:177] <- "Spatial weights, Gower D"

ss$type <- factor(ss$type, levels=c("Structural", "Behavioral", "Spatial weights, Centroid or knn",
	"Spatial weights, Gower D"))

ss$variable <- as.character(ss$variable)
want  <- c("NY.GDP.PCAP.KD.l1", "..AA.ZF....l1", "ST.INT.RCPT.CD.l1", 
	"DEMOC.l1", "Oil.Futures.l1", "FP.CPI.TOTL.l1", "MS.MIL.XPND.GD.ZS.l1",
	"i.conf.GOVtGOV.l1", "W.knn4.std.protest.tALL.l1", "BX.KLT.DINV.CD.WD.l1",
	"W.gower.pol.protest.tALL.l1", "SP.DYN.LE00.FE.IN.l1")
highlight <- ss[match(want, ss$variable), ]
highlight$label <- c("GDP p.c.", "Exch. rate to SDR", 
	"Tourism receipts", "Polity DEM", "Oil Futures", "CPI", "Mil. exp.", 
	"Intra-gov. conflict", "knn4 protest", "FDI", "Gower Polity protest",
	"Life expt., female")

y.min <- max(1, round_any(min(ss$total), 100, floor))
p <- ggplot(data=ss, aes(x=total, y=between.p)) +
	geom_abline(intercept=0.5, slope=0, size=0.5, col="cyan", alpha=0.5) +
	geom_point(aes(color=type), size=2, alpha=0.5) +
	scale_x_log10(limits=c(y.min, 2e+27), minor_breaks=c(0:27), breaks=c(10^seq(from=0, to=27, by=5)), labels=expression(10^0, 10^5, 10^10, 10^15, 10^20, 10^25)) +
	scale_y_continuous(limits=c(0, 1)) +
	labs(x="Total variance", y="Between/Total", col="") +
	theme_bw() +
	theme(plot.title = element_text(vjust=2), legend.position="bottom",
		legend.direction="vertical", legend.margin=unit(-0.6,"cm")) +
	guides(col = guide_legend(ncol = 2)) 
p <- p + geom_text(data=highlight, aes(x=total, y=between.p, label=label),
	alpha=0.7, size=4) 
ggsave(file="graphics/var_all.png", p, width=5, height=5.5, dpi=300)

# Add iso-lines showing variables with equal within variance

iso_within <- function(x0) { 
  # For x0, returns points accross between/total raio with equal
  # within variance
  # x0 = constant within
  # y  = between/total ratio, y coordinate
  # total = updated total, x coordinate
  Ny = 1000
  Nx = length(x0)
  res <- data.frame(x0=vector("numeric", Nx*Ny), 
  	x=vector("numeric", Nx*Ny), y=vector("numeric", Nx*Ny))

  y <- seq(0, 0.999, length.out=Ny)
  for (i in seq_along(x0)) {
  	idx0 <- (i-1) * Ny + 1
  	idx1 <- i * Ny
    res[idx0:idx1, ] <- data.frame(
    	x0 = x0[i], 
    	x = x0[i] + y * x0[i] / (1 - y),
    	y = y)
  }
  res
}

xy <- iso_within(10^seq(0, 25, by=5))
p + geom_line(data=xy, aes(x=x, y=y, group=factor(x0)),
	colour="gray50", size=0.2)

##-----------------------------------------------------------------------------
##
##    	Figure 3: Plot risk versus hazard over time for Thailand
##		For this, we load the model estimates file, the models are
##		re-estimated below.
##
##------------------------------------------------------------------------------


# These are created later on from the ensemble and theme forecasts
load("data/all_preds.rda")


dpi <- 300
png("graphics/thailand.png", height=3*dpi, width=4*dpi, pointsize=15)
cex <- 2
one.country <- all.preds[all.preds$country=="Thailand", ]
n.thm.mdls <- 7  # How many theme models are there?

# Pretty lables for x-axis
x.at <- c("2001-03-01", "2004-01-01", "2007-01-01", "2010-01-01", 
	"2013-01-01", "2014-09-01")
x.at <- as.Date(x.at)
x.labs <- format(x.at, "%Y-%m")

layout(matrix(c(1,2), 2, 1), heights=c(3, 2))

par(mar=c(3, 3.5, 2, 2), cex=cex, las=1)

# Top plot with conditional hazard

plot(one.country$date, one.country$ch.Ensemble, type="l", col="darkblue", 
	lwd=cex*2, bty="n", xlab="", ylab="", cex.axis=1, xaxt="n", ylim=c(0, 0.08))
theme.mdls.ch <- one.country[, paste0("ch.m", 1:n.thm.mdls)]
matplot(one.country$date, theme.mdls.ch, type="l", add=TRUE, 
	col=alpha("gray50", 0.3), lty=1)

# Fancify plot
grid(nx=NA, ny=NULL)
title(main="Conditional Hazard", adj=0)
axis(1, at=x.at, labels=x.labs)
rect(test.end+1, 0, max(all.preds$date), 0.08, density=3)  # mark forecast period

# Mark ILC that occurred
abline(v=one.country$date[one.country$irr.t==1], col="darkred", lwd=cex*2)
text(as.Date("2006-11-01"), y=0.04, adj=c(0, 0), labels=c("IRC occured in\nSeptember 2009"))

# ILC in forecast period on 22 May 2014
abline(v=as.Date("2014-05-22"), col="darkred", lwd=cex*2)

# Bottom plot with conditional risk
plot(one.country$date, one.country$risk.Ensemble, type="l", col="darkblue", 
	lwd=cex*2, bty="n", xlab="", ylab="", ylim=c(0, 1), xaxt="n")
theme.mdls.cr <- one.country[, paste0("risk.m", 1:n.thm.mdls)]
matplot(one.country$date, theme.mdls.cr, type="l", add=TRUE, 
	col=alpha("gray50", 0.3), lty=1)	

# Fancify plot
title(main="Conditional Risk", adj=0)
grid(nx=NA, ny=NULL)
axis(1, at=x.at, labels=x.labs)
rect(test.end+1, 0, max(all.preds$date), 1, density=3)  # mark forecast period

dev.off()


##------------------------------------------------------------------------------
##
##    	Estimate thematic models
##
##------------------------------------------------------------------------------


##------------------------------------------------------------------------------
##
##		Supplemental Table 1: Leader characteristics
##
##------------------------------------------------------------------------------

model1 <- spduration::spdur(
	"duration ~ 1 + log10(i.matl.conf.DIStGOV.l1+1) + 
		log10(i.matl.coop.GOVtGOV.l1+1)",
	"atrisk ~ 1 + ldr.irregular + ldr.foreign + log10(mths.in.power+1)",
	data=train)

modelFitStats(model2)
xtable(table.spdur(model1), caption="Leader characteristics model estimates",
	label="tab:chatter")

##------------------------------------------------------------------------------
##
##		Supplemental Table 2: Public discontent
##
##------------------------------------------------------------------------------


model2 <- spduration::spdur(
	"duration ~ 1 + log10(i.verb.coop.GOVtGOV.l1+1) + 
		log10(i.verb.conf.GOVtDIS.l1+1) + log10(i.verb.conf.DIStGOV.l1+1) + 
  		log10(i.protest.tGOV.l1+1)",
  	"atrisk ~ 1 + IT.NET.USER.P2.l1 + IT.CEL.SETS.P2.l1 + log10(exclpop.l1+1)",
  	data=train)

modelFitStats(model2)
xtable(table.spdur(model2), caption="Public interactions model estimates",
	label="tab:chatter")


##------------------------------------------------------------------------------
##
##		Supplemental Table 3: Global instability models
##   
##------------------------------------------------------------------------------

model3 <- spduration::spdur(
	"duration ~ 1 + W.knn4.std.ins.l.count.both.l1 + SP.DYN.LE00.FE.IN.l1 + 
		PARCOMP.l1 ",
	"atrisk ~ 1 + log10(NY.GDP.PCAP.KD.l1) + exclpop.l1",
	data=train)

modelFitStats(model3)
xtable(table.spdur(model3), caption="Global instability model estimates",
	label="tab:pitf5")


##------------------------------------------------------------------------------
##
##		Supplemental Table 4: Protest model
##
##------------------------------------------------------------------------------

model4 <- spduration::spdur(
	"duration ~ 1 + eth.rel.l.count.l1 + reb.l.count.both.l1 + protest.tALL.l1 + 
	W.gower.pol.reb.l.count.both.l1",
	"atrisk ~ 1+ dom.cris.i.count.l1 + log10(MS.MIL.XPND.GD.ZS.l1)",
	data=train)

modelFitStats(model4)
xtable(table.spdur(model4), caption="Protest model estimates")


##------------------------------------------------------------------------------
##
##		Supplemental Table 5: Contagion model
##
##------------------------------------------------------------------------------

model5 <- spduration::spdur(
	"duration ~ 1 + log10(W.centdist.std.opp_resistance.l1+1) +
  		log10(W.centdist.std.repression.l1+1)",
	"atrisk ~ 1 + Amnesty.l1 + log10(ProxElection.l1+1) +
  		log10(opp_resistance.l1+1) + log10(SP.POP.TOTL.l1)",
  	data=train)


modelFitStats(model5)
xtable(table.spdur(model5), caption="Contagion model estimates",
	label="tab:contagion")


##------------------------------------------------------------------------------
##
##		Supplemental Table 6: Internal conflict
##
##------------------------------------------------------------------------------

model6 <- spduration::spdur(
  "duration ~ 1 +   log(intratension.l1 *i.protest.tGOV.l1*IT.CEL.SETS.P2.l1 + 1 )+
  log(intratension.l1+1) +
  log(i.protest.tGOV.l1+1) +
  log(IT.CEL.SETS.P2.l1+1) ",
  "atrisk ~ 1 + log10(NY.GDP.PCAP.KD.l1)  +ProxElection.l1 + AUTOC.l1 ",
  data=train)

modelFitStats(model6)
xtable(table.spdur(model6), caption="Internal conflict model estimates",
	label="tab:internal")


##------------------------------------------------------------------------------
##
##		Supplemental Table 7: Financial stability
##
##------------------------------------------------------------------------------


model7 <- spduration::spdur(
 "duration ~  log(1+ifs__cpi.i) + log(1+ifsinrsv.i)",
 "atrisk ~ 1 + log10(NY.GDP.PCAP.KD.l1) + Amnesty.l1 + log10(ProxElection.l1+1) +
   	log10(opp_resistance.l1+1) + log10(SP.POP.TOTL.l1)",
 data=train)
summary(model7)
modelFitStats(model7)
xtable(table.spdur(model7), caption="Financial Stability model estimates",
	label="tab:fin")

# Optionally, save model estimates. Can also load Duke model estimates if 
# you want to skip the model estimation above.
save(model1, model2, model3, model4, model5, model6, model7, 
	file="data/model_estimates.rda")


##------------------------------------------------------------------------------
##
##    	Calibrate ensemble model
##
##------------------------------------------------------------------------------


#load("data/model_estimates.rda")

# Set number of models and their names. We will need this several times
# for the code below.
n.models <- 7
model.names <- c("Leader char.", "Public disc.", "Global Inst.", "Protest", 
	"Contagion", "Internal conflict", "Financial")


##		Calculate theme model predictions
pr.out  <- data.frame(observed=test$failure)
for (i in 1:n.models) {
	a.model <- get(paste0("model", i))
	p.name  <- paste0("pred.", i)
	print(p.name)
	pr.out[, p.name] <- predict(a.model, data=test, stat="conditional hazard")
	# fix exact 0's to slightly above 0, otherwise EBMA will not work
	pr.out[, p.name] <- replace(pr.out[, p.name], pr.out[, p.name]<=0, 1e-19)
}

# Subset calibration and test period predictions
pr.calib <- pr.out[test$date>=calib.start & test$date<=calib.end, ]
pr.test  <- pr.out[test$date>=test.start  & test$date<=test.end,  ]

# Create ensemble data object
ensemble.data <- makeForecastData(
	.predCalibration 	= pr.calib[, 2:(n.models+1)],
	.outcomeCalibration = pr.calib[, "observed"],
	.predTest 			= pr.test[, 2:(n.models+1)],
	.outcomeTest 		= pr.test[, "observed"],
	.modelNames=model.names
	)

# Optionally, save or load ensemble data object
#save(ensemble.data, file="data/ensemble_data.rda")


# Calibrate ensemble model
#load("data/ensemble_data.rda")
ensemble <- calibrateEnsemble(ensemble.data, model="logit", maxIter=25000, 
	exp=2, const=0.0001)

# This table will show the fit statitics for the ensemble and component models,
# as well as the model weights. The constant and slope parameters are 
# optionally used in EBMA to aggregate theme model predictions, however we 
# do not use them here.
summary(ensemble)

# Optionally, save or load fitted ensemble object
#save(ensemble.data, ensemble, file="data/ensemble.rda")
#load("data/ensemble.rda")



##------------------------------------------------------------------------------
##
##    	Table 1: Ensemble model, monthly observations
##
##------------------------------------------------------------------------------


##
##		Calculate AUC, accuracy, prec., recall, ROC for train/calib/test
##		for all models and ensemble
##

# Create data frames with model predictions and observed data
in.sample.df <- rbind(train, test[with(test, date>=calib.start & date<=calib.end), ])
in.sample.preds <- predData(in.sample.df)
oo.sample.preds <- predData(test[with(test, date>=test.start & date<=test.end), ])


##
##		Accuracy in country-month data
##

maxF <- function(pred, obs, data=NULL) {
	# Wrapper for ROCR that returns maximum F score as numeric
	if (!is.null(data)) {
		pred <- eval(substitute(pred), envir=data)
		obs  <- eval(substitute(obs), envir=data)
	}
	rocr.df <- prediction(pred, obs)
	perf <- performance(rocr.df, "f")
	res  <- max(perf@y.values[[1]], na.rm=T)
	return(res)
}

auc2 <- function(pred, obs) {
	# Wrapper for ROCR that returns AUC as numeric
	df <- prediction(pred, obs)
	perf <- performance(df, "auc")
	res  <- perf@y.values[[1]]
	return(res)
}

cm.fit.table <- data.frame(Model=c("Ensemble", model.names), W=NA, AUC=NA, F=NA,
	AUC=NA, F=NA)
cm.fit.table[, "W"] <- c(NA, ensemble@modelWeights)
cm.fit.table[, 3] <- apply(in.sample.preds[3:10], 2, auc2, obs=in.sample.preds$obs)
cm.fit.table[, 4] <- apply(in.sample.preds[3:10], 2, maxF, obs=in.sample.preds$obs)
cm.fit.table[, 5] <- apply(oo.sample.preds[3:10], 2, auc2, obs=oo.sample.preds$obs)
cm.fit.table[, 6] <- apply(oo.sample.preds[3:10], 2, maxF, obs=oo.sample.preds$obs)

# Print for latex
print(xtable(cm.fit.table, digits=2, label="tab:ensemble", 
	caption="Ensemble model"), include.rownames=FALSE)

##
##		How is fit in the calibration data only?
##

calib.preds <- predData(test[with(test, date>=calib.start & date<=calib.end), ])
apply(calib.preds[3:10], 2, auc2, obs=calib.preds$obs)
apply(calib.preds[3:10], 2, maxF, obs=calib.preds$obs)


##------------------------------------------------------------------------------
##
##    	Figure 5: Thematic andensemble in-sample prediction correlations
##
##------------------------------------------------------------------------------

load("data/all_preds.rda")

cor.dat <- as.matrix(cor(all.preds[all.preds$date<=test.end, c(5:11, 4)]))
cor.dat <- round(cor.dat, digits=2)
cor.dat <- melt(cor.dat)

labels <- list("Leader char."="ch.m1", "Public discontent"="ch.m2", 
	"Global instability"="ch.m3", "Protest"="ch.m4", "Contagion"="ch.m5", 
	"Internal conflict"="ch.m6", "Financial"="ch.m7", "Ensemble"="ch.Ensemble")
levels(cor.dat$Var1) <- rev(labels)
levels(cor.dat$Var2) <- labels

# Source: http://www.peterhaschke.com/r/2013/04/23/CorrelationMatrix.html
p <- ggplot(cor.dat, aes(Var2, Var1, fill = value)) + 
	geom_tile() + 
	geom_text(aes(Var2, Var1, label = value), color = "#073642", size = 3) +
	scale_fill_gradient(name="Cor", low = "#fdf6e3", high = "steelblue",
    breaks=seq(0, 1, by = 0.2), limits = c(0, 1)) +
    labs(x = "", y = "") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(filename="graphics/pred_cor.png", plot=p, width=5, height=4, units="in",
	dpi=300)


##------------------------------------------------------------------------------
##
##		Figure 6: ROC curves
##
##------------------------------------------------------------------------------


rocData <- function(pred, obs, data=NULL) {
	if (!is.null(data)) {
		pred <- eval(substitute(pred), envir=data)
		obs  <- eval(substitute(obs), envir=data)
	}
	rocr.df <- prediction(pred, obs)
	rocr.pr <- performance(rocr.df, "tpr", "fpr")
	xy <- data.frame(recall=rocr.pr@x.values[[1]], precision=rocr.pr@y.values[[1]])
	return(xy)
}


xy <- data.frame(Partition="In sample", rocData(Ensemble, obs, in.sample.preds))
xy <- rbind(xy,
	data.frame(Partition="Test", rocData(Ensemble, obs, oo.sample.preds)))

p <- ggplot(data=xy, aes(x=recall, y=precision, col=Partition)) + 
	scale_x_continuous(expand = c(0.01, 0.01)) + 
	scale_y_continuous(expand = c(0.01, 0.01)) + 
	geom_line(show_guide=TRUE) +
	geom_abline(slope=1, color="gray", alpha=0.5) +
	labs(x="FPR", y="TPR") + 
	theme_bw() +
	theme(legend.position=c(0.75,0.25))

ggsave(filename="graphics/roc.png", plot=p, width=3.4, height=3, units="in",
	dpi=300)



##------------------------------------------------------------------------------
##
##		Figure 7: Precision-recall curves
##
##------------------------------------------------------------------------------


pr.data <- function(pred, obs, data=NULL) {
	if (!is.null(data)) {
		pred <- eval(substitute(pred), envir=data)
		obs  <- eval(substitute(obs), envir=data)
	}
	rocr.df <- prediction(pred, obs)
	rocr.pr <- performance(rocr.df, "prec", "rec")
	xy <- data.frame(recall=rocr.pr@x.values[[1]], precision=rocr.pr@y.values[[1]])
	return(xy)
}


xy <- data.frame(Partition="In sample", pr.data(Ensemble, obs, in.sample.preds))
xy <- rbind(xy,
	data.frame(Partition="Test", pr.data(Ensemble, obs, oo.sample.preds)))

pr <- ggplot(data=xy, aes(x=recall, y=precision, col=Partition)) + 
	geom_line(show_guide=TRUE) + labs(x="Recall", y="Precision") +
	scale_y_log10() + theme_bw() + theme(legend.position=c(0.75,0.75))

ggsave(filename="graphics/rec_prec_plot.png", plot=pr, width=3.4, height=3, units="in",
	dpi=300)

xy.themes <- lapply(in.sample.preds[, 4:10], pr.data, obs=in.sample.preds$obs)
xy.themes <- ldply(xy.themes, .id="model")

pr + geom_line(data=xy.themes, aes(group=model), colour="gray", alpha=0.5, show_guide=FALSE)



##------------------------------------------------------------------------------
##
##    	Table 2: Top 20 forecasts
##
##------------------------------------------------------------------------------


# How many months out to predict?
n.ahead <- 6

##		Create forecasts

# Ensemble forecasts for n.ahead months
fcast <- ensembleForecast(n.ahead, n.models=7)
fcast <- fcast[order(fcast[, 1], decreasing=TRUE), ]
head(fcast)

# Add cumulative total
total <- data.frame(
	gwcode=test[format(test$date, "%Y-%m")==format(test.end, "%Y-%m"), ]$ccode,
	total=(1 - apply(1 - ensembleForecast(n.ahead, n.models=7), 1, prod))
	)

# List of top 10 forecasts
top.list <- head(total[order(total$total, decreasing=TRUE), ], 10)
top.list <- data.frame(Country=rownames(top.list), Probability=top.list$total)
print(xtable(top.list, caption="Top 10 forecasts for IRC between April and September 2014 (6 months) using March 2014 data", 
	label="tab:forecast"), include.rownames=FALSE)


##------------------------------------------------------------------------------
##
##    	Figure 8: Map of forecasts
##
##------------------------------------------------------------------------------


# Function for plotting on world map
source("R/utilities/worldMap.R")

# Scale by 100 for proper coloring
total$scaled <- total$total*100

# Plot
dpi <- 300
png("graphics/map_6months.png", width=3*dpi, height=1.26*dpi, pointsize=20)
worldMap("scaled", "gwcode", total, legend.title="Percent")
dev.off()


##
##		Done
##