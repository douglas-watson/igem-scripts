#!/usr/bin/env R -f
#
# nadine-plots.R - Nadine's plots
#
# Author:
# 	Douglas Watson <douglas@watsons.ch>
#
# Date:
# 	Continued on 20th September 2011
# 
# License:
# 	GNU GPL v3.0
#
##########################################

library('ggplot2')
library('reshape')
library('chron')

################################
# Nadine's experiments
################################

od = read.table('data/16.09-nadine-OD.dat', header=TRUE, sep='\t')
rfu = read.table('data/16.09-nadine-RFU.dat', header=TRUE, sep='\t')
samplenames = read.table('config/16.09-nadine.titles', header=TRUE, sep='\t')

# Melt and merge, convert time, insert sample concentration and name
od.m = melt(od, id=c('Time'), variable_name='Code')
rfu.m = melt(rfu, id=c('Time'), variable_name='Code')
both.m = data.frame(Time=od.m$Time, Code=od.m$Code, OD=od.m$value,
	RFU=rfu.m$value)
both.m$Time = as.POSIXct(sapply(both.m$Time, toString), format="%H:%M:%S")
both.m = merge(both.m, samplenames)

# Normalise RFU, calculate average over the last hour (keep repeats separate)
both.m$NormRFU = with(both.m, RFU/OD)

# calculate last hour averages, and repeat averages
averages = ddply(both.m, .variables=c('Time', 'SampleName', 'Concentration'),
	.fun=summarise, avg=mean(NormRFU), std=sd(NormRFU)/sqrt(2))
both.m$Hours = as.integer(format(both.m$Time, "%H"))
last_hour_averages = ddply(subset(both.m, Hours == max(Hours) - 1), 
		.variables=c('SampleName', 'Concentration', 'Repeat'), .fun=summarise,
		avg=mean(NormRFU), std=sd(NormRFU)/sqrt(6))

# Now separate experiments
# Current numbering: (output of levels(samplenames$SampleName))
# [1] "J61002 Plac-RFP + IPTG"      
# [2] "J61002 Ptet-RFP + ATC"        
# [3] "J6-pLac-RFP in BL21 + IPTG"      
# [4] "pSB3K1-pConst-TetR + J6-pTet-RFP + ATC"                      
# [5] "pSB3K1 Pconst-TetR Ptet-LacI + J61002 Plac-RFP + ATC"   
# [6] "pSB3K1 Pconst-TetR Ptet-LacI + J61002 Plac-RFP + IPTG"         

print(levels(samplenames$SampleName))
for ( i in seq(1, 6) ) {
	# Assign variables expN, expN.av, and expN.lh
	assign(paste('exp', i, sep=''), subset(both.m, 
			SampleName == levels(SampleName)[i]))
	assign(paste('exp', i, '.av', sep=''), subset(averages, 
			SampleName == levels(SampleName)[i]))
	assign(paste('exp', i, '.lh', sep=''), subset(last_hour_averages, 
			SampleName == levels(SampleName)[i]))
}

# Dynamics plots:
dynplot = function(all_data, av_data, colourlabel, clpalette='YlOrRd', 
	smooth=FALSE) {
	maintitle = toString(unique(all_data$SampleName))
	all_data$Concentration = factor(all_data$Concentration)
	av_data$Concentration = factor(av_data$Concentration)
	dyn = qplot(Time, NormRFU, data = all_data, geom = 'point',
		colour = Concentration,
		xlab = "Time", ylab = "Normalised RFU [a.u.]", main=maintitle) +
		#geom_smooth(method=loess, n=10) +
		scale_colour_brewer(colourlabel, palette = clpalette) + 
		scale_x_datetime(minor = '10 min', format = '%H h')
	if ( smooth == TRUE ) {
		dyn = dyn + geom_smooth()
	}
	return(dyn)
}

dyn_exp1 = dynplot(exp1, exp1.av, "Conc. [µM]", smooth=TRUE)
dyn_exp1

exp1$Concentration = factor(exp1$Concentration)
qplot(OD, RFU, data = exp1, colour= Concentration, log='y') +
	scale_colour_brewer(palette="YlOrRd")
qplot(Time, OD, data = exp1, colour=Concentration) +
	scale_colour_brewer(palette="YlOrRd")
max(exp1$OD)

# Dynamics and dose response plot for exp2:
dosplot = function(lh_data, xunit='[µM]') {
	maintitle = toString(unique(lh_data$SampleName))
	xbreaks = unique(lh_data$Concentration)[-2]
	p = qplot(Concentration, avg, data=lh_data, geom='point',
		colour=factor(Repeat),
		xlab = paste("Conc.", xunit), log = 'x', 
		ylab="Normalised RFU at saturation [a.u]", main=maintitle )  +
		geom_linerange(aes(ymin=avg-1.96*std, ymax=avg+1.96*std)) +
		scale_x_continuous(breaks = xbreaks) +
		scale_colour_hue("Repeat")
	return(p)
}

# Make both plots for all experiments
dyn_exp1 = dynplot(exp1, exp1.av, "Conc. [µM]", smooth=TRUE)
dyn_exp2 = dynplot(exp2, exp2.av, "Conc. [ng / ml]", clpalette="Spectral")
dyn_exp3 = dynplot(exp3, exp3.av, "Conc. [µM]")
dyn_exp4 = dynplot(exp4, exp4.av, "Conc. [ng / ml]", clpalette="Spectral",
	smooth=TRUE)
dyn_exp5 = dynplot(exp5, exp5.av, "Conc. [ng / ml]",
	smooth=TRUE)
dyn_exp6 = dynplot(exp6, exp6.av, "Conc. [µM]",
	smooth=TRUE)

dos_exp1 = dosplot(exp1.lh)
dos_exp2 = dosplot(exp2.lh, xunit='[ng/ml]')
dos_exp3 = dosplot(exp3.lh)
dos_exp4 = dosplot(exp4.lh, xunit='[ng/ml]')
dos_exp5 = dosplot(exp5.lh, xunit='[ng/ml]')
dos_exp6 = dosplot(exp6.lh)

# And save them
for ( i in seq(1, 6) ) { 
	ggsave(paste("plots/nadine-exp", i, "-dynamics.pdf", sep=''), 
		plot=get(paste('dyn_exp', i, sep='')), width=16, height=10)
	ggsave(paste("plots/nadine-exp", i, "-doseresponse.pdf", sep=''), 
		plot=get(paste('dos_exp', i, sep='')), width=16, height=10)
}

######################################
# TetR mutants
######################################

od = read.table('data/20.09-tetrmutants-OD.dat', header=TRUE, sep='\t')
rfu = read.table('data/20.09-tetrmutants-RFU.dat', header=TRUE, sep='\t')
samplenames = read.table('config/20.09-tetR.titles', header=TRUE, sep='\t')

# Melt and merge, convert time, insert sample concentration and name
od.m = melt(od, id=c('Time'), variable_name='Code')
rfu.m = melt(rfu, id=c('Time'), variable_name='Code')
both.m = data.frame(Time=od.m$Time, Code=od.m$Code, OD=od.m$value,
	RFU=rfu.m$value)
both.m$Time = as.POSIXct(sapply(both.m$Time, toString), format="%H:%M:%S")
both.m = merge(both.m, samplenames)

# Normalise RFU, 
both.m$NormRFU = with(both.m, RFU/OD)

# calculate last hour averages, and repeat averages
averages = ddply(both.m, .variables=c('Time', 'SampleName', 'Concentration'),
	.fun=summarise, avg=mean(NormRFU), std=sd(NormRFU)/sqrt(length(NormRFU)))
both.m$Hours = as.integer(format(both.m$Time, "%H"))
last_hour_averages = ddply(subset(both.m, Hours == max(Hours) - 1), 
		.variables=c('SampleName', 'Concentration', 'Repeat'), .fun=summarise,
		avg=mean(NormRFU), std=sd(NormRFU)/sqrt(length(NormRFU)))


# Separate by mutants: 
# [1] "E37AW43ST141A" "J36F43S"       "P39K"          "P39QY42M"     
# [5] "V36F"          "Y42FK108E"    
print(levels(samplenames$SampleName))
for ( i in seq(1, 6) ) {
	# Assign variables tetRN, tetRN.av, and tetRN.lh
	assign(paste('tetR', i, sep=''), subset(both.m, 
			SampleName == levels(SampleName)[i]))
	assign(paste('tetR', i, '.av', sep=''), subset(averages, 
			SampleName == levels(SampleName)[i]))
	assign(paste('tetR', i, '.lh', sep=''), subset(last_hour_averages, 
			SampleName == levels(SampleName)[i]))
}

factor

# Induction (or 'dynamics' plots)
dyn_tetR1 = dynplot(tetR1, tetR1.av, "Conc. [ng/ml]", clpalette="Spectral")
dyn_tetR2 = dynplot(tetR2, tetR2.av, "Conc. [ng/ml]", clpalette="Spectral")
dyn_tetR3 = dynplot(tetR3, tetR3.av, "Conc. [ng/ml]", clpalette="Spectral")
dyn_tetR4 = dynplot(tetR4, tetR4.av, "Conc. [ng/ml]", clpalette="Spectral")
dyn_tetR5 = dynplot(tetR5, tetR5.av, "Conc. [ng/ml]", clpalette="Spectral")
dyn_tetR6 = dynplot(tetR6, tetR6.av, "Conc. [ng/ml]", clpalette="Spectral")

dosplot2 = function(lh_data, xunit='[µM]') {
	maintitle = toString(unique(lh_data$SampleName))
	xbreaks = unique(lh_data$Concentration)[-2]
	p = qplot(Concentration, avg, data=lh_data, geom='line',
		xlab = paste("Conc.", xunit), log = 'x', 
		ylab="Normalised RFU at saturation [a.u]", main=maintitle,
		ylim = c(0, max(lh_data$avg)))  +
		geom_linerange(aes(ymin=avg-1.96*std, ymax=avg+1.96*std)) +
		scale_x_continuous(breaks = xbreaks) +
		scale_colour_hue("Repeat")
	return(p)
}

dos_tetR1 = dosplot2(tetR1.lh, xunit='[ng/ml]')
dos_tetR2 = dosplot2(tetR2.lh, xunit='[ng/ml]')
dos_tetR3 = dosplot2(tetR3.lh, xunit='[ng/ml]')
dos_tetR4 = dosplot2(tetR4.lh, xunit='[ng/ml]')
dos_tetR5 = dosplot2(tetR5.lh, xunit='[ng/ml]')
dos_tetR6 = dosplot2(tetR6.lh, xunit='[ng/ml]')

# Compare induction for each mutant at 0 ATC
ATC0 = subset(both.m, Concentration == '0')
p_atc0 = qplot(Time, NormRFU, data = ATC0, colour = SampleName,
		xlab = 'Time', ylab = 'Normalised fluorescence [a.u.]',
		main = 'TetR mutants, no ATC') +
	scale_x_datetime(minor = '10 min', format = '%H h') +
	scale_colour_brewer("Mutant", palette="Accent")

# Same but for 2000 ATC
ATC2000 = subset(both.m, Concentration == '2000')
p_atc2000 = qplot(Time, NormRFU, data = ATC2000, colour = SampleName,
		xlab = 'Time', ylab = 'Normalised fluorescence [a.u.]',
		main = 'TetR mutants, 2000 ng/ml ATC') +
	scale_x_datetime(minor = '10 min', format = '%H h') +
	scale_colour_brewer("Mutant", palette="Accent")

######################################
# PNG exports
######################################

# The first plots I made
for ( i in seq(1, 6) ) {
	ggsave(paste("plots/nadine-exp", i, "-induction.png", sep=''), 
		plot=get(paste('dyn_exp', i, sep='')),
		width=12, height=7.5)
	ggsave(paste("plots/nadine-exp", i, "-doseresponse.png", sep=''), 
		plot=get(paste('dos_exp', i, sep='')),
		width=12, height=7.5)
}

# TetR mutants

for ( i in seq(1, 6) ) {
 	ggsave(paste("plots/tetR-", get(paste("tetR", i, sep=''))$SampleName[1],
		"-induction.png", sep=""),
		plot=get(paste('dyn_tetR', i, sep='')),
		width=8, height=5)
 	ggsave(paste("plots/tetR-", get(paste("tetR", i, sep=''))$SampleName[1],
		"-doseresponse.png", sep=""),
		plot=get(paste('dos_tetR', i, sep='')),
		width=8, height=5)
}

ggsave("plots/tetR-ATC0-induction.png", plot=p_atc0, width=8, height=5)
ggsave("plots/tetR-ATC2000-induction.png", plot=p_atc2000, width=8, height=5)
