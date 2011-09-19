#!/usr/bin/env R -f
#
# plate-reader.R - Analyse data from plate reader experiments, for optical
# 		density and fluorescence
#
# Author:
# 	Douglas Watson <douglas@watsons.ch>
#
# Date:
# 	Started on 9th September 2011
# 
# License:
# 	GNU GPL v3.0
#
##########################################

library('ggplot2')
library('reshape')
library('chron')

###############################
# Load data
###############################
od = read.table('data/11.09-lysis-final-OD.dat', header=TRUE, sep='\t')
samplenames = read.table('config/11.09-lysis-final.titles', header=TRUE, 
	sep='\t')

# Convert time to POSIX datetime objects. Introduces fake date so it works well.
od$Time = as.POSIXct(sapply(od$Time, toString), format="%H:%M:%S")

od.m = melt(od, id.vars=c('Time'), variable_name='Code')

# Insert treatment names, concentration, and repeat number
merged = merge(od.m, samplenames)

# Keep only 1:10 samples (best results with these and kick out the outlier
merged = merged[grep('1:10', merged$SampleName),]
merged = merged[(merged$Code != "A9"),]

# Average repeats of a same treatment
averaged = ddply(merged, .variables=c('Time', 'SampleName', 'Concentration'), 
	.fun=summarise, avg=mean(value), std=sd(value)/sqrt(3))

################################
# Lysis dynamics
################################

# We are only interested in the T7-C2 (1:10) sample (the one that lyses best)
only_tyc2 = merged[merged$SampleName == 'T7-C2 (1:10)',]
only_tyc2_av = averaged[averaged$SampleName == 'T7-C2 (1:10)',]
only_tyc2$Concentration = factor(only_tyc2$Concentration)
only_tyc2_av$Concentration = factor(only_tyc2_av$Concentration)


# Plot dynamics: points for the separate treatments, and a line representing
# the average of each treatment.
p_lysis_dynamics = qplot(Time, avg, data = only_tyc2_av,
	colour=Concentration, geom='line', 
	xlab="Time", ylab="Optical Density [a.u.]",
	main="T7-C2 (1:10) - Lysis dynamics") +
	scale_x_datetime(major='1 hours', minor='10 min', format="%H h") +
	scale_colour_discrete() +
	# geom_errorbar(aes(ymin=avg-1.96*std, ymax=avg+1.96*std)) +
	geom_point(aes(y=value), data=only_tyc2) +
	scale_colour_brewer("IPTG Conc. [µM]", palette="YlOrRd")
p_lysis_dynamics

ggsave("plots/lysis_dynamics.pdf", width=16, height=10)


################################
# Dose response
################################


# Average over the six measurements in the last hour. Keep repeats separate.
merged$Hours = as.integer(format(merged$Time, "%H"))
last_hour = merged[merged$Hours == max(merged$Hours),]
dose_response = ddply(last_hour, .variables=c('SampleName', 'Concentration', 
		'Repeat'), .fun=summarise, hourly_avg=mean(value), std=sd(value)/sqrt(6))

# Average the repeats
dose_response_avg = ddply(last_hour, .variables=c('SampleName', 
		'Concentration'), 
	.fun=summarise, avg=mean(value), std=sd(value)/sqrt(6*3))

# Again, points for each repeat, then a line as general average
dose_response$SampleName = factor(dose_response$SampleName)
dose_response_avg$SampleName = factor(dose_response_avg$SampleName)
p_dose_response = qplot(Concentration, hourly_avg, data=dose_response, 
	colour=SampleName, geom='point', xlab="IPTG concentration [µM]", 
	ylab="Optical Density [a.u.]") +
	scale_colour_hue("Promoter type") +
	geom_line(aes(y=avg), data=dose_response_avg)
ggsave("plots/dose-response.pdf", width=16, height=10)

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
dynplot = function(all_data, av_data, colourlabel, clpalette='YlOrRd') {
	maintitle = toString(unique(all_data$SampleName))
	all_data$Concentration = factor(all_data$Concentration)
	av_data$Concentration = factor(av_data$Concentration)
	dyn = qplot(Time, NormRFU, data = all_data, geom = 'point',
		colour = Concentration,
		xlab = "Time", ylab = "Normalised RFU [a.u.]", main=maintitle) +
		geom_smooth() +
		scale_colour_brewer(colourlabel, palette = clpalette) + 
		scale_x_datetime(minor = '10 min', format = '%H h')
	return(dyn)
}

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
dyn_exp1 = dynplot(exp1, exp1.av, "Conc. [µM]")
dyn_exp2 = dynplot(exp2, exp3.av, "Conc. [ng / ml]", clpalette="Spectral")
dyn_exp3 = dynplot(exp3, exp2.av, "Conc. [µM]")
dyn_exp4 = dynplot(exp4, exp3.av, "Conc. [ng / ml]", clpalette="Spectral")
dyn_exp5 = dynplot(exp5, exp4.av, "Conc. [ng / ml]")
dyn_exp6 = dynplot(exp6, exp5.av, "Conc. [µM]")

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

###############################
# Non-random variants
###############################

od = read.table('data/10.09-rfp-nonrand-OD.dat', header=TRUE, sep='\t')
rfu = read.table('data/10.09-rfp-nonrand-RFU.dat', header=TRUE, sep='\t')
samplenames = read.table('config/10.09-rfp-nonrand.titles', header=TRUE, 
	sep='\t')

# Transform time to a time object
rfu.m = melt(rfu, id.vars=c('Time'), variable_name='Code')
od.m = melt(od, id.vars=c('Time'), variable_name='Code')
normalised.m = data.frame(
	Time=as.POSIXct(sapply(rfu$Time, toString), format="%H:%M:%S"),
	Code=rfu.m$Code, norm_rfu=rfu.m$value/od.m$value)

# Merge with sample names (by code), to get name, strength, etc.
normalised.m = merge(normalised.m, samplenames)

# Extract only the last hour (saturation), and average over those points
normalised.m$Hours = as.integer(format(normalised.m$Time, "%H"))
last_hour = normalised.m[normalised.m$Hours == max(normalised.m$Hours),]
last_hour$Concentration = factor(last_hour$Concentration)
last_hour_avg = ddply(last_hour, 
	.variables=c('SampleName', 'Variant', 'Strength', 'Concentration'), 
	.fun=summarise, avg=mean(norm_rfu), std=sd(norm_rfu)/sqrt(18))

# Order the labels by affinity (RFU at 250 uM IPTG) [NO LONGER USED]
# last_250 = last_hour_avg[last_hour_avg$Concentration == 250, ]
# ordered_breaks <- last_250$SampleName[order(-last_250$avg)]
# last_hour_avg$SampleName = ordered(last_hour_avg$SampleName,
	# levels=ordered_breaks)


# Bar chart, ordered by expected affinity. Separated by concentration and lac
# variant
w = 10
dodge <- position_dodge(width=w)
last_hour_avg$Strength = as.integer(last_hour_avg$Strength)
p_non_random = qplot(Strength, avg, data=last_hour_avg, 
	fill=Concentration, facets=Variant ~ .,
	geom='bar', position=dodge, stat='identity', width=w,
	xlab="Expected promoter efficiency", ylab="Normalised RFU at Saturation [a.u]",
	main="Designed T7 Variants - Dose Response") + 
	geom_errorbar(aes(ymin=avg-1.96*std, ymax=avg+1.96*std),
		position=dodge, width=w/5) +
	scale_fill_hue("IPTG conc. [µM]") +
	scale_x_continuous(breaks=c(14, 30, 54, 80, 100, 110))
ggsave("plots/non-random.pdf", width=16, height=10)

## Compare expression ratio (with and without IPTG)
last_hour_clean = last_hour[c('Variant', 'Strength', 'Concentration', 
		'Repeat', 'norm_rfu')] 
last_hour.m = melt(last_hour_clean, 
	id=c('Variant', 'Strength', 'Concentration', 'Repeat'))
last_hour.m$Concentration = paste("IPTG", last_hour.m$Concentration, 
	sep='')
last_hour.c = cast(last_hour.m, 
	Variant + Strength + Repeat ~ Concentration)
last_hour.c$InductionRatio = with(last_hour.c, IPTG250/IPTG0)

ratio_averages = ddply(last_hour.c, .variables=c('Variant', 'Strength'),
	.fun=summarise, avg=mean(InductionRatio), std=sd(InductionRatio)/sqrt(3))

ratio_averages$Variant = ordered(ratio_averages$Variant,
	levels=c("Without lac operator", "With lac operator"))
p_induction_ratio = qplot(Strength, avg, data=ratio_averages, fill=Variant,
	geom='bar', stat='identity', position=dodge, width=w,
	xlab="Expected promoter efficiency [% of WT]", 
	ylab="Induction Ratio (250 µM IPTG / 0 µM IPTG)",
	main="Induction ratio comparison for designed T7 variants") +
	geom_errorbar(aes(ymin=avg-1.96*std, ymax=avg+1.96*std), 
		position=dodge, width=w/5) +
	scale_x_continuous(breaks=ratio_averages$Strength) +
	scale_fill_hue("Promoter type")
ggsave("plots/induction_ratio.pdf", width=16, height=10)

###############################
# Random variants
###############################

rand_od = read.table('data/11.09-rfp-rand-OD.dat', header=TRUE, sep='\t')
rand_rfu = read.table('data/11.09-rfp-rand-RFU.dat', header=TRUE, sep='\t')
# rand_names = read.table('config/11.09-rfp-rand.titles', header=TRUE, 
	# sep='\t')

rand_rfu.m = melt(rand_rfu, id.vars=c('Time'), variable_name='Code')
rand_od.m = melt(rand_od, id.vars=c('Time'), variable_name='Code')
rand_norm.m = data.frame(
	Time=as.POSIXct(sapply(rand_rfu$Time, toString), format="%H:%M:%S"),
	Code=rand_rfu.m$Code, norm_rfu=rand_rfu.m$value/rand_od.m$value)

# Compare designed and random mutants on a boxplot
normalised.m$type = "Designed"
rand_norm.m$type = "Random"
both = merge(normalised.m, rand_norm.m, all=T)

p_var_comparison = qplot(type, norm_rfu, data = both, geom = 'boxplot', 
	binwidth=50, position='dodge', ylab="Normalised RFU [a.u.]", xlab=NULL, 
	  main="Variability in promoter efficiency")
ggsave("plots/variability_comparison.pdf", width=3, height=7)

######################################
# DNA recovery experiment
######################################

rfp = read.table('data/17.09-dnarecov-RFP.dat', header=TRUE, sep='\t')
gfp = read.table('data/17.09-dnarecov-GFP.dat', header=TRUE, sep='\t')

rfp.m = melt(rfp, id=c('Hour'), variable_name='Sample')
gfp.m = melt(gfp, id=c('Hour'), variable_name='Sample')

# Rename 'value' column to RFU or GFU
names(rfp.m)[3] = 'RFU'
names(gfp.m)[3] = 'GFU'
 
# merge and melt again, clean up sample name
dna_recovery = merge(rfp.m, gfp.m)
dna_recovery$Sample = factor(dna_recovery$Sample, 
	labels=c('C2-RFP & C11-GFP', 'C2-GFP & C11-RFP'))
dna_recovery.m = melt(dna_recovery, id=c('Hour', 'Sample'), 
	variable_name='Measure')

# Plot comparison of two output colours, for both samples

qplot(Hour, value, data = dna_recovery.m, geom = 'point',
	colour = Sample,
	xlab = 'Time [hours]', ylab = 'Fluorescence [a.u.]') +
	facet_grid(Measure ~ ., scales = 'free_y')

######################################
# PNG exports
######################################

# Reminder: here are the plots I have
# p_dose_response	
# p_induction_ratio	
# p_lysis_dynamics	
# p_non_random	
# p_nonrand	
# p_nonrand_big	
# p_nonrand_sub	
# p_rand	
# p_rand_big	
# p_rand_sub	
# p_var_comparison	

ggsave("plots/dose_response.png", plot=p_dose_response, 
	width=12, height=7.5)
ggsave("plots/induction_ratio.png", plot=p_induction_ratio, 
	width=12, height=7.5)
ggsave("plots/lysis_dynamics.png", plot=p_lysis_dynamics,
	width=12, height=7.5)
ggsave("plots/non_random_dose_response.png", plot=p_non_random,
	width=12, height=7.5)
ggsave("plots/varability_comparison.png", plot=p_var_comparison, 
	width=4, height=10)

# Nadine's stuff

for ( i in seq(1, 6) ) {
	ggsave(paste("plots/nadine-exp", i, "-induction.png", sep=''), 
		plot=get(paste('dyn_exp', i, sep='')),
		width=12, height=7.5)
	ggsave(paste("plots/nadine-exp", i, "-doseresponse.png", sep=''), 
		plot=get(paste('dos_exp', i, sep='')),
		width=12, height=7.5)
}
