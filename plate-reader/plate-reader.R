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
	ylab="Normalised RFU at saturation [a.u.]") +
	scale_colour_hue("Promoter type") +
	geom_line(aes(y=avg), data=dose_response_avg)
ggsave("plots/dose-response.pdf", width=16, height=10)


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
	main="Designed Mutants - Dose Response") + 
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
	main="Induction ratio comparison for designed mutants") +
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
