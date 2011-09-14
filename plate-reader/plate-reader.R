#!/usr/bin/env R
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
qplot(Time, avg, data = only_tyc2_av,
	colour=Concentration, geom='line', 
	xlab="Time [hour:min]", ylab="Optical Density [a.u.]",
	main="T7-C2 (1:10) - Lysis dynamics") +
	scale_x_datetime(major='1 hours', minor='10 min') +
	scale_colour_discrete() +
	# geom_errorbar(aes(ymin=avg-1.96*std, ymax=avg+1.96*std)) +
	geom_point(aes(y=value), data=only_tyc2) +
	scale_colour_hue("Concentration [µM]")

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
qplot(Concentration, hourly_avg, data=dose_response, colour=SampleName, 
	geom='point', xlab="Concentration [µM]", 
	ylab="Last hour's averaged RFU [a.u.]") +
	geom_line(aes(y=avg), data=dose_response_avg)
ggsave("plots/dose-response.pdf", width=16, height=10)


###############################
# Non-random variants
###############################

od = read.table('data/10.09-rfp-nonrand-OD.dat', header=TRUE, sep='\t')
rfu = read.table('data/10.09-rfp-nonrand-RFU.dat', header=TRUE, sep='\t')
samplenames = read.table('config/10.09-rfp-nonrand.titles', header=TRUE, 
	sep='\t')

rfu.m = melt(rfu, id.vars=c('Time'), variable_name='Code')
od.m = melt(od, id.vars=c('Time'), variable_name='Code')
normalised.m = data.frame(
	Time=as.POSIXct(sapply(rfu$Time, toString), format="%H:%M:%S"),
	Code=rfu.m$Code, norm_rfu=rfu.m$value/od.m$value)

normalised.m = merge(normalised.m, samplenames)

normalised.m$Hours = as.integer(format(normalised.m$Time, "%H"))
last_hour = normalised.m[normalised.m$Hours == max(normalised.m$Hours),]

last_hour$Concentration = factor(last_hour$Concentration)
last_hour_avg = ddply(last_hour, 
	.variables=c('SampleName', 'Concentration'), 
	.fun=summarise, avg=mean(norm_rfu), std=sd(norm_rfu)/sqrt(18))

dodge <- position_dodge(width=0.9)
last_250 = last_hour_avg[last_hour_avg$Concentration == 250, ]
ordered_breaks <- last_250$SampleName[order(-last_250$avg)]

last_hour_avg$SampleName = ordered(last_hour_avg$SampleName,
	levels=ordered_breaks)

qplot(SampleName, avg, data=last_hour_avg, geom='bar',
	fill=Concentration, position='dodge',
	xlab="Sample", ylab="Normalised fluorescence [a.u]") +
	geom_linerange(aes(ymin=avg-1.96*std, ymax=avg+1.96*std),
			position=dodge) +
	scale_fill_hue("Concentration [µM]")
ggsave("plots/non-random.pdf", width=16, height=10)

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

qplot(norm_rfu, data=rand_norm.m, geom='histogram',
	binwidth=10)
ggsave("plots/rand_variability.pdf", width=16, height=10)

qplot(norm_rfu, data=normalised.m, geom='histogram',
	binwidth=10)
ggsave("plots/nonrand_variability.pdf", width=16, height=10)
