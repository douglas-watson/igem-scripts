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
args = commandArgs(trailingOnly = TRUE)

fluo = read.table('data/05-09-cotransfs-fluo.dat', header=TRUE, sep='\t')
od = read.table('data/05-09-cotransfs-OD.dat', header=TRUE, sep='\t')
samplenames = read.table('config/samplenames.conf', header=TRUE, sep='\t')

# Rename column names to sample names, convert times to time objects
names(fluo)[-1] = sapply(samplenames$SampleName, toString)
names(od)[-1] = sapply(samplenames$SampleName, toString)
# Convert time to POSIX datetime objects. Introduces fake date so it works well.
fluo$Time = as.POSIXct(sapply(fluo$Time, toString), format="%H:%M:%S")
od$Time = as.POSIXct(sapply(od$Time, toString), format="%H:%M:%S")

# Add an hours column
fluo$Hours = as.integer(format(fluo$Time, "%H"))
od$Hours = as.integer(format(od$Time, "%H"))

###############################
# Helper Functions
###############################

colour_comparison_plot = function(molten_data, samples, maintitle,
	outfile=NULL, ylabel="Fluorescence") {
	# Plot in plots/`outfile`.pdf, on a single figure, the values in
	# `molten_data`, only for the samples in `samples`. Samples are
	# distinguished by their colour.
	#
	# If unspecified, `outfile` defaults to the main title.

	if ( is.null(outfile) )
		outfile = maintitle

	subdata = subset(molten_data, SampleName %in% samples)
	# For proper colouring, the SampleName factor must be redefined:
	subdata$SampleName = factor(subdata$SampleName)

	p = qplot(Time, avg, data=subdata, colour=SampleName, 
		main=maintitle, ylab=ylabel) +
		geom_errorbar(aes(ymin=avg - std, ymax=avg + std)) +
		scale_x_datetime(major='2 hours', minor='10 min', format='%H h',
			limits=c(min(subdata$Time), max(subdata$Time))) +
		scale_colour_hue('Sample Name')
	ggsave(paste('plots/', outfile, '.pdf', sep=''), width=16, height=10)

	return(p)
}

###############################
# Calculations and Plotting
###############################

fluo.m = melt(fluo, id.vars=c('Time', 'Hours'), 
	variable_name='SampleName')
od.m = melt(od, id.vars=c('Time', 'Hours'), 
	variable_name='SampleName')

combined.m = data.frame(Time=fluo.m$Time, Hours=fluo.m$Hours,
	SampleName=fluo.m$SampleName, fluo=fluo.m$value, od=od.m$value)


# Plot raw fluorescence versus time, facetted by sample type
p = qplot(Time, value, data = fluo.m) + 
	facet_wrap(~ SampleName, ncol=12) +
	scale_x_datetime(major='4 hours', format="%Hh")
# ggsave("plots/facet_plot_F.pdf", width=24, height=15)

# And the same for optical density
n = qplot(Time, value, data = od, ylab = "Optical Density") +
	facet_wrap(~ SampleName, ncol=12)
# ggsave("plots/facet_plot_OD.pdf", width=24, height=14)

# Calculate average fluorescence per hours, normalised by OD
normalised = ddply(combined.m, .variables=c('SampleName', 'Hours'), 
	.fun=summarise, avg=mean(fluo)/mean(exp(od)), std=sd(fluo/exp(od)))
normalised$Time = as.POSIXct(sapply(normalised$Hours, toString), format="%H")

colour_comparison_plot(normalised,
	maintitle = "Some Plot",
	ylabel = "Normalised Fluorescence",
	samples = c('HUGH-Port',
				'CLINTON-Rizzuto',
				'CHRISTIAN-Boose',
				'JULIO-Worster',
				'ERIK-Klippel',
				'TYRONE-Mccrum',
				'CLAYTON-Sheperd'),
	)


