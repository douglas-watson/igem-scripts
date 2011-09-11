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

###############################
# Load data
###############################
fluo = read.table('data/05-09-cotransfs-fluo.dat', header=TRUE, sep='\t')
samplenames = read.table('config/samplenames.conf', header=TRUE, sep='\t')
titles = read.table('config/titles.conf', header=TRUE, sep='\t')

# Rename column names to sample names, convert times to time objects
names(fluo)[-1] = sapply(samplenames$SampleName, toString)
# Convert time to POSIX datetime objects. Introduces fake date so it works well.
fluo$Time = as.POSIXct(sapply(fluo$Time, toString), format="%H:%M:%S")

###############################
# Helper Functions
###############################

comparison_plot = function(data, melted_data) {
	# Make plot with title data$PlotTitle containing the fluorescence intensity
	# data for the corresponding data$SampleName, from the fluo.m data frame.
	# Values for each sample are plotted as a smooth curve, distinguished by
	# their colour.

	maintitle = data$PlotTitle[1]
	subdata = subset(melted_data, SampleName %in% data$SampleName)

	# For proper colouring, the SampleName factor must be redefined:
	subdata$SampleName = factor(subdata$SampleName)
	p = qplot(Time, value, data=subdata, colour=SampleName, geom='smooth',
		main=maintitle, ylab='Fluorescence') +	
		scale_x_datetime(major='2 hours', format='%H h') +
		scale_colour_hue('Sample Name', breaks=data$SampleName, 
			labels=sub('-', ' ', data$SampleName))
	ggsave(paste('plots/', maintitle, '.pdf', sep=''), width=16, height=10)
	return(p)
}

###############################
# Plotting
###############################

fluo.m = melt(fluo, id.vars='Time', variable_name='SampleName')

# Plot fluorescence versus time, facetted by sample type
p = qplot(Time, value, data = fluo.m) + 
	facet_wrap(~ SampleName, ncol=12) +
	scale_x_datetime(major='4 hours', format="%Hh")


# Comparison plots, as specified in config/titles.conf
d_ply(titles, .variables=('PlotTitle'), .fun=comparison_plot, melted_data=fluo.m, .print=FALSE)

