library(igraph) # R famous graph library
library(data.table) # A really fast data frame
library(ggplot2) # Nice plots in R
library(gridExtra) # To plot in a Grid
library(microbenchmark) #useful for profiling an ensamble of functions
library(cheddar) # R foodweb library
library(sqldf) # Data Fram SQL

# Loading our functions
source('network_reading.R')
source('network_aggregation.R')
source('network_sampling.R')
source('videos.R')
source('misc.R')

