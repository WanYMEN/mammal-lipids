## For normalizing the raw data matrix from CAMERA

rm(list=ls())

library(preprocessCore)
library(ggplot2)
library(reshape2)

setwd("D://work//skoltech//lipid//writing//GitHub//data")

## for milk ##

DATA <- read.csv("milk_FA.rawdata.csv", header= TRUE, row.names = 1)

data <- normalize.quantiles(as.matrix(DATA), copy=FALSE)

write.csv(data, "milk_FA.normalized.csv", quote=F, row.names = TRUE)

## for milk end##


## for brain ##

DATA <- read.csv("brain_FA.rawdata.csv", header= TRUE, row.names = 1)

data <- normalize.quantiles(as.matrix(DATA), copy=FALSE)

write.csv(data, "brain_FA.normalized.csv", quote=F, row.names = TRUE)

## for brain end##
