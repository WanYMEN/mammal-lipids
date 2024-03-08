rm(list=ls())

library(ggplot2)
library(gridExtra)
library(ggrepel)

setwd("D://work//skoltech//lipid//writing//GitHub//data")

info.milk <- read.csv("brain_milk_FA.info.milk.csv", row.names = 1)
info.PFC <- read.csv("brain_milk_FA.info.PFC.csv", row.names = 1)
info.CB <- read.csv("brain_milk_FA.info.CB.csv", row.names = 1)

milkSpecies <- c("HS", "HSm", "HSs", "MM", "MF", "CA", "SS")
brainSpecies <- c("HS", "MM", "CA", "SS")

data.milk <- read.csv("brain_milk_FA.normalized.milk.csv", row.names = 1)
data.milk[1:2, 1:5]
data.milk <- log2(data.milk[, rownames(info.milk)])

bothFA <- rownames(data.milk)

data.PFC <- read.csv("brain_milk_FA.normalized.PFC.csv", row.names = 1)
data.PFC <- log2(data.PFC[bothFA, rownames(info.PFC)])
data.CB <- read.csv("brain_milk_FA.normalized.CB.csv", row.names = 1)
data.CB <- log2(data.CB[bothFA, rownames(info.CB)])


## brain and milk intensity correlation across species
intensityDots <- function(data.brain, info.brain, region, corMethods)
{
all_cor <- c()
for(i in 1:length(brainSpecies)) {
    mean.brain <- rowMeans(data.brain[, rownames(info.brain)[info.brain$species==brainSpecies[i]]], na.rm = TRUE)
    milkCounterpart.1 <- milkSpecies[substr(milkSpecies, 1, 1)==substr(brainSpecies[i], 1, 1)]
        
    for(j in 1:length(milkCounterpart.1)) {
        if(milkCounterpart.1[j]=="HS") {
          mean.milk <- rowMeans(data.milk[, rownames(info.milk)[info.milk$species %in% c("HSm", "HSs")]], na.rm = TRUE)
        } else {
          mean.milk <- rowMeans(data.milk[, rownames(info.milk)[info.milk$species==milkCounterpart.1[j]]], na.rm = TRUE)
        }
        cr <- cor.test(mean.milk, mean.brain, method = corMethods)
        coef <- formatC(cr$estimate, digits = 4, format = "f")
        cr.pv <- formatC(cr$p.value, digits = 6, format = "f")
        all_cor.tmp <- c(region, paste("milk", milkCounterpart.1[j], sep="."), paste(region, brainSpecies[i], sep="."), cr$estimate, cr$p.value)
        all_cor <- rbind(all_cor, all_cor.tmp)
    } 
}
return(all_cor)
}

plotList.PFC.pes <- intensityDots(data.PFC, info.PFC, "PFC", "pearson")
plotList.CB.pes <- intensityDots(data.CB, info.CB, "CB", "pearson")

all_cor_value <- rbind(plotList.PFC.pes, plotList.CB.pes)
colnames(all_cor_value) <- c("Region", "Milk_species", "Brain_species", "Coefficient", "P.value")

write.table(all_cor_value, "brain_milk_FA.pearsonCorrelation.txt", row.names = FALSE, quote = FALSE, sep = "\t")


