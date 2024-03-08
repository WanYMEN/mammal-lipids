rm(list=ls())

library(dplyr)
library(stringr)
library(car)

setwd("D://work//skoltech//lipid//writing//GitHub//data")

data <- read.csv("brain_FA.normalized.csv", row.names = 1)

info <- read.csv("brain_FA.info.csv", row.names = 1)
info <- info[info$species != "QC", ] 

level.s <- c("HS", "PT", "MM", "CA", "SS")
level.g <- c("PR", "CA", "SS")
level.p <- c("HS", "PT", "MM")

info$species <- factor(info$species, levels = level.s, order = TRUE)
info$group <- factor(info$group, levels = level.g, order = TRUE)

str(info)

info$scaledDay <- info$age
info[info$species=="PT", ]$scaledDay <- (info[info$species=="PT", ]$scaledDay)*2
info[info$species=="MM", ]$scaledDay <- (info[info$species=="MM", ]$scaledDay)*2.5
info[info$species=="CA", ]$scaledDay <- (info[info$species=="CA", ]$scaledDay)*5
info[info$species=="SS", ]$scaledDay <- (info[info$species=="SS", ]$scaledDay)*5

info$logScaledDay <- log2(info$scaledDay + 1)

info.PFC <- info[info$region=="PFC", ]
data.PFC <- data[, rownames(info.PFC)]

info.CB <- info[info$region=="CB", ]
data.CB <- data[, rownames(info.CB)]


adjustedIntensity <- function(dat, inf, item, Flevel)
{
cf <- dat

for(i in 1:nrow(dat)){
    dt.tmp <- data.frame(cbind(as.numeric(as.character(dat[i, ])), as.character(inf[, "sex"]), inf[, "logScaledDay"], as.character(inf[, item])))
    colnames(dt.tmp) <- c("value", "sex", "age", "species")
    dt.tmp$value <- log2(as.numeric(as.character(dt.tmp$value)))
    dt.tmp$age <- as.numeric(as.character(dt.tmp$age))
    dt.tmp$sex <- as.numeric(as.factor(dt.tmp$sex)) - 1
    dt.tmp$species <- factor(dt.tmp$species, levels = Flevel, order = TRUE)
    dt.tmp$species <- as.numeric(dt.tmp$species) - 1

    lmmod1 <- lm(value ~ sex + age + species, data = dt.tmp)
    Coef <- coef(lmmod1)
    cf[i, ] <- dt.tmp$value - Coef["sex"]*dt.tmp$sex - Coef["age"]*dt.tmp$age
    
}
return(cf)
}


##PFC
## species
adj.s.PFC <- data.frame(adjustedIntensity(data.PFC, info.PFC, "species", level.s))
adj.s.PFC <- 2^adj.s.PFC
write.csv(adj.s.PFC, "brain_FA.sexAgeAdjust.species.PFC.csv", row.names = TRUE, quote = FALSE)

## primates
subInfo.PFC <- info.PFC[info.PFC$group == "PR", ]
subData.PFC <- data.PFC[, rownames(subInfo.PFC)]

adj.p.PFC <- data.frame(adjustedIntensity(subData.PFC, subInfo.PFC, "species", level.p))
adj.p.PFC <- 2^adj.p.PFC
write.csv(adj.p.PFC, "brain_FA.sexAgeAdjust.primate.PFC.csv", row.names = TRUE, quote = FALSE)


##CB
## species
adj.s.CB <- data.frame(adjustedIntensity(data.CB, info.CB, "species", level.s))
adj.s.CB <- 2^adj.s.CB
write.csv(adj.s.CB, "brain_FA.sexAgeAdjust.species.CB.csv", row.names = TRUE, quote = FALSE)

## primates
subInfo.CB <- info.CB[info.CB$group == "PR", ]
subData.CB <- data.CB[, rownames(subInfo.CB)]

adj.p.CB <- data.frame(adjustedIntensity(subData.CB, subInfo.CB, "species", level.p))
adj.p.CB <- 2^adj.p.CB
write.csv(adj.p.CB, "brain_FA.sexAgeAdjust.primate.CB.csv", row.names = TRUE, quote = FALSE)


