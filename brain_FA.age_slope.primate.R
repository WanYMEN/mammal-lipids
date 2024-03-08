rm(list=ls())

library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)


setwd("D://work//skoltech//lipid//writing//GitHub//data")

data <- read.csv("brain_FA.normalized.csv", row.names = 1)
data <- log2(data)

info <- read.csv("brain_FA.info.csv", row.names = 1)

allSpecies <- c("HS","PT", "MM")
allRegion <- c("PFC", "CB")

info <- info[info$species %in% allSpecies, ]
info$species <- factor(info$species, levels=allSpecies, order = TRUE)

info$age <- as.numeric(as.character(info$age))

info$scaledDay <- info$age
info[info$species=="PT", ]$scaledDay <- (info[info$species=="PT", ]$scaledDay)*2
info[info$species=="MM", ]$scaledDay <- (info[info$species=="MM", ]$scaledDay)*2.5
info$logScaledDay <- log2(info$scaledDay + 1)

colors <- c("HS" = "#FF0000", "PT"= "#A52A2A", "MM" = "#FF8C00")

info <- info[info$age <1000, ]

info.PFC <- info[info$region=="PFC", ]
data.PFC <- data[, rownames(info.PFC)]

info.CB <- info[info$region=="CB", ]
data.CB <- data[, rownames(info.CB)]


## lines plots and slopes for intensity across age 
ageSlopes <- function(ageData, ageInfo, allSpecies, type, Xlab)
{
m=0
slopes <- c()
Pvalues <- c()
slopes.2_5 <- c()
slopes.97_5 <- c()
for(i in 1:nrow(ageData)) {
    m <- m+1
    ints.tmp <- as.numeric(as.character(ageData[i, ]))
    dt.tmp <- data.frame(cbind(ints.tmp, ageInfo[, c(type, "species")]))
    colnames(dt.tmp) <- c("intensity", "age", "species")

    slope.sub <- c()
    pvalue.sub <- c()
    slopes.2_5.sub <- c()
    slopes.97_5.sub <- c()
    for(j in 1:length(allSpecies)) {
       dt.sub <- dt.tmp[dt.tmp$species==allSpecies[j], ]
       lm.tmp <- lm(dt.sub$intensity ~ dt.sub$age)
       slope.tmp <- as.numeric(coef(lm.tmp)[2])
       slope.sub <- cbind(slope.sub, slope.tmp)
       slopes.2_5.tmp <- confint(lm.tmp, 'dt.sub$age', level = 0.95)[1, 1]
       slopes.97_5.tmp <- confint(lm.tmp, 'dt.sub$age', level = 0.95)[1, 2]
       slopes.2_5.sub <- cbind(slopes.2_5.sub, slopes.2_5.tmp)
       slopes.97_5.sub <- cbind(slopes.97_5.sub, slopes.97_5.tmp)
       
       sm <- summary(lm.tmp)
       pvalue.sub <- cbind(pvalue.sub, formatC(sm$coefficients[2, "Pr(>|t|)"], digits = 6, format = "f"))
    }
    slopes <- rbind(slopes, slope.sub)
    Pvalues <- rbind(Pvalues, pvalue.sub)
    slopes.2_5 <- rbind(slopes.2_5, slopes.2_5.sub)
    slopes.97_5 <- rbind(slopes.97_5, slopes.97_5.sub)
}

colnames(slopes) <- allSpecies
rownames(slopes) <- rownames(ageData)

colnames(Pvalues) <- allSpecies
rownames(Pvalues) <- rownames(ageData)

colnames(slopes.2_5) <- paste(allSpecies, "2_5", sep=".")
rownames(slopes.2_5) <- rownames(ageData)

colnames(slopes.97_5) <- paste(allSpecies, "97_5", sep=".")
rownames(slopes.97_5) <- rownames(ageData)

return(list(slopes, Pvalues, slopes.2_5, slopes.97_5))

}


logScaledDay.PFC <- ageSlopes(data.PFC, info.PFC, allSpecies, "logScaledDay", "log2(scaledDays+1)")
slope_CI.PFC <- cbind(logScaledDay.PFC[[1]], logScaledDay.PFC[[3]], logScaledDay.PFC[[4]])
write.table(slope_CI.PFC, "brain_FA.age_slope.PFC.txt", row.names = TRUE, quote = FALSE, sep = "\t")

logScaledDay.CB <- ageSlopes(data.CB, info.CB, allSpecies, "logScaledDay", "log2(scaledDays+1)")
slope_CI.CB <- cbind(logScaledDay.CB[[1]], logScaledDay.CB[[3]], logScaledDay.CB[[4]])
write.table(slope_CI.CB, "brain_FA.age_slope.CB.txt", row.names = TRUE, quote = FALSE, sep = "\t")

