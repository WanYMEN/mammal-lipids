rm(list=ls())

library(ggplot2)
library(reshape2)
library(stringr)
library(gridExtra)
library(ggrepel)
library(plyr)

setwd("D://work//skoltech//lipid//writing//GitHub//data")

DATA.milk <- read.csv("milk_FA.normalized.csv", row.names = 1)
info.milk <- read.csv("milk_FA.info.csv", row.names = 1)

info.shg <- read.csv("milk_FA.info.humanShangHai.csv", row.names = 1)
info.msk <- read.csv("milk_FA.info.humanMoscow.csv", row.names = 1)

info.msk <- info.msk[info.msk$LactationStage<500, ]

info.milk.hm <- rbind(info.shg[, c("LactationStage", "Parity", "N_C", "Sex")], info.msk[, c("LactationStage", "Parity", "N_C", "Sex")])
info.milk.hm$LS <- as.numeric(formatC(info.milk.hm$LactationStage/7, digits = 0, format = "d")) + 1
info.milk.hm[info.milk.hm$LS >4, ]$LS <- ">4"

LSgroup <- c("1", "2", "3", "4", ">4")

info.milk.hm$LS <- factor(info.milk.hm$LS, levels = LSgroup, order =TRUE)

group.s <-  c("MM", "CA", "SS")

FC.PFC <- read.csv("brain_milk_FA.species_FC.PFC.csv", row.names = 1)
FC.CB <- read.csv("brain_milk_FA.species_FC.CB.csv", row.names = 1)

data.milk <- DATA.milk[rownames(FC.PFC), ]

LSgroup.corTest <- function(FC.brain, dt.milk.ot, corMethods, groupSpecies){
LScorTest <- c()
for(i in 1:length(LSgroup)){
    nm.tmp <- rownames(info.milk.hm)[info.milk.hm$LS==LSgroup[i]]
    dt.tmp <- data.milk[, nm.tmp]

    LS.FC <- c()
    LSmean_tmp <- c()
    for(j in 1:nrow(data.milk)){
        mean.dt <- mean(as.numeric(dt.tmp[j, ]), na.rm = TRUE)
        mean.milk <- mean(as.numeric(dt.milk.ot[j, ]), na.rm = TRUE)
        FC.tmp <- log2(mean.dt/mean.milk)
        LS.FC <- c(LS.FC, FC.tmp)
    }

   cr <- cor.test(LS.FC, FC.brain, method = corMethods)
   coef <- formatC(cr$estimate, digits = 4, format = "f")
   cr.pv <- formatC(cr$p.value, digits = 6, format = "f")
   sub <- cbind(coef, cr.pv)
   LScorTest <- rbind(LScorTest, sub)
}
return(LScorTest)
}


for(g in 1:length(group.s)){
    inf.milk.ot <- info.milk[info.milk$species==group.s[g], ]
    FC.PFC.ot <- FC.PFC[, grepl(paste("PFC.HS", group.s[g], sep="."), colnames(FC.PFC))]
    FC.CB.ot <- FC.CB[rownames(FC.PFC), grepl(paste("CB.HS", group.s[g], sep="."), colnames(FC.CB))]
    data.milk.ot <- data.milk[, rownames(inf.milk.ot)]

    PFC_cor <- LSgroup.corTest(FC.PFC.ot, data.milk.ot, "pearson", group.s[g])
    CB_cor <- LSgroup.corTest(FC.CB.ot, data.milk.ot, "pearson", group.s[g])
    
    corTest.FC.pes <- data.frame(cbind(PFC_cor, CB_cor), stringsAsFactors = FALSE)
    colnames(corTest.FC.pes) <- c("coef.PFC", "pvalue.PFC", "coef.CB", "pvalue.CB")
    rownames(corTest.FC.pes) <- LSgroup

    write.csv(corTest.FC.pes, paste("brain_milk_FA.species_FC_LS_pearson.HS_", group.s[g], ".csv", sep=""), row.names = TRUE, quote = FALSE)
}
