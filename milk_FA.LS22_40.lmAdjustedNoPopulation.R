rm(list=ls())

library(ggplot2)
library(reshape2)
library(car)
library(gridExtra)

setwd("D://work//skoltech//lipid//writing//GitHub//data")

data <- read.csv("milk_FA.normalized.csv", row.names = 1)

info.shg <- read.csv("milk_FA.info.humanShangHai.csv", row.names = 1)
info.msk <- read.csv("milk_FA.info.humanMoscow.csv", row.names = 1)

info.shg.filt <- info.shg[info.shg$LactationStage >=22 & info.shg$LactationStage<=40, ]
info.msk.filt <- info.msk[info.msk$LactationStage <=40, ]
info.msk.filt <- info.msk.filt[info.msk.filt$Parity %in% c("1", "2"), ]

info.hm <- rbind(info.shg.filt[, c("LactationStage", "Parity", "N_C", "Sex")], info.msk.filt[, c("LactationStage", "Parity", "N_C", "Sex")])

info.hm$Population <- factor(c(rep("Shanghai", nrow(info.shg.filt)), rep("Moscow", nrow(info.msk.filt))))
info.hm$Parity <- factor(info.hm$Parity)
info.hm$N_C <- factor(info.hm$N_C)
info.hm$Sex <- factor(info.hm$Sex)

info.hm$species <- c(rep("HSs", nrow(info.shg.filt)), rep("HSm", nrow(info.msk.filt)))
info.hm$species <- factor(info.hm$species, levels = c("HSm", "HSs"), order = TRUE)

all.species <- unique(as.character(info.hm$species))

logdata.hm <- log2(data[, rownames(info.hm)])


milkFactorAttribution <- function(dat, inf)
{
cf <- dat
cf[] <- NA
result <- c()
for(i in 1:nrow(dat)) {
    ints <- as.numeric(as.character(dat[i,]))
    allDat <- cbind(ints, inf)
    
    lmmod1 <- lm(ints ~ log2(LactationStage) + Parity + N_C + Sex + Population, data = allDat)
    aov <- Anova(lmmod1, type = 2)
    coef <- coef(lmmod1)
    names(coef) <- c("Intercept", "LS", "Parity", "N_C", "Sex", "Population")
    cf[i, ] <- ints - coef["LS"]*log2(inf$LactationStage) - coef["Parity"]*as.numeric(inf$Parity) - coef["N_C"]*as.numeric(inf$N_C) - coef["Sex"]*as.numeric(inf$Sex)
    
    sq <- aov[1:5, "Sum Sq"]/sum(aov[, "Sum Sq"])
    pr <- aov[1:5, "Pr(>F)"]
    fdr <- rep("NA", 5)
    sub <- c(as.character(sq), as.character(pr), fdr)
    result <- rbind(result, sub)
}

colnames(result) <- c("LactationStage.R2", "Parity.R2", "N_C.R2", "Sex.R2", "Population.R2", "LactationStage.Pvalue", "Parity.Pvalue", "N_C.Pvalue", "Sex.Pvalue", "Population.Pvalue", "LactationStage.FDR", "Parity.FDR", "N_C.FDR", "Sex.FDR", "Population.FDR")
result[, 11:15] <- apply(result[, 6:10], 2, function(x) {x[is.na(x)] <- 1; p.adjust(x, "fdr")})
result <- data.frame(apply(result, 2, function(x) as.numeric(as.character(x))))
rownames(result) <- rownames(dat)

result$Population.sig <- ifelse(result$Population.FDR <0.05, "yes", "no")

return(list(result, cf))
}


attrib <- milkFactorAttribution(logdata.hm, info.hm)

write.csv(data.frame(attrib[2]), "milk_FA.LS22_40.lmAdjustedNoPopulation.csv", row.names = TRUE, quote = FALSE)


specificity <- data.frame(attrib[1])

proportionPercentage <- function(dat)
{
pro.sp <- c()
for(i in 1:length(all.species)) {
    sampn.sp <- rownames(info.hm[info.hm$species==all.species[i], ])
    dt.sp <- dat[, sampn.sp]
    pro.tmp <- rowMeans(dt.sp, na.rm = TRUE)
    pro.sp <- cbind(pro.sp, pro.tmp)
}

colnames(pro.sp) <- all.species
pro.sp[is.nan(pro.sp)] <- 0
mean.sp <- pro.sp

return(pro.sp)
}


spMe.data <- data.frame(proportionPercentage(data))
specificity$HSm <- spMe.data[rownames(specificity), "HSm"] 
specificity$HSs <- spMe.data[rownames(specificity), "HSs"] 
specificity$pop.sig <- apply(specificity, 1, function(x) if(x["Population.sig"]=="yes") {ifelse(x["HSm"] >x["HSs"], "HSm", "HSs")} else {NA})

write.csv(specificity, "milk_FA.specific.humans.csv", row.names = TRUE, quote = FALSE)
