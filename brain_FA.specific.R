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


FAspecific <- function(dat, inf, item, Flevel)
{

result <- matrix(NA, nrow(dat), 6)
colnames(result) <- c("maxSigSpecies", "mean.tg", "mean.ot", "FC", "FCup1.5", "aovPV.species")
rownames(result) <- rownames(dat)

for(i in 1:nrow(dat)){
    dt.tmp <- data.frame(cbind(as.numeric(as.character(dat[i, ])), as.character(inf[, "sex"]), inf[, "logScaledDay"], as.character(inf[, item])))
    colnames(dt.tmp) <- c("value", "sex", "age", "species")
    dt.tmp$value <- as.numeric(as.character(dt.tmp$value))
    dt.tmp$age <- as.numeric(as.character(dt.tmp$age))
    dt.tmp$sex <- factor(dt.tmp$sex)
    dt.tmp$species <- as.character(dt.tmp$species)

    sp.me <- aggregate(. ~ species, dt.tmp, mean)
    rownames(sp.me) <- sp.me$species
    max.sp <- rownames(sp.me)[sp.me$value==max(sp.me$value)]

    if(length(max.sp)>1){
       m=0
       for(j in 1:length(max.sp)){
          dt.tmp$species.1 <- dt.tmp$species
          dt.tmp[dt.tmp$species.1!=max.sp[j], ]$species.1 <- "others"
          me.tg <- mean(dt.tmp[dt.tmp$species.1==max.sp[j], ]$value, na.rm=TRUE)
          me.ot <- mean(dt.tmp[dt.tmp$species.1!=max.sp[j], ]$value, na.rm=TRUE)
          FC.tmp <- me.tg/me.ot
          ot.sp <- Flevel[Flevel!=max.sp[j]]
 
          n=0
          for(k in 1:length(ot.sp)){
              if(!(ot.sp[k] %in% unique(as.character(sp.me$species)))){
                  n = n + 1
              } else {
                  if((me.tg/sp.me[sp.me$species==ot.sp[k], ]$value) >=1.5){ n = n + 1 }
              }
          }
           
          if(FC.tmp >=1.5 & n==length(ot.sp)){
             max.sp.tmp <- max.sp[j]
             m = m + 1
          }
       }
       if(m==1){ max.sp <- max.sp.tmp }
    }

    if(length(max.sp)==1){
       dt.tmp$species.1 <- dt.tmp$species
       dt.tmp[dt.tmp$species.1!=max.sp, ]$species.1 <- "others"
       
       me.tg <- mean(dt.tmp[dt.tmp$species.1==max.sp, ]$value, na.rm=TRUE)
       me.ot <- mean(dt.tmp[dt.tmp$species.1!=max.sp, ]$value, na.rm=TRUE)
       FC.tmp <- me.tg/me.ot
       dt.tmp$species.1 <- factor(dt.tmp$species.1, levels = c(max.sp, "others"), order = TRUE)
       ot.sp <- Flevel[Flevel!=max.sp]
       
       n=0
       for(k in 1:length(ot.sp)){
           if(!(ot.sp[k] %in% unique(as.character(sp.me$species)))){
              n = n + 1
           } else {
              if((me.tg/sp.me[sp.me$species==ot.sp[k], ]$value) >=1.5){ n = n + 1 }
           }
       }
       
       AV <- Anova(lm(log2(value) ~ sex + age + species.1, data = dt.tmp), type = 2)
       AV.p <- AV["species.1", "Pr(>F)"]
       result[i, ] <- c(max.sp, log2(me.tg), log2(me.ot), FC.tmp, n, AV.p)
    } else { result[i, ] <- c(NA, NA, NA, NA, NA, NA) }

}

return(result)

}


## PFC ##
## species ##
sp.s.PFC <- data.frame(FAspecific(data.PFC, info.PFC, "species", level.s))
sp.s.PFC <- sp.s.PFC[!is.na(sp.s.PFC$maxSigSpecies) & as.numeric(as.character(sp.s.PFC$FC)) >=1.5 & sp.s.PFC$FCup1.5==(length(level.s) - 1), ]
sp.s.PFC$aovQV.species <- p.adjust(as.numeric(as.character(sp.s.PFC$aovPV.species)), method = "BH")

sig.s.PFC <- sp.s.PFC[sp.s.PFC$aovQV.species <0.05, ]
write.csv(sig.s.PFC, "brain_FA.specific.species.PFC.csv", row.names = TRUE, quote = FALSE)

## group ##
sp.g.PFC <- data.frame(FAspecific(data.PFC, info.PFC, "group", level.g))
sp.g.PFC <- sp.g.PFC[!is.na(sp.g.PFC$maxSigSpecies) & as.numeric(as.character(sp.g.PFC$FC)) >=1.5 & sp.g.PFC$FCup1.5==(length(level.g) - 1), ]
sp.g.PFC$aovQV.species <- p.adjust(as.numeric(as.character(sp.g.PFC$aovPV.species)), method = "BH")

sig.g.PFC <- sp.g.PFC[sp.g.PFC$aovQV.species <0.05, ]
write.csv(sig.g.PFC, "brain_FA.specific.group.PFC.csv", row.names = TRUE, quote = FALSE)

## primates ##
subInfo.PFC <- info.PFC[info.PFC$group == "PR", ]
subData.PFC <- data.PFC[, rownames(subInfo.PFC)]

sp.p.PFC <- data.frame(FAspecific(subData.PFC, subInfo.PFC, "species", level.p))
sp.p.PFC <- sp.p.PFC[!is.na(sp.p.PFC$maxSigSpecies) & as.numeric(as.character(sp.p.PFC$FC)) >=1.5 & sp.p.PFC$FCup1.5==(length(level.p) - 1), ]
sp.p.PFC$aovQV.species <- p.adjust(as.numeric(as.character(sp.p.PFC$aovPV.species)), method = "BH")

sig.p.PFC <- sp.p.PFC[sp.p.PFC$aovQV.species <0.05, ]
write.csv(sig.p.PFC, "brain_FA.specific.primate.PFC.csv", row.names = TRUE, quote = FALSE)



## CB ##
## species ##
sp.s.CB <- data.frame(FAspecific(data.CB, info.CB, "species", level.s))
sp.s.CB <- sp.s.CB[!is.na(sp.s.CB$maxSigSpecies) & as.numeric(as.character(sp.s.CB$FC)) >=1.5 & sp.s.CB$FCup1.5==(length(level.s) - 1), ]
sp.s.CB$aovQV.species <- p.adjust(as.numeric(as.character(sp.s.CB$aovPV.species)), method = "BH")

sig.s.CB <- sp.s.CB[sp.s.CB$aovQV.species <0.05, ]
write.csv(sig.s.CB, "brain_FA.specific.species.CB.csv", row.names = TRUE, quote = FALSE)

## group ##
sp.g.CB <- data.frame(FAspecific(data.CB, info.CB, "group", level.g))
sp.g.CB <- sp.g.CB[!is.na(sp.g.CB$maxSigSpecies) & as.numeric(as.character(sp.g.CB$FC)) >=1.5 & sp.g.CB$FCup1.5==(length(level.g) - 1), ]
sp.g.CB$aovQV.species <- p.adjust(as.numeric(as.character(sp.g.CB$aovPV.species)), method = "BH")

sig.g.CB <- sp.g.CB[sp.g.CB$aovQV.species <0.05, ]
write.csv(sig.g.CB, "brain_FA.specific.group.CB.csv", row.names = TRUE, quote = FALSE)

## primates ##
subInfo.CB <- info.CB[info.CB$group == "PR", ]
subData.CB <- data.CB[, rownames(subInfo.CB)]

sp.p.CB <- data.frame(FAspecific(subData.CB, subInfo.CB, "species", level.p))
sp.p.CB <- sp.p.CB[!is.na(sp.p.CB$maxSigSpecies) & as.numeric(as.character(sp.p.CB$FC)) >=1.5 & sp.p.CB$FCup1.5==(length(level.p) - 1), ]
sp.p.CB$aovQV.species <- p.adjust(as.numeric(as.character(sp.p.CB$aovPV.species)), method = "BH")

sig.p.CB <- sp.p.CB[sp.p.CB$aovQV.species <0.05, ]
write.csv(sig.p.CB, "brain_FA.specific.primate.CB.csv", row.names = TRUE, quote = FALSE)


