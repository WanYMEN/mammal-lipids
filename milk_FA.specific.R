
## it takes a long time to run this script, thus it is better to run on a super cluster ##

rm(list=ls())

library(dplyr)
library(stringr)

setwd("D://work//skoltech//lipid//writing//GitHub//data")
setwd("/home/wangym/work/skoltech/lipid/writing/GitHub/data")

data <- read.csv("milk_FA.normalized.csv", row.names = 1)
info <- read.csv("milk_FA.info.csv", row.names = 1)

info <- info[!(info$group %in% c("QC", "BF")), ] 

level.s <- c("HS", "MM", "MF", "BT", "BG", "CA", "SS")
level.g <- c("PR", "BO", "CA", "SS")
level.p <- c("HS", "MM", "MF")

info$species <- str_replace_all(info$species, "HSm|HSs", "HS")
info$species <- factor(info$species, levels = level.s, order = TRUE)
info$group <- factor(info$group, levels = level.g, order = TRUE)

data <- data[, rownames(info)]


##
FAspecific <- function(dat, inf, item, Flevel)
{

result <- matrix(NA, nrow(dat), 6)
colnames(result) <- c("maxSigSpecies", "mean.tg", "mean.ot", "FC", "FCup1.5", "tTest.pvalue")
rownames(result) <- rownames(dat)

for(i in 1:nrow(dat)){
    dt.tmp <- data.frame(cbind(as.numeric(as.character(dat[i, ])), as.character(inf[, item])))
    colnames(dt.tmp) <- c("value", "species")
    dt.tmp$value <- as.numeric(as.character(dt.tmp$value))
    dt.tmp$species <- factor(dt.tmp$species, levels = Flevel, order = TRUE)

    sp.me <- aggregate(. ~ species, dt.tmp, mean, na.rm=TRUE)
    rownames(sp.me) <- sp.me$species
    max.sp <- rownames(sp.me)[sp.me$value==max(sp.me$value)]

    if(length(max.sp)>1){
       m=0
       for(j in 1:length(max.sp)){
           me.tg <- mean(dt.tmp[dt.tmp$species==max.sp[j], ]$value, na.rm=TRUE)
           me.ot <- mean(dt.tmp[dt.tmp$species!=max.sp[j], ]$value, na.rm=TRUE)
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
       dt.tg <- dt.tmp[dt.tmp$species==max.sp, "value"]
       dt.ot <- dt.tmp[dt.tmp$species!=max.sp, "value"]
       
       me.tg <- mean(dt.tg, na.rm=TRUE)
       me.ot <- mean(dt.ot, na.rm=TRUE)
       FC.tmp <- me.tg/me.ot
       ot.sp <- Flevel[Flevel!=max.sp]
       
       n=0
       for(k in 1:length(ot.sp)){
           if(!(ot.sp[k] %in% unique(as.character(sp.me$species)))){
              n = n + 1
           } else {
              if((me.tg/sp.me[sp.me$species==ot.sp[k], ]$value) >=1.5){ n = n + 1 }
           }
       }
       
       tTest.p <- tryCatch({t.test(log2(dt.tg), log2(dt.ot), paired=FALSE)$p.value}, error = function(e) {ifelse(grep("not enough 'y' observations", as.character(e)), 0, 1)})
       result[i, ] <- c(max.sp, log2(me.tg), log2(me.ot), FC.tmp, n, tTest.p)
    } else { result[i, ] <- c(NA, NA, NA, NA, NA, NA) }

}

return(result)

}


##
FAspecificImputation <- function(mat, dat, inf, item)
{
dat <- dat[rownames(mat), ]
mat$impute1000 = 0
sigFA <- nrow(mat[mat$tTest.FDR <0.05, ])

for(m in 1:1000){
    value.tmp <- matrix(NA, nrow(dat), 2)
    rownames(value.tmp) <- rownames(dat)
    for(i in 1:nrow(dat)){
        tg <- nrow(inf[inf[, item]==as.character((mat[i,]$maxSigSpecies)), ])
        smp.tg <- sample(rownames(inf), tg)
        smp.ot <- rownames(inf)[!(rownames(inf) %in% smp.tg)]
        
        rd.tg <- dat[i, smp.tg]
        rd.ot <- dat[i, smp.ot]
        
        value.tmp[i, 1] <- tryCatch({t.test(log2(rd.tg), log2(rd.ot), paired=FALSE)$p.value}, error = function(e) {ifelse(grep("not enough 'y' observations", as.character(e)), 0, 1)})
    }
    value.tmp[, 2] <- p.adjust(as.numeric(as.character(value.tmp[, 1])), method = "BH")
    sig.tmp <- ifelse(value.tmp[, 2] <0.05, 1, 0)
    mat$impute1000 = mat$impute1000 + sig.tmp
}
return(mat)
}


## species
species.sp <- data.frame(FAspecific(data, info, "species", level.s))
species.sp <- species.sp[!is.na(species.sp$maxSigSpecies) & as.numeric(as.character(species.sp$FC)) >=1.5 & species.sp$FCup1.5==(length(level.s) - 1), ]
species.sp$tTest.FDR <- p.adjust(as.numeric(as.character(species.sp$tTest.pvalue)), method = "BH")

imputeSpecies <- FAspecificImputation(species.sp, data, info, "species")

impt.s <- data.frame(imputeSpecies)
impt.s$iRatio <- impt.s$impute1000/1000
impt.s$sig <- impt.s$maxSigSpecies
impt.s[impt.s$tTest.FDR >=0.05 | impt.s$iRatio >=0.05, ]

write.csv(impt.s, "milk_FA.specific.species.csv", row.names = TRUE, quote = FALSE)


## group
group.sp <- data.frame(FAspecific(data, info, "group", level.g))
group.sp <- group.sp[!is.na(group.sp$maxSigSpecies) & as.numeric(as.character(group.sp$FC)) >=1.5 & group.sp$FCup1.5==(length(level.g) - 1), ]
group.sp$tTest.FDR <- p.adjust(as.numeric(as.character(group.sp$tTest.pvalue)), method = "BH")

imputeGroup <- FAspecificImputation(group.sp, data, info, "group")

impt.g <- data.frame(imputeGroup)
impt.g$iRatio <- impt.g$impute1000/1000
impt.g$sig <- impt.g$maxSigSpecies
impt.g[impt.g$tTest.FDR >=0.05 | impt.g$iRatio >=0.05, ]

write.csv(impt.g, "milk_FA.specific.group.csv", row.names = TRUE, quote = FALSE)


## primates
subInfo.p <- info[info$group == "PR", ] 
subData.p <- data[, rownames(subInfo.p)]

primate.sp <- data.frame(FAspecific(subData.p, subInfo.p, "species", level.p))
primate.sp <- primate.sp[!is.na(primate.sp$maxSigSpecies) & as.numeric(as.character(primate.sp$FC)) >=1.5 & primate.sp$FCup1.5==(length(level.p) - 1), ]
primate.sp$tTest.FDR <- p.adjust(as.numeric(as.character(primate.sp$tTest.pvalue)), method = "BH")

imputePrimate <- FAspecificImputation(primate.sp, subData.p, subInfo.p, "species")

impt.p <- data.frame(imputePrimate)
impt.p$iRatio <- impt.p$impute1000/1000
impt.p$sig <- impt.p$maxSigSpecies
impt.p[impt.p$tTest.FDR >=0.05 | impt.p$iRatio >=0.05, ]

write.csv(impt.p, "milk_FA.specific.primate.csv", row.names = TRUE, quote = FALSE)

