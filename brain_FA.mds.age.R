rm(list=ls())

library(ggplot2)
library(dplyr)
library(stringr)

setwd("D://work//skoltech//lipid//writing//GitHub//data")

data <- read.csv("brain_FA.normalized.csv", row.names = 1)
data <- log2(data)

info <- read.csv("brain_FA.info.csv", row.names = 1)
info <- info[info$region!="QC", ]
info$age <- as.numeric(as.character(info$age))

allSpecies <- c("HS","PT", "MM", "CA", "SS")
allRegion <- c("PFC", "CB")

info$species <- factor(info$species, levels = allSpecies, order = TRUE)

info$scaledAge <- info$age
info[info$species=="HS", ]$scaledAge <- info[info$species=="HS", ]$age/365 ## 1 year
info[info$species=="PT", ]$scaledAge <- info[info$species=="PT", ]$age/45 ## 1.5 month
info[info$species=="MM", ]$scaledAge <- info[info$species=="MM", ]$age/365 ## 1 year
info[info$species=="CA", ]$scaledAge <- info[info$species=="CA", ]$age/75 ## 2.5 month
info[info$species=="SS", ]$scaledAge <- info[info$species=="SS", ]$age/60 ## 2 month
info[info$scaledAge>=1, ]$scaledAge <- 1

colors <- c("HS" = "#FF0000", "PT"= "#A52A2A", "MM" = "#FF8C00", "CA" = "#00FFFF", "SS" = "#6495ED")


MDS <- function(dat, inf, lab, region, dim, width, height, x_lim, x_bre, y_lim, y_bre)
{
mds <- t(dat) %>% scale(center = TRUE, scale = TRUE) %>% dist() %>% cmdscale(k=3, eig = TRUE)
mds.p <- data.frame(mds$points)
colnames(mds.p) <- c("Dim1", "Dim2", "Dim3")
mds.p$species <- inf$species
mds.p$age <- inf$scaledAge
expl.var <- round(mds$eig/sum(mds$eig)*100)

if(dim=="12"){
    pdf(paste("brain_FA.mds", lab, region, "age.Dim12.pdf", sep="."), w=width, h=height)
    print(ggplot(mds.p, aes(x=Dim1, y=Dim2, color = species, size = age)) + geom_point(alpha = 0.5) + scale_colour_manual(name = "species", values = colors) + scale_size_continuous(name = "age/max", limits=c(0, 1), range = c(1, 3), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '>=1')) + xlab(paste("Dim1 ", "(", expl.var[1], "%)", sep="")) + ylab(paste("Dim2 ", "(", expl.var[2], "%)", sep="")) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"), axis.title = element_text(size=9, color = "black"), axis.text = element_text(size=7, color = "black"), axis.line = element_line(colour = "black", size = 0.6), axis.ticks = element_line(colour = "black", size = 0.6), legend.spacing.y = unit(1, 'mm'), legend.text=element_text(size=rel(0.7)), legend.key.size = unit(1, "mm")) + scale_x_continuous(limits= x_lim, expand = c(0, 0), breaks = x_bre) + scale_y_continuous(limits=y_lim, expand = c(0, 0), breaks = y_bre))
    dev.off()
} else if(dim=="13"){
   pdf(paste("brain_FA.mds", lab, region, "age.Dim13.pdf", sep="."), w=width, h=height)
   print(ggplot(mds.p, aes(x=Dim1, y=Dim3, color = species, size = age)) + geom_point(alpha = 0.5) + scale_colour_manual(name = "species", values = colors) + scale_size_continuous(name = "age/max", limits=c(0, 1), range = c(1, 3), breaks = c(0, 0.5, 1), labels = c('0', '0.5', '>=1')) + xlab(paste("Dim1 ", "(", expl.var[1], "%)", sep="")) + ylab(paste("Dim3 ", "(", expl.var[3], "%)", sep="")) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"), axis.title = element_text(size=9, color = "black"), axis.text = element_text(size=7, color = "black"), axis.line = element_line(colour = "black", size = 0.6), axis.ticks = element_line(colour = "black", size = 0.6), legend.spacing.y = unit(1, 'mm'), legend.text=element_text(size=rel(0.7)), legend.key.size = unit(1, "mm")) + scale_x_continuous(limits= x_lim, expand = c(0, 0), breaks = x_bre) + scale_y_continuous(limits=y_lim, expand = c(0, 0), breaks = y_bre))
   dev.off()
}
}


## species ## 
##PFC
info.PFC <- info[info$region=="PFC", ]
data.PFC <- data[, rownames(info.PFC)]

MDS(data.PFC, info.PFC, "species", "PFC", "13", 2.6, 1.55, c(-6, 11), seq(-6, 11, 8.5), c(-6.8, 6.4), seq(-6.8, 6.4, 6.6))


##CB
info.CB <- info[info$region=="CB", ]
data.CB <- data[, rownames(info.CB)]

MDS(data.CB, info.CB, "species", "CB", "12", 2.6, 1.55, c(-6.2, 7.8), seq(-6.2, 7.8, 7), c(-6.2, 7.8), seq(-6.2, 7.8, 7))


## primate ##
#PFC
subInfo.PFC <- info.PFC[info.PFC$group == "PR", ] 
subData.PFC <- data.PFC[, rownames(subInfo.PFC)]

MDS(subData.PFC, subInfo.PFC, "primate", "PFC", "12", 2.6, 1.55, c(-5, 14), seq(-5, 14, 9.5), c(-9, 8), seq(-9, 8, 8.5))

#CB
subInfo.CB <- info.CB[info.CB$group == "PR", ] 
subData.CB <- data.CB[, rownames(subInfo.CB)]

MDS(subData.CB, subInfo.CB, "primate", "CB", "12", 2.6, 1.55, c(-6, 9), seq(-6, 9, 7.5), c(-6, 9), seq(-6, 9, 7.5))

