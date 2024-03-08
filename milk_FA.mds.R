rm(list=ls())

library(ggplot2)
library(dplyr)
library(stringr)

setwd("D://work//skoltech//lipid//writing//GitHub//data")

data <- read.csv("milk_FA.normalized.csv", row.names = 1)
data <- log2(data)

info <- read.csv("milk_FA.info.csv", row.names = 1)

info$species <- factor(info$species, levels = c("HSm", "HSs", "MM", "MF", "BT", "BG", "CA", "SS", "BF", "QC"), order = TRUE)
info$group <- factor(info$group, levels = c("PR", "BO", "CA", "SS", "BF", "QC"), order = TRUE)

info$species.1 <- info$species
info$species.1 <- str_replace_all(info$species.1, "HSm|HSs", "HS")
info$species.1 <- factor(info$species.1, levels = c("HS", "MM", "MF", "BT", "BG", "CA", "SS", "BF", "QC"), order = TRUE)

rownames(info)[which(!(rownames(info) %in% colnames(data)))]
colnames(data)[which(!(colnames(data) %in% rownames(info)))]

colors.s <- c("HS" = "#FF0000", "HSm" = "#CD0000", "HSs" = "#EE4000", "MM" = "#FF8C00", "MF" = "#FFD700", "BT" = "#00FF00", "BG" = "#32CD32", "CA" = "#00FFFF", "SS" = "#6495ED", "BF" = "#B452CD", "QC" = "#696969")

data <- data[, rownames(info)]

MDS <- function(dat, inf, clct, lab, COLOR, x_lim, x_bre, y_lim, y_bre, wid, len)
{
mds <- t(dat) %>% scale(center = TRUE, scale = TRUE) %>% dist() %>% cmdscale(eig = TRUE)
mds.p <- data.frame(cbind(mds$points, inf[, clct]))
colnames(mds.p) <- c("Dim1", "Dim2", clct)
expl.var <- round(mds$eig/sum(mds$eig)*100)

pdf(paste("milk_FA.mds.", lab, ".pdf", sep=""), w=wid, h=len)
print(ggplot(mds.p, aes(x=Dim1, y=Dim2, color = inf[, clct])) + geom_point(size=0.8) + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.border = element_rect(colour = "black"), axis.title = element_text(size=9, color = "black"), axis.text = element_text(size=7, color = "black"), axis.line = element_line(colour = "black", size = 0.6), axis.ticks = element_line(colour = "black", size = 0.6), legend.spacing.y = unit(1, 'mm'), legend.text=element_text(size=rel(0.5)), legend.key.size = unit(3, "mm")) + xlab(paste("Dim1 ", "(", expl.var[1], "%)", sep="")) + ylab(paste("Dim2 ", "(", expl.var[2], "%)", sep="")) + scale_colour_manual(name = clct, values = COLOR) + scale_x_continuous(limits= x_lim, expand = c(0, 0), breaks = x_bre) + scale_y_continuous(limits=y_lim, expand = c(0, 0), breaks = y_bre))
dev.off()
}


## all
subInfo.1 <- info[info$group != "QC", ] 
subData.1 <- data[, rownames(subInfo.1)]

MDS(subData.1, subInfo.1, "species.1", "all", colors.s, c(-13, 17), seq(-13, 17, 15), c(-8, 14), seq(-8, 14, 11), 2.7, 1.6)

## phyla
subInfo.2 <- info[!(info$group %in% c("BF", "QC")), ] 
subData.2 <- data[, rownames(subInfo.2)]

colors.g <- c("PR" = "#FF6347", "BO" = "#3CB371", "CA" = "#00FFFF", "SS" = "#6495ED")

MDS(subData.2, subInfo.2, "group", "phyla", colors.g,  c(-13, 17), seq(-13, 17, 15), c(-8, 14), seq(-8, 14, 11), 2.45, 1.6)

## primates
subInfo.p1 <- info[info$group == "PR", ] 
subData.p1 <- data[, rownames(subInfo.p1)]

MDS(subData.p1, subInfo.p1, "species.1", "primates", colors.s,  c(-10, 20), seq(-10, 20, 15), c(-17, 9), seq(-17, 9, 13), 2.76, 1.6)


## only two human population
subInfo.HS <- info[info$species %in% c("HSm", "HSs"), ] 
subData.HS <- data[, rownames(subInfo.HS)]

MDS(subData.HS, subInfo.HS, "species", "humans", colors.s, c(-10, 20), seq(-10, 20, 15), c(-15, 11), seq(-15, 11, 13), 2.64, 1.6)


