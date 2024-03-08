rm(list=ls())

library(ggplot2)
library(reshape2)
library(stringr)
library(gridExtra)
library(cowplot)

setwd("D://work//skoltech//lipid//writing//GitHub//data")

data <- read.csv("brain_FA.normalized.csv", row.names = 1)

info <- read.csv("brain_FA.info.csv", row.names = 1)
info <- info[info$species != "QC", ]

level.s <- c("HS", "PT", "MM", "CA", "SS")
level.g <- c("PR", "CA", "SS")
level.p <- c("HS", "PT", "MM")

info$species <- factor(info$species, levels = level.s, order = TRUE)
info$group <- factor(info$group, levels = level.g, order = TRUE)

info.PFC <- info[info$region=="PFC", ]
info.CB <- info[info$region=="CB", ]

colors.s <- c("HS" = "#FF0000", "PT"= "#A52A2A", "MM" = "#FF8C00", "CA" = "#00FFFF", "SS" = "#6495ED")
colors.g <- c("PR" = "#FF6347", "CA" = "#00FFFF", "SS" = "#6495ED")

all.FAs <- rownames(data)


sm.FA <- data.frame(FA = all.FAs, clength = unlist(strsplit(all.FAs, ":"))[2*(1:length(all.FAs))-1], dbond = unlist(strsplit(all.FAs, ":"))[2*(1:length(all.FAs))])
rownames(sm.FA) <- all.FAs

sm.FA$clength <- as.numeric(as.character(sm.FA$clength))
sm.FA$dbond <- as.numeric(as.character(sm.FA$dbond))

## SFA:saturated; even-unsaturated:EUFA; odd-unsaturated: OUFA; longChain-unsaturated: LUFA
sm.FA$group <- as.character(apply(sm.FA, 1, function(x) {if(x[3]==0) {"SFA"} else if(as.numeric(x[2]) <26 & (as.numeric(x[2]) %% 2 == 0)) {"EUFA"} else if(as.numeric(x[2]) <26 & (as.numeric(x[2]) %% 2 != 0)) {"OUFA"} else {"LUFA"}}))

sm.FA$group <- factor(sm.FA$group , levels = c("SFA", "EUFA", "OUFA", "LUFA"), order = TRUE)


## species percentage of proportion
PtProp <- function(dat, inf, item, label, Flevel, region, colors)
{
pro.sp <- c()
for(i in 1:length(Flevel)) {
    sampn.sp <- rownames(inf[inf[, item]==Flevel[i], ])
    dt.sp <- dat[, sampn.sp]
    pro.tmp <- rowMeans(dt.sp, na.rm = TRUE)
    pro.sp <- cbind(pro.sp, pro.tmp)
}

colnames(pro.sp) <- Flevel
pro.sp <- apply(pro.sp, 2, function(x) x/sum(x))
pro.sp <- apply(pro.sp, 1, function(x) (x/sum(x))*100)
mt.pro.sp <- melt(pro.sp)

mt.pro.sp$Var1 <- factor(mt.pro.sp$Var1, levels = Flevel, order = TRUE)

mt.pro.sp$group <- sm.FA[as.character(mt.pro.sp$Var2), "group"]
mt.pro.sp$group <- factor(mt.pro.sp$group, levels = c("SFA", "EUFA", "OUFA", "LUFA"), order = TRUE)

pdf(paste("brain_FA", region, "percentage_bar", label, "pdf", sep="."), w=3.3, h=1.25)
print(ggplot() + geom_bar(data = mt.pro.sp, aes(x=Var2, y = value, fill=Var1), stat = "identity", width=1) + ylab("PercentageOfProportion") + xlab("") + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_text(size=7, color = "black"), axis.ticks.y = element_line(colour = "black", size = 0.6), axis.title.y = element_text(size=8), panel.spacing = unit(0, "lines"), strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm"), size=6), panel.border = element_rect(size = 0.3, colour = "black"), strip.background = element_rect(size = 0.35, colour = "black"), plot.margin = unit(c(0.1, 0.5, -1.1, 0), "lines"), legend.position="none") + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(name = "species", values = colors) + facet_grid(~ group, scales="free", space="free_x"))
dev.off()
}
##


data.PFC.s <- read.csv("brain_FA.sexAgeAdjust.species.PFC.csv", row.names = 1)
data.CB.s <- read.csv("brain_FA.sexAgeAdjust.species.CB.csv", row.names = 1)
data.PFC.s <- data.PFC.s[rownames(data), rownames(info.PFC)]
data.CB.s <- data.CB.s[rownames(data), rownames(info.CB)]

PtProp(data.PFC.s, info.PFC, "species", "species", level.s, "PFC", colors.s)
PtProp(data.CB.s, info.CB, "species", "species", level.s, "CB", colors.s)


subInfo.PFC <- info.PFC[info.PFC$group == "PR", ] 
subInfo.CB <- info.CB[info.CB$group == "PR", ] 

data.PFC.p <- read.csv("brain_FA.sexAgeAdjust.primate.PFC.csv", row.names = 1)
data.CB.p <- read.csv("brain_FA.sexAgeAdjust.primate.CB.csv", row.names = 1)

data.PFC.p <- data.PFC.p[rownames(data), rownames(subInfo.PFC)]
data.CB.p <- data.CB.p[rownames(data), rownames(subInfo.CB)]

PtProp(data.PFC.p, subInfo.PFC, "species", "primate", level.p, "PFC", colors.s)
PtProp(data.CB.p, subInfo.CB, "species", "primate", level.p, "CB", colors.s)

