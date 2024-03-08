rm(list=ls())

setwd("D://work//skoltech//lipid//writing//GitHub//data")

library(ggplot2)
library(reshape2)
library(stringr)
library(gridExtra)

DATA <- read.csv("milk_FA.normalized.csv", row.names = 1)

info <- read.csv("milk_FA.info.csv", row.names = 1)
info <- info[!(info$species %in% c("BF", "QC")), ]

info$species <- str_replace_all(info$species, "HSm|HSs", "HS")
info$species <- factor(info$species, levels = c("HS", "MM", "MF", "BT", "BG", "CA", "SS"), order = TRUE)
info <- info[order(info$species), ]

all.species <- unique(as.character(info$species))
all.FAs <- rownames(DATA)

data <- DATA[, rownames(info)]

colors <- c("HS" = "#FF0000", "MM" = "#FF8C00", "MF" = "#FFD700", "BT" = "#00FF00", "BG" = "#32CD32", "CA" = "#00FFFF", "SS" = "#6495ED")


sm.FA <- data.frame(FA = all.FAs, clength = unlist(strsplit(all.FAs, ":"))[2*(1:length(all.FAs))-1], dbond = unlist(strsplit(all.FAs, ":"))[2*(1:length(all.FAs))])
rownames(sm.FA) <- all.FAs

sm.FA$FA <- factor(sm.FA$FA , levels = all.FAs, order = TRUE)

sm.FA$clength <- as.numeric(as.character(sm.FA$clength))
sm.FA$dbond <- as.numeric(as.character(sm.FA$dbond))

## SFA:saturated; even-unsaturated:EUFA; odd-unsaturated: OUFA; poly-unsaturated: PUFA

sm.FA$group <- as.character(apply(sm.FA, 1, function(x) {if(x[3]==0) {"SFA"} else if(as.numeric(x[3]) <5 & (as.numeric(x[2]) %% 2 == 0)) {"EUFA"} else if(as.numeric(x[3]) <5 & (as.numeric(x[2]) %% 2 != 0)) {"OUFA"} else {"PUFA"}}))

sm.FA$group <- factor(sm.FA$group , levels = c("SFA", "EUFA", "OUFA", "PUFA"), order = TRUE)


## species percentage of proportion
pro.sp <- c()
for(i in 1:length(all.species)) {
    sampn.sp <- rownames(info[info$species==all.species[i], ])
    dt.sp <- data[, sampn.sp]
    pro.tmp <- rowMeans(dt.sp, na.rm = TRUE)
    pro.sp <- cbind(pro.sp, pro.tmp)
}

colnames(pro.sp) <- all.species
pro.sp[is.nan(pro.sp)] <- 0
pro.sp <- apply(pro.sp, 2, function(x) x/sum(x))
head(pro.sp)
pro.sp <- apply(pro.sp, 1, function(x) (x/sum(x))*100)
mt.pro.sp <- melt(pro.sp)
head(mt.pro.sp)

mt.pro.sp$Var1 <- factor(mt.pro.sp$Var1, levels = c("HS", "MM", "MF", "BT", "BG", "CA", "SS"), order = TRUE)
mt.pro.sp$Var2 <- factor(mt.pro.sp$Var2, levels = all.FAs, order = TRUE)

mt.pro.sp$group <- sm.FA[as.character(mt.pro.sp$Var2), "group"]
mt.pro.sp$group <- factor(mt.pro.sp$group, levels = c("SFA", "EUFA", "OUFA", "PUFA"), order = TRUE)

pdf("milk_FA.percentage.species.pdf",w=3, h=1.4)
ggplot() + geom_bar(data = mt.pro.sp, aes(x=Var2, y = value, fill=Var1), stat = "identity", width=1) + ylab("PercentageOfProportion") + xlab("") + theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size=9), axis.text.y = element_text(size=7, color = "black"), axis.ticks.y = element_line(colour = "black", size = 0.6), panel.spacing = unit(0, "lines"), strip.text.x = element_text(margin = margin(0.5,0,0.5,0, "mm"), size=6, colour = "black"), panel.border = element_rect(size = 0.3, colour = "black"), strip.background = element_rect(size = 0.3, colour = "black"), plot.margin = unit(c(0, .5, -1.1, 0), "lines"), legend.position="none") + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(name = "species", values = colors) + facet_grid(~ group, scales="free", space="free_x")
dev.off()

