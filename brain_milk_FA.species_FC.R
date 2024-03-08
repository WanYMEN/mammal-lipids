rm(list=ls())

library(ggplot2)
library(reshape2)
library(stringr)
library(gridExtra)
library(ggrepel)

setwd("D://work//skoltech//lipid//writing//GitHub//data")

DATA.milk <- read.csv("milk_FA.normalized.csv", row.names = 1)
INFO.milk <- read.csv("milk_FA.info.csv", row.names = 1)

info.H <-  read.csv("milk_FA.info.human.csv", row.names = 1)
info.H <- info.H[info.H$LactationStage<=15, ]
info.H$species <- str_replace_all(info.H$species, "HSm|HSs", "HS")

info.ot <- INFO.milk[INFO.milk$species %in% c("MM", "CA", "SS"), ]
info.ot$LactationStage <- NA

info.milk <- rbind(info.H, info.ot)
info.milk$species <- factor(info.milk$species, levels = c("HS", "MM", "CA", "SS"), order = TRUE)
info.milk <- info.milk[order(info.milk$species), ]


DATA.brain <- read.csv("brain_FA.normalized.csv", row.names = 1)
info.brain <- read.csv("brain_FA.info.csv", row.names = 1)

info.brain <- info.brain[info.brain$species %in% c("HS", "MM", "CA", "SS"), ]
info.brain$species <- factor(info.brain$species, levels = c("HS", "MM", "CA", "SS"), order = TRUE)
info.brain <- info.brain[order(info.brain$species), ]

info.PFC <- info.brain[info.brain$region=="PFC", ]
info.CB <- info.brain[info.brain$region=="CB", ]

bothFA <- intersect(rownames(DATA.milk), rownames(DATA.brain))

data.milk <- DATA.milk[bothFA, rownames(info.milk)]
data.PFC <- DATA.brain[bothFA, rownames(info.PFC)]
data.CB <- DATA.brain[bothFA, rownames(info.CB)]

allSpecies <- c("HS", "MM", "CA", "SS")


## brain and milk comparison of Fold Change between species
FoldChangesLS2Week <- function(data.brain, info.brain, brainRegion)
{
result <- c()
for(j in 1:length(bothFA)) { 
    logFC <- c()
    ncd=0
    sub <- c()
    cn <- c()

    for(k in 1:length(allSpecies)) { 
        ncd=ncd+1
        ints.brain.1 <- as.numeric(as.character(data.brain[j, rownames(info.brain)[info.brain$species==allSpecies[k]]]))
        ints.milk.1 <- as.numeric(as.character(data.milk[j, rownames(info.milk)[info.milk$species==allSpecies[k]]]))

        if(k<length(allSpecies)) {
          for(n in (k+1):length(allSpecies)) {
              ints.milk.2 <- as.numeric(as.character(data.milk[j, rownames(info.milk)[info.milk$species==allSpecies[n]]]))
              logFC[ncd] <- log2(mean(ints.milk.1, na.rm = TRUE)/mean(ints.milk.2, na.rm = TRUE))
              cn <- cbind(cn, paste("log2FC.milk(", allSpecies[k], "/", allSpecies[n], ")", sep=""))
              sub <- cbind(sub, logFC[ncd])
              ncd=ncd+1
              
              ints.brain.2 <- as.numeric(as.character(data.brain[j, rownames(info.brain)[info.brain$species==allSpecies[n]]]))
              logFC[ncd] <- log2(mean(ints.brain.1, na.rm = TRUE)/mean(ints.brain.2, na.rm = TRUE))
              cn <- cbind(cn, paste("log2FC.", brainRegion, "(", allSpecies[k], "/", allSpecies[n], ")", sep=""))
              sub <- cbind(sub, logFC[ncd])
              ncd=ncd+1
          }
        }
   }
 result <- rbind(result, sub)
}

colnames(result) <- as.character(cn)
rownames(result) <- bothFA
return(result)

}

result.PFC <- FoldChangesLS2Week(data.PFC, info.PFC, "PFC")
result.CB <- FoldChangesLS2Week(data.CB, info.CB, "CB")

write.csv(result.PFC, "brain_milk_FA.species_FC.PFC.csv", row.names = TRUE, quote = FALSE)
write.csv(result.CB, "brain_milk_FA.species_FC.CB.csv", row.names = TRUE, quote = FALSE)

