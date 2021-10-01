library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyverse)
library(readr)
summaryfolder = "summary"

read_counts <- function(type) {
  fname = file.path(summaryfolder,sprintf("%s_%s.tsv",type,"counts"))
  tbl <- read_tsv(fname)
  return(tbl)
}


types <- c("CAZY","MEROPS","Pfam")

d <- lapply(types,read_counts)

CAZY <- pivot_longer(data = d[[1]], 
                     cols = -c(1:1),
                     names_to = "Species", 
                     values_to = "Abundance")

MEROPS <- pivot_longer(data = d[[2]], 
                     cols = -c(1:1),
                     names_to = "Species", 
                     values_to = "Abundance")

Pfam <- pivot_longer(data = d[[3]], 
                       cols = -c(1:2),
                       names_to = "Species", 
                       values_to = "Abundance")

CAZY.heatmap <- ggplot(data = CAZY, mapping = aes(x = Species,
                                                  y = DOMAIN,
                                                  fill = Abundance)) +
  geom_tile() +
  xlab(label = "Species") + theme(axis.text.x=element_text(angle = -90, hjust = 0))

CAZY.heatmap
