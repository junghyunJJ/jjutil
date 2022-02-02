# https://rdrr.io/cran/RapidoPGS/man/EUR_ld.blocks38.html
library(tidyverse)
load("EUR_ld.blocks38.RData")
ld_blodk_eur38 <- data.frame(EUR_ld.blocks38)[,1:3]
colnames(ld_blodk_eur38) <- c("CHR", "START", "STOP")
ld_blodk_eur38$CHR <- sub("chr","",ld_blodk_eur38$CHR)
ld_blodk_eur38 <- as_tibble(ld_blodk_eur38)
save(ld_blodk_eur38, file = "RHOGE_ld_blodk_eur38.RData")

colnames(ld_blodk_eur38) <- c("chrom", "start", "stop")
write.table(ld_blodk_eur38, "focus_grch38.eur.loci.bed", row.names = F, quote = F, sep = " ")

load("EUR_ld.blocks.RData")
ld_blodk_eur19 <- data.frame(EUR_ld.blocks)[,1:3]
colnames(ld_blodk_eur19) <- c("CHR", "START", "STOP")
ld_blodk_eur19$CHR <- sub("chr","",ld_blodk_eur19$CHR)
ld_blodk_eur19 <- as_tibble(ld_blodk_eur19)
save(ld_blodk_eur19, file = "RHOGE_ld_blodk_eur19.RData")

colnames(ld_blodk_eur19) <- c("chrom", "start", "stop")
write.table(ld_blodk_eur19, "focus_grch37.eur.loci.bed", row.names = F, quote = F, sep = " ")
