library(stringr)
library(tidyverse)

FullIntronList <- read.table("Misc/FullIntronList.bed", sep='\t', col.names=c("chr", "start", "stop", "Intron", "score", "strand", "type"))
Donors <- read.table("Misc/MotifScores/Donors.txt", sep='\t', col.names = c("Intron", "donor_seq", "DonorScore"))
Acceptors <- read.table("Misc/MotifScores/Acceptors.txt", sep='\t', col.names = c("Intron", "acceptor_seq", "AcceptorScore"))
Branchpoints <- read.table("Misc/MotifScores/Branchpoints.bestPerIntron.txt", sep='\t', col.names=c("seqid", "eqez", "ss_dist", "bp_seq", "bp_scr", "y_cont", "ppt_off", "ppt_len", "ppt_scr", "svm_scr")) %>%
  mutate(Intron=str_sub(seqid, end=-4))

FullIntronList %>%
  left_join(Donors, by="Intron") %>%
  left_join(Acceptors, by="Intron") %>%
  left_join(Branchpoints, by="Intron") %>%
  dplyr::select(-c("chr", "start", "stop", "score", "seqid")) %>%
  write.table(file="../output/IntronFeatures.txt", sep='\t', row.names = F, quote=F)
