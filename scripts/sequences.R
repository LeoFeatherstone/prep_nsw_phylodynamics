#############################################################
# Write sequences to separate FASTA files based on clusters #
#############################################################

library(readxl)
library(ape)

ann <- read_excel("sequences/large_clusters.xlsx")
seq <- read.FASTA("sequences/clusters_combined.fasta")

for (c in unique(ann$Cluster)) {
  ids <- ann$New_name_large_cluster[ann$Cluster == c]
  sel <- seq[names(seq) %in% ids]
  if (length(sel) > 0) {
    write.FASTA(sel, file = paste0("sequences/", c, ".fasta"))
  }
}
