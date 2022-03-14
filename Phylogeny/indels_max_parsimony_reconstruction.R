library(dplyr)
library(castor)
library(ape)
all_indels <- read.table("/slowhome/AI/fedonin.gg/dbykova/big_dataset/combined_mq40_keep_complex_std_names_filtered_with_DR_indels.tsv", sep=" ", header=F)
tree_labeled <- read.tree("/slowhome/AI/fedonin.gg/dbykova/9drugs/9drugs.snps_pruned.nwk")
indels <- all_indels[all_indels$V1 %in% tree_labeled$tip.label,]
indels <- slice(indels, match(tree_labeled$tip.label, indels$V1))
x <- indels[,2]
anc_rec <- asr_max_parsimony(tree_labeled, x+1)
df <- data.frame(anc_rec$ancestral_likelihoods)
anc_indels <- apply(df, 1, which.max) - 1
result <- cbind(tree_labeled$node.label, anc_indels)
for (j in 3:ncol(indels)){
x <- indels[,j]
anc_rec <- asr_max_parsimony(tree_labeled, x+1)
df <- data.frame(anc_rec$ancestral_likelihoods)
anc_indels <- apply(df, 1, which.max) - 1
result <- cbind(result, anc_indels)
}
write.table(result, file="/slowhome/AI/fedonin.gg/dbykova/9drugs/9drugs_rec.indels", sep="\t", row.names=F, col.names=F, quote=F, eol="\n")


