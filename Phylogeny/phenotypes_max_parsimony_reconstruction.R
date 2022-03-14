library(ape)
library(castor)
drugs <- c('Isoniazid', 'Rifampicin', 'Pyrazinamide', 'Ethambutol', 'Streptomycin', 
           'Amikacin', 'Capreomycin', 'Ofloxacin', 'Moxifloxacin')
for (i in drugs){
  phenotypes <- read.table(paste0("input_data/", i, ".pheno"), sep="\t", header=F)
  tree_labeled <- read.tree(file=paste0("../../snp_small_trees/", i, ".snps_snp_tree.treefile"))
  phenotypes <- phenotypes[phenotypes$V1 %in% tree_labeled$tip.label,]
  x <- phenotypes$V2
  names(x) <- phenotypes$V1
  x_ordered <- x[order(factor(names(x), levels = tree_labeled$tip.label))]
  anc_rec <- asr_max_parsimony(tree_labeled, x_ordered+1)
  df <- data.frame(anc_rec$ancestral_likelihoods)
  anc_pheno <- apply(df, 1, which.max) - 1
  result <- cbind(tree_labeled$node.label, anc_pheno)
  write.table(result, file=paste0("mp_anc_rec/", i, ".mp.pheno"), sep="\t", row.names=F, col.names=F, quote=F, eol="\n")
}
