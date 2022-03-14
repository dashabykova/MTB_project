from ete3 import Tree
import sys

#tree to prune
path_to_tree = "/slowhome/AI/fedonin.gg/dbykova/big_dataset/iqtree_snp/combined_mq40_keep_complex_std_names_filtered_with_DR_snp.fasta.treefile"
t = Tree(path_to_tree, format=1)
#tips to keep
tips = "/slowhome/AI/fedonin.gg/dbykova/9drugs/9drugs.tips"
tips_list = [l.strip() for l in open(tips).readlines()]
t.prune(tips_list)
t.write(format=1, outfile="%s_pruned.nwk" % path_to_tree)