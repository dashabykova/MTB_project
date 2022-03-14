import sys
from Bio import Phylo

def parent_childs_dict(tree):
    parents = {}
    for clade in tree.find_clades(order="level"):
        parents[str(clade)] = []
        for child in clade:
            parents[str(clade)].append(child.name)
    return parents

def compare(child_var1, child_var2, parent_variants):
    child1 = open(child_var1).readlines()
    child2 = open(child_var2).readlines()
    parent_var = open(parent_variants).readlines()
    diff = (set(child1) & set(child2)) - (set(child1) & set(child2) & set(parent_var))
    if len(diff) > 0:
        snps = 0
        indels = 0
        parent_name = parent_variants.split('/')[-1][:parent_variants.split('/')[-1].rfind('.')]
        for v in diff:
            pos, ref, alt = v.strip().split('\t')
            if len(ref) == len(alt):
                snps += 1
            else:
                indels += 1
        print("{0} has {1} snps and {2} indels missing in both children".format(parent_name, snps, indels))
        return 1
    return 0
    
if __name__ == "__main__":
    path_to_tree = sys.argv[1]
    path_to_anc_variants = sys.argv[2]
    path_to_leaves = sys.argv[3]
    tree = Phylo.read(path_to_tree, 'newick')
    parents = parent_childs_dict(tree)
    for parent, childs in parents.items():
        if len(childs) == 2:
            if (not childs[0].startswith('Node')) and (not childs[1].startswith('Node')):
                is_diff = compare(path_to_leaves + childs[0] + '.variants', path_to_leaves + childs[1] + '.variants', 
                                  path_to_anc_variants + parent + '.variants')
                if is_diff:
                    print(childs[0], childs[1], parent)