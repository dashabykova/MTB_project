import sys
import numpy as np
from os import makedirs, listdir
from os.path import exists

if __name__ == "__main__":
    #drug = sys.argv[1]
    #out_path = '/slowhome/AI/fedonin.gg/dbykova/indels_parsimony_reconstructed/' + drug + '/'
    out_path = '/export/data/bykova/MTB/data/9drugs/indels_parsed/'
    makedirs(out_path, exist_ok=True)
    data_path = '/export/data/bykova/MTB/data/'
    index_indel = {}
    i = 0
    with open(data_path + 'combined_mq40_keep_complex_std_names_filtered_with_DR.indel_list') as indels:
        for line in indels:
            index_indel[i] = line
            i += 1
    node_indels = {}
    with open(data_path + '9drugs/9drugs_rec.indels') as f:
        for line in f:
            temp = line.strip().split('\t')
            if len(temp[1:]) == i:
                node_indels[temp[0]] = list(np.nonzero([int(j) for j in temp[1:]])[0])
            else:
                print('Node ' + temp[0] + ' has wrong number of indels')
                if (temp[0] in ['0', '1']) & (len(temp) == i):
                    node_indels['empty_node'] = list(np.nonzero([int(j) for j in temp])[0])
    for n in node_indels:
        with open(out_path + n + '.indels', 'w') as outf:
            for index in node_indels[n]:
                outf.write(index_indel[index])