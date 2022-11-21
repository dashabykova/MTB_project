## Epistatic interactions in Mycobacterium tuberculosis genes 
Here will be scripts for the research project

* Phylogeny: ancestor states reconstruction and reconstruction results parsing
* Downstream: functional analysis of genes found in phenotype-association tests (to be added)
* Annotation: variant calling results annotation

### R libraries used:
dplyr
castor
ape

### Python libraries used:
ete3 
biopython 
pandas

### Annotation
All annotation scripts can be used independently. Annotation.py is also a part of the pipeline decribed in https://github.com/Reshetnikoff/m.tuberculosis-research-code. Annotation.py is a preliminary data processing script. Nucleotide sequences of each protein-coding gene is translated and aligned to the reference protein sequences with extreme care to take into account possible frameshifts, start or stop codon loss, etc. 
The script takes 3 arguments as input: (1) a file with nucleotide variants to annotate, (2) folder to write an output file, (3) threshold of mutated gene ORF length: gene is considered to be "broken" if the length of its ORF is below the threshold (0.5 is recommended). The script needs H37Rv reference fasta sequence, CDS annotation and non-coding regions annotation (all these files can be found in Annotation/data, here Mycobrowser v.3 annotation is used). Also, the script needs blast to be installed and available through the command line. The script creates isolate_name_result.tsv file which contains information about all mutations relative to genes. An example of the output can be found in Examples/ folder. 

**Example of use:**
```
python3 Annotation.py SAMN03649057.variants test/ 0.5
```

* Here SAMN03648483.variants is a file with all nucleotide replacements (including insertions and deletions) in an isolate SAMN03648483 genome relative to H37Rv reference genome (all ".variants" files used in the analysis can be found in https://github.com/Reshetnikoff/m.tuberculosis-research-code/tree/main/db/nucl_data.tar.gz),
* test/ is an output folder,
* 0.5 is a "broken gene" threshold. 

**What is the script doing**

For an isolate we take the nucleotide sequence of each protein-coding gene along with its 180 bp upstream and downstream regions which differs from the reference sequence and run blastx to compare this nucleotide sequence with the reference protein sequence of a given gene. Next, the first start codon upstream of the first position of the best blastx hit among all reading frames is found (including alternative start codons, GTG and TTG). The isolate’s nucleotide sequence is translated from this start codon to the first stop codon within the same reading frame. The protein sequence obtained is aligned to the reference protein sequence with Needleman–Wunsch algorithm. The gene is labeled as “broken” if it is shorter than {threshold} the reference protein length or 30% longer than the reference protein. If the gene is broken, no mutation is considered in this gene. Otherwise, all mutations inside the gene are annotated according to position in protein based on the alignment. SNVs, insertions and deletions at a distance less than 100 bp upstream of a gene are annotated as upstream relative to this gene. Nucleotide variants in noncoding genes are annotated using the reference genome gene boundaries. 

**Output file description:**

Annotation file contains six tab-separated columns: gene name, a mutation position according to corresponding protein sequence, reference amino acid, alternative amino acid, mutation type and codon column (for snps). Mutation type can be either 'snp' for single amino acid substitution, 'syn' for synonymous substitutions, 'ins' for an insertion and 'del' for a deletion. For insertions there is a dash ('-') symbol for reference amino acid sequence, for deletions there is no reference or alternative amino acids, just a deletion length (in amino acids) instead. Also, there can be additional labels in annotation such as 'right_extension' or 'left_extension' if a protein sequence is prolonged, and 'right_clip' and 'left_clip' if a protein sequence is shortened due to start or stop codon mutation. If a protein sequence is shortened by more than a {threshold} or prolonged by more than 30%, the gene is considered to be 'broken'.

Please note that for the analysis performed in https://github.com/Reshetnikoff/m.tuberculosis-research-code during the matrices preparation synonymous substitutions are ignored, the codon column is also ignored, as well as 'right_extension', 'right_clip', 'left_extension' and 'left_clip' features.




