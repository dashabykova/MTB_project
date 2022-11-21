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
All annotation scipts can be used independently. Annotation.py is also a part of the pipeline decribed in https://github.com/Reshetnikoff/m.tuberculosis-research-code. Annotation.py is a preliminary data processing script. Nucleotide sequences of each protein-coding gene is translated and aligned to the reference protein sequences with extreme care to take into account possible frameshifts, start or stop codon loss, etc. 
The script takes 4 arguments as input: (1) file with nucleotide variants to annotate, (2) folder to write an output file, (3) threshold of mutated gene ORF length: gene is considered to be "broken" if the length of its ORF is below the threshold (0.5 is recommended), (4) folder for temporary files to put into. 

Example of use:

python3 Annotation.py SAMN03648483.variants test/ 0.5 blastx_tempfiles/

* Here SAMN03648483.variants is a file with all nucleotide replacements (including insertions and deletions) in an isolate SAMN03648483 genome relative to H37Rv reference genome (all ".variants" files used in the analysis can be found in https://github.com/Reshetnikoff/m.tuberculosis-research-code/tree/main/db/nucl_data.tar.gz),
* test/ is an output folder,
* 0.5 is a "broken gene" threshold,
* blastx_tempfiles/ is a folder to store temporary blastx files, in the end they will be deleted. 




