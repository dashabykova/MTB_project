## Epistatic interactions in Mycobacterium tuberculosis genes 
Here will be scripts for the research project

- Phylogeny: ancestor states reconstruction and reconstruction results parsing
- Downstream: functional analysis of genes found in phenotype-association tests (to be added)
- Annotation: variant calling results annotation

### R libraries used:
dplyr
castor
ape

### Python libraries used:
ete3 
biopython 
pandas

### Annotation
All annotation scipts can be used independently. They are also the part of the pipeline decribed in https://github.com/Reshetnikoff/m.tuberculosis-research-code. Nucleotide variants to annotate are stored in https://github.com/Reshetnikoff/m.tuberculosis-research-code/tree/main/db/nucl_data.tar.gz. To obtain the annotation for a particular isolate, one has to run:

