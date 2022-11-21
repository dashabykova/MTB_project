#!/usr/bin/env python
# coding: utf-8

import sys
import pandas as pd
from Bio import Seq, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
import numpy as np
import os
import subprocess
import warnings
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

def first_start_from_left(seq):
    i = len(seq)
    while i >= 3:
        if (seq[i - 3:i] == 'ATG') | (seq[i - 3:i] == 'GTG') | (seq[i - 3:i] == 'TTG'): 
            return i - 3
        else:
            i -= 3
    return -1

def get_codon(protein_pos, nuc_seq):
    return nuc_seq[protein_pos*3 - 3:protein_pos*3]

def analyze_alignment(alignment, alt_seq, ref_seq):
    protein_ref = alignment[0].seqA
    protein_alt = alignment[0].seqB
    mutation_dict = {}
    al_start = alignment[0].start
    al_end = alignment[0].end
    indel_length = 0
    identity = 0
    gap_number_r = protein_ref[:al_start].count('-')
    gap_number_a = protein_alt[:al_start].count('-')
    for i in range(al_start, al_end):
        if (protein_alt[i] == '-'):
            indel_length += 1
            gap_number_a += 1
            if indel_length == 1:
                indel_start = i
            mtype = 'del'
        elif (protein_ref[i] == '-'):
            gap_number_r += 1
            indel_length += 1
            if indel_length == 1:
                indel_start = i - 1
            mtype = 'ins'
        else:
            if indel_length > 0:
                mutation_dict[indel_start] = (mtype, indel_length)
                indel_length = 0
            alt_codon = get_codon(i - gap_number_a + 1, alt_seq)
            ref_codon = get_codon(i - gap_number_r + 1, ref_seq)
            if protein_ref[i] != protein_alt[i]:
                mutation_dict[i] = ('snp', alt_codon)
            else:
                identity += 1
                if alt_codon != ref_codon:
                    mutation_dict[i] = ('syn', alt_codon)
    return identity/(al_end - al_start + 1), mutation_dict 

def translate_and_align(translated_old, query, blast_start, blast_end, frame, ref_seq, break_threshold=0.7):
    if frame < 0:
        myseq = Seq.reverse_complement(query)
        true_start = first_start_from_left(myseq[:len(query) - blast_end + 3])
    else:
        myseq = query
        true_start = first_start_from_left(myseq[:blast_start + 2])
    if true_start != -1:
        translated_new = Seq.translate(myseq[true_start:], table='Bacterial', to_stop=True)
        new_length = len(translated_new) * 3 + 3
        if (len(translated_new) <= len(translated_old) * break_threshold) | (len(translated_new) >= len(translated_old) * 1.3):
            return 'broken'
        else:
            alignment = pairwise2.align.localms(translated_old, translated_new, 4, -2, -7, -1)
            old_length = len(translated_old) * 3 + 3
            alt_seq = myseq[true_start:true_start + new_length]
            if frame < 0:
                ref_seq = Seq.reverse_complement(ref_seq)
            identity, mutation_dict = analyze_alignment(alignment, alt_seq, ref_seq)
            if identity < break_threshold:
                return 'broken'
            return [alignment, mutation_dict, alt_seq]
    else:
        return 'broken'
        
def get_mtype(ref, alt):
    if ref.isalpha() & (len(ref) == 1) & alt.isalpha() & (len(alt) == 1):
        return 'snp'
    elif (len(ref) < len(alt)) | (ref == '-'):
        return 'ins'
    else:
        return 'del'
    
def analyze_mutations(value, alt_pos, alignment, frame, pos_mutation, mutation_dict, 
                      output, alt_seq, ref_seq, altered_positions, upstream=100):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    name, start, end = value
    al_start = alignment[0].start
    al_end = alignment[0].end
    protein_ref = alignment[0].seqA
    protein_alt = alignment[0].seqB
    left, right, left_ext = False, False, False
    np_alt_pos = np.asarray(alt_pos)
    #process alignment (nonsynonymous mutations inside protein)
    if frame > 0:
        np_alt_pos = np_alt_pos - start
    else:
        np_alt_pos = end - np_alt_pos
    if (al_start != 0) & (np_alt_pos[(np_alt_pos < al_start * 3 + upstream) & (np_alt_pos > -al_start * 3)].size > 0):
        left_ref = protein_ref[:al_start]
        left_alt = protein_alt[:al_start]
        left = True
        if left_ref[0] == '-':
            left_ext = True
            mtype = 'left_extension'
            output.write('\t'.join([name, str(al_start),
                                     left_alt, mtype]) + '\n')
        else:
            mtype = 'left_clip'
            left_clip_length = len(left_ref)
            output.write('\t'.join([name, str(al_start), mtype]) + '\n')
    if (al_end < len(alignment[0].seqA)):
        right_ref = protein_ref[al_end:]
        right_alt = protein_alt[al_end:]
        right = True
        if right_ref[-1] == '-':
            mtype = 'right_extension'
            output.write('\t'.join([name, str(al_end + 1 - protein_ref[:al_end].count('-')),
                                     right_alt, mtype]) + '\n')
        else:
            mtype = 'right_clip'
            right_clip_length = len(right_ref)
            output.write('\t'.join([name, str(al_end + 1 - protein_ref[:al_end].count('-')), 
                         str(right_clip_length), mtype]) + '\n')
    for m in mutation_dict:
        pos = m - protein_ref[:m].count('-') + 1
        if mutation_dict[m][0] == 'del':
            output.write('\t'.join([name, str(pos), str(mutation_dict[m][1]), 'del']) + '\n')
        elif mutation_dict[m][0] == 'ins':
            output.write('\t'.join([name, str(pos), '-', protein_alt[m:m + mutation_dict[m][1]], 
                                    'ins']) + '\n')
        else:
            mtype, codon = mutation_dict[m]
            output.write('\t'.join([name, str(pos), protein_ref[m], protein_alt[m], mtype, codon]) + '\n')
    #process non-coding and synonymous mutations
    for p in alt_pos:
        if (p >= start) and (p <= end):
            altered_positions.append(p)
        else:
            if frame > 0:
                from_start =  p - start - 1
            else:
                from_start = end - p
            if (from_start < 0) and (from_start > -upstream):
                ref = pos_mutation[p-1][0]
                alt = pos_mutation[p-1][1]
                mtype = get_mtype(ref, alt)
                output.write('\t'.join([name, str(from_start), ref, 
                                            alt, mtype]) + '\n')
                altered_positions.append(p)

if __name__ == "__main__":
    #reference sequence
    nuc_ref = SeqIO.read('data/h37rv.fasta', 'fasta')
    #MTB coding genes annotation
    cds = pd.read_csv('data/h37rv_new_proteins.tsv', sep='\t')
    #file with variant calling result for an isolate (the list of variants)
    input_f = sys.argv[1]
    #output folder
    outfolder = sys.argv[2]
    #threshold for identity
    threshold = float(sys.argv[3])
    #folder for blastx temporary files
    blast_folder = 'blastx_tempfiles/'
    os.makedirs(blast_folder, exist_ok=True)
    organism = input_f.split('/')[-1][:-9]
    outname = outfolder + organism + '_result.tsv'
    nonref = pd.read_csv(input_f, sep='\t', header=None, names=['pos', 'ref', 'alt'])
    #mutations dictionary
    positions = np.asarray(nonref.pos)
    pos_mutation = {}
    for i in range(len(nonref)):
        mutation = nonref.iloc[i]
        pos = mutation['pos']
        alt = mutation['alt']
        ref = mutation['ref']
        pos_mutation[pos - 1] = [ref, alt]

    cds.dropna(subset=['protein'], inplace=True)

    flank = 180
    nonref_dict = {}
    #get nonreference nucleotide sequence with flanking regions and reference protein sequence for each protein coding gene
    altered_positions = []
    for i in range(len(cds)):
        seq = cds.iloc[i]
        name = seq['name']
        gene_synonym = seq['gene_synonym']
        start = seq['start']
        end = seq['end']
        ref_seq_p = seq['protein'].replace(' ', '')
        ref_seq = str(nuc_ref.seq[start:end])
        strand = seq['strand']
        locus_tag = seq['locus_tag']
        changed = list(positions[(positions >= max(start - flank, 0)) & (positions <= min(end + flank, len(nuc_ref.seq)))])
        if changed:
            changed.sort()
            #altered_positions.extend(changed)
            nonref_seq = str(nuc_ref.seq[max(start - flank, 0):changed[0] - 1]) + pos_mutation[changed[0] - 1][1]
            len_mut = len(pos_mutation[changed[0] - 1][0])
            for i in range(1, len(changed)):
                nonref_seq = nonref_seq + str(nuc_ref.seq[changed[i - 1] + len_mut - 1:changed[i] - 1]) + pos_mutation[changed[i] - 1][1]
                len_mut = len(pos_mutation[changed[i] - 1][0])
            nonref_seq = nonref_seq + str(nuc_ref.seq[changed[-1] + len_mut - 1:min(end + flank, len(nuc_ref.seq))])
            if str(name) == 'nan':
                if str(gene_synonym) == 'nan':
                    name = locus_tag
                else:
                    name = gene_synonym
            nonref_dict[(name, start, end)] = (changed, nonref_seq, ref_seq_p, ref_seq)

    #working with nonref_dict, doing blastx (nonreference nucleotide vs reference protein sequence) and processing result
    blast_temp = blast_folder + str(threshold) + organism
    for value in nonref_dict:
        name, start, end = value
        with open(blast_temp + 'query.fasta', 'w') as temp:
            temp.write('>query\n')
            query = nonref_dict[value][1]
            temp.write(query)
        with open(blast_temp + 'subject.fasta', 'w') as subj_temp:
            subj_temp.write('>subject\n')
            subj = nonref_dict[value][2]
            subj_temp.write(subj)
        subprocess.run(['blastx','-out', blast_temp + 'res.xml', '-outfmt', '5', '-query', blast_temp + 'query.fasta', 
                        '-evalue', '0.001', '-subject', blast_temp + 'subject.fasta', '-word_size', '2',
                       '-query_gencode', '11'])
        with open(blast_temp + "res.xml") as result_handle:
            blast_record = NCBIXML.read(result_handle)
            if blast_record.alignments:
                best_al = blast_record.alignments[0]
            else:
                with open(outname, 'a') as broken:   
                    broken.write(name + '\tbroken\n')
                continue
            blast_start = best_al.hsps[0].query_start
            blast_end = best_al.hsps[0].query_end
            frame = best_al.hsps[0].frame[0]
        ref_seq = nonref_dict[value][3]
        translation_res = translate_and_align(subj, query, blast_start, blast_end, frame, ref_seq, break_threshold=threshold)
        if translation_res == 'broken':
            with open(outname, 'a') as broken:
                broken.write(name + '\tbroken\n')
        else:
            alignment, mutation_dict, alt_nuc_seq = translation_res
            alt_pos = nonref_dict[value][0]
            output = open(outname, 'a')
            analyze_mutations(value, alt_pos, alignment, frame, pos_mutation, 
                              mutation_dict, output, alt_nuc_seq, ref_seq, altered_positions)
            output.close()
    os.remove(blast_temp + 'subject.fasta')
    os.remove(blast_temp + 'query.fasta')
    os.remove(blast_temp + 'res.xml')

    output.close()

    #non-coding genes
    output = open(outname, 'a')
    #non-coding genes annotated
    nc_genes = pd.read_csv('data/h37rv_noncoding.tsv', sep='\t')
    for i in range(len(nc_genes)):
        seq = nc_genes.iloc[i]
        gene_type = seq['type']
        name = seq['name'] 
        gene_synonym = seq['gene_synonym']
        start = seq['start']
        end = seq['end']
        strand = seq['strand']
        locus_tag = seq['locus_tag']
        changed = list(positions[(positions >= start) & (positions <= end)])
        if changed:
            altered_positions.extend(changed)
            if str(name) == 'nan':
                if str(gene_synonym) == 'nan':
                    name = locus_tag
                else:
                    name = gene_synonym
            for i in changed:
                if strand == 1:
                    from_start = i - start + 1
                else:
                    from_start = end - i + 1
                mtype = get_mtype(pos_mutation[i - 1][0], pos_mutation[i - 1][1])
                output.write('\t'.join([name, str(from_start), pos_mutation[i - 1][0], 
                                        pos_mutation[i - 1][1], mtype]) + '\n')
    output.close()
    
    #all variants left
    for key in pos_mutation:
        if key + 1 not in altered_positions:
            with open(outname, 'a') as output:
                mtype = get_mtype(pos_mutation[key][0], pos_mutation[key][1])
                output.write('\t'.join(['-', str(key), pos_mutation[key][0], 
                                        pos_mutation[key][1], mtype]) + '\n')

