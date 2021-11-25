# -*- coding: utf-8 -*-

import numpy as np
from collections import OrderedDict
#import argparse
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import os

GFF_File = 'Brassica_napus.AST_PRJEB5043_v1.44.sorted.gff3'
VCF_FILE = 'raw.g5mac3dplm.recode.vcf'
SAMPLE_NAMES_FILE = 'sample_names.txt'

GENE_TO_LOOK = ['BnaA02g00370D', 'BnaA03g02820D', 'BnaA03g13630D', 'BnaA10g22080D', 'BnaC02g00490D', 'BnaC03g04170D', 'BnaC03g16530D', 'BnaC09g46500D', 'BnaC09g46540D', 'BnaA02g06490D', 'BnaA04g15260D', 'BnaC04g53290D', 'BnaC04g53290D', 'BnaC04g53290D', 'BnaA05g05010D', 'BnaCnng45490D', 'BnaAnng19140D', 'BnaA03g37880D', 'BnaA03g39820D', 'BnaA06g24000D', 'BnaA04g13710D']
#GENE_TO_LOOK = ['BnaC02g00490D']

def gff_parse(gff_file, genes_to_look):
    
    Genes = {}
    genes_found = 0
    with open(gff_file, "r") as in_f:
        line = in_f.readline()
        while line:
            line = line.rstrip('\n')
            if line.startswith('#') != True:
                contig = line.split('\t')[0]
                locus_type = line.split('\t')[2]
                if locus_type ==  'gene':
                    gene_name = line.split("\t")[8].split(";")[1].lstrip('Name=')
                    if gene_name in genes_to_look:
                        if contig not in Genes.keys():
                            Genes[contig] = {}
                        (start_pos, stop_pos) = (line.split("\t")[3], line.split("\t")[4])
                        Genes[contig][gene_name] = (start_pos, stop_pos)
                        genes_found += 1
            line = in_f.readline()
    
    return Genes, genes_found


def vcf_parse(vcf_file, Genes):
    # complete_allele_patterns: {SLXXXX:{LKXXXX:{'0001','10101'}, LKXX:..}, SLXX..}
    complete_allele_patterns = {}
    # super_contigs: { LKXXXX: OrderedDict{POS1: (ref, alt), POS2: (ref, alt)..}, LKXXX..}
    #order of posttion values in dict is crucial for this - so LKXXX dicts are all ordered.  
    super_contigs = {}

    with open(vcf_file, 'r') as in_f:
        line = in_f.readline()
        while line:
            line = line.rstrip('\n')
            if line.startswith('#CHROM'):
                samples = line.split('\t')[9:]
                print("Reading the {} file, {} samples found.".format(vcf_file, len(samples)))
                for indv in samples:
                    complete_allele_patterns[indv] = {}
            if line.startswith("#") != True:
                locus = line.split('\t')[0]
                if locus in list(Genes.keys()):
                    super_contigs[locus] = OrderedDict()
                    for indv in complete_allele_patterns.keys():
                        complete_allele_patterns[indv][locus] = ['','']
            line = in_f.readline()

    with open(vcf_file, 'r') as in_f:
        line = in_f.readline()
        while line:
            line = line.rstrip('\n')
            if line.startswith('#') != True:
                line_array = line.split('\t')
                locus = line_array[0]
                position = line_array[1]
                ref = line_array[3]
                alt = line_array[4]
                variant_type = line_array[7].split(';')[-1].lstrip('TYPE=')
                if locus in Genes.keys():
                    super_contigs[locus][position] = (position, variant_type, ref, alt)
                    
                    #for genotype
                    all_samples_genotypes = []
                    for data in line_array[9:]:
                        G = data.split(':')[0]
                        if G == ".":
                            all_samples_genotypes.append("./.")        
                        else:
                            all_samples_genotypes.append(G)
                    
                    assert len(all_samples_genotypes) == len(complete_allele_patterns.keys())
                    for i in range(len(all_samples_genotypes)):
                        
                        complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][0] = complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][0] + all_samples_genotypes[i].split("/")[0]
                        complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][1] = complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][1] + all_samples_genotypes[i].split("/")[1]
                        #else:
                        #    print("ERROR in Genotype representation for {}:{} {}".format(locus, position, all_samples_genotypes[i]))
            line = in_f.readline()


    return complete_allele_patterns, super_contigs


def alleles_slicer(start_stop, allele_seqs, positions, gene, locus, sample):
        
    """
        Takes whole contig allele sequences and outputs sliced sequences.
        
        input:
            start_stop: a tuple of ints 
            allele_seqs: a tuple of two alleles patterns. a array/tuple of strings.
            positions: a list of positions for every allele in haplotype sequence.
        
        returns:
            a tupple of sliced allele seqs
    """
    
    try:
        assert len(allele_seqs[0]) == len(positions)
        assert len(allele_seqs[0]) == len(allele_seqs[1])
    except:
        AssertionError
        print("Lengths of position array {} and both allele strings {},{} do not match. {}:{}:{}".format(len(positions),len(allele_seqs[0]), len(allele_seqs[1]), gene,locus,sample))
    pos_array = []
    sliced_allele_0 = ''
    sliced_allele_1 = ''
    for i in range(len(positions)):
        if int(positions[i]) >= int(start_stop[0]) and int(positions[i]) <= int(start_stop[1]):
            pos_array.append(positions[i])
            sliced_allele_0 = sliced_allele_0 + allele_seqs[0][i]    
            sliced_allele_1 = sliced_allele_1 + allele_seqs[1][i]
    
    return (pos_array, [sliced_allele_0, sliced_allele_1])


def final_haplotype_pruning(Genes, complete_allele_patterns ,super_contigs):

    final_data = {}
    for super_locus in super_contigs.keys():
        super_locus_dict = {}
        position_array = list(super_contigs[super_locus].keys())
        for gene in Genes[super_locus].keys():
            current_gene_dict = {}
            start_stop = Genes[super_locus][gene]
            for sample in complete_allele_patterns.keys():
                current_gene_dict[sample] = {} #fill in the innermost dict first
                
                raw_allele_seqs = complete_allele_patterns[sample][super_locus]
                pos_array, sliced_allele_seqs = alleles_slicer(start_stop, raw_allele_seqs, position_array, gene, super_locus, sample)
                
                #marker format : Locus:POS:var_type:Ref allele:alt_allele1,alt_allele2
                marker_array = [super_locus+':'+pos+':'+super_contigs[super_locus][pos][1]+':'+super_contigs[super_locus][pos][2] + ':'+super_contigs[super_locus][pos][3] for pos in pos_array]
                current_gene_dict[sample] = [marker_array, sliced_allele_seqs[0], sliced_allele_seqs[1]]
            final_data[gene] = current_gene_dict

    return final_data


def distance_calculator(seqA_array, seqB_array):
    """
    Return distance between two loci based on two set of allele patterns

    Parameters
    ----------
    seqA_array : List of two strings
        eg. ['101001', '101101']
    seqB_array : List of two strings
        eg. ['101001', '101101']
        

    Returns
    -------
    float
        distance between two samples, averaged. 
        Scoring scheme: +1 if both allele patterns differ
                        +0.5 if only one is different 
                         
    """
    try:
        assert len(seqA_array[0]) == len(seqB_array[0])
        assert len(seqA_array[1]) == len(seqB_array[1])
    except:
        AssertionError
        print("Sequence length mismatch.")
    pattern0 = zip(seqA_array[0], seqB_array[0])
    difference0 = sum([1 for allele in pattern0 if allele[0] != allele[1]])
    pattern1 = zip(seqA_array[1], seqB_array[1])
    difference1 = sum([1 for allele in pattern1 if allele[0] != allele[1]])

    return (difference0 + difference1)/2


def distance_calculator_eu(seqA_array, seqB_array):
    """
    Return distance between two loci based on two set of allele patterns

    Parameters
    ----------
    seqA_array : List of two strings
        eg. ['101001', '101101']
    seqB_array : List of two strings
        eg. ['101001', '101101']
        

    Returns
    -------
    float
        distance between two samples, averaged. 
        Scoring scheme: +1 if both allele patterns differ
                        +0.5 if only one is different 
                         
    """
    try:
        assert len(seqA_array[0]) == len(seqB_array[0])
        assert len(seqA_array[1]) == len(seqB_array[1])
    except:
        AssertionError
        print("Sequence length mismatch.")
    vectorA_1 = np.array([int(i) for i in seqA_array[0]], dtype=int)
    vectorA_2 = np.array([int(i) for i in seqA_array[1]], dtype=int)

    vectorB_1 = np.array([int(i) for i in seqB_array[0]], dtype=int)
    vectorB_2 = np.array([int(i) for i in seqB_array[1]], dtype=int)
    
    eu_dist_1 = int(np.sqrt(np.sum(np.square(vectorA_1 - vectorB_1))))
    eu_dist_2 = int(np.sqrt(np.sum(np.square(vectorA_2 - vectorB_2))))
    
    return eu_dist_1 + eu_dist_2


def supplementary_output(final_data, gene, sample_names_dict, output_file_name = None):
    if output_file_name == None:
        output_file_name = gene+"_haps.csv"
    else:
        output_file_name = os.path.join(output_file_name, gene+"_haps.csv")
        
    with open(output_file_name, "w") as out_f:
        #first line
        contig_line = " , ," + ','.join([s.split(":")[0] for s in final_data[gene][list(final_data[gene].keys())[0]][0]]) + "\n"
        zeroth_line = "POS, ,"+ ','.join([s.split(":")[1] for s in final_data[gene][list(final_data[gene].keys())[0]][0]]) + "\n"
        first_line = 'TYPE, ,'+ ','.join([s.split(":")[2].replace(',',';') for s in final_data[gene][list(final_data[gene].keys())[0]][0]]) + "\n"
        second_line = 'REF, ,'+ ','.join([s.split(":")[3] for s in final_data[gene][list(final_data[gene].keys())[0]][0]]) + "\n\n"
        out_f.write(contig_line+zeroth_line+first_line+second_line)
        for sample in final_data[gene]:
            haplotype1 = final_data[gene][sample][1]
            haplotype2 = final_data[gene][sample][2]
            seq1 = []
            seq2 = []
            for i in range(len(haplotype1)):
                ref = [final_data[gene][sample][0][i].split(':')[3]]
                alt = final_data[gene][sample][0][i].split(':')[4].split(",")
                ref_alt = ref + alt
                    
                if haplotype1[i] != '.' and haplotype2[i] != '.':
                    seq1.append(ref_alt[int(haplotype1[i])]) if int(haplotype1[i]) != 0 else seq1.append('*')
                    if haplotype1[i] == haplotype2[i]:
                        seq2.append('-')
                    else:
                        seq2.append(ref_alt[int(haplotype2[i])]) if int(haplotype2[i]) != 0 else seq2.append('*')
                        
                else:
                    if haplotype1[i] == '.':
                        seq1.append('.')
                        if haplotype2[i] == '.':
                            seq2.append('.')
                        else:
                            seq2.append(ref_alt[int(haplotype2[i])]) if int(haplotype2[i]) != 0 else seq2.append('*')
                    else:
                        seq1.append(ref_alt[int(haplotype1[i])]) if int(haplotype1[i]) != 0 else seq1.append('*')
                        if haplotype2[i] == '.':
                            seq2.append('.')
                        else:
                            seq2.append(ref_alt[int(haplotype2[i])]) if int(haplotype2[i]) != 0 else seq2.append('*')
            
            sample_line = sample_names_dict[sample]+","+ sample+","+','.join(seq1)+"\n"+' , '+','+','.join(seq2)+"\n\n"
            
            out_f.write(sample_line)        
                    
def create_dendrogram(linkage_matrix, gene, sample_line_names, output_file_name=None):
    
    if output_file_name == None:
        output_file_name = gene+".jpg"
    else:
        output_file_name = os.path.join(output_file_name, gene+".jpg")
        
    
    plt.figure(figsize=(15,7))
    plt.box(False)
    dendrogram(linkage_matrix, labels=sample_line_names)
    plt.yticks(np.arange(0,max(plt.yticks()[0]), 1))                    
    plt.grid(color='lightgrey', which='major', axis='y', linestyle='-')
    plt.title(gene)
    plt.savefig(output_file_name, dpi=300)
    
#def cut_off_input_file()


    
def main():
    Genes, num_genes_found = gff_parse(GFF_File, GENE_TO_LOOK)
    print("File successfully read, {} out of {} genes found.".format(num_genes_found, len(GENE_TO_LOOK)))
    complete_haplotypes, variant_data = vcf_parse(VCF_FILE, Genes)
    print("Data from VCF file loaded.")
    final_data = final_haplotype_pruning(Genes, complete_haplotypes, variant_data)

    samples = list(final_data[GENE_TO_LOOK[0]].keys())
    with open(SAMPLE_NAMES_FILE, "r") as in_f:
            lines = in_f.readlines()

    sample_names_dict = {}
    for line in lines[1:]:
        sample_names_dict[line.split("\t")[0]] = line.split("\t")[1].replace("\n",'')
    sample_line_names = [sample_names_dict[s] for s in samples]
    
    for bna_gene in final_data.keys():
        print("Working on {}".format(bna_gene))
        samples = list(final_data[bna_gene].keys())
        Distance_Matrix = np.zeros((len(samples), len(samples)))

        for vertical in range(Distance_Matrix.shape[0]):
            for horizontal in range(Distance_Matrix.shape[1]):
                Distance_Matrix[vertical][horizontal] = distance_calculator(final_data[bna_gene][samples[vertical]][1:],final_data[bna_gene][samples[horizontal]][1:])

        dists = squareform(Distance_Matrix)
        linkage_matrix = linkage(dists, "single")
        
        
        
        if 'Haplotypes' not in os.listdir():
            os.mkdir('Haplotypes')
        if bna_gene not in os.listdir('Haplotypes'):
            os.mkdir(os.path.join('Haplotypes', bna_gene))
            
        supplementary_output(final_data, bna_gene,sample_names_dict, os.path.join('Haplotypes', bna_gene))
        create_dendrogram(linkage_matrix, bna_gene, sample_line_names, os.path.join('Haplotypes', bna_gene))
        print("Finished gene {}\n".format(bna_gene))
        
    print("Dendrograms created. Provide the desired cut off values.\n")
    input("Program paused. Press any key to continue..")
        
        
        
    return final_data, Genes, linkage_matrix, sample_line_names, Distance_Matrix
