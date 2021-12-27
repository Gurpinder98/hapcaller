import numpy as np

def gff_parse(gff_file, genes_to_look):
    """
    Reads GFF file and returns a dictionary of following structure:
    {
        contig: {
            gene : (start_position, stop position)
            ..
        }
        ..
    }
    """
    genes = {}
    with open(gff_file, "r") as in_f:
        for line in in_f:
            line = line.rstrip('\n')
            if line.startswith('#') != True:
                contig = line.split('\t')[0]
                locus_type = line.split('\t')[2]
                if locus_type ==  'gene':
                    gene_name = line.split("\t")[8].split(";")[1].lstrip('Name=')
                    if gene_name in genes_to_look:
                        if contig not in genes.keys():
                            genes[contig] = {}
                        (start_pos, stop_pos) = (line.split("\t")[3], line.split("\t")[4])
                        genes[contig][gene_name] = (start_pos, stop_pos)    
    return genes


def vcf_parse(vcf_file, genes):
    """
    reads the vcf file and outputs a dictionary of the following structure:

    complete_allele_patterns: {
        sample: {
            contig: ('pattern1', 'pattern2', [information about each position])
        }
    }
    """
    # complete_allele_patterns: {SLXXXX:{LKXXXX:{'0001','10101'}, LKXX:..}, SLXX..}
    complete_allele_patterns = {}


    with open(vcf_file, 'r') as in_f:
        
        for line in in_f:
            line = line.rstrip('\n')

            # getting information about number of samples from the "#CHROM" line. 
            if line.startswith('#CHROM'):
                samples = line.split('\t')[9:] #9 is specific to the VCF file. 
                print(f"Reading the {vcf_file} file, {len(samples)} samples found.")
                
                for indv in samples:
                    complete_allele_patterns[indv] = {}
            
            if line.startswith('#') != True:
                line_array = line.split('\t')
                locus = line_array[0]
                position = line_array[1]
                ref = line_array[3]
                alt = line_array[4]
                variant_type = line_array[7].split(';')[-1].lstrip('TYPE=')
                
                # if locus is in genes dict passed on from gff parse
                if locus in genes.keys():
                    
                    all_samples_genotypes = []
                    for data in line_array[9:]:
                        G = data.split(':')[0]
                        if G == ".":
                            all_samples_genotypes.append("./.")        
                        else:
                            all_samples_genotypes.append(G)
                    
                    assert len(all_samples_genotypes) == len(complete_allele_patterns.keys())
                    
                    

                    for indv in complete_allele_patterns:
                        complete_allele_patterns[indv][locus] = ['','']
                    for i in range(len(all_samples_genotypes)):
                        
                        complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][0] = complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][0] + all_samples_genotypes[i].split("/")[0]
                        complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][1] = complete_allele_patterns[list(complete_allele_patterns.keys())[i]][locus][1] + all_samples_genotypes[i].split("/")[1]


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