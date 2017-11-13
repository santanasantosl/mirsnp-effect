from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import os
import re
import subprocess
import random
import multiprocessing as mp
import math
import sys


def filter_vcf(parsed_mirdsnp, reference_genome_dir, outfile):

    all_snps = list()
    for mirna in parsed_mirdsnp.keys():
        all_snps = all_snps + parsed_mirdsnp[mirna]

    all_snps = list(set(all_snps))
    dbsnp_vcf = os.path.join(reference_genome_dir, 'dbsnp_138.hg19.vcf')
    dbsnp_vcf_handle = open(dbsnp_vcf, 'r')

    filtered_vcf = os.path.join(reference_genome_dir, outfile)
    filtered_vcf_handle = open(filtered_vcf, 'w')

    for line in dbsnp_vcf_handle:
        line = line.strip()
        if not re.search('^#', line):
            fields = line.split('\t')
            rs_number = fields[2]
            if rs_number in all_snps:
                filtered_vcf_handle.write(line + '\n')
    filtered_vcf_handle.close()


def generate_fasta_files(annotated_variants, vcf_list, reference_genome_dict , kmer_length):

    reference_sequences = dict()
    mutated_sequences = dict()
    alleles = dict()

    for line in vcf_list:
        if not re.search('^#', line):
            #print line
            fields = line.split('\t')
            chr = fields[0]
            pos = fields[1]
            rs_num = fields[2]
            ref = fields[3]
            alt = fields[4]

            alleles[rs_num] = dict()
            alleles[rs_num]['ref'] = ref
            alleles[rs_num]['alt'] = alt

            #print annotated_variants.keys()
            if annotated_variants.has_key(rs_num):
                snp_strand = annotated_variants[rs_num][0]
                gene_strand = annotated_variants[rs_num][3]

                if len(ref) == 1 and len(alt) == 1:
                    if reference_genome_dict.has_key(chr):
                            reference_genome_seq = str(reference_genome_dict[chr].seq)
                            #reference_genome_seq = Seq(reference_genome_seq, generic_dna)
                            #reference_genome_seq = str(reference_genome_seq.transcribe().reverse_complement())
                            start_pos = int(pos) - kmer_length
                            end_pos = int(pos) + kmer_length
                            #print rs_num + '\t' + str(start_pos) + '\t' + str(end_pos)
                            reference_sequence = reference_genome_seq[start_pos - 1 :end_pos]
                            reference_sequence = reference_sequence[0:kmer_length] + ref + reference_sequence[kmer_length + 1:]
                            mutated_sequence = ""

                            ####### RSNUM AND GENE RELATIONSHIP IS NEEDED ###########
                            if snp_strand == '+':
                                mutated_sequence = reference_sequence[0:kmer_length] + alt + reference_sequence[kmer_length + 1:]

                                if gene_strand == '-':
                                    mutated_sequence = Seq(mutated_sequence, generic_dna)
                                    mutated_sequence = str(mutated_sequence.transcribe().reverse_complement())

                                    reference_sequence = Seq(reference_sequence, generic_dna)
                                    reference_sequence = str(reference_sequence.transcribe().reverse_complement())

                            elif snp_strand == '-':
                                reference_sequence = Seq(reference_sequence, generic_dna)
                                reference_sequence = str(reference_sequence.transcribe().reverse_complement())
                                mutated_sequence = reference_sequence[0:kmer_length] + alt + reference_sequence[
                                                                                             kmer_length + 1:]

                                if gene_strand == '+':
                                    mutated_sequence = Seq(mutated_sequence, generic_dna)
                                    mutated_sequence = str(mutated_sequence.transcribe().reverse_complement())

                                    reference_sequence = Seq(reference_sequence, generic_dna)
                                    reference_sequence = str(reference_sequence.transcribe().reverse_complement())

                            for i in range(1, len(reference_sequence) - kmer_length):
                                current_reference_sequence = reference_sequence[i:i+kmer_length]
                                if len(current_reference_sequence) == kmer_length:

                                    if not reference_sequences.has_key(rs_num):
                                        reference_sequences[rs_num] = list()
                                    reference_sequences[rs_num].append(current_reference_sequence)

                            for i in range(1, len(mutated_sequence) - kmer_length):
                                current_mutated_sequence = mutated_sequence[i:i+kmer_length]
                                if len(current_mutated_sequence) == kmer_length:
                                    if not mutated_sequences.has_key(rs_num):
                                        mutated_sequences[rs_num] = list()
                                    mutated_sequences[rs_num].append(current_mutated_sequence)
                else:
                    print 'NON FATAL ERROR: variant ' + rs_num + ' will be excluded from analysis.'
                    print 'Reference or alternative alleles have length longer that 1.'

    return reference_sequences, mutated_sequences, alleles


def run_rnahybrid(rnahybrid, target, query):

    random_int = str(random.randint(1, 10000000000))
    target_tmp = "target." + str(random_int) + ".tmp"
    query_tmp = "query." + str(random_int) + ".tmp"
    output_tmp = "output." + str(random_int) + ".tmp"

    position = 0
    mfe = 0

    #rnahybrid_cmd = rnahybrid + " -s 3utr_human " + str(target) + " " + str(query)
    target_tmp_handle = open(target_tmp, "w")
    target_tmp_handle.write('>TARGET\n' + str(target))
    target_tmp_handle.close()

    query_tmp_handle = open(query_tmp, "w")
    query_tmp_handle.write('>QUERY\n' + str(query))
    query_tmp_handle.close()

    rnahybrid_cmd = rnahybrid + " -s 3utr_human -t " + target_tmp + " -q " + query_tmp + " > " + output_tmp
    #print rnahybrid_cmd
    subprocess.call(rnahybrid_cmd, shell=True)

    mfe, position = parse_rnahybrid_out(output_tmp)

    os.remove(target_tmp)
    os.remove(query_tmp)
    os.remove(output_tmp)
    return mfe, position


def get_energy_pvalues(mirna, reference_mfe, alternative_mfe, mirna_length, distribution_dir):

    distribution_dir = os.path.join(distribution_dir, str(mirna_length) + 'mer')

    contents = os.listdir(distribution_dir)

    reference_pvalue = "NA"
    alternative_pvalue = "NA"
    log_ratio = "NA"
    pvalue_of_log_ratio = "NA"
    #energy_difference = "NA"
    #pvalue_energy_difference = "NA"

    background_dist_file_status = 0
    ratio_dist_file_status = 0
    for content in contents:
        #print content
        if re.search("^background.dist." + mirna.replace('-', '_') + '.txt', content):
            background_dist_file_status = 1
        elif re.search("^ratio.dist." + mirna.replace('-', '_') + '.txt', content):
            ratio_dist_file_status = 1

    if background_dist_file_status == 1 and ratio_dist_file_status == 1:
        outfile = "tmp_" + str(random.randint(1, 10000000000)) + ".txt"

        calculate_pvalues_bin = "calculate_binding_energy_pvalues.R "
        get_pvalues_cmd = "Rscript " + calculate_pvalues_bin +\
                          " --mirna=" + str(mirna) +\
                          " --reference=" + str(reference_mfe) +\
                          " --distribution_dir=" + str(distribution_dir) +\
                          " --alternative=" + str(alternative_mfe) +\
                          " --out=" + outfile
                         # " --dist_dir=" + str(distribution_dir)


        subprocess.call(get_pvalues_cmd, shell=True)
        outfile_handle = open(outfile, 'r')
        counter = 1
        for line in outfile_handle:
            line = line.strip()
            if counter > 1:
                fields = line.split('\t')
                reference_pvalue = fields[1]
                alternative_pvalue = fields[3]
                log_ratio = fields[4]
                pvalue_of_log_ratio = fields[5]
                #energy_difference = fields[6]
                #pvalue_energy_difference = fields[7]

            counter += 1
        outfile_handle.close()
        os.remove(outfile)

    #return reference_pvalue, alternative_pvalue, log_ratio, pvalue_of_log_ratio, energy_difference, pvalue_energy_difference
    return reference_pvalue, alternative_pvalue, log_ratio, pvalue_of_log_ratio


def parse_validated_snps(mirsnp_file):

    parsed_mirsnp = dict()

    mirsnp_file_handle = open(mirsnp_file, 'r')
    for line in mirsnp_file_handle:
        line = line.strip()
        fields = line.split('\t')
        mirna = fields[0]
        rsnum = fields[1]
        gene = fields[2]
        if not parsed_mirsnp.has_key(mirna):
            parsed_mirsnp[mirna] = list()
        parsed_mirsnp[mirna].append(rsnum)

    return parsed_mirsnp


def parse_rnahybrid_out(rnahyrbid_out):

    mfe = 0
    position = 0

    rnahyrbid_out_handle  = open(rnahyrbid_out, 'r')
    for line in rnahyrbid_out_handle:
        line = line.strip()
        if re.search('^mfe', line):
            mfe = line.split(':')[1].split(' ')[1]

        elif re.search('^position ', line):
            position = line.split(' ')[2]

    rnahyrbid_out_handle.close()
    return mfe, position


def filter_vcf_given_list_of_snps(list_of_snps, reference_genome_dir, outfile):

    dbsnp_vcf = os.path.join(reference_genome_dir, 'dbsnp_138.hg19.vcf')

    outfile = os.path.join(reference_genome_dir, outfile.replace('.vcf', ''))

    vcftools_bin = "vcftools"
    filter_vcf = vcftools_bin + " --vcf " + dbsnp_vcf + " --snps " + list_of_snps + " --recode --out " + outfile
    subprocess.call(filter_vcf, shell=True)

    os.rename(outfile + '.recode.vcf', outfile + '.vcf')


def run_prediction_snp(rnahybrid, snps, alleles,
                       reference_sequences_dict, mutated_sequences_dict,
                       mirna_names, mirna_dict, distribution_dir):

    output_file = "snp_energy_" + str(random.randint(1, 10000000000)) + ".txt"
    print "Running predictions and saving them to file " + output_file
    output_file_handle = open(output_file, 'w')
    output_file_handle.write('RS_NUM\tMIRNA\tALLELE1\tALLELE1_MFE\tALLELE1_PVALUE\tALLELE2\tALLELE2_MFE\tALLELE2_PVALUE\tPVALUE_RATIO\tRATIO_PVALUE\tPOSITION_ALLELE1\tPOSIITON_ALLELE2\n')

    for mirna in mirna_names:
        for current_snp in snps:

            reference_mfe = 0
            reference_position = 0

            mutated_mfe = 0
            mutated_position = 0

            if mirna_dict.has_key(mirna):
                mirna_seq = str(mirna_dict[mirna].seq)
                mirna_length = len(mirna_seq)

                reference_sequences = reference_sequences_dict[mirna_length]
                mutated_sequences = mutated_sequences_dict[mirna_length]

                ref_allele = 'NA'
                alt_allele = 'NA'

                if alleles[mirna_length].has_key(current_snp):

                    ref_allele = alleles[mirna_length][current_snp]['ref']
                    alt_allele = alleles[mirna_length][current_snp]['alt']

                if reference_sequences.has_key(current_snp):
                    ## REFERENCE SEQUENCES
                    for current_reference_sequence in reference_sequences[current_snp]:
                        current_ref_mfe, current_ref_pos = run_rnahybrid(rnahybrid, current_reference_sequence, mirna_seq)
                        #print current_ref_mfe
                        if float(current_ref_mfe) < float(reference_mfe):
                            reference_mfe = current_ref_mfe
                            reference_position = current_ref_pos

                    ## MUTATED SEQUENCES
                    for current_mutated_sequence in mutated_sequences[current_snp]:
                        current_mut_mfe, current_mut_pos = run_rnahybrid(rnahybrid, current_mutated_sequence, mirna_seq)
                        if float(current_mut_mfe) < float(mutated_mfe):
                            mutated_mfe = current_mut_mfe
                            mutated_position = current_mut_pos

                    reference_pvalue, alternative_pvalue, \
                    log_ratio, pvalue_of_log_ratio = get_energy_pvalues(mirna,
                                                                  reference_mfe,
                                                                  mutated_mfe,
                                                                  mirna_length,
                                                                  distribution_dir)

                    output_file_handle.write(str(current_snp) + '\t' + str(mirna)
                                             + '\t' + str(ref_allele) + '\t' + str(reference_mfe) + '\t' + str(reference_pvalue)
                                             + '\t' + str(alt_allele) + '\t' + str(mutated_mfe) + '\t' + str(alternative_pvalue)
                                             + '\t' + str(log_ratio) + '\t' + str(pvalue_of_log_ratio)
                                             + '\t' + str(reference_position) + '\t' + str(mutated_position) + '\n')

    output_file_handle.close()


def run_prediction_for_list_of_snps_mp(annotated_variants_fn,
                                       list_of_snps_fn, reference_genome_dir, mirna_names_fn,
                                       mirna_fasta, rnahybrid, number_of_threads, distribution_dir, debug=False):

    snps = list()
    annotated_variants = dict()

    list_of_snps_fh = open(list_of_snps_fn, 'r')
    for line in list_of_snps_fh:
        line = line.strip()
        snps.append(line)
    list_of_snps_fh.close()

    mirna_names = list()
    mirna_names_fh = open(mirna_names_fn, 'r')
    for line in mirna_names_fh:
        line = line.strip()
        mirna_names.append(line)
    mirna_names_fh.close()


    ## Load reference genome into dict
    print "Loading genome into memory.."
    reference_genome = os.path.join(reference_genome_dir, 'main_hg19.fasta')
    reference_genome_dict = SeqIO.to_dict(SeqIO.parse(open(reference_genome, 'rU'), 'fasta'))

    ## Load mirna fasta into dict
    print "Loading mirna sequences into memory.."
    mirna_fasta_dict = SeqIO.to_dict(SeqIO.parse(open(mirna_fasta, 'rU'), 'fasta'))

    ## Load annotated variants into memory
    print "Loading variant strand into memory.."
    annotated_variants_handle = open(annotated_variants_fn, 'r')
    for line in annotated_variants_handle:
        line = line.strip()

        if not re.search('^RS_NUM', line):
            fields = line.split('\t')
            snp = fields[0]
            annotated_variants[snp] = fields[1:]

    #print annotated_variants
    annotated_variants_handle.close()

    ## Filter VCF
    random_int = str(random.randint(1, 10000000000))
    random_int = str(3542265467)

    filtered_vcf = "tmp_" + random_int + ".vcf"
    print "Creating filtered VCF file (" + filtered_vcf + ").."
    #filter_vcf_given_list_of_snps(list_of_snps_fn, reference_genome_dir, filtered_vcf)
    filtered_vcf_handle = open(os.path.join(reference_genome_dir, filtered_vcf), 'r')

    print "Loading VCF file into memory.."
    filtered_vcf_list = list()

    for line in filtered_vcf_handle:
        line = line.strip()
        filtered_vcf_list.append(line)

    filtered_vcf_handle.close()

    # Split mirna names into chunks in order to allow multiprocessing
    print "Splitting mirna into chunks.."
    number_of_elements_per_chunk = int(math.ceil(len(mirna_names)/float(number_of_threads)))
    split_mirnas_names_list = [mirna_names[i:i+number_of_elements_per_chunk] for i in xrange(0, len(mirna_names), number_of_elements_per_chunk)]

    print "Generating reference and mutated sequences.."
    reference_sequences = dict()
    mutated_sequences = dict()
    alleles = dict()

    kmer_lengths = [7,8,16,17,18,19,20,21,22,23,24,25,26,27,28]

    for kmer_length in kmer_lengths:
        print "Processing kmer length " + str(kmer_length)
        reference_sequences[kmer_length], \
        mutated_sequences[kmer_length], \
        alleles[kmer_length] = generate_fasta_files(annotated_variants, filtered_vcf_list,
                                                    reference_genome_dict, kmer_length)
    ## Choose the function depending on the debug status
    if not debug:

        processes = [mp.Process(target=run_prediction_snp,
                                args=(rnahybrid, snps, alleles, reference_sequences, mutated_sequences,
                                      split_mirnas_names_list[i],
                                      mirna_fasta_dict, distribution_dir)) for i in range(0, len(split_mirnas_names_list))]
        for p in processes:
            p.start()
        for p in processes:
            p.join()

    else:

        processes = [mp.Process(target=run_prediction_snp_debug,
                                args=(rnahybrid, snps, alleles, reference_sequences, mutated_sequences,
                                      split_mirnas_names_list[i],
                                      mirna_fasta_dict, distribution_dir)) for i in range(0, len(split_mirnas_names_list))]
        for p in processes:
            p.start()
        for p in processes:
            p.join()


def run_prediction_snp_debug(rnahybrid, snps, alleles,
                             reference_sequences_dict, mutated_sequences_dict,
                             mirna_names, mirna_dict, distribution_dir):

    output_file = "snp_energy_" + str(random.randint(1, 10000000000)) + ".txt"
    print "Running predictions and saving them to file " + output_file
    output_file_handle = open(output_file, 'w')
    output_file_handle.write('RS_NUM\tMIRNA\tALLELE1_MFE\tALLELE1_PVALUE\tALLELE2_MFE\tALLELE2_PVALUE\tPVALUE_RATIO\tRATIO_PVALUE\tENERGY_DIFFERENCE\tPVALUE_ENERGY_DIFFERENCE\tPOSITION_ALLELE1\tPOSIITON_ALLELE2\n')

    for mirna in mirna_names:

        for current_snp in snps:

            reference_mfe = 0
            reference_position = 0

            mutated_mfe = 0
            mutated_position = 0

            if mirna_dict.has_key(mirna):

                mirna_seq = mirna_dict[mirna].seq
                mirna_length = len(mirna_seq)

                ref_allele = 'NA'
                alt_allele = 'NA'

                if alleles[mirna_length].has_key(current_snp):
                    ref_allele = alleles[mirna_length][current_snp]['ref']
                    alt_allele = alleles[mirna_length][current_snp]['alt']

                reference_sequences = reference_sequences_dict[mirna_length]
                mutated_sequences = mutated_sequences_dict[mirna_length]

                if reference_sequences.has_key(current_snp):

                    ## REFERENCE SEQUENCES
                    reference_counter = 0
                    reference_min_energy_counter = 0
                    for current_reference_sequence in reference_sequences[current_snp]:
                        current_ref_mfe, \
                        current_ref_pos = run_rnahybrid(rnahybrid,
                                                        current_reference_sequence,
                                                        mirna_seq)

                        print mirna + '\t' + mirna_seq + '\t' + str(current_snp) + '\t' + str(ref_allele) + \
                              '\t' + str(current_ref_mfe) + '\t' + str(reference_counter) + \
                              '\t' + current_reference_sequence + '\t' + 'REF'

                        if float(current_ref_mfe) < float(reference_mfe):
                            reference_mfe = current_ref_mfe
                            reference_position = current_ref_pos
                            reference_min_energy_counter = reference_counter

                        reference_counter += 1

                    current_mut_mfe, \
                    current_mut_pos = run_rnahybrid(rnahybrid,
                                                    mutated_sequences[current_snp][reference_min_energy_counter],
                                                    mirna_seq)

                    print "MINIMUM_ENERGY:" + str(mirna) + '\t' + str(mirna_seq) + '\t' + str(current_snp) + \
                          '\t' + str(reference_mfe) + '\t' + str(reference_min_energy_counter) + '\t' + \
                          str(reference_sequences[current_snp][reference_min_energy_counter]) + '\t' \
                          + str(current_mut_mfe) + '\t' + \
                          mutated_sequences[current_snp][reference_min_energy_counter] + '\t' + 'REF'

                    reference_pvalue, alternative_pvalue, \
                    log_ratio, \
                    pvalue_of_log_ratio = get_energy_pvalues(mirna,
                                                             reference_mfe,
                                                             current_mut_mfe,
                                                             mirna_length,
                                                             distribution_dir)

                    output_file_handle.write(str(current_snp) + '\t' + str(mirna)
                                             + '\t' + str(reference_mfe) + '\t' + str(reference_pvalue)
                                             + '\t' + str(current_mut_mfe) + '\t' + str(alternative_pvalue)
                                             + '\t' + str(log_ratio) + '\t' + str(pvalue_of_log_ratio)
                                             + '\t' + str(current_mut_pos) + '\t' + 'REF' + '\n')

                    ## MUTATED SEQUENCES
                    mutated_counter = 0
                    mutated_min_energy_counter = 0
                    for current_mutated_sequence in mutated_sequences[current_snp]:
                        current_mut_mfe, current_mut_pos = run_rnahybrid(rnahybrid, current_mutated_sequence, mirna_seq)

                        print mirna + '\t' + str(mirna_seq) + '\t' + str(current_snp) + '\t' + str(alt_allele) +\
                              '\t' + str(current_mut_mfe) + '\t' + str(mutated_counter) + '\t' + \
                              current_mutated_sequence + '\t' + 'MUT'

                        if float(current_mut_mfe) < float(mutated_mfe):
                            mutated_mfe = current_mut_mfe
                            mutated_position = current_mut_pos

                        mutated_counter += 1

                    current_ref_mfe, \
                    current_ref_pos = run_rnahybrid(rnahybrid,
                                                    reference_sequences[current_snp][mutated_min_energy_counter],
                                                    mirna_seq)

                    print "MINIMUM_ENERGY:" + str(mirna) + '\t' + str(mirna_seq) + '\t' + str(current_snp) + '\t' + str(
                        mutated_mfe) + '\t' + str(mutated_min_energy_counter) + '\t' + \
                          str(mutated_sequences[current_snp][mutated_min_energy_counter]) + '\t' + \
                          str(current_ref_mfe) + '\t' + reference_sequences[current_snp][mutated_min_energy_counter] + \
                          '\t' + 'MUT'

                    reference_pvalue, alternative_pvalue, \
                    log_ratio, pvalue_of_log_ratio = get_energy_pvalues(
                        mirna, mutated_mfe, current_ref_mfe, mirna_length, distribution_dir)

                    output_file_handle.write(str(current_snp) + '\t' + str(mirna) + '\t'
                                             + str(current_ref_mfe) + '\t' + str(reference_pvalue)
                                             + '\t' + str(mutated_mfe) + '\t' + str(alternative_pvalue)
                                             + '\t' + str(log_ratio) + '\t' + str(pvalue_of_log_ratio)
                                             + '\t' + str(current_ref_pos) + '\t' + 'MUT' + '\n')
    output_file_handle.close()


def run_prediction_for_list_of_snps_sp(rnahybrid, snps,
                                       reference_sequences, mutated_sequences, mirnas, mirna_fasta_dict):
    print "Running predictions..."
    run_prediction_snp(rnahybrid, snps, reference_sequences, mutated_sequences, mirnas, mirna_fasta_dict)


def create_variant_gene_annotation_file(snps_3utr_intersection_file, refseq_id_to_symbol, variants_strand, outfile):

    refseq_id_to_symbol_dict = dict()
    refseq_id_to_strand_dict = dict()

    print "Loading gene information..."
    refseq_id_to_symbol_handle = open(refseq_id_to_symbol, 'r')

    for line in refseq_id_to_symbol_handle:
        line = line.strip()
        fields = line.split('\t')
        refseq = fields[0]
        strand = fields[2]
        gene = fields[3]
        refseq_id_to_symbol_dict[refseq] = gene
        refseq_id_to_strand_dict[refseq] = strand
        #print refseq + '\t' + gene

    refseq_id_to_symbol_handle.close()

    print "Loading variants information..."
    variants_strand_dict = dict()
    variants_strand_handle = open(variants_strand, 'r')

    for line in variants_strand_handle:
        line = line.strip()
        fields = line.split('\t')
        variants_strand_dict[fields[0]] = fields[1]
    variants_strand_handle.close()

    print "Generating output file..."
    outfile_handle = open(outfile, 'w')
    outfile_handle.write('RS_NUM\tVARIANT_STRAND\tREFSEQ_ID\tGENE_SYMBOL\tGENE_STRAND\n')

    snps_3utr_intersection_file_handle = open(snps_3utr_intersection_file, 'r')
    for line in snps_3utr_intersection_file_handle:
        line = line.strip()
        fields = line.split('\t')
        rsid = fields[3]
        annotations = fields[7].split(';')
        for annotation in annotations:
            refseq_id = annotation.split('_')[0] + '_' + annotation.split('_')[1]
            if refseq_id_to_symbol_dict.has_key(refseq_id):
                outfile_handle.write(rsid + '\t' + variants_strand_dict[rsid] + '\t' + refseq_id + '\t' + refseq_id_to_symbol_dict[refseq_id] + '\t' + refseq_id_to_strand_dict[refseq_id] + '\n')
    outfile_handle.close()


def read_input_file(input_fn, annotation_fn):

    # Check if input file exist
    if not os.path.exists(input_fn):
        sys.exit('Input file ' + str(input_fn) + ' does not exist!')

    # Open input file and parse it accordingly
    input_handle = open(input_fn, 'r')

    # Stores rsid and coordinates in different lists
    rsid_list = list()
    coordinates_list = list()

    for line in input_handle:
        line = line.strip()
        # If line begins with rsid stores it in the rsid list
        if re.search('^rs', line):
            rsid_list.append(line)
        # If line correspond to genomic coordinates stores it in the coordinates list
        # Line should be in the format CHR\tPOSITION\tREF\tALT\tSTRAND\tANNOTATION
        # Annotation field is optional
        else:
            fields = line.split('\t')
            if not len(fields) >= 5:
                sys.exit('ERROR: ' + line + ' has an invalid format. Please, check for valid formats in the README file')

            chrom = fields[0]
            pos = fields[1]
            ref = fields[2]
            alt = fields[3]
            strand = fields[4]
            annotation = "."

            if len(fields) == 6:
                annotation = fields[5]

            coordinates_list.append(chrom + '\t' + pos + '\t' + ref + '\t' + alt + '\t' + strand + '\t' + annotation)

    input_handle.close()

    if len(rsid_list) > 0:
        if not os.path.exists(annotation_fn):
            sys.exit('Annotation file ' + str(annotation_fn) + ' does not exist!')

        annotation_handle = open(annotation_fn, 'r')
        for line in annotation_handle:
            line = line.strip()
            fields = line.split('\t')


        annotation_handle.close()


if __name__ == "__main__":

    if len(sys.argv) != 7:
        sys.exit("usage: python generate_background_distribution_ratio.py <RNAHYBRID> <SNP_LIST> <REFERENCE_GENOME_DIR> <MIRNA_LIST> <MIRNA_FASTA> <ANNOTATED_VARIANTS>")

    rnahybrid = sys.argv[1]
    snp_list = sys.argv[2]
    reference_genome_dir = sys.argv[3]
    mirna_list_name = sys.argv[4]
    mirna_fasta = sys.argv[5]
    annotated_variants = sys.argv[6]
    # Number of threads
    number_of_threads = 1
    # Distribution directory location
    distribution_dir = "background_dist"

    run_prediction_for_list_of_snps_mp(annotated_variants,
                                       snp_list,
                                       reference_genome_dir,
                                       mirna_list_name,
                                       mirna_fasta,
                                       rnahybrid,
                                       number_of_threads,
                                       distribution_dir)
