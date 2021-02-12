from functions.GTtable import *
import pandas as pd
import codecs
import allel
import re

def VCFtoTable(vcf_file):
    # extract data from VCF file
    callset = allel.read_vcf(vcf_file, fields='*', alt_number=1)
    
    ##### VCF DataFrame #####
    vcf_dataframe = pd.DataFrame()
    
    # define fields for dataframe by trying to extract each one
    # if fails, be able to continue
    try:
        samples = callset['samples']
    except:
        print("Could not find SAMPLES field.")
    try:
        chrom = callset['variants/CHROM']
        vcf_dataframe['CHROM'] = chrom
    except:
        print("Could not extract CHROM field.")
    try:
        pos = callset['variants/POS']
        vcf_dataframe['POS'] = pos
    except:
        print("Could not extract POS field.")
    try:
        ref = callset['variants/REF']
        vcf_dataframe['REF'] = ref
    except:
        print("Could not extract REF field.")
    try:
        alt = callset['variants/ALT']
        vcf_dataframe['ALT'] = alt
    except:
        print("Could not extract ALT field.")
    try:
        qual = callset['variants/QUAL']
        vcf_dataframe['QUAL'] = qual
    except:
        print("Could not extract QUAL field.")
    try:
        dp_info = callset['variants/DP']
        vcf_dataframe['DP'] = dp_info
    except:
        print("Could not extract INFO/DP field.")
    try:
        dp_format = callset['calldata/DP']
        dp_df = pd.DataFrame()
        for sample in range(len(samples)):
            column_name = "%s_DP" % samples[sample]
            vcf_dataframe[column_name] = dp_format[:, sample]
    except:
        print("Could not extract FORMAT/DP field.")    
    try:
        gt = allel.GenotypeArray(callset['calldata/GT']).to_gt()
        gt_df = pd.DataFrame()
        for sample in range(len(samples)):
            column_name = "%s_GT" % samples[sample]
            vcf_dataframe[column_name] = gt[:, sample]
            vcf_dataframe[column_name] = vcf_dataframe[column_name].map(lambda gt: codecs.decode(gt, 'UTF-8')) 
    except:
        print("Could not extract GT field")
    try:
        ad = allel.GenotypeArray(callset['calldata/AD']).to_gt()
        ad_df = pd.DataFrame()
        for sample in range(len(samples)):
            column_name = "%s_AD" % samples[sample]
            vcf_dataframe[column_name] = ad[:, sample]
            vcf_dataframe[column_name] = vcf_dataframe[column_name].map(lambda ad: codecs.decode(ad, 'UTF-8'))
    except:
        print("Could not extract AD field.")
    try:
        an = callset['variants/AN']
        vcf_dataframe['AN'] = an
    except:
        print("Could not extract AN field.")
    try:
        typ = callset['variants/TYPE']
        vcf_dataframe['TYPE'] = typ
    except:
        print("Could not extract TYPE field.")
        

    ##### Chromosome Length DataFrame #####    
    # extract chromosome names from vcf_dataframe
    chromosomes = vcf_dataframe.CHROM.unique()

    # extract VCF headers
    headers = allel.read_vcf_headers(vcf_file)
    
    # extract chromosome lengths from 'headers' by using 'contigs' keyword
    contigs = []
    for line in headers[0]:
        if 'contig' in line:
            for chromosome in chromosomes:
                if chromosome in line:
                    contigs.append(line)
                    break
                        
    # create dictionary of chromosome's name and length
    chrom_lengths = {}
    keyword_id = 'ID=(.*?),'
    keyword_len = 'length=(.+?)>'
    for contig in contigs:
        ID = re.search(keyword_id, contig)
        for chromosome in chromosomes:
            if chromosome == ID.group(1):
                length = re.search(keyword_len, contig)
                chrom_lengths[chromosome] = int(length.group(1))
                
    # convert dictionary to dataframe
    chrom_len = pd.DataFrame(chrom_lengths, index=['LEN']).T
    chrom_len.index.name = 'CHROM'
    
    return samples, vcf_dataframe, chrom_len
