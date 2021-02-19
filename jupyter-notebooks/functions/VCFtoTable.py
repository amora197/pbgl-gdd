from functions.AttributeGuide import *
from functions.CTtable import *
from natsort import natsorted
import pandas as pd
import codecs
import allel
import re

# Function to extract info from VCF file and returns 4 things:
# 1) attributes - dictionary of available attributes in VCF
# 2) samples - list of samples in VCF 
# 3) vcf_dataframe - pandas dataframe with attributes' data per column
# 4) chrom_len - pandas dataframe with chromosome names and lengths

def vcf_to_table(vcf_file):
    
    ########## Headers Extraction ##########  
    # extract VCF headers
    # headers contains a list of lists/dictionaries
    # headers[0] - list of all the header lines
    # headers[1] - dict of FILTER field
    # headers[2] - dict of INFO fields
    # headers[3] - dict of FORMAT fields
    # headers[4] - list of samples
    headers = allel.read_vcf_headers(vcf_file)
    
    # headers included in VCF, such as CHROM, POS, etc...
    for line in headers[0]:
        if line[:2] != '##':
            default_fields = np.transpose(line.replace('#', '').split('\t'))

    mandatory_fields = []
    for field in default_fields:
        if 'FILTER' in field:
            continue
        if (field == 'INFO') or (field == 'FORMAT'):
            break
        mandatory_fields.append(field)
    
    
    ########## Attribute Descriptions ##########
    attributes = { 'FORMAT': {},
                   'INFO':   {} }
    
    for info in sorted(headers[2]):
        attributes['INFO'][info] = {}
        attributes['INFO'][info]['Description'] = headers[2][info]['Description']
        attributes['INFO'][info]['Type'] = headers[2][info]['Type']

    for info in sorted(headers[3]):
        attributes['FORMAT'][info] = {}
        attributes['FORMAT'][info]['Description'] = headers[3][info]['Description']
        attributes['FORMAT'][info]['Type'] = headers[3][info]['Type']
    
    # print attributes
    attribute_guide(attributes)
    
    
    ########## Attributes Choosing ##########
    # prompt user to choose attributes to include
    input_instructions = """
    ########## Attribute Choosing ##########

    The following will be included by default: 
    %s,
    ['GT' (per sample), 'DP' (overall and per sample)]

    Enter the fields to include in the prompt below separated
    by commas as such:

    >AD, TYPE, AN, LEN, ...

    To leave the default fields only, leave empty and 
    press ENTER/RETURN.

    >""" % mandatory_fields
    
    raw_input = input(input_instructions)
    fields_to_include = raw_input.replace(' ', '').split(',')
    
    
    ########## VCF Data Extraction ##########
    # extract data from VCF file
    print("\nExtracting data from VCF...\n")
    callset = allel.read_vcf(vcf_file, fields='*', alt_number=1)
    samples = callset['samples']
    
    
    ########## VCF DataFrame ##########
    print("\nCreating dataframe with VCF data...\n")
    vcf_dataframe = pd.DataFrame()
    
    # add mandatory fields to dataframe
    for field in mandatory_fields:
        vcf_dataframe[field] = callset['variants/%s' % field]

    vcf_dataframe['DP'] = callset['variants/DP']
    
    gt = allel.GenotypeArray(callset['calldata/GT']).to_gt()
    for sample in range(len(samples)):
        column_name = "GT_%s" % samples[sample]
        vcf_dataframe[column_name] = gt[:, sample]
        vcf_dataframe[column_name] = vcf_dataframe[column_name].map(lambda gt: codecs.decode(gt, 'UTF-8'))
    
    dp_format = callset['calldata/DP']
    for sample in range(len(samples)):
        column_name = "DP_%s" % samples[sample]
        vcf_dataframe[column_name] = dp_format[:, sample]
    
    # include fields chosen by user
    for field in fields_to_include:
        if field == 'AD':
            ad_format = allel.GenotypeArray(callset['calldata/AD']).to_gt()
            for sample in range(len(samples)):
                column_name = "AD_%s" % samples[sample]
                vcf_dataframe[column_name] = ad_format[:, sample]
                vcf_dataframe[column_name] = vcf_dataframe[column_name].map(lambda ad: codecs.decode(ad, 'UTF-8'))

        elif field == 'GL':
            gl_format = [[[] for gl in range(3)] for sample in samples]
            for sample in range(len(samples)):
                for row in range(len(callset['calldata/GL'])):
                    for gl in range(3):
                        gl_format[sample][gl].append(callset['calldata/GL'][row][sample][gl])

            for sample in range(3):
                column_name = "GL_%s_%i" % (samples[sample], gl)
                vcf_dataframe[column_name] = np.transpose(gl_format[sample][gl])

        elif (field in attributes['FORMAT']) and (field in attributes['INFO']):
            vcf_dataframe[field] = callset['variants/%s' % field]

            field_format = callset['calldata/%s' % field]
            for sample in range(len(samples)):
                column_name = "%s_%s" % (field, samples[sample])
                vcf_dataframe[column_name] = field_format[:, sample]

        elif field in attributes['FORMAT']:
            field_format = callset['calldata/%s' % field]
            for sample in range(len(samples)):
                column_name = "%s_%s" % (field, samples[sample])
                vcf_dataframe[column_name] = field_format[:, sample]

        elif field in attributes['INFO']:
            vcf_dataframe[field] = callset['variants/%s' % field]
    
        
    
    ########## Chromosome Length DataFrame ########## 
    print("\nCreating dataframe with chromosome names and lenghts...\n")
    # extract unique chromosome names from vcf_dataframe
    chromosomes = natsorted(vcf_dataframe.CHROM.unique())
    
    # extract chromosome lengths from 'headers[0]' by using 'contigs' keyword
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
        # search for keyword_id inside contig line
        ID = re.search(keyword_id, contig)
        for chromosome in chromosomes:
            if chromosome == ID.group(1):
                # search for keyword_len inside contig line
                length = re.search(keyword_len, contig)
                chrom_lengths[chromosome] = int(length.group(1))
                
    # convert dictionary to dataframe
    chrom_len = pd.DataFrame(chrom_lengths, index=['LEN']).T
    chrom_len.index.name = 'CHROM'
    
    ########## Return Variables and Dataframes ##########
    print("\nDONE!\n")
    return attributes, samples, vcf_dataframe, chrom_len