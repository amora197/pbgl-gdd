from functions.AttributeGuide import *
import allel


# function to extract attributes from VCF headers and returns dictionary of 
# available fields to choose from in vcf_to_table() function.

def extract_attributes(vcf_file):
    # extract VCF headers
    # headers contains a list of lists/dictionaries
    # headers[0] - list of all the header lines
    # headers[1] - dict of FILTER field
    # headers[2] - dict of INFO fields
    # headers[3] - dict of FORMAT fields
    # headers[4] - list of samples

    headers = allel.read_vcf_headers(vcf_file)
    
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
    
    return attributes
