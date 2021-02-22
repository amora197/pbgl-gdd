import re
import pandas as pd


# function to extract specific VCF data by attributes
def extract_df_data(vcf_df, filter_expression):
    # make copy of vcf_df
    vcf_df_filtered = vcf_df
    
    # make copy of string with conditions and make item list per condition
    filters_raw = filter_expression
    filters_raw = filters_raw.replace(' ', '').split(',')
    
    # keywords to search for comparatos
    comparators = r'((?<!=)==(?!=)|<=|>=|!=|<|>)'
    
    # split each condition into parts using comparators as key
    filters = []
    for filtre in filters_raw:
        action = re.split(comparators, filtre)
        filters.append(action)
    
    # loop per condition while also checking for field duplicates
    attributes = []
    for filtre in filters:
        attribute, comparator, value = filtre
        
        if attribute in attributes:
            continue
        attributes.append(attribute)
        
        indexes = []
        for line in filters:
            if attribute in line:
                indexes.append(filters.index(line))
                
        filter_list = []
        for index in indexes:
            attribute, comparator, value = filters[index]
            value = to_int(value)
            
            if "GT" in attribute:
                filter_list.append("(vcf_df_filtered['%s'] %s '%s')" % (attribute, comparator, value))
            elif type(value) == int:
                filter_list.append("(vcf_df_filtered['%s'] %s %i)" % (attribute, comparator, value))
            else:
                filter_list.append("(vcf_df_filtered['%s'] %s '%s')" % (attribute, comparator, value))
                
            if index == indexes[-1]:
                break
                
            if comparator == '==':
                filter_list.append("|")
            else:
                filter_list.append("&")
            
        filter_list = ' '.join([str(elem) for elem in filter_list])
        vcf_df_filtered = vcf_df_filtered[eval(filter_list)]
    
    # reset indexes
    vcf_df_filtered.reset_index(inplace=True, drop=True)
    
    return vcf_df_filtered

# function to try to convert value to integer
def to_int(condition):
    try:
        # convert it to integer
        return int(condition)
    except:
        # return as string
        return condition 
