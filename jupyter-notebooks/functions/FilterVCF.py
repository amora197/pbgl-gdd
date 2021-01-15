import re
import pandas as pd


# function to filter VCF attributes by values
def filter_vcf(vcf_df, filter_list):
    vcf_df_filtered = vcf_df
    
    filters_raw = filter_list.replace(' ', '').split(',')
    comparators = r'((?<!=)==(?!=)|<=|>=|!=|<|>)'
    filters = []
    for filter in filters_raw:
        action = re.split(comparators, filter)
        filters.append(action)
    
    for filter in filters:
        attribute, comparator, value = filter
        value = to_int(value)

        if comparator == '==':
            to_filter = vcf_df_filtered[attribute] == value
        elif comparator == '!=':
            to_filter = vcf_df_filtered[attribute] != value
        elif comparator == '<':
            to_filter = vcf_df_filtered[attribute] < value
        elif comparator == '<=':
            to_filter = vcf_df_filtered[attribute] <= value
        elif comparator == '>':
            to_filter = vcf_df_filtered[attribute] > value
        elif comparator == '>=':
            to_filter = vcf_df_filtered[attribute] >= value

        vcf_df_filtered = vcf_df_filtered[to_filter]
    
    # reset indexes
    vcf_df_filtered.reset_index(inplace=True, drop=True)

    return vcf_df_filtered

# function to try to convert value to integer
def to_int(condition):
    try:
        return int(condition)
    except:
        return condition 
