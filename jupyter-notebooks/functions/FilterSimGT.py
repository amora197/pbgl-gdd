import pandas as pd

# function to filter out variant positions POS from VCF dataframe (vcf_df)
# where the samples have the same genotypes, such as 0/0 and 1/1 

def filter_sim_gt(samples, vcf_dataframe, genotypes):
    vcf_df = vcf_dataframe
    gt_samples = ['GT_%s' % sample for sample in samples]

    # loop through each genotype in genotypes list
    for genotype in genotypes: 
        # list storing conditions
        filter_list = []
        for sample in gt_samples:
            filter_list.append("(vcf_df['%s'] == '%s')" % (sample, genotype))
            if sample == gt_samples[-1]:
                break
            filter_list.append("&")
        
        # convert list into one long string
        gt_filter = ' '.join([str(elem) for elem in filter_list])
        
        # filter vcf by by converting gt_filter from string to executable code
        vcf_df = vcf_df.drop(vcf_df[eval(gt_filter)].index)

    # reset dataframe indexes for plotting functions
    vcf_df.reset_index(inplace=True, drop=True)

    return vcf_df