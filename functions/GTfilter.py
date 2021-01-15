import pandas as pd

def filter_similar_gt(samples, vcf_df, genotype):
    
    gt_samples = ['%s_GT' % sample for sample in samples]
    sampleA = gt_samples[0]
    sampleB = gt_samples[1]
    
    filt = (vcf_df[sampleA] == genotype) & (vcf_df[sampleB] == genotype)
    vcf_df = vcf_df.drop(vcf_df[filt].index)
    
    # reset indexes
    vcf_df.reset_index(inplace=True, drop=True)
    
    return vcf_df
