import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def stats(samples, vcf_df):
    # function to plot histograms
    def plot_variant_hist(df, chromosome, attribute, bins=50, MSTD=False, xmin=0, xmax=0):
        x = df[attribute][:]
        fig, ax = plt.subplots(figsize=(12, 9), dpi=500, facecolor='w', edgecolor='whitesmoke')
        sns.despine(ax=ax, offset=10)
        ax.hist(x, bins=bins, rwidth=1.0)
        if MSTD:
            xMax = df[attribute][:].max()
            mean = np.mean(x)
            stdDev = np.std(x)
            plt.axvline(x=mean-stdDev, ls='--', color='g', alpha=0.75)
            plt.axvline(x=mean+stdDev, ls='--', color='g', alpha=0.75)
            ax.text(xMax/2, xMax/4, 'Mean: %i\nStdDev: %i' % (mean, stdDev))
        ax.set_xlabel(attribute)
        ax.set_ylabel('Amount')
        if (xmin!=0) and (xmax!=0):
            plt.xlim(xmin, xmax)
        if chromosome=='all':
            ax.set_title('%s-distribution-%s-all-chromosomes' % (attribute, samples_names))
            plt.savefig('%s-distribution-%s-all-chromosomes.png' % (attribute, samples_names))
        else:
            ax.set_title('%s-distribution-chromosome-%s-%s' % (attribute, chromosome, samples_names))
            plt.savefig('%s-distribution-chromosome-%s-%s.png' % (attribute, chromosome, samples_names))
        plt.close('all')
        
    samples_names = '-'.join(map(str, samples))

    chromosomes = chrom_len.index
    lengths = chrom_len.LEN

    gt_samples = ['%s_GT' % sample for sample in samples]
    dp_samples = ['%s_DP' % sample for sample in samples]

    attributes = ['POS', 'QUAL']
    
    for gt_sample in gt_samples:
        attributes.append(gt_sample)

    for dp_sample in dp_samples:
        attributes.append(dp_sample)
    
    plot_variant_hist(vcf_df, 'all', 'QUAL', bins=75, MSTD=True)
    for dp in dp_samples:
        plot_variant_hist(vcf_df, 'all', dp, bins=75, MSTD=True)
        
    for chrom in chromosomes:
        df_chrom = vcf_df[vcf_df.CHROM == chrom]
        for attribute in attributes:
            if (attribute == 'QUAL') or (attribute in dp_samples):
                plot_variant_hist(df_chrom, chrom, attribute, bins=100, MSTD=True)
            else:
                plot_variant_hist(df_chrom, chrom, attribute, bins=50)
