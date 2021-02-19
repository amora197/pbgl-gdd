import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# function to plot histograms
def variant_hist(samples, df, chromosome, attribute, bins=50, MSTD=False, xmin=0, xmax=0):
    fig, ax = plt.subplots(figsize=(12, 9), dpi=500, facecolor='w', edgecolor='whitesmoke')
    sns.despine(ax=ax, offset=10)
    
    samples = '-'.join(map(str, samples))
    x = df[attribute][:]
    
    if xmin == 0:
        remainder = xmax % bins
    else:
        remainder = (xmax - xmin) % bins
        
    if (xmax!=0):
        if remainder == 0:
            ax.hist(x, bins=bins, range=(xmin, xmax), align='mid')
        else:
            ax.hist(x, bins=bins, range=(xmin, xmax-remainder+bins), align='mid')
        plt.xlim(xmin, xmax)
    else:
        ax.hist(x, bins=bins, align='mid')

    bottom, top = plt.ylim()
    
    if (xmax != 0):
        right = xmax
        if (xmin != 0):
            left = xmin
            plt.xlim(left, right)
        else:
            plt.xlim(0, right)
    else:
        left, right = plt.xlim()

    if MSTD:
        mean = np.mean(x)
        stdDev = np.std(x)
        plt.axvline(x=mean-stdDev, ls='--', color='g', alpha=0.75)
        plt.axvline(x=mean+stdDev, ls='--', color='g', alpha=0.75)
        ax.text(right/2, top/2, 'Mean: %i\nStdDev: %i' % (mean, stdDev), fontsize=16)

    ax.set_xlabel(attribute)
    ax.set_ylabel('Amount')

    if chromosome=='all':
        ax.set_title('%s-distribution-%s-all-chromosomes' % (attribute, samples))
        plt.savefig('%s-distribution-%s-all-chromosomes.png' % (attribute, samples))
        plt.savefig('%s-distribution-%s-all-chromosomes.svg' % (attribute, samples))
    else:
        ax.set_title('%s-distribution-chromosome-%s-%s' % (attribute, chromosome, samples))
        plt.savefig('%s-distribution-chromosome-%s-%s.png' % (attribute, chromosome, samples))
        plt.savefig('%s-distribution-chromosome-%s-%s.svg' % (attribute, chromosome, samples))
