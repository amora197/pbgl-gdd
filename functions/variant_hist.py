import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# function to plot histograms
def plot_variant_hist(samples, df, chromosome, attribute, bins=50, MSTD=False, xmin=0, xmax=0):
    samples = '-'.join(map(str, samples))
    x = df[attribute][:]
    fig, ax = plt.subplots(figsize=(12, 9), dpi=500, facecolor='w', edgecolor='whitesmoke')
    sns.despine(ax=ax, offset=10)
    
    if (xmax!=0):
        ax.hist(x, bins=bins, range=(xmin,xmax), rwidth=1.0)
    else:
        ax.hist(x, bins=bins, rwidth=1.0)

    left, right = plt.xlim()
    bottom, top = plt.ylim()

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
