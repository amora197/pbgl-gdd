from functions.VCFtoTable import *
import matplotlib.pyplot as plt
import numpy as np

'''
Function that plots a figure per chromosome, showing
the genotypes of parent and mutant per variant position.
'''

def GTplots(samples, vcf_dataframe, chrom_len):
    samples_names =  '-'.join(map(str, samples))
    
    # extract chromosomes and their respective lengths
    chromosomes = chrom_len.index
    lengths = chrom_len.LEN
    max_length = lengths.max()
    
    # list the samples ID_GT from vcf_dataframe and the genotypes to search for
    gt_samples = ['%s_GT' % sample for sample in samples]
    genotypes = ['0/0', '0/1', '1/1']

    # store the variant positions per GT per sample per chromosome in list of list of list
    # positions has list of chromosomes         [chrom1, chrom2, ..., chromN]
    # chromN has list of samples per chromosome [sample1, sample, ..., sampleN]
    # sampleN has list of genotypes per sample  [pos-of-0/0, pos-of-0/1, pos-of-1/1]
    positions = [ [ [] for sample in samples] for chromosome in chromosomes]
    for chromosome in range(len(chromosomes)):
        for sample in range(len(samples)):
            for gt in genotypes:
                # look for specific GT (0/0, 0/1, or 1/1)
                filt = (vcf_dataframe[gt_samples[sample]]==gt) & (vcf_dataframe.CHROM==chromosomes[chromosome])
                positions[chromosome][sample].append(vcf_dataframe.POS[filt].to_numpy())

    # loop per chromosome
    number_of_samples = len(samples)
    for chromosome in range(len(chromosomes)):
        plt.figure(num=chromosome, figsize=(21,4), dpi=500, facecolor='w', edgecolor='whitesmoke')

        # loop per sample
        for sample in range(number_of_samples):
            plt.subplot(number_of_samples, 1, sample + 1)
            ax = plt.gca()
            ax.axes.yaxis.set_ticks([])

            if sample == 0:
                plt.title('plot-samples-%s-chromosome-%s' % (samples_names, chromosomes[chromosome]))
                ax.axes.xaxis.set_ticks([])

            plt.xlim([-max_length*0.01, max_length*1.01])
            plt.ylabel(samples[sample])

            # plot per GT
            for genotype in range(len(genotypes)):
                # plot 0/0's
                if genotype == 0:
                    x = positions[chromosome][sample][genotype]
                    y = np.zeros(len(positions[chromosome][sample][genotype]))
                    plt.scatter(x, y, marker='|', s=5500, c='b', linewidth=0.025)
                # plot 0/1's
                elif genotype == 1:
                    x = positions[chromosome][sample][genotype]
                    y = np.zeros(len(positions[chromosome][sample][genotype]))
                    plt.scatter(x, y, marker='|', s=5500, c='g', linewidth=0.025)
                # plot 1/1's
                elif genotype == 2:
                    x = positions[chromosome][sample][genotype]
                    y = np.zeros(len(positions[chromosome][sample][genotype]))
                    plt.scatter(x, y, marker='|', s=5500, c='y', linewidth=0.025)
                else:
                    continue
            
            if (sample + 1) == number_of_samples:
                plt.xlabel('Position (bp)')
                
                # add legend
                leg1 = plt.legend(bbox_to_anchor=(1.01,1.0), loc="center left", borderaxespad=0,
                                  title='GT',frameon=False, prop={'weight': 'extra bold'},
                                  labels=['0/0','0/1','1/1'], fontsize='large'
                                 )
                leg2 = plt.legend(bbox_to_anchor=(1.01,1.0), loc="center left", borderaxespad=0,
                                  title=' ', frameon=False, markerfirst=False, prop={'weight': 'extra bold'},
                                  labels=['|','|','|'], labelcolor=['b','g','y'],  fontsize='large'
                                 )
                ax.add_artist(leg1)

        plt.subplots_adjust(hspace=0)
        plt.savefig('plot-samples-%s-chromosome-%s.png' % (samples_names, chromosomes[chromosome]))
        plt.savefig('plot-samples-%s-chromosome-%s.svg' % (samples_names, chromosomes[chromosome]))
        plt.show()
