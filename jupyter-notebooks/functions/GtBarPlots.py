import matplotlib.pyplot as plt
import numpy as np

'''
Function to bar-plot genotype qunatity per chromosome
for two samples per variant position.
'''

def gt_bar_plots(samples, vcf_df, chrom_len, window_size):
    samples_names = '-'.join(map(str, samples))
    
    # extract chromosomes and their respective lengths
    chromosomes = chrom_len.index
    lengths = chrom_len.LEN
    max_length = lengths.max()
    
    # list the samples ID_GT from vcf_df and the genotypes to search for
    gt_samples = ['GT_%s' % sample for sample in samples]
    genotypes = ['0/0', '0/1', '1/1']
    
    # extract CHROM, POS, and samples' GT from vcf_df into new dataframe
    df = vcf_df.filter(['CHROM', 'POS', gt_samples[0], gt_samples[1]])
    
    # loop per chromosome
    for chrom in range(len(chromosomes)):
        plt.figure(num=chrom, figsize=(21,4), dpi=500, facecolor='w', edgecolor='whitesmoke')
        
        if lengths[chrom] < (lengths.max()*0.05):
            continue
        
        # loop per sample
        for sample in range(len(samples)):
            # extract info using current chromosome
            filt = (df.CHROM == chromosomes[chrom])
            df_chrom = df[filt]
            df_chrom = df_chrom.filter(['POS', gt_samples[sample]])
            
            # setup windows
            bins = np.arange(0, lengths[chrom], window_size)
            
            # use window midpoints as x coordinate
            ind = (bins[1:] + bins[:-1]) / 2
            ind = np.append(ind, ind[-1]+window_size)
            
            # count genotypes
            gt00 = []
            gt01 = []
            gt11 = []
            gtNN = []
            
            for x in range(len(bins)):
                if bins[x] == bins[-1]:
                    pos_filter = (df_chrom.POS >= bins[x])
                else:
                    pos_filter = (df_chrom.POS >= bins[x]) & (df_chrom.POS < bins[x+1])

                for gt in genotypes:
                    count = 0
                    myfilter = (df_chrom[gt_samples[sample]] == gt)
                    gt_count = df_chrom.loc[myfilter & pos_filter, 'POS'].count()
                    if gt == '0/0':
                        gt00.append(gt_count)
                    elif gt == '0/1':
                        gt01.append(gt_count)
                    elif gt == '1/1':
                        gt11.append(gt_count)

                gt_total = df_chrom.loc[pos_filter, 'POS'].count()
                gtNN.append(gt_total - gt00[x] - gt01[x] - gt11[x])

            gt00 = np.array(gt00)
            gt01 = np.array(gt01)
            gt11 = np.array(gt11)
            gtNN = np.array(gtNN)
            
            # plot
            plt.subplot(len(samples), 1, sample+1)
            plt.ylabel('%s' % samples[sample])
            ax = plt.gca()
            
            if sample == 1:
                gt00 *= -1
                gt01 *= -1
                gt11 *= -1
                gtNN *= -1
                
            width = window_size
            
            p1 = plt.bar(ind, gtNN, width, color='silver')
            p2 = plt.bar(ind, gt00, width, bottom=gtNN, color='b')
            p3 = plt.bar(ind, gt01, width, bottom=gtNN+gt00, color='g')
            p4 = plt.bar(ind, gt11, width, bottom=gtNN+gt00+gt01, color='y')
                        
            if sample == 0:
                plt.title('bar-plot-samples-%s-chromosome-%s' % (samples_names, chromosomes[chrom]))
                plt.legend((p1[0], p2[0], p3[0], p4[0]), ('Other', '0/0', '0/1', '1/1'), ncol=4, framealpha=0)
                ax.axes.xaxis.set_ticks([])
                
            if (sample + 1) == len(samples):
                plt.xlabel('Position (bp)')
                ticks =  ax.get_yticks()
                ax.set_yticklabels([int(abs(tick)) for tick in ticks])

            plt.xlim([-max_length*0.01, max_length*1.01])
                
        plt.subplots_adjust(hspace=0)
        plt.savefig('bar-plot-samples-%s-chromosome-%s.png' % (samples_names, chromosomes[chrom]))
        plt.savefig('bar-plot-samples-%s-chromosome-%s.svg' % (samples_names, chromosomes[chrom]))
        plt.show()
