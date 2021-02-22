import matplotlib.pyplot as plt
import numpy as np

'''
Plot all the genotypes per chromosomes per sample in one figure.
'''

def gt_plot(samples, vcf_dataframe, chrom_len, linethickness=0.02):
    
    samples_names =  '-'.join(map(str, samples))

    chromosomes = chrom_len.index
    lengths = chrom_len.LEN
    max_length = lengths.max()

    gt_samples = ['GT_%s' % sample for sample in samples]
    genotypes = ['0/0', '0/1', '1/1']

    positions = [ [ [] for sample in samples] for chromosome in chromosomes]
    for chromosome in range(len(chromosomes)):
        for sample in range(len(samples)):
            for gt in genotypes:
                filt = (vcf_dataframe[gt_samples[sample]]==gt) & (vcf_dataframe.CHROM==chromosomes[chromosome])
                positions[chromosome][sample].append(vcf_dataframe.POS[filt].to_numpy())

    subplot_index = 1
    plt.figure(figsize=(11,8.5), dpi=500, facecolor='w', edgecolor='whitesmoke')
    plt.suptitle('gt-plot-%s' % samples_names, fontsize=20, fontweight='black', verticalalignment='bottom')
    for chromosome in range(len(chromosomes)):
        for sample in range(len(gt_samples)):            
            plt.subplot(len(chromosomes)*len(samples), 1, subplot_index)
            ax = plt.gca()
            ax.axes.yaxis.set_ticks([])

            if subplot_index != (len(chromosomes)*len(samples)):
                ax.axes.xaxis.set_ticks([])

            plt.xlim([-max_length*0.01, max_length*1.01])
            plt.ylabel('%s_%s' % (samples[sample], chromosomes[chromosome]), rotation=0, labelpad=65, 
                       verticalalignment='center', fontsize=7)

            for genotype in range(len(genotypes)):
                if genotype == 0:
                    x = positions[chromosome][sample][genotype]
                    y = np.zeros(len(positions[chromosome][sample][genotype]))
                    plt.scatter(x, y, marker='|', s=5500, c='b', linewidth=linethickness)
                elif genotype == 1:
                    x = positions[chromosome][sample][genotype]
                    y = np.zeros(len(positions[chromosome][sample][genotype]))
                    plt.scatter(x, y, marker='|', s=5500, c='g', linewidth=linethickness)
                elif genotype == 2:
                    x = positions[chromosome][sample][genotype]
                    y = np.zeros(len(positions[chromosome][sample][genotype]))
                    plt.scatter(x, y, marker='|', s=5500, c='y', linewidth=linethickness)
                else:
                    continue
            subplot_index += 1

    leg1 = plt.legend(bbox_to_anchor=(1.01,len(chromosomes)*len(samples)*0.5), loc="center left", borderaxespad=0,
                      frameon=False, 
                      prop={'size': 20},
                      labels=['0/0','0/1','1/1'], fontsize='large')
    leg2 = plt.legend(bbox_to_anchor=(1.01,len(chromosomes)*len(samples)*0.5), loc="center left", borderaxespad=0,
                      frameon=False, markerfirst=False, 
                      prop={'weight': 'extra bold', 'size': 20},
                      labels=['|','|','|'], labelcolor=['b','g','y'],  fontsize='large')

    ax.add_artist(leg1)
    plt.subplots_adjust(hspace=0)
    plt.xlabel('Position (bp)', fontsize=20)
    fig_name_png = 'gt-plot-%s.png'  % samples_names
    plt.savefig(fig_name_png, bbox_inches='tight')
    fig_name_svg = 'gt-plot-%s.svg'  % samples_names
    plt.savefig(fig_name_svg, bbox_inches='tight')
    plt.show()
