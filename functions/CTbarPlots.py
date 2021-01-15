import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def CTbarPlots(samples, vcf_df, chrom_len, window_size):
    
    ''' 
    Dictionary to create contingency table
                                     Mutant (SampleB)
                         +---------+---------+---------+-------+
                         |  0/0    |   0/1   |   1/1   | other |
               +---------+---------+---------+---------+-------+
               |   0/0   |   a     |    b    |    c    |       |
    Parent     +---------+---------+---------+---------+       +
   (SampleA)   |   0/1   |   d     |    e    |    f    |       |
               +---------+---------+---------+---------+       +
               |   1/1   |   g     |    h    |    i    |       |
               +---------+---------+---------+---------+       +
               |  other  |                                 j   |
               +---------+---------+---------+---------+-------+
    '''
    
    samples_names =  '-'.join(map(str, samples))
    chromosomes = chrom_len.index
    lengths = chrom_len.LEN
    
    gt_samples = ['%s_GT' % sample for sample in samples]
    genotypes = ['0/0', '0/1', '1/1']
    
    df = vcf_df.filter(['CHROM', 'POS', gt_samples[0], gt_samples[1]])
    
    for chrom in range(len(chromosomes)):
        if lengths[chrom] < (lengths.max()*0.05):
            continue
        
        chrom_filter = (df.CHROM == chromosomes[chrom])
        df_chrom = df[chrom_filter]
        
        a = []
        b = []
        c = []
        d = []
        e = []
        f = []
        g = []
        h = []
        i = []
        j = []
        
        # setup windows
        bins = np.arange(0, lengths[chrom], window_size)
        
        # use window midpoints as x coordinate
        ind = (bins[1:] + bins[:-1]) / 2
        ind = np.append(ind, ind[-1]+window_size)
        
        for x in range(len(bins)):
            if bins[x] == bins[-1]:
                pos_filter = (df_chrom.POS >= bins[x])
            else:
                pos_filter = (df_chrom.POS >= bins[x]) & (df_chrom.POS < bins[x+1])

            df_bin = df_chrom[pos_filter]

            sampleA = df_bin[gt_samples[0]]
            sampleB = df_bin[gt_samples[1]]

            for sampleA_gt in genotypes:
                for sampleB_gt in genotypes:     
                    gt_count = df_bin.loc[(sampleA == sampleA_gt) & (sampleB == sampleB_gt), 'POS'].count()

                    if sampleA_gt == '0/0':
                        if sampleB_gt == '0/0':
                            a.append(gt_count)
                        elif sampleB_gt == '0/1':
                            b.append(gt_count)
                        elif sampleB_gt == '1/1':
                            c.append(gt_count)
                    elif sampleA_gt == '0/1':
                        if sampleB_gt == '0/0':
                            d.append(gt_count)
                        elif sampleB_gt == '0/1':
                            e.append(gt_count)
                        elif sampleB_gt == '1/1':
                            f.append(gt_count)
                    elif sampleA_gt == '1/1':
                        if sampleB_gt == '0/0':
                            g.append(gt_count)
                        elif sampleB_gt == '0/1':
                            h.append(gt_count)
                        elif sampleB_gt == '1/1':
                            i.append(gt_count)

            gt_total = df_bin['POS'].count()
            gt_other = gt_total - a[x] - b[x] - c[x] - d[x] - e[x] - f[x] - g[x] - h[x] - i[x]
            j.append(gt_other)
            
        a = np.array(a)
        b = np.array(b)
        c = np.array(c)
        d = np.array(d)
        e = np.array(e)
        f = np.array(f)
        g = np.array(g)
        h = np.array(h)
        i = np.array(i)
        j = np.array(j)
        
        # plot
        plt.figure(num=chrom, figsize=(21,4), dpi=500, facecolor='w', edgecolor='whitesmoke')
        plt.ylabel('Quantity')
        plt.xlabel('Position (bp)')
        plt.title('contingency-bar-plot-samples-%s-chromosome-%s' % (samples_names, chromosomes[chrom]))

        width = window_size

        p1 = plt.bar(ind,  a, width, color='tab:blue')
        p2 = plt.bar(ind,  b, width, color='tab:orange', bottom=a)
        p3 = plt.bar(ind,  c, width, color='tab:green',  bottom=a+b)
        p4 = plt.bar(ind,  d, width, color='tab:red',    bottom=a+b+c)
        p5 = plt.bar(ind,  e, width, color='tab:purple', bottom=a+b+c+d)
        p6 = plt.bar(ind,  f, width, color='tab:brown',  bottom=a+b+c+d+e)
        p7 = plt.bar(ind,  g, width, color='tab:pink',   bottom=a+b+c+d+e+f)
        p8 = plt.bar(ind,  h, width, color='tab:gray',   bottom=a+b+c+d+e+f+g)
        p9 = plt.bar(ind,  i, width, color='tab:olive',  bottom=a+b+c+d+e+f+g+h)
        p10 = plt.bar(ind, j, width, color='tab:cyan',   bottom=a+b+c+d+e+f+g+h+i)

        plt.legend((p1[0], p2[0], p3[0], p4[0], p5[0], p6[0], p7[0], p8[0], p9[0], p10[0]), 
                   ('a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'),
                    ncol=10, framealpha=0)
        
        plt.savefig('contingency-bar-plot-samples-%s-chromosome-%s.png' % (samples_names, chromosomes[chrom]))
        plt.savefig('contingency-bar-plot-samples-%s-chromosome-%s.svg' % (samples_names, chromosomes[chrom]))
        plt.show()
