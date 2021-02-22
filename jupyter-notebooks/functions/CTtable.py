import dataframe_image as dfi
import pandas as pd
import numpy as np

def ct_table(samples, vcf_dataframe, chromosome):
    
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

    number_of_samples = len(samples)
    samples_names =  '-'.join(map(str, samples))
    gt_samples = ['GT_%s' % sample for sample in samples]

    for sample in range(number_of_samples):
        
        # check if it's last sample (ie, no comparison possible)
        if gt_samples[sample] == gt_samples[-1]:
            break
        
        sampleA = vcf_dataframe[gt_samples[sample]]
        sampleB = vcf_dataframe[gt_samples[sample+1]]
        
        gt_dict = {
                    'a': 0,
                    'b': 0, 
                    'c': 0, 
                    'd': 0,
                    'e': 0,
                    'f': 0,
                    'g': 0, 
                    'h': 0,
                    'i': 0,
                    'j': 0
                  }

        for gt in range(len(sampleA)):
            if sampleA[gt] == '0/0':
                if sampleB[gt] == '0/0':
                    gt_dict['a'] += 1
                elif sampleB[gt] == '0/1':
                    gt_dict['b'] += 1
                elif sampleB[gt] == '1/1':
                    gt_dict['c'] += 1
                else:
                    gt_dict['j'] +=1
            elif sampleA[gt] == '0/1':
                if sampleB[gt] == '0/0':
                    gt_dict['d'] += 1
                elif sampleB[gt] == '0/1':
                    gt_dict['e'] += 1
                elif sampleB[gt] == '1/1':
                    gt_dict['f'] += 1
                else:
                    gt_dict['j'] +=1
            elif sampleA[gt] == '1/1':
                if sampleB[gt] == '0/0':
                    gt_dict['g'] += 1
                elif sampleB[gt] == '0/1':
                    gt_dict['h'] += 1
                elif sampleB[gt] == '1/1':
                    gt_dict['i'] += 1
                else:
                    gt_dict['j'] +=1
            else:
                gt_dict['j'] +=1

        # convert dictionary values into np.array            
        genotypes = np.array(list(gt_dict.values()))

        # create contingency table
        contingency_table = pd.DataFrame(data=[[genotypes[0],genotypes[1], genotypes[2], genotypes[9]],
                                               [genotypes[3],genotypes[4], genotypes[5], genotypes[9]],
                                               [genotypes[6],genotypes[7], genotypes[8], genotypes[9]],
                                               [genotypes[9],genotypes[9], genotypes[9], genotypes[9]]],
                                         index=[[gt_samples[sample], gt_samples[sample], 
                                                 gt_samples[sample], gt_samples[sample]], 
                                                ['0/0', '0/1', '1/1', 'other']],
                                         columns=[[gt_samples[sample+1], gt_samples[sample+1], 
                                                   gt_samples[sample+1], gt_samples[sample+1]],
                                                ['0/0', '0/1', '1/1', 'other']])

        # save contingency table as image
        contingency_table.dfi.export('CT-filtered-%s-%s.png' % (samples_names, chromosome))
    
    print('\nContingency Table - Chromosome %s\n' % chromosome)
    print(contingency_table)
    
    return contingency_table
