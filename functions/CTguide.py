import pandas as pd

def ct_guide():
    mutant = 'Mutant'
    progenitor = 'Progenitor'
    
    contingency_table_guide = pd.DataFrame(data=[['a', 'b', 'c', ' '],
                                                 ['d', 'e', 'f', ' '],
                                                 ['g', 'h', 'i', ' '],
                                                 [' ', ' ', ' ', 'j']],
                                          index=[[progenitor, progenitor, progenitor, progenitor], 
                                                 ['0/0', '0/1', '1/1', 'other']],
                                          columns=[[mutant, mutant, mutant, mutant],
                                                   ['0/0', '0/1', '1/1', 'other']])
    
    return contingency_table_guide
