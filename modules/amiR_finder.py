
'''
Created on 2 Mar 2016

@author: steve
'''
def optimal_amiRNA(dna_sequence):
    """
    Id any seqs of 21 nt with T at nucleotide 1 and A or T at nucleotide 10
    Input is the target (not the complement)
    optimal_RNAs --> [(amiRNA,*amiRNA)]
    """

    dna_sequence.upper()
    dna_sequence = complement(dna_sequence)
    pos = 0
    optimal_amiRNAs=[]
    while pos <= (len(dna_sequence)-21):
        if dna_sequence[pos] =='T':
            if dna_sequence[pos+9] == 'T' or dna_sequence[pos+9] == 'A':
                comp = complement(dna_sequence[pos:pos+21])
                optimal_amiRNAs.append((dna_sequence[pos:pos+21],comp))
        pos+=1
    return optimal_amiRNAs 

def complement(seq):
    """Provides the complement in the 5' - 3' direction

    Assumption: reference consists of A, G, C, T only

    complement(str) --> str
    """
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return (''.join(d[c] if c in d else c for c in reversed(seq)))    