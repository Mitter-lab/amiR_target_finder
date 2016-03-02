'''
Created on 2 Mar 2016

@author: steve
'''
#Copyright Stephen Fletcher

import get_ref
"""
Set conserved region sequence length (ie. 21nt)
Script returns all sequences of that length that
appear in all provided refernce sequences, 
along with their positions in each.  Additionally,
longest_com_seq function returns the longest conserved 
sequence (even when less that the set length), and the
most conserved sequence of the set length.

Default window length = 21
"""

def get_seq_lib(ref_seq, window):

    seq_lib = {}

    pos = 0


    for header,seq in ref_seq.iteritems():
        pos = 0
        while pos < (len(seq)-window+1):  
            sub_seq = seq[pos:pos+window]
            #if new k-mer
            if sub_seq not in seq_lib:
                seq_lib[sub_seq]={header:[pos]}
            #existing k-mer in existing header seq --> append pos
            elif header in seq_lib[sub_seq]:
                seq_lib[sub_seq][header].append(pos)
            #existing k-mer in new header seq --> header = [pos]
            elif header not in seq_lib[sub_seq]:
                seq_lib[sub_seq][header] = [pos]
            else:
                pass
            pos +=1
    return seq_lib

def get_sub_seqs(ref_file, window):
    """
    Returns sub-sequences of length win that are present in all reference
    sequences as a list of tuples containing dictionaries
    
    cons_seq = [(subseq, {header,[position,]}),]
     
    """
    ref_seq = get_ref.get_ref_f_strand(ref_file)
    cons_seq = []
    seq_lib=get_seq_lib(ref_seq, window)
    
    for sub_seq, occurences in seq_lib.iteritems():
        if len (occurences) == len (ref_seq):
            cons_seq.append((sub_seq,occurences))
    return cons_seq


def get_most_sub_seqs(ref_file, window=21):
    """
    Returns the most common sub-sequences of length win
    """
    
    
    seq_present={}
    ref_seq = get_ref.get_ref_f_strand(ref_file)
    seq_lib=get_seq_lib(ref_seq, window)
    
    cons_seq=[(0,{1:2})]
 
    set_length = 1
    for sub_seq, occurences in seq_lib.iteritems():
 
        if len(occurences) > len(cons_seq[0][1]):
 
            cons_seq=[(sub_seq,occurences)]
            set_length = len(occurences)
            for key in occurences.keys(): 
                seq_present[key]=1
        elif len(occurences) == set_length and set_length > 1:
            cons_seq.append((sub_seq,occurences)) 
            for key in occurences.keys(): 
                seq_present[key]=1
 
 
    print 'Most conserved k-mer present in {0} reference sequences'\
    .format(len(seq_present))
    return cons_seq


def longest_com_seq(ref_seq, win = 21):
    max_len = 0 #max lenghth of conserved seq
    stop_iterating = 30 #round to stop iterating
    
    for i in range(stop_iterating):
        cons_seqs = get_sub_seqs(ref_seq,i)
        if len(cons_seqs) != 0:
            max_len = i
        else:
            break
    refs_covered=[]
    best_cons_seqs = get_most_sub_seqs(ref_seq,win)
    print '\nMax conserved seq. len. = {0} nt'.format(max_len)
    print '\n{0} seqs identified of length {1} are present in {2} reference sequences:'\
    .format(len(best_cons_seqs ), win, len(best_cons_seqs[0][1]))
    for seq in best_cons_seqs:
        print "\nSeq = {0}\n".format(seq[0])
        for key,value in seq[1].iteritems():
            print "Reference = {0} at position {1}".format(key,value)
            if key not in refs_covered:
 
                refs_covered.append(key)


def complement(sequence):
    """Provides the complement in the 5' - 3' direction

    Assumption: reference consists of A, G, C, T only

    complement(str) --> str
    """
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(d[c] if c in d else c for c in reversed(sequence))

