'''
Created on 2 Mar 2016

@author: steve
'''
#Copyright Stephen Fletcher

import get_ref
import amiR_finder
import copy
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

def get_sub_seqs(ref_seq, window):
    """
    Returns sub-sequences of length win that are present in all reference
    sequences as a list of tuples containing dictionaries
    
    cons_seq = [(subseq, {header,[position,]}),]
     
    """

    cons_seq = []
    seq_lib=get_seq_lib(ref_seq, window)
    
    for sub_seq, occurences in seq_lib.iteritems():
        if len (occurences) == len (ref_seq):
            cons_seq.append((sub_seq,occurences))
    return cons_seq


def get_most_sub_seqs(ref_seq, window=21):
    """
    Returns the most common sub-sequences of length win
    """
    
    
    seq_present={}
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
 
    if cons_seq==[(0,{1:2})]:
        rand_header_seq = ref_seq.popitem()
        cons_seq = [(rand_header_seq[1][:window],{rand_header_seq[0]:0})]
    return cons_seq


def get_best_amiRs(ref_seq, amiRs, win, recur_depth=1000):
    """

    """

    if len(ref_seq)==0:
        return amiRs
        
    elif len(ref_seq)==1:
        for key, value in ref_seq.iteritems():
            targ_amiRNAs = amiR_finder.optimal_amiRNA(value)
            if targ_amiRNAs != []:
                amiRs.append(targ_amiRNAs[0])#select only the first optimal amiRNA
                return amiRs
            else:
                print "no amiR found in {0}".format(key)
                return amiRs
    else:
        cons_seqs = get_most_sub_seqs(ref_seq,win)
        for i in cons_seqs:
            targ_amiRNAs = amiR_finder.optimal_amiRNA(i[0])
            if targ_amiRNAs != []:
                amiRs.append(targ_amiRNAs[0])#select only the first optimal amiRNA
                for key in i[1].keys():
                    if key in ref_seq: 
                        del ref_seq[key]
                    return get_best_amiRs(ref_seq, amiRs, win, recur_depth)
        return False


def best_amiR(ref_file, win=21, max_targ = 23):
    ref_seq = get_ref.get_ref_f_strand(ref_file)
    while win < max_targ:
        amiRs = get_best_amiRs(ref_seq, [], win)
        if amiRs is False:
            win+=1
        else:
            print amiRs
            break
    if win == max_targ: print "No set of amiRs found that cover all sequences"
    
    

print best_amiR('/Users/steve/seq/Neena_lab/tospovirus/ref/tospo_Hanu_13_3_15/TSWV-M.txt')



# def longest_com_seq(ref_file, win = 21):
#     
#     
#     ref_seq = get_ref.get_ref_f_strand(ref_file)
#     
#     max_len = 0 #max lenghth of conserved seq
#     stop_iterating = 30 #round to stop iterating
#     
#     for i in range(stop_iterating):
#         cons_seqs = get_sub_seqs(ref_seq,i)
#         if len(cons_seqs) != 0:
#             max_len = i
#         else:
#             break
#     cons_seqs=cons_seqs = get_sub_seqs(ref_seq,max_len)
#     refs_covered=[]
#     best_cons_seqs = get_most_sub_seqs(ref_seq,win)
#     print '\nMax conserved seq. len. = {0} nt'.format(max_len)
#     print '\nConserved seqs:\n'
#     for i in cons_seqs:
#         print "Seq = {0}".format(i[0])
#     print '\n{0} seqs identified of length {1} are present in {2} reference sequences:'\
#     .format(len(best_cons_seqs ), win, len(best_cons_seqs[0][1]))
#     for seq in best_cons_seqs:
#         print "\nSeq = {0}\n".format(seq[0])
#         for key,value in seq[1].iteritems():
#             print "Reference = {0} at position {1}".format(key,value)
#             if key not in refs_covered:
#                 refs_covered.append(key)


def complement(sequence):
    """Provides the complement in the 5' - 3' direction

    Assumption: reference consists of A, G, C, T only

    complement(str) --> str
    """
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(d[c] if c in d else c for c in reversed(sequence))

