'''
Created on 2 Mar 2016

@author: steve
'''
#Copyright Stephen Fletcher

import get_ref
import amiR_finder

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
    #TODO: remove position recording
    
    
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
        #need to find real amiR here
        rand_header_seq = ref_seq.popitem()
        #TODO: this is really inefficient - must fix
        cons_seq = [(amiR_finder.optimal_amiRNA(rand_header_seq[1])[0][1],
                     {rand_header_seq[0]:0})]
    return cons_seq


def get_best_amiRs(ref_seq, amiRs, win):
    """
    Recursive function that attempts to find the least amount of amiRNAs for a 
    given window size that target all sequences.  Returns Fasle if this can't be 
    done, or returns a Set of (amiR/*amiR) tuples
    """

    if len(ref_seq)==0:
        return amiRs
        
    elif len(ref_seq)==1:
        for key, value in ref_seq.iteritems():
            targ_amiRNAs = amiR_finder.optimal_amiRNA(value)
            if targ_amiRNAs != []:
                amiRs.add(targ_amiRNAs[0])#select only the first optimal amiRNA
                return amiRs
            else:
                print "no amiR found in {0}".format(key)
                return amiRs
    else:
        cons_seqs = get_most_sub_seqs(ref_seq,win)
        for i in cons_seqs:
            targ_amiRNAs = amiR_finder.optimal_amiRNA(i[0])
            if targ_amiRNAs !=[]: break

        if targ_amiRNAs != []:
            amiRs.add(targ_amiRNAs[0])#select only the first optimal amiRNA
            for key in cons_seqs[0][1].keys():
                if key in ref_seq:
                    del ref_seq[key]
            return get_best_amiRs(ref_seq, amiRs, win)
        return False


def check_results(amiRNA_list, ref_file):
    """
    Confirm the is an amiRNA for all headers in ref_seq
    """
    
    ref_dict = get_ref.get_ref_f_strand(ref_file)

    headers_included=set()
    for amiRNA in amiRNA_list:
        for header, seq in ref_dict.iteritems():
            pos = 0
            while pos < (len(seq)-len(amiRNA)+1):  
                sub_seq = seq[pos:pos+len(amiRNA)]

                if sub_seq == amiRNA:
                    
                    headers_included.add(header)
                pos+=1
                
    return len(ref_dict.keys()) == len(headers_included)


def best_amiR(ref_file, win=21, max_targ = 30):
    """
    Iterates from win to max_targ and finds the minimal amiRNA set
    to target all sequences
    """
    #TODO: if can't find anything - do something
    
    ref_seq = get_ref.get_ref_f_strand(ref_file)
    min_amiRNAs = 10000000
    print '-'*50
    print '{0} reference sequences loaded'.format(len(ref_seq))
    while win < max_targ:
        amiRs=set()
        ref_seq = get_ref.get_ref_f_strand(ref_file) #fix this
        answer = get_best_amiRs(ref_seq, amiRs, win)
        print "win = {0}".format(win)
        if answer is False: print "No amiRNA found"
        else: 
            print "no of amiRs = {0}".format(len(answer))
            if len(answer)<min_amiRNAs:
                min_amiRNAs = len(answer)
                best_set = answer
        win+=1
    print "\nRecommended amiRNAs/amiRNA*s are:\n"
    for i in best_set: print i 
    test_set = []
    for i in best_set: test_set.append(i[1])
    print "\nTest all seqs targeted = {0}".format(check_results(test_set, 
                                                                ref_file))
    print '-'*50  
    
    

def complement(sequence):
    """Provides the complement in the 5' - 3' direction

    Assumption: reference consists of A, G, C, T only

    complement(str) --> str
    """
    d = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(d[c] if c in d else c for c in reversed(sequence))

