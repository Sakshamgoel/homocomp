"""
Created on Wed Jan 18 16:21:06 2023

@author: Drew Beckwith, Saksham Goel

"""

import numpy as np
import pyabpoa as pa


np.random.seed(seed=1)

def generate_sequence(n):
    return("".join(np.random.choice(["A","C","G","T"], n)))

def mutate(str, snp_rate, indel_rate):
    x = [c for c in str]
    i = 0
    while i < len(x):
        if np.random.random() < snp_rate:
            x[i] = generate_sequence(1)
        if np.random.random() < indel_rate:
            length = np.random.geometric(0.5)
            if np.random.random() < 0.5: # insertion
                x[i] = x[i] + generate_sequence(length)
            else:
                for j in range(i,i+length):
                    if j < len(x):
                        x[j] = ""
                        i += 1
        i += 1
    return("".join(x))

"""
homoCompress

This function preforms homopolymer compression on a given
sequence.


@param seq: The sequence to preform homopolyer compression on

@return a tuple with the first element being an array with the 
    string values of the homopolymer compression and the second
    array in the tuple is an array with the count of each base at
    that location in the array.

"""
def homoCompress(seq):
    compSeq = []
    countSeq = []
    i = 0
    while i < len(seq):
        starter = seq[i]
        count = 1
        p = i + 1
        while p < len(seq) and seq[i] == seq[p]:
            p = p + 1
            i = i + 1
            count = count + 1
        compSeq.append(starter)
        countSeq.append(count)
        i = i + 1
    return ("".join(compSeq), countSeq)


def compressedPartialOrderAllignment(seqs):
    seqsReadyForPoa = []
    seqCounts = []
    for seq in seqs:
        homoCompTuple = homoCompress(seq)
        seqReady = homoCompTuple[0]
        seqSingleCount = homoCompTuple[1]
        seqsReadyForPoa.append(seqReady)
        seqCounts.append(seqSingleCount)
        
    print(seqsReadyForPoa)
    print(seqCounts)
    a = pa.msa_aligner()
    res=a.msa(seqsReadyForPoa, out_cons=True, out_msa=True, out_pog='pog.png', incr_fn='') # perform multiple sequence alignment 
                                                                # generate a figure of alignment graph to pog.png

    for seq in res.cons_seq:
        print(seq) 
    
seqs=[
    'CCGAAGA',
    'CCGAACTCGA',
    'CCCGGAAGA',
    'CCGAAGA'
    ]
compressedPartialOrderAllignment(seqs)
    
    
    
s1 = generate_sequence(100)

s2 = mutate(s1, 0.1, 0.1)
#print(s1, end = '')
#print('')
#print(s2)
#print('')
#print(homoCompress(s1))
