"""
Created on Wed Jan 18 16:21:06 2023

@author: Drew Beckwith, Saksham Goel

"""

import numpy as np
import pyabpoa as pa
from skbio.alignment import global_pairwise_align_nucleotide
from skbio import TabularMSA, DNA


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
        
    #print(seqsReadyForPoa)
    #print(seqCounts)
    a = pa.msa_aligner()
    res=a.msa(seqsReadyForPoa, out_cons=True, out_msa=True, out_pog='pog.png', incr_fn='') # perform multiple sequence alignment 
                                                                # generate a figure of alignment graph to pog.png
    return(res.cons_seq[0])
    #for seq in res.cons_seq:
       # print(seq)
       
       
def pairwiseAlign(seqs, consenSeq):
    for seq in seqs:
        alignment = global_pairwise_align_nucleotide(DNA(consenSeq), DNA(seqs[seq][0]))
        seqs[seq] = seqs[seq] + ([str(alignment[0][0]), str(alignment[0][1])],)
        print(alignment[0][0])
        print(alignment[0][1])
        print("\n\n")
    return seqs
    

def expansionMean(toExpandWith, consenSeq):
    print(consenSeq)
    print("\n\n\n")
    print("TO ALLIGN WITH")
    print(toExpandWith)
    #x = 0
    #for base in consenSeq:
        #count = 0
        #toDivide = 0
        #for seq in toExpandWith:
            #print(toExpandWith[seq][2][x])
           
                
        
        
    


seqs=[]   
s1 = generate_sequence(100)
s2 = mutate(s1, 0.1, 0.1)
s3 = mutate(s1, 0.1, 0.1)

seqs.append(s1)
seqs.append(s2)
seqs.append(s3)

poaCompressed = compressedPartialOrderAllignment(seqs)

print("### Sequences Before Compression ###")
for seq in seqs:
    print(seq)
print("\n\n")

print("### Sequences After Compression ###")
homecompSeqs = {}
homoId = 0
for seq in seqs:
    print(homoCompress(seq)[0])
    print(homoCompress(seq)[1])
    homecompSeqs[homoId] = homoCompress(seq)
    homoId = homoId + 1
print("\n\n")


print("### Consensus Sequence Before Expansion ###")
print(poaCompressed)
print("\n\n")

print(seqs)
print("### GLOBAL PAIRWISE ALIGN ###")
print("\n\n\n")
readyForExpansion = pairwiseAlign(homecompSeqs, poaCompressed)

expansionMean(readyForExpansion, poaCompressed)


