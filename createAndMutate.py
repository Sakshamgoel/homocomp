"""
Created on Wed Jan 18 16:21:06 2023

@author: Drew Beckwith, Saksham Goel

"""

import numpy as np
import pyabpoa as pa
from skbio.alignment import global_pairwise_align_nucleotide
from skbio import TabularMSA, DNA
import statistics
from scipy.stats import poisson

import pandas as pd
import math


from plotnine import *
from plotnine.data import *
#%matplotlib inline


#np.random.seed(seed=2) #test for bottom string dash
#np.random.seed(seed=3) #test for top string dash/very bad allignment
#np.random.seed(seed=4) ##test for top string dash/mismatch



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


def mutatePolymerComp(str, snp_rate, indel_rate):
    x = [c for c in str]
    i = 0
    finalStr = ""
    while i < len(x):
        starter = str[i]
        poisMean = 1
        p = i + 1
        while p < len(str) and str[i] == str[p]:
            p = p + 1
            i = i + 1
            poisMean = poisMean + 1
           
        charNumToAdd = poisson.rvs(mu=poisMean, size=1)
        for y in range(0, charNumToAdd[0]):
            finalStr = finalStr + starter
        i = i + 1
    print(str)
    print(finalStr)
    return(finalStr)


def globalPariwiseAllignment(seqsToAllign):
    if len(seqsToAllign) !=2:
        return "Unable to allign more or less than 2 allignments"
    else:
        topSeq = seqsToAllign[0]
        botSeq = seqsToAllign[1]
        matrix = np.empty((len(botSeq) + 1, len(topSeq) + 1))
        #print("botSeq  ", len(botSeq) + 1)
        matrix[0,0] = 0
        for x in range(1, len(botSeq)+1):
            matrix[x,0] = math.floor(matrix[x-1,0] - 1)
        for x in range(1, len(topSeq)+1):
            matrix[0,x] = int(matrix[0,x-1] - 1)
            
        for x in range(1, len(botSeq)+1):
            for y in range(1, len(topSeq)+1):
                if (topSeq[y-1] == botSeq[x-1]):
                    num1 = matrix[x-1, y-1] + 1
                    num2 = matrix[x-1,y-1]-1
                    num3 = matrix[x-1, y]-1
                    num4 = matrix[x, y-1] -1
                    matrix[x, y] = max(num1,num2,num3,num4)
                else:
                    num2 = matrix[x-1,y-1]-1
                    num3 = matrix[x-1, y]-1
                    num4 = matrix[x, y-1] -1
                    matrix[x, y] = max(num2,num3,num4)
           
        np.set_printoptions(suppress = True)
        backR = len(botSeq)
        backC = len(topSeq)
        allignSeqBot = ""
        allignSeqTop = ""
        while backR > 0 or backC > 0:
            maxNum = []
            maxNum.append(matrix[backR-1, backC-1])
            maxNum.append(matrix[backR-1, backC])
            maxNum.append(matrix[backR, backC-1])
            #print(matrix[backX, backY-1])
            maxToGo = max(maxNum)
            if (maxToGo == matrix[backR-1, backC-1]):
                allignSeqBot = botSeq[backR-1] + allignSeqBot 
                allignSeqTop = topSeq[backC-1] + allignSeqTop 
                backR-=1
                backC-=1
            elif (maxToGo == matrix[backR, backC-1]):
                allignSeqBot = "-" + allignSeqBot 
                allignSeqTop = topSeq[backC-1] + allignSeqTop
                backC-=1
            else:
                allignSeqBot = botSeq[backR-1] + allignSeqBot 
                allignSeqTop = "-" + allignSeqTop 
                backR-=1
        returnTup = ([allignSeqTop, allignSeqBot], matrix[len(botSeq), len(topSeq)])
        #print(allignSeqTop)
        #print(allignSeqBot)
        return returnTup
#seqs = ["ATCAGTGTAT", "TCTGTGTAT"]
#seqs = ["GCAGTC", "GACTC"]
#globalPariwiseAllignment(seqs)
#print()


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
       
def partialOrderAllignment(seqs):

    a = pa.msa_aligner()
    res=a.msa(seqs, out_cons=True, out_msa=True, out_pog='pog.png', incr_fn='') # perform multiple sequence alignment 
                                                                # generate a figure of alignment graph to pog.png
    return(res.cons_seq[0])
    #for seq in res.cons_seq:
       # print(seq)
       
       
def pairwiseAlign(homecompSeqs, consenSeq):
    for id in homecompSeqs:
        #alignment = global_pairwise_align_nucleotide(DNA(consenSeq), DNA(homecompSeqs[id][0]))
        alignment = globalPariwiseAllignment([consenSeq, homecompSeqs[id][0]])
        homecompSeqs[id] = homecompSeqs[id] + ([str(alignment[0][0]), str(alignment[0][1])],)
        print(alignment[0][0])
        print(alignment[0][1])
        print(alignment[1])
        print("\n")
    return homecompSeqs
    

def expansionMean(toExpandWith, consenSeq):

    locDict = {}
    for align in toExpandWith:
        locArray = []
        topSeq = toExpandWith[align][2][0]
        botSeq = toExpandWith[align][2][1]
        consensusLoc = 0
        allignLoc = 0
        while consensusLoc < len(topSeq):
            if topSeq[consensusLoc] != "-": 
                if (topSeq[consensusLoc] == botSeq[consensusLoc]): #match
                    locArray.append(allignLoc)
                    allignLoc+=1
                elif botSeq[consensusLoc] == "-": #deletion
                    locArray.append(-1)
                else: #mismatch
                    locArray.append(-1)
                    allignLoc+=1
            else:
                allignLoc+=1
                
            
            consensusLoc += 1  
        locDict[align] = locArray
        
    finalStr = ""
    finalArr = []
    for base in range(len(consenSeq)):
        if (base == 14):
            numToDivide = 0
            count = 0
            for key in locDict:
                if locDict[key][base] != -1:
                    #print("locDict[key][base]")
                    #print(locDict[key][base])
                    #print("toExpandWith[key][1][locDict[key][base]]")
                 
                    #print(toExpandWith[key][1][locDict[key][base]])
                    numToDivide = numToDivide + toExpandWith[key][1][locDict[key][base]]
                    count += 1
        else:
            numToDivide = 0
            count = 0
            for key in locDict:
                if locDict[key][base] != -1:
                    numToDivide = numToDivide + toExpandWith[key][1][locDict[key][base]]
                    count += 1
        finalBaseLength = 0
        if (count > 0): 
            finalBaseLength = round(numToDivide/count)
            
        charCount = 0
        while charCount < finalBaseLength:
            finalStr = finalStr + consenSeq[base]
            charCount+=1
        
                
    return finalStr


def expansionMedian(toExpandWith, consenSeq):

    locDict = {}
    for align in toExpandWith:
        locArray = []
        topSeq = toExpandWith[align][2][0]
        botSeq = toExpandWith[align][2][1]
        consensusLoc = 0
        allignLoc = 0
        while consensusLoc < len(topSeq):
            if topSeq[consensusLoc] != "-": 
                if (topSeq[consensusLoc] == botSeq[consensusLoc]): #match
                    locArray.append(allignLoc)
                    allignLoc+=1
                elif botSeq[consensusLoc] == "-": #deletion
                    locArray.append(-1)
                else: #mismatch
                    locArray.append(-1)
                    allignLoc+=1
            else:
                allignLoc+=1
                
            
            consensusLoc += 1  
        locDict[align] = locArray
    finalStr = ""
    finalArr = []
    for base in range(len(consenSeq)):
        medianArr = []
        for key in locDict:
            if locDict[key][base] != -1:
                medianArr.append(toExpandWith[key][1][locDict[key][base]])
        finalBaseLength = 0
        if (len(medianArr) > 0):
            finalBaseLength = statistics.median(medianArr)
        charCount = 0
        while charCount < finalBaseLength:
            finalStr = finalStr + consenSeq[base]
            charCount+=1
        
                
    return finalStr


def expansionMode(toExpandWith, consenSeq):

    locDict = {}
    for align in toExpandWith:
        locArray = []
        topSeq = toExpandWith[align][2][0]
        botSeq = toExpandWith[align][2][1]
        consensusLoc = 0
        allignLoc = 0
        while consensusLoc < len(topSeq):
            if topSeq[consensusLoc] != "-": 
                if (topSeq[consensusLoc] == botSeq[consensusLoc]): #match
                    locArray.append(allignLoc)
                    allignLoc+=1
                elif botSeq[consensusLoc] == "-": #deletion
                    locArray.append(-1)
                else: #mismatch
                    locArray.append(-1)
                    allignLoc+=1
            else:
                allignLoc+=1
                
            
            consensusLoc += 1  
        locDict[align] = locArray
    finalStr = ""

    for base in range(len(consenSeq)):
        modeArr = []
        for key in locDict:
            if locDict[key][base] != -1:
                modeArr.append(toExpandWith[key][1][locDict[key][base]])
        finalBaseLength = 0
        if (len(modeArr) > 0):
            finalBaseLength = statistics.median(modeArr)
        charCount = 0
        while charCount < finalBaseLength:
            finalStr = finalStr + consenSeq[base]
            charCount+=1
        
                
    return finalStr

def calcAvgScore(seqs, allignmentToTest):
    avgToDivide = 0
    for seq in seqs:
        #alignment = global_pairwise_align_nucleotide(DNA(allignmentToTest), DNA(seq))
        alignment = globalPariwiseAllignment([allignmentToTest, seq])
        avgToDivide += alignment[1]
    finalAverage = avgToDivide/len(seqs)
    return finalAverage

    
    
def singleRun(howManySeqs, length):
    seqs=[]
    smallbase = generate_sequence(length) 
    for x in range(howManySeqs):
        s1 = mutatePolymerComp(smallbase, 0.1, 0.1)
        seqs.append(s1)
        
    
    
    
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
    print("\n\n")
    print("### GLOBAL PAIRWISE ALIGN ###")
    readyForExpansion = pairwiseAlign(homecompSeqs, poaCompressed)
    
    
    print("\n\n")
    print("### MEAN EXPANSION FINAL STRING ###")
    finalMean = expansionMean(readyForExpansion, poaCompressed)
    print(finalMean)
    print("\n\n")
    
    print("### MEDIAN EXPANSION FINAL STRING ###")
    finalMedian = expansionMedian(readyForExpansion, poaCompressed)
    print(finalMedian)
    print("\n\n")
    
    print("### MODE EXPANSION FINAL STRING ###")
    finalMode = expansionMode(readyForExpansion, poaCompressed)
    print(finalMode)
    print("\n\n")
    
    
    print("### Average Alignment Score Compressed Mean ###")
    print(seqs)
    avgMean = calcAvgScore(seqs, finalMean)
    print(avgMean)
    
    print("### Average Alignment Score Compressed Median ###")
    print(seqs)
    avgMedian = calcAvgScore(seqs, finalMedian)
    print(avgMedian)
    
    
    print("### Average Alignment Score Compressed Mode ###")
    print(seqs)
    avgModes = calcAvgScore(seqs, finalMode)
    print(avgModes)
    
    
    
    print("\n\n")
    print("### Partial Order Allignment Without Compression ###")
    poaNonCompressed = partialOrderAllignment(seqs)
    print(poaNonCompressed)
    print(finalMedian)
    print(smallbase)
    print("\n\n")
    print("### Average Alignment Score Without Compression ###")
    print(seqs)
    avgMeanNonComp = calcAvgScore(seqs, poaNonCompressed)
    print(avgMeanNonComp)
    
    

def multipleRunwithGraph(seqLengthToCheck, howManyTimesEach):
    data = {"comp": [],
            "seqLength": [],
            "value": []
            }
    
    #toCheck = [100]
    for checkNum in seqLengthToCheck:
        meanToAppend = 0
        modeToAppend = 0
        medianToAppend = 0
        nonCompToAppend = 0
        for x in range(howManyTimesEach):
            # List that stores all the mutated sequences
            seqs = []
            
            # base sequence on which mutations are done
            base = generate_sequence(checkNum)
            
            # Generating mutations on the base sequence
            s1 = mutate(base, 0.1, 0.1)
            s2 = mutate(base, 0.1, 0.1)
            s3 = mutate(base, 0.1, 0.1)
            
            seqs.append(s1)
            seqs.append(s2)
            seqs.append(s3)
            
            # Partial Order Align the compressed mutated sequences
            # The result is the consensus sequence of the compressed mutated sequences
            poaCompressed = compressedPartialOrderAllignment(seqs)
            
            # Data type to store the sequences
            # Key: the ID of the sequence
            # Value: A tuple, where the first value of tuple is a string of unique basepairs of that sequence
            # and the second value of tuple is an array of integers, where each integer corresponds to the 
            # frequency of the unique basepairs in the first value.
            homecompSeqs = {}
            homoId = 0
            for seq in seqs:
                homecompSeqs[homoId] = homoCompress(seq)
                homoId = homoId + 1
            
            # Pairwise Align adds the pairwise aligned sequences to the tuple in homocompSeqs
            readyForExpansion = pairwiseAlign(homecompSeqs, poaCompressed)
            finalMean = expansionMean(readyForExpansion, poaCompressed)
            finalMedian = expansionMedian(readyForExpansion, poaCompressed)
            finalMode = expansionMode(readyForExpansion, poaCompressed)
            
            avgMean = calcAvgScore(seqs, finalMean)
            avgMedian = calcAvgScore(seqs, finalMedian)
            avgMode = calcAvgScore(seqs, finalMode)
            
            poaNonCompressed = partialOrderAllignment(seqs)
            avgMeanNonComp = calcAvgScore(seqs, poaNonCompressed)
            
            meanToAppend += avgMean
            medianToAppend += avgMedian 
            modeToAppend += avgMode
            nonCompToAppend += avgMeanNonComp
        finalMean = meanToAppend/howManyTimesEach
        finalMedian = medianToAppend/howManyTimesEach
        finalMode = modeToAppend/howManyTimesEach
        nonCompToAppend = nonCompToAppend/howManyTimesEach
        data["comp"].append("Mean")
        data["comp"].append("Median")
        data["comp"].append("Mode")
        data["comp"].append("NonCompressed")
        data["seqLength"].append(checkNum)
        data["seqLength"].append(checkNum)
        data["seqLength"].append(checkNum)
        data["seqLength"].append(checkNum)
        data["value"].append(finalMean)
        data["value"].append(finalMedian)
        data["value"].append(finalMode)
        data["value"].append(nonCompToAppend)
    
    
    
    
    
        
    df = pd.DataFrame(data)
    print(ggplot(df) + geom_bar(aes(x="seqLength", y="value", fill="comp"), stat = "identity", position = "dodge2"))


singleRun(3, 100)
#multipleRunwithGraph([100, 200], 5)









