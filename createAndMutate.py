"""
Created on Wed Jan 18 16:21:06 2023

@author: Drew Beckwith, Saksham Goel

"""
import numpy as np

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

s1 = generate_sequence(100)

s2 = mutate(s1, 0.1, 0.1)
print(s1, end = '')
print('')
print(s2)
