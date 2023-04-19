#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 19 15:04:04 2023

@author: drewbeckwith
"""
import numpy as np


def generate_sequence(n):
    return("".join(np.random.choice(["A","C","G","T"], n)))



def main():
    randomSeq = generate_sequence(500)
    f = open("/Users/drewbeckwith/documents/github/homocomp/reads10/ref.fasta", "w")
    f.write(">ref\n" + randomSeq)
    
    f.close()

if __name__ == "__main__":
    main()