#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def SkewArray(Genome):

    Genome = '_' + Genome
    skewArray = [0 for i in Genome]

    n = len(Genome)
    for i in range(n):

        if Genome[i] in ('A', 'T'):
            skewArray[i] = skewArray[i-1]
        elif Genome[i] in ('G'):
            skewArray[i] = skewArray[i-1]+1
        elif Genome[i] in ('C'):
            skewArray[i] = skewArray[i-1]-1

    return skewArray

text = 'CATGGGCATCGGCCATACGCC'
text = 'AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT'
res = SkewArray(Genome=text)
print('SkewArray: {}'.format(" ".join(str(x) for x in res)))