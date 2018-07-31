#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.wk4_gibbs_sampling import *

def test_normalize():
    probas = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}
    probas = {'A': 0.22, 'C': 0.54, 'G': 0.58, 'T': 0.36, 'X': 0.3}
    res = Normalize(probas)
    print(res)

def test_weighted_dice():
    probas = {'ccgc': 0.25, 'cggc': 0.25, 'ggcg': 0.25, 'gcgt': 0.25, 'cgtt': 0.25, 'gtta': 0.25, 'ttag': 0.25}
    probas = {'AA': 0.2, 'TT': 0.2, 'CC': 0.1, 'GG': 0.1, 'AT': 0.4}
    res = WeightedDie(probas)
    print('\n {}'.format(res))

def test_ProfileGeneratedString():

    text = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
    k = 5
    profile = {
        'A': [0.2, 0.2, 0.3, 0.2, 0.3],
        'C': [0.4, 0.3, 0.1, 0.5, 0.1],
        'G': [0.3, 0.3, 0.5, 0.2, 0.4],
        'T': [0.1, 0.2, 0.1, 0.1, 0.2]}

    res = ProfileGeneratedString(text, profile, k)
    print(res)

def test_GibbsSampler():
    t = 5  # nbr of string
    k = 8
    N = 100
    Dna = [
        'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA',
        'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',
        'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',
        'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',
        'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']

    res = GibbsSampler(Dna, k, t, N)
    [print(i) for i in res]