# -*- coding: utf-8 -*-
from wk3_motifs import *

motifs=(
'AACGTA',
'CCCGTT',
'CACCTT',
'GGATTA',
'TTCCGG')

def test_count():
    """get the counts of each of the DNA bases"""
    res = get_count_matrix(Motifs=motifs)
    print('\nCount Matrix:\n')
    [print(k,v) for k,v in res.items()]


def test_profile_matrix():
    """subroutine: get_count_matrix"""
    res = get_profile_matrix(Motifs=motifs)
    print('\nProfile Matrix:\n')
    [print(k,v) for k,v in res.items()]


def test_consensus_motif():
    """subroutine: get_count_matrix"""
    res = get_consensus_motif(Motifs=motifs)
    print('\nConsensus Motif: {}'.format(res))


def test_consensus_motif_score():
    """subroutines: get_consensus_motifs (which needs get_count_matrix)"""
    res = get_consensus_motif_score(Motifs=motifs)
    print('\nConsensus Motif get_consensus_motif_score: {}'.format(res))


def test_string_probability_from_profile_matrix():

    text = 'ACGGGGATTACC'
    profile={
        'A': [0.2, 0.2, 0.0, 0.0, 0.0, 0.0, 0.9, 0.1, 0.1, 0.1, 0.3, 0.0],
        'C': [0.1, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1, 0.2, 0.4, 0.6],
        'G': [0.0, 0.0, 1.0, 1.0, 0.9, 0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
        'T': [0.7, 0.2, 0.0, 0.0, 0.1, 0.1, 0.0, 0.5, 0.8, 0.7, 0.3, 0.4]}
    print('\nstring probability for {}: {}'.format(text,
                                                   get_string_probability_from_profile_matrix(Text=text,
                                                                                              profile=profile)))


def test_get_most_probable_kmer_from_profile_matrix():
    """
    applies sliding k-window across input text
    at each slide, gets probaility for cuurrent window and stores it
    at the end the max proba is used to select the most probable k-mer

    subroutine: get_string_probability_from_profile_matrix
    """

    text = 'ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT'
    k = 5
    profile = {
        'A': [0.2, 0.2, 0.3, 0.2, 0.3],
        'C': [0.4, 0.3, 0.1, 0.5, 0.1],
        'G': [0.3, 0.3, 0.5, 0.2, 0.4],
        'T': [0.1, 0.2, 0.1, 0.1, 0.2]}

    print('\nProfileMostProbableKmer: {}'.format(get_most_probable_kmer_from_profile_matrix(text=text, k=k, profile=profile)))


def test_GreedyMotifSearch():
    t = 5 # nbr of strings
    k = 3
    dna = [
    'GGCGTTCAGGCA',
    'AAGAATCAGTCA',
    'CAAGGAGTTCGC',
    'CACGTCAATCAC',
    'CAATAATATTCG']
    print('\nGreedyMotifSearch results: {}'.format(GreedyMotifSearch(Dna=dna, k=k, nbr_strings=t)))