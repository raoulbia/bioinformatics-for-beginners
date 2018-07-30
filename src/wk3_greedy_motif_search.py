#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, nbr_strings):
    BestMotifs = []
    for i in range(0, nbr_strings):
        BestMotifs.append(Dna[i][0:k])
    # print(BestMotifs)
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, nbr_strings):
            P = get_profile_matrix(Motifs[0:j])
            Motifs.append(get_most_probable_kmer_from_profile_matrix(Dna[j], k, P))
        if get_consensus_motif_score(Motifs) < get_consensus_motif_score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

