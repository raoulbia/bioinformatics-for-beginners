#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
Recall that Profile(Motifs) returns a profile matrix that we can use to get (k-mer) motifs.

In general, we can begin from a collection of randomly chosen k-mers Motifs in Dna,
construct Profile(Motifs), and use this profile to generate a new collection of k-mers:

    Motifs(Profile(Motifs), Dna).

our hope is that Motifs(Profile(Motifs), Dna) has a better score than the
original collection of k-mers Motifs. We can then form the profile matrix of these k-mers,

    Profile(Motifs(Profile(Motifs), Dna))

and use it to form the most probable k-mers,

    Motifs(Profile(Motifs(Profile(Motifs), Dna)), Dna).

etc ...

"""

# first, import the random package
import random
# Next, copy your RandomizedMotifSearch function (along with all required subroutines)
# from Motifs.py below this line
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M

    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs

def RandomMotifs(Dna, k, t):
    """Important Note: an independent random starting position should be generated for each Dna string"""
    m = len(Dna[0])
    # print('\nm: {} k: {} m-k+1:{}'.format(m,k, m-k))

    random_motifs = []
    for i in range(t):
        r = random.randint(1, m - k)
        # print(r)
        # if r > m - k :
        #     r = r - k
        #     print('  {}'.format(r+k))
        # print(Dna[i][r:r+k])

        random_motifs.append(Dna[i][r:r + k])
    return random_motifs


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)

    profileMatrix = Count(Motifs=Motifs)
    for k, v in profileMatrix.items():
        v[:] = [i / t + 4 for i in v] # +4 to account for initialzing count matrix with 1's
    return profileMatrix

def Count(Motifs):
    countMatrix = {}
    k = len(Motifs[0])
    # print(Motifs)
    # init dict of lists (one list for each DNA base)
    for symbol in "ACGT":
        countMatrix[symbol] = []
        for j in range(k):
            countMatrix[symbol].append(1.0)
    # populate dict of lists
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            countMatrix[symbol][j] += 1
    return countMatrix

def Motifs(Profile, Dna):
    k = len(Profile['A'])
    motifs = []
    for i in range(len(Dna)):
        motifs.append(most_probable_kmer(Dna[i], k, Profile))
    return motifs

def most_probable_kmer(text, k, profile):
    n = len(text)
    most_probable = {}
    for i in range(n - k + 1):
        pattern = text[i:i + k]
        most_probable[pattern] = Pr(pattern, profile)
    return max(most_probable, key=most_probable.get)

def Pr(Text, profile):
    p = 1
    pos = 0 # init pos
    for i in Text:
        p *= profile[i][pos]
        # move one position to the right
        pos += 1
    return p

def Score(Motifs):
    # my approach
    # score = 0
    # consensus = Consensus(Motifs=Motifs)
    # motifs_by_col = [''.join(x) for x in zip(*Motifs)]
    # n = len(motifs_by_col)
    # for i in range(n):
    #     #score = sum the number of symbols in the j-th column of Motifs that do not match
    #     count = [1 for x in motifs_by_col[i] if not x == consensus[i]]
    #     score += sum(count) * (1 / n ** 4)  # multiply by (1/n**4) to account for intializing the count matrix with ones
    # return score

    # Jorge Iván Fuentes Rosado
    score = 0
    consensus = Consensus(Motifs)
    for row in Motifs:
        for c1, r1 in zip(consensus, row):
            if c1 != r1:
                score += 1
    return score


def Consensus(Motifs):
    k = len(Motifs[0])
    count_matrix = Count(Motifs)
    consensus = ""
    for i in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count_matrix[symbol][i] > m:
                m = count_matrix[symbol][i]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

# Copy the ten strings occurring in the hyperlinked DosR dataset below.
Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC",
"CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG",
"ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC",
"GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC",
"GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG",
"CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA",
"GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA",
"GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG",
"GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG",
"TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

# set t equal to the number of strings in Dna, k equal to 15, and N equal to 100.
t = 10  # nbr of string
k = 15
N = 100


# Call RandomizedMotifSearch(Dna, k, t) N times, storing the best-scoring set of motifs
# resulting from this algorithm in a variable called BestMotifs
BestMotifs_Iter = []
for i in range(N):
    BestMotifs_Iter.append(RandomizedMotifSearch(Dna=Dna, k=k, t=t))
# print(BestMotifs_Iter)
best_score = 1
best_motifs = []
for motifs in BestMotifs_Iter:
    score = Score(motifs)
    if score < best_score:
        best_score = score
        best_motifs = motifs
[print(i) for i in best_motifs]
print(best_score)



""" for reference
Solution by Jorge Iván Fuentes Rosado


import random
def RandomMotifs(Dna, k, t):
    ranMotifs=[]
    for i in range(t):
        j=random.randint(0, len(Dna[0])-k)
        ranMotifs.append(Dna[i][j:k+j])
    return ranMotifs


def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {} # output variable
    profile = CountWithPseudocounts(Motifs)
    dt=t+len(profile)
    for i in 'ACTG':
        for j in range(k):
            profile[i][j] = profile[i][j]/dt
    return profile

def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {}
    for symbol in 'ACGT':
        count[symbol]=[]
        for j in range(0,len(Motifs[0])):
            count[symbol].append(1.0)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            symbol = Motifs[i][j]
            count[symbol][j] += 1.0
    return count

def Motifs(Profile, Dna):
    k=len(Profile['A'])
    motifs=[]
    for i in range(len(Dna)):
        motifs.append(ProfileMostProbablePattern(Dna[i], k, Profile))
    return motifs

def ProfileMostProbablePattern(Text, k, Profile): #finds the most likely Kmer in a string that satisfy a profile matrix.
    a =-1
    ProbPattern=""
    for i in range(len(Text)-k+1):
        if Pr(Text[i:i+k], Profile) > a:
            a=Pr(Text[i:i+k], Profile)
            ProbPattern=Text[i:i+k]
    return ProbPattern

def Pr(Text, Profile):# finds probability of DNA sequence Klenght given  Profile.
    pr=1
    for i in range(len(Text)):
        pr*=Profile[Text[i]][i]
    return pr

def Score(Motifs):
    H_D = 0
    p = Consensus(Motifs)
    for i in range(len(Motifs)):
        for j in range(len(Motifs[0])):
            if Motifs[i][j] != p[j]:
                H_D += 1
    return H_D

def Consensus(Motifs):
    count=CountWithPseudocounts(Motifs)
    consensus=""
    for j in range(k):
        m=0
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m=count[symbol][j]
                frequentSymbol=symbol
        consensus+=frequentSymbol
    return consensus

def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs




# Copy the ten strings occurring in the hyperlinked DosR dataset below.
Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]


def RepeatedRandomizedMotifSearch(Dna, k, t):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(N):
        Motifs = RandomizedMotifSearch(Dna, k, t)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs
# set t equal to the number of strings in Dna, k equal to 15, and N equal to 100.
t=len(Dna)
k=15
N=100
# Call RandomizedMotifSearch(Dna, k, t) N times, storing the best-scoring set of motifs
# resulting from this algorithm in a variable called BestMotifs
BestMotifs=RepeatedRandomizedMotifSearch(Dna, k, t)
# Print the BestMotifs variable
print(BestMotifs)
# Print Score(BestMotifs)
print(Score(BestMotifs))

"""