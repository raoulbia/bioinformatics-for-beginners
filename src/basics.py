

# oriC of Vibrio cholerae
ori = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"

### Sliding Window pattern search
def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count
pattern = "TGATCA"
print('number of occurrences of pattern {}: {}\n'.format(pattern, PatternCount(ori, pattern)))


### The Frequent Words Problem
"""
Find the most frequent k-mers in a string.
Input: A string Text and an integer k.
Output: All most frequent k-mers in Text.
"""
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if Pattern in freq:
            freq[Pattern] += 1
        else:
            freq[Pattern] = 1
    return freq
print('Frequency Map: {}\n'.format(FrequencyMap(Text="CGATATATCCATAG", k=3)))

def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words
k = 3
seq = "CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA"
print('Most Frequent {}-mer patterns in input text {} : {}\n'.format(k, seq, FrequentWords(Text=seq, k=k)))


### Reverse & Complement
seq = "GATTACA"

def Reverse(Pattern):
    reversePattern = ''
    for char in Pattern:
        reversePattern = char + reversePattern
    return reversePattern
print('Reverse: {}\n'.format(Reverse(Pattern=seq)))

def Commplement(Pattern):
    dic = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complementPattern = ''.join(dic[x] for x in Pattern)
    return complementPattern
print('Complement: {}\n'.format(Commplement(Pattern=seq)))

def ReverseComplement(Pattern):
    return Reverse(Commplement(Pattern=Pattern))
print('Reverse Complement of {}: {}\n'.format(seq, ReverseComplement(Pattern=seq)))


### Pattern Matching / Pattern Count
"""
Pattern Matching Problem:  Find all occurrences of a pattern in a string.
Input: Strings Pattern and Genome.
Output: All starting positions in Genome where Pattern appears as a substring.
"""
def PatternMatching(Pattern, Genome):
    positions = []
    n = len(Genome)
    k = len(Pattern)
    for i in range(n-k+1):
        if Genome[i:i+k] == Pattern:
            positions.append(i)
    return positions

Pattern = 'CGCG'
Genome='CGCGATACGTTACATACATGATAGACCGCGCGCGATCATATCGCGATTATC'
print('{} occurs in the genome {} at the following indices:\n{}\n'.format(Pattern, Genome, PatternMatching(Pattern=Pattern, Genome=Genome)))