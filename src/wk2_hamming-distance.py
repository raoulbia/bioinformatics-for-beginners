#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def HammingDistance(p, q):
    mismatch = []
    [mismatch.append(x) for x, y in zip(p, q) if x != y]
    return len(mismatch)

# Test 0  Sample Dataset (your code is not run on this dataset)
# Res = 3
p ='GGGCCGTTGGT'
q = 'GGACCGTTGAC'
res = HammingDistance(p, q)
print('Hamming Distance: {}. Expected: 3'.format(res))

# Test 1  Make sure you return a positive number
# Res = 4
p = 'AAAA'
q = 'TTTT'

res = HammingDistance(p, q)
print('Hamming Distance: {}. Expected: 4'.format(res))

# Test 2 Check if you're accidentally returning Edit Distance instead of Hamming Distance
# Res = 8
p = 'ACGTACGT'
q = 'TACGTACG'

res = HammingDistance(p, q)
print('Hamming Distance: {}. Expected: 8'.format(res))

# Test 3 Check if you're accidentally counting the number of matches instead of the number of mismatches
# Res = 6
p = 'ACGTACGT'
q = 'CCCCCCCC'

res = HammingDistance(p, q)
print('Hamming Distance: {}. Expected: 6'.format(res))

# Test 4 Example where the two strings have no matches
# Res = 8
p = 'ACGTACGT'
q = 'TGCATGCA'

res = HammingDistance(p, q)
print('Hamming Distance: {}. Expected: 8'.format(res))

# Test 5 Check if you're missing the first index of the strings
# Res = 15
p = 'GATAGCAGCTTCTGAACTGGTTACCTGCCGTGAGTAAATTAAAATTTTATTGACTTAGGTCACTAAATACT'
q = 'AATAGCAGCTTCTCAACTGGTTACCTCGTATGAGTAAATTAGGTCATTATTGACTCAGGTCACTAACGTCT'

res = HammingDistance(p, q)
print('Hamming Distance: {}. Expected: 15'.format(res))

# Test 6 Check if you're missing the last index of the strings
# Res = 28
p = 'AGAAACAGACCGCTATGTTCAACGATTTGTTTTATCTCGTCACCGGGATATTGCGGCCACTCATCGGTCAGTTGATTACGCAGGGCGTAAATCGCCAGAATCAGGCTG'
q = 'AGAAACCCACCGCTAAAAACAACGATTTGCGTAGTCAGGTCACCGGGATATTGCGGCCACTAAGGCCTTGGATGATTACGCAGAACGTATTGACCCAGAATCAGGCTC'

res = HammingDistance(p, q)
print('Hamming Distance: {}. Expected: 28'.format(res))

# Test 7  Full dataset
# Res = 103
p = 'AGCGTGCCGAAATATGCCGCCAGACCTGCTGCGGTGGCCTCGCCGACTTCACGGATGCCAAGTGCATAGAGGAAGCGAGCAAAGGTGGTTTCTTTCGCTTTATCCAGCGCGTTAACCACGTTCTGTGCCGACTTT'
q = 'AGAAACAGACCGCTATGTTCAACGATTTGTTTTATCTCGTCACCGGGATATTGCGGCCACTCATCGGTCAGTTGATTACGCAGGGCGTAAATCGCCAGAATCAGGCTGAAACCCTACGGACAGGTTTACGAACCT'

res = HammingDistance(p, q)
print('Hamming Distance: {}. Expected: 103'.format(res))