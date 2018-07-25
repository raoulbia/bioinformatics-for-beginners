#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def HammingDistance(p, q):
    mismatch = []
    [mismatch.append(x) for x, y in zip(p, q) if x != y]
    return len(mismatch)


# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
# (similar to PatternMatching functio)
def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    n = len(Text)
    k = len(Pattern)
    for i in range(n - k + 1):
        p = Text[i:i + k]
        q = Pattern
        nbr_mismatches = HammingDistance(p, q)
        if nbr_mismatches <= d:
            positions.append(i)
    return positions

# more test cases: http://bioinformaticsalgorithms.com/data/testdatasets/week2/06%20-%20ApproximatePatternMatching.txt
text = 'CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT'
pattern ='ATTCTGGA'
d = 4
print('{} occurs at the following indices:\n{}\n'.format(pattern,
                                                         ApproximatePatternMatching(Text=text,
                                                                                    Pattern=pattern,
                                                                                    d=d)))


def ApproximatePatternCount(Text, Pattern, d):
    count = 0
    n = len(Text)
    k = len(Pattern)
    for i in range(n - k + 1):
        p = Text[i:i + k]
        q = Pattern
        nbr_mismatches = HammingDistance(p, q)
        if nbr_mismatches <= d:
            count = count + 1
    return count

# more test cases: http://bioinformaticsalgorithms.com/data/testdatasets/week2/07%20-%20ApproximatePatternCount.txt
text = 'TTTAGAGCCTTCAGAGG'
pattern ='GAGG'
d = 2
print('{} with {} mismatch allowed occurs {} times'.format(pattern,
                                                           d,
                                                           ApproximatePatternCount(Text=text,
                                                                                    Pattern=pattern,
                                                                                    d=d)))
