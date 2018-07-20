

""" Exercise from Section 1.5
How many combined occurrences of "ATGATCAAG" and "CTTGATCAT" can you find in the
ori region of the Thermotoga petrophila bacterium
"""

Text = 'AACTCTATACCTCCTTTTTGTCGAATTTGTGTGATTTATAGAGAAAATCTTATTAACTGAAACTAAAATGGTAGGTTTGGTGGTAGGTTTTGTGTACATTTTGTAGTATCTGATTTTTAATTACATACCGTATATTGTATTAAATTGACGAACAATTGCATGGAATTGAATATATGCAAAACAAACCTACCACCAAACTCTGTATTGACCATTTTAGGACAACTTCAGGGTGGTAGGTTTCTGAAGCTCTCATCAATAGACTATTTTAGTCTTTACAAACAATATTACCGTTCAGATTCAAGATTCTACAACGCTGTTTTAATGGGCGTTGCAGAAAACTTACCACCTAAAATCCAGTATCCAAGCCGATTTCAGAGAAACCTACCACTTACCTACCACTTACCTACCACCCGGGTGGTAAGTTGCAGACATTATTAAAAACCTCATCAGAAGCTTGTTCAAAAATTTCAATACTCGAAACCTACCACCTGCGTCCCCTATTATTTACTACTACTAATAATAGCAGTATAATTGATCTGA'


def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

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

pattern = "ATGATCAAG"
reverse_pattern = "CTTGATCAT"
count_1 = PatternCount(Text=Text, Pattern=pattern)
count_2 = PatternCount(Text=Text, Pattern=reverse_pattern)
print('number of occurrences of pattern {}: {}\n'.format(pattern, count_1))
print('number of occurrences of pattern {}: {}\n'.format(reverse_pattern, count_2))


# find the following six 9-mers appear in this region three or more times
frequency_map = FrequencyMap(Text=Text, k=len(pattern))
print('Frequency Map: {}\n'.format(frequency_map))
res = [x for x in frequency_map if frequency_map[x] >= 3]
print('Patterns that appear 3 or more times: {}\n'.format(res))


#### The ClumpFinding Problem
"""
slide a window of fixed length L along the genome, looking for a region where a k-mer appears several times in short succession
typically L = 500 which reflects the typical length of ori in bacterial genomes

More formally, given integers L and t, a k-mer Pattern forms an (L, t)-clump inside a (longer) string Genome
if there is an interval of Genome of length L in which this k-mer appears at least t times.

Pattern is not known in advance!!! we must find it!

Output: All k-mers forming (L, t)-clumps in genome

"""

print('\nThe ClumpFinding Problem')
print('------------------------')

genome = str.lower('gatcagcataagggtccCTGCAATGCATGACAAGCCTGCAGTtgttttac')  # len 50 chars
k = 4 # 4-mer
L = 25 # # genome interva / slice width
t = 3  # appears at least 3 times

def sliceGenome(genome, L):
    regions = []
    for i in range(len(genome)):
        slice = genome[i:i+L]
        if len(slice) == 25:
            regions.append(slice)
        else:
            break
    return regions


# get region slices
regions  = sliceGenome(genome=genome, L=L)

# find clumps
clumps = {}
for region in regions:
    frequency_map = FrequencyMap(Text=region, k=k)
    for key in frequency_map:
        if frequency_map[key] >= t:
            clumps[key] = frequency_map[key]
print('({},{})-clumps: {} in genome\n'.format(L, t, clumps))


