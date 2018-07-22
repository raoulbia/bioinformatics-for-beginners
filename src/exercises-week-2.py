def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count

def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    print('ExtendedGenome: {}'.format(ExtendedGenome))
    for i in range(n):
        current_slice = ExtendedGenome[i:i + (n // 2)]
        count = PatternCount(current_slice, symbol)
        print("nbr of {}'s in {}: {}".format(symbol, current_slice, count))
        array[i] = count
    return array

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(Genome[0:n//2], symbol)
    print('\nGenome:{} \nsymbol count for first half {}: {}'.format(Genome, Genome[0:n//2], array))

    # must start at 1 otherwise the first loop will do array[0-1]
    for i in range(1, n):

        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1

        print(ExtendedGenome[i-1: i+(n//2)-1],
              ExtendedGenome[i-1],
              ExtendedGenome[i+(n//2)-1], array)
    return array

text = 'AAAAGGGG'
amino_acid = 'A'
res_1 = SymbolArray(Genome=text, symbol=amino_acid)
print('number of occurrences of {} that we encounter in each window of '
      'ExtendedGenome: {}'.format(amino_acid, res_1))

res_2 = FasterSymbolArray(Genome=text, symbol=amino_acid)
print('number of occurrences of {} that we encounter in each window of '
      'ExtendedGenome: {}'.format(amino_acid, res_2))