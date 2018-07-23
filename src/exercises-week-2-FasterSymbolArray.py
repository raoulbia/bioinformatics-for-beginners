#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import decimal
from timeit import default_timer as timer
ts = timer()

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
    print('Extended Genome: {}'.format(ExtendedGenome))
    for i in range(n):
        current_slice = ExtendedGenome[i:i + (n // 2)]
        count = PatternCount(current_slice, symbol)
        array[i] = count
    return array

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(Genome[0:n//2], symbol)

    # must start at 1 otherwise the first loop will do array[0-1]
    for i in range(1, n):

        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1

    return array

text = 'AAAAGGGG'
amino_acid = 'A'

start = timer()
res_1 = SymbolArray(Genome=text, symbol=amino_acid)
print('number of occurrences of {} that we encounter in each window of '
      'ExtendedGenome: {}'.format(amino_acid, res_1))
time = decimal.Decimal(timer() - start)
print('timer: {}'.format(round(time,6)))

start = timer()
res_2 = FasterSymbolArray(Genome=text, symbol=amino_acid)
print('Fast: number of occurrences of {} that we encounter in each window of '
      'ExtendedGenome: {}'.format(amino_acid, res_2))
time = decimal.Decimal(timer() - start)
print('timer: {}'.format(round(time,6)))