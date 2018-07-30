#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from src.wk4_gibbs_sampling import *

def test_noralize():
    probas = {'A': 0.1, 'C': 0.1, 'G': 0.1, 'T': 0.1}
    res = Normalize(probas)
    print(res)