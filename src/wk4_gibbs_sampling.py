#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def  Normalize(Probabilities):

    total = sum(Probabilities.values())

    for k, v in Probabilities.items():
        Probabilities[k] = v / total
    return Probabilities