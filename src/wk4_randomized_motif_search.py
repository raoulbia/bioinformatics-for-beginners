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

