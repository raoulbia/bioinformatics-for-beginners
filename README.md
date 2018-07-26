#### Bioinformatics for Beginners

Repository for scripts used in the Coursera course 
[Biology Meets Programming: Bioinformatics for Beginners](https://www.coursera.org/learn/bioinformatics/home/welcome).

#### Glossary

**Motifs**

Are sets of k-mers. they represent potential binding sites i.e. the hidden messages in DNA.

**Consensus Motif**

A consensus motif is selected from a set of motifs by looking for the most frequently occurring 
nucleotides in each column of the motif matrix (ties are broken arbitrarily).

If we select Motifs correctly from the collection of upstream regions,
then the corresponding consensus motif provides a candidate regulatory motif for these regions.
    
**Profile Matrix**

A profile matrix is constructed from a collection of k-mers Motifs in Dna. It is based on 
frequency of occurrence of ACGT bases in a set of motifs (k-mers).

As we slide across a text string (DNA) using a window of length k, the profile matrix is used 
to  compute a probability score of the current k-long slice using the profile matrix. At the end of the 
sliding (for a given text string), the k-mer with the highest probability is selected. This 
process is repeated for each text string.
 

#### Notes

certain genes are able to control the transcription of other genes because the regulatory proteins that
they encode are transcription factors, or master regulatory proteins that turn other genes on and off.

A transcription factor regulates a gene by binding to a specific short DNA interval called a regulatory motif,
or transcription factor binding site, in the gene’s upstream region.

A gene’s upstream region is a 600-1000 nucleotide-long region preceding the start of the gene.

For example, CCA1 expresses master regulatory proteins that bind to "AAAAAATCT" in the upstream region of
many genes regulated by CCA1.

this is a good blog post explaining how the greedy algorithm works:
http://www.mrgraeme.co.uk/greedy-motif-search/

The greedy motif search process relies on the following nested subroutines:

```
# initial Motifs = first k letters substring from DNA[0]
|--- get_profile_matrix(Motifs) 
|    |-- get_count_matrix(Motifs)

# intital profile matrix is used to get a most-probable-k-mer for each DNA string
|--- get_most_probable_kmer_from_profile_matrix((Text, profile))
|    |-- get_string_probability_from_profile_matrix

# at each iteration we need to get a motif score and see if current score is improved
|--- get_consensus_motif_score(Motifs)
|    |--get_consensus_motif(Motifs)
|       |-- get_count_matrix(Motifs)
```

