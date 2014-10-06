Fractionation
=============

A long and probably ill organized script used to analyze data from CoGe

This script was sued to generate some statisitcs and figures for the recent PNAS submission. All the code used in that analysis is here, although there is much that was not included in the final draft of the MS. Unfortunately, the script is long and over-burderned with analyses and could really do with being split up into smaller, more easily digeastable chunks.

The script requires a few bits of data to get going, those being:

[1] CoGe output "final syntenic gene output" file. This can be generated for at CoGe.org. this can be used to parse to format of [2] and [3] below.
[2] A matrix with rowname representing reference genome genes, and columns representing syntenic blocks, "-" added where there is no gene in the corresponding block, and the gene name written where there is a match.
[3] The same matrix as [2] expcet with the chromosome to which the "cotton" gene is mapped to.
[4] KaKs data generated using CoGe
[4] siRNA data, in this case the data are from Lei Gong
[5] TE annotation data, from the D5 genomne release
[6] The nearest TE to each gene, generated with BEDtools
[7] The GC content also generated with BEDtools
[8] Expression data generated used in REnny-Byfield et al. 2014 GBE.

An
