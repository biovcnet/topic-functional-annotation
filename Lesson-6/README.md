# Binder for BVCN Functional Annotation lesson 6

Initially forked from [here](https://github.com/binder-examples/conda). Thank you to the awesome [binder](https://mybinder.org/) team!

[![Binder](https://mybinder.org/badge_logo.svg)](https://gesis.mybinder.org/binder/v2/gh/Arkadiy-Garber/bvcn-binder-kegg-koala/master?urlpath=lab)

Part of the [Bioinformatics Virtual Coordination Network](https://biovcnet.github.io/) :)

Vist the [KEGG-Decoder GitHub](https://github.com/bjtully/BioData/tree/master/KEGGDecoder) page for more information


## Walkthrough

There are 4 example files from the various KEGG KOALA outputs, `blastkoala.txt`, `ghostkoala.txt`, `kofamkoala.txt`, and `NORP_subset.txt`. Each one has the data formatted in a specific way that may (or may not) work with KEGG-Decoder.

`blastkoala.txt` is a single genome formated correctly

```
BAFMBKGE_00001  K03596
BAFMBKGE_00002  K03100
BAFMBKGE_00003  K03685
```

As a single genome, KEGG-Decoder can only produce a `static` visualization output - `interactive` visualization requires at least 3 genomes.

```
KEGG-decoder -i blastkoala.txt -o blastkoala.ko -v static
```

The `blastkoala.svg` file contains a heatmap comparing the genome results to the KEGG-Decoder pathways. The more `red` the cell of the heatmap, the more complete the pathway

`ghostkoala.txt` has a single genome that is formatted incorrectly

```
NC_015736.1_1   K01703
NC_015736.1_2   K01704
```

KEGG-Decoder can only accept gene IDs with a single `_`. KEGG-Decoder uses all of the string prior to the first `_` as the name to associate with each row of the output and heatmap. Identical strings would result in the merging of multiple genomes together

```
sed 's/NC_/NC/g' ghostkoala.txt > ghostkoala-mod.txt
KEGG-decoder -i ghostkoala-mod.txt -o ghostkoala.ko -v static
```

`kofamkoala.txt` is a complex table that requires some pre-processing

```
sed 's/ \+ /\t/g' kofamkoala.txt | cut -f1,2 | sed 's/\* //g' |  grep -v "#" > kofamkoala-mod.txt
KEGG-decoder -i kofamkoala-mod.txt -o kofamkoala.ko -v static
```

The genome name for the output will be `WP` as it is the string prior to the first `_`.

`NORP_subset.txt` contains KOALA results for 4 genomes. The `NORP.ko` file will be created after each KEGG-Decoder run, but will contain identical contents. It is the numerical values used to the create the heatmap. Because there are 3+ genomes, we can make an `interative` visualization as well as a `static` visualization.

```
KEGG-decoder -i NORP_subset.txt -o NORP.ko -v static
KEGG-decoder -i NORP_subset.txt -o NORP.ko -v interactive
```
