# Hands-on Tutorial - Prodigal

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/biovcnet/topic-functional-annotation/master?urlpath=lab)

Using Binder to populate a terminal environment that contains the protein prediction tool - Prodigal.

## Walkthrough

Move into the Lesson 1 directory

```
cd Lesson-1
```

'ls' the content

```
ls
```

A single FASTA file is present `MED921-Poseidoniales.fasta`

Generate the Prodigal help command

```
prodigal -h
```

We are going to compare the results running of Prodigal with the `single` and `meta` parameter options associated with the `-p` flag

```
prodigal -i MED921-Poseidoniales.fasta -a MED921-proteins-single.faa -d MED921-genes-single.fna -o prodigal-MED921-single.out -m -p single
```

and

```
prodigal -i MED921-Poseidoniales.fasta -a MED921-proteins-meta.faa -d MED921-genes-meta.fna -o prodigal-MED921-meta.out -m -p meta
```

Check the contents of the output file with `less`. You can refer the [Prodigal](https://github.com/hyattpd/prodigal/wiki/understanding-the-prodigal-output) help page to determine what the various outputs mean.

Compare the total number of predicted proteins

```
grep ">" MED921-proteins-single.faa | wc -l
```

and 

```
grep ">" MED921-proteins-meta.faa | wc -l
```

Are these differences meaningful? 
