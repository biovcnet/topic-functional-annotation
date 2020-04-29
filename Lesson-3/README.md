# Hands-on Tutorial - HMMER

[![Binder](http://mybinder.org/badge_logo.svg)](http://mybinder.org/v2/gh/biovcnet/functional-annotation-lesson-3-binder/master?urlpath=lab)

Using Binder to populate a terminal environment that contains the alignment tools [HMMER](http://hmmer.org/).

---

## Organism and Genes of Interest

[Bugula neritina](https://en.wikipedia.org/wiki/Bryozoa), "moss animals" 

---

## Build a position-sensitive model of a gene of interest

Use MyD88: innate immune signal transduction adaptor

- protein sequences collected from UniProt

- [What is MyD88?](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4109969/)


---

## Compare MyD88 proteins against Bugula neritina genome

We'll compare the ability to detect this protein using BLAST and HMMER

```bash
cd example_data/
ls
```

We have a set of 16 MyD88 proteins `MyD88.faa` and the CANU-filtered Bugula neritina proteome `Bugula.pep`

Make a BLAST database

```bash
makeblastdb -in Bugula.pep -dbtype prot
```

BLAST MyD88 proteins against the B. neritina genome to see if we get a result.

```bash
blastp -query MyD88.faa -db Bugula.pep -outfmt 6 -evalue 1e-5 -out blastp.MyD88.Bugula.pep.outfmt6
```
This should be a rather short blast operation, because we only have a few sequences.

***Did we get any results at all?***

```bash
less blastp.MyD88.Bugula.pep.outfmt6
```

### Let's try a probabilistic approach with HMMER (HMMER: biosequence analysis using profile hidden Markov models)

HMMER is a program that uses hidden Markov profiles to predict the probability that two proteins are related.

The first step is to give it a set of data for calculating probabilities.  

There are two common sources for reference data sets for HMMER:

1.) precomputed hmm profiles in public databases (e.g. Pfam)

2.) manually created multiple sequence alignments

We are going to do approach #2 first.

We have a set of data already that we can create our own hmm profile from, the MyD88 set of data.

## Align sequences

First we have to align the sequences. It is important to eyeball them at some point to make sure that they are not composed of non-overlapping, too long or too short sequences. It is better to be more conservative in the number of sequences rather than trying to include **everything**.  Some databases contain mislabeled sequences that can confound your search for proteins.

```bash
muscle -in MyD88.faa -out MyD88.msa
```

Take a look inside the file to see if the alignments make sense.  Kick out any sequences that are making the sequences less cohesive.

```bash
less -S MyD88.msa
```

Then, we are going to use hmmbuild to analyze the occurrence of each peptide in relation to the position in the aligned proteins.

```bash
hmmbuild MyD88.hmm MyD88.msa
```

The final step to create a mathematical matrix that shows the probability of each transition based on the data provided.

```bash
hmmpress MyD88.hmm
ls
```

The MyD88.htm is now a searchable hmm profile that we can use with our genome to find proteins which are probably related.

```bash
hmmscan --domtblout Bugula.MyD88.domtblout MyD88.hmm Bugula.pep
```

Look at the results in the Bugula.MyD88.domtblout

```bash
less Bugula.MyD88.domtblout
```

Are there any potential matches to our hmm model?

---

## Predict multiple genes simultaneously

Building HMM models for specific genes of interest can be time consuming. In this example, the sequences were collected and screened prior to the start of the tutorial. There were numerous HMM databases to switch to which can be used to annotate multiple proteins at once. Specific HMM databases will be covered in future lessons, but one example is [Pfam](https://pfam.xfam.org/)

(myBinder does not want to download the full Pfam database, so we'll skip this part of the example)

`hmmsearch --domtblout Bugula.Pfam.domtblout Pfam-A.hmm Bugula.pep`


Things to consider:

    1.) How much does the reference database influence our ability to annotate a genome?
    2.) Which approach was more sensitive - local alignments with BLAST or hmm profile search with HMMER?
    3.) What roles do these approaches play in a robust annotation pipeline?
