# Hands-on Tutorial - HMMER

[![Binder](http://mybinder.org/badge_logo.svg)](http://mybinder.org/v2/gh/biovcnet/functional-annotation-lesson-3-binder/master?urlpath=lab)

Using Binder to populate a terminal environment that contains the alignment tools [HMMER](http://hmmer.org/).

---

## Organism and Genes of Interest

We recently finished annotating the Bryozoan Bugula neritina, which is hundreds of millions of years old, and very difficult to relate to other known genes in invertebrates.  

Here is a picture of what a Bryozoan looks like.  They colonize rocks in the marine environment.

![Bugula neritina](https://upload.wikimedia.org/wikipedia/commons/d/d4/Bugula_neritina_%28YPM_IZ_101969%29_002.jpeg) 

---

## Build a position-sensitive model of a gene of interest

Use MyD88: innate immune signal transduction adaptor

Innate immunity is a pathogenic defense mechanism in invertebrates.  The genes can be highly conserved. MyD88 is a known signal transduction adaptor for toll-like receptors that triggers a pathogen-associated molecular pattern  of response in other organisms.  Innate immunity genes will allow organisms to recognize their invaders and fend them off or allow them to cohabitate.

Bryozoans have a population of bacterial symbionts, much like the relationship between zooxanthellae and coral.  These microbial organisms are important for fighting cancer, as they produce a unique set of compounds known as [bryostatins.](https://www.ncbi.nlm.nih.gov/pubmed/24033119)

Understanding how the bryozoan allows these symbionts to colonize instead of fending them off is an interesting puzzle to solve.

The starting point in this tutorial is a set of protein models of MyD88 from other organisms.

- protein sequences collected from UniProt (we will cover this in an upcoming tutorial)

- [Here is a paper that describes more about innate immunity genes such as MyD88.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4109969/)


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

To exit from the "less" command, push lower case "q" on your keyboard.

#### Spoiler Alert on what you should see

<details><summary>CLICK ME</summary>
<p>
#### The file is empty
</p>
</details>

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

There are two key parts of the alignment that suggest it might not be great quality - the string a gaps in most sequences at the start and end of the alignment. We can removed those sequences and re-run muscle.

```bash
rm MyD88.msa
nano MyD88.faa
```

Then find the sequences `tr|M4I212|M4I212_CRAGI` and `tr|A0A2B4RVP2|A0A2B4RVP2_STYPI` and use `Ctrl+k` to cut the lines (header + sequence) for both sequences. Save and close. Re-run the multiple sequence alignment.

```bash
muscle -in MyD88.faa -out MyD88.msa
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

Compare hmmscan to hmmsearch

```bash
hmmsearch --tblout Bugula.MyD88.tblout MyD88.hmm Bugula.pep
less Bugula.MyD88.tblout
```

Are there any potential matches to our hmm model?

---

## Predict multiple genes simultaneously

Building HMM models for specific genes of interest can be time consuming. In this example, the sequences were collected and screened prior to the start of the tutorial. There are numerous HMM databases to switch to which can be used to annotate multiple proteins at once. Specific HMM databases will be covered in future lessons, but one example is [Pfam](https://pfam.xfam.org/)

(myBinder does not want to download the full Pfam database, so we'll skip this part of the example)

`hmmsearch --domtblout Bugula.Pfam.domtblout Pfam-A.hmm Bugula.pep`


Things to consider:

    1.) How much does the reference database influence our ability to annotate a genome?
    2.) Which approach was more sensitive - local alignments with BLAST or hmm profile search with HMMER?
    3.) What roles do these approaches play in a robust annotation pipeline?
