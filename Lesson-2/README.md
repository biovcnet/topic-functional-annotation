# Hands-on Tutorial - BLAST & DIAMOND

[![Binder](http://mybinder.org/badge_logo.svg)](http://mybinder.org/v2/gh/biovcnet/functional-annotation-lesson-2-binder/master?urlpath=lab)

Using Binder to populate a terminal environment that contains the alignment tools [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) and [DIAMOND](https://diamond.readthedocs.io/en/latest/).

## Set-up

Create a directory and move our data into that directory.

```bash
mkdir blast-example
mv *faa ./blast-example
cd blast-example
```
We will be working with 2 files:

* `Cyanobacteria-TMED155.orfs.faa`
* `carbon-fixation-markers.faa`

`Cyanobacteria-TMED155.orfs.faa` contains 1,585 predicted proteins from a cyanobacteria metagenome-assembled genome (MAG) reconstructed from the Mediterranean Sea

As determined by CheckM - the MAG TMED155 is 64.63% complete with 1.29% redundancy.

`carbon-fixation-markers.faa` contains 326 proteins chosen to represent key marker genes in the 5 established carbon fixation pathways - Calvin-Benson-Bassham, reverse TCA, Wood-Ljungdahl, 3-hydroxypropionate bicycle, & 3-hydroxypropionate/4-hydroxybutyrate cycle.

## BLAST Walkthrough

Create a BLAST index of the 'subject' sequences. In this case, the subject are the sequences we are comparing to the genome.

```bash
makeblastdb -in carbon-fixation-markers.faa -dbtype prot
```

Creates 3 index files that end in `*phr`, `*pin`, `*psq`

To compare the TMED155 proteins against the carbon fixation gene database, run a `BLASTP` command.

```bash
blastp -query Cyanobacteria-TMED155.orfs.faa -db carbon-fixation-markers.faa -out BLAST_fulloutput.txt -evalue 1e-20 -num_descriptions 5 -num_alignments 5
```

Let's check the contents of the output

```bash
less +11093 BLAST_fulloutput.txt
```

This format is very similar to the standard BLAST web interface, which might be easy for a person to understand, but is terrible for a machine.

We will use the tabular out format option to force BLAST it give us a tab-delimited table. And we will select which values we want to see and in what order. The full list of options can be found [here](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6).

```bash
blastp -query Cyanobacteria-TMED155.orfs.faa -db carbon-fixation-markers.faa -out BLAST_output.tab -evalue 1e-20 -max_target_seqs 10 -outfmt '6 qseqid qstart qend sseqid slen sstart send bitscore pident evalue'
less BLAST_output.tab
```

The output contains a lot of poor BLAST HSP (high scoring pairs). We can discern this by the low percent identity values (2nd to last column). We can use Unix to only select for matches with high percent identity.

```bash
awk '{if ($9>=50) print }' BLAST_output.tab
```

And sort that output

```bash
awk '{if ($9>=50) print }' BLAST_output.tab | sort -nrk 9,9
```

Unfortunately, the BLAST output only captures the subject ID, but the descriptions of the `carbon-fixation-markers.faa` contains functional annotations.

We can grab all of them at once with something like...

```bash
awk '{if ($9>=50) print }' BLAST_output.tab | cut -f1,4 | sort -k 2 > TMED155-BLAST-matches.tmp
cut -f2 TMED155-BLAST-matches.tmp > gene-matches.tmp
grep -f gene-matches.tmp carbon-fixation-markers.faa | sed 's/>//' | sort > gene-matches-descriptions.tmp
paste TMED155-BLAST-matches.tmp gene-matches-descriptions.tmp | sort > TMED155-BLAST-matches.tsv
rm *.tmp
```

## DIAMOND Walkthrough

We will follow a similar series of steps to perform a DIAMOND BLAST search.

Make a DIAMOND index of the subject sequences

```bash
diamond makedb --in carbon-fixation-markers.faa -d carbon-fixation-markers
```

And then perform a DIAMOND BLASTP search - in `fast` mode

```bash
diamond blastp -p 1 -q Cyanobacteria-TMED155.orfs.faa -d carbon-fixation-markers.dmnd -o TMED155-diamond-fast.tab --max-target-seqs 15 -f 6 qseqid qstart qend sseqid slen sstart send bitscore pident evalue
```
Perform a DIAMOND BLASTP search - in `more-sensitive` mode

```bash
diamond blastp --more-sensitive -p 1 -q Cyanobacteria-TMED155.orfs.faa -d carbon-fixation-markers.dmnd -o TMED155-diamond-more-sensitive.tab --max-target-seqs 15 -f 6 qseqid qstart qend sseqid slen sstart send bitscore pident evalue
```

Notice the longer compute time and then compare the outputs

```bash
wc -l TMED155-diamond-fast.tab
wc -l TMED155-diamond-more-sensitive.tab

awk '{if ($9>=50) print }' TMED155-diamond-fast.tab | sort -nrk 9,9 | wc -l
awk '{if ($9>=50) print }' TMED155-diamond-more-sensitive.tab | sort -nrk 9,9 | wc -l
```

`more-sensitive` detects more significant matches than `fast`. But for high percent identity matches (>50%) `more-sensitive` and `fast` have the same number of matches.

Compare the annotations as above:

```bash
awk '{if ($9>=50) print }' TMED155-diamond-more-sensitive.tab | cut -f1,4 | sort -k 2 > TMED155-diamond-matches.tmp
cut -f2 TMED155-diamond-matches.tmp > gene-matches.tmp
grep -f gene-matches.tmp carbon-fixation-markers.faa | sed 's/>//' | sort > gene-matches-descriptions.tmp
paste TMED155-diamond-matches.tmp gene-matches-descriptions.tmp | sort > TMED155-diamond-matches.tsv
rm *.tmp
```

Interestingly, DIAMOND BLAST detected an other signifcant match in TMED155 compared to BLAST for protein `120126_4` with similarity to `formate dehydrogenase beta subunit`. This does not change our interpretation of the carbon fixation potential for TMED155 as formate dehydrogenase can be part of multiple cellular pathways.

