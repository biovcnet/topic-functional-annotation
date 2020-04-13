# TOPIC: Functional Annotation
A repository for the lessons and tutorials for the Functional Annotation TOPIC channel of the [BVCN](https://biovcnet.github.io/)


## Prerequisites
* [A working Unix environment](https://github.com/biovcnet/biovcnet.github.io/wiki/1.-Setting-up-a-local-Linux-(or-Unix)-environment)
* [Experience with the command line](https://github.com/biovcnet/biovcnet.github.io/wiki/2.-Using-the-Command-line)

# Overview
This BVCN topic will cover:

* how to predict opening reading frames for all three domains of life
* the methodologies used to assign functions to proteins
* the intricacies of the specific tools and databases that can used for varying levels of specificity

# Lesson 1
### Title: How to predict open reading frames on DNA?
Goals

* Tools for predicting open reading frames for Bacteria and Archaea
* Tools for predicting open reading frames for Eukaryotes
* Converting open reading frames to predicted proteins

[Watch the lesson](https://youtu.be/uGjjN-q7N2E)

[Watch the tutorial](https://youtu.be/on2fZveY8sU)

[Follow the tutorial](https://github.com/biovcnet/topic-functional-annotation/blob/master/Lesson-1/README.md)

[Access the presentation](https://github.com/biovcnet/topic-functional-annotation/blob/master/Lesson-1/Lesson-1.pdf)

### Lesson 1 FAQ
Q. When you talk about multiples N and use the command -m, that N are genomes or contings/clusters?

A. The Ns refer to a run of base pairs in a DNA sequence for which it is unclear which bp is actually there. A sequence like this would be fine ATGCANGTCAGA, but sequence like ATGCCCNNNNNNNNNNNNNNGTACAC would be thrown out. The reason being that some assemblers will stitch together sequences using Ns to represent a gap and that gap could be hundreds of bps long. If Prodigal ran across that gap it runs the potential of creating a chimeric protein.

Q. The -p function when you decide to use single or meta, why you would use one over another?

A. [From the Prodgial GitHub issues](https://github.com/hyattpd/Prodigal/issues/57#issuecomment-536608100): You always want "-p single" when you know everything belongs to one genome. So, yes, "-p single" (the default, don't need to specify this) on binned contigs. "-p meta" should be used on short sequences (when insufficient data exists to train on), and on any mixed samples (metagenomic assemblies, etc.)

Q. The input is raw data - such as a sequence obtained but yet processed? 

A. The input would be any bacterial or archaeal genome that you have the untranslated DNA for but are looking to generate the proteins.

# Lesson 2
### Title: Orthology based functional annotation
Goals

* Difference between homolog and ortholog?
* Tools for orthology based functional annotation
    * BLAST+
    * DIAMOND
* Databases
* Interpreting similarity based results

[Watch the tutorial](https://www.youtube.com/watch?v=oHg5SJYRHA0)

Follow the tutorial here

Access the presentation [here](https://github.com/biovcnet/topic-functional-annotation/blob/master/Lesson-2/Presentation1.pdf)

