# Binder for BVCN MagicLamp lesson

Initially forked from [here](https://github.com/binder-examples/conda). Thank you to the awesome [binder](https://mybinder.org/) team!

[![Binder](https://mybinder.org/badge_logo.svg)](https://gesis.mybinder.org/binder/v2/gh/Arkadiy-Garber/bvcn-binder-magiclamp/master?urlpath=lab)

Part of the [Bioinformatics Virtual Coordination Network](https://biovcnet.github.io/) :)


## Walkthrough

Enter the MagicCave

    cd MagicCave/

print the MagicLamp help menu

    MagicLamp.py help

print WspGenie help menu

    MagicLamp.py WspGenie -h

run WspGenie on test dataset

    MagicLamp.py WspGenie -bin_dir test_dataset/ -bin_ext fna -out wspgenie_out


go into the wspgenie output directory and check out the output file

    cd wspgenie_out/
    less -S wspgenie-summary.csv

check out the gene predictions

    cd ORF_calls/
    cd ../../

mv ORF calls to the main directory

    mv wspgenie_out/ORF_calls/ ./

print LithoGenie help menu

    MagicLamp.py LithoGenie -h

run LithoGenie on ORF calls

    MagicLamp.py LithoGenie -bin_dir ORF_calls/ -bin_ext faa --orfs -out lithogenie_out

check out the output

    cd lithogenie_out/
    less -S lithogenie-summary.csv
    less lithogenie.ALL.heatmap.csv
    cd ../

re-run LithoGenie to create a .heatmap.csv for an element-of-interest

    MagicLamp.py LithoGenie -bin_dir ORF_calls/ -bin_ext faa --orfs -out lithogenie_out --skip -cat sulfur
    # answer 'y' to the question
    MagicLamp.py LithoGenie -bin_dir ORF_calls/ -bin_ext faa --orfs -out lithogenie_out --skip -cat iron

check out the new results

    cd lithogenie_out/
    less lithogenie.sulfur.heatmap.csv
    less lithogenie.iron.heatmap.csv

print the HmmGenie help menu

    MagicLamp.py HmmGenie -h

run HmmGenie with a set of HMMs for gas vesicle formation

    MagicLamp.py HmmGenie -hmm_dir MagicCave/hmms/gas/ -hmm_ext hmm -bin_dir test_dataset/ -bin_ext fna -out gas_out

check out the results and re-run HmmGenie with more stringent parameters

    MagicLamp.py HmmGenie -hmm_dir MagicCave/hmms/gas/ -hmm_ext hmm -bin_dir test_dataset/ -bin_ext fna -out gas_out -clu 5

check out the results

    cd gas_out/
    less -S genie-summary.csv










# Binder for BVCN Functional Annotation lesson

Initially forked from [here](https://github.com/binder-examples/conda). Thank you to the awesome [binder](https://mybinder.org/) team!

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Arkadiy-Garber/bvcn-binder-hmmgenie/master?urlpath=lab)

Part of the [Bioinformatics Virtual Coordination Network](https://biovcnet.github.io/) :)

## Walkthrough

Run phmmer on the Cyc1 fasta file

    phmmer -A Cyc1.refseq.msa --tblout Cyc1.refseq.tblout -E 1E-20 Cyc1.faa ../refseq_db/refseq_nr.sample.faa

Build HMM file from MSA (multiple sequence alignment) file, using hmmbuild

    hmmbuild Cyc1.hmm Cyc1.refseq.msa

Query the Cyc1 HMM file against refseq database sample

    hmmsearch --tblout Cyc1.hmm.refseq.tblout Cyc1.hmm ../refseq_db/refseq_nr.sample.faa

Examine the output file. What do the bit scores look like for likely false positives

    less Cyc1.hmm.refseq.tblout

Move into directory containing MtrA FASTA file, and create an alignment using Muscle.

    muscle -in MtrA.faa -out MtrA.fa

Build HMM file from MSA (multiple sequence alignment) file, using hmmbuild

    hmmbuild MtrA.hmm MtrA.fa

Query the MtrA HMM file against refseq database sample

    hmmsearch --tblout MtrA.hmm.nr.tblout MtrA.hmm ../refseq_db/refseq_nr.sample.faa

Examine the output file. What do the bit scores look like for likely false positives

    less MtrA.hmm.nr.tblout

Move the HMM files into a single directory

    mv MtrA.hmm ../HMMs/
    mv Cyc1.hmm ../HMMs/

Check out the Pfam-derived HMM and bitscores.txt file

    less Catalase.hmm

Run HmmGenie (MagicLamp) on test dataset using the new HMM collection

    MagicLamp.py HmmGenie -hmm_dir HMMs/ -hmm_ext hmm -bin_dir test_data/ -bin_ext txt -out hmmgenie_out -eval 1E-1
    MagicLamp.py HmmGenie -hmm_dir HMMs/ -hmm_ext hmm -bin_dir test_data/ -bin_ext txt -out hmmgenie_out -bit HMMs/bitscores.txt
    MagicLamp.py HmmGenie -hmm_dir HMMs/ -hmm_ext hmm -bin_dir test_data/ -bin_ext txt -out hmmgenie_out -bit HMMs/bitscores.txt -clu 2




