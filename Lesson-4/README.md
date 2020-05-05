# Binder for BVCN FeGenie lesson

Initially forked from [here](https://github.com/binder-examples/conda). Thank you to the awesome [binder](https://mybinder.org/) team!

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Arkadiy-Garber/bvcn-binder-FeGenie/master?urlpath=lab)

Part of the [Bioinformatics Virtual Coordination Network](https://biovcnet.github.io/) :)



## Walkthrough

Enter the main FeGenie directory

    cd FeGenie

print the FeGenie help menu

    FeGenie -h

run FeGenie on test dataset

    FeGenie.py -bin_dir genomes/ -bin_ext fna -out fegenie_out

Go into the output directory and check out the output files

    cd fegenie_out
    less FeGenie-geneSummary-clusters.csv

run FeGenie on gene calls

    FeGenie.py -bin_dir ORFs/ -bin_ext faa -out fegenie_out --orfs

run FeGenie on gene calls, and use reference database (RefSeq sub-sample) for cross-validation

    FeGenie.py -bin_dir ORFs/ -bin_ext faa -out fegenie_out --orfs -ref refseq_db/refseq_nr.sample.faa

