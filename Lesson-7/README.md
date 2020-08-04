# Binder for BVCN AntiSMASH lesson

Initially forked from [here](https://github.com/binder-examples/conda). Thank you to the awesome [binder](https://mybinder.org/) team!

[![Binder](https://mybinder.org/badge_logo.svg)](https://gesis.mybinder.org/binder/v2/gh/Arkadiy-Garber/bvcn-binder-antismash/master?urlpath=lab)

Part of the [Bioinformatics Virtual Coordination Network](https://biovcnet.github.io/) :)


## Walkthrough

Put bioinf_tools into PATH

    export PATH="${PATH}:/home/jovyan/MagicCave:/home/jovyan/bioinf_tools"

Download antismash databases

    download-antismash-databases

Print the antismash help menu

    antismash -h

Run the simplest version of the command

    antismash Sodalis_praecaptivus-fixed.gbk

Run the simplest version of the command with an assembly

    antismash --genefinding-tool prodigal Sodalis_praecaptivus-fixed.fasta

 Convert GenBank output file to FASTA format

    bit-genbank-to-fasta -i NZ_CP006569.1.region001.gbk -o NZ_CP006569.1.region001.fasta
    bit-genbank-to-AA-seqs -i NZ_CP006569.1.region001.gbk -f NZ_CP006569.1.region001.faa
    
Print full help menu

    antismash --help-showall

Run antismash with additional analyses

    antismash --genefinding-tool prodigal --fullhmmer --pfam2go --asf --cb-knownclusters --clusterhmmer --cf-create-clusters Sodalis_praecaptivus-fixed.fasta



