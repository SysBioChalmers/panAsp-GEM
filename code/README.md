# Aspergillus pan-genome reconstruction workflow

In this project, we reconstructed a pan-genome for Aspergillus, based on previously established genome-scale metabolic models (GEMs) for
(1) Aspergillus niger, (2) Aspergillus oryzae and (3) Aspergillus fumigatus. The different steps of the workflow are detailed below.

## Obtain template GEMs and fasta sequences

First, we obtain fasta protein sequences and GEMs for Aspergillus niger, oryzae, and fumigatus. 

- The GEMs were obtained as detailed in `model/templates/README.md`
- The fasta sequence files were obtained and processed as detailed in `data/README.md`

## BPGA/USEARCH analysis

The first step in our workflow for generating a pan-genome for Aspergillus is running the BPGA algorithm by 
[Chaudhari et al. (2016)](https://www.nature.com/articles/srep24373). The following steps are taken:

- BPGA was downloaded from https://sourceforge.net/projects/bpgatool/files/
- BPGA step 1 was run using all sequences from `data/fasta/` with option 4 "Use any protein fasta files"
- BPGA step 2 was run with option 1 "Use usearch clustering algorithm" and a sequence identity cutoff of 0.5
- From the output, the files `u_clusters.txt` and `INPUT_all_original.fasta` were retained and stored in
the `data/genome/` folder.

Note that running BPGA with these settings is equivalent to running [USEARCH](http://drive5.com/usearch/)
with default settings and a sequence identity cutoff of 0.5.

## genes.rmd
## construct_panAspGEM.mlx
## metabolites_niger.rmd
## metabolites_niger_oryzae.rmd
## metabolites_oryzae.rmd
## metabolites_Ani_Afu.rmd
## curate_panAsp.mlx

## (curate_panaAsp.rmd and metabolites_panAspv1_OtherAsp.Rmd)



