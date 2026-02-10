# Folder allModels

Contains an *xml* file for the GEM of each of the 187 oryzae strains that were included in this study.

# Folder biolog

Fodler containing all data (raw, processed, visualizations) of the BioLog growth experiment, the results of our FBA analyses, and a comparison between our simulations and the experimental data.

# Folder clades

Contains the metadata that were gathered for all 187 oryzae strains in this study. Additionally, contains coordinates for each of these strains in reduced dimension space based on the reaction content of the different strains.

# Folder fasta

Fasta sequences were obtained and modified as follows:

## A. fumigatus

Download: https://nribf21.nrib.go.jp/CAoGD/search.cgi?prj=01909&sobj=cds&sgnm=5 > Download Data > Protein Sequences > Download
Modify: sed 's/|.*//' "$input_file" > "$output_file" # linux, not windows

## A. oryzae

Download: https://nribf21.nrib.go.jp/CAoGD/search.cgi?prj=01909&sobj=cds&sgnm=446 > Download Data > Protein Sequences > Download
Modify: sed 's/|.*//' "$input_file" > "$output_file" # linux, not windows

## A. niger

Download https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Aspni7# > Mycocosm > Annotation > 
All models, Filtered and Not > Proteins > Aspni7_all_proteins_20131226.aa.fasta.gz
Modify1: For consistency, change filetype from .fasta to .faa
Modify2: grep "\S" "$input_file" > "$output_file" # Remove empty lines
Modify3: awk -F"|" '{ print (NF==1) ? $1 : ">ANIG_"$3 }' "$input_file" > "$output_file" # format identifiers
Modify4: sed -e 's/\*$//' "$input_file" > "$output_file" # remove trailing asterisk

## A. oryzae isolates
Cleaned data obtained from Casper. I just changed the extension from .fas to .faa
Modify1: ren *.* *.faa # probably only works on windows

# Folder figures

Contains the figures that are code-generated and included in our manuscript or in its supplmentary materials.

# Folder genome

- INPUT_all_original.fasta: Output from BPGA run. This file contains all the (protein) sequences used as input for BPGA. Used as input for genes.Rmd
- u_clusters.txt: Output from BPGA run. This file contains the clustering information for each (protein) sequence in the BPGA input. Used as input for genes.Rmd
- BPGA2ortho_GEM_custom.csv: Output from genes.Rmd. This file links the original protein sequence identifiers to BPGA clusters. Used as input for construct_drafts.m

# Folder growthProfiler

All data and visualizations related to the growth experiment performed with the growth profiler.

# Folder intermediate

Folder containing helper files that take a while to construct. Some scripts will read input from this folder in stead of running time consuming steps. Not intended for end users. 

