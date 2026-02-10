# Aspergillus pan-genome reconstruction workflow

In this project, we reconstructed a pan-genome for Aspergillus, based on previously established genome-scale metabolic models (GEMs) for (1) Aspergillus niger, (2) Aspergillus oryzae and (3) Aspergillus fumigatus. The different steps of the workflow are detailed below. Note that large files may be unavailable from this GitHub repository. For these files, we refer to the associated Zenodo repository.

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

## Constructing the pan-oryzae GEM

Constructing the pan-oryzae GEM involves the scripts:

1. A set of MATLAB helper functions
	- cleanGrRules_local.m
	- getGenesFromGrRules.m
	- getModelFromOrthology_local.m
	- id2clust.m
	- replaceGrRules_local.m
2. m1_construct_panAspGEM.mlx

The script m1_construct_panAspGEM.mlx requires inputs generated in the R scripts:

3. r1_genes.rmd
4. r2_metabolites_oryzae.rmd
5. r3_metabolites_Ani_Afu.rmd

These inputs can directly be obtained from the folder `data/intermediate`. The output of this script is the file `data/intermediate/panAsp_v1_187_50_1385.xml`, which serves as a first draft of the pan-oryzae GEM.

## Performing FBA simulations to refine the pan-oryzae GEM

Refining the pan-oryzae GEM is initially achieved by performing FBA simulations, and identifying gaps and errors in the GEM by comparing FBA predictions with observed oryzae growth profiles obtained from BioLog experiments.

Refining the GEM and performing FBA simulations is performed in the scripts:

1. m2_curate_panAspGEM.mlx: initial curation to fix reaction reversibility and handling of protons and water in reactions, in addition to fixing a few clear errors.
2. p1_map_kegg_chebi.ipyb: needed for p2_constructEnsemble.ipynb
3. helpers.py: needed for p2_constructEnsemble.ipynb
4. p2_constructEnsemble.ipynb: script to wrangle the models of the different oryzae strains into a single data object for more efficient downstream analysis
5. r4_BioLog_prepare.rmd: In this script, we prepare the outputs of the BioLog phenotyping experiment to serve as input for the FBA analyses performed in *simulations.ipynb*. In addition, we define a set of core metabolites based on the BioLog experimental data.
6. p3_simulations.ipynb: In this script, we perform FBA analyses for eight different oryzae isolates that were used in *in vivo* experiments where growth on different sources of carbon and nitrogen was assessed. The goal of this script is to further curate/gapfill these GEMs, such that FBA simulations on growth versus no growth align with the observed growth phenotypes of the experiment.

## Downstream analyses and visualizations

1. m3_visualizeModels.mlx: visualizes the 187 strain-specific GEMs in a tSNE plot, used for manuscript Figure 2.
2. r5_Visualizations_for_manuscript.rmd: Visualizes summary statistics displayed in Figures 1B and 1D. Also generates heatmaps for the BioLog experiments displayed in Figures 3 and S3.
3. p4_compareModels.ipynb: computes summary statistics for the pan-oryzae GEM, used throughout the manuscript text and in Figures 1C and 2A.

