## pAo: a consensus genome-scale metabolic model for Aspergillus oryzae

[![Version](https://badge.fury.io/gh/{{organization or username}}%2F{{repository name}}.svg)](https://badge.fury.io/gh/sysbiochalmers/yeast-gem)  
[![Zenodo](https://zenodo.org/badge/{{Zenodo ID}}.svg)](https://zenodo.org/badge/latestdoi/{{Zenodo ID}})  
[![Gitter chat](https://badges.gitter.im/{{organization or username}}/{{repository name}}.svg)](https://gitter.im/{{organization or username}}/{{repository name}})


#### Description

Aspergillus oryzae (koji mold) is a key microorganism in traditional food fermentations including soy sauce, sake, and miso, and is important in novel culinary applications and modern biotechnology, such as sustainable meat alternatives and enzyme production. Despite its industrial importance, until recently, the most recent genome-scale metabolic model (GEM) for A. oryzae dated back to 2008 and was limited to a single strain (RIB40). Here, we present pAo, a pan-GEM for A. oryzae, integrating genomic data from 187 strains to capture species-wide metabolic diversity. Our model comprises 2,018 reactions (a 52% increase over the RIB40-based model) and includes previously overlooked pathways, such as cytochrome P450-mediated xenobiotic metabolism and extended amino acid metabolism. Using this pan-GEM, we derived strain-specific GEMs and validated them through high-throughput phenotypic screening on 290 substrates. Growth experiments on industrial carbon sources (glucose, glycerol, maltose, and xylose) revealed significant metabolic diversity across strains. This resource enables informed strain selection for biotechnological applications and provides a foundation for future metabolic engineering in A. oryzae.

#### Citation

Gilis, J., van der Luijt, C.R.B., Feller, M., Sanchez-Giron Barba, C., Sommer, M.O.A., Jahn, L.J., Kerkhoven, E.J. (2026). pAo: a consensus genome-scale metabolic model for Aspergillus oryzae capturing intra-species diversity.

#### Keywords

> Keywords are be separated by semicolons.
> The `Model source` field contains the source(s) of the current model, eg existing GEMs. If possible, use the Markdown format to add the URL with the DOI. The (NCBI) taxonomy ID should be provided in the [format from identifiers.org](https://registry.identifiers.org/registry/taxonomy). For the genome identifier, please provide the ENA/GenBank/RefSeq identifier via *identifiers.org*, or from other sources such as PATRIC or KBase.  

**Utilisation:** {{ experimental data reconstruction; _in silico_ strain design; model template }}  
**Field:** {{ metabolic-network reconstruction }}  
**Type of model:** {{ reconstruction; curated }}  
**Model source:** {{ TODO }}  
**Taxonomic name:** {{ _Aspergillus oryzae_ }}  
**Taxonomy ID:** {{ [taxonomy:5062](https://identifiers.org/taxonomy:5062) }}  
**Genome ID:** {{ [insdc.gca:GCA_009687165.1](https://identifiers.org/insdc.gca:GCA_009687165.1)  }}  
**Metabolic system:** {{ full metabolism }}  
**Tissue:**  
**Bioreactor:**    
**Cell type:**  
**Cell line:**  
**Strain:** {{ 187 strains }}  
**Condition:** {{ aerobic; glucose-limited; nitrogen-limited }}  

### Repository structure

There are three main folders in this project. Each of these folders contains a README file with more details.

- code: contains all the source code required to reproduce all aspects of this project (data preprocessing, data analysis, data visualization) as well as the virtual environment used for running the python analyses.
- data: contains raw data, intermediate data files, and some final outputs (simulation results, figures).
- model: contains the final pan-oryzae model in several different formats.

### Usage of the final pan-oryzae model

The `model/` folder contains the consensus pan-model `pAo.xml`, the eight experimentally validated
strain models `pAo_<strain>.xml`, and two pickled collections (see `model/README.md` for details).

The pan-model can be loaded into MATLAB (RAVEN/COBRA) with:

```
pAo = importModel("./pAo.xml");
```

It can be loaded into python using:

```
data_dir = Path("./model")
data_dir = data_dir.resolve()
model_path = data_dir / "pAo.xml"
panOryzae = read_sbml_model(str(model_path.resolve()), skip_validation=True)
```

The strain-specific collections are distributed as pickled lists of cobra models:

```
from pickle import load
# the 8 validated strains + template + pan-model (with experimental gap-fills)
with open("./model/pAo_strain-GEMs_validated.pickle", 'rb') as infile:
    validated = load(infile)

# all 187 automated strain reconstructions (curated, no gap-fills)
with open("./model/pAo_strain-GEMs_all-187.pickle", 'rb') as infile:
    all187 = load(infile)
```

