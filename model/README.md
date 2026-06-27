# pAo model files

The consensus genome-scale metabolic model of *Aspergillus oryzae* (**pAo**) and its
strain-specific derivatives. Every model here has been put through the shared, FBA-neutral
model-quality curation in [`code/curate_models.py`](../code/curate_models.py) — the energy-generating-cycle
bound fix, KEGG chemical-formula corrections, and whole-H₂O/H⁺ mass-balancing — and the set is
assembled by [`code/finalize_models.py`](../code/finalize_models.py).

## Files

| File | Contents |
|---|---|
| `pAo.xml` | The consensus **pan-model** as SBML. Curated; energy-generating-cycle free. |
| `pAo_<strain>.xml` | The **8 experimentally validated** strain models as SBML (RIB40, NRRL_2217, NRRL_3483, NRRL_3488, NRRL_5589, NRRL_5592, NRRL_35890, NRRL_471), each carrying its experimental gap-fills. |
| `pAo_strain-GEMs_validated.pickle` | A pickled `list` of cobra models: the 8 validated strains **+ `template` + `Pan_oryzae`** (10 members), with gap-fills. This is what the strain-comparison analyses load. |
| `pAo_strain-GEMs_all-187.pickle` | A pickled `list` of all **187** automated strain reconstructions, curated but **without** experimental gap-fills. Large (~167 MB) — distributed via Zenodo, not committed to git. |

Load a pickle with:

```python
from pickle import load
with open("pAo_strain-GEMs_validated.pickle", "rb") as f:
    models = load(f)        # list of cobra.Model
```

## The 8 validated strains appear twice — by design

The eight validated strains are present in **both** collections, with **different** content:

- in `pAo_strain-GEMs_validated.pickle` / `pAo_<strain>.xml` they carry the **experimental
  gap-fills** added during BioLog-guided curation, so they reproduce the observed growth phenotypes;
- in `pAo_strain-GEMs_all-187.pickle` they are the **automated, gene-content-only**
  reconstructions (no gap-fills), so that all 187 strains are treated identically and the
  pan-genome reaction statistics are not biased by the eight hand-curated members.

## Regeneration

`finalize_models.py` rebuilds every file above from two intermediates (kept out of git, on Zenodo):

- `data/intermediate/strain-GEMs_automated_187.pickle` — the 187 automated reconstructions (from `p2`);
- `data/intermediate/strain-GEMs_validated_gapfilled.pickle` — the 10 gap-filled members (from `p3`).

It needs only cobra + gurobi (no medusa / jupyter):

```bash
python code/finalize_models.py
```

The curation is FBA-neutral: formulae and charges are metadata, and the bound fix only removes the
three spurious energy-generating cycles, so `slim_optimize()` is unchanged for every model.
