"""Derive the 187-strain pan-genome reaction classification used to label reactions as
core / accessory in the strain-comparison analysis (consistent with the manuscript, Figure 2).

A reaction is classified by the fraction of the 187 strain-specific GEMs that contain it:
  strict_core (100%) | softcore (>95%) | shell (5-95%) | cloud (<5%, >1 strain) | strain_specific (1)
  class = core      if present in >95% of strains, else accessory.

Run once against the 187-strain gemList pickle to (re)generate
`data/genome/reaction_classification_187.csv`. The pickle is large and not committed; the derived
CSV is the committed, static artifact the pipeline reads. The 8 curated strain models alone are NOT
sufficient — they collapse the >95% threshold — so the full 187-strain collection is required here.

Usage:
    python code/derive_reaction_classification.py [path/to/strain-GEMs_automated_187.pickle]
"""
import os, sys
from pickle import load
import numpy as np, pandas as pd

DEFAULT_PICKLE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              'data', 'intermediate', 'strain-GEMs_automated_187.pickle')
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'data', 'genome', 'reaction_classification_187.csv')


def main(pickle_path):
    with open(pickle_path, 'rb') as f:
        all_models = load(f)
    # keep only the strain-specific models (drop the template and the pan-model)
    strains = [m for m in all_models
               if m.id not in {'template', 'Pan_oryzae', 'panAsp', 'pAo'} and not m.id.lower().startswith('pan')]
    n = len(strains)
    if n < 100:
        sys.exit(f"Expected ~187 strain models, found {n} in {pickle_path}. "
                 f"Use the full 187-strain collection, not the 8-strain/ensemble subset.")
    all_rxns = sorted({r.id for m in strains for r in m.reactions})
    pres = pd.DataFrame({m.id: {r.id: 1 for r in m.reactions} for m in strains},
                        index=all_rxns).fillna(0).astype(np.int8)
    cnt = pres.sum(axis=1); frac = cnt / n
    out = pd.DataFrame({'reaction': all_rxns, 'n_strains': cnt.values, 'fraction': frac.values})
    out['category'] = np.select(
        [out.fraction == 1, out.fraction > 0.95, out.fraction > 0.05, out.n_strains > 1],
        ['strict_core', 'softcore', 'shell', 'cloud'], default='strain_specific')
    out['class'] = np.where(out.fraction > 0.95, 'core', 'accessory')
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    out.to_csv(OUT, index=False)
    cats = out['category'].value_counts().to_dict()
    print(f"{n} strains; {len(out)} reactions -> {OUT}")
    print(f"  core (>95%): {(out['class'] == 'core').sum()};  accessory (<=95%): {(out['class'] == 'accessory').sum()}")
    print(f"  categories: {cats}")


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_PICKLE
    if not os.path.exists(path):
        sys.exit(f"Pickle not found: {path}\nPass the path to strain-GEMs_automated_187.pickle as an argument.")
    main(path)
