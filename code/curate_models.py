"""Reusable, FBA-neutral model-quality curation for the pAo GEMs.

`curate(model)` applies, in place:
  1. the energy-generating-cycle (EGC) bound fix (where the reactions are present),
  2. chemical-formula corrections from data/genome/metabolite_formula_curation.csv (KEGG), and
  3. mass-balancing of reactions that are off by whole H2O or H+ molecules (missing water in
     hydrolyses, missing protons in redox steps), skipping reactions that touch a metabolite with
     no parseable formula (generic pseudo-metabolites that have no single formula).

None of this changes FBA predictions: formulae/charges are metadata, and the bound fix only removes
the three spurious cycles. Used by both p3_simulations.ipynb (pan model export) and
finalize_models.py (strain collections).
"""
import os
import re
from collections import Counter

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
FORMULA_CSV = os.path.join(_HERE, "..", "data", "genome", "metabolite_formula_curation.csv")

# Restrict each reaction to its thermodynamically feasible direction to remove the EGCs.
EGC_FIX = [("OtherAsp_R04962", "lower_bound", 0.0), ("OtherAsp_R01708", "upper_bound", 0.0),
           ("r1736", "lower_bound", 0.0), ("r1737", "upper_bound", 0.0)]

_GENERIC = re.compile(r"[RX*()]")          # generic/polymeric formula tokens (R-groups, (..)n)
_BIOMASS = {"r1897", "r2359", "r2358"}     # unbalanced by definition
_base = lambda mid: re.sub(r"\[[a-z]\]$", "", mid)   # strip compartment suffix -> KEGG base id


def _parseable(met):
    f = met.formula
    if not f or str(f).strip() in ("", "nan", "None") or _GENERIC.search(f):
        return False
    try:
        return bool(met.elements)
    except Exception:
        return False


def _imbalance(rxn):
    if not all(_parseable(x) for x in rxn.metabolites):
        return None
    return {k: v for k, v in rxn.check_mass_balance().items() if k != "charge" and abs(v) > 1e-6}


def curate(model, formula_csv=FORMULA_CSV):
    """Apply the EGC, formula and water/proton curations to `model` in place; return it."""
    # (1) EGC bound fix (only for reactions the model actually contains)
    for rid, attr, val in EGC_FIX:
        if model.reactions.has_id(rid):
            setattr(model.reactions.get_by_id(rid), attr, val)

    # (2) chemical-formula corrections (KEGG)
    fix = pd.read_csv(formula_csv).set_index("metabolite_base_id")["formula"].to_dict()
    for met in model.metabolites:
        if _base(met.id) in fix:
            met.formula = fix[_base(met.id)]

    # (3) balance whole-H2O / whole-H+ imbalances
    bnd = {r.id for r in model.boundary}
    for r in [x for x in model.reactions if x.id not in bnd and x.id not in _BIOMASS]:
        d = _imbalance(r)
        if not d:
            continue
        comp = Counter(x.compartment for x in r.metabolites).most_common(1)[0][0]
        if set(d) == {"H", "O"} and abs(d["H"]) == 2 * abs(d["O"]) and (d["H"] > 0) == (d["O"] > 0):
            wid = f"C00001[{comp}]" if model.metabolites.has_id(f"C00001[{comp}]") else "C00001[c]"
            if model.metabolites.has_id(wid):
                r.add_metabolites({model.metabolites.get_by_id(wid): -d["O"]})
        elif set(d) == {"H"}:
            pid = f"C00080[{comp}]" if model.metabolites.has_id(f"C00080[{comp}]") else "C00080[c]"
            if model.metabolites.has_id(pid):
                r.add_metabolites({model.metabolites.get_by_id(pid): -d["H"]})
    return model
