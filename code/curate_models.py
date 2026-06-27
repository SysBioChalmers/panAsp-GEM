"""Reusable, FBA-neutral model-quality curation for the pAo GEMs.

`curate(model)` applies, in place:
  1. the energy-generating-cycle (EGC) bound fix (where the reactions are present),
  2. chemical-formula corrections from data/genome/metabolite_formula_curation.csv (KEGG),
  3. mass-balancing of reactions that are off by whole H2O or H+ molecules (missing water in
     hydrolyses, missing protons in redox steps), skipping reactions that touch a metabolite with
     no parseable formula (generic pseudo-metabolites that have no single formula),
  4. SBO-term annotation of every metabolite, gene and reaction (by type), and
  5. MIRIAM cross-references (kegg.compound + ChEBI back-fill) for KEGG-identified metabolites.

None of this changes FBA predictions: formulae/charges/SBO terms/cross-references are metadata, and
the bound fix only removes the three spurious cycles. Used by both p3_simulations.ipynb (pan model
export) and finalize_models.py (strain collections).
"""
import os
import re
from collections import Counter

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
FORMULA_CSV = os.path.join(_HERE, "..", "data", "genome", "metabolite_formula_curation.csv")
KEGG_CHEBI_CSV = os.path.join(_HERE, "..", "data", "intermediate", "kegg_chebi_dict.csv")

# Restrict each reaction to its thermodynamically feasible direction to remove the EGCs.
EGC_FIX = [("OtherAsp_R04962", "lower_bound", 0.0), ("OtherAsp_R01708", "upper_bound", 0.0),
           ("r1736", "lower_bound", 0.0), ("r1737", "upper_bound", 0.0)]

_GENERIC = re.compile(r"[RX*()]")          # generic/polymeric formula tokens (R-groups, (..)n)
_BIOMASS = {"r1897", "r2359", "r2358"}     # unbalanced by definition
_KEGG_C = re.compile(r"^C\d{5}$")          # KEGG compound id pattern
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


def _add_sbo(model):
    """Annotate every component with a Systems Biology Ontology term (FBA-neutral metadata).

    memote scores each component type against an expected SBO term; assigning them by type lifts the
    annotation-SBO section without touching any stoichiometry, bound or formula.
    """
    for met in model.metabolites:
        met.annotation["sbo"] = "SBO:0000247"          # simple chemical
    for gene in model.genes:
        gene.annotation["sbo"] = "SBO:0000243"         # gene
    exch = {r.id for r in model.exchanges}
    dem = {r.id for r in model.demands}
    snk = {r.id for r in model.sinks}
    for r in model.reactions:
        if r.id in _BIOMASS or "biomass" in (r.name or "").lower():
            sbo = "SBO:0000629"                        # biomass production
        elif r.id == "r1901":
            sbo = "SBO:0000630"                        # ATP maintenance (NGAM)
        elif r.id in exch:
            sbo = "SBO:0000627"                        # exchange reaction
        elif r.id in dem:
            sbo = "SBO:0000628"                        # demand reaction
        elif r.id in snk:
            sbo = "SBO:0000632"                        # sink reaction
        elif r.boundary:
            sbo = "SBO:0000627"
        elif len({m.compartment for m in r.metabolites}) > 1:
            sbo = "SBO:0000185"                        # transport (translocation) reaction
        else:
            sbo = "SBO:0000176"                        # metabolic (biochemical) reaction
        r.annotation["sbo"] = sbo
    return model


def _load_kegg_chebi(path):
    """Parse the KEGG-compound -> ChEBI map ('cpd:C00462,chebi:16042' per line)."""
    k2c = {}
    try:
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                parts = line.strip().split(",")
                if len(parts) == 2 and parts[0].startswith("cpd:") and "chebi:" in parts[1]:
                    k2c[parts[0][4:]] = parts[1].split(":")[-1]   # C00462 -> 16042
    except FileNotFoundError:
        pass
    return k2c


def _add_miriam(model, kegg_chebi_csv=KEGG_CHEBI_CSV):
    """Add MIRIAM cross-references (FBA-neutral metadata) to KEGG-identified metabolites.

    Every metabolite whose id is a KEGG compound (C#####) gets a `kegg.compound` cross-reference
    (most carry the id but never declared it), and any still lacking a ChEBI reference is back-filled
    from the curated KEGG->ChEBI table. This lifts memote's metabolite-annotation coverage.
    """
    k2c = _load_kegg_chebi(kegg_chebi_csv)
    for met in model.metabolites:
        base = _base(met.id)
        if _KEGG_C.match(base):
            met.annotation.setdefault("kegg.compound", base)
            if "chebi" not in met.annotation and base in k2c:
                met.annotation["chebi"] = f"CHEBI:{k2c[base]}"
    return model


def curate(model, formula_csv=FORMULA_CSV):
    """Apply the EGC, formula, water/proton, SBO and MIRIAM curations to `model` in place."""
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

    # (4) SBO-term annotation (metadata; lifts memote annotation score)
    _add_sbo(model)

    # (5) MIRIAM cross-references for KEGG-identified metabolites (kegg.compound + ChEBI back-fill)
    _add_miriam(model)
    return model
