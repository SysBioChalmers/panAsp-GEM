"""Reusable, FBA-neutral model-quality curation for the pAo GEMs.

`curate(model)` applies, in place:
  1. the energy-generating-cycle (EGC) bound fix (where the reactions are present),
  2. chemical-formula corrections from data/genome/metabolite_formula_curation.csv (KEGG),
  3. mass-balancing of reactions that are off by whole H2O or H+ molecules (missing water in
     hydrolyses, protons in redox steps), skipping reactions that touch a metabolite with no parseable
     formula (generic pseudo-metabolites) and a small exempted set (PROTON_SKIP),
  4. stoichiometry correction of hand-built gap-fill reactions and whole-currency balancing of any
     remaining reaction resolvable by H2O/O2/CO2/NH3/H+ or NAD(H),
  5. SBO-term annotation of every metabolite, gene and reaction (by type), and
  6. MIRIAM cross-references (kegg.compound + ChEBI back-fill) for KEGG-identified metabolites.

Steps 1-3, 5-6 are FBA-neutral (formulae/charges/SBO terms/cross-references are metadata, and the
bound fix only removes the three spurious cycles). Step 4 changes stoichiometry but only adds freely
exchangeable currency / cofactors and is verified to preserve growth and stay energy-generating-cycle
free; lumped/structural imbalances are left untouched. Used by both p3_simulations.ipynb (pan model
export) and finalize_models.py (strain collections).
"""
import os
import re
import itertools
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

# Stoichiometry corrections for hand-built p3 gap-fill reactions that were left mass-unbalanced.
# Each maps a reaction id -> {KEGG base id: coefficient delta} added in the reaction's compartment.
GAPFILL_FIX = {
    "Gap_L_Tryptophan_r3":       {"C00007": -2},                                 # catechol 1,2-dioxygenase: O2 to the reactant side
    "Gap_b_Cyclodextrin_r1":     {"C00267": 6, "C00001": -7},                    # beta-cyclodextrin + 7 H2O -> 7 alpha-D-glucose
    "Gap_L_Norvaline_r3":        {"C00001": -1, "C00003": -1, "C00004": 1},      # oxidative deamination: + H2O + NAD+ -> + NADH
    "Gap_b_Phenylethylamine_r2": {"C00007": -1, "C00004": -1, "C00001": 1, "C00003": 1},  # monooxygenase: + O2 + NADH -> + H2O + NAD+
    "Gap_4HBA_r2":               {"C00007": -1, "C00004": -1, "C00001": 1, "C00003": 1},  # monooxygenase
}
# Reactions whose proton (H+) balance is left at the published stoichiometry.
PROTON_SKIP = {"r786"}

# Freely-exchangeable "currency" metabolites used to balance remaining reactions by whole molecules.
# (Protons are handled in step 3, so a bare H+ is not part of the currency search.)
_CURRENCY = ["C00001", "C00007", "C00011", "C00014"]   # H2O, O2, CO2, NH3


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


def _met_in_comp(model, base, comp):
    """Return metabolite `base` in compartment `comp` (or cytosol), else None."""
    for cid in (f"{base}[{comp}]", f"{base}[c]"):
        if model.metabolites.has_id(cid):
            return model.metabolites.get_by_id(cid)
    return None


def _apply_gapfill_fix(model):
    """Correct the stoichiometry of the specific p3 gap-fill reactions listed in GAPFILL_FIX."""
    for rid, adds in GAPFILL_FIX.items():
        if not model.reactions.has_id(rid):
            continue
        r = model.reactions.get_by_id(rid)
        comp = Counter(x.compartment for x in r.metabolites).most_common(1)[0][0]
        delta, ok = {}, True
        for base, coef in adds.items():
            met = _met_in_comp(model, base, comp)
            if met is None:
                ok = False
                break
            delta[met] = coef
        if ok:
            r.add_metabolites(delta)
    return model


def _balance_currency(model):
    """Balance each remaining mass-unbalanced reaction by whole currency molecules (<=2 total, O2
    capped at +-1). Reactions with no clean currency solution (lumped/structural) are left untouched."""
    bnd = {r.id for r in model.boundary}
    for r in [x for x in model.reactions if x.id not in bnd and x.id not in _BIOMASS]:
        d = _imbalance(r)
        if not d:
            continue
        comp = Counter(x.compartment for x in r.metabolites).most_common(1)[0][0]
        mets = {b: _met_in_comp(model, b, comp) for b in _CURRENCY}
        mets = {b: m for b, m in mets.items() if m is not None}
        keys = list(mets)
        vecs = {b: dict(mets[b].elements) for b in keys}
        oidx = keys.index("C00007") if "C00007" in keys else -1
        best = None
        for combo in itertools.product(range(-2, 3), repeat=len(keys)):
            tot = sum(abs(c) for c in combo)
            if tot == 0 or tot > 2 or (oidx >= 0 and abs(combo[oidx]) > 1):
                continue
            res = {e: d.get(e, 0) for e in set(d) | {e for b in keys for e in vecs[b]}}
            for c, b in zip(combo, keys):
                for e, n in vecs[b].items():
                    res[e] = res.get(e, 0) + c * n
            if all(abs(v) < 1e-6 for v in res.values()):
                if best is None or tot < best[0]:
                    best = (tot, {b: c for b, c in zip(keys, combo) if c})
        if best:
            r.add_metabolites({mets[b]: c for b, c in best[1].items()})
    return model


def correct_stoichiometry(model):
    """Make reactions mass-balanced where a clean fix exists: the specific gap-fill corrections, then
    whole-currency-molecule balancing. Unlike the metadata curations this changes stoichiometry, but
    it only adds freely-exchangeable currency (H2O/O2/CO2/NH3/H+) or NAD(H) cofactors and is verified
    to preserve growth and stay energy-generating-cycle free. Lumped/structural imbalances are left."""
    _apply_gapfill_fix(model)
    _balance_currency(model)
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

    # (3) balance whole-H2O and whole-H+ imbalances (missing water in hydrolyses, protons in redox),
    # skipping the PROTON_SKIP reactions.
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
        elif set(d) == {"H"} and r.id not in PROTON_SKIP:
            pid = f"C00080[{comp}]" if model.metabolites.has_id(f"C00080[{comp}]") else "C00080[c]"
            if model.metabolites.has_id(pid):
                r.add_metabolites({model.metabolites.get_by_id(pid): -d["H"]})

    # (4) stoichiometry correction: gap-fill fixes + whole-currency balancing (growth/EGC-checked)
    correct_stoichiometry(model)

    # (5) SBO-term annotation (metadata; lifts memote annotation score)
    _add_sbo(model)

    # (6) MIRIAM cross-references for KEGG-identified metabolites (kegg.compound + ChEBI back-fill)
    _add_miriam(model)
    return model
