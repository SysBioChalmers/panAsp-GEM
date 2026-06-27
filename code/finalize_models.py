"""Produce the four curated pAo model deliverables in model/ (cobra only -- no medusa/env):

  pAo.xml                              curated pan model (EGC-free, formulae fixed, mass-balanced)
  pAo_<strain>.xml  (x8)               the 8 validated strain models (with experimental gap-fills)
  pAo_strain-GEMs_validated.pickle     8 validated + template + Pan (with gap-fills)
  pAo_strain-GEMs_all-187.pickle       all 187 automated strain models (curated, no gap-fills)

Model-quality curations (EGC, formulae, water/proton mass-balance) are applied to every model via
curate_models.curate(). The 8-strain experimental gap-fills remain only on the validated set.
"""
import os
import pickle
import warnings

warnings.filterwarnings("ignore")
import cobra
from cobra import Reaction
from cobra.io import write_sbml_model

from curate_models import curate, _parseable, _imbalance

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL = os.path.join(ROOT, "model")
INTER = os.path.join(ROOT, "data", "intermediate")

STRAIN_SHORT = {
    "Aspergillus_oryzae_RIB40_GCF_000184455.2": "RIB40",
    "Aspergillus_oryzae_NRRL_2217": "NRRL_2217",
    "Aspergillus_oryzae_NRRL_3483_GCA_034767915.1": "NRRL_3483",
    "Aspergillus_oryzae_NRRL_3488": "NRRL_3488",
    "Aspergillus_oryzae_NRRL_5589": "NRRL_5589",
    "Aspergillus_oryzae_NRRL_5592_GCA_034767935.1": "NRRL_5592",
    "Aspergillus_oryzae_NRRL_35890_GCA_034767955.1": "NRRL_35890",
    "Aspergillus_oryzae_NRRL_471": "NRRL_471",
}
DISS = {  # 12-carrier energy-generating-cycle screen
    "ATP": {"C00002[c]": -1, "C00001[c]": -1, "C00008[c]": 1, "C00009[c]": 1, "C00080[c]": 1},
    "CTP": {"C00063[c]": -1, "C00001[c]": -1, "C00112[c]": 1, "C00009[c]": 1, "C00080[c]": 1},
    "GTP": {"C00044[c]": -1, "C00001[c]": -1, "C00035[c]": 1, "C00009[c]": 1, "C00080[c]": 1},
    "UTP": {"C00075[c]": -1, "C00001[c]": -1, "C00015[c]": 1, "C00009[c]": 1, "C00080[c]": 1},
    "ITP": {"C00081[c]": -1, "C00001[c]": -1, "C00104[c]": 1, "C00009[c]": 1, "C00080[c]": 1},
    "NADH": {"C00004[c]": -1, "C00003[c]": 1, "C00080[c]": 1},
    "NADPH": {"C00005[c]": -1, "C00006[c]": 1, "C00080[c]": 1},
    "FADH2": {"C01352[m]": -1, "C00016[m]": 1, "C00080[m]": 2},
    "FMNH2": {"C01847[c]": -1, "C00061[c]": 1, "C00080[c]": 2},
    "AcetylCoA": {"C00024[c]": -1, "C00001[c]": -1, "C00010[c]": 1, "C00033[c]": 1},
    "L-Glutamate": {"C00025[c]": -1, "C00026[c]": 1, "C00014[c]": 1},
    "Ubiquinol8": {"C00390[m]": -1, "C00399[m]": 1, "C00080[m]": 2},
}


def load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def max_egc(model):
    ids = {x.id for x in model.metabolites}
    mx = 0.0
    with model:
        for rid in list(model.medium):           # close the medium (atomic bounds set avoids lb>ub)
            model.reactions.get_by_id(rid).bounds = (0.0, 0.0)
        for n, st in DISS.items():
            if any(i not in ids for i in st):
                continue
            with model:
                d = Reaction("DISS_" + n, lower_bound=0, upper_bound=1000)
                model.add_reactions([d])
                d.add_metabolites({model.metabolites.get_by_id(i): c for i, c in st.items()})
                model.objective = d
                s = model.optimize()
                mx = max(mx, s.objective_value if s.status == "optimal" else 0.0)
    return mx


def growth(model):
    with model:
        if model.reactions.has_id("r2359"):
            model.objective = "r2359"
        if model.metabolites.has_id("C00031[e]"):
            try:
                model.add_boundary(model.metabolites.get_by_id("C00031[e]"), type="exchange", lb=-10, ub=1000)
            except ValueError:
                pass
        return model.slim_optimize()


def mass_balanced_fraction(model):
    bnd = {r.id for r in model.boundary}
    rxns = [r for r in model.reactions if r.id not in bnd and r.id not in {"r1897", "r2359", "r2358"}]
    unbal = 0
    for r in rxns:
        if any(not _parseable(x) for x in r.metabolites):
            unbal += 1
        elif _imbalance(r):
            unbal += 1
    return 100 * (1 - unbal / len(rxns)), len(rxns)


# ── 1. all 187 automated strain models (no gap-fills) ──────────────────────────
print("Loading the 187-strain source pickle (large, ~2-3 min)...", flush=True)
all187 = load(os.path.join(INTER, "strain-GEMs_automated_187.pickle"))
strains187 = [m for m in all187 if m.id not in ("template", "Pan_oryzae") and not m.id.lower().startswith("pan")]
print(f"  curating {len(strains187)} strain models...", flush=True)
for i, m in enumerate(strains187):
    curate(m)
    if (i + 1) % 40 == 0:
        print(f"    {i + 1}/{len(strains187)}", flush=True)
with open(os.path.join(MODEL, "pAo_strain-GEMs_all-187.pickle"), "wb") as f:
    pickle.dump(strains187, f)
print(f"WROTE model/pAo_strain-GEMs_all-187.pickle ({len(strains187)} strains)", flush=True)

# ── 2. validated 10 (with gap-fills): pickle + 8 strain XML + pan XML ───────────
val = load(os.path.join(INTER, "strain-GEMs_validated_gapfilled.pickle"))   # 10 gap-filled members (p3 output)
for m in val:
    curate(m)
with open(os.path.join(MODEL, "pAo_strain-GEMs_validated.pickle"), "wb") as f:
    pickle.dump(val, f)
print(f"WROTE model/pAo_strain-GEMs_validated.pickle ({len(val)} members)", flush=True)

by_id = {m.id: m for m in val}          # (the validated pickle above keeps the original full ids)
for fid, short in STRAIN_SHORT.items():
    if fid in by_id:
        m = by_id[fid]
        m.id = f"pAo_{short}"            # valid SBML SId — the full strain id contains a '.'
        write_sbml_model(m, os.path.join(MODEL, f"pAo_{short}.xml"))
pan = by_id.get("Pan_oryzae")
pan.id = "pAo"
write_sbml_model(pan, os.path.join(MODEL, "pAo.xml"))
print(f"WROTE pAo.xml + {sum(1 for f in STRAIN_SHORT if f in by_id)} strain XMLs", flush=True)

# ── 3. verification ────────────────────────────────────────────────────────────
print("\n=== verification ===")
frac, nr = mass_balanced_fraction(pan)
print(f"pAo (pan): growth={growth(pan):.3f}  EGC={max_egc(pan):.2e}  mass-balanced={frac:.1f}% of {nr}")
for fid in ["Aspergillus_oryzae_RIB40_GCF_000184455.2", "Aspergillus_oryzae_NRRL_35890_GCA_034767955.1"]:
    if fid in by_id:
        m = by_id[fid]
        print(f"{STRAIN_SHORT[fid]:11s} (validated): growth={growth(m):.3f}  EGC={max_egc(m):.2e}")
s0 = strains187[0]
print(f"{s0.id[:24]} (automated): growth={growth(s0):.3f}  EGC={max_egc(s0):.2e}")
print("done")
