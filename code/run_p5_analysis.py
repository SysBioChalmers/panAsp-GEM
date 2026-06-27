"""
Analysis script corresponding to p5_strain_comparison.ipynb.
Run from the code/ directory or project root.

Addresses reviewer points Q5 (EGCs), Q12 (flux variation), Q13 (accessory
reactions, NGAM), Q14 (metabolic costs, subsystem profiles) on the eight
curated A. oryzae strain models plus the Pan_oryzae pan-model.

Methodological choices (documented for the rebuttal):
  * Carbon sources are fed at EQUAL CARBON (C-mol) supply, not equal mmol, so
    growth/cost differences reflect metabolism rather than carbon content.
  * A non-zero non-growth ATP maintenance (NGAM) is imposed (the shipped model
    has lb=0). Growth-associated maintenance (GAM) is reported.
  * Flux variability analysis uses loopless=True so reported ranges are
    thermodynamically attainable (free of internal loop-law violations).
"""
import os, sys, warnings, pickle
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from cobra import Reaction, Metabolite
from cobra.flux_analysis import pfba, flux_variability_analysis, single_reaction_deletion
from pickle import load

warnings.filterwarnings('ignore')
sns.set_theme(style='whitegrid', context='notebook')

# ── paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT    = os.path.dirname(SCRIPT_DIR)
OUT_DIR    = os.path.join(PROJECT, 'data', 'intermediate')
os.makedirs(OUT_DIR, exist_ok=True)

# ── constants ─────────────────────────────────────────────────────────────────
TARGET_STRAINS = [
    'Aspergillus_oryzae_RIB40_GCF_000184455.2',
    'Aspergillus_oryzae_NRRL_2217',
    'Aspergillus_oryzae_NRRL_3483_GCA_034767915.1',
    'Aspergillus_oryzae_NRRL_3488',
    'Aspergillus_oryzae_NRRL_5589',
    'Aspergillus_oryzae_NRRL_5592_GCA_034767935.1',
    'Aspergillus_oryzae_NRRL_35890_GCA_034767955.1',
    'Aspergillus_oryzae_NRRL_471',
]
STRAIN_SHORT = {
    'Aspergillus_oryzae_RIB40_GCF_000184455.2':     'RIB40',
    'Aspergillus_oryzae_NRRL_2217':                  'NRRL_2217',
    'Aspergillus_oryzae_NRRL_3483_GCA_034767915.1':  'NRRL_3483',
    'Aspergillus_oryzae_NRRL_3488':                  'NRRL_3488',
    'Aspergillus_oryzae_NRRL_5589':                  'NRRL_5589',
    'Aspergillus_oryzae_NRRL_5592_GCA_034767935.1':  'NRRL_5592',
    'Aspergillus_oryzae_NRRL_35890_GCA_034767955.1': 'NRRL_35890',
    'Aspergillus_oryzae_NRRL_471':                   'NRRL_471',
}
CARBON_SOURCES = {
    'Glucose':  'C00031[e]',
    'Glycerol': 'C00116[e]',
    'Maltose':  'C00208[e]',
    'Xylose':   'C00181[e]',
}
# Carbon atoms per molecule — used to normalise uptake to equal C-mol supply.
CARBON_ATOMS = {'Glucose': 6, 'Glycerol': 3, 'Maltose': 12, 'Xylose': 5}
CS_ORDER = ['Glucose', 'Glycerol', 'Maltose', 'Xylose']
STRAINS_ORDERED = [STRAIN_SHORT[s] for s in TARGET_STRAINS]

BIOMASS_RXN = 'r2359'        # growth objective (drains Biomass[c])
BIOMASS_ASM = 'r1897'        # biomass-formation reaction (holds GAM ATP term)
NGAM_RXN    = 'r1901'        # non-growth ATP maintenance
O2_EXCHANGE = 'r2202'        # O2 exchange (uptake reported as negative flux)
ATP_C       = 'C00002[c]'    # cytosolic ATP

# Carbon uptake normalised to equal C-mol across substrates, anchored to the
# documented specific glucose uptake rate qS = 1.12 mmol glucose/gDW/h used for the
# biomass measurement (manuscript Suppl. Table S5; A1560, glucose/ammonia, mu=0.1/h).
# Glucose -> 1.12 mmol/gDW/h (6 C); the other substrates are fed at the same carbon
# flux: glycerol 2.24, maltose 0.56, xylose 1.344 mmol/gDW/h.
GLUCOSE_QS       = 1.12                                  # mmol glucose/gDW/h (paper)
CARBON_CMOL_RATE = GLUCOSE_QS * CARBON_ATOMS['Glucose']  # = 6.72 mmol C/gDW/h
# Non-growth ATP maintenance. The shipped model has r1901 lb=0 (no maintenance);
# we impose a literature-consistent value for filamentous fungi (~1 mmol
# ATP/gDW/h; cf. A. niger 1.5-3.6). Section 6 sweeps this to show robustness.
NGAM_VALUE = 1.0
# Growth-associated maintenance, read from r1897 (reported, not changed).
GAM_ATP    = 49.0
# Numerical floor: differences below this (relative) are at the LP solver
# tolerance and are reported as "not resolved", not as biology.
NOISE_REL  = 1e-4

# Pan-model EGC curation (already applied in the pickle via p3_simulations.ipynb
# cell 112). Recorded here so Section 1 can demonstrate causality by reverting.
EGC_FIX = {
    # reaction_id: (attribute, fixed_value, reverted_value)
    'OtherAsp_R04962': ('lower_bound', 0.0,  -1000.0),
    'OtherAsp_R01708': ('upper_bound', 0.0,   1000.0),
    'r1736':           ('lower_bound', 0.0,  -1000.0),
    'r1737':           ('upper_bound', 0.0,   1000.0),
}


# ── load models ────────────────────────────────────────────────────────────────
print("=== Loading models ===")
pickle_path = os.path.join(PROJECT, 'model', 'pAo_strain-GEMs_validated.pickle')
with open(pickle_path, 'rb') as f:
    all_models = load(f)
models = {m.id: m for m in all_models if m.id in TARGET_STRAINS}
print(f"Loaded {len(models)} / {len(TARGET_STRAINS)} target strain models")
for sid in TARGET_STRAINS:
    status = 'OK' if sid in models else 'MISSING'
    print(f"  [{status}] {STRAIN_SHORT[sid]}")

# Pan_oryzae pan-model (EGC testing only; not used for growth predictions).
pan_model = next((m for m in all_models if m.id == 'Pan_oryzae'), None)
if pan_model is not None:
    pan_model.objective = pan_model.reactions.get_by_id(BIOMASS_RXN)
    print(f"  [OK] Pan_oryzae  (rxns={len(pan_model.reactions)}, mets={len(pan_model.metabolites)})")
else:
    print("  [MISSING] Pan_oryzae")

# Set biomass (r2359) as default objective (pickle ships with r1901 as default).
for model in models.values():
    model.objective = model.reactions.get_by_id(BIOMASS_RXN)

# Gapfill r766 (5-proFAR isomerase, histidine biosynthesis) for NRRL_35890.
_NRRL35890 = 'Aspergillus_oryzae_NRRL_35890_GCA_034767955.1'
if _NRRL35890 in models:
    m35 = models[_NRRL35890]
    if 'r766' not in {r.id for r in m35.reactions}:
        _ref = next(m for sid, m in models.items()
                    if sid != _NRRL35890 and 'r766' in {r.id for r in m.reactions})
        _ref_rxn = _ref.reactions.get_by_id('r766')
        _new_rxn = Reaction('r766', name=_ref_rxn.name,
                            lower_bound=_ref_rxn.lower_bound,
                            upper_bound=_ref_rxn.upper_bound)
        _new_rxn.subsystem = _ref_rxn.subsystem
        m35.add_reactions([_new_rxn])
        _m35_met_ids = {m.id for m in m35.metabolites}
        _met_dict = {}
        for _met, _coef in _ref_rxn.metabolites.items():
            if _met.id in _m35_met_ids:
                _met_dict[m35.metabolites.get_by_id(_met.id)] = _coef
            else:
                _new_met = Metabolite(id=_met.id, name=_met.name,
                                     formula=_met.formula, compartment=_met.compartment)
                m35.add_metabolites([_new_met])
                _met_dict[_new_met] = _coef
        _new_rxn.add_metabolites(_met_dict)
        print(f"  Applied r766 gapfill for NRRL_35890 (histidine biosynthesis)")

# Report the energy parameters (GAM/NGAM) the model actually uses.
_ref = models[TARGET_STRAINS[0]]
_asm = _ref.reactions.get_by_id(BIOMASS_ASM)
_gam = _asm.metabolites.get(_ref.metabolites.get_by_id(ATP_C), None)
print(f"\nEnergy parameters:")
print(f"  GAM (ATP in {BIOMASS_ASM} biomass formation): {_gam} mmol ATP/gDW")
print(f"  NGAM ({NGAM_RXN}) shipped lb = "
      f"{_ref.reactions.get_by_id(NGAM_RXN).lower_bound}; imposed = {NGAM_VALUE} mmol ATP/gDW/h")


# ── helper ─────────────────────────────────────────────────────────────────────
def setup_condition(model, cs_name, carbon_met_id, ngam=NGAM_VALUE):
    """Configure a growth condition inside a `with model:` block:
       (1) carbon uptake normalised to equal C-mol (CARBON_CMOL_RATE),
       (2) non-growth ATP maintenance (NGAM) lower bound.
    Returns the substrate uptake rate used (mmol substrate/gDW/h)."""
    uptake = CARBON_CMOL_RATE / CARBON_ATOMS[cs_name]
    if carbon_met_id in {m.id for m in model.metabolites}:
        model.add_boundary(model.metabolites.get_by_id(carbon_met_id),
                           type='exchange', lb=-uptake, ub=1000)
    else:
        print(f"  Warning: {carbon_met_id} not in {model.id}")
    if ngam is not None and NGAM_RXN in {r.id for r in model.reactions}:
        model.reactions.get_by_id(NGAM_RXN).lower_bound = ngam
    return uptake


def atp_turnover(model, fluxes):
    """Total ATP regenerated per hour in a flux distribution (mmol ATP/gDW/h):
    sum of ATP *production* across all reactions touching cytosolic ATP.
    A defensible energetic-cost proxy (equals ATP consumption at steady state)."""
    if ATP_C not in {m.id for m in model.metabolites}:
        return float('nan')
    atp = model.metabolites.get_by_id(ATP_C)
    total = 0.0
    for rxn in atp.reactions:
        v = fluxes.get(rxn.id, 0.0)
        produced = rxn.metabolites[atp] * v
        if produced > 0:
            total += produced
    return total


# ═══════════════════════════════════════════════════════════════════════════════
# 1. EGC detection (Q5) — with before/after-curation demonstration on Pan_oryzae
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== 1. EGC detection ===")
DISSIPATION_SPECS = {
    'ATP':         {'C00002[c]': -1, 'C00001[c]': -1, 'C00008[c]': 1, 'C00009[c]': 1, 'C00080[c]': 1},
    'CTP':         {'C00063[c]': -1, 'C00001[c]': -1, 'C00112[c]': 1, 'C00009[c]': 1, 'C00080[c]': 1},
    'GTP':         {'C00044[c]': -1, 'C00001[c]': -1, 'C00035[c]': 1, 'C00009[c]': 1, 'C00080[c]': 1},
    'UTP':         {'C00075[c]': -1, 'C00001[c]': -1, 'C00015[c]': 1, 'C00009[c]': 1, 'C00080[c]': 1},
    'ITP':         {'C00081[c]': -1, 'C00001[c]': -1, 'C00104[c]': 1, 'C00009[c]': 1, 'C00080[c]': 1},
    'NADH':        {'C00004[c]': -1, 'C00003[c]': 1,  'C00080[c]': 1},
    'NADPH':       {'C00005[c]': -1, 'C00006[c]': 1,  'C00080[c]': 1},
    'FADH2':       {'C01352[m]': -1, 'C00016[m]': 1,  'C00080[m]': 2},
    'FMNH2':       {'C01847[c]': -1, 'C00061[c]': 1,  'C00080[c]': 2},
    'AcetylCoA':   {'C00024[c]': -1, 'C00001[c]': -1, 'C00010[c]': 1, 'C00033[c]': 1},
    'L-Glutamate': {'C00025[c]': -1, 'C00026[c]': 1,  'C00014[c]': 1},
    'Ubiquinol8':  {'C00390[m]': -1, 'C00399[m]': 1,  'C00080[m]': 2},
}

def test_egc(model, specs=DISSIPATION_SPECS, tol=1e-6, verbose=False):
    met_ids = {m.id for m in model.metabolites}
    results = {}
    with model:
        model.medium = {r_id: 0.0 for r_id in model.medium}
        for name, stoich in specs.items():
            if any(mid not in met_ids for mid in stoich):
                results[name] = float('nan')
                continue
            with model:
                diss = Reaction(f'DISS_{name}', lower_bound=0.0, upper_bound=1000.0)
                model.add_reactions([diss])
                diss.add_metabolites(
                    {model.metabolites.get_by_id(mid): c for mid, c in stoich.items()})
                model.objective = diss
                sol = model.optimize()
                flux = sol.objective_value if sol.status == 'optimal' else 0.0
                results[name] = flux
                if verbose and flux > tol:
                    print(f"    EGC [{name}] flux={flux:.1f}")
    return results

egc_rows = []
for strain_id, model in models.items():
    print(f"  {STRAIN_SHORT[strain_id]}", end=' ', flush=True)
    row = test_egc(model)
    row['model'] = STRAIN_SHORT[strain_id]
    egc_rows.append(row)
    print('done')

# Pan_oryzae: after curation (current pickle) should be all-zero.
egc_before = {}
if pan_model is not None:
    print("  Pan_oryzae (after fix)", end=' ', flush=True)
    row = test_egc(pan_model)
    row['model'] = 'Pan_oryzae'
    egc_rows.append(row)
    print('done')

    # Demonstrate causality: revert the four curated bounds -> EGCs return.
    print("  Pan_oryzae (before fix, bounds reverted)", end=' ', flush=True)
    with pan_model:
        for rid, (attr, _fixed, reverted) in EGC_FIX.items():
            setattr(pan_model.reactions.get_by_id(rid), attr, reverted)
        egc_before = test_egc(pan_model)
    print('done')

egc_df = pd.DataFrame(egc_rows).set_index('model')
numeric_cols = egc_df.select_dtypes(include='number')
max_flux = numeric_cols.max().max()
print(f"\nMax dissipation flux (all {len(egc_df)} models x 12 carriers): {max_flux:.2e}")
print("EGC result:", "NO EGCs detected" if max_flux < 1e-6 else "EGCs DETECTED")
if egc_before:
    nz = {k: round(v, 0) for k, v in egc_before.items() if v == v and abs(v) > 1e-6}
    print(f"Pan_oryzae BEFORE the four-reaction curation: EGCs on {nz}")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. FBA growth predictions (Q12) — carbon-normalised, with biomass yield
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== 2. FBA growth predictions (equal C-mol uptake, NGAM imposed) ===")
fba_rows = []
for strain_id, model in models.items():
    for cs_name, cs_met in CARBON_SOURCES.items():
        with model:
            uptake = setup_condition(model, cs_name, cs_met)
            model.objective = model.reactions.get_by_id(BIOMASS_RXN)
            sol = model.optimize()
            growth = sol.objective_value if sol.status == 'optimal' else 0.0
            o2 = sol.fluxes.get(O2_EXCHANGE, float('nan')) if sol.status == 'optimal' else float('nan')
        # Biomass yield per C-mol of substrate (uptake-rate independent).
        yield_cmol = growth / CARBON_CMOL_RATE if growth else 0.0
        dbl_h = (np.log(2) / growth) if growth > 1e-9 else float('nan')
        fba_rows.append({
            'strain': strain_id, 'strain_short': STRAIN_SHORT[strain_id],
            'carbon_source': cs_name, 'uptake_mmol': uptake,
            'growth_rate': growth, 'yield_per_cmol': yield_cmol,
            'doubling_h': dbl_h, 'o2_flux': o2, 'status': sol.status,
        })

growth_df = pd.DataFrame(fba_rows)
pivot_growth = growth_df.pivot(index='strain_short', columns='carbon_source', values='growth_rate')
pivot_yield  = growth_df.pivot(index='strain_short', columns='carbon_source', values='yield_per_cmol')
print("FBA growth rates (h-1), equal C-mol supply:")
print(pivot_growth[CS_ORDER].to_string())

# Noise-floor report: is the inter-strain spread above solver tolerance?
print("\nInter-strain growth spread (is it resolved above solver tolerance?):")
for cs in CS_ORDER:
    col = pivot_growth[cs]
    rng = col.max() - col.min()
    rel = rng / col.mean() if col.mean() else 0
    flag = 'RESOLVED' if rel > NOISE_REL else f'at noise floor (<{NOISE_REL:.0e})'
    print(f"  {cs:9s} mean={col.mean():.4f}  spread={rng:.2e}  rel={rel:.1e}  [{flag}]")
print(f"\nDoubling times (h): glucose ~{growth_df.query('carbon_source==\"Glucose\"')['doubling_h'].mean():.2f}  "
      f"(absolute rates depend on the {CARBON_CMOL_RATE:.0f} mmol-C/gDW/h uptake assumption; "
      f"yields are uptake-independent)")

# Compare FBA with BioLog experimental data (phenotype concordance)
biolog_ordered = None
biolog_path = os.path.join(PROJECT, 'data', 'biolog', 'PM_all.csv')
try:
    pm = pd.read_csv(biolog_path, sep=';')
    cs_map = {'D_Glucose': 'Glucose', 'Glycerol': 'Glycerol', 'D_Maltose': 'Maltose', 'D_Xylose': 'Xylose'}
    pm_sub = pm[pm['Medium'].isin(cs_map) & pm['Isolate'].isin(STRAIN_SHORT.values())].copy()
    pm_sub['carbon_source'] = pm_sub['Medium'].map(cs_map)

    def _frac_true(series):
        n_true  = (series == True).sum()   # noqa: E712
        n_false = (series == False).sum()  # noqa: E712
        total = n_true + n_false
        return n_true / total if total > 0 else float('nan')

    def _n_rep(series):
        return int((series == True).sum() + (series == False).sum())  # noqa: E712

    biolog_consensus = (
        pm_sub.groupby(['Isolate', 'carbon_source'])['Observed']
        .apply(_frac_true).reset_index(name='frac_true'))
    biolog_n = (pm_sub.groupby(['Isolate', 'carbon_source'])['Observed']
                .apply(_n_rep).reset_index(name='n_rep'))
    biolog_pivot = biolog_consensus.pivot(index='Isolate', columns='carbon_source', values='frac_true')
    biolog_ordered = biolog_pivot.reindex(
        index=sorted(STRAIN_SHORT.values()), columns=CS_ORDER)
    print("\nBioLog consensus growth (fraction of replicates positive; "
          f"median replicates/condition = {int(biolog_n['n_rep'].median())}):")
    print(biolog_ordered.to_string())

    # Concordance — compare ONLY where BioLog data exists (many conditions untested).
    bio = biolog_ordered.reindex(index=STRAINS_ORDERED, columns=CS_ORDER)
    model_grows = (pivot_growth.reindex(index=STRAINS_ORDERED, columns=CS_ORDER)
                   > 0.01 * pivot_growth.values.max())
    tested       = bio.notna()
    biolog_grows = (bio > 0.5)
    n_tested = int(tested.values.sum())
    n_agree  = int(((model_grows == biolog_grows) & tested).values.sum())
    fp_mask  = (model_grows & (~biolog_grows) & tested)   # model grows, wet lab does not
    n_fp     = int(fp_mask.values.sum())
    fp_list  = [f"{STRAINS_ORDERED[i]}/{CS_ORDER[j]}"
                for i, j in zip(*np.where(fp_mask.values))]
    print(f"\nModel/BioLog concordance over TESTED conditions only "
          f"({n_tested}/{model_grows.size} have data; median 1 replicate):")
    print(f"  agree: {n_agree}/{n_tested}   model false-positives: {n_fp}"
          + (f"  ({', '.join(fp_list)})" if fp_list else ""))
    growth_df['biolog_frac'] = growth_df.apply(
        lambda r: biolog_ordered.loc[r['strain_short'], r['carbon_source']]
        if (biolog_ordered is not None and r['strain_short'] in biolog_ordered.index) else np.nan, axis=1)
except Exception as e:
    print(f"  BioLog comparison skipped: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# 3. pFBA metabolic costs (Q12, Q14) — yield + ATP turnover + network flux load
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== 3. pFBA metabolic costs ===")
pfba_rows = []
for strain_id, model in models.items():
    for cs_name, cs_met in CARBON_SOURCES.items():
        with model:
            setup_condition(model, cs_name, cs_met)
            try:
                sol = pfba(model)
                total_flux = sol.fluxes.abs().sum()
                growth = sol.fluxes.get(BIOMASS_RXN, 0.0)
                atp = atp_turnover(model, sol.fluxes)
                o2 = abs(sol.fluxes.get(O2_EXCHANGE, float('nan')))
                status = sol.status
            except Exception as e:
                print(f"  pFBA failed {STRAIN_SHORT[strain_id]}/{cs_name}: {e}")
                total_flux, growth, atp, o2, status = (float('nan'),)*4 + ('failed',)
        pfba_rows.append({
            'strain': strain_id, 'strain_short': STRAIN_SHORT[strain_id],
            'carbon_source': cs_name, 'total_flux': total_flux, 'growth_rate': growth,
            'atp_turnover': atp, 'o2_uptake': o2,
            'yield_per_cmol': (growth / CARBON_CMOL_RATE) if growth else 0.0,
            'atp_per_biomass': (atp / growth) if growth else float('nan'),
            'o2_per_biomass':  (o2 / growth) if growth else float('nan'),
            'status': status,
        })

pfba_df = pd.DataFrame(pfba_rows)
print("pFBA biomass yield per C-mol (gDW per mmol C):")
print(pfba_df.pivot(index='strain_short', columns='carbon_source',
                    values='yield_per_cmol').reindex(columns=CS_ORDER).round(5).to_string())
print("\npFBA ATP turnover per biomass (mmol ATP / gDW):")
print(pfba_df.pivot(index='strain_short', columns='carbon_source',
                    values='atp_per_biomass').reindex(columns=CS_ORDER).round(1).to_string())


# ═══════════════════════════════════════════════════════════════════════════════
# 4. FVA flux variability (Q12, Q14) — LOOPLESS (thermodynamically attainable)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== 4. FVA flux variability (loopless=True) ===")
fva_data = {}
for strain_id, model in models.items():
    fva_data[strain_id] = {}
    for cs_name, cs_met in CARBON_SOURCES.items():
        with model:
            setup_condition(model, cs_name, cs_met)
            sol = model.optimize()
            if sol.status != 'optimal' or sol.objective_value < 1e-6:
                print(f"  Skip FVA {STRAIN_SHORT[strain_id]}/{cs_name}: no growth")
                continue
            print(f"  FVA {STRAIN_SHORT[strain_id]}/{cs_name} (growth={sol.objective_value:.4f})...",
                  end=' ', flush=True)
            fva = flux_variability_analysis(model, fraction_of_optimum=0.9,
                                            processes=1, loopless=True)
            fva['range'] = fva['maximum'] - fva['minimum']
            fva_data[strain_id][cs_name] = fva
            print('done')
print("FVA complete.")

fva_rows = []
for strain_id, cs_dict in fva_data.items():
    for cs_name, fva in cs_dict.items():
        tmp = fva.copy()
        tmp['strain'] = strain_id
        tmp['strain_short'] = STRAIN_SHORT[strain_id]
        tmp['carbon_source'] = cs_name
        tmp.index.name = 'reaction'
        fva_rows.append(tmp.reset_index())
fva_long = pd.concat(fva_rows, ignore_index=True) if fva_rows else pd.DataFrame()
print(f"FVA long-format rows: {len(fva_long)}")

# Standard-vs-loopless comparison on RIB40 (demonstrates loop removal; Q5/Q12).
print("\n  Standard-vs-loopless FVA comparison (RIB40)...", end=' ', flush=True)
loop_compare_rows = []
rib = models[TARGET_STRAINS[0]]
for cs_name, cs_met in CARBON_SOURCES.items():
    with rib:
        setup_condition(rib, cs_name, cs_met)
        if rib.optimize().objective_value < 1e-6:
            continue
        std = flux_variability_analysis(rib, fraction_of_optimum=0.9, processes=1, loopless=False)
        loo = flux_variability_analysis(rib, fraction_of_optimum=0.9, processes=1, loopless=True)
    std_r = std['maximum'] - std['minimum']
    loo_r = loo['maximum'] - loo['minimum']
    for rid in std_r.index:
        loop_compare_rows.append({
            'reaction': rid, 'carbon_source': cs_name,
            'range_standard': std_r[rid], 'range_loopless': loo_r[rid],
            'loop_inflation': std_r[rid] - loo_r[rid],
        })
loop_compare = pd.DataFrame(loop_compare_rows)
if not loop_compare.empty:
    n_loop = int((loop_compare['loop_inflation'] > 1.0).sum())
    print(f"done. {n_loop} reaction-conditions had loop-inflated standard ranges (>1 unit).")


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Accessory reactions (Q13) — flux-active on LOOPLESS ranges
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== 5. Accessory reactions (loopless flux activity) ===")
all_rxn_ids = set()
for m in models.values():
    all_rxn_ids.update(r.id for r in m.reactions)
presence = pd.DataFrame(
    {STRAIN_SHORT[sid]: {r.id: 1 for r in m.reactions} for sid, m in models.items()},
    index=sorted(all_rxn_ids)).fillna(0).astype(int)
n_strains = len(models)
rxn_counts = presence.sum(axis=1)
# Core/accessory follow the 187-strain pan-genome definition (core = present in >95% of the 187
# strain-specific models), to stay consistent with the manuscript (Figure 2). The classification is
# a static table derived once from strain-GEMs_automated_187.pickle (code/derive_reaction_classification.py).
_cls187 = pd.read_csv(os.path.join(PROJECT, 'data', 'genome', 'reaction_classification_187.csv'))
CORE187 = set(_cls187.loc[_cls187['class'] == 'core', 'reaction'])
ACC187  = set(_cls187.loc[_cls187['class'] == 'accessory', 'reaction'])
RXN_N187 = dict(zip(_cls187['reaction'], _cls187['n_strains']))
core_rxns      = sorted(set(all_rxn_ids) & CORE187)
accessory_rxns = sorted(set(all_rxn_ids) & ACC187)
n_unclassified = len(set(all_rxn_ids) - CORE187 - ACC187)
print(f"Core reactions (>95% of 187 strains): {len(core_rxns)}")
print(f"Accessory reactions (<=95% of 187 strains): {len(accessory_rxns)}")
print(f"Reactions in the 8 models absent from the 187 classification (v3 gap-fills etc.): {n_unclassified}")

acc_summary = pd.DataFrame()
if not fva_long.empty:
    rxn_info = {}
    for m in models.values():
        for r in m.reactions:
            rxn_info.setdefault(r.id, {'name': r.name, 'subsystem': r.subsystem})
    acc_fva = fva_long[fva_long['reaction'].isin(accessory_rxns)].copy()
    active_acc = acc_fva[acc_fva['range'] > 1e-6]
    if not active_acc.empty:
        acc_summary = (active_acc.groupby(['reaction', 'carbon_source'])
                       .agg(n_active_strains=('strain', 'nunique'), max_range=('range', 'max'))
                       .reset_index())
        acc_summary['n_strains_present'] = acc_summary['reaction'].map(RXN_N187)  # of 187
        acc_summary['name']      = acc_summary['reaction'].map(lambda x: rxn_info.get(x, {}).get('name', ''))
        acc_summary['subsystem'] = acc_summary['reaction'].map(lambda x: rxn_info.get(x, {}).get('subsystem', ''))
    n_active = acc_summary['reaction'].nunique() if not acc_summary.empty else 0
    print(f"Flux-active accessory reactions (loopless range > 0 on >=1 substrate): "
          f"{n_active} / {len(accessory_rxns)}")

# 5b. Accessory single-deletion growth impact — the direct "impact on growth" test
# (FVA range > 0 only shows a reaction *can* carry flux, not that it affects growth).
print("\n=== 5b. Accessory single-deletion growth impact (glucose) ===")
acc_set = set(accessory_rxns)
acc_del_rows = []
for strain_id, model in models.items():
    with model:
        setup_condition(model, 'Glucose', CARBON_SOURCES['Glucose'])
        base = model.slim_optimize()
        for r in [x for x in model.reactions if x.id in acc_set]:
            lb, ub = r.lower_bound, r.upper_bound
            r.lower_bound, r.upper_bound = 0.0, 0.0
            ko = model.slim_optimize()
            r.lower_bound, r.upper_bound = lb, ub
            ret = (ko / base) if (base and ko == ko) else float('nan')
            acc_del_rows.append({
                'strain': STRAIN_SHORT[strain_id], 'reaction': r.id,
                'base_growth': base, 'ko_growth': ko, 'growth_retained': ret,
                'pct_drop': 100 * (1 - ret) if ret == ret else float('nan')})
acc_del_df = pd.DataFrame(acc_del_rows)
if not acc_del_df.empty:
    n_impact = int((acc_del_df['pct_drop'] > 1).sum())
    print(f"Single deletions evaluated: {len(acc_del_df)} (strain x accessory present)")
    print(f"  reducing growth >1%: {n_impact};  max single-deletion drop: "
          f"{acc_del_df['pct_drop'].max():.3f}%")
    print("  -> no accessory reaction individually changes the growth optimum"
          if n_impact == 0 else "  -> some accessory reactions affect growth (see CSV)")


# ═══════════════════════════════════════════════════════════════════════════════
# 6. GAM/NGAM sensitivity (Q13)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== 6. GAM/NGAM sensitivity ===")
ref_model = models[TARGET_STRAINS[0]]
print(f"GAM (ATP in {BIOMASS_ASM}): {GAM_ATP} mmol ATP/gDW")
print(f"NGAM ({NGAM_RXN}) imposed value: {NGAM_VALUE} mmol ATP/gDW/h "
      f"(shipped default 0; literature range ~1-3.6 for filamentous fungi)")

ngam_sweep_rows = []
for ngam in np.linspace(0, 10, 21):
    for cs_name, cs_met in CARBON_SOURCES.items():
        with ref_model:
            setup_condition(ref_model, cs_name, cs_met, ngam=ngam)
            ref_model.objective = ref_model.reactions.get_by_id(BIOMASS_RXN)
            sol = ref_model.optimize()
            growth = sol.objective_value if sol.status == 'optimal' else 0.0
        ngam_sweep_rows.append({'ngam': ngam, 'carbon_source': cs_name, 'growth_rate': growth})
ngam_df = pd.DataFrame(ngam_sweep_rows)
# slope per substrate
for cs in CS_ORDER:
    g = ngam_df[ngam_df['carbon_source'] == cs].sort_values('ngam')
    if len(g) > 1:
        slope = np.polyfit(g['ngam'], g['growth_rate'], 1)[0]
        print(f"  {cs:9s} d(growth)/d(NGAM) = {slope:.4f} h-1 per mmol ATP/gDW/h")


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Subsystem activity (Q14)
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== 7. Subsystem activity ===")
sub_rows = []
for strain_id, model in models.items():
    for cs_name, cs_met in CARBON_SOURCES.items():
        with model:
            setup_condition(model, cs_name, cs_met)
            try:
                sol = pfba(model)
            except Exception as e:
                print(f"  pFBA failed {STRAIN_SHORT[strain_id]}/{cs_name}: {e}")
                continue
            for rxn in model.reactions:
                flux = abs(sol.fluxes.get(rxn.id, 0))
                if flux < 1e-9:
                    continue
                sub_rows.append({
                    'strain': strain_id, 'strain_short': STRAIN_SHORT[strain_id],
                    'carbon_source': cs_name, 'reaction': rxn.id,
                    'subsystem': rxn.subsystem if rxn.subsystem else 'Unknown',
                    'abs_flux': flux})
sub_df = pd.DataFrame(sub_rows)
print(f"Non-zero flux entries: {len(sub_df)}")

sub_agg = (sub_df.groupby(['strain_short', 'carbon_source', 'subsystem'])['abs_flux']
           .sum().reset_index())
total_by_cond = (sub_df.groupby(['strain_short', 'carbon_source'])['abs_flux']
                 .sum().rename('total_flux').reset_index())
sub_agg = sub_agg.merge(total_by_cond, on=['strain_short', 'carbon_source'])
sub_agg['fraction'] = sub_agg['abs_flux'] / sub_agg['total_flux']

# Unknown-subsystem flux fraction (completeness report; M-sub)
unk = sub_agg[sub_agg['subsystem'] == 'Unknown'].groupby('carbon_source')['fraction'].mean()
print("Mean 'Unknown'-subsystem flux fraction by substrate:")
for cs in CS_ORDER:
    print(f"  {cs:9s} {unk.get(cs, 0)*100:.1f}%")

max_frac = sub_agg.groupby('subsystem')['fraction'].max()
active_subsystems = max_frac[max_frac > 0.005].sort_values(ascending=False).index.tolist()
print(f"Active subsystems (>0.5%): {len(active_subsystems)}")


# ═══════════════════════════════════════════════════════════════════════════════
# 8. Export CSVs
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== 8. Export ===")
pfba_econ = pfba_df.copy()
pfba_econ['flux_per_growth'] = pfba_econ['total_flux'] / pfba_econ['growth_rate'].replace(0, np.nan)

exports = [
    (egc_df,      'strain_comparison_EGC.csv',         True),
    (growth_df,   'strain_comparison_FBA.csv',         False),
    (pfba_df,     'strain_comparison_pFBA.csv',        False),
    (pfba_econ,   'strain_comparison_pFBA_efficiency.csv', False),
    (fva_long,    'strain_comparison_FVA.csv',         False),
    (loop_compare,'strain_comparison_FVA_loopless_vs_standard.csv', False),
    (sub_agg,     'strain_comparison_subsystems.csv',  False),
    (ngam_df,     'strain_comparison_NGAM.csv',        False),
    (acc_summary, 'strain_comparison_accessory.csv',   False),
    (acc_del_df,  'strain_comparison_accessory_deletion.csv', False),
]
if biolog_ordered is not None:
    exports.append((biolog_ordered, 'strain_comparison_biolog.csv', True))
if egc_before:
    eb = pd.DataFrame([egc_before], index=['Pan_oryzae_before_fix'])
    exports.append((eb, 'strain_comparison_EGC_before_fix.csv', True))
for df, fname, use_index in exports:
    path = os.path.join(OUT_DIR, fname)
    if df is not None and not df.empty:
        df.to_csv(path, index=use_index)
        print(f"  Saved: {fname}")
    else:
        print(f"  Skipped (empty): {fname}")


# ═══════════════════════════════════════════════════════════════════════════════
# 9. Visualizations
# ═══════════════════════════════════════════════════════════════════════════════
print("\n=== 9. Visualizations ===")

rxn_meta = {}
for _m in models.values():
    for _r in _m.reactions:
        rxn_meta.setdefault(_r.id, {'name': _r.name or _r.id,
                                    'subsystem': _r.subsystem or 'Unknown'})
rxn_meta_df = pd.DataFrame.from_dict(rxn_meta, orient='index')
rxn_meta_df.index.name = 'reaction_id'
rxn_meta_df.to_csv(os.path.join(OUT_DIR, 'strain_comparison_reaction_info.csv'))

def _rxn_label(rxn_id, max_name=42):
    name = rxn_meta.get(rxn_id, {}).get('name', '')
    return f"{rxn_id}  {name[:max_name]}" if name and name != rxn_id else rxn_id

def _save(fig, stem):
    fig.savefig(os.path.join(OUT_DIR, stem + '.png'), dpi=150, bbox_inches='tight')
    fig.savefig(os.path.join(OUT_DIR, stem + '.pdf'), bbox_inches='tight')
    plt.close(fig)


# ── 9a. Figures: all-reaction loopless-FVA heatmap (supplementary) and the four-panel
# manuscript figure (A growth, B pFBA optimum, C leave-one-out, D substrate consistency).
# These are the only two figures produced; earlier exploratory figures were removed. ────────
print("\n  Manuscript + supplementary figures...")

# Single-reaction deletion on RIB40 (panel C: core vs accessory essentiality).
_rib = models[TARGET_STRAINS[0]]
with _rib:
    setup_condition(_rib, 'Glucose', CARBON_SOURCES['Glucose'])
    _base = _rib.slim_optimize()
    _del = single_reaction_deletion(_rib, processes=1)
_del = _del.reset_index(drop=True)
_del['rid'] = _del['ids'].apply(lambda s: list(s)[0])
_del['ret'] = (100 * _del['growth'] / _base).fillna(0).clip(lower=0)
_core_ret = _del[_del['rid'].map(lambda r: r in CORE187)]['ret'].values
_acc_ret = (acc_del_df['growth_retained'] * 100).fillna(0).clip(lower=0).values if not acc_del_df.empty else np.array([])
_n_core_ess = int((_core_ret < 1).sum()); _n_acc_ess = int((_acc_ret < 99).sum())

# Substrate palette + the Figure-5 strain palette/order (so panel A reads next to Figure 5).
_CS_COL = {'Glucose': '#2c3e50', 'Glycerol': '#e67e22', 'Maltose': '#16a085', 'Xylose': '#2980b9'}
_STRAIN_COL = {'NRRL_2217': '#9cc078', 'NRRL_3483': '#90c0e4', 'NRRL_3488': '#e49c9c',
               'NRRL_35890': '#f0d884', 'NRRL_471': '#cca8cc', 'NRRL_5589': '#f0b478',
               'NRRL_5592': '#d8b49c', 'RIB40': '#90cccc'}
_FIG5_ORDER = ['NRRL_2217', 'NRRL_3483', 'NRRL_3488', 'NRRL_35890', 'NRRL_471',
               'NRRL_5589', 'NRRL_5592', 'RIB40']
# ggplot theme_bw aesthetic (white panel, light grid, thin grey border) matching Figure 5.
_GGPLOT_BW = {
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'axes.edgecolor': '#777777',
    'axes.linewidth': 0.8, 'axes.grid': True, 'grid.color': '#dcdcdc', 'grid.linewidth': 0.7,
    'axes.axisbelow': True, 'axes.spines.top': True, 'axes.spines.right': True,
    'axes.spines.left': True, 'axes.spines.bottom': True, 'xtick.color': '#4d4d4d',
    'ytick.color': '#4d4d4d', 'text.color': '#1a1a1a', 'axes.labelcolor': '#1a1a1a'}

def _strain_legend_handles():
    from matplotlib.lines import Line2D
    return [Line2D([0], [0], marker='s', ls='', ms=9, color=_STRAIN_COL[s], label=s.replace('_', ' '))
            for s in _FIG5_ORDER]

def _draw_growth(ax):
    """Panel A: FBA-predicted maximum growth rate per strain, grouped by carbon source and
    coloured as in Figure 5. Bars within a substrate are near-identical (strain CV < 0.1%) — the
    models reach the same growth rate, in contrast to the experimental variation in Figure 5."""
    gv = (growth_df.pivot(index='strain_short', columns='carbon_source', values='growth_rate')
          .reindex(index=_FIG5_ORDER, columns=CS_ORDER))
    n = len(_FIG5_ORDER); width = 0.105
    for j, s in enumerate(_FIG5_ORDER):
        xs = [i + (j - (n - 1) / 2) * width for i in range(len(CS_ORDER))]
        ax.bar(xs, [gv.loc[s, c] for c in CS_ORDER], width=width, color=_STRAIN_COL[s],
               edgecolor='white', linewidth=0.4, zorder=3)
    ax.set_xticks(range(len(CS_ORDER))); ax.set_xticklabels(CS_ORDER)
    ax.set_ylabel('Predicted growth rate (h$^{-1}$)'); ax.set_ylim(0, gv.values.max() * 1.55)
    ax.grid(axis='x', linewidth=0)
    ax.legend(handles=_strain_legend_handles(), loc='upper center', ncol=4, fontsize=7,
              frameon=True, framealpha=.9, columnspacing=.8, handletextpad=.3, borderpad=.4)

def _draw_pfba(ax, metric='flux_per_growth', ylabel='Total pFBA flux per unit growth\n(mmol gDW$^{-1}$ h$^{-1}$ per h$^{-1}$)'):
    """Each substrate is a tight cluster of 8 strain points: the parsimonious optimum is
    set by the carbon source, not by accessory-genome differences between strains."""
    rng = np.random.default_rng(0)
    for xi, cs in enumerate(CS_ORDER):
        d = pfba_econ[pfba_econ['carbon_source'] == cs]
        ax.scatter(xi + (rng.random(len(d)) - .5) * .18, d[metric], s=55, c=_CS_COL[cs],
                   alpha=.85, edgecolors='white', lw=.5, zorder=3)
    ax.set_xticks(range(len(CS_ORDER))); ax.set_xticklabels(CS_ORDER); ax.set_xlim(-.5, len(CS_ORDER) - .5)
    ax.set_ylabel(ylabel); ax.margins(y=.16)

def _draw_loo(ax):
    rng = np.random.default_rng(0)
    ax.scatter(0 + (rng.random(len(_core_ret)) - .5) * .5, _core_ret, s=20, c='#34495e', alpha=.5, edgecolors='none')
    ax.scatter(1 + (rng.random(len(_acc_ret)) - .5) * .5, _acc_ret, s=20, c='#8e44ad', alpha=.5, edgecolors='none')
    ax.axhline(100, color='#27ae60', lw=1.2, ls='--')
    ax.set_xticks([0, 1]); ax.set_xticklabels([f'core (>95% of 187)\n(RIB40, n={len(_core_ret)})', f'accessory (≤95%)\n(8 strains, n={len(_acc_ret)})'])
    ax.set_ylabel('Growth retained after\nsingle-reaction deletion (%)'); ax.set_ylim(-5, 108); ax.set_xlim(-.6, 1.6)

def _draw_consist(ax, ntop=7):
    """Subsystems on x, Δ-vs-glucose on y; each of 8 strains is a point (jittered)."""
    pv = sub_agg.pivot_table(index=['strain_short', 'subsystem'], columns='carbon_source', values='fraction').reset_index()
    for c in ['Glycerol', 'Maltose', 'Xylose']:
        pv[c + '_d'] = pv[c] - pv['Glucose']
    grp = pv.groupby('subsystem')[['Glycerol_d', 'Maltose_d', 'Xylose_d']]
    ranked = grp.mean().abs().max(axis=1).sort_values(ascending=False)
    ranked = ranked[~ranked.index.isin(['Unknown', ''])]      # drop the uninterpretable bucket
    ts = ranked.head(ntop).index.tolist()
    sd = grp.std().loc[ts].max().max() * 100
    cm = {'Glycerol': '#e67e22', 'Maltose': '#16a085', 'Xylose': '#2980b9'}; off = {'Glycerol': -.24, 'Maltose': 0, 'Xylose': .24}
    rng = np.random.default_rng(1)
    for xi, ss in enumerate(ts):
        d = pv[pv['subsystem'] == ss]
        for sn in ['Glycerol', 'Maltose', 'Xylose']:
            yv = d[sn + '_d'].values * 100; xj = xi + off[sn] + (rng.random(len(yv)) - .5) * .13
            ax.scatter(xj, yv, s=50, c=cm[sn], alpha=.85, edgecolors='white', lw=.5)
    ax.axhline(0, color='#888', lw=1)
    ax.set_xticks(range(len(ts))); ax.set_xticklabels([s[:18] for s in ts], rotation=32, ha='right', fontsize=8)
    ax.set_ylabel('Δ subsystem flux share\nvs glucose (percentage points)')
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], marker='o', ls='', color=cm[s], label=s) for s in cm],
              fontsize=8.5, loc='upper right', title='substrate')

# Supplementary heatmap data: per-strain loopless FVA range (glucose), by subsystem.
_sg = fva_long[fva_long['carbon_source'] == 'Glucose']
_Rfull = _sg.pivot(index='reaction', columns='strain_short', values='range').reindex(columns=STRAINS_ORDERED).fillna(0.0)
_R = _Rfull[_Rfull.max(axis=1) > 1e-6]                        # reactions that carry flux somewhere
_subof = pd.Series({r: rxn_meta.get(r, {}).get('subsystem', 'Unknown') for r in _R.index})
_ord = pd.DataFrame({'s': _subof, 'm': _R.mean(axis=1)}).sort_values(['s', 'm'], ascending=[True, False]).index
_A = _R.reindex(_ord); _Asub = _subof.reindex(_ord)
_TICKS = [0, 10, 100, 1000, 2000]

def _draw_A(ax):
    """Log-scaled FVA-range heatmap of all flux-carrying reactions, grouped by subsystem."""
    im = ax.imshow(np.log10(_A.values + 1), cmap='viridis', aspect='auto', interpolation='nearest',
                   vmin=0, vmax=np.log10(2001))
    ax.set_xticks(range(len(STRAINS_ORDERED))); ax.set_xticklabels(STRAINS_ORDERED, rotation=40, ha='right', fontsize=8)
    ax.set_yticks([])
    b = 0
    for s, grp in _Asub.groupby(_Asub, sort=False):
        n = len(grp); ax.axhline(b - .5, color='white', lw=.6)
        if n >= 7:
            ax.text(-0.8, b + n / 2, str(s)[:30], va='center', ha='right', fontsize=6.5)
        b += n
    fa = [i for i, s in enumerate(_Asub) if s == 'Fatty acid biosynthesis']
    if fa:
        ax.annotate('fatty-acid\nbiosynthesis', xy=(7.6, (fa[0] + fa[-1]) / 2), xytext=(9.4, (fa[0] + fa[-1]) / 2),
                    fontsize=8, va='center', color='#c0392b', arrowprops=dict(arrowstyle='-[', color='#c0392b', lw=1.4))
    return im

def _cbar(fig, im, ax):
    cb = fig.colorbar(im, ax=ax, shrink=.4, label='flux range (mmol/gDW/h)')
    cb.set_ticks([np.log10(t + 1) for t in _TICKS]); cb.set_ticklabels(_TICKS)

# Supplementary figure: all-reaction loopless-FVA range heatmap.
if not fva_long.empty:
    fig, ax = plt.subplots(figsize=(7.2, 10)); imv = _draw_A(ax)
    ax.set_title(f'FVA flux range is near-identical across strains for most reactions\n'
                 f'{len(_A)} reactions with non-zero range (rows, by subsystem)', fontsize=10.5)
    _cbar(fig, imv, ax); _save(fig, 'fig_fva_allrxn_heatmap')

# Manuscript figure (Figure-5 ggplot aesthetic): A growth | B pFBA / C leave-one-out / D consistency.
if not pfba_econ.empty and not acc_del_df.empty and not growth_df.empty:
    with plt.rc_context(_GGPLOT_BW):
        fig = plt.figure(figsize=(14, 11))
        outer = fig.add_gridspec(2, 1, hspace=0.34)
        top = outer[0].subgridspec(1, 2, width_ratios=[.52, .48], wspace=.2)
        bot = outer[1].subgridspec(1, 2, width_ratios=[.34, .66], wspace=.22)
        axA = fig.add_subplot(top[0]); _draw_growth(axA)
        axA.set_title('A', loc='left', fontsize=17, fontweight='bold')
        axB = fig.add_subplot(top[1]); _draw_pfba(axB)
        axB.set_title('B', loc='left', fontsize=17, fontweight='bold')
        axC = fig.add_subplot(bot[0]); _draw_loo(axC)
        axC.set_title('C', loc='left', fontsize=17, fontweight='bold')
        axD = fig.add_subplot(bot[1]); _draw_consist(axD)
        axD.set_title('D', loc='left', fontsize=17, fontweight='bold')
        _save(fig, 'fig_mechanistic_manuscript')
    print("  Saved: fig_fva_allrxn_heatmap, fig_mechanistic_manuscript")

print("\n=== DONE ===")
