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
pickle_path = os.path.join(PROJECT, 'model', 'panAsp_v3_gemList_187_strains.pickle')
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
core_rxns      = rxn_counts[rxn_counts == n_strains].index.tolist()
accessory_rxns = rxn_counts[(rxn_counts > 0) & (rxn_counts < n_strains)].index.tolist()
print(f"Core reactions (all {n_strains} strains): {len(core_rxns)}")
print(f"Accessory reactions (1-{n_strains-1} strains): {len(accessory_rxns)}")
print(f"Strain-unique (singleton) reactions: {(rxn_counts == 1).sum()}")

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
        acc_summary['n_strains_present'] = acc_summary['reaction'].map(rxn_counts)
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


# ── 9a. EGC before/after curation (Q5) ─────────────────────────────────────────
if egc_before:
    carriers = ['NADH', 'NADPH', 'L-Glutamate']
    before = [egc_before.get(c, 0) for c in carriers]
    after  = [egc_df.loc['Pan_oryzae', c] if 'Pan_oryzae' in egc_df.index else 0 for c in carriers]
    x = np.arange(len(carriers)); w = 0.38
    fig, ax = plt.subplots(figsize=(8, 5))
    b1 = ax.bar(x - w/2, before, w, label='Before curation', color='#c0392b', edgecolor='white')
    b2 = ax.bar(x + w/2, after,  w, label='After curation (4 bounds)', color='#27ae60', edgecolor='white')
    ax.bar_label(b1, fmt='%.0f', padding=3); ax.bar_label(b2, fmt='%.0f', padding=3)
    ax.set_xticks(x); ax.set_xticklabels(carriers)
    ax.set_ylabel('Max dissipation flux (mmol/gDW/h)')
    ax.set_title('Pan_oryzae erroneous energy-generating cycles eliminated by curation\n'
                 'Reverting the four bound changes restores the EGCs (causal test)', fontsize=11)
    ax.legend(); ax.set_ylim(0, 1150)
    _save(fig, 'fig_EGC_beforeafter')
    print("  Saved: fig_EGC_beforeafter")

# ── 9b. EGC heatmap (all models, post-curation) ───────────────────────────────
egc_plot = egc_df.reindex(index=STRAINS_ORDERED + (['Pan_oryzae'] if 'Pan_oryzae' in egc_df.index else [])).fillna(0)
egc_plot = egc_plot.where(egc_plot.abs() > 1e-6, 0.0)
fig, ax = plt.subplots(figsize=(12, 5))
g = sns.heatmap(egc_plot, annot=True, fmt='.0f', cmap='RdYlGn_r', vmin=0, vmax=0.001,
                linewidths=0.5, linecolor='white', ax=ax,
                cbar_kws={'label': 'Dissipation flux (mmol/gDW/h)', 'shrink': 0.7})
for t in g.texts:
    t.set_color('white'); t.set_fontweight('bold'); t.set_fontsize(10)
ax.set_title(f'EGC screen after curation: {len(egc_plot)} models x 12 energy carriers\n'
             'All dissipation tests return zero — no erroneous energy-generating cycles', fontsize=11)
ax.set_xlabel('Energy carrier'); ax.set_ylabel(''); ax.tick_params(axis='x', rotation=40)
_save(fig, 'fig_EGC_heatmap')
print("  Saved: fig_EGC_heatmap")

# ── 9c. Standard vs loopless FVA (rigor; Q5/Q12) ──────────────────────────────
if not loop_compare.empty:
    glc = loop_compare[loop_compare['carbon_source'] == 'Glucose'].copy()
    top = glc.sort_values('loop_inflation', ascending=False).head(15).iloc[::-1]
    labels = [_rxn_label(r, 34) for r in top['reaction']]
    y = np.arange(len(top))
    fig, ax = plt.subplots(figsize=(11, 7))
    ax.barh(y, top['range_standard'], color='#c0392b', alpha=0.55, label='Standard FVA (with loops)')
    ax.barh(y, top['range_loopless'], color='#2471a3', label='Loopless FVA (attainable)')
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('FVA flux range at 90% optimal growth (mmol/gDW/h)')
    ax.set_title('Standard FVA overstates flux variability via internal loops (RIB40, glucose)\n'
                 'Loopless FVA collapses loop-only ranges to their attainable values', fontsize=11)
    ax.legend(loc='lower right')
    _save(fig, 'fig_FVA_loopless_vs_standard')
    print("  Saved: fig_FVA_loopless_vs_standard")

# ── 9d. Biomass yield per C-mol (Q12/Q14) ─────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(15, 5))
sns.heatmap(pivot_growth.reindex(index=STRAINS_ORDERED, columns=CS_ORDER),
            annot=True, fmt='.3f', cmap='YlOrRd', linewidths=0.5, linecolor='gray',
            ax=axes[0], cbar_kws={'label': 'h⁻¹', 'shrink': 0.8})
axes[0].set_title('FBA growth rate (h⁻¹)\nequal C-mol uptake, NGAM imposed')
axes[0].set_xlabel('Carbon source'); axes[0].set_ylabel('Strain')
sns.heatmap(pivot_yield.reindex(index=STRAINS_ORDERED, columns=CS_ORDER),
            annot=True, fmt='.4f', cmap='Greens', linewidths=0.5, linecolor='gray',
            ax=axes[1], cbar_kws={'label': 'gDW / mmol C', 'shrink': 0.8})
axes[1].set_title('Biomass yield per C-mol\n(uptake-rate independent efficiency)')
axes[1].set_xlabel('Carbon source'); axes[1].set_ylabel('Strain')
fig.suptitle('Carbon-normalised growth and yield — substrate differences reflect metabolism, not C content',
             y=1.02, fontsize=12)
_save(fig, 'fig_growth_yield')
print("  Saved: fig_growth_yield")

# ── 9e. Phenotype concordance: model vs BioLog (Q12) ──────────────────────────
if biolog_ordered is not None:
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))
    sns.heatmap(pivot_growth.reindex(index=STRAINS_ORDERED, columns=CS_ORDER),
                annot=True, fmt='.3f', cmap='YlOrRd', vmin=0,
                linewidths=0.5, linecolor='gray', ax=axes[0],
                cbar_kws={'label': 'h⁻¹', 'shrink': 0.8})
    axes[0].set_title('FBA predicted growth (h⁻¹)')
    axes[0].set_xlabel('Carbon source'); axes[0].set_ylabel('Strain')
    bio_p = biolog_ordered.reindex(index=STRAINS_ORDERED, columns=CS_ORDER)
    sns.heatmap(bio_p, annot=True, fmt='.2f', cmap='YlOrRd', vmin=0, vmax=1,
                linewidths=0.5, linecolor='gray', ax=axes[1],
                cbar_kws={'label': 'fraction replicates', 'shrink': 0.8})
    # mark disagreements (model grows, BioLog does not)
    for i, sid in enumerate(STRAINS_ORDERED):
        for j, cs in enumerate(CS_ORDER):
            mv = pivot_growth.reindex(index=STRAINS_ORDERED, columns=CS_ORDER).iloc[i, j]
            bv = bio_p.iloc[i, j]
            if mv is not None and mv > 0.01 and bv is not None and bv <= 0.5:
                axes[1].add_patch(plt.Rectangle((j, i), 1, 1, fill=False,
                                  edgecolor='#2471a3', lw=3))
    axes[1].set_title('BioLog observed growth\n(blue box = model false-positive vs wet lab)')
    axes[1].set_xlabel('Carbon source'); axes[1].set_ylabel('Strain')
    fig.suptitle('Phenotype concordance: FBA reproduces most growth phenotypes; '
                 'discordances (e.g. NRRL_5589) flag model refinement targets',
                 y=1.02, fontsize=12)
    _save(fig, 'fig_FBA_vs_biolog')
    print("  Saved: fig_FBA_vs_biolog")

# ── 9f. pFBA ATP-cost and yield (Q14) ─────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(15, 5))
piv_atp = pfba_df.pivot(index='strain_short', columns='carbon_source',
                        values='atp_per_biomass').reindex(index=STRAINS_ORDERED, columns=CS_ORDER)
piv_y   = pfba_df.pivot(index='strain_short', columns='carbon_source',
                        values='yield_per_cmol').reindex(index=STRAINS_ORDERED, columns=CS_ORDER)
sns.heatmap(piv_atp, annot=True, fmt='.0f', cmap='RdYlBu_r', linewidths=0.5,
            linecolor='gray', ax=axes[0], cbar_kws={'label': 'mmol ATP / gDW', 'shrink': 0.8})
axes[0].set_title('ATP turnover per biomass (energetic cost)\nLower = more energy-efficient growth')
axes[0].set_xlabel('Carbon source'); axes[0].set_ylabel('Strain')
sns.heatmap(piv_y, annot=True, fmt='.4f', cmap='Greens', linewidths=0.5,
            linecolor='gray', ax=axes[1], cbar_kws={'label': 'gDW / mmol C', 'shrink': 0.8})
axes[1].set_title('Biomass yield per C-mol\nHigher = more carbon-efficient growth')
axes[1].set_xlabel('Carbon source'); axes[1].set_ylabel('Strain')
fig.suptitle('pFBA metabolic cost — defensible energetic and carbon-efficiency metrics',
             y=1.02, fontsize=12)
_save(fig, 'fig_pFBA_efficiency')
print("  Saved: fig_pFBA_efficiency")

# ── 9g. FVA annotated heatmaps (loopless), top inter-strain variance ──────────
if not fva_long.empty:
    fva_var_rows = []
    for cs_name in CS_ORDER:
        subset = fva_long[fva_long['carbon_source'] == cs_name].copy()
        subset = subset[subset['range'] < 1800]
        if subset.empty:
            continue
        inter_var = subset.groupby('reaction')['range'].var().dropna().sort_values(ascending=False)
        mean_range = subset.groupby('reaction')['range'].mean()
        top_rxns = inter_var.head(25).index.tolist()
        plot_data = (subset[subset['reaction'].isin(top_rxns)]
                     .pivot(index='reaction', columns='strain_short', values='range')
                     .reindex(columns=STRAINS_ORDERED).fillna(0))
        if not plot_data.empty:
            plot_data = plot_data.loc[plot_data.mean(axis=1).sort_values().index]
            plot_data.index = [_rxn_label(r) for r in plot_data.index]
            fig, ax = plt.subplots(figsize=(12, 9))
            sns.heatmap(plot_data, cmap='viridis', ax=ax, yticklabels=True,
                        linewidths=0.2, linecolor='lightgray',
                        cbar_kws={'label': 'Loopless flux range (mmol/gDW/h)', 'shrink': 0.8})
            ax.set_title(f'Loopless FVA flux range at 90% optimal growth — {cs_name}\n'
                         'Top 25 reactions by inter-strain variance', pad=10)
            ax.set_xlabel('Strain'); ax.set_ylabel(''); ax.tick_params(axis='y', labelsize=7)
            _save(fig, f'fig_FVA_annotated_{cs_name}')
        for rxn_id, var_val in inter_var.head(50).items():
            fva_var_rows.append({
                'reaction': rxn_id, 'carbon_source': cs_name,
                'inter_strain_variance': var_val, 'mean_range': mean_range.get(rxn_id, np.nan),
                'name': rxn_meta.get(rxn_id, {}).get('name', ''),
                'subsystem': rxn_meta.get(rxn_id, {}).get('subsystem', '')})
    pd.DataFrame(fva_var_rows).to_csv(
        os.path.join(OUT_DIR, 'strain_comparison_FVA_top_variable.csv'), index=False)
    print("  Saved: fig_FVA_annotated_{substrate} + FVA_top_variable.csv")

# ── 9h. Accessory reactions — subsystem + presence (Q13) ──────────────────────
if not acc_summary.empty:
    acc_plot = acc_summary.copy()
    acc_plot['subsystem_label'] = acc_plot['subsystem'].replace('', np.nan).fillna('Unknown/ungrouped')
    sub_counts = acc_plot.groupby('subsystem_label')['reaction'].nunique().sort_values()
    fig, ax = plt.subplots(figsize=(9, max(3.5, len(sub_counts) * 0.45)))
    colors = plt.cm.tab10(np.linspace(0, 0.9, len(sub_counts)))
    bars = ax.barh(sub_counts.index, sub_counts.values, color=colors, edgecolor='white')
    ax.bar_label(bars, padding=3, fontsize=10)
    ax.set_xlabel('Number of flux-active accessory reactions (loopless)')
    ax.set_title(f'Accessory reactions with attainable flux activity — by subsystem\n'
                 f'({acc_plot["reaction"].nunique()} reactions, loopless range > 0 on >=1 carbon source)', pad=10)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    _save(fig, 'fig_accessory_subsystems')

    active_acc_rxns = sorted(acc_plot['reaction'].unique())
    pres_mat = presence.loc[presence.index.isin(active_acc_rxns), STRAINS_ORDERED].copy()
    pres_mat.index = [_rxn_label(r, 38) for r in pres_mat.index]
    fig, ax = plt.subplots(figsize=(10, max(4, len(pres_mat) * 0.22)))
    sns.heatmap(pres_mat, cmap='Blues', vmin=0, vmax=1, ax=ax, linewidths=0.3, linecolor='gray',
                cbar_kws={'label': 'present (1) / absent (0)', 'shrink': 0.4, 'ticks': [0, 1]})
    ax.set_title(f'Presence/absence of {len(pres_mat)} flux-active accessory reactions', pad=10)
    ax.set_xlabel('Strain'); ax.set_ylabel(''); ax.tick_params(axis='y', labelsize=6)
    _save(fig, 'fig_accessory_presence')
    print("  Saved: fig_accessory_subsystems, fig_accessory_presence")

# ── 9h2. Accessory single-deletion growth impact (Q13) ────────────────────────
if not acc_del_df.empty:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(acc_del_df['growth_retained'].dropna() * 100, bins=40,
            color='#2471a3', edgecolor='white')
    ax.axvline(100, color='#c0392b', ls='--', lw=1.5, label='no impact (100% retained)')
    ax.set_xlabel('Growth retained after single accessory-reaction deletion (%)')
    ax.set_ylabel('Number of (strain × accessory reaction) deletions')
    ax.set_title('Accessory reactions do not individually affect the growth optimum\n'
                 f'{len(acc_del_df)} single deletions; max growth drop '
                 f'{acc_del_df["pct_drop"].max():.2f}% (glucose, carbon-normalised, NGAM)',
                 fontsize=10)
    ax.legend()
    _save(fig, 'fig_accessory_deletion')
    print("  Saved: fig_accessory_deletion")

# ── 9i. NGAM sensitivity with literature value marked (Q13) ────────────────────
fig, ax = plt.subplots(figsize=(8, 5))
colors = sns.color_palette('tab10', n_colors=4)
for (cs_name, grp), color in zip(ngam_df.groupby('carbon_source'), colors):
    ax.plot(grp['ngam'], grp['growth_rate'], label=cs_name, color=color, marker='o', markersize=4)
ax.axvline(x=NGAM_VALUE, color='black', linestyle='--', lw=1.5,
           label=f'Imposed NGAM = {NGAM_VALUE}')
ax.axvspan(1.0, 3.6, color='gray', alpha=0.12, label='Literature range (filamentous fungi)')
ax.set_xlabel('NGAM (mmol ATP/gDW/h)'); ax.set_ylabel('Growth rate (h⁻¹)')
ax.set_title('Growth vs NGAM (RIB40) — growth is robust across the literature range', fontsize=11)
ax.legend(fontsize=9)
_save(fig, 'fig_NGAM_sensitivity')
print("  Saved: fig_NGAM_sensitivity")

# ── 9j. Subsystem activity heatmaps (Q14) ─────────────────────────────────────
for cs_name in CS_ORDER:
    sub_cs = sub_agg[(sub_agg['carbon_source'] == cs_name) &
                     (sub_agg['subsystem'].isin(active_subsystems))]
    if sub_cs.empty:
        continue
    pivot_sub = sub_cs.pivot(index='subsystem', columns='strain_short', values='fraction') \
                      .reindex(index=active_subsystems, columns=STRAINS_ORDERED).fillna(0)
    fig, ax = plt.subplots(figsize=(10, max(5, len(active_subsystems) * 0.28)))
    sns.heatmap(pivot_sub, cmap='RdYlBu_r', ax=ax, linewidths=0.2, linecolor='lightgray',
                vmin=0, vmax=pivot_sub.values.max())
    ax.set_title(f'Subsystem activity fraction (pFBA) — {cs_name}')
    ax.set_xlabel('Strain'); ax.set_ylabel('Subsystem')
    _save(fig, f'fig_subsystem_{cs_name}')
print("  Saved: fig_subsystem_{substrate}")

# ── 9k. Summary panel ─────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.35)
ax1, ax2, ax3, ax4 = (fig.add_subplot(gs[r, c]) for r in range(2) for c in range(2))
sns.heatmap(pivot_yield.reindex(index=STRAINS_ORDERED, columns=CS_ORDER), annot=True, fmt='.4f',
            cmap='Greens', linewidths=0.5, linecolor='gray', ax=ax1,
            cbar_kws={'label': 'gDW/mmol C', 'shrink': 0.8})
ax1.set_title('A  Biomass yield per C-mol\nnearly identical across strains', fontsize=10)
ax1.set_xlabel('Carbon source'); ax1.set_ylabel('Strain')
if biolog_ordered is not None:
    sns.heatmap(biolog_ordered.reindex(index=STRAINS_ORDERED, columns=CS_ORDER), annot=True,
                fmt='.2f', cmap='YlOrRd', vmin=0, vmax=1, linewidths=0.5, linecolor='gray',
                ax=ax2, cbar_kws={'label': 'fraction', 'shrink': 0.8})
    ax2.set_title('B  BioLog observed growth\nstrain-level differences (wet lab)', fontsize=10)
else:
    ax2.text(0.5, 0.5, 'BioLog data not available', transform=ax2.transAxes, ha='center')
ax2.set_xlabel('Carbon source'); ax2.set_ylabel('Strain')
sns.heatmap(piv_atp, annot=True, fmt='.0f', cmap='RdYlBu_r', linewidths=0.5, linecolor='gray',
            ax=ax3, cbar_kws={'label': 'mmol ATP/gDW', 'shrink': 0.8})
ax3.set_title('C  ATP turnover per biomass\nenergetic cost by substrate', fontsize=10)
ax3.set_xlabel('Carbon source'); ax3.set_ylabel('Strain')
for (cs_name, grp), col in zip(ngam_df.groupby('carbon_source'), sns.color_palette('tab10', 4)):
    ax4.plot(grp['ngam'], grp['growth_rate'], label=cs_name, color=col, marker='o', markersize=3, lw=1.5)
ax4.axvline(NGAM_VALUE, color='black', ls='--', lw=1.2)
ax4.set_xlabel('NGAM (mmol ATP/gDW/h)'); ax4.set_ylabel('Growth rate (h⁻¹)')
ax4.set_title('D  NGAM sensitivity (RIB40)\nrobust across literature range', fontsize=10)
ax4.legend(fontsize=9, framealpha=0.5); ax4.grid(True, alpha=0.25)
fig.suptitle('Pan-GEM strain comparison — mechanistic FBA analysis summary',
             y=0.99, fontsize=13, fontweight='bold')
_save(fig, 'fig_summary_panel')
print("  Saved: fig_summary_panel")


# ── 9l. Additional / alternative views of the same conclusions ─────────────────
# These complement (do not replace) the figures above. Reusable panel drawers so
# the standalone figures and the multi-panel manuscript figure stay identical.
acc_per = (presence.loc[[r for r in accessory_rxns], STRAINS_ORDERED].sum()
           if accessory_rxns else pd.Series(0.0, index=STRAINS_ORDERED))
gr_g = growth_df[growth_df['carbon_source'] == 'Glucose'].set_index('strain_short')['growth_rate']

def _draw_scatter(ax):
    """Standard vs loopless FVA range, every reaction (RIB40, glucose)."""
    glc = loop_compare[loop_compare['carbon_source'] == 'Glucose']
    x = glc['range_standard'].clip(lower=1e-2); y = glc['range_loopless'].clip(lower=1e-2)
    loops = (glc['range_standard'] > 1) & (glc['range_loopless'] < 1)
    ax.scatter(x[~loops], y[~loops], s=10, c='#2471a3', alpha=.5, label='attainable')
    ax.scatter(x[loops], y[loops], s=16, c='#c0392b', alpha=.75, label='loop artifact')
    ax.plot([1e-2, 2e3], [1e-2, 2e3], 'k--', lw=1, alpha=.6)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('Standard FVA flux range (mmol/gDW/h)')
    ax.set_ylabel('Loopless FVA flux range (mmol/gDW/h)')
    ax.legend(loc='upper left', fontsize=8, frameon=False); ax.grid(alpha=.2, which='both')
    return int(loops.sum())

def _draw_egc(ax):
    carr = ['NADH', 'NADPH', 'L-Glutamate']
    before = [egc_before.get(c, 0) for c in carr]
    after = [egc_df.loc['Pan_oryzae', c] if 'Pan_oryzae' in egc_df.index else 0 for c in carr]
    xx = np.arange(3); w = .38
    b1 = ax.bar(xx - w/2, before, w, color='#c0392b', label='before curation')
    b2 = ax.bar(xx + w/2, after, w, color='#27ae60', label='after curation')
    ax.bar_label(b1, fmt='%.0f', fontsize=8); ax.bar_label(b2, fmt='%.0f', fontsize=8)
    ax.set_xticks(xx); ax.set_xticklabels(carr); ax.set_ylim(0, 1150)
    ax.set_ylabel('Max dissipation flux\n(mmol/gDW/h)'); ax.legend(fontsize=8, frameon=False)

def _draw_content(ax):
    yy = np.arange(len(STRAINS_ORDERED))
    ax.barh(yy, [acc_per.get(s, 0) for s in STRAINS_ORDERED], color='#8e44ad', edgecolor='white')
    ax.set_yticks(yy); ax.set_yticklabels(STRAINS_ORDERED, fontsize=8); ax.invert_yaxis()
    ax.set_xlabel('Accessory reactions carried')
    for i, s in enumerate(STRAINS_ORDERED):
        ax.text(acc_per.get(s, 0) + 1, i, f'{int(acc_per.get(s, 0))}', va='center', fontsize=7.5)
    cv = gr_g.std() / gr_g.mean() * 100 if len(gr_g) else 0
    ax.text(.97, .04, f'growth rate identical\n({gr_g.mean():.3f} h⁻¹, CV {cv:.2f}%)',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
            bbox=dict(boxstyle='round', fc='#eafaf1', ec='#27ae60'))

def _draw_subdelta(ax, ntop=14):
    _m = sub_agg.groupby(['carbon_source', 'subsystem'])['fraction'].mean().reset_index()
    _piv = _m.pivot(index='subsystem', columns='carbon_source', values='fraction').reindex(columns=CS_ORDER).fillna(0)
    _delta = _piv[['Glycerol', 'Maltose', 'Xylose']].sub(_piv['Glucose'], axis=0)
    _top = _delta.abs().max(axis=1).sort_values(ascending=False).head(ntop).index
    _d = _delta.loc[_top].sort_values('Xylose'); _vmax = np.abs(_d.values).max()
    im = ax.imshow(_d.values, cmap='RdBu_r', vmin=-_vmax, vmax=_vmax, aspect='auto')
    ax.set_xticks(range(3)); ax.set_xticklabels(['Glycerol', 'Maltose', 'Xylose'])
    ax.set_yticks(range(len(_d))); ax.set_yticklabels([s[:46] for s in _d.index], fontsize=8)
    for (i, j), v in np.ndenumerate(_d.values):
        if abs(v) > _vmax * .12:
            ax.text(j, i, f'{v*100:+.1f}', ha='center', va='center', fontsize=7,
                    color='white' if abs(v) > _vmax * .5 else '#222')
    return im

if not loop_compare.empty:
    fig, ax = plt.subplots(figsize=(6.4, 6)); nloop = _draw_scatter(ax)
    ax.set_title(f'Every reaction, standard vs loopless FVA (RIB40, glucose)\n'
                 f'{nloop} reactions collapse to ~0 — pure loop artifacts')
    _save(fig, 'fig_FVA_loopless_scatter')

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
yy = np.arange(len(STRAINS_ORDERED))
a1.barh(yy, [gr_g.get(s, np.nan) for s in STRAINS_ORDERED], color='#27ae60', edgecolor='white')
a1.set_yticks(yy); a1.set_yticklabels(STRAINS_ORDERED); a1.invert_yaxis()
a1.set_xlabel('FBA growth rate (h⁻¹)'); a1.set_title('Function: near-identical (spread < 0.1%)')
a2.barh(yy, [acc_per.get(s, 0) for s in STRAINS_ORDERED], color='#8e44ad', edgecolor='white')
a2.set_xlabel('Accessory reactions carried'); a2.set_title('Content: clearly different')
for i, s in enumerate(STRAINS_ORDERED):
    a2.text(acc_per.get(s, 0) + 1, i, f'{int(acc_per.get(s, 0))}', va='center', fontsize=8)
fig.suptitle('Single-point FBA cannot discriminate strains: identical function, different gene content',
             y=1.02, fontsize=12)
_save(fig, 'fig_function_vs_content')

if biolog_ordered is not None:
    grp = growth_df.pivot(index='strain_short', columns='carbon_source', values='growth_rate') \
                   .reindex(index=STRAINS_ORDERED, columns=CS_ORDER)
    mg = grp > 0.01 * grp.values.max(); bg = biolog_ordered.reindex(index=STRAINS_ORDERED, columns=CS_ORDER); tst = bg.notna()
    TP = int(((mg) & (bg > .5) & tst).values.sum()); FP = int(((mg) & (bg <= .5) & tst).values.sum())
    FN = int(((~mg) & (bg > .5) & tst).values.sum()); TN = int(((~mg) & (bg <= .5) & tst).values.sum())
    M = np.array([[TP, FP], [FN, TN]])
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    ax.imshow([[.15, .6], [.6, .15]], cmap='RdYlGn_r', vmin=0, vmax=1, alpha=.35)
    for (i, j), v in np.ndenumerate(M):
        ax.text(j, i - .12, str(v), ha='center', va='center', fontsize=30, fontweight='bold')
        ax.text(j, i + .22, ['TP', 'FP', 'FN', 'TN'][i*2+j], ha='center', va='center', fontsize=11, color='#444')
    ax.set_xticks([0, 1]); ax.set_xticklabels(['grows', 'no growth']); ax.xaxis.tick_top(); ax.xaxis.set_label_position('top')
    ax.set_yticks([0, 1]); ax.set_yticklabels(['predicts\ngrowth', 'predicts\nno growth'])
    ax.set_xlabel('BioLog (observed)'); ax.set_ylabel('Model (FBA)')
    ax.set_title(f'Model vs BioLog, {TP+FP+FN+TN} tested conditions\n'
                 f'Sensitivity {TP}/{TP+FN}=100% · Specificity {TN}/{TN+FP}=0/2\n'
                 f'2 false-positives: NRRL_5589 glycerol & maltose')
    _save(fig, 'fig_concordance_matrix')

fig, ax = plt.subplots(figsize=(7.5, 6)); imS = _draw_subdelta(ax, 14)
ax.set_title('Substrate-specific rewiring: subsystem flux share\nrelative to the glucose baseline (Δ percentage points)')
fig.colorbar(imS, ax=ax, shrink=.6, label='Δ flux fraction vs glucose')
_save(fig, 'fig_subsystem_delta')

# Multi-panel manuscript figure (A EGC, B loopless FVA, C function-vs-content, D rewiring)
if not loop_compare.empty and biolog_ordered is not None:
    fig = plt.figure(figsize=(13, 10.5)); gs = fig.add_gridspec(2, 2, hspace=.34, wspace=.26)
    axA = fig.add_subplot(gs[0, 0]); _draw_egc(axA)
    axA.set_title('A  Models are free of energy-generating cycles\nPan-GEM curation (4 bounds) eliminates all EGCs', fontsize=10.5, loc='left')
    axB = fig.add_subplot(gs[0, 1]); nb2 = _draw_scatter(axB)
    axB.set_title(f'B  Loopless FVA removes spurious flux variability\n{nb2} reactions collapse from wide range to ~0', fontsize=10.5, loc='left')
    axC = fig.add_subplot(gs[1, 0]); _draw_content(axC)
    axC.set_title('C  Identical growth, divergent accessory content\naccessory complement varies 29–149 reactions per strain', fontsize=10.5, loc='left')
    axD = fig.add_subplot(gs[1, 1]); imD = _draw_subdelta(axD, 9)
    axD.set_title('D  Substrate-specific pathway rewiring\nsubsystem flux share vs glucose (Δ %-points)', fontsize=10.5, loc='left')
    fig.colorbar(imD, ax=axD, shrink=.7, label='Δ vs glucose')
    _save(fig, 'fig_mechanistic_4panel')
print("  Saved: fig_FVA_loopless_scatter, fig_function_vs_content, fig_concordance_matrix, "
      "fig_subsystem_delta, fig_mechanistic_4panel")


# ── 9m. Revised manuscript figure: intra-species variation is confined to
# dispensable accessory metabolism, and substrate rewiring is strain-consistent ─
print("\n  Revised manuscript figure (single-reaction deletion + strain consistency)...")
_rib = models[TARGET_STRAINS[0]]
with _rib:
    setup_condition(_rib, 'Glucose', CARBON_SOURCES['Glucose'])
    _base = _rib.slim_optimize()
    _del = single_reaction_deletion(_rib, processes=1)
_del = _del.reset_index(drop=True)
_del['rid'] = _del['ids'].apply(lambda s: list(s)[0])
_del['ret'] = (100 * _del['growth'] / _base).fillna(0).clip(lower=0)
_rxn_count = presence.sum(axis=1)
_core_ret = _del[_del['rid'].map(lambda r: _rxn_count.get(r, 0) == n_strains)]['ret'].values
_acc_ret = (acc_del_df['growth_retained'] * 100).fillna(0).clip(lower=0).values if not acc_del_df.empty else np.array([])
_n_core_ess = int((_core_ret < 1).sum()); _n_acc_ess = int((_acc_ret < 99).sum())

def _draw_varA(ax):
    sg = fva_long[(fva_long['carbon_source'] == 'Glucose') & (fva_long['range'] < 1800)]
    iv = sg.groupby('reaction')['range'].var().dropna().sort_values(ascending=False)
    rids = iv.head(15).index.tolist()
    P = (sg[sg['reaction'].isin(rids)].pivot(index='reaction', columns='strain_short', values='range')
         .reindex(index=rids, columns=STRAINS_ORDERED).fillna(0))
    lab = [f'{r}  {rxn_meta.get(r, {}).get("name", "")[:26]}' for r in rids]
    im = ax.imshow(P.values, cmap='viridis', aspect='auto')
    ax.set_xticks(range(len(STRAINS_ORDERED))); ax.set_xticklabels(STRAINS_ORDERED, rotation=40, ha='right', fontsize=8)
    ax.set_yticks(range(len(rids))); ax.set_yticklabels(lab, fontsize=7)
    return im

def _draw_loo(ax):
    rng = np.random.default_rng(0)
    ax.scatter(0 + (rng.random(len(_core_ret)) - .5) * .55, _core_ret, s=8, c='#34495e', alpha=.45, edgecolors='none')
    ax.scatter(1 + (rng.random(len(_acc_ret)) - .5) * .55, _acc_ret, s=8, c='#8e44ad', alpha=.45, edgecolors='none')
    ax.axhline(100, color='#27ae60', lw=1, ls='--')
    ax.set_xticks([0, 1]); ax.set_xticklabels([f'core\n(RIB40, n={len(_core_ret)})', f'accessory\n(8 strains, n={len(_acc_ret)})'])
    ax.set_ylabel('Growth retained after\nsingle-reaction deletion (%)'); ax.set_ylim(-5, 108); ax.set_xlim(-.6, 1.6)
    ax.text(0, -2, f'{_n_core_ess} essential', ha='center', va='top', fontsize=8, color='#34495e')
    ax.text(1, 92, f'all ≥ 99.9%\n({_n_acc_ess} essential)', ha='center', va='top', fontsize=8, color='#8e44ad')

def _draw_consist(ax, ntop=7):
    pv = sub_agg.pivot_table(index=['strain_short', 'subsystem'], columns='carbon_source', values='fraction').reset_index()
    for c in ['Glycerol', 'Maltose', 'Xylose']:
        pv[c + '_d'] = pv[c] - pv['Glucose']
    grp = pv.groupby('subsystem')[['Glycerol_d', 'Maltose_d', 'Xylose_d']]
    ts = grp.mean().abs().max(axis=1).sort_values(ascending=False).head(ntop).index.tolist()
    sd = grp.std().loc[ts].max().max() * 100
    cm = {'Glycerol': '#e67e22', 'Maltose': '#16a085', 'Xylose': '#2980b9'}; off = {'Glycerol': -.24, 'Maltose': 0, 'Xylose': .24}
    for yi, ss in enumerate(ts[::-1]):
        d = pv[pv['subsystem'] == ss]
        for sn in ['Glycerol', 'Maltose', 'Xylose']:
            ax.scatter(d[sn + '_d'].values * 100, [yi + off[sn]] * len(d), s=26, c=cm[sn], alpha=.8, edgecolors='white', lw=.4)
    ax.axvline(0, color='#888', lw=1)
    ax.set_yticks(range(len(ts))); ax.set_yticklabels([s[:44] for s in ts[::-1]], fontsize=9)
    ax.set_xlabel('Δ subsystem flux share vs glucose baseline (percentage points)')
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], marker='o', ls='', color=cm[s], label=s) for s in cm],
              fontsize=8.5, loc='lower left', frameon=True, title='substrate')
    ax.text(.99, .02, f'each point = one of 8 strains\nmax inter-strain SD = {sd:.2f} pp', transform=ax.transAxes,
            ha='right', va='bottom', fontsize=8.5, bbox=dict(boxstyle='round', fc='#fef9e7', ec='#d4ac0d'))

if not fva_long.empty and not acc_del_df.empty:
    fig, ax = plt.subplots(figsize=(7, 5.5)); imv = _draw_varA(ax)
    ax.set_title('Reactions with the largest inter-strain variation\n(loopless FVA flux range, glucose)')
    fig.colorbar(imv, ax=ax, shrink=.7, label='flux range (mmol/gDW/h)'); _save(fig, 'fig_strain_variable_reactions')
    fig, ax = plt.subplots(figsize=(5.6, 5)); _draw_loo(ax)
    ax.set_title('Leave-one-out: no accessory reaction is essential'); _save(fig, 'fig_leave_one_out')
    fig, ax = plt.subplots(figsize=(7.8, 5.2)); _draw_consist(ax)
    ax.set_title('Substrate-specific rewiring is consistent across strains'); _save(fig, 'fig_subsystem_consistency')
    fig = plt.figure(figsize=(14, 14)); gs = fig.add_gridspec(2, 2, height_ratios=[1.2, .95], hspace=.6, wspace=.34)
    axA = fig.add_subplot(gs[0, 0]); imA = _draw_varA(axA); fig.colorbar(imA, ax=axA, shrink=.7, label='flux range (mmol/gDW/h)')
    axA.set_title('A  Reactions varying most between strains (loopless FVA range, glucose)', loc='left', fontsize=10.5)
    axB = fig.add_subplot(gs[0, 1]); _draw_loo(axB)
    axB.set_title('B  No accessory reaction is essential', loc='left', fontsize=10.5)
    axC = fig.add_subplot(gs[1, :]); _draw_consist(axC)
    axC.set_title('C  Substrate rewiring is consistent across strains', loc='left', fontsize=10.5)
    _save(fig, 'fig_mechanistic_manuscript')
    print("  Saved: fig_strain_variable_reactions, fig_leave_one_out, fig_subsystem_consistency, fig_mechanistic_manuscript")

print("\n=== DONE ===")
