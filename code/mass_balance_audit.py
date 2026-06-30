"""
Mass-balance audit of the panAsp v3 model collection.

Classifies every internal reaction as balanced / imbalanced and, for the
imbalanced ones, separates:
  * imbalances explained by the two metabolite-annotation bugs
    (C00080 proton annotated 'H2'/charge 0; C00007 O2 with no formula), and
  * genuine stoichiometric errors that remain after those annotations are
    corrected in-memory (extra/missing protons, or carbon-skeleton C/N/O/P/S).

Outputs:
  data/intermediate/mass_balance_audit.csv      (per-reaction detail)
  data/intermediate/fig_mass_balance_audit.png  (summary)
The audit does NOT modify the saved model — annotation fixes are in-memory only,
to reveal the true picture.
"""
import os, warnings
from pickle import load
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(PROJECT, 'data', 'intermediate')
PICKLE  = os.path.join(PROJECT, 'model', 'pAo_strain-GEMs_validated.pickle')

with open(PICKLE, 'rb') as f:
    all_models = load(f)
pan = next(m for m in all_models if m.id == 'Pan_oryzae')


def _src(rid):
    if rid.startswith('OtherAsp'): return 'OtherAsp (heterologous)'
    if rid.startswith('Afu'):      return 'Afu (A. fumigatus)'
    if rid.startswith('Gap'):      return 'Gapfill'
    if rid and rid[0] == 'r' and rid[1:].split('_')[0].isdigit(): return 'native (template)'
    return 'other'


# ── 1. metabolites with no usable formula (block balance assessment) ───────────
noform = [m for m in pan.metabolites if not m.formula]
noform_ids = {m.id for m in noform}
print(f"Metabolites with no formula: {len(noform)}")
base_no = {}
for m in noform:
    base = m.id.split('[')[0]
    base_no[base] = base_no.get(base, 0) + 1
for base, n in sorted(base_no.items(), key=lambda x: -x[1])[:15]:
    nm = next((x.name for x in noform if x.id.split('[')[0] == base), '')
    print(f"   {base:10s} x{n:3d}  {nm}")

# ── 2. in-memory annotation fixes (proton + O2) to reveal the true picture ─────
for m in pan.metabolites:
    base = m.id.split('[')[0]
    if base == 'C00080':            # proton: 'H2'/0  ->  'H'/+1
        m.formula, m.charge = 'H', 1
    elif base == 'C00007':          # molecular O2: no formula -> 'O2'
        m.formula = 'O2'
PROTON = 'C00080'
O2     = 'C00007'

# refresh no-formula set after fixes
noform_ids = {m.id for m in pan.metabolites if not m.formula}

# ── 3. classify every internal reaction ───────────────────────────────────────
rows = []
for r in pan.reactions:
    if r.boundary:
        continue
    bases = {m.id.split('[')[0] for m in r.metabolites}
    has_noform = any((m.id in noform_ids) for m in r.metabolites)
    try:
        bal = r.check_mass_balance()
    except ValueError:
        rows.append({'reaction': r.id, 'source': _src(r.id),
                     'subsystem': r.subsystem or 'Unknown',
                     'status': 'unassessable (no-formula metabolite)',
                     'imbalance': '', 'category': 'unassessable',
                     'has_proton': PROTON in bases, 'has_O2': O2 in bases})
        continue
    if not bal:
        continue
    keys = {k for k in bal if abs(bal[k]) > 1e-9}
    proton_like = keys <= {'H', 'charge'}
    o2_like     = keys <= {'O', 'H', 'charge'} and O2 in bases
    if proton_like:
        cat = 'proton (H/charge) stoichiometry'
    elif keys <= {'O'} or o2_like:
        cat = 'oxygen'
    elif keys & {'C', 'N', 'P', 'S'}:
        cat = 'carbon-skeleton (C/N/P/S)'
    else:
        cat = 'other'
    rows.append({'reaction': r.id, 'source': _src(r.id),
                 'subsystem': r.subsystem or 'Unknown',
                 'status': 'imbalanced',
                 'imbalance': '; '.join(f'{k}:{bal[k]:+g}' for k in sorted(keys)),
                 'category': cat,
                 'has_proton': PROTON in bases, 'has_O2': O2 in bases})

audit = pd.DataFrame(rows)
audit.to_csv(os.path.join(OUT_DIR, 'mass_balance_audit.csv'), index=False)

n_internal = sum(1 for r in pan.reactions if not r.boundary)
n_flagged  = len(audit)
print(f"\nInternal reactions: {n_internal}")
print(f"Flagged (imbalanced or unassessable) AFTER fixing proton+O2 annotations: {n_flagged}")
print(f"  balanced: {n_internal - n_flagged} ({100*(n_internal-n_flagged)/n_internal:.0f}%)")

print("\nBy category:")
print(audit['category'].value_counts().to_string())
print("\nImbalanced (assessable) by source:")
print(audit[audit['status'] == 'imbalanced']['source'].value_counts().to_string())
print("\nProton-stoichiometry errors by source:")
ps = audit[audit['category'] == 'proton (H/charge) stoichiometry']
print(ps['source'].value_counts().to_string())
print("\nTop subsystems among genuine carbon-skeleton errors:")
cs = audit[audit['category'] == 'carbon-skeleton (C/N/P/S)']
print(cs['subsystem'].value_counts().head(10).to_string())

# ── 4. figure ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(15, 6))
cat_order = ['proton (H/charge) stoichiometry', 'carbon-skeleton (C/N/P/S)',
             'oxygen', 'other', 'unassessable']
cat_counts = audit['category'].value_counts().reindex(cat_order).fillna(0)
colors = ['#e67e22', '#c0392b', '#2980b9', '#7f8c8d', '#bdc3c7']
b = axes[0].barh(range(len(cat_counts)), cat_counts.values, color=colors, edgecolor='white')
axes[0].set_yticks(range(len(cat_counts)))
axes[0].set_yticklabels([c.replace(' ', '\n', 1) for c in cat_counts.index], fontsize=9)
axes[0].bar_label(b, fmt='%.0f', padding=3)
axes[0].set_xlabel('Number of internal reactions')
axes[0].set_title(f'Imbalance category after correcting the proton + O₂ annotations\n'
                  f'({n_internal - n_flagged}/{n_internal} balanced; {n_flagged} flagged)', fontsize=10)
axes[0].invert_yaxis()

# stacked by source for the imbalanced (assessable) reactions
imb = audit[audit['status'] == 'imbalanced']
piv = imb.pivot_table(index='source', columns='category', values='reaction',
                      aggfunc='count', fill_value=0)
piv = piv.reindex(columns=[c for c in cat_order if c in piv.columns])
piv.plot(kind='barh', stacked=True, ax=axes[1],
         color={'proton (H/charge) stoichiometry': '#e67e22',
                'carbon-skeleton (C/N/P/S)': '#c0392b', 'oxygen': '#2980b9', 'other': '#7f8c8d'})
axes[1].set_xlabel('Number of imbalanced reactions')
axes[1].set_ylabel('')
axes[1].set_title('Imbalanced reactions by source and category', fontsize=10)
axes[1].legend(fontsize=8, title='')
fig.suptitle('Mass-balance audit — Pan_oryzae (proton/O₂ annotations corrected in-memory)',
             y=1.02, fontsize=12, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, 'fig_mass_balance_audit.png'), dpi=150, bbox_inches='tight')
fig.savefig(os.path.join(OUT_DIR, 'fig_mass_balance_audit.pdf'), bbox_inches='tight')
plt.close(fig)
print("\nSaved: mass_balance_audit.csv, fig_mass_balance_audit.{png,pdf}")
