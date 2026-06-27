"""
Rebuild Section 9 (Enhanced visualizations) of p5_strain_comparison.ipynb with
embedded outputs. Run AFTER run_p5_analysis.py has produced all CSVs/PNGs.

The cells carry both runnable code and pre-computed outputs (figures embedded as
base64 display_data, summary tables as stream text) so the notebook renders the
results without a 30-min FVA re-run, while remaining reproducible if executed.
"""
import json, base64
import pandas as pd
import numpy as np
from pathlib import Path

PROJECT = Path(__file__).parent.parent
OUT_DIR = PROJECT / 'data' / 'intermediate'
NB_PATH = PROJECT / 'code' / 'p5_strain_comparison.ipynb'
CS_ORDER = ['Glucose', 'Glycerol', 'Maltose', 'Xylose']


# ── helpers ────────────────────────────────────────────────────────────────────
def _b64(fname):
    p = OUT_DIR / fname
    return base64.b64encode(p.read_bytes()).decode('utf-8') if p.exists() else None

def _stream(text_lines):
    return {'output_type': 'stream', 'name': 'stdout',
            'text': text_lines if isinstance(text_lines, list) else [text_lines]}

def _png(fname, desc='<Figure>'):
    b = _b64(fname)
    return {'output_type': 'display_data', 'metadata': {},
            'data': {'image/png': b, 'text/plain': [desc]}} if b else None

def _code(lines, outputs=None):
    src = lines if isinstance(lines, str) else '\n'.join(lines)
    return {'cell_type': 'code', 'execution_count': None,
            'metadata': {'section9': True}, 'source': src, 'outputs': outputs or []}

def _md(text):
    return {'cell_type': 'markdown', 'metadata': {'section9': True}, 'source': text}

def _read_csv(fname, **kw):
    p = OUT_DIR / fname
    try:
        return pd.read_csv(p, **kw)
    except Exception:
        return None


# ── build cells ────────────────────────────────────────────────────────────────
cells = []

cells.append(_md(
    '## 9. Figures\n\n'
    'Two figures are produced by `run_p5_analysis.py` from the CSVs in `data/intermediate/`: the '
    'four-panel manuscript figure and the supplementary all-reaction loopless-FVA heatmap. '
    'Methodological choices addressing the reviewer: carbon sources fed at **equal C-mol** supply, '
    'non-growth ATP maintenance (**NGAM**) imposed, and flux variability analysis run **loopless** so '
    'reported ranges are thermodynamically attainable.'))

cells.append(_code(
    "import pandas as pd, numpy as np\n"
    "from IPython.display import Image, display\n"
    "OUT='../data/intermediate/'\n"
    "print('Loaded')",
    outputs=[_stream(['Loaded\n'])]))

# 9a manuscript figure (4-panel)
cells.append(_md(
    '### 9a. Manuscript figure\n\n'
    'The proposed manuscript figure (`fig_mechanistic_manuscript`), styled to match Figure 5 (ggplot '
    'aesthetic, shared strain colours), shows that the strain models do not differ appreciably in the '
    'solution space probed by the tested growth conditions: **(A)** FBA predicts a near-identical '
    'maximum growth rate for every strain on each carbon source (strain CV < 0.1%), in direct contrast '
    'to the experimental variation in Figure 5; **(B)** the parsimonious (pFBA) optimum is set by the '
    'carbon source, not the strain (strain CV < 0.3%); **(C)** single-reaction deletion — 169 essential '
    'core reactions vs 0/634 accessory deletions affecting growth, so the accessory genome is '
    'dispensable for the optimum; **(D)** substrate-specific subsystem rewiring with the 8 strains '
    'overlapping tightly (max inter-strain SD ≈ 3 percentage points).'))
cells.append(_code(
    "display(Image(OUT+'fig_mechanistic_manuscript.png'))",
    outputs=[o for o in [_png('fig_mechanistic_manuscript.png', '<manuscript figure (4-panel)>')] if o]))

# 9b supplementary all-reaction loopless-FVA heatmap
cells.append(_md(
    '### 9b. Supplementary figure — all-reaction loopless-FVA range heatmap\n\n'
    'The loopless flux range of every flux-carrying reaction (827 reactions with a non-zero range), '
    'grouped by subsystem and shown on a log scale: the eight strains are visually indistinguishable '
    'for the large majority of reactions, with the few differences concentrated in fatty-acid '
    'biosynthesis and mitochondrial transport — reactions that are inactive at the parsimonious '
    'optimum, which is why the growth and pFBA predictions remain near-identical across strains.'))
cells.append(_code(
    "display(Image(OUT+'fig_fva_allrxn_heatmap.png'))",
    outputs=[o for o in [_png('fig_fva_allrxn_heatmap.png', '<all-reaction FVA heatmap>')] if o]))

# ── write notebook ─────────────────────────────────────────────────────────────
nb = json.loads(NB_PATH.read_text(encoding='utf-8'))

def _is_sec9(cell):
    if cell.get('metadata', {}).get('section9'):       # cells we created (idempotent)
        return True
    src = ''.join(cell.get('source', ''))
    # legacy heuristics (cells from earlier versions without the metadata tag)
    if cell.get('cell_type') == 'code' and ('OUT+' in src or 'OUT +' in src or 'Image(OUT' in src):
        return True   # Section 9 code cells use the bare `OUT` var; sections 0-8 use OUT_DIR
    return ('## 9.' in src or '### 9' in src or 'fig_EGC' in src or
            'fig_summary_panel' in src or 'fig_FVA_loopless' in src or
            'fig_growth_yield' in src or
            src.strip().startswith("import pandas as pd, numpy"))

existing = [i for i, c in enumerate(nb['cells']) if _is_sec9(c)]
for i in sorted(existing, reverse=True):
    nb['cells'].pop(i)
if existing:
    print(f"Removed {len(existing)} pre-existing Section 9 cells")

nb['cells'].extend(cells)
NB_PATH.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding='utf-8')
print(f"Wrote {len(cells)} Section 9 cells to {NB_PATH.name}")
n_png = sum(1 for c in cells for o in c.get('outputs', [])
            if o.get('output_type') == 'display_data')
print(f"Embedded {n_png} figures")
