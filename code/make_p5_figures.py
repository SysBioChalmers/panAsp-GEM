"""Stand-alone builder for the two manuscript figures, decoupled from the heavy analysis.

Reads the plot-ready data from the CSVs in data/intermediate/ (written by run_p5_analysis.py) so the
figure LAYOUT can be iterated in seconds without re-running FBA/FVA. The one model-dependent input
(panel C: RIB40 single-reaction deletion) is computed once and cached to
strain_comparison_core_deletion.csv.

Produces, in data/intermediate/:
  fig_mechanistic_manuscript.png / .svg   (A growth | B pFBA / C leave-one-out / D consistency)
  fig_fva_allrxn_heatmap.png   / .svg     (all-reaction loopless-FVA range, by subsystem)

Run:  python code/make_p5_figures.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(SCRIPT_DIR)
OUT = os.path.join(PROJECT, "data", "intermediate")

# ── tweakable layout knobs (iterate here) ──────────────────────────────────────
L = dict(
    # panel D (subsystem consistency) x-axis labels — wrapped to 2 lines, no truncation
    D_ntop=7, D_rotation=40, D_fontsize=8,
    # heatmap — narrow heatmap, full labels on the left, black boundary ticks sticking out left
    HM_figsize=(4.6, 14.5), HM_label_fontsize=5.0, HM_label_min_rows=4,
    HM_sep_color="black", HM_sep_lw=0.4, HM_tick=0.6,
    SAVE_PAD=0.18,
)


def _wrap2(s):
    """Split a long label into two lines at the space-or-slash nearest the middle."""
    if len(s) <= 16:
        return s
    mid = len(s) / 2
    cands = [i + 1 for i, c in enumerate(s) if c in " /" and 3 < i + 1 < len(s) - 3]
    if not cands:
        return s
    j = min(cands, key=lambda i: abs(i - mid))
    return s[:j].rstrip() + "\n" + s[j:].strip()

CS_ORDER = ["Glucose", "Glycerol", "Maltose", "Xylose"]
STRAINS_ORDERED = ["RIB40", "NRRL_2217", "NRRL_3483", "NRRL_3488",
                   "NRRL_5589", "NRRL_5592", "NRRL_35890", "NRRL_471"]
_CS_COL = {"Glucose": "#2c3e50", "Glycerol": "#e67e22", "Maltose": "#16a085", "Xylose": "#2980b9"}
_STRAIN_COL = {"NRRL_2217": "#9cc078", "NRRL_3483": "#90c0e4", "NRRL_3488": "#e49c9c",
               "NRRL_35890": "#f0d884", "NRRL_471": "#cca8cc", "NRRL_5589": "#f0b478",
               "NRRL_5592": "#d8b49c", "RIB40": "#90cccc"}
_FIG5_ORDER = ["NRRL_2217", "NRRL_3483", "NRRL_3488", "NRRL_35890", "NRRL_471",
               "NRRL_5589", "NRRL_5592", "RIB40"]
_GGPLOT_BW = {
    "figure.facecolor": "white", "axes.facecolor": "white", "axes.edgecolor": "#777777",
    "axes.linewidth": 0.8, "axes.grid": True, "grid.color": "#dcdcdc", "grid.linewidth": 0.7,
    "axes.axisbelow": True, "xtick.color": "#4d4d4d", "ytick.color": "#4d4d4d",
    "text.color": "#1a1a1a", "axes.labelcolor": "#1a1a1a"}


def _csv(name):
    return pd.read_csv(os.path.join(OUT, name))


# ── panel C cache: RIB40 single-reaction deletion (the only model-dependent input) ──────────────
def _ensure_core_deletion():
    path = os.path.join(OUT, "strain_comparison_core_deletion.csv")
    if os.path.exists(path):
        return pd.read_csv(path)
    print("  computing panel-C RIB40 single-reaction deletion (one-time, cached)...", flush=True)
    import warnings, pickle
    warnings.filterwarnings("ignore")
    from cobra.flux_analysis import single_reaction_deletion
    RIB40 = "Aspergillus_oryzae_RIB40_GCF_000184455.2"
    with open(os.path.join(PROJECT, "model", "pAo_strain-GEMs_validated.pickle"), "rb") as f:
        rib = {m.id: m for m in pickle.load(f)}[RIB40]
    cls = _csv_path("../data/genome/reaction_classification_187.csv")
    core = set(cls.loc[cls["class"] == "core", "reaction"])
    with rib:
        rib.objective = rib.reactions.get_by_id("r2359")
        rib.add_boundary(rib.metabolites.get_by_id("C00031[e]"), type="exchange", lb=-1.12, ub=1000)
        rib.reactions.get_by_id("r1901").lower_bound = 1.0
        base = rib.slim_optimize()
        d = single_reaction_deletion(rib, processes=1).reset_index(drop=True)
    d["rid"] = d["ids"].apply(lambda s: list(s)[0])
    d["retained"] = (100 * d["growth"] / base).fillna(0).clip(lower=0)
    out = d[d["rid"].map(lambda r: r in core)][["rid", "retained"]]
    out.to_csv(path, index=False)
    print(f"  cached {len(out)} core-reaction deletions -> {os.path.basename(path)}", flush=True)
    return out


def _csv_path(rel):
    return pd.read_csv(os.path.join(SCRIPT_DIR, rel))


# ── load all figure data ────────────────────────────────────────────────────────
growth_df = _csv("strain_comparison_FBA.csv")
pfba_econ = _csv("strain_comparison_pFBA_efficiency.csv")
sub_agg = _csv("strain_comparison_subsystems.csv")
fva_long = _csv("strain_comparison_FVA.csv")
rxn_info = _csv("strain_comparison_reaction_info.csv").set_index("reaction_id")
acc_del = _csv("strain_comparison_accessory_deletion.csv")
core_del = _ensure_core_deletion()

_core_ret = core_del["retained"].values
_acc_ret = (acc_del["growth_retained"] * 100).fillna(0).clip(lower=0).values


def _save(fig, stem):
    fig.savefig(os.path.join(OUT, stem + ".png"), dpi=150, bbox_inches="tight", pad_inches=L["SAVE_PAD"])
    try:
        fig.savefig(os.path.join(OUT, stem + ".svg"), bbox_inches="tight", pad_inches=L["SAVE_PAD"])
        msg = f"  saved {stem}.png / .svg"
    except PermissionError:
        msg = f"  saved {stem}.png  (SVG locked/open in a viewer — skipped)"
    plt.close(fig)
    print(msg, flush=True)


def _strain_legend_handles():
    from matplotlib.lines import Line2D
    return [Line2D([0], [0], marker="s", ls="", ms=9, color=_STRAIN_COL[s], label=s.replace("_", " "))
            for s in _FIG5_ORDER]


def _draw_growth(ax):
    gv = (growth_df.pivot(index="strain_short", columns="carbon_source", values="growth_rate")
          .reindex(index=_FIG5_ORDER, columns=CS_ORDER))
    n = len(_FIG5_ORDER); width = 0.105
    for j, s in enumerate(_FIG5_ORDER):
        xs = [i + (j - (n - 1) / 2) * width for i in range(len(CS_ORDER))]
        ax.bar(xs, [gv.loc[s, c] for c in CS_ORDER], width=width, color=_STRAIN_COL[s],
               edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_xticks(range(len(CS_ORDER))); ax.set_xticklabels(CS_ORDER)
    ax.set_ylabel("Predicted growth rate (h$^{-1}$)"); ax.set_ylim(0, gv.values.max() * 1.55)
    ax.grid(axis="x", linewidth=0)
    ax.legend(handles=_strain_legend_handles(), loc="upper center", ncol=4, fontsize=7,
              frameon=True, framealpha=.9, columnspacing=.8, handletextpad=.3, borderpad=.4)


def _draw_pfba(ax):
    rng = np.random.default_rng(0)
    for xi, cs in enumerate(CS_ORDER):
        d = pfba_econ[pfba_econ["carbon_source"] == cs]
        ax.scatter(xi + (rng.random(len(d)) - .5) * .18, d["flux_per_growth"], s=55, c=_CS_COL[cs],
                   alpha=.85, edgecolors="white", lw=.5, zorder=3)
    ax.set_xticks(range(len(CS_ORDER))); ax.set_xticklabels(CS_ORDER); ax.set_xlim(-.5, len(CS_ORDER) - .5)
    ax.set_ylabel("Total pFBA flux per unit growth\n(mmol gDW$^{-1}$ h$^{-1}$ per h$^{-1}$)"); ax.margins(y=.16)


def _draw_loo(ax):
    rng = np.random.default_rng(0)
    ax.scatter(0 + (rng.random(len(_core_ret)) - .5) * .5, _core_ret, s=20, c="#34495e", alpha=.5, edgecolors="none")
    ax.scatter(1 + (rng.random(len(_acc_ret)) - .5) * .5, _acc_ret, s=20, c="#8e44ad", alpha=.5, edgecolors="none")
    ax.axhline(100, color="#27ae60", lw=1.2, ls="--")
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"core (>95% of 187)\n(RIB40, n={len(_core_ret)})",
                        f"accessory (≤95%)\n(8 strains, n={len(_acc_ret)})"])
    ax.set_ylabel("Growth retained after\nsingle-reaction deletion (%)"); ax.set_ylim(-5, 108); ax.set_xlim(-.6, 1.6)


def _draw_consist(ax, ntop=None):
    ntop = ntop or L["D_ntop"]
    pv = sub_agg.pivot_table(index=["strain_short", "subsystem"], columns="carbon_source", values="fraction").reset_index()
    for c in ["Glycerol", "Maltose", "Xylose"]:
        pv[c + "_d"] = pv[c] - pv["Glucose"]
    grp = pv.groupby("subsystem")[["Glycerol_d", "Maltose_d", "Xylose_d"]]
    ranked = grp.mean().abs().max(axis=1).sort_values(ascending=False)
    ranked = ranked[~ranked.index.isin(["Unknown", ""])]
    ts = ranked.head(ntop).index.tolist()
    cm = {"Glycerol": "#e67e22", "Maltose": "#16a085", "Xylose": "#2980b9"}
    off = {"Glycerol": -.24, "Maltose": 0, "Xylose": .24}
    rng = np.random.default_rng(1)
    for xi, ss in enumerate(ts):
        d = pv[pv["subsystem"] == ss]
        for sn in ["Glycerol", "Maltose", "Xylose"]:
            yv = d[sn + "_d"].values * 100; xj = xi + off[sn] + (rng.random(len(yv)) - .5) * .13
            ax.scatter(xj, yv, s=50, c=cm[sn], alpha=.85, edgecolors="white", lw=.5)
    ax.axhline(0, color="#888", lw=1)
    ax.set_xticks(range(len(ts)))
    ax.set_xticklabels([_wrap2(s) for s in ts], rotation=L["D_rotation"], ha="right",
                       rotation_mode="anchor", fontsize=L["D_fontsize"])
    ax.set_ylabel("Δ subsystem flux share\nvs glucose (percentage points)")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0], [0], marker="o", ls="", color=cm[s], label=s) for s in cm],
              fontsize=8.5, loc="upper right", title="substrate")


# ── heatmap data ────────────────────────────────────────────────────────────────
_sg = fva_long[fva_long["carbon_source"] == "Glucose"]
_Rfull = _sg.pivot(index="reaction", columns="strain_short", values="range").reindex(columns=STRAINS_ORDERED).fillna(0.0)
_R = _Rfull[_Rfull.max(axis=1) > 1e-6]
_subof = pd.Series({r: (rxn_info["subsystem"].get(r, "Unknown") if r in rxn_info.index else "Unknown") for r in _R.index}).fillna("Unknown")
_ord = pd.DataFrame({"s": _subof, "m": _R.mean(axis=1)}).sort_values(["s", "m"], ascending=[True, False]).index
_A = _R.reindex(_ord); _Asub = _subof.reindex(_ord)
_TICKS = [0, 10, 100, 1000, 2000]


def _draw_A(ax):
    im = ax.imshow(np.log10(_A.values + 1), cmap="viridis", aspect="auto", interpolation="nearest",
                   vmin=0, vmax=np.log10(2001))
    ax.set_xticks(range(len(STRAINS_ORDERED))); ax.set_xticklabels(STRAINS_ORDERED, rotation=40, ha="right", fontsize=7)
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    tick = L["HM_tick"]
    b, bounds = 0, [0]
    for s, grp in _Asub.groupby(_Asub, sort=False):
        n = len(grp)
        if n >= L["HM_label_min_rows"]:
            ax.text(-0.5 - tick - 0.15, b + n / 2, str(s), va="center", ha="right",
                    fontsize=L["HM_label_fontsize"])
        b += n; bounds.append(b)
    # thin black subsystem-boundary ticks, sticking out to the left of the heatmap (no interior lines)
    for bb in bounds:
        ax.plot([-0.5, -0.5 - tick], [bb - 0.5, bb - 0.5],
                color=L["HM_sep_color"], lw=L["HM_sep_lw"], clip_on=False)
    ax.set_xlim(-0.5, len(STRAINS_ORDERED) - 0.5)
    ax.set_ylim(len(_A) - .5, -.5)
    return im


def _cbar(fig, im, ax):
    cb = fig.colorbar(im, ax=ax, shrink=.4, label="flux range (mmol/gDW/h)")
    cb.set_ticks([np.log10(t + 1) for t in _TICKS]); cb.set_ticklabels(_TICKS)


def main():
    print("Building figures from CSVs in data/intermediate/ ...")
    # ── supplementary heatmap ──
    fig, ax = plt.subplots(figsize=L["HM_figsize"]); imv = _draw_A(ax)
    ax.set_title(f"FVA flux range is near-identical across strains for most reactions\n"
                 f"{len(_A)} reactions with non-zero range (rows, by subsystem)", fontsize=10.5)
    _cbar(fig, imv, ax); _save(fig, "fig_fva_allrxn_heatmap")

    # ── four-panel manuscript figure ──
    with plt.rc_context(_GGPLOT_BW):
        fig = plt.figure(figsize=(14, 11))
        outer = fig.add_gridspec(2, 1, hspace=0.42)
        top = outer[0].subgridspec(1, 2, width_ratios=[.52, .48], wspace=.2)
        bot = outer[1].subgridspec(1, 2, width_ratios=[.34, .66], wspace=.22)
        for cell, fn, lab in [(top[0], _draw_growth, "A"), (top[1], _draw_pfba, "B"),
                              (bot[0], _draw_loo, "C"), (bot[1], _draw_consist, "D")]:
            ax = fig.add_subplot(cell); fn(ax)
            ax.set_title(lab, loc="left", fontsize=17, fontweight="bold")
        _save(fig, "fig_mechanistic_manuscript")
    print("done")


if __name__ == "__main__":
    main()
