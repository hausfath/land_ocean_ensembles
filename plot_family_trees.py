"""Render the LSAT and SST family trees as left-to-right dendrograms.

Outputs:
    family_tree_land.png
    family_tree_ocean.png

Trees are encoded as nested Python lists (strings are leaves, lists are
internal branches). Top-level method families are shown as labelled bands;
each family carries equal probability at the root (P = 1/N).
"""
from pathlib import Path

import matplotlib.patches as mp
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent

BLUE = "#1f4e8c"
DARK = "#0e2949"
EDGE = "#5b6673"
NATIVE_FC, NATIVE_EC = "#e9f3ff", BLUE
PSEUDO_FC, PSEUDO_EC = "#fff7dc", "#b8912a"
BAND_FC = "#f2f4f7"


# ---------- tree definitions ---------------------------------------------
# Each top-level branch is one method family (equal P=1/5 at root).
# Within the CRU-lineage family the two leaves share P=1/10 each.
LAND_TREE = [
    "berkeley_earth",          # HOMOG / SCALPEL family
    "noaa_land",               # PAIRWISE PHA family
    ["crutem5", "glosatlat"],  # CRU-lineage family
    "dclsat",                  # DYNAMICAL-CONSTRAINT family
    "c_lsat",                  # CMA-HOMOGENIZATION family
]
LAND_FAMILIES = ["Homogenization / scalpel", "Pairwise homogenization",
                 "CRU lineage", "Dynamical constraint", "CMA homogenization"]

OCEAN_TREE = [
    "hadsst4",
    "ersstv6",
    ["cobe_sst3", "cobe_sst3_donor"],  # COBE family: two SST3 siblings (native + HadSST4-donor)
    "dcent_sst",
]
OCEAN_FAMILIES = ["HadSST", "ERSST", "COBE", "Dynamical constraint"]

LAND_LEAVES = {
    "berkeley_earth": ("Berkeley Earth Hi-res",       "1750–2026", "10 native members"),
    "noaa_land":      ("NOAAGlobalTemp Land v6.1",    "1850–2026", "DCLSAT pseudo-ensemble (m 1–100)"),
    "crutem5":        ("CRUTEM 5.1.0.0",              "1850–2026", "200 σ-synthesized members"),
    "glosatlat":      ("GloSATLAT 1.0.0.0",           "1781–2021", "200 σ-synthesized members"),
    "dclsat":         ("DCLSAT-I (DCENT-I v1.1.0.0)", "1850–2025", "200 native members (infilled)"),
    "c_lsat":         ("C-LSAT 2.1 (CMA)",            "1850–2025", "DCLSAT pseudo-ensemble (m 101–200)"),
}

OCEAN_LEAVES = {
    "hadsst4":         ("HadSST 4.2.0.0",       "1850–2026", "200 native members + σ noise"),
    "ersstv6":         ("ERSSTv6",              "1850–2025", "1000 native members"),
    "cobe_sst3":       ("COBE-SST3 (native)",   "1870–2024", "300 native perturbation members"),
    "cobe_sst3_donor": ("COBE-SST3 (donor)",    "1850–2024", "HadSST4-wrapped pseudo-ensemble"),
    "dcent_sst":       ("DCENT-I SST v1.1.0.0", "1850–2025", "200 native members (infilled)"),
}


# ---------- layout / render ----------------------------------------------
def _iter_leaves(node):
    if isinstance(node, str):
        yield node
    else:
        for c in node:
            yield from _iter_leaves(c)


def render(tree, leaf_info, families, root_lines, out_path):
    leaves = list(_iter_leaves(tree))
    n = len(leaves)
    ROW = 1.0
    leaf_y = {name: i * ROW for i, name in enumerate(leaves)}

    def max_depth(node):
        return 0 if isinstance(node, str) else 1 + max(max_depth(c) for c in node)

    depth = max_depth(tree)

    ROOT_W, ROOT_X0 = 3.7, 0.0
    SPINE_X = ROOT_X0 + ROOT_W + 0.55
    COL_W = 0.75
    LEAF_X = SPINE_X + depth * COL_W + 0.35
    BOX_W, BOX_H = 5.6, 0.74

    node_pos, edges = [], []

    def layout(node, d):
        if isinstance(node, str):
            y = leaf_y[node]
            node_pos.append((LEAF_X, y, "leaf", node))
            return y
        ys = [layout(c, d + 1) for c in node]
        y = 0.5 * (min(ys) + max(ys))
        x = SPINE_X + d * COL_W
        node_pos.append((x, y, "node", None))
        for c, cy in zip(node, ys):
            cx = LEAF_X if isinstance(c, str) else SPINE_X + (d + 1) * COL_W
            edges.append((x, y, cx, cy))
        return y

    root_y = layout(tree, 0)

    fig, ax = plt.subplots(figsize=(11.5, 0.62 * n + 1.6))

    fam_bounds, i = [], 0
    for branch in tree:
        k = len(list(_iter_leaves(branch)))
        fam_bounds.append((i, i + k - 1))
        i += k
    for fi, ((i0, i1), fam) in enumerate(zip(fam_bounds, families)):
        y0, y1 = i0 * ROW - 0.44, i1 * ROW + 0.44
        if fi % 2 == 0:
            ax.add_patch(mp.FancyBboxPatch(
                (SPINE_X - 0.30, y0), LEAF_X + BOX_W + 0.35 - (SPINE_X - 0.30), y1 - y0,
                boxstyle="round,pad=0.0,rounding_size=0.10",
                fc=BAND_FC, ec="none", zorder=0))
        ax.text(LEAF_X + BOX_W + 0.28, 0.5 * (y0 + y1), fam,
                ha="left", va="center", fontsize=11, color=EDGE,
                fontweight="bold", style="italic")

    for x1, y1, x2, y2 in edges:
        ax.plot([x1, x1], [y1, y2], color=EDGE, lw=1.5, solid_capstyle="round", zorder=1)
        ax.plot([x1, x2], [y2, y2], color=EDGE, lw=1.5, solid_capstyle="round", zorder=1)

    for x, y, kind, _ in node_pos:
        if kind == "node":
            ax.add_patch(mp.Circle((x, y), 0.07, color=BLUE, zorder=3))

    for x, y, kind, lbl in node_pos:
        if kind != "leaf":
            continue
        title, span, note = leaf_info[lbl]
        pseudo = "pseudo" in note.lower() or "synthesized" in note.lower()
        ax.add_patch(mp.FancyBboxPatch(
            (LEAF_X, y - BOX_H / 2), BOX_W, BOX_H,
            boxstyle="round,pad=0.0,rounding_size=0.14",
            fc=PSEUDO_FC if pseudo else NATIVE_FC,
            ec=PSEUDO_EC if pseudo else NATIVE_EC, lw=1.3, zorder=2))
        ax.text(LEAF_X + 0.18, y - 0.085, title, ha="left", va="center",
                fontsize=12.5, fontweight="bold", color="#111", zorder=4)
        ax.text(LEAF_X + 0.18, y + 0.21, f"{span}    {note}", ha="left", va="center",
                fontsize=10.5, color="#444", zorder=4)

    RH = 1.35
    ax.add_patch(mp.FancyBboxPatch(
        (ROOT_X0, root_y - RH / 2), ROOT_W, RH,
        boxstyle="round,pad=0.0,rounding_size=0.18",
        fc=BLUE, ec=DARK, lw=1.4, zorder=2))
    cx = ROOT_X0 + ROOT_W / 2
    ax.text(cx, root_y - 0.33, root_lines[0], ha="center", va="center",
            fontsize=13.5, fontweight="bold", color="white", zorder=3)
    ax.text(cx, root_y + 0.06, root_lines[1], ha="center", va="center",
            fontsize=11, color="#dce9fb", zorder=3)
    ax.text(cx, root_y + 0.40, f"P = 1/{len(tree)} per family", ha="center",
            va="center", fontsize=11, color="#dce9fb", zorder=3)
    ax.annotate("", xy=(ROOT_X0 + ROOT_W - 0.02, root_y), xytext=(SPINE_X, root_y),
                arrowprops=dict(arrowstyle="<|-", color=BLUE, lw=2.4,
                                mutation_scale=24), zorder=1)

    ax.legend(handles=[
        mp.Patch(facecolor=NATIVE_FC, edgecolor=NATIVE_EC, label="Native ensemble"),
        mp.Patch(facecolor=PSEUDO_FC, edgecolor=PSEUDO_EC,
                 label="Pseudo-ensemble (donor or σ-synthesized members)"),
    ], loc="upper left", bbox_to_anchor=(0.0, -0.015), frameon=False,
        fontsize=11.5, ncol=2, handlelength=1.6, bbox_transform=ax.transAxes)

    ax.set_xlim(-0.3, LEAF_X + BOX_W + 3.4)
    ax.set_ylim(-0.85, (n - 1) * ROW + 0.95)
    ax.invert_yaxis()
    ax.set_axis_off()

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_path}  ({n} leaves, max depth {depth})")


def main():
    render(LAND_TREE, LAND_LEAVES, LAND_FAMILIES,
           ("Land LSAT ensemble", "10,000 members · 1850–2025"),
           ROOT / "family_tree_land.png")
    render(OCEAN_TREE, OCEAN_LEAVES, OCEAN_FAMILIES,
           ("Ocean SST ensemble", "10,000 members · 1850–2025"),
           ROOT / "family_tree_ocean.png")


if __name__ == "__main__":
    main()
