"""Render the LSAT and SST family trees as left-to-right dendrograms.

Outputs:
    family_tree_land.png
    family_tree_ocean.png

Trees are encoded as nested Python lists (strings are leaves, lists are
internal branches). Probabilities are not drawn explicitly — the tree
topology + the methodology MDs are sufficient documentation.
"""
from pathlib import Path

import matplotlib.patches as mp
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent


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

OCEAN_TREE = [
    "hadsst4",
    "ersstv6",
    "cobe_sst2",
    "dcent_sst",
]


LAND_LEAVES = {
    "berkeley_earth": ("Berkeley Earth Highres",  "1750–2026", "10 native members"),
    "noaa_land":      ("NOAAGlobalTemp Land v6.1","1850–2026", "+ DCLSAT (m 1-100) pseudo-ensemble"),
    "crutem5":        ("CRUTEM 5.1.0.0",          "1850–2026", "+ synth 200 from σ components"),
    "glosatlat":      ("GloSATLAT 1.0.0.0",       "1781–2021", "+ synth 200 from σ components"),
    "dclsat":         ("DCLSAT-I (DCENT-I v1.1.0.0)","1850–2025","200 native members (infilled)"),
    "c_lsat":         ("C-LSAT 2.1 (CMA)",        "1850–2025", "+ DCLSAT (m 101-200) pseudo-ensemble"),
}

OCEAN_LEAVES = {
    "hadsst4":   ("HadSST 4.2.0.0",            "1850–2026", "200 native bias members + σ noise"),
    "ersstv6":   ("ERSSTv6",                   "1850–2024", "1000 native members + 2025 frozen-offset"),
    "cobe_sst2": ("COBE-SST 2",                "1850–2026", "+ HadSST4 pseudo-ensemble"),
    "dcent_sst": ("DCENT-I SST (v1.1.0.0)",    "1850–2025", "200 native members (infilled)"),
}


# ---------- layout helpers -----------------------------------------------
def walk_leaves(node, acc):
    if isinstance(node, str):
        acc.append(node); return
    for c in node:
        walk_leaves(c, acc)


def render(tree, leaf_info, root_title, root_subtitle, fig_title, out_path,
           figsize):
    leaves = []
    walk_leaves(tree, leaves)

    ROW_GAP = 1.6
    leaf_y = {name: i * ROW_GAP for i, name in enumerate(leaves)}

    node_positions = []
    edges = []

    def layout(node, depth):
        if isinstance(node, str):
            y = leaf_y[node]
            node_positions.append((depth, y, "leaf", node))
            return y
        child_ys = [layout(c, depth + 1) for c in node]
        y = sum(child_ys) / len(child_ys)
        node_positions.append((depth, y, "node", None))
        for c, cy in zip(node, child_ys):
            edges.append((depth, y, depth + 1, cy, isinstance(c, str)))
        return y

    root_y = layout(tree, 0)
    max_depth = max(d for d, *_ in node_positions)
    n_leaves = len(leaves)

    LABEL_X = max_depth + 0.4
    ROOT_OUT_X = -3.5

    fig, ax = plt.subplots(figsize=figsize)

    for x1, y1, x2, y2, is_leaf_edge in edges:
        ax.plot([x1, x1], [y1, y2], color="#444", lw=1.1)
        end_x = LABEL_X if is_leaf_edge else x2
        ax.plot([x1, end_x], [y2, y2], color="#444", lw=1.1)

    for d, y, kind, _ in node_positions:
        if kind == "node":
            ax.add_patch(mp.Circle((d, y), 0.085, color="#1f4e8c", zorder=3))

    for d, y, kind, lbl in node_positions:
        if kind != "leaf":
            continue
        title, span, note = leaf_info.get(lbl, ("?", "?", ""))
        is_pseudo = "pseudo" in note.lower() or "synth" in note.lower()
        box_color = "#fffbe6" if is_pseudo else "#e9f3ff"
        edge_color = "#cfa01a" if is_pseudo else "#1f4e8c"
        ax.annotate(title, xy=(LABEL_X, y), xytext=(LABEL_X + 0.1, y - 0.12),
                    ha="left", va="bottom", fontsize=10.5, fontweight="bold",
                    family="sans-serif", zorder=4)
        ax.annotate(f"{span}   {note}",
                    xy=(LABEL_X, y), xytext=(LABEL_X + 0.1, y + 0.08),
                    ha="left", va="top", fontsize=9.0, color="#333",
                    family="sans-serif", zorder=4)
        ax.add_patch(mp.FancyBboxPatch(
            (LABEL_X, y - 0.55), 6.4, 1.10,
            boxstyle="round,pad=0.0,rounding_size=0.18",
            fc=box_color, ec=edge_color, lw=1.0, zorder=1,
        ))

    out_w, out_h = 4.5, 1.5
    out_right = ROOT_OUT_X + out_w / 2 + 0.05
    out_left = ROOT_OUT_X - out_w / 2 + 0.05

    ax.annotate("", xy=(out_right, root_y),
                xytext=(0, root_y),
                arrowprops=dict(arrowstyle="-|>", color="#1f4e8c",
                                lw=2.6, mutation_scale=22),
                zorder=2)
    ax.add_patch(mp.FancyBboxPatch(
        (out_left, root_y - out_h / 2), out_w, out_h,
        boxstyle="round,pad=0.0,rounding_size=0.22",
        fc="#1f4e8c", ec="#0e2949", lw=1.2, zorder=2,
    ))
    cx = (out_left + out_right) / 2
    ax.text(cx, root_y - 0.40, root_title,
            ha="center", va="center",
            fontsize=12, fontweight="bold", color="white", zorder=3)
    ax.text(cx, root_y + 0.10, "10,000 members",
            ha="center", va="center",
            fontsize=9.5, color="#e9f3ff", zorder=3)
    ax.text(cx, root_y + 0.40, root_subtitle,
            ha="center", va="center",
            fontsize=9.5, color="#e9f3ff", zorder=3)

    ax.text((out_left + LABEL_X + 6) / 2, -1.0, fig_title,
            ha="center", va="center", fontsize=12, color="#1f4e8c",
            fontweight="bold")

    ax.set_xlim(ROOT_OUT_X - out_w / 2 - 0.4, LABEL_X + 7.1)
    ax.set_ylim(-2.0, (n_leaves - 1) * ROW_GAP + 1.2)
    ax.invert_yaxis()
    ax.set_axis_off()

    legend_elems = [
        mp.Patch(facecolor="#e9f3ff", edgecolor="#1f4e8c",
                 label="Native-ensemble dataset"),
        mp.Patch(facecolor="#fffbe6", edgecolor="#cfa01a",
                 label="Pseudo-ensemble (regular dataset + donor or σ-synthesized perturbations)"),
    ]
    ax.legend(handles=legend_elems, loc="lower center",
              bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=10, ncol=2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"wrote {out_path}")
    print(f"  {n_leaves} leaves, max depth {max_depth}")
    plt.close(fig)


def main():
    render(LAND_TREE, LAND_LEAVES,
           root_title="Land LSAT ensemble",
           root_subtitle="1850–2025  (family-tree weighted)",
           fig_title=("Land surface air temperature family tree\n"
                      "(6 LSAT products, 5 equal-weight method families)"),
           out_path=ROOT / "family_tree_land.png",
           figsize=(13, 8))

    render(OCEAN_TREE, OCEAN_LEAVES,
           root_title="Ocean SST ensemble",
           root_subtitle="1850–2025  (family-tree weighted)",
           fig_title=("Sea-surface temperature family tree\n"
                      "(4 SST products, equal-weight method families)"),
           out_path=ROOT / "family_tree_ocean.png",
           figsize=(13, 6))


if __name__ == "__main__":
    main()
