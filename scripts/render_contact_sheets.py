#!/usr/bin/env python3
"""Contact sheets of the schematic panel of every structure, grouped by family.

A reviewer asked for the gallery as a package rather than 33 separate files. Pasting 33 full figures into one
image makes the text unreadable, so each sheet carries only the LEFT panel of each figure -- the plain outline
schematic and its section -- which is the part under review, cropped from the rendered PNG at full resolution and
laid out with the structure name and its one-line description. The full figures stay available individually.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/render_contact_sheets.py
Outputs: analysis/gallery/_sheet_<family>.png  and  analysis/gallery/_sheet_index.md
"""
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "/anvil/scratch/x-rywang/Au_Cl/scripts")
from render_gallery import CJK, FAMILY_ZH  # noqa: E402

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
GAL = f"{ROOT}/analysis/gallery"
ORDER = ["flat Au(111)", "point defect", "reconstruction-related", "strip step", "vicinal step face",
         "kink / edge rearrangement", "single-layer island", "single-layer pit", "composite"]
LEFT_FRAC = 0.31          # the schematic column occupies the leftmost ~31% of each rendered figure


def main():
    gal = json.load(open(f"{GAL}/gallery.json"))
    by = {}
    for sid, m in gal.items():
        by.setdefault(m["family"], []).append(sid)
    lines = ["# 结构简笔示意 · 分族拼图", "",
             "每张拼图取各结构图的最左一栏（简笔俯视轮廓 + 侧面剖线），原分辨率裁切，未缩放文字。",
             "完整图（含按配位数着色的俯视图与侧视图）仍是每个结构一个 PNG。","注意：hcp 配准标记（蓝色 ×）、加上去的 Au（绿环）与移走的位点（紫色虚圆）画在完整图的俯视图上，""拼图这一栏看不到——例如 A1-fcc 与 A1-hcp 的简笔轮廓本来就相同，区别只在配准。", ""]
    made = 0
    for fam in ORDER:
        sids = sorted(by.get(fam, []))
        if not sids: continue
        imgs = []
        for sid in sids:
            p = f"{GAL}/{sid}.png"
            if not os.path.exists(p): continue
            a = plt.imread(p)
            # also stop above the legend strip, which spans the full width and bleeds into the crop
            imgs.append((sid, a[:int(0.88 * a.shape[0]), :int(LEFT_FRAC * a.shape[1])]))
        if not imgs: continue
        ncol = 2 if len(imgs) > 1 else 1
        nrow = int(np.ceil(len(imgs) / ncol))
        hmax = max(im.shape[0] for _, im in imgs); wmax = max(im.shape[1] for _, im in imgs)
        fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * wmax / 190, nrow * (hmax / 190 + 0.42)), dpi=190)
        axes = np.atleast_1d(axes).ravel()
        for ax in axes: ax.axis("off")
        for ax, (sid, im) in zip(axes, imgs):
            ax.imshow(im); ax.axis("off")
            # the CJK fallback face has no U+00C5; the sheet caption is plain text, so spell the unit
            note = gal[sid].get("schematic_note", "").replace("\u00c5", " 埃")
            ax.set_title(f"{sid}", fontsize=11, color="#14181c", pad=2)
            ax.text(0.5, -0.02, note, transform=ax.transAxes, ha="center", va="top",
                    fontproperties=CJK, fontsize=7.6, color="#5a616a", wrap=True)
        fig.suptitle(f"{FAMILY_ZH.get(fam, fam)}　·　{len(imgs)} 个结构", fontproperties=CJK,
                     fontsize=15, y=0.995, color="#14181c")
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        # full-family slug: "single-layer island" and "single-layer pit" collided on the first word and
        # overwrote each other's sheet
        slug = fam.replace("/", "-").replace("(", "").replace(")", "").replace(" ", "-").lower()
        out = f"{GAL}/_sheet_{slug}.png"
        fig.savefig(out, facecolor="white"); plt.close(fig)
        made += 1
        lines.append(f"- **{FAMILY_ZH.get(fam, fam)}**（{len(imgs)}）：`{os.path.basename(out)}` — "
                     + "、".join(sids))
    open(f"{GAL}/_sheet_index.md", "w").write("\n".join(lines) + "\n")
    print(f"{made} contact sheets -> {GAL}/_sheet_*.png")


if __name__ == "__main__":
    main()
