#!/usr/bin/env python3
"""One page with all candidate cells as small top views, grouped by class, each clickable to a full-size image.

Each thumbnail is the reconstructed periodic cell (as cut, before the DFT z-shift) seen from above, atoms as discs
coloured by layer (same palette as the lineage figure), the cell outline in blue; the enlarged view is the same
drawing at full resolution with a caption (state id, class/split, size, atoms, requested U, repair, CN stratum).

Usage (from Au_Cl/):  rough_sampling_v1/pyrun_rs.sh rough_sampling_v1/render_candidates_page.py [--manifest rough200/rough200_manifest.json] [--out figures]
Output: rough_sampling_v1/<out>/candidates/<state_id>.png and <out>/candidates.html
"""
import argparse
import html
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from ase.io import read  # noqa: E402
from matplotlib.patches import Circle, Polygon  # noqa: E402

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
R = f"{ROOT}/rough_sampling_v1"
sys.path.insert(0, R)
from extract_cells import layers_from_z  # noqa: E402
from render_rough_lineage import LAYER_COL, RAD  # noqa: E402

CLASS_NAME = {"A": "A — roughened strips", "B": "B — islands", "C": "C — pits", "D": "D — strips with atoms moved to the step foot"}


def draw(at, path, title):
    c = at.get_cell().array; P = at.get_positions(); lay = layers_from_z(at)
    fig, ax = plt.subplots(figsize=(6.4, 6.0), dpi=150)
    for i in np.argsort(P[:, 2]):
        ax.add_patch(Circle(P[i, :2], RAD, fc=LAYER_COL.get(int(lay[i]), "#000"), ec="k", lw=0.35))
    # a faded ring of periodic images so the seams can be judged
    for si in (-1, 0, 1):
        for sj in (-1, 0, 1):
            if (si, sj) == (0, 0): continue
            sh = si * c[0][:2] + sj * c[1][:2]
            for i in np.argsort(P[:, 2]):
                q = P[i, :2] + sh
                if -4 < q[0] - 0 < c[0][0] + c[1][0] + 4 and -4 < q[1] < c[1][1] + 4:
                    ax.add_patch(Circle(q, RAD, fc=LAYER_COL.get(int(lay[i]), "#000"), ec="none", alpha=0.25))
    poly = np.array([[0, 0], c[0][:2], c[0][:2] + c[1][:2], c[1][:2]])
    ax.add_patch(Polygon(poly, closed=True, fill=False, ec="#0033cc", lw=2.5))
    ax.set_xlim(-4, c[0][0] + c[1][0] + 4); ax.set_ylim(-4, c[1][1] + 4); ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_title(title, fontsize=11)
    fig.tight_layout(pad=0.3); fig.savefig(path); plt.close(fig)
    # quantise to a 64-colour palette: the alpha-blended image ring otherwise makes ~350 KB per file; this gives ~60 KB
    from PIL import Image
    im = Image.open(path).convert("RGB").quantize(colors=64, method=Image.Quantize.MEDIANCUT); im.save(path, optimize=True)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--manifest", default="rough200/rough200_manifest.json"); ap.add_argument("--out", default="figures")
    a = ap.parse_args()
    man = json.load(open(f"{R}/{a.manifest}")); S = man["states"]
    outdir = f"{R}/{a.out}/candidates"; os.makedirs(outdir, exist_ok=True)
    for k, s in enumerate(S):
        png = f"{outdir}/{s['state_id']}.png"
        if not os.path.exists(png):
            at = read(f"{s['cell_dir']}/cell.extxyz")
            draw(at, png, f"{s['cell_id']}\n{s['cell']}, {s['n_atoms']} Au, U = {s['U_V']:+.2f} V")
        if (k + 1) % 25 == 0: print(f"  {k+1}/{len(S)}", flush=True)
    # page
    status = "FROZEN" if man.get("frozen") else "CANDIDATE (not frozen; nothing submitted)"
    cards = []
    for cls in "ABCD":
        rows = [s for s in S if s["cls"] == cls]
        cards.append(f'<h2>{html.escape(CLASS_NAME[cls])} <span class="n">{len(rows)} states</span></h2><div class="grid">')
        for s in sorted(rows, key=lambda s: (s["split"], s["cell"], s["state_id"])):
            rep = "no repair" if s["repair"].startswith("no repair") else "seam band minimised"
            cap = (f"{s['state_id']} | class {s['cls']} / {s['split']} | {s['cell']}, {s['n_atoms']} Au, k {s['kpoints']} | U = {s['U_V']:+.2f} V "
                   f"(TARGETMU {s['TARGETMU_eV']:.4f} eV) | centre CN {s['cn']} ({s['cn_stratum']}) | {rep} | parent {s['parent_id']} @ {s.get('time_ps') or 0:.0f} ps"
                   + (f" | {s['pair']}" if s.get("pair") else ""))
            cards.append(f'<figure><a href="candidates/{html.escape(s["state_id"])}.png" data-cap="{html.escape(cap)}">'
                         f'<img loading="lazy" src="candidates/{html.escape(s["state_id"])}.png" alt="{html.escape(s["state_id"])}"></a>'
                         f'<figcaption>{html.escape(s["cell"])} · {s["n_atoms"]} Au · U {s["U_V"]:+.2f} V · {html.escape(s["split"])}</figcaption></figure>')
        cards.append("</div>")
    page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Rough candidates ({len(S)})</title>
<style>
:root{{--bg:#fff;--fg:#111;--mut:#666;--line:#ddd}} @media(prefers-color-scheme:dark){{:root:not([data-theme=light]){{--bg:#141414;--fg:#eee;--mut:#aaa;--line:#333}}}}
:root[data-theme=dark]{{--bg:#141414;--fg:#eee;--mut:#aaa;--line:#333}}
body{{margin:0;padding:16px;background:var(--bg);color:var(--fg);font:14px/1.45 "DejaVu Sans",system-ui,sans-serif}}
h1{{font-size:20px;margin:0 0 4px}} h2{{font-size:16px;margin:22px 0 8px;border-bottom:1px solid var(--line);padding-bottom:4px}} .n{{color:var(--mut);font-weight:normal}}
p{{color:var(--mut);margin:0 0 6px;max-width:1100px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:10px}}
figure{{margin:0}} figure img{{width:100%;height:auto;display:block;border:1px solid var(--line);border-radius:4px;cursor:zoom-in;background:#fff}}
figcaption{{font-size:11px;color:var(--mut);margin-top:2px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}}
#lb{{display:none;position:fixed;inset:0;background:rgba(0,0,0,.88);z-index:9;flex-direction:column;align-items:center;justify-content:center;padding:16px;cursor:zoom-out}}
#lb img{{max-width:96vw;max-height:84vh;background:#fff;border-radius:6px}} #lb div{{color:#eee;font-size:13px;margin-top:10px;text-align:center;max-width:96vw}}
.leg span{{display:inline-block;width:12px;height:12px;border-radius:50%;border:1px solid #333;vertical-align:middle;margin:0 4px 0 10px}}
</style></head><body>
<h1>Rough-sampling candidate cells — {len(S)} states, status: {html.escape(status)}</h1>
<p>Top views of the reconstructed periodic cells (as cut from the MD frames, before the DFT z-shift); faded discs are the periodic images, the blue outline is the cell.
Click a thumbnail to enlarge. Grouped by parent class; within a class ordered by split (test, train, val), size, id.
Source: <code>rough200/rough200_manifest.json</code> (seed {man['seed']}); per-state table: <code>rough200/rough200_cells.csv</code>.</p>
<p class="leg">Layer colours: <span style="background:{LAYER_COL[0]}"></span>0 (fixed) <span style="background:{LAYER_COL[1]}"></span>1 (fixed)
<span style="background:{LAYER_COL[2]}"></span>2 (pit floor) <span style="background:{LAYER_COL[3]}"></span>3 (terrace) <span style="background:{LAYER_COL[4]}"></span>4 (adatom level)</p>
{''.join(cards)}
<div id="lb" onclick="this.style.display='none'"><img id="lbimg" alt=""><div id="lbcap"></div></div>
<script>
document.querySelectorAll('.grid a').forEach(a=>a.addEventListener('click',e=>{{e.preventDefault();const lb=document.getElementById('lb');
document.getElementById('lbimg').src=a.getAttribute('href');document.getElementById('lbcap').textContent=a.dataset.cap;lb.style.display='flex';}}));
document.addEventListener('keydown',e=>{{if(e.key==='Escape')document.getElementById('lb').style.display='none';}});
</script></body></html>"""
    open(f"{R}/{a.out}/candidates.html", "w").write(page)
    print(f"wrote {R}/{a.out}/candidates.html with {len(S)} thumbnails")


if __name__ == "__main__":
    main()
