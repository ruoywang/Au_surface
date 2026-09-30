#!/usr/bin/env python3
"""Build the standalone HTML report from analysis/web/data.json (inlined, because the artifact CSP blocks fetch).

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_web_report.py
Output: analysis/web/index.html  (publish with the gallery and map PNGs as supporting files)
"""
import json
import os

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
WEB = f"{ROOT}/analysis/web"

HTML = r"""<title>Au(111) Charging Atlas</title>
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Spectral:ital,wght@0,400;0,500;0,600;1,400&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{
  --bg:#f6f4ef; --surface:#fffefb; --surface-2:#f0ece4; --ink:#15181c; --ink-2:#565d65; --muted:#8a9199;
  --line:#e0dbd1; --line-2:#cdc6b9;
  --gold:#96650d; --gold-soft:#e8d9b4; --teal:#1d6a85; --teal-soft:#cfe3ea; --red:#a8382a; --amber:#c4841f;
  --terrace:#b99a3e; --sub:#79828c;
  --shadow:0 1px 2px rgba(30,24,10,.05),0 6px 20px -12px rgba(30,24,10,.18);
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
  --sans:"IBM Plex Sans",system-ui,-apple-system,Segoe UI,sans-serif;
  --serif:"Spectral",Georgia,"Times New Roman",serif;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --bg:#121417; --surface:#191c21; --surface-2:#21252b; --ink:#e9ebee; --ink-2:#a6adb5; --muted:#79818a;
  --line:#2a2e35; --line-2:#3a3f47;
  --gold:#d9a743; --gold-soft:#4a3c1c; --teal:#5fb2cd; --teal-soft:#1d3b46; --red:#d9705f; --amber:#dda23e;
  --terrace:#c9ab52; --sub:#8b949e;
  --shadow:0 1px 2px rgba(0,0,0,.3),0 8px 24px -14px rgba(0,0,0,.7);
}}
:root[data-theme="dark"]{
  --bg:#121417; --surface:#191c21; --surface-2:#21252b; --ink:#e9ebee; --ink-2:#a6adb5; --muted:#79818a;
  --line:#2a2e35; --line-2:#3a3f47;
  --gold:#d9a743; --gold-soft:#4a3c1c; --teal:#5fb2cd; --teal-soft:#1d3b46; --red:#d9705f; --amber:#dda23e;
  --terrace:#c9ab52; --sub:#8b949e;
  --shadow:0 1px 2px rgba(0,0,0,.3),0 8px 24px -14px rgba(0,0,0,.7);
}
*{box-sizing:border-box}
body{background:var(--bg);color:var(--ink);font-family:var(--sans);font-size:15px;line-height:1.62;margin:0;
  -webkit-font-smoothing:antialiased}
.wrap{max-width:1140px;margin:0 auto;padding-inline:20px;padding-block:0}
h1,h2,h3{font-family:var(--serif);font-weight:500;text-wrap:balance;margin:0}
h1{font-size:clamp(30px,5.4vw,50px);line-height:1.1;letter-spacing:-.015em}
h2{font-size:clamp(21px,3vw,28px);line-height:1.2;margin-bottom:.35em}
h3{font-size:17px;line-height:1.3}
p{margin:0 0 .9em}
a{color:var(--teal)}
.lede{font-size:clamp(16px,2.1vw,18.5px);color:var(--ink-2);max-width:66ch}
.eyebrow{font-family:var(--mono);font-size:11px;letter-spacing:.13em;text-transform:uppercase;color:var(--gold);
  margin-bottom:.5em}
.mono{font-family:var(--mono);font-variant-numeric:tabular-nums}
small,.small{font-size:12.5px;color:var(--ink-2);line-height:1.5}
.caption{font-size:12.5px;color:var(--muted);line-height:1.5;margin-top:.5em}

header.top{border-bottom:1px solid var(--line);background:
  radial-gradient(1200px 340px at 12% -40%,var(--gold-soft),transparent 65%),var(--surface)}
header.top .wrap{padding-block:46px 34px}
.kicker{display:flex;flex-wrap:wrap;gap:8px;margin-top:20px}
.chip{font-family:var(--mono);font-size:11.5px;border:1px solid var(--line-2);border-radius:999px;
  padding:3px 11px;color:var(--ink-2);background:var(--surface)}

nav.sticky{position:sticky;top:env(safe-area-inset-top,0px);z-index:40;background:color-mix(in srgb,var(--bg) 92%,transparent);
  backdrop-filter:blur(10px);border-bottom:1px solid var(--line)}
nav.sticky .wrap{display:flex;gap:2px;overflow-x:auto;padding-block:7px;scrollbar-width:none}
nav.sticky .wrap::-webkit-scrollbar{display:none}
nav.sticky a{font-family:var(--mono);font-size:11.5px;color:var(--ink-2);text-decoration:none;white-space:nowrap;
  padding:5px 11px;border-radius:6px}
nav.sticky a:hover{background:var(--surface-2);color:var(--ink)}

section{padding-block:52px;border-bottom:1px solid var(--line)}
section:last-of-type{border-bottom:0}
.finding{display:grid;grid-template-columns:minmax(0,1fr);gap:22px}
.claim{font-family:var(--serif);font-size:clamp(18px,2.4vw,22px);line-height:1.38;color:var(--ink);
  border-left:3px solid var(--gold);padding-left:16px;max-width:60ch;margin:4px 0 2px}
.card{background:var(--surface);border:1px solid var(--line);border-radius:10px;box-shadow:var(--shadow)}
.pad{padding:18px}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:18px}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:1px;background:var(--line);
  border:1px solid var(--line);border-radius:10px;overflow:hidden}
.stat{background:var(--surface);padding:14px 16px}
.stat .v{font-family:var(--mono);font-size:22px;font-weight:500;letter-spacing:-.02em;color:var(--ink)}
.stat .l{font-size:11.5px;color:var(--ink-2);margin-top:2px;line-height:1.35}
.tablewrap{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--surface)}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{text-align:left;padding:8px 12px;border-bottom:1px solid var(--line);white-space:nowrap}
th{font-family:var(--mono);font-size:10.5px;letter-spacing:.07em;text-transform:uppercase;color:var(--muted);
  font-weight:500;position:sticky;top:0;background:var(--surface)}
td.num{font-family:var(--mono);font-variant-numeric:tabular-nums;text-align:right}
tr:last-child td{border-bottom:0}
tbody tr:hover{background:var(--surface-2)}

.legend{display:flex;flex-wrap:wrap;gap:10px 16px;font-size:12px;color:var(--ink-2);margin-top:10px}
.legend span{display:inline-flex;align-items:center;gap:6px}
.sw{width:11px;height:11px;border-radius:3px;display:inline-block}
.filters{display:flex;flex-wrap:wrap;gap:6px;margin:16px 0 20px}
.filters button{font-family:var(--mono);font-size:11.5px;border:1px solid var(--line-2);background:var(--surface);
  color:var(--ink-2);border-radius:999px;padding:4px 12px;cursor:pointer}
.filters button[aria-pressed="true"]{background:var(--ink);color:var(--bg);border-color:var(--ink)}
.filters button:focus-visible{outline:2px solid var(--teal);outline-offset:2px}

.gal{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:20px}
.gcard{background:var(--surface);border:1px solid var(--line);border-radius:10px;overflow:hidden;box-shadow:var(--shadow)}
.gcard img{display:block;width:100%;height:auto;background:#fff}
.ghead{display:flex;justify-content:space-between;align-items:baseline;gap:10px;padding:12px 15px 0}
.ghead h3{font-size:16px}
.gmeta{display:flex;flex-wrap:wrap;gap:5px;padding:9px 15px 14px}
.tag{font-family:var(--mono);font-size:10.5px;padding:2px 7px;border-radius:5px;background:var(--surface-2);
  color:var(--ink-2);border:1px solid var(--line)}
.tag.k{background:color-mix(in srgb,var(--red) 16%,transparent);border-color:transparent;color:var(--red)}
.tag.e{background:color-mix(in srgb,var(--amber) 18%,transparent);border-color:transparent;color:var(--amber)}
svg{display:block;max-width:100%;height:auto}
svg text{font-family:var(--mono);fill:var(--ink-2)}
.chartbox{overflow-x:auto}
.note{background:var(--surface-2);border:1px solid var(--line);border-left:3px solid var(--teal);border-radius:8px;
  padding:13px 16px;font-size:13px;color:var(--ink-2);line-height:1.55}
.note b{color:var(--ink);font-weight:600}
footer{padding-block:38px 56px;color:var(--muted);font-size:12.5px}
@media (max-width:560px){.gal{grid-template-columns:1fr}section{padding-block:38px}}
</style>

<header class="top"><div class="wrap">
  <div class="eyebrow">Constant-potential DFT · VASP + VASPsol++ · 1 M implicit electrolyte</div>
  <h1>Au(111) Charging Atlas</h1>
  <p class="lede" style="margin-top:14px">How surface morphology changes the metal's ability to hold charge, and how
  much of that atomic-scale difference actually reaches the ions in the liquid. Thirty-three Au(111) structures,
  107 fixed geometries, 291 converged constant-potential states at three electron chemical potentials.</p>
  <div class="kicker" id="kicker"></div>
</div></header>

<nav class="sticky"><div class="wrap">
  <a href="#pzc">1 · PZC, not capacitance</a>
  <a href="#which">2 · Which morphologies</a>
  <a href="#where">3 · Where the anions go</a>
  <a href="#filter">4 · The electrolyte filter</a>
  <a href="#omega">5 · Potential and stability</a>
  <a href="#gallery">Structure gallery</a>
  <a href="#method">Method</a>
</div></nav>

<div class="wrap">

<section id="pzc"><div class="finding">
  <div><div class="eyebrow">Finding 1</div><h2>Morphology moves the potential of zero charge; it barely moves the capacitance</h2></div>
  <p class="claim">Across all 104 geometries whose sampled points bracket &sigma;&nbsp;=&nbsp;0, the zero-charge
  potential spans 254&nbsp;mV while the secant capacitance spans only 15%. At a fixed potential the charge difference
  between two morphologies comes from the PZC offset about eight times over.</p>
  <div class="stats" id="s1"></div>
  <div class="card pad"><div class="chartbox"><div id="c_sigma"></div></div>
    <div class="caption">Every geometry's surface charge density against the internal potential
    U&nbsp;=&nbsp;&mu;<sub>0</sub>&nbsp;&minus;&nbsp;&mu;<sub>e</sub>, from the actual converged &mu;<sub>e</sub> and
    N<sub>e</sub>. The curves are nearly parallel: the family-to-family difference is a horizontal shift, not a change
    of slope.</div>
    <div class="legend" id="leg1"></div></div>
  <div class="note"><b>Why the split matters.</b> A defect that carries more positive charge at a given potential is
  not necessarily more polarisable. Writing &sigma;(U)&nbsp;=&nbsp;C&nbsp;(U&nbsp;&minus;&nbsp;U<sub>pzc</sub>),
  the spread of &sigma; across geometries divides into &lang;C&rang;&thinsp;&times;&thinsp;spread(U<sub>pzc</sub>) and
  spread(C)&thinsp;&times;&thinsp;|U&nbsp;&minus;&nbsp;&lang;U<sub>pzc</sub>&rang;|. The first term dominates
  everywhere in the sampled window.</div>
  <div class="tablewrap"><table id="t_decomp"></table></div>
  <div class="caption">Cross-cell comparisons carry a numerical offset: the audit records two <em>flat</em> cells
  differing by 72&nbsp;meV in neutral &mu;<sub>e</sub> from cell shape, k-mesh and PREC alone. Inside one cell that
  offset is common to every row, so the within-cell blocks are morphology only.</div>
</div></section>

<section id="which"><div class="finding">
  <div><div class="eyebrow">Finding 2</div><h2>Protrusions lower the zero-charge potential, vacancies raise it</h2></div>
  <p class="claim">Measured against the flat terrace in the <em>same</em> cell, adatoms and step edges shift the PZC
  down by 29&ndash;97&nbsp;mV, while vacancies shift it up by 8&ndash;13&nbsp;mV. The stacking-fault surfaces barely
  move at all. The sign follows the classic work-function argument: a protrusion smooths the electron spill-out and
  lowers the work function; a depression does the opposite.</p>
  <div class="card pad"><div class="chartbox"><div id="c_dpzc"></div></div>
    <div class="caption">Ideal configurations only, each against the flat member of its own cell. Bars show the PZC
    shift; the dot marks the capacitance change over the same comparison.</div></div>
  <div class="note"><b>Capacitance follows, weakly.</b> Added atoms raise the capacitance more than removed ones:
  three adatoms give +15%, one adatom +8%, a step +8%, while a single vacancy gives +1% and the reconstructions are
  flat to within 1%. Roughness that protrudes into the electrolyte both lowers the PZC and increases the capacitance,
  but the PZC effect is the one that changes the charge at a given potential.</div>
</div></section>

<section id="where"><div class="finding">
  <div><div class="eyebrow">Finding 3</div><h2>The anion hot spot is not above the defect</h2></div>
  <p class="claim">Resolving the anion excess by the coordination number of the nearest surface atom, the
  under-coordinated sites hold <em>less</em> anion enrichment than the surrounding terrace in 36 of 41 geometries
  &mdash; by 1.7% at the median and up to 6.1%. Every exception is a vacancy rim, which is a depression rather than
  a protrusion. Both integration conventions give the same ordering.</p>
  <div class="card pad"><div class="chartbox"><div id="c_regions"></div></div>
    <div class="caption">Enrichment K<sub>&Omega;</sub> = &int;n<sub>&minus;</sub> / (n<sub>b</sub>&int;S<sub>ion</sub>)
    in the boundary-relative window, at the most positive sampled potential of each structure. Bars are grouped by
    coordination class; the dashed line is that structure's cell average.</div>
    <div class="legend" id="leg3"></div></div>
  <div class="note"><b>This is the charging result seen from the liquid side.</b> An adatom lowers the local
  zero-charge potential, so at a fixed electrode potential that patch is relatively <em>less</em> positive than the
  terrace around it and attracts fewer anions. A vacancy raises it, and the vacancy rims are exactly the geometries
  that break the trend. Findings 2 and 3 are one piece of physics measured from the two sides of the interface.
  <br><br>K and &Gamma; answer different questions and are both reported: a tight pocket can reach a high
  concentration while holding very few extra ions, and a broad region can be barely enriched yet dominate the total
  excess. Area shares are given so the two never get confused.</div>
  <div id="maps"></div>
</div></section>

<section id="filter"><div class="finding">
  <div><div class="eyebrow">Finding 4</div><h2>The electrolyte is a low-pass filter with a half-transmission wavelength near 8&nbsp;&Aring;</h2></div>
  <p class="claim">Comparing every Fourier mode of the metal's charge response with the same mode of the anion
  response, transmission falls from full at long wavelength to well under 1% at the atomic scale. Half the amplitude
  is gone by &lambda;&nbsp;&asymp;&nbsp;7.7&nbsp;&Aring;, about two and a half nearest-neighbour spacings, and only
  16% of the metal's lateral contrast survives into the anion map.</p>
  <div class="stats" id="s4"></div>
  <div class="card pad"><div class="chartbox"><div id="c_tf"></div></div>
    <div class="caption">Normalised transmission |F<sub>ion</sub>(k)|&thinsp;/&thinsp;|F<sub>metal</sub>(k)| against
    lateral wavelength, pooled over every geometry with both end potentials. The ratio compares two different
    quantities, so each curve is normalised to its own longest-wavelength bin and only the shape is meaningful. Band
    is the interquartile range.</div></div>
  <div class="note"><b>What survives and what does not.</b> The Au&ndash;Au spacing is 2.94&nbsp;&Aring;: atomic
  corrugation of the metal charge reaches the ions at under 2% of its relative amplitude. Defect-scale features
  &mdash; a kink repeat of 10&nbsp;&Aring;, an island of comparable width &mdash; pass at 35% or more. A linearised
  Poisson&ndash;Boltzmann slab with the 1&nbsp;M Debye length gives the same qualitative k-filtering, but the measured
  curve is steeper at short wavelength, which is where the 4&nbsp;&Aring; ion-exclusion cavity, not the Debye
  screening, sets the standoff. The reference is a guide, not a fit.</div>
</div></section>

<section id="omega"><div class="finding">
  <div><div class="eyebrow">Finding 5</div><h2>Over &plusmn;0.2&nbsp;V the potential does not re-rank these morphologies</h2></div>
  <p class="claim">Integrating the electron-number curves gives the potential-induced change in relative grand
  potential without ever subtracting two total energies. Across the sampled window the largest cross-morphology value
  is 29&nbsp;meV &mdash; the same size as the scatter among sampled configurations of a single morphology.</p>
  <div class="card pad"><div class="chartbox"><div id="c_omega"></div></div>
    <div class="caption">&Delta;(&Omega;<sub>A</sub>&nbsp;&minus;&nbsp;&Omega;<sub>B</sub>) across U&nbsp;=&nbsp;&plusmn;0.2&nbsp;V
    for same-composition, same-cell pairs of different structures. Negative means A is relatively stabilised as the
    potential moves positive. The shaded band is the median-to-maximum spread among pairs drawn from one morphology's
    own sampled configurations.</div></div>
  <div class="note"><b>What this is and is not.</b> It is the potential-<em>induced</em> part only. It does not say
  which geometry is more stable at the reference potential, which needs a consistent energy baseline the dataset has
  not yet fixed, and it is never applied across different Au counts, which would need a reservoir term. The
  &plusmn;0.5&nbsp;V extension widens the window by 2.5&times; and will scale these values with it.</div>
</div></section>

<section id="gallery"><div class="finding">
  <div><div class="eyebrow">Structure gallery</div><h2>All thirty-three structures, coloured by coordination</h2></div>
  <p class="lede">Each render is the geometry that was actually computed, not a re-generated idealisation. Atoms carry
  the same coordination classes used to cut the anion regions above, so the pictures and the maps read against each
  other directly.</p>
  <div class="filters" id="filters"></div>
  <div class="gal" id="gal"></div>
</div></section>

<section id="method"><div class="finding">
  <div><div class="eyebrow">Method</div><h2>Conventions, and what these numbers do not cover</h2></div>
  <div class="grid2">
    <div class="card pad"><h3>Potential and charge</h3><p class="small" style="margin-top:8px">
      U&nbsp;=&nbsp;&mu;<sub>0</sub>&nbsp;&minus;&nbsp;&mu;<sub>e</sub> with &mu;<sub>0</sub>&nbsp;=&nbsp;&minus;4.9071&nbsp;eV,
      an internal reference: not a potential versus RHE and not a per-structure PZC.
      &sigma;&nbsp;=&nbsp;&minus;e&thinsp;(N<sub>e</sub>&nbsp;&minus;&nbsp;N<sub>e</sub><sup>0</sup>)&thinsp;/&thinsp;A<sub>proj</sub>.
      Every value uses the converged &mu;<sub>e</sub> and N<sub>e</sub> read from the run, never the target, because
      the convergence criterion lets a run stop up to 10&nbsp;meV from its target.</p></div>
    <div class="card pad"><h3>Ion densities</h3><p class="small" style="margin-top:8px">
      Reconstructed from the converged potential and the accessibility mask through the model's own constitutive
      relation, with the finite-size saturation term. No explicit ions are added and no chloride is present: this is
      the non-specific anion response of a continuum electrolyte, not a chemisorption preference.</p></div>
    <div class="card pad"><h3>Integration windows</h3><p class="small" style="margin-top:8px">
      Absolute: everything below z&nbsp;=&nbsp;31&nbsp;&Aring;, which partitions the cell exactly. Boundary-relative:
      10.4&nbsp;&Aring; upward from each column's own accessibility boundary, so a raised island keeps the same
      fraction of the decaying tail as the terrace beside it. Both are reported and the conclusions hold under each.</p></div>
    <div class="card pad"><h3>Known limits</h3><p class="small" style="margin-top:8px">
      Forces carry an egg-box error of order 0.02&nbsp;eV/&Aring; per atom from the production grid settings, so no
      conclusion here rests on small force differences. Perturbed and collectively deformed configurations are a
      designed sampling, not a thermal ensemble, and are used for sensitivity ranges rather than averages. Path images
      are sparse and are not a minimum-energy path.</p></div>
  </div>
</div></section>

</div>
<footer><div class="wrap">
  Generated from <span class="mono">dataset_v1/states.json</span> and the per-column field reductions.
  <span id="stamp"></span>
</div></footer>

<script type="application/json" id="D">__DATA__</script>
<script>
const D = JSON.parse(document.getElementById("D").textContent);
const FAM = ["flat Au(111)","point defect","reconstruction-related","strip step","vicinal step face",
             "kink / edge rearrangement","single-layer island","single-layer pit","composite"];
const FC = {"flat Au(111)":"#4c78a8","point defect":"#e08a2e","reconstruction-related":"#54a24b",
  "strip step":"#9c6bb0","vicinal step face":"#cf4f45","kink / edge rearrangement":"#3fa0a8",
  "single-layer island":"#c8a33a","single-layer pit":"#8a6a4e","composite":"#8d8f94"};
const CSS = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const fmt = (x,n=2) => (x===null||x===undefined||!isFinite(x)) ? "–" : x.toFixed(n);
const el = (t,a={},...k) => {const e=document.createElementNS(t==="svg"||SVGT.has(t)?"http://www.w3.org/2000/svg":"http://www.w3.org/1999/xhtml",t);
  for(const [p,v] of Object.entries(a)) if(v!==null&&v!==undefined) e.setAttribute(p,v);
  for(const c of k.flat()) if(c!==null&&c!==undefined) e.append(c.nodeType?c:document.createTextNode(c)); return e;};
const SVGT = new Set(["g","rect","circle","line","path","text","polyline","polygon","svg","defs","clipPath","tspan"]);

function frame(w,h,m){const s=el("svg",{viewBox:`0 0 ${w} ${h}`,width:w,height:h,role:"img"});
  return {s,w,h,m,X:null,Y:null};}
function axes(F,xd,yd,xl,yl,xt,yt){
  const {s,w,h,m}=F; const X=x=>m.l+(x-xd[0])/(xd[1]-xd[0])*(w-m.l-m.r);
  const Y=y=>h-m.b-(y-yd[0])/(yd[1]-yd[0])*(h-m.t-m.b); F.X=X;F.Y=Y;
  const line=CSS("--line"), mut=CSS("--muted");
  for(const t of xt){s.append(el("line",{x1:X(t),x2:X(t),y1:m.t,y2:h-m.b,stroke:line,"stroke-width":1}));
    s.append(el("text",{x:X(t),y:h-m.b+15,"text-anchor":"middle","font-size":10,fill:mut},String(t)));}
  for(const t of yt){s.append(el("line",{x1:m.l,x2:w-m.r,y1:Y(t),y2:Y(t),stroke:line,"stroke-width":1}));
    s.append(el("text",{x:m.l-7,y:Y(t)+3.4,"text-anchor":"end","font-size":10,fill:mut},String(t)));}
  s.append(el("text",{x:(m.l+w-m.r)/2,y:h-4,"text-anchor":"middle","font-size":10.5,fill:CSS("--ink-2")},xl));
  s.append(el("text",{x:12,y:(m.t+h-m.b)/2,"text-anchor":"middle","font-size":10.5,fill:CSS("--ink-2"),
    transform:`rotate(-90 12 ${(m.t+h-m.b)/2})`},yl));
  return F;}

/* ---- kicker + stamp ---- */
const NG=Object.keys(D.geometries).length, NS=Object.keys(D.structures).length;
const NPT=Object.values(D.geometries).reduce((a,g)=>a+g.pts.length,0);
document.getElementById("kicker").append(
  ...[[NS,"structures"],[NG,"distinct geometries"],[NPT,"converged states"],
      ["−3","potentials: −5.1071 / −4.9071 / −4.7071 eV"]]
     .map(([v,l])=>el("span",{class:"chip"},`${v} ${l}`)));
document.getElementById("stamp").textContent = " · " + Object.keys(D.transmission).length +
  " geometries carry the paired-potential field analysis.";

/* ---- Finding 1 stats ---- */
const zs=Object.values(D.geometries).filter(g=>g.z!==null).map(g=>g.z);
const cs=Object.values(D.geometries).filter(g=>g.C!==null).map(g=>g.C);
const dAll=D.decomposition["all geometries"];
document.getElementById("s1").append(...[
  [`${(1000*(Math.max(...zs)-Math.min(...zs))).toFixed(0)} mV`,"span of the zero-charge potential across all geometries"],
  [`${(100*(Math.max(...cs)-Math.min(...cs))/median(cs)).toFixed(0)}%`,"span of the secant capacitance over the same set"],
  [`${median(cs).toFixed(1)} µF/cm²`,"median secant capacitance"],
  [`${dAll["U=+0.2V"].ratio_pzc_over_capacitance.toFixed(1)}×`,"charge spread from the PZC shift, against the capacitance term, at U = +0.2 V"],
].map(([v,l])=>el("div",{class:"stat"},el("div",{class:"v"},v),el("div",{class:"l"},l))));
function median(a){const b=[...a].sort((x,y)=>x-y);const n=b.length;return n%2?b[(n-1)/2]:(b[n/2-1]+b[n/2])/2;}

/* ---- sigma(U) ---- */
(function(){
  const F=frame(860,420,{l:54,r:14,t:14,b:40});
  axes(F,[-0.56,0.56],[-7,7],"U = μ₀ − μₑ   (V, internal reference)",
       "σ   (µC/cm²)",[-0.5,-0.25,0,0.25,0.5],[-6,-4,-2,0,2,4,6]);
  F.s.append(el("line",{x1:F.X(-0.56),x2:F.X(0.56),y1:F.Y(0),y2:F.Y(0),stroke:CSS("--line-2"),"stroke-width":1.4}));
  for(const g of Object.values(D.geometries)){
    const pts=g.pts.map(p=>`${F.X(p[0]).toFixed(1)},${F.Y(p[1]).toFixed(1)}`).join(" ");
    F.s.append(el("polyline",{points:pts,fill:"none",stroke:FC[g.f]||"#888","stroke-width":1.1,
      "stroke-opacity":.62,"stroke-linecap":"round"}));}
  document.getElementById("c_sigma").append(F.s);
  document.getElementById("leg1").append(...FAM.map(f=>el("span",{},
    el("i",{class:"sw",style:`background:${FC[f]}`}),f)));
})();

/* ---- decomposition table ---- */
(function(){
  const t=document.getElementById("t_decomp");
  t.append(el("thead",{},el("tr",{},...["set","n","spread U_pzc (mV)","spread C (%)","σ spread from PZC","from C","ratio"]
    .map(h=>el("th",{},h)))));
  const tb=el("tbody");
  for(const [k,v] of Object.entries(D.decomposition)){
    const w=v["U=+0.2V"];
    tb.append(el("tr",{},el("td",{},k),el("td",{class:"num"},String(v.n)),
      el("td",{class:"num"},fmt(v.U_pzc_spread_mV,0)),el("td",{class:"num"},fmt(v.C_spread_pct,1)),
      el("td",{class:"num"},fmt(w.from_pzc_shift_uC_per_cm2,2)),el("td",{class:"num"},fmt(w.from_capacitance_uC_per_cm2,2)),
      el("td",{class:"num"},fmt(w.ratio_pzc_over_capacitance,1)+"×")));}
  t.append(tb);
})();

/* ---- dU_pzc bars ---- */
(function(){
  const rows=Object.entries(D.structures).filter(([k,v])=>v.dz!==null&&v.dz!==undefined)
    .sort((a,b)=>a[1].dz-b[1].dz);
  if(!rows.length){document.getElementById("c_dpzc").append(el("p",{class:"small"},"no in-cell flat reference available"));return;}
  const H=30, w=860, m={l:238,r:60,t:26,b:34}, h=m.t+m.b+rows.length*H;
  const F=frame(w,h,m);
  const lo=Math.min(-100,...rows.map(r=>r[1].dz)), hi=Math.max(30,...rows.map(r=>r[1].dz));
  const X=x=>m.l+(x-lo)/(hi-lo)*(w-m.l-m.r); F.X=X;
  F.s.append(el("line",{x1:X(0),x2:X(0),y1:m.t-6,y2:h-m.b,stroke:CSS("--line-2"),"stroke-width":1.4}));
  for(const t of [-100,-75,-50,-25,0,25]){
    F.s.append(el("line",{x1:X(t),x2:X(t),y1:m.t-6,y2:h-m.b,stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:X(t),y:h-m.b+14,"text-anchor":"middle","font-size":10,fill:CSS("--muted")},String(t)));}
  rows.forEach(([k,v],i)=>{
    const y=m.t+i*H+H/2;
    F.s.append(el("rect",{x:Math.min(X(0),X(v.dz)),y:y-7,width:Math.abs(X(v.dz)-X(0)),height:14,rx:2,
      fill:FC[v.f]||"#888","fill-opacity":.85}));
    F.s.append(el("text",{x:m.l-9,y:y+3.6,"text-anchor":"end","font-size":11,fill:CSS("--ink")},k));
    if(v.dC!==null&&v.dC!==undefined){
      F.s.append(el("circle",{cx:w-m.r+18,cy:y,r:4.4,fill:CSS("--surface"),stroke:CSS("--gold"),"stroke-width":1.6}));
      F.s.append(el("text",{x:w-m.r+27,y:y+3.4,"font-size":9.5,fill:CSS("--gold")},(v.dC>0?"+":"")+v.dC.toFixed(0)+"%"));}
  });
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-3,"text-anchor":"middle","font-size":10.5,fill:CSS("--ink-2")},
    "ΔU_pzc against the flat terrace of the same cell   (mV)"));
  F.s.append(el("text",{x:w-m.r+14,y:m.t-12,"font-size":9.5,fill:CSS("--gold")},"ΔC"));
  document.getElementById("c_dpzc").append(F.s);
})();

/* ---- region K by coordination class ---- */
(function(){
  const CLS=[["kink/adatom",CSS("--red")],["edge/rim",CSS("--amber")],["terrace",CSS("--terrace")],
             ["sub-surface/foot",CSS("--sub")]];
  const rows=Object.entries(D.regions).sort((a,b)=>(a[1].regions.all.K)-(b[1].regions.all.K));
  if(!rows.length){document.getElementById("c_regions").append(el("p",{class:"small"},"field reduction still running"));return;}
  const RH=46,w=860,m={l:196,r:96,t:22,b:40},h=m.t+m.b+rows.length*RH;
  const vals=rows.flatMap(([k,v])=>Object.values(v.regions).map(r=>r.K));
  const lo=Math.min(...vals)-0.01, hi=Math.max(...vals)+0.01;
  const X=x=>m.l+(x-lo)/(hi-lo)*(w-m.l-m.r);
  const F=frame(w,h,m);
  for(let t=Math.ceil(lo*20)/20;t<=hi;t+=0.05){
    F.s.append(el("line",{x1:X(t),x2:X(t),y1:m.t-4,y2:h-m.b,stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:X(t),y:h-m.b+14,"text-anchor":"middle","font-size":10,fill:CSS("--muted")},t.toFixed(2)));}
  F.s.append(el("line",{x1:X(1),x2:X(1),y1:m.t-4,y2:h-m.b,stroke:CSS("--line-2"),"stroke-width":1.4}));
  rows.forEach(([k,v],i)=>{
    const y0=m.t+i*RH, present=CLS.filter(([c])=>v.regions[c]);
    const bh=Math.min(8,(RH-16)/Math.max(present.length,1));
    F.s.append(el("text",{x:m.l-9,y:y0+RH/2-2,"text-anchor":"end","font-size":11,fill:CSS("--ink")},k));
    F.s.append(el("text",{x:m.l-9,y:y0+RH/2+10,"text-anchor":"end","font-size":9,fill:CSS("--muted")},
      "U = "+(v.U>0?"+":"")+v.U.toFixed(2)+" V"));
    present.forEach(([c,col],j)=>{
      const r=v.regions[c], y=y0+8+j*(bh+1.6);
      F.s.append(el("rect",{x:X(lo),y:y,width:Math.max(0,X(r.K)-X(lo)),height:bh,rx:1.5,fill:col,"fill-opacity":.9}));
      F.s.append(el("text",{x:X(r.K)+5,y:y+bh-0.6,"font-size":8.6,fill:CSS("--ink-2")},
        r.K.toFixed(3)+"  ("+(100*r.a).toFixed(0)+"% area)"));});
    const a=v.regions.all;
    F.s.append(el("line",{x1:X(a.K),x2:X(a.K),y1:y0+5,y2:y0+RH-6,stroke:CSS("--ink"),"stroke-width":1.2,
      "stroke-dasharray":"3 2"}));
  });
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-3,"text-anchor":"middle","font-size":10.5,fill:CSS("--ink-2")},
    "KΩ  =  anion concentration in that region, relative to bulk"));
  document.getElementById("c_regions").append(F.s);
  document.getElementById("leg3").append(...CLS.map(([c,col])=>el("span",{},el("i",{class:"sw",style:`background:${col}`}),c)),
    el("span",{},el("i",{class:"sw",style:`background:${CSS("--ink")}`}),"cell average"));
})();

/* ---- transfer function ---- */
(function(){
  const curves=Object.values(D.transmission).filter(v=>v.tf&&v.tf.b.length);
  if(!curves.length){document.getElementById("c_tf").append(el("p",{class:"small"},"field reduction still running"));return;}
  const halves=curves.map(v=>v.tf.half).filter(x=>x);
  document.getElementById("s4").append(...[
    [`${median(halves).toFixed(1)} Å`,"wavelength at which half the amplitude is lost (median over geometries)"],
    [`${curves.length}`,"geometries with both end potentials analysed"],
    ["2.94 Å","Au–Au nearest-neighbour spacing, for scale"],
    [`${(100*median(Object.values(D.transmission).map(v=>v.ratio).filter(x=>x))).toFixed(0)}%`,
     "of the metal's lateral contrast that survives into the anion map"],
  ].map(([v,l])=>el("div",{class:"stat"},el("div",{class:"v"},v),el("div",{class:"l"},l))));
  const bins=new Map();
  for(const v of curves) for(const [lam,T] of v.tf.b){
    const key=Math.round(Math.log(lam)*2.2); if(!bins.has(key)) bins.set(key,[]); bins.get(key).push([lam,T]);}
  const pooled=[...bins.entries()].map(([k,a])=>{
    const L=median(a.map(x=>x[0])), Ts=a.map(x=>x[1]).sort((p,q)=>p-q);
    const q=f=>Ts[Math.min(Ts.length-1,Math.floor(f*Ts.length))];
    return {lam:L,T:median(Ts),lo:q(.25),hi:q(.75),n:a.length};}).sort((a,b)=>a.lam-b.lam);
  const w=860,h=400,m={l:56,r:16,t:16,b:42};
  const F=frame(w,h,m);
  const lx=x=>Math.log10(x), xd=[lx(0.7),lx(24)];
  const X=x=>m.l+(lx(x)-xd[0])/(xd[1]-xd[0])*(w-m.l-m.r);
  const Y=y=>h-m.b-(Math.log10(Math.max(y,2e-3))-Math.log10(2e-3))/(0-Math.log10(2e-3))*(h-m.t-m.b);
  for(const t of [1,2,3,5,10,20]){
    F.s.append(el("line",{x1:X(t),x2:X(t),y1:m.t,y2:h-m.b,stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:X(t),y:h-m.b+15,"text-anchor":"middle","font-size":10,fill:CSS("--muted")},String(t)));}
  for(const t of [0.002,0.01,0.05,0.2,1]){
    F.s.append(el("line",{x1:m.l,x2:w-m.r,y1:Y(t),y2:Y(t),stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:m.l-7,y:Y(t)+3.4,"text-anchor":"end","font-size":10,fill:CSS("--muted")},
      t>=1?"1":String(t)));}
  F.s.append(el("polygon",{points:pooled.map(p=>`${X(p.lam).toFixed(1)},${Y(p.hi).toFixed(1)}`).join(" ")+" "+
    pooled.slice().reverse().map(p=>`${X(p.lam).toFixed(1)},${Y(p.lo).toFixed(1)}`).join(" "),
    fill:CSS("--teal"),"fill-opacity":.16}));
  F.s.append(el("polyline",{points:pooled.map(p=>`${X(p.lam).toFixed(1)},${Y(p.T).toFixed(1)}`).join(" "),
    fill:"none",stroke:CSS("--teal"),"stroke-width":2.2,"stroke-linejoin":"round"}));
  for(const p of pooled) F.s.append(el("circle",{cx:X(p.lam),cy:Y(p.T),r:3.2,fill:CSS("--teal")}));
  const hm=median(halves);
  F.s.append(el("line",{x1:X(hm),x2:X(hm),y1:m.t,y2:h-m.b,stroke:CSS("--gold"),"stroke-width":1.5,"stroke-dasharray":"5 3"}));
  F.s.append(el("text",{x:X(hm)+6,y:m.t+13,"font-size":10,fill:CSS("--gold")},`λ½ ≈ ${hm.toFixed(1)} Å`));
  F.s.append(el("line",{x1:X(2.94),x2:X(2.94),y1:m.t,y2:h-m.b,stroke:CSS("--red"),"stroke-width":1.2,"stroke-dasharray":"2 3"}));
  F.s.append(el("text",{x:X(2.94)+5,y:h-m.b-8,"font-size":9.5,fill:CSS("--red")},"Au–Au 2.94 Å"));
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-4,"text-anchor":"middle","font-size":10.5,fill:CSS("--ink-2")},
    "lateral wavelength λ   (Å, log scale)"));
  F.s.append(el("text",{x:13,y:(m.t+h-m.b)/2,"text-anchor":"middle","font-size":10.5,fill:CSS("--ink-2"),
    transform:`rotate(-90 13 ${(m.t+h-m.b)/2})`},"normalised transmission (log)"));
  document.getElementById("c_tf").append(F.s);
})();

/* ---- dOmega ---- */
(function(){
  const seen=new Set(), rows=[];
  for(const p of D.pairs){
    const a=p.a.split("__"), b=p.b.split("__");
    if(a[1]!=="ideal"&&a[1]!=="relaxed") continue;
    if(b[1]!=="ideal"&&b[1]!=="relaxed") continue;
    const key=[a[0],b[0]].sort().join("|"); if(seen.has(key)) continue; seen.add(key);
    rows.push({l:`${a[0]} ${a[1]}  vs  ${b[0]} ${b[1]}`,v:p.d_relative_Omega_eV,z:p.dU_pzc_mV,n:p.n_atoms});}
  rows.sort((x,y)=>Math.abs(y.v)-Math.abs(x.v));
  const R=rows.slice(0,14);
  if(!R.length){document.getElementById("c_omega").append(el("p",{class:"small"},"no cross-structure pairs"));return;}
  const H=30,w=860,m={l:300,r:70,t:30,b:36},h=m.t+m.b+R.length*H;
  const lim=Math.max(0.035,...R.map(r=>Math.abs(r.v)))*1.12;
  const X=x=>m.l+(x+lim)/(2*lim)*(w-m.l-m.r);
  const F=frame(w,h,m);
  const sc=D.same_structure_pair_scale;
  F.s.append(el("rect",{x:X(-sc.max),y:m.t-8,width:X(sc.max)-X(-sc.max),height:h-m.b-m.t+8,
    fill:CSS("--ink-2"),"fill-opacity":.09}));
  F.s.append(el("text",{x:X(0),y:m.t-13,"text-anchor":"middle","font-size":9.5,fill:CSS("--ink-2")},
    `within one morphology: |ΔΩ| up to ${(1000*sc.max).toFixed(0)} meV`));
  for(const t of [-0.03,-0.02,-0.01,0,0.01,0.02,0.03]){
    if(Math.abs(t)>lim) continue;
    F.s.append(el("line",{x1:X(t),x2:X(t),y1:m.t-8,y2:h-m.b,stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:X(t),y:h-m.b+14,"text-anchor":"middle","font-size":10,fill:CSS("--muted")},
      (1000*t).toFixed(0)));}
  F.s.append(el("line",{x1:X(0),x2:X(0),y1:m.t-8,y2:h-m.b,stroke:CSS("--line-2"),"stroke-width":1.4}));
  R.forEach((r,i)=>{const y=m.t+i*H+H/2;
    F.s.append(el("rect",{x:Math.min(X(0),X(r.v)),y:y-8,width:Math.abs(X(r.v)-X(0)),height:16,rx:2,
      fill:r.v<0?CSS("--teal"):CSS("--gold"),"fill-opacity":.85}));
    F.s.append(el("text",{x:m.l-9,y:y+3.6,"text-anchor":"end","font-size":10.5,fill:CSS("--ink")},r.l));
    F.s.append(el("text",{x:w-m.r+8,y:y+3.6,"font-size":9.5,fill:CSS("--muted")},
      r.z===null?"":`ΔU₀ ${r.z>0?"+":""}${r.z.toFixed(0)} mV`));});
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-3,"text-anchor":"middle","font-size":10.5,fill:CSS("--ink-2")},
    "Δ(Ωₐ − Ωᵦ) across U = ±0.2 V   (meV)"));
  document.getElementById("c_omega").append(F.s);
})();

/* ---- gallery ---- */
(function(){
  const fams=FAM.filter(f=>Object.values(D.structures).some(s=>s.f===f));
  const bar=document.getElementById("filters"), gal=document.getElementById("gal");
  let active="all";
  const mk=(id,label)=>{const b=el("button",{type:"button","aria-pressed":String(active===id),"data-f":id},label);
    b.addEventListener("click",()=>{active=id;[...bar.children].forEach(c=>c.setAttribute("aria-pressed",
      String(c.dataset.f===id)));render();}); return b;};
  bar.append(mk("all",`all ${Object.keys(D.structures).length}`),
    ...fams.map(f=>mk(f,`${f} (${Object.values(D.structures).filter(s=>s.f===f).length})`)));
  function render(){
    gal.textContent="";
    const items=Object.entries(D.structures).filter(([k,v])=>active==="all"||v.f===active)
      .sort((a,b)=>FAM.indexOf(a[1].f)-FAM.indexOf(b[1].f)||a[0].localeCompare(b[0]));
    for(const [k,v] of items){
      const tags=[el("span",{class:"tag"},`${v.n} Au`),el("span",{class:"tag"},`A ${v.A} Å²`),
        el("span",{class:"tag"},v.f)];
      if(v.cn.kink) tags.push(el("span",{class:"tag k"},`${v.cn.kink} CN≤6`));
      if(v.cn.edge) tags.push(el("span",{class:"tag e"},`${v.cn.edge} CN 7–8`));
      if(v.C!==null&&v.C!==undefined) tags.push(el("span",{class:"tag"},`C ${v.C.toFixed(1)} µF/cm²`));
      if(v.z!==null&&v.z!==undefined) tags.push(el("span",{class:"tag"},
        `U₀ ${v.z>0?"+":""}${(1000*v.z).toFixed(0)} mV`));
      if(v.dz!==null&&v.dz!==undefined) tags.push(el("span",{class:"tag"},
        `ΔU₀ ${v.dz>0?"+":""}${v.dz.toFixed(0)} mV vs flat`));
      gal.append(el("figure",{class:"gcard",style:"margin:0"},
        el("div",{class:"ghead"},el("h3",{},k),el("span",{class:"small mono"},v.cfg)),
        el("img",{src:`gallery/${k}.png`,alt:`${k}: top and side view, atoms coloured by coordination number`,
          loading:"lazy"}),
        el("div",{class:"gmeta"},...tags)));}
  }
  render();
})();

/* ---- per-structure maps ---- */
(function(){
  const host=document.getElementById("maps"); if(!D.mapped||!D.mapped.length) return;
  host.append(el("h3",{style:"margin:26px 0 4px"},"Anion maps, structure by structure"),
    el("p",{class:"small",style:"margin:0 0 14px;max-width:66ch"},
      "Per-column anion excess at the most positive sampled potential, its change across the sampled window, the "+
      "metal's own change over the same step, and the coordination map the regions are cut from."));
  const g=el("div",{class:"gal"});
  for(const k of D.mapped) g.append(el("figure",{class:"gcard",style:"margin:0"},
    el("div",{class:"ghead"},el("h3",{},k)),
    el("img",{src:`maps/${k}.png`,alt:`${k}: anion excess, its potential-driven change, the metal charge change and the coordination map`,
      loading:"lazy"})));
  host.append(g);
})();
</script>
"""


def main():
    data = json.load(open(f"{WEB}/data.json"))
    mapped = sorted(f[:-4] for f in os.listdir(f"{ROOT}/analysis/maps")) if os.path.exists(f"{ROOT}/analysis/maps") else []
    data["mapped"] = [m for m in mapped if m in data["structures"]]
    out = HTML.replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/"))
    open(f"{WEB}/index.html", "w").write(out)
    print(f"{os.path.getsize(f'{WEB}/index.html')/1024:.0f} kB -> {WEB}/index.html "
          f"({len(data['structures'])} structures, {len(data['mapped'])} maps)")


if __name__ == "__main__":
    main()
