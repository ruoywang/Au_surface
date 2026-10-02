#!/usr/bin/env python3
"""Build the standalone HTML report (Chinese) from analysis/web/data.json, inlined because the artifact CSP blocks fetch.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_web_report.py
Output: analysis/web/index.html  (publish with the gallery and map PNGs as supporting files)
"""
import json
import os

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
WEB = f"{ROOT}/analysis/web"

HTML = r"""<title>Au surface morphology and charging response</title>
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Noto+Serif:wght@400;500;600&family=Inter:wght@300;400;500;700&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{
  --bg:#f6f4ef; --surface:#fffefb; --surface-2:#f0ece4; --ink:#15181c; --ink-2:#565d65; --muted:#8a9199;
  --line:#e0dbd1; --line-2:#cdc6b9;
  --gold:#96650d; --gold-soft:#e8d9b4; --teal:#1d6a85; --red:#a8382a; --amber:#c4841f;
  --terrace:#b99a3e; --sub:#79828c;
  --shadow:0 1px 2px rgba(30,24,10,.05),0 6px 20px -12px rgba(30,24,10,.18);
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
  --sans:"Inter",system-ui,-apple-system,"Segoe UI",Helvetica,Arial,sans-serif;
  --serif:"Noto Serif",Georgia,"Times New Roman",serif;
}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){
  --bg:#121417; --surface:#191c21; --surface-2:#21252b; --ink:#e9ebee; --ink-2:#a6adb5; --muted:#79818a;
  --line:#2a2e35; --line-2:#3a3f47;
  --gold:#d9a743; --gold-soft:#4a3c1c; --teal:#5fb2cd; --red:#d9705f; --amber:#dda23e;
  --terrace:#c9ab52; --sub:#8b949e;
  --shadow:0 1px 2px rgba(0,0,0,.3),0 8px 24px -14px rgba(0,0,0,.7);
}}
:root[data-theme="dark"]{
  --bg:#121417; --surface:#191c21; --surface-2:#21252b; --ink:#e9ebee; --ink-2:#a6adb5; --muted:#79818a;
  --line:#2a2e35; --line-2:#3a3f47;
  --gold:#d9a743; --gold-soft:#4a3c1c; --teal:#5fb2cd; --red:#d9705f; --amber:#dda23e;
  --terrace:#c9ab52; --sub:#8b949e;
  --shadow:0 1px 2px rgba(0,0,0,.3),0 8px 24px -14px rgba(0,0,0,.7);
}
*{box-sizing:border-box}
body{background:var(--bg);color:var(--ink);font-family:var(--sans);font-size:15.5px;line-height:1.78;margin:0;
  font-weight:400;-webkit-font-smoothing:antialiased}
.wrap{max-width:1180px;margin:0 auto;padding-inline:20px;padding-block:0}
/* Figures break out of the 1180px reading column. A 12-inch-wide render shown at 560px puts its 9pt labels
   below 7 CSS px, which is not readable at all; at ~1700px the same labels land near 18px. Text stays in the
   narrow column because long measures are hard to read. */
.bleed{width:min(1780px,96vw);margin-left:50%;transform:translateX(-50%)}

h1,h2,h3{font-family:var(--serif);font-weight:600;text-wrap:balance;margin:0;letter-spacing:.01em}
h1{font-size:clamp(28px,5vw,46px);line-height:1.22}
h2{font-size:clamp(20px,2.9vw,27px);line-height:1.35;margin-bottom:.4em}
h3{font-size:17px;line-height:1.5}
p{margin:0 0 .9em}
a{color:var(--teal)}
.lede{font-size:clamp(15.5px,2vw,17.5px);color:var(--ink-2);max-width:40em;line-height:1.85}
.eyebrow{font-family:var(--mono);font-size:11px;letter-spacing:.13em;text-transform:uppercase;color:var(--gold);
  margin-bottom:.6em}
.mono{font-family:var(--mono);font-variant-numeric:tabular-nums}
small,.small{font-size:13px;color:var(--ink-2);line-height:1.75}
.caption{font-size:13px;color:var(--muted);line-height:1.7;margin-top:.6em}

header.top{border-bottom:1px solid var(--line);background:
  radial-gradient(1200px 340px at 12% -40%,var(--gold-soft),transparent 65%),var(--surface)}
header.top .wrap{padding-block:46px 34px}
.kicker{display:flex;flex-wrap:wrap;gap:8px;margin-top:20px}
.chip{font-family:var(--mono);font-size:11.5px;border:1px solid var(--line-2);border-radius:999px;
  padding:3px 11px;color:var(--ink-2);background:var(--surface)}

nav.sticky{position:sticky;top:env(safe-area-inset-top,0px);z-index:40;
  background:color-mix(in srgb,var(--bg) 92%,transparent);backdrop-filter:blur(10px);border-bottom:1px solid var(--line)}
nav.sticky .wrap{display:flex;gap:2px;overflow-x:auto;padding-block:7px;scrollbar-width:none}
nav.sticky .wrap::-webkit-scrollbar{display:none}
nav.sticky a{font-size:12.5px;color:var(--ink-2);text-decoration:none;white-space:nowrap;padding:5px 11px;border-radius:6px}
nav.sticky a:hover{background:var(--surface-2);color:var(--ink)}

section{padding-block:54px;border-bottom:1px solid var(--line)}
section:last-of-type{border-bottom:0}
.finding{display:grid;grid-template-columns:minmax(0,1fr);gap:22px}
.claim{font-family:var(--serif);font-size:clamp(17px,2.2vw,20px);line-height:1.8;color:var(--ink);
  border-left:3px solid var(--gold);padding-left:17px;max-width:40em;margin:4px 0 2px}
.card{background:var(--surface);border:1px solid var(--line);border-radius:10px;box-shadow:var(--shadow)}
.pad{padding:18px}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(310px,1fr));gap:18px}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(158px,1fr));gap:1px;background:var(--line);
  border:1px solid var(--line);border-radius:10px;overflow:hidden}
.stat{background:var(--surface);padding:14px 16px}
.stat .v{font-family:var(--mono);font-size:22px;font-weight:500;letter-spacing:-.02em;color:var(--ink)}
.stat .l{font-size:12px;color:var(--ink-2);margin-top:3px;line-height:1.6}
.tablewrap{overflow-x:auto;border:1px solid var(--line);border-radius:10px;background:var(--surface)}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{text-align:left;padding:8px 12px;border-bottom:1px solid var(--line);white-space:nowrap}
th{font-size:11.5px;color:var(--muted);font-weight:500;position:sticky;top:0;background:var(--surface)}
td.num{font-family:var(--mono);font-variant-numeric:tabular-nums;text-align:right}
tr:last-child td{border-bottom:0}
tbody tr:hover{background:var(--surface-2)}

.legend{display:flex;flex-wrap:wrap;gap:10px 18px;font-size:12.5px;color:var(--ink-2);margin-top:12px}
.legend span{display:inline-flex;align-items:center;gap:6px}
.sw{width:11px;height:11px;border-radius:3px;display:inline-block;flex:none}
.filters{display:flex;flex-wrap:wrap;gap:6px;margin:16px 0 20px}
.filters button{font-size:12.5px;border:1px solid var(--line-2);background:var(--surface);
  color:var(--ink-2);border-radius:999px;padding:4px 13px;cursor:pointer;font-family:inherit}
.filters button[aria-pressed="true"]{background:var(--ink);color:var(--bg);border-color:var(--ink)}
.filters button:focus-visible{outline:2px solid var(--teal);outline-offset:2px}

/* ONE figure per row. Two-up put every structure figure at ~560px, where its labels were illegible. */
.gal{display:grid;grid-template-columns:1fr;gap:26px}
.gcard{background:var(--surface);border:1px solid var(--line);border-radius:10px;overflow:hidden;box-shadow:var(--shadow)}
.gcard button.zoom{display:block;width:100%;padding:0;border:0;background:#fff;cursor:zoom-in;line-height:0}
.gcard button.zoom:focus-visible{outline:2px solid var(--teal);outline-offset:-2px}
.gcard img{display:block;width:100%;height:auto}
.ghead{display:flex;justify-content:space-between;align-items:baseline;gap:10px;padding:13px 16px 0}
.ghead h3{font-size:16.5px;font-family:var(--mono);font-weight:500}
.gmeta{display:flex;flex-wrap:wrap;gap:5px;padding:10px 16px 15px}
.tag{font-family:var(--mono);font-size:10.5px;padding:2px 7px;border-radius:5px;background:var(--surface-2);
  color:var(--ink-2);border:1px solid var(--line)}
.tag.k{background:color-mix(in srgb,var(--red) 16%,transparent);border-color:transparent;color:var(--red)}
.tag.e{background:color-mix(in srgb,var(--amber) 18%,transparent);border-color:transparent;color:var(--amber)}
.hint{font-size:12px;color:var(--muted);padding:0 16px 12px}
svg{display:block;max-width:100%;height:auto}
/* A chart scaled down to phone width would put its 10-11px labels below legibility, so below 720px the charts
   keep their drawn size and the box scrolls sideways instead. */
.chartbox{overflow-x:auto;-webkit-overflow-scrolling:touch}
@media (max-width:720px){
  .chartbox svg{max-width:none}
  .chartbox::after{content:"\2190 scroll sideways for the full chart";display:block;font-size:11.5px;color:var(--muted);
    padding:6px 2px 0;position:sticky;left:0}
}
.note{background:var(--surface-2);border:1px solid var(--line);border-left:3px solid var(--teal);border-radius:8px;
  padding:14px 17px;font-size:14px;color:var(--ink-2);line-height:1.8}
.note b{color:var(--ink);font-weight:500}
dialog.lb{border:0;padding:0;background:transparent;max-width:98vw;max-height:96vh;overflow:visible}
dialog.lb::backdrop{background:rgba(12,14,17,.86)}
dialog.lb img{display:block;max-width:98vw;max-height:88vh;width:auto;height:auto;background:#fff;border-radius:8px}
dialog.lb .bar{display:flex;justify-content:space-between;align-items:center;gap:14px;color:#e8eaec;
  font-size:13px;padding:9px 4px 0}
dialog.lb button{background:transparent;border:1px solid rgba(232,234,236,.35);color:#e8eaec;border-radius:6px;
  padding:4px 12px;font-size:12.5px;cursor:pointer;font-family:inherit}
footer{padding-block:40px 60px;color:var(--muted);font-size:12.5px}
@media (max-width:560px){section{padding-block:38px}}
</style>

<header class="top"><div class="wrap">
  <div class="eyebrow">Constant-potential DFT · VASP + VASPsol++ · 1 M implicit electrolyte</div>
  <h1>Au surface morphology and charging response</h1>
  <p class="lede" style="margin-top:14px">Au(111) terraces and the related vicinal faces Au(211), (221), (332)
  and (554): how surface morphology changes the metal's ability to hold charge, and how much of that
  atomic-scale difference actually reaches the ions in the liquid.</p>
  <div class="kicker" id="kicker"></div>
  <div class="note" style="margin-top:20px;max-width:none"><b>Scope (every chart on this page uses the same
  one).</b> <span id="scope"></span></div>
</div></header>

<nav class="sticky"><div class="wrap">
  <a href="#origin">Where the states come from</a>
  <a href="#pzc">1 · What moves is the zero-charge point</a>
  <a href="#which">2 · Which morphologies, which way</a>
  <a href="#where">3 · Where the anions actually gather</a>
  <a href="#filter">4 · The electrolyte as a low-pass filter</a>
  <a href="#omega">5 · Potential and relative stability</a>
  <a href="#gallery">Structure gallery</a>
  <a href="#maps">Anion maps</a>
  <a href="#method">Method and limits</a>
  <a href="#revlog">Revision log</a>
</div></nav>

<div class="wrap">

<section id="origin"><div class="finding">
  <div><div class="eyebrow">Sampling</div><h2>How one built geometry becomes many computed states</h2></div>
  <p class="lede">The <span data-n="ngeom2"></span> geometries are <b>not</b> <span data-n="ngeom2"></span>
  separately designed defects. 33 of them are built by hand (see the gallery below); every other one is a
  relaxation of those at the common reference potential, or a random displacement, collective deformation or
  path image generated from that relaxed geometry.</p>
  <div class="card pad bleed" id="originfig"></div>
  <div class="note"><b>What this sampling is, and what it is not.</b> The structure classes are built from
  crystallography and geometric operations; the relaxations supply a local reference state; the perturbations
  and path images cover off-equilibrium configurations. They are <b>not</b> a statistic of how often such
  configurations occur in molecular dynamics, and <b>not</b> the set of stable morphologies DFT would find on
  its own. Each geometry then carries several potentials, but not every geometry has all five.</div>
</div></section>

<section id="pzc"><div class="finding">
  <div><div class="eyebrow">Finding 1</div><h2>At a fixed potential, the charge differences show up mostly as
  a shift of the zero-charge point</h2></div>
  <p class="claim">Over the geometries whose sampled points bracket &sigma;&nbsp;=&nbsp;0, the zero-charge
  point spans about <span data-n="pzcspread"></span> while the secant capacitance spans about
  <span data-n="cspread"></span>. Within this sample and this potential window, the charge scale that goes
  with the shift in zero-charge point is roughly an order of magnitude larger than the one that goes with the
  change in capacitance.</p>
  <div class="stats" id="s1"></div>
  <div class="card pad bleed"><div class="chartbox"><div id="c_sigma"></div></div>
    <div class="caption">Surface charge density against the internal potential
    U&nbsp;=&nbsp;&mu;<sub>0</sub>&nbsp;&minus;&nbsp;&mu;<sub>e</sub> for every geometry, taken from the
    <em>actually converged</em> &mu;<sub>e</sub> and N<sub>e</sub>. All five potentials are now in scope. The
    curves are close to parallel, so family to family the difference is mostly a sideways shift, and they are
    also close to straight: over the 79 geometries with all five points the worst departure from a line is
    0.081&nbsp;&micro;C/cm² against a range of 7.32, about 1%. One secant capacitance is therefore meaningful
    across the whole window, which the &plusmn;0.2&nbsp;V data could not establish. This is a statement about
    the <b>total</b> surface charge of a fixed geometry at these sampled points only; it does not imply that
    forces, the local potential, the bound charge or the ion concentration are linear in U.</div>
    <div class="legend" id="leg1"></div></div>
  <div class="note"><b>What this ratio is, and what it is not.</b> A defect carrying more positive charge at a
  given potential does not mean it is easier to polarise. Expanding
  &sigma;<sub>i</sub>&nbsp;=&nbsp;C<sub>i</sub>&thinsp;(U&nbsp;&minus;&nbsp;Z<sub>i</sub>) about the reference
  gives &delta;&sigma;<sub>i</sub>&nbsp;=&nbsp;&minus;C<sub>0</sub>&delta;Z<sub>i</sub>&nbsp;+&nbsp;(U&minus;Z<sub>0</sub>)&delta;C<sub>i</sub>&nbsp;&minus;&nbsp;&delta;C<sub>i</sub>&delta;Z<sub>i</sub>:
  besides the first two terms there is a cross term, and across geometries &delta;C and &delta;Z are themselves
  correlated. The table below gives only the <b>scale</b> of the first two terms, estimated from their
  <b>ranges</b>, and their ratio. It is <b>not a variance decomposition and must not be read as "term X
  explains Y per cent"</b>.
  <br><br>Nor is the capacitance "almost unchanged": A3, three adatoms, is about 15% above the flat member of
  its own cell. The accurate statement is that within this potential window the charge difference from the
  capacitance change is smaller than the one from the shift in zero-charge point.</div>
  <div class="tablewrap"><table id="t_decomp"></table></div>
  <div class="caption">Comparisons across cells carry a numerical offset: the audit recorded two <em>flat</em>
  cells differing by 72&nbsp;meV in neutral &mu;<sub>e</sub> from cell shape, k-mesh and PREC alone. Inside one
  cell that offset is common to every row, which is why grouping by cell isolates the morphology effect
  better, though it does not mean every numerical error cancels. The full set of <span data-n="ngeom"></span>
  geometries also contains deliberate perturbations, collective deformations and different cells and sampling
  settings, so the overall spread should not be attributed entirely to "the defect type itself".</div>
</div></section>

<section id="which"><div class="finding">
  <div><div class="eyebrow">Finding 2</div><h2>In these matched pairs, protrusions push the zero-charge point
  negative and vacancies push it positive</h2></div>
  <p class="claim">Against the flat terrace <em>in the same cell</em>: the adatoms and steps in this dataset
  lower the zero-charge point by 29&ndash;97&nbsp;mV, the three vacancy configurations raise it by
  8&ndash;13&nbsp;mV, and the two stacking reconstructions barely move it. The sign agrees with the classical
  work-function argument, in which a protrusion smooths the electron spill-out and lowers the work function
  while a depression does the opposite. This describes these configurations; it is not a general rule for all
  protrusions and depressions.</p>
  <div class="card pad"><div class="chartbox"><div id="c_dpzc"></div></div>
    <div class="caption">Ideal configurations only, each against the flat member of its own cell. Bar length
    is the shift in zero-charge point; the dot on the right marks the capacitance change for the same
    pair.</div></div>
  <div class="note"><b>The capacitance moves the same way, but far less.</b> Adding atoms matters more than
  removing them: three adatoms +15%, a single adatom +8%, a step +8%, while a single vacancy is only +1% and
  both reconstructions stay within 1%. Roughness that protrudes into the electrolyte both lowers the
  zero-charge point and raises the capacitance.
  <br><br><b>Mind the sign.</b> From &sigma;&nbsp;=&nbsp;C&thinsp;(U&nbsp;&minus;&nbsp;U<sub>pzc</sub>), at
  fixed U, &Delta;&sigma;&nbsp;=&nbsp;&minus;C&thinsp;&Delta;U<sub>pzc</sub>: <b>a negative shift of the
  zero-charge point means more positive charge overall at the same potential</b>. The data agree. A1-fcc sits
  65&nbsp;mV below T-4x4 in the same cell, and at U&nbsp;&asymp;&nbsp;+0.2&nbsp;V its &sigma; is +3.18 against
  +2.22&nbsp;&micro;C/cm².</div>
</div></section>

<section id="where"><div class="finding">
  <div><div class="eyebrow">Finding 3</div><h2>Charging more overall and gathering more ions somewhere are not
  the same thing</h2></div>
  <p class="claim">Partitioning the anion excess by the coordination number of the nearest surface atom, the
  region-averaged enrichment over the under-coordinated sites (CN&nbsp;&le;&nbsp;8, area-weighted) is
  <em>lower</em> than over the terrace of the same structure: <span data-n="uc_main"></span>. Each geometry is
  taken at its most positive sampled potential, ideal and relaxed geometries only. The exceptions are listed
  below, and the two integration conventions give the same ordering.</p>
  <div class="note" style="margin-top:0"><b>Exceptions.</b> <span data-n="uc_exc"></span>
  This count used to be written into the prose as 37/41 with four exceptions; it is now recomputed from
  <code>regions.json</code> at build time, with its definition (CN&nbsp;&le;&nbsp;8 area-weighted, most
  positive sampled potential, ideal and relaxed geometries only) written out alongside it.</div>
  <p class="claim" style="border-left-color:var(--teal)">The per-column spatial correlation supports the same
  point: correlating the metal's positive-charge gain (&minus;&Delta;n<sub>e</sub>) against the anion gain
  (&Delta;&Gamma;<sub>&minus;</sub>) column by column, <span data-n="corr"></span>. This describes the
  <b>spatial covariation</b> over the whole map, that the two undulations are broadly opposite in sign. It
  does not give peak positions, and does not guarantee that every structure's maxima are displaced; for that,
  see the individual maps below.</p>
  <div class="card pad"><div class="chartbox"><div id="c_regions"></div></div>
    <div class="caption">Enrichment ratio in the boundary-relative window,
    K<sub>&Omega;</sub>&nbsp;=&nbsp;&int;n<sub>&minus;</sub>&thinsp;/&thinsp;(n<sub>b</sub>&int;S<sub>ion</sub>),
    at each structure's most positive sampled potential. Bars are grouped by coordination class; the dashed
    line is that structure's whole-cell mean. This chart draws one representative geometry per structure,
    ideal or relaxed, so 33 rows; the count quoted above covers both the ideal and the relaxed geometry of the
    same structure.</div>
    <div class="legend" id="leg3"></div></div>
  <div class="note"><b>Overall charging and the local ion distribution do not map onto each other simply.</b>
  Even when the zero-charge point shifts negative and the metal carries more total positive charge at the same
  potential, the liquid region above the under-coordinated atoms need not have a higher region-averaged anion
  concentration. The local distribution is set jointly by the self-consistent potential, the ion
  accessibility and the geometric partition.
  <br><br>K and &Gamma; answer different questions and both are given: a small region can be highly
  concentrated yet hold very few extra ions, while a large region can be barely enriched and still supply most
  of the total excess. The area fraction of each region is shown alongside to keep the two apart. A region
  average can only say that the mean concentration there is lower; where the hot spots are is supported by the
  per-column correlation above and by the maps below.</div>
</div></section>

<section id="filter"><div class="finding">
  <div><div class="eyebrow">Finding 4</div><h2>Short-wavelength spatial structure is strongly damped in the
  ion response</h2></div>
  <p class="claim">Comparing each Fourier mode of the metal's charge response with the same mode of the anion
  response, the relative spectral amplitude falls monotonically as the wavelength gets <b>shorter</b>: at the
  atomic scale (&lambda;&nbsp;&asymp;&nbsp;1&nbsp;&Aring;) only a fraction of a per cent of the reference
  survives, while at the defect scale (&lambda;&nbsp;&asymp;&nbsp;9&nbsp;&Aring;) more than 30% does. This is
  consistent with a picture of spatial screening and smoothing.</p>
  <div class="stats" id="s4"></div>
  <div class="card pad"><div class="chartbox"><div id="c_tf"></div></div>
    <div class="caption">Relative spectral amplitude
    |F<sub>ion</sub>(k)|&thinsp;/&thinsp;|F<sub>metal</sub>(k)| against in-plane wavelength, pooled over every
    geometry that has both end potentials. Electron areal density and ion areal density share dimensions but
    are different physical quantities, so their ratio has no absolute scale; each curve is therefore
    normalised to <b>its own longest-wavelength bin</b>. That is an analysis choice, not something the units
    force, so the 1.0 at the long-wavelength end is a definition and <b>does not mean 100% physical
    transmission</b>; different cells also admit different smallest wavevectors, so the reference bin is not
    identical across curves. Only the shape of the curve is meaningful. The shaded band is the interquartile
    range.</div></div>
  <div class="note"><b>Both numbers are indicators, not material constants.</b> Each depends on the
  normalisation and binning described above and should be read by its definition:
  <br>· <b><span data-n="half"></span></b>: the wavelength at which the relative spectral amplitude falls to
  half of its own reference. This is an empirical indicator <em>under the present Fourier binning, filtering
  and per-geometry normalisation</em>, <b>not a universal resolution of the electrolyte</b>.
  <br>· <b><span data-n="ratio"></span></b>: the ratio of the two normalised contrasts
  (max&minus;min)/&lang;|f|&rang;. It is <b>not</b> the percentage of charge, of Fourier modes or of
  information that survives.
  <br><br><b>The short-wavelength damping is not attributed to any single mechanism.</b> The model contains a
  non-local cavity, a non-linear dielectric and the ion response, all coupled. The observed fall of amplitude
  with wavelength is consistent with spatial screening, but one ratio curve at one parameter set cannot
  establish that the short-wavelength end is set by the 4&nbsp;&Aring; cavity rather than by Debye screening.
  The linearised Poisson&ndash;Boltzmann slab is drawn as a trend reference, not as a fit.</div>
</div></section>

<section id="omega"><div class="finding">
  <div><div class="eyebrow">Finding 5</div><h2>What the potential contributes to the relative grand potential
  of same-composition configurations</h2></div>
  <p class="claim">Integrating the electron-number curve gives the potential-<em>induced</em> change in
  relative grand potential without ever subtracting two total energies. Over the full &plusmn;0.5&nbsp;V
  window the potential-induced change in relative grand potential reaches tens to over a hundred meV between
  configurations (largest cross-morphology value <span data-n="dom"></span>), and the charging contribution
  from rearrangements <em>within</em> one morphology can be comparable to or larger than the cross-morphology
  contribution (up to <span data-n="domsame"></span>). <b>Without a baseline for the relative grand potential
  at the reference potential, the full stability ranking cannot be judged.</b> The two numbers above are
  potential-induced changes for different configuration pairs, not full relative grand potentials and not
  error bars, so comparing them says nothing about whether a ranking changes: if two morphologies differ by
  10&nbsp;meV at the reference potential, 89&nbsp;meV re-ranks them; if by 1&nbsp;eV, it does not.</p>
  <div class="card pad"><div class="chartbox"><div id="c_omega"></div></div>
    <div class="caption">Pairs of different structures with the same composition in the same cell,
    D&nbsp;=&nbsp;+&int;<sub>&mu;₀&minus;w</sub><sup>&mu;₀+w</sup>&thinsp;[N<sub>A</sub>&minus;N<sub>B</sub>]&thinsp;d&mu;,
    w&nbsp;=&nbsp;0.2&nbsp;V. A negative value means A is relatively stabilised as the potential moves
    positive. The grey band is the spread over pairs drawn from the sampled configurations of one morphology.
    Sign: because &part;&Omega;/&part;&mu;<sub>e</sub>&nbsp;=&nbsp;&minus;N<sub>e</sub> and
    U&nbsp;=&nbsp;&mu;₀&minus;&mu;<sub>e</sub>, U&nbsp;=&nbsp;+w is the <em>lower</em> &mu;.</div></div>
  <div class="note"><b>What this is, and what it is not.</b> It is only the potential-<em>induced</em> part.
  It does not say which configuration is more stable at the reference potential; that needs a common energy
  baseline the dataset has not fixed, and without it there is no way to tell whether the ranking changes. If
  two configurations start a few meV apart, tens of meV are enough to re-rank them; if they start 1&nbsp;eV
  apart, they are not. It also never crosses different Au atom counts, which would require a reservoir term.
  <br><br>These integrals are now taken over the actual &plusmn;0.5&nbsp;V curves, not scaled up from the
  narrower window as an earlier draft warned against.</div>
</div></section>

<section id="gallery"><div class="finding">
  <div><div class="eyebrow">Structure gallery</div><h2>All thirty-three structures, coloured by coordination
  number</h2></div>
  <p class="lede">Each figure has three panels. On the left is a <b>plain outline with no atoms drawn</b>, so
  the shape of the defect is readable at a glance; in the middle, a top view coloured by coordination number;
  on the right, a side view. All three are generated from the geometry that was actually computed, not from a
  regenerated idealisation and not by hand, so the cartoon cannot drift from the real structure.
  <b>Click any figure to enlarge.</b></p>
  <div class="note" style="max-width:none"><b>What the marks mean.</b>
  <br>· <b>Colour</b> = coordination number (Au neighbours within 3.4&nbsp;&Aring;). The environment named in
  brackets is only what that coordination number <em>usually</em> means; coordination is not an identity.
  <br>· <b>Brightness and outline</b> = whether the atom belongs to the upper surface, by the same rule the
  region analysis uses: three or more higher neighbours within 2.35&nbsp;&Aring; counts as buried.
  <br>· <b>Green ring</b> = an added Au; <b>purple dashed circle</b> = a removed site. Drawn only where the
  structure family is itself defined by adding or removing atoms, so that a step's upper terrace or a vicinal
  step edge is never mislabelled as an adatom.
  <br>· <b>Blue &times;</b> = hcp-like registry (an atom directly two layers below); <b>brown open square</b> =
  transition band or domain wall (registry between fcc-like and hcp-like). The thresholds are 0.35 and 0.75 of
  the fcc in-plane offset. This is a <b>threshold-dependent display classification</b>, not the definition of
  a domain; it need not match the count from the different threshold used at build time, and the threshold
  should not be tuned to reproduce that count. The reference fcc offset comes from the <em>same-layer</em>
  nearest-neighbour spacing (R2's compressed layer measures 2.753&nbsp;&Aring;, not the 4.40&nbsp;&Aring; a
  cross-layer projection gives). These marks are drawn <b>only on structures whose normal is (111)</b>: the
  test uses the global z axis and the (111) interlayer spacing, which means nothing on
  Au(211)/(221)/(332)/(554). A few unexplainable marks once appeared there and have been withdrawn.
  <br>· <b>A&ndash;A&prime;</b> (and <b>B&ndash;B&prime;</b> where needed) = where the section was actually
  taken. The panel title says whether it is a real section or a projected envelope along the line of sight.
  The section is referenced to the <b>terrace level</b>, so protrusions sit above the baseline and depressions
  below it.
  <br>· <b>Heavy brown dashes</b> = the parent structure's original edge. Drawn only on structures modified
  from an existing one, and taken from the build-time atom mapping.
  <br>· For a finite feature (island, pit, point defect, composite) the atom side view draws only a
  <b>narrow band</b> around A&ndash;A&prime;. A whole-cell projection stacks the top-layer atoms in front of
  and behind a pit on top of it, and the pit looks filled in. The band half-width is a <b>true perpendicular
  distance</b>; it used to be computed as a fractional-coordinate difference times the other cell vector's
  length, which in a 60&deg; cell is narrower than the caption by 1/sin60, so a stated
  &plusmn;3.23&nbsp;&Aring; was really &plusmn;2.80&nbsp;&Aring;.
  <br><br><b>The outline follows atom connectivity, not a height partition, and its scale is the same-layer
  spacing.</b> A height partition is not the outline of a structure: a lower atom can win territory between
  two upper ones and cut a connected cluster in two. The outline is now the union of discs around the feature
  atoms, of radius 0.62&nbsp;&times;&nbsp;the same-layer nearest-neighbour spacing. The words <b>same
  layer</b> are what matters: the previous version took the nearest neighbour in projection, and an adatom is
  only a<sub>0</sub>/&radic;3&nbsp;=&nbsp;1.70&nbsp;&Aring; laterally from the terrace atoms beneath it, so on
  A3, C1 and the small pits more than half the atoms had such a cross-layer partner as their nearest and the
  median came out 1.70 instead of 2.94&nbsp;&Aring;. Discs 2.10&nbsp;&Aring; across cannot bridge a
  2.94&nbsp;&Aring; bond, so A3's three-atom cluster was drawn as three circles and Pit-7-compact's floor as
  twelve. Checked structure by structure: the number of disc components now <b>equals</b> the number of bonded
  clusters, every time.
  <br><br><b>Periodic copies are translations.</b> The distance field is computed once on the base cell and
  translated. Previously each displayed copy re-called a distance function whose image search reached only
  &plusmn;1 cell, so a copy two cells from the centre lost the feature entirely: the three copies of Step-8x2
  covered 53.18% / 51.11% / 0.40% of their area, which no periodic structure should do.
  <br><br><b>The section corresponds to A&ndash;A&prime; point by point.</b> The profile used to be rolled
  while the A and A&prime; labels stayed put, so the two came apart: Pit-19-8x8 samples terrace 16 &rarr; pit
  floor 75 &rarr; terrace 8 along its marked line, and the roll turned that into pit 38 &rarr; terrace 24
  &rarr; pit 37, which is exactly why a pit looked like a central bump. Centring a feature now shifts the
  whole display together, the origin A, the plan line, the profile and the atom band, rather than one curve.
  The atom side view also uses the along-line and across-line coordinates rather than projecting on the
  global x axis.
  <br><br><b>One periodic image per atom.</b> The along-line coordinate and the perpendicular distance used to
  be folded independently: the perpendicular one was wrapped into the nearest band while the along-line one
  kept the atom's original position. In a sheared cell the across vector has a component along the line, so
  that is not allowed. With <b>a</b>&nbsp;=&nbsp;(10,&nbsp;0) and <b>b</b>&nbsp;=&nbsp;(5,&nbsp;8.6603), the
  same atom written as <b>r</b> and as <b>r</b>+<b>b</b> came out at 9.5 and 4.5&nbsp;&Aring; along the line,
  half a period apart from nothing but a change of notation. The image <b>r</b>&minus;m<b>b</b> is now chosen
  from the perpendicular distance and both coordinates are read off it.
  <br><br><b>The drawing has a regression test, and the test itself is checked by a mutation test.</b>
  <code>scripts/check_gallery_drawing.py</code> asserts five things for all 33 structures: the disc outline
  has as many connected components as the feature atoms' bond graph (with the bond cutoff taken from the
  shortest Au&ndash;Au distance in the cell, <b>not</b> from the spacing function under test); every periodic
  copy is a translation of the base array with an identical level-set coverage; the profile equals the field
  re-sampled independently at A&nbsp;+&nbsp;s&middot;t&#770;; the band width and the side-view coordinates
  match an independent image search; and the view is unchanged by an in-plane lattice translation. It
  <b>exits non-zero</b> on failure. <code>scripts/check_drawing_mutations.py</code> puts each of the four bugs
  back one at a time and requires the regression to report it. All four are caught; the disc-radius mutation
  raises eight connectivity failures, among them A3 (1 bonded cluster, 3 circles) and Pit-7-compact (1 and
  12).
  <br><br><b>No "ion-accessible boundary" is drawn.</b> That region is set by SION and varies over terraces,
  islands and pits, which is one of this study's own results, so a horizontal line at a fixed height would
  contradict it. Only the electrolyte side is labelled.
  <br><br><b>Marks specific to certain families.</b>
  <br>· <b>The four vicinal faces</b> (Au(211)/(221)/(332)/(554)): the plan-view gradient is only a height
  distribution and cannot stand in for a step diagram, so the section is taken along the lattice direction of
  greatest relief and marks the local (111) terrace and the riser at the periodic seam <b>qualitatively</b>.
  The only number given is the <b>nominal angle between the macroscopic (hkl) and (111)</b>, computed from the
  plane normals as arccos[(h+k+l)/(&radic;3&middot;&radic;(h²+k²+l²))], and <b>not measured from this
  figure</b>. Terrace width and step height are <b>not</b> printed: those would need a chosen atom row and a
  stated measurement direction, which this figure does not do.
  <br>&nbsp;&nbsp;(A "tilt" used to be printed here, computed as ptp(profile)/(L&minus;grid step), the
  apparent slope of the height-partition field along one lattice direction over nearly a whole period. When
  the cut happens to be the steepest descent and the terrace is wide it comes close, as for Au554 at 5.6&deg;
  against a nominal 5.77&deg;, but it is not the same quantity and is not generally close: Au211 gave
  25.8&deg; against a nominal 19.47&deg;. Withdrawn.)
  <br>· <b>Kinks and the edge-detachment end state</b>: the parent's original straight edge is overlaid in
  brown dashes and the changed atoms circled. Kink-edge1/2's parent is the same strip without a kink,
  recovered by the build's own rule (four full rows of three atoms plus one atom alone in a fifth row, which
  is the kink); Step-8x2_edge-vacancy_plus_foot-adatom is diffed atom by atom against Step-8x2, giving 1 added
  and 1 removed.
  <br>· <b>C1</b>: the seven island atoms sit at exactly the <b>same height</b> as the step's upper terrace,
  so no "higher than most atoms" rule can see them. They come from the build-time atom mapping instead,
  diffing against the parent Step-8x4: 7 added, 0 removed. The section line is also pulled onto the row that
  passes through them.
  <br>· <b>C2</b>: its island and pit are not on the same lattice line, and a single straight line along a
  lattice vector can only clip the rim of one of them, which draws one pit as three shallow ones. It gets
  <b>two</b> sections, A&ndash;A&prime; through the island and B&ndash;B&prime; through the pit, each with its
  own atom band.
  <br><br><b>The outline is an aid to recognising the morphology; do not read boundary positions, widths or
  peak positions off it.</b> The morphology class comes from the frozen plan's structure family plus periodic
  connectivity, not from the area fraction of the high and low regions.</div>
  <div class="filters" id="filters"></div>
  <div class="gal bleed" id="gal"></div>
</div></section>

<section id="maps"><div class="finding">
  <div><div class="eyebrow">Anion maps</div><h2>The anion distribution, structure by structure</h2></div>
  <p class="lede">Four panels each: the per-column anion excess at the most positive sampled potential, how it
  changes between the two ends of the sampled window, the metal's own positive-charge change over the same
  potential step (&minus;&Delta;n<sub>e</sub>), and the coordination-number map the regions are cut from.
  Tiled periodically so the repeat is visible. <b>Click to enlarge.</b></p>
  <div class="note" style="max-width:none"><b>How to read the colour scales.</b> Each panel has <b>its
  own</b> range, so <b>amplitudes cannot be compared between panels</b>, only the spatial distribution within
  one. Single-signed quantities use a sequential scale (dark = large); quantities that cross zero use a
  diverging scale centred on zero. Units: &Gamma;<sub>&minus;</sub> and &Delta;&Gamma;<sub>&minus;</sub> are
  ions per projected area (&Aring;<sup>&minus;2</sup>); the metal panel shows the change in <b>positive</b>
  charge, &minus;&Delta;n<sub>e</sub> (e/&Aring;<sup>2</sup>, positive = electrons lost); coordination number
  is a discrete scale.</div>
  <div class="gal bleed" id="mapgrid"></div>
</div></section>

<section id="method"><div class="finding">
  <div><div class="eyebrow">Method</div><h2>Conventions, and what these numbers do not cover</h2></div>
  <div class="grid2">
    <div class="card pad"><h3>Potential and charge</h3><p class="small" style="margin-top:8px">
      U&nbsp;=&nbsp;&mu;<sub>0</sub>&nbsp;&minus;&nbsp;&mu;<sub>e</sub> with
      &mu;<sub>0</sub>&nbsp;=&nbsp;&minus;4.9071&nbsp;eV is an internal reference for this project. It is
      neither a potential against RHE nor each structure's own potential of zero charge.
      &sigma;&nbsp;=&nbsp;&minus;e&thinsp;(N<sub>e</sub>&nbsp;&minus;&nbsp;N<sub>e</sub><sup>0</sup>)&thinsp;/&thinsp;A<sub>proj</sub>.
      Every number comes from the &mu;<sub>e</sub> and N<sub>e</sub> the run actually converged to, never from
      the target, because the convergence criterion lets a run stop up to 10&nbsp;meV from it.</p></div>
    <div class="card pad"><h3>Ion density</h3><p class="small" style="margin-top:8px">
      Reconstructed from the converged potential and the ion-accessibility mask through the model's own
      constitutive relation, including the finite-size saturation term. No explicit ion is added and there is
      no chlorine in the system: this is the non-specific anion response of a continuum electrolyte, not a
      chemisorption preference of Cl.</p></div>
    <div class="card pad"><h3>Integration windows</h3><p class="small" style="margin-top:8px">
      The absolute window takes everything below z&nbsp;&lt;&nbsp;31&nbsp;&Aring;, so the columns sum exactly
      to the whole cell. The boundary-relative window takes 10.4&nbsp;&Aring; upward from each column's own
      accessibility boundary, so a raised island and the terrace beside it keep the same fraction of the decay
      tail. Both are reported, and the findings above hold under either.</p></div>
    <div class="card pad"><h3>How the regions are cut (a definition that was corrected)</h3>
      <p class="small" style="margin-top:8px">
      Each column is assigned to the nearest <b>un-buried</b> surface atom and classified by that atom's
      coordination number. The first version treated every atom within 3&nbsp;&Aring; of the top as a
      candidate, but the (111) interlayer spacing is only 2.4&nbsp;&Aring;, so the second layer qualified too
      and it sits directly below the hollow sites: on a perfectly flat T-4x4, 48% of the columns were assigned
      to a CN&nbsp;12 second-layer atom, splitting one flat terrace into "terrace" and "sub-surface". With the
      un-buried rule a flat face is 100% terrace again. The conclusions above hold before and after the fix
      (before: 36/41, median &minus;1.7%).</p></div>
    <div class="card pad"><h3>How the scope is handled</h3><p class="small" style="margin-top:8px">
      Every cross-structure number uses all five potentials (&minus;5.4071 / &minus;5.1071 / &minus;4.9071 /
      &minus;4.7071 / &minus;4.4071&nbsp;eV). The &plusmn;0.5&nbsp;V extension completed on 2026-10-01, 214 of
      214 tasks with no failures, which is why the scope was widened; while it was incomplete the page used
      the three base potentials only, so that a geometry sampled to &plusmn;0.5&nbsp;V was never compared
      against one sampled to &plusmn;0.2&nbsp;V. Widening was checked first and changes very little: the
      spread of the zero-charge point is identical at 254.3&nbsp;mV and the median secant capacitance moves
      from 11.84 to 11.59&nbsp;&micro;C/cm&middot;V&#8315;&sup1;. It improves coverage, because <b>both outer
      points exist for all 107 geometries</b> while the middle three do not (105 / 85 / 101 of 107): the
      zero-charge point is now obtained for 107 of 107 rather than 104. Quantities that need more points are
      still computed only where those points exist, and each chart states how many geometries contributed.
      <br><br>One state is a partial loss.
      <span class="mono">C2-island+pit__ideal__mu-5.4071</span> was killed at its walltime while writing the
      bound-charge grid: its SCF had closed, so N<sub>e</sub> and &mu;<sub>e</sub> are valid and it is used in
      the charging and region analyses, but it has no metal-to-ion transmission pair, which is why 106 of 107
      geometries appear in that chart. The snapshot at the top gives the counts and build time this page used.</p></div>
    <div class="card pad"><h3>Revision log</h3><p class="small" style="margin-top:8px">
      This page has been through several rounds of external review and was changed accordingly. The body keeps
      only the currently correct explanation; what changed and why is in the expandable
      <a href="#revlog">revision log</a> at the foot. The raw data and the scripts are traceable in the
      repository.</p></div>
    <div class="card pad"><h3>Known limits</h3><p class="small" style="margin-top:8px">
      At production settings the forces carry an egg-box error of about 0.02&nbsp;eV/&Aring; per atom, so no
      conclusion here rests on a small force difference. The perturbed and collectively deformed
      configurations are a designed sampling, not a thermodynamic ensemble; they are used to give a
      sensitivity range and are never averaged over. The path images are sparse and are not a minimum energy
      path.</p></div>
  </div>
</div></section>

<section id="revlog"><div class="finding">
  <div><div class="eyebrow">Revision log</div><h2>What was changed, and why</h2></div>
  <p class="lede">The body keeps only the currently correct explanation. This records the substantive changes
  made after each round of external review, so they can be traced; every one was checked against the data or
  a unit test before it was made.</p>
  <details open><summary style="cursor:pointer;font-weight:500;padding:10px 0">Latest round (the
  &plusmn;0.5&nbsp;V window)</summary>
  <ul class="small" style="line-height:1.9;max-width:62em">
    <li><b>The scope was widened from &plusmn;0.2 to &plusmn;0.5&nbsp;V.</b> The extension completed on
      2026-10-01, 214 of 214 tasks with no failures, 325 node-h. What it changed was checked before the
      switch: the spread of the zero-charge point is identical at 254.3&nbsp;mV, the median secant capacitance
      moves 11.84&nbsp;&rarr;&nbsp;11.59, its spread 14.3&nbsp;&rarr;&nbsp;15.4%, and coverage improves from
      104 to 107 of 107 geometries because both outer points exist everywhere while the middle three do
      not.</li>
    <li><b>Two results the narrow window could not give.</b> &sigma;(U) is straight across the full
      1&nbsp;V: over the 79 geometries with all five points the worst departure from a line is
      0.081&nbsp;&micro;C/cm² against a range of 7.32, about 1%, so one secant capacitance is well defined
      over the window. That is the total surface charge of a fixed geometry only; nothing is implied about
      forces, local potentials, bound charge or ion concentrations, and the multi-potential comparisons in
      the data stay in place. The grand-potential numbers also grew with the window, to a cross-morphology maximum of
      89&nbsp;meV and a within-morphology maximum of 138&nbsp;meV (29 and 27 at &plusmn;0.2&nbsp;V). An
      earlier draft of this entry said the comparison made the "not enough to re-rank" conclusion firmer;
      <b>that was wrong</b> and was corrected on review. Both are potential-induced changes for different
      configuration pairs, not full relative grand potentials and not error bars, so their ratio bears on
      nothing about the ranking. Finding 5 now states only what the numbers show.</li>
    <li><b>The under-coordination comparison, now at about +0.5&nbsp;V.</b> This is a difference of
      <em>region enrichment factors</em> K at each geometry's most positive sampled potential, now
      +0.5&nbsp;V rather than +0.2&nbsp;V; it is not the response ratio
      &Delta;K<sub>low</sub>/&Delta;K<sub>terrace</sub> and should not be read as one. At +0.5&nbsp;V the
      under-coordinated regions enrich less than the terrace in <span data-n="uc_main"></span>, against a
      median 1.4% at +0.2&nbsp;V. The exceptions drop from six to five, and the largest of them, the dilute
      8&times;8 island at +3.6%, is no longer one.</li>
    <li><b>One state is a partial loss, recorded rather than hidden.</b>
      <span class="mono">C2-island+pit__ideal__mu-5.4071</span> was killed at its walltime while writing the
      bound-charge grid (534 of 713&nbsp;MB; 29&thinsp;370&thinsp;510 of 39&thinsp;200&thinsp;000 values), and
      the ion grid was never written. Its SCF had closed, so N<sub>e</sub> and &mu;<sub>e</sub> are valid and
      it is used in the charging and region analyses; it has no transmission pair, which is why that chart
      covers 106 of 107 geometries. Re-running it would cost about 20 node-h for one transmission point, so it
      was not re-run.</li>
  </ul></details>
  <details><summary style="cursor:pointer;font-weight:500;padding:10px 0">Previous round (sampling, and
  figure legibility)</summary>
  <ul class="small" style="line-height:1.9;max-width:62em">
    <li><b>An opening figure of the sampling tree.</b> The body used to start straight at Finding 1, leaving
      no way to know where the 107 geometries came from, and inviting the reading that they are 107 separately
      designed defects. The figure shows one thing only: how 33 hand-built structures become 107 geometries
      through relaxation, random displacement, collective deformation and path images. It also names the five
      potentials explicitly, because "constant potential" reads as a single setting otherwise.</li>
    <li><b>Every figure is now one per row and wider than the text column.</b> The structure gallery and the
      maps were laid out two to a row at about 560&nbsp;px each; a 9&nbsp;pt label on a 12-inch canvas then
      renders at about 7 CSS px, which is unreadable. One per row at up to 1780&nbsp;px puts the same label
      near 18&nbsp;px. The maps went from 1&times;4 to 2&times;2 (canvas 17.2&nbsp;&rarr;&nbsp;10 inches) with
      their type enlarged as well, and the SVG charts now scale with their container.</li>
    <li><b>The under-coordinated versus terrace count is computed, not hardcoded.</b> The body used to state
      "37 of 41 geometries, four exceptions, all vacancy rims, exceeding by only 0.1&ndash;0.3%" as a literal,
      with no record of its definition, and it does not reproduce from <code>regions.json</code>. It is now
      recomputed at build time with the definition written out: <span data-n="uc_main"></span>. The exceptions
      are <span data-n="uc_exc"></span>The old wording missed two of them and understated the largest
      excess.</li>
    <li><b>The "terrace tilt" on the vicinal sections was withdrawn.</b> What it printed was
      ptp(profile)/(L&minus;grid step), the apparent slope of the height-partition field along one lattice
      direction over nearly a whole period. Au554 giving 5.6&deg; against a nominal 5.77&deg; was a
      coincidence; Au211 gave 25.8&deg; against a nominal 19.47&deg;. Only the nominal angle from the plane
      normals is shown now, labelled as not measured from the figure.</li>
  </ul></details>
  <details><summary style="cursor:pointer;font-weight:500;padding:10px 0">Second round (what the figures
  assert)</summary>
  <ul class="small" style="line-height:1.9;max-width:62em">
    <li><b>The drawn "ion-accessible boundary" was removed.</b> A horizontal dashed line used to be drawn at a
      fixed 4.2&nbsp;&Aring; above the highest atom in the atom side view and at a fixed 3.2&nbsp;&Aring; in
      the outline, two different offsets, neither read from SION. The accessible region varies over terraces,
      islands and pits, which is one of this study's results, so representing it by a line at fixed height
      contradicts the work. Only the electrolyte side is labelled now.</li>
    <li><b>The morphology class is no longer decided by area fraction.</b> The old rule called anything with
      30&ndash;70% high area a step, so Pit-19-8x8, a finite shallow pit, was labelled "two-level terrace
      (step)". It now comes from the frozen plan's structure family plus periodic connectivity: a strip step
      spans the cell, an island or a pit closes on itself.</li>
    <li><b>A section no longer pretends to be one.</b> When a cut through the feature centre was unavailable
      the code fell back to a projected envelope. The title now says honestly which of the two it is and the
      plan view shows where A&ndash;A&prime; actually lies; the solid body below the section is marked as
      substrate indication only.</li>
    <li><b>Information that height cannot show was added.</b> hcp registry marks (blue &times;) were
      introduced, since R1-hcp-terminated and a flat fcc face, and A1-hcp and A1-fcc, had been indistinguishable
      in the figures; so were the "added Au" and "removed site" marks, gated by structure family so that a
      step's upper terrace is never mislabelled as an adatom.</li>
    <li><b>Coordination counts are for the upper surface only.</b> A whole-slab count includes the back face,
      so a four-layer flat slab read "terrace 32" when only 16 atoms face the electrolyte.</li>
    <li>The Fourier paragraph had its direction reversed (it falls as the wavelength <em>shortens</em>);
      "different dimensions" became "different physical quantities"; the per-column correlation is no longer
      described as locating peaks; "the offset within a cell cancels exactly" became "isolates the morphology
      effect better, without guaranteeing full cancellation"; and "all three base potentials are complete for
      every geometry" was corrected, since some reference structures only ever had one or two.</li>
  </ul></details>
  <details><summary style="cursor:pointer;font-weight:500;padding:10px 0">First round (science and
  code)</summary>
  <ul class="small" style="line-height:1.9;max-width:62em">
    <li><b>The grand-potential integral had the wrong sign.</b> Since
      &part;&Omega;/&part;&mu;&nbsp;=&nbsp;&minus;N and U&nbsp;=&nbsp;&mu;₀&minus;&mu;, U&nbsp;=&nbsp;+w is
      the lower &mu; and the integral must not be negated again. Unit test: with
      N<sub>A</sub>&minus;N<sub>B</sub>&nbsp;=&nbsp;1 over a &plusmn;0.2&nbsp;V window the answer is
      +0.4&nbsp;eV; the old code gave &minus;0.4&nbsp;eV, which named the opposite geometry as the one a
      positive potential stabilises.</li>
    <li><b>Finding 3's mechanism was backwards, and it confused global with local.</b> From
      &sigma;&nbsp;=&nbsp;C(U&minus;U<sub>pzc</sub>), a negative shift of the zero-charge point means
      <em>more</em> positive at the same potential, not less; and the zero-charge point here is a scalar for
      the whole periodic cell and cannot be assigned to one region. Retracted and replaced with what the data
      do support.</li>
    <li><b>The metal&ndash;ion phase comparison had a sign error.</b> Charging positive removes electrons and
      adds anions, so comparing the two directly put a spurious &pi; phase between patterns that are in
      register. It now compares &minus;&Delta;n<sub>e</sub> against
      &Delta;&Gamma;<sub>&minus;</sub>.</li>
    <li><b>The statistical scope was being contaminated by the unfinished extension.</b> The charging analysis
      had been ingesting &plusmn;0.5&nbsp;V states as they finished, 68 of the 107 geometries already carried
      five points, while the page still said three potentials. Cross-structure numbers now use the three base
      potentials only, with the extension drawn separately.</li>
    <li><b>"Eight tenths of it comes from the zero-charge point" does not hold.</b> A ratio of ranges is a
      scale comparison, not a variance decomposition: the expansion has a cross term and &delta;C and
      &delta;Z are correlated. "The capacitance hardly changes" went too, since A3 is 15% above the flat
      member of its cell.</li>
    <li>Withdrawn: that the potential does not re-rank stability (no reference-potential baseline), and that
      &plusmn;0.5&nbsp;V would scale these values by 2.5. The half-amplitude wavelength and the contrast ratio
      were downgraded to definition-dependent indicators, and the short-wavelength damping is no longer
      attributed solely to the 4&nbsp;&Aring; cavity.</li>
    <li><b>The region partition split a flat terrace in two.</b> Columns were assigned to the nearest atom
      within 3&nbsp;&Aring; of the top, but the interlayer spacing is only 2.4&nbsp;&Aring; and the second
      layer sits directly under the hollow sites: a perfectly flat T-4x4 had 48% of its columns labelled
      CN&nbsp;12. Assigning to un-buried atoms restores a flat face to 100% terrace and moved the count from
      36/41 to 37/41.</li>
    <li>The &sigma; chart had its y axis hardcoded to &plusmn;7&nbsp;&micro;C/cm² while the data reach 7.32,
      so the curves were clipped; it is derived from the data now.</li>
  </ul></details>
</div></section>

</div>
<footer><div class="wrap">
  Built from <span class="mono">dataset_v1/states.json</span> and the per-column field reduction. <span id="stamp"></span>
</div></footer>

<dialog class="lb" id="lb">
  <div class="bar"><span id="lbcap"></span><button type="button" id="lbclose">Close (Esc)</button></div>
  <img id="lbimg" alt="">
</dialog>

<script type="application/json" id="D">__DATA__</script>
<script>
const D = JSON.parse(document.getElementById("D").textContent);
const FAM = ["flat Au(111)","point defect","reconstruction-related","strip step","vicinal step face",
             "kink / edge rearrangement","single-layer island","single-layer pit","composite"];
const FZH = {"flat Au(111)":"flat Au(111)","point defect":"point defect",
  "reconstruction-related":"reconstruction-related","strip step":"strip step",
  "vicinal step face":"vicinal step face","kink / edge rearrangement":"kink / edge rearrangement",
  "single-layer island":"single-layer island","single-layer pit":"single-layer pit","composite":"composite"};
const CZH = {"kink/adatom":"kink / adatom","edge/rim":"edge / rim","terrace":"terrace",
  "sub-surface/foot":"step foot / sub-surface"};
const FC = {"flat Au(111)":"#4c78a8","point defect":"#e08a2e","reconstruction-related":"#54a24b",
  "strip step":"#9c6bb0","vicinal step face":"#cf4f45","kink / edge rearrangement":"#3fa0a8",
  "single-layer island":"#c8a33a","single-layer pit":"#8a6a4e","composite":"#8d8f94"};
const CSS = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const fmt = (x,n=2) => (x===null||x===undefined||!isFinite(x)) ? "–" : x.toFixed(n);
const SVGT = new Set(["g","rect","circle","line","path","text","polyline","polygon","svg","tspan"]);
const el = (t,a={},...k) => {const e=document.createElementNS(SVGT.has(t)?"http://www.w3.org/2000/svg":"http://www.w3.org/1999/xhtml",t);
  for(const [p,v] of Object.entries(a)) if(v!==null&&v!==undefined) e.setAttribute(p,v);
  for(const c of k.flat()) if(c!==null&&c!==undefined) e.append(c.nodeType?c:document.createTextNode(c)); return e;};
function median(a){const b=[...a].sort((x,y)=>x-y);const n=b.length;return n%2?b[(n-1)/2]:(b[n/2-1]+b[n/2])/2;}
/* every number quoted in the prose is written here from the same data object the charts read, so the text can
   never drift from the figures the way a hardcoded "291 states at three potentials" did. */
function fill(key,text){for(const e of document.querySelectorAll(`[data-n="${key}"]`)) e.textContent=text;}
function frame(w,h,m){return {s:el("svg",{viewBox:`0 0 ${w} ${h}`,role:"img",
  style:`width:100%;min-width:${w}px;height:auto`}),w,h,m};}

/* ---- lightbox ---- */
const LB=document.getElementById("lb"), LBI=document.getElementById("lbimg"), LBC=document.getElementById("lbcap");
document.getElementById("lbclose").addEventListener("click",()=>LB.close());
LB.addEventListener("click",e=>{if(e.target===LB||e.target===LBI) LB.close();});
function zoomable(src,cap,alt){
  const b=el("button",{type:"button",class:"zoom","aria-label":`enlarge ${cap}`},
    el("img",{src:src,alt:alt,loading:"lazy"}));
  b.addEventListener("click",()=>{LBI.src=src;LBI.alt=alt;LBC.textContent=cap;
    if(LB.showModal) LB.showModal(); else LB.setAttribute("open","");});
  return b;}

/* ---- kicker + data snapshot, all generated from the same object the charts read ---- */
const SN=D.snapshot||{};
const NG=Object.keys(D.geometries).length, NS=Object.keys(D.structures).length;
const NPT=Object.values(D.geometries).reduce((a,g)=>a+g.pts.length,0);
document.getElementById("kicker").append(
  ...[[NS,"structures"],[NG,"distinct geometries"],[NPT,"electronic states"],["5","potentials"],["\u00b10.5 V","window, complete"]]
     .map(([v,l])=>el("span",{class:"chip"},`${v} ${l}`)));
document.getElementById("scope").append(
  `All five potentials, \u22125.4071 / \u22125.1071 / \u22124.9071 / \u22124.7071 / \u22124.4071 eV: `
  +`${NPT} states, ${NG} geometries, ${NS} structures. The \u00b10.5 V extension completed on 2026-10-01 `
  +`(214 of 214 tasks, no failures), so every number on this page now uses the full \u00b10.5 V window. `
  +`The region and spatial analysis covers ${SN.n_region_states||"?"} states `
  +`(${SN.n_region_structures||"?"} structures); ${Object.keys(D.transmission).length} of ${NG} geometries `
  +`admit a metal-to-ion comparison.`,
  el("br"),
  el("span",{class:"mono",style:"font-size:12px"},
    `built ${SN.built||"?"}${SN.commit?"  \u00b7  commit "+SN.commit:""}`));
document.getElementById("stamp").textContent = ` \u00b7 snapshot ${SN.built||""}${SN.commit?" / "+SN.commit:""}`;

/* ---- Finding 1 stats ---- */
const zs=Object.values(D.geometries).filter(g=>g.z!==null).map(g=>g.z);
const cs=Object.values(D.geometries).filter(g=>g.C!==null).map(g=>g.C);
const dAll=D.decomposition["all geometries"];
const PZCSPREAD=1000*(Math.max(...zs)-Math.min(...zs)), CSPREAD=100*(Math.max(...cs)-Math.min(...cs))/median(cs);
document.getElementById("s1").append(...[
  [`${PZCSPREAD.toFixed(0)} mV`,`zero-charge point spread (${zs.length} geometries bracketed by their samples)`],
  [`${CSPREAD.toFixed(0)}%`,"secant-capacitance spread, same geometries"],
  [`${median(cs).toFixed(1)} µF/cm²`,"median secant capacitance"],
  [`${dAll["U=+0.2V"].scale_ratio.toFixed(1)}\u00d7`,"ratio of the two scales at U = +0.2 V (not a share)"],
].map(([v,l])=>el("div",{class:"stat"},el("div",{class:"v"},v),el("div",{class:"l"},l))));
fill("pzcspread",`${PZCSPREAD.toFixed(0)} mV`);
fill("cspread",`${CSPREAD.toFixed(0)}%`);
fill("ngeom",String(NG));
{const p=D.pairs.map(x=>Math.abs(x.d_relative_Omega_eV));
 fill("dom",`${(1000*Math.max(...p)).toFixed(0)} meV`);
 fill("domsame",`${(1000*D.same_structure_pair_scale.max).toFixed(0)} meV`);}
{const T=Object.values(D.transmission);
 const rs=T.map(v=>v.r).filter(x=>x!==null&&x!==undefined);
 const neg=rs.filter(x=>x<0).length;
 fill("corr",`${neg} of ${rs.length} geometries are negative, median r = ${median(rs).toFixed(2)}`);
 const hs=T.map(v=>v.tf&&v.tf.half).filter(x=>x);
 fill("half",`half-amplitude wavelength ${median(hs).toFixed(1)} Å`);
 fill("ratio",`contrast ratio ${(100*median(T.map(v=>v.ratio).filter(x=>x))).toFixed(0)}%`);
}
{const U=D.undercoord;
 if(U){
  const txt=`${U.n_lower} of ${U.n} geometries carrying both region types, median `
    +`${Math.abs(U.median_pct).toFixed(1)}% lower, at most ${Math.abs(U.min_pct).toFixed(1)}% lower `
    +`(${U.min_structure})`;
  fill("uc_main",txt);
  const v=U.exceptions.filter(e=>/^V[0-9]/.test(e.structure)).length;
  fill("uc_exc",`${U.exceptions.length} in all: `
    +U.exceptions.map(e=>`${e.structure} (${e.config}, ${e.pct.toFixed(1)}% higher)`).join(", ")
    +`. ${v} are vacancy structures; the others are the small island and pit in the dilute 8\u00d78 cell. `);}}
fill("ngeom2",String(NG));

/* ---- the opening figure ---- */
(function(){
  const n=document.getElementById("originfig"); if(!n) return;
  const src=`${D.imgbase}gallery/_lineage.png`;
  const cap="33 hand-built structures become 107 geometries through relaxation, random displacement, collective deformation and path images; each geometry then carries several potentials.";
  n.append(zoomable(src,cap,"branching diagram from the hand-built geometries to the 107 computed ones"),
           el("div",{class:"caption"},cap), el("div",{class:"hint"},"click to enlarge"));
})();

/* ---- sigma(U) ---- */
(function(){
  const F=frame(880,420,{l:56,r:14,t:14,b:42}), m=F.m, {w,h}=F;
  // derive the axis limits from the data: hardcoding +-7 clipped the curves once the +-0.5 V states arrived
  const allS=Object.values(D.geometries).flatMap(g=>g.pts.map(p=>p[1]));
  const allU=Object.values(D.geometries).flatMap(g=>g.pts.map(p=>p[0]));
  const SY=Math.ceil(Math.max(...allS.map(Math.abs))+0.4), SX=Math.max(...allU.map(Math.abs))+0.05;
  const X=x=>m.l+(x+SX)/(2*SX)*(w-m.l-m.r), Y=y=>h-m.b-(y+SY)/(2*SY)*(h-m.t-m.b);
  for(const t of [-0.5,-0.25,0,0.25,0.5]){
    F.s.append(el("line",{x1:X(t),x2:X(t),y1:m.t,y2:h-m.b,stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:X(t),y:h-m.b+16,"text-anchor":"middle","font-size":11,fill:CSS("--muted")},String(t)));}
  for(const t of [-6,-4,-2,0,2,4,6]){
    F.s.append(el("line",{x1:m.l,x2:w-m.r,y1:Y(t),y2:Y(t),stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:m.l-8,y:Y(t)+3.8,"text-anchor":"end","font-size":11,fill:CSS("--muted")},String(t)));}
  F.s.append(el("line",{x1:X(-SX),x2:X(SX),y1:Y(0),y2:Y(0),stroke:CSS("--line-2"),"stroke-width":1.4}));
  // all five potentials are in scope since the extension completed, so there is no faint overlay any more
  // and no gold scope boundary at +-0.2 V; the sampled potentials themselves are marked instead
  for(const g of Object.values(D.geometries))
    F.s.append(el("polyline",{points:g.pts.map(p=>`${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join(" "),
      fill:"none",stroke:FC[g.f]||"#888","stroke-width":1.2,"stroke-opacity":.7,"stroke-linecap":"round"}));
  for(const u of [-0.5,-0.2,0,0.2,0.5]) F.s.append(el("line",{x1:X(u),x2:X(u),y1:m.t,y2:h-m.b,
    stroke:CSS("--gold"),"stroke-width":1,"stroke-dasharray":"4 3","stroke-opacity":.45}));
  F.s.append(el("text",{x:X(0.5)-5,y:m.t+13,"text-anchor":"end","font-size":10,fill:CSS("--gold")},
    "the five sampled potentials"));
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-5,"text-anchor":"middle","font-size":11,fill:CSS("--ink-2")},
    "U = μ₀ − μₑ   (V)"));
  F.s.append(el("text",{x:14,y:(m.t+h-m.b)/2,"text-anchor":"middle","font-size":11,fill:CSS("--ink-2"),
    transform:`rotate(-90 14 ${(m.t+h-m.b)/2})`},"σ   (µC/cm²)"));
  document.getElementById("c_sigma").append(F.s);
  document.getElementById("leg1").append(...FAM.map(f=>el("span",{},el("i",{class:"sw",style:`background:${FC[f]}`}),FZH[f])));
})();

/* ---- decomposition table ---- */
(function(){
  const t=document.getElementById("t_decomp");
  t.append(el("thead",{},el("tr",{},...["Set","n","U0 spread (mV)","C spread (%)","\u03c3 scale: U0 term","C term","ratio"]
    .map(h=>el("th",{},h)))));
  const tb=el("tbody");
  const zh=k=>k==="all geometries"?"all geometries":k.replace(/^within cell A=(\d+) A\^2 \((.*)\)$/,"within one cell, A=$1 Å² ($2)");
  for(const [k,v] of Object.entries(D.decomposition)){
    const w=v["U=+0.2V"];
    tb.append(el("tr",{},el("td",{},zh(k)),el("td",{class:"num"},String(v.n)),
      el("td",{class:"num"},fmt(v.U_pzc_spread_mV,0)),el("td",{class:"num"},fmt(v.C_spread_pct,1)),
      el("td",{class:"num"},fmt(w.from_pzc_shift_uC_per_cm2,2)),el("td",{class:"num"},fmt(w.from_capacitance_uC_per_cm2,2)),
      el("td",{class:"num"},fmt(w.ratio_pzc_over_capacitance,1)+"\u00d7")));}
  t.append(tb);
})();

/* ---- dU_pzc bars ---- */
(function(){
  const rows=Object.entries(D.structures).filter(([k,v])=>v.dz!==null&&v.dz!==undefined).sort((a,b)=>a[1].dz-b[1].dz);
  if(!rows.length) return;
  const H=31,w=880,m={l:250,r:66,t:28,b:36},h=m.t+m.b+rows.length*H;
  const lo=Math.min(-100,...rows.map(r=>r[1].dz)), hi=Math.max(30,...rows.map(r=>r[1].dz));
  const X=x=>m.l+(x-lo)/(hi-lo)*(w-m.l-m.r);
  const F=frame(w,h,m);
  for(const t of [-100,-75,-50,-25,0,25]){
    F.s.append(el("line",{x1:X(t),x2:X(t),y1:m.t-6,y2:h-m.b,stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:X(t),y:h-m.b+15,"text-anchor":"middle","font-size":11,fill:CSS("--muted")},String(t)));}
  F.s.append(el("line",{x1:X(0),x2:X(0),y1:m.t-6,y2:h-m.b,stroke:CSS("--line-2"),"stroke-width":1.4}));
  rows.forEach(([k,v],i)=>{const y=m.t+i*H+H/2;
    F.s.append(el("rect",{x:Math.min(X(0),X(v.dz)),y:y-7,width:Math.abs(X(v.dz)-X(0)),height:14,rx:2,
      fill:FC[v.f]||"#888","fill-opacity":.85}));
    F.s.append(el("text",{x:m.l-10,y:y+4,"text-anchor":"end","font-size":11.5,fill:CSS("--ink")},k));
    if(v.dC!==null&&v.dC!==undefined){
      F.s.append(el("circle",{cx:w-m.r+16,cy:y,r:4.4,fill:CSS("--surface"),stroke:CSS("--gold"),"stroke-width":1.6}));
      F.s.append(el("text",{x:w-m.r+25,y:y+3.6,"font-size":10,fill:CSS("--gold")},(v.dC>0?"+":"")+v.dC.toFixed(0)+"%"));}});
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-3,"text-anchor":"middle","font-size":11,fill:CSS("--ink-2")},
    "shift of the zero-charge point vs the flat terrace in the same cell (mV)"));
  F.s.append(el("text",{x:w-m.r+12,y:m.t-13,"font-size":10,fill:CSS("--gold")},"capacitance change"));
  document.getElementById("c_dpzc").append(F.s);
})();

/* ---- region K ---- */
(function(){
  const CLS=[["kink/adatom",CSS("--red")],["edge/rim",CSS("--amber")],["terrace",CSS("--terrace")],
             ["sub-surface/foot",CSS("--sub")]];
  const rows=Object.entries(D.regions).sort((a,b)=>a[1].regions.all.K-b[1].regions.all.K);
  if(!rows.length) return;
  const RH=48,w=880,m={l:210,r:100,t:22,b:40},h=m.t+m.b+rows.length*RH;
  const vals=rows.flatMap(([k,v])=>Object.values(v.regions).map(r=>r.K));
  const lo=Math.min(...vals)-0.008, hi=Math.max(...vals)+0.008;
  const X=x=>m.l+(x-lo)/(hi-lo)*(w-m.l-m.r);
  const F=frame(w,h,m);
  for(let t=Math.ceil(lo*20)/20;t<=hi;t+=0.05){
    F.s.append(el("line",{x1:X(t),x2:X(t),y1:m.t-4,y2:h-m.b,stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:X(t),y:h-m.b+15,"text-anchor":"middle","font-size":11,fill:CSS("--muted")},t.toFixed(2)));}
  // the axis starts near the data, not at zero, so mark the bulk baseline explicitly rather than let the bar
  // lengths imply a ratio they do not carry
  if(lo<1&&hi>1){
    F.s.append(el("line",{x1:X(1),x2:X(1),y1:m.t-10,y2:h-m.b,stroke:CSS("--teal"),"stroke-width":1.6}));
    F.s.append(el("text",{x:X(1)+5,y:m.t-2,"font-size":10.5,fill:CSS("--teal")},"bulk concentration K = 1"));}
  else{
    F.s.append(el("text",{x:m.l,y:m.t-2,"font-size":10.5,fill:CSS("--teal")},
      `note: the axis starts at ${lo.toFixed(2)}, so the bulk baseline K = 1 is off-scale`));}
  rows.forEach(([k,v],i)=>{
    const y0=m.t+i*RH, present=CLS.filter(([c])=>v.regions[c]);
    const bh=Math.min(8,(RH-18)/Math.max(present.length,1));
    F.s.append(el("text",{x:m.l-10,y:y0+RH/2-3,"text-anchor":"end","font-size":11.5,fill:CSS("--ink")},k));
    F.s.append(el("text",{x:m.l-10,y:y0+RH/2+11,"text-anchor":"end","font-size":9.5,fill:CSS("--muted")},
      "U = "+(v.U>0?"+":"")+v.U.toFixed(2)+" V"));
    present.forEach(([c,col],j)=>{
      const r=v.regions[c], y=y0+9+j*(bh+1.8);
      F.s.append(el("rect",{x:X(lo),y:y,width:Math.max(0,X(r.K)-X(lo)),height:bh,rx:1.5,fill:col,"fill-opacity":.9}));
      F.s.append(el("text",{x:X(r.K)+5,y:y+bh-0.4,"font-size":9,fill:CSS("--ink-2")},
        r.K.toFixed(3)+"  (area "+(100*r.a).toFixed(0)+"%)"));});
    F.s.append(el("line",{x1:X(v.regions.all.K),x2:X(v.regions.all.K),y1:y0+5,y2:y0+RH-7,stroke:CSS("--ink"),
      "stroke-width":1.2,"stroke-dasharray":"3 2"}));});
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-3,"text-anchor":"middle","font-size":11,fill:CSS("--ink-2")},
    "K\u03a9 = anion concentration in that region, relative to bulk"));
  document.getElementById("c_regions").append(F.s);
  document.getElementById("leg3").append(...CLS.map(([c,col])=>el("span",{},el("i",{class:"sw",style:`background:${col}`}),CZH[c])),
    el("span",{},el("i",{class:"sw",style:`background:${CSS("--ink")}`}),"whole-cell mean"));
})();

/* ---- transfer function ---- */
(function(){
  const curves=Object.values(D.transmission).filter(v=>v.tf&&v.tf.b.length);
  if(!curves.length) return;
  const halves=curves.map(v=>v.tf.half).filter(x=>x);
  document.getElementById("s4").append(...[
    [`${median(halves).toFixed(1)} Å`,"wavelength at half amplitude (median over geometries)"],
    [`${curves.length}`,"geometries with both end potentials, included here"],
    ["2.94 Å","Au\u2013Au nearest-neighbour spacing, as a scale reference"],
    [`${(100*median(Object.values(D.transmission).map(v=>v.ratio).filter(x=>x))).toFixed(0)}%`,
     "fraction of the metal's lateral contrast surviving into the anion map"],
  ].map(([v,l])=>el("div",{class:"stat"},el("div",{class:"v"},v),el("div",{class:"l"},l))));
  const bins=new Map();
  for(const v of curves) for(const [lam,T] of v.tf.b){
    const key=Math.round(Math.log(lam)*2.2); if(!bins.has(key)) bins.set(key,[]); bins.get(key).push(T);}
  const pooled=[...bins.entries()].map(([k,a])=>{
    const s=[...a].sort((p,q)=>p-q), q=f=>s[Math.min(s.length-1,Math.floor(f*s.length))];
    return {lam:Math.exp(k/2.2),T:median(s),lo:q(.25),hi:q(.75)};}).sort((a,b)=>a.lam-b.lam);
  const w=880,h=400,m={l:60,r:16,t:16,b:44}, F=frame(w,h,m);
  const l10=Math.log10, x0=l10(0.7), x1=l10(24);
  const X=x=>m.l+(l10(x)-x0)/(x1-x0)*(w-m.l-m.r);
  const Y=y=>h-m.b-(l10(Math.max(y,2e-3))-l10(2e-3))/(0-l10(2e-3))*(h-m.t-m.b);
  for(const t of [1,2,3,5,10,20]){
    F.s.append(el("line",{x1:X(t),x2:X(t),y1:m.t,y2:h-m.b,stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:X(t),y:h-m.b+16,"text-anchor":"middle","font-size":11,fill:CSS("--muted")},String(t)));}
  for(const t of [0.002,0.01,0.05,0.2,1]){
    F.s.append(el("line",{x1:m.l,x2:w-m.r,y1:Y(t),y2:Y(t),stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:m.l-8,y:Y(t)+3.8,"text-anchor":"end","font-size":11,fill:CSS("--muted")},t>=1?"1":String(t)));}
  F.s.append(el("polygon",{points:pooled.map(p=>`${X(p.lam).toFixed(1)},${Y(p.hi).toFixed(1)}`).join(" ")+" "+
    pooled.slice().reverse().map(p=>`${X(p.lam).toFixed(1)},${Y(p.lo).toFixed(1)}`).join(" "),
    fill:CSS("--teal"),"fill-opacity":.16}));
  F.s.append(el("polyline",{points:pooled.map(p=>`${X(p.lam).toFixed(1)},${Y(p.T).toFixed(1)}`).join(" "),
    fill:"none",stroke:CSS("--teal"),"stroke-width":2.2,"stroke-linejoin":"round"}));
  for(const p of pooled) F.s.append(el("circle",{cx:X(p.lam),cy:Y(p.T),r:3.2,fill:CSS("--teal")}));
  const hm=median(halves);
  F.s.append(el("line",{x1:X(hm),x2:X(hm),y1:m.t,y2:h-m.b,stroke:CSS("--gold"),"stroke-width":1.5,"stroke-dasharray":"5 3"}));
  F.s.append(el("text",{x:X(hm)+6,y:m.t+14,"font-size":11,fill:CSS("--gold")},`half transmission \u2248 ${hm.toFixed(1)} Å`));
  F.s.append(el("line",{x1:X(2.94),x2:X(2.94),y1:m.t,y2:h-m.b,stroke:CSS("--red"),"stroke-width":1.2,"stroke-dasharray":"2 3"}));
  F.s.append(el("text",{x:X(2.94)+5,y:h-m.b-9,"font-size":10,fill:CSS("--red")},"Au–Au 2.94 Å"));
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-5,"text-anchor":"middle","font-size":11,fill:CSS("--ink-2")},
    "in-plane wavelength \u03bb  (Å, log axis)"));
  F.s.append(el("text",{x:14,y:(m.t+h-m.b)/2,"text-anchor":"middle","font-size":11,fill:CSS("--ink-2"),
    transform:`rotate(-90 14 ${(m.t+h-m.b)/2})`},"normalised transmission (log)"));
  document.getElementById("c_tf").append(F.s);
})();

/* ---- dOmega ---- */
(function(){
  const seen=new Set(), rows=[];
  for(const p of D.pairs){
    const a=p.a.split("__"), b=p.b.split("__");
    if(!["ideal","relaxed"].includes(a[1])||!["ideal","relaxed"].includes(b[1])) continue;
    const key=[a[0],b[0]].sort().join("|"); if(seen.has(key)) continue; seen.add(key);
    rows.push({l:`${a[0]} ${a[1]}  vs  ${b[0]} ${b[1]}`,v:p.d_relative_Omega_eV,z:p.dU_pzc_mV});}
  rows.sort((x,y)=>Math.abs(y.v)-Math.abs(x.v));
  const R=rows.slice(0,14); if(!R.length) return;
  const H=31,w=880,m={l:310,r:74,t:32,b:38},h=m.t+m.b+R.length*H;
  const lim=Math.max(0.035,...R.map(r=>Math.abs(r.v)))*1.12;
  const X=x=>m.l+(x+lim)/(2*lim)*(w-m.l-m.r);
  const F=frame(w,h,m), sc=D.same_structure_pair_scale;
  F.s.append(el("rect",{x:X(-sc.max),y:m.t-8,width:X(sc.max)-X(-sc.max),height:h-m.b-m.t+8,
    fill:CSS("--ink-2"),"fill-opacity":.09}));
  F.s.append(el("text",{x:X(0),y:m.t-14,"text-anchor":"middle","font-size":10.5,fill:CSS("--ink-2")},
    `within one morphology: |\u0394\u03a9| at most ${(1000*sc.max).toFixed(0)} meV`));
  for(const t of [-0.03,-0.02,-0.01,0,0.01,0.02,0.03]){
    if(Math.abs(t)>lim) continue;
    F.s.append(el("line",{x1:X(t),x2:X(t),y1:m.t-8,y2:h-m.b,stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:X(t),y:h-m.b+15,"text-anchor":"middle","font-size":11,fill:CSS("--muted")},(1000*t).toFixed(0)));}
  F.s.append(el("line",{x1:X(0),x2:X(0),y1:m.t-8,y2:h-m.b,stroke:CSS("--line-2"),"stroke-width":1.4}));
  R.forEach((r,i)=>{const y=m.t+i*H+H/2;
    F.s.append(el("rect",{x:Math.min(X(0),X(r.v)),y:y-8,width:Math.abs(X(r.v)-X(0)),height:16,rx:2,
      fill:r.v<0?CSS("--teal"):CSS("--gold"),"fill-opacity":.85}));
    F.s.append(el("text",{x:m.l-10,y:y+4,"text-anchor":"end","font-size":11,fill:CSS("--ink")},r.l));
    F.s.append(el("text",{x:w-m.r+8,y:y+4,"font-size":10,fill:CSS("--muted")},
      r.z===null?"":`ΔU₀ ${r.z>0?"+":""}${r.z.toFixed(0)} mV`));});
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-3,"text-anchor":"middle","font-size":11,fill:CSS("--ink-2")},
    "\u0394(\u03a9\u2090 \u2212 \u03a9\u1d66) between U = \u00b10.2 V   (meV)"));
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
    ...fams.map(f=>mk(f,`${FZH[f]} (${Object.values(D.structures).filter(s=>s.f===f).length})`)));
  function render(){
    gal.textContent="";
    for(const [k,v] of Object.entries(D.structures).filter(([k,v])=>active==="all"||v.f===active)
        .sort((a,b)=>FAM.indexOf(a[1].f)-FAM.indexOf(b[1].f)||a[0].localeCompare(b[0]))){
      const tags=[el("span",{class:"tag"},`${v.n} Au`),el("span",{class:"tag"},`area ${v.A} Å²`),
        el("span",{class:"tag"},FZH[v.f]||v.f)];
      if(v.cn.kink) tags.push(el("span",{class:"tag k"},`${v.cn.kink} \u00d7 CN\u22646`));
      if(v.cn.edge) tags.push(el("span",{class:"tag e"},`${v.cn.edge} \u00d7 CN 7\u20138`));
      if(v.C!=null) tags.push(el("span",{class:"tag"},`C ${v.C.toFixed(1)} µF/cm²`));
      if(v.z!=null) tags.push(el("span",{class:"tag"},`U₀ ${v.z>0?"+":""}${(1000*v.z).toFixed(0)} mV`));
      if(v.dz!=null) tags.push(el("span",{class:"tag"},`ΔU₀ ${v.dz>0?"+":""}${v.dz.toFixed(0)} mV`));
      gal.append(el("figure",{class:"gcard",style:"margin:0"},
        el("div",{class:"ghead"},el("h3",{},k),el("span",{class:"small"},v.cfg==="ideal"?"ideal":"relaxed")),
        zoomable(`${D.imgbase}gallery/${k}.png`,`${k} \u00b7 outline, top view and side view`,
          `${k}: outline, top view and side view, atoms coloured by coordination number`),
        v.sch?el("div",{class:"hint",style:"padding-top:10px;padding-bottom:0"},"outline: "+v.sch):null,
        el("div",{class:"gmeta"},...tags)));}
  }
  render();
})();

/* ---- maps ---- */
(function(){
  const host=document.getElementById("mapgrid"); if(!D.mapped||!D.mapped.length) return;
  for(const k of D.mapped){
    const v=D.structures[k]||{};
    host.append(el("figure",{class:"gcard",style:"margin:0"},
      el("div",{class:"ghead"},el("h3",{},k),el("span",{class:"small"},FZH[v.f]||"")),
      zoomable(`${D.imgbase}maps/${k}.png`,`${k} \u00b7 anion maps`,
        `${k}: anion excess, its change with potential, the metal's electron-count change, and the `
        +`coordination map`),
      el("div",{class:"hint"},"click to enlarge")));}
})();
</script>
"""


# The Artifact platform wraps the file in its own document skeleton; GitHub Pages does not, so a standalone build
# has to supply it. These are the parts of that skeleton the page actually relies on.
SKELETON_HEAD = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<meta name="color-scheme" content="light dark">
<meta name="description" content="Constant-potential DFT charging response of 33 Au(111) morphologies at three electron chemical potentials, and how much of it reaches the ions.">
<style>
:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}
body{margin:0}
img{max-width:100%}
[hidden]{display:none!important}
</style>
"""
SKELETON_TAIL = "\n</body>\n</html>\n"


def undercoordinated_vs_terrace():
    """Do under-coordinated columns enrich anions less than the terrace of the SAME structure?

    Computed here instead of being written into the prose. The page used to state 37/41 with four exceptions as
    a literal; recomputing it from regions.json with the definition spelled out gives 37 of 43 with six, four of
    which are the vacancies. The definition, which the old literal did not record: the under-coordinated region
    is the area-weighted mean over the CN <= 8 classes (kink/adatom and edge/rim), compared with the terrace
    class, at each geometry's MOST POSITIVE sampled potential, over ideal and relaxed geometries only."""
    import numpy as _np
    S = json.load(open(f"{ROOT}/analysis/spatial/regions.json"))["states"]
    best = {}
    for v in S.values():
        if v["config"] not in ("ideal", "relaxed"): continue
        g = v["geometry_id"]
        if g not in best or v["U"] > best[g]["U"]: best[g] = v
    rows = []
    for v in best.values():
        R = v["regions"]
        sel = [c for c in ("kink/adatom", "edge/rim") if c in R]
        if not sel or "terrace" not in R: continue
        a = sum(R[c]["area_fraction"] for c in sel)
        K = sum(R[c]["K_rel"] * R[c]["area_fraction"] for c in sel) / a
        rows.append((v["structure_id"], v["config"], 100.0 * (K / R["terrace"]["K_rel"] - 1.0)))
    rel = [r[2] for r in rows]
    lo = [r for r in rows if r[2] < 0]
    ex = [r for r in rows if r[2] >= 0]
    imin = int(_np.argmin(rel))
    return dict(n=len(rows), n_lower=len(lo), median_pct=float(_np.median(rel)),
                min_pct=float(min(rel)), min_structure=rows[imin][0],
                exceptions=[{"structure": r[0], "config": r[1], "pct": r[2]} for r in ex])


def build(data, standalone=False, imgbase="", out_path=None):
    data = dict(data, imgbase=imgbase)
    page = HTML.replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/"))
    if standalone:
        k = page.index("</style>") + len("</style>")
        page = SKELETON_HEAD + page[:k] + "\n</head>\n<body>\n" + page[k:] + SKELETON_TAIL
    open(out_path, "w").write(page)
    return os.path.getsize(out_path)


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--standalone", action="store_true", help="emit a complete HTML document (for GitHub Pages)")
    ap.add_argument("--imgbase", default="", help="prefix for the gallery/ and maps/ image paths")
    ap.add_argument("--out", default=None)
    ap.add_argument("--commit", default=None, help="short sha of the commit the page is published from")
    a = ap.parse_args()
    data = json.load(open(f"{WEB}/data.json"))
    if a.commit: data.setdefault("snapshot", {})["commit"] = a.commit
    # the lead figures' numbers come from the same file the figures were drawn from, so the prose beside a
    # figure cannot drift from the figure
    data["undercoord"] = undercoordinated_vs_terrace()
    mapped = sorted(f[:-4] for f in os.listdir(f"{ROOT}/analysis/maps")) if os.path.exists(f"{ROOT}/analysis/maps") else []
    data["mapped"] = [m for m in mapped if m in data["structures"]]
    out = a.out or f"{WEB}/index.html"
    n = build(data, standalone=a.standalone, imgbase=a.imgbase, out_path=out)
    print(f"{n/1024:.0f} kB -> {out}  ({len(data['structures'])} structures, {len(data['mapped'])} maps"
          f"{', standalone, imgbase=' + repr(a.imgbase) if a.standalone else ''})")


if __name__ == "__main__":
    main()
