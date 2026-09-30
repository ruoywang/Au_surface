#!/usr/bin/env python3
"""Build the standalone HTML report (Chinese) from analysis/web/data.json, inlined because the artifact CSP blocks fetch.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_web_report.py
Output: analysis/web/index.html  (publish with the gallery and map PNGs as supporting files)
"""
import json
import os

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
WEB = f"{ROOT}/analysis/web"

HTML = r"""<title>Au 表面形貌与充电响应</title>
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Noto+Serif+SC:wght@400;500;600&family=Noto+Sans+SC:wght@300;400;500;700&family=IBM+Plex+Mono:wght@400;500&display=swap">
<style>
:root{
  --bg:#f6f4ef; --surface:#fffefb; --surface-2:#f0ece4; --ink:#15181c; --ink-2:#565d65; --muted:#8a9199;
  --line:#e0dbd1; --line-2:#cdc6b9;
  --gold:#96650d; --gold-soft:#e8d9b4; --teal:#1d6a85; --red:#a8382a; --amber:#c4841f;
  --terrace:#b99a3e; --sub:#79828c;
  --shadow:0 1px 2px rgba(30,24,10,.05),0 6px 20px -12px rgba(30,24,10,.18);
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;
  --sans:"Noto Sans SC","PingFang SC","Hiragino Sans GB","Microsoft YaHei",system-ui,sans-serif;
  --serif:"Noto Serif SC","Songti SC","SimSun",Georgia,serif;
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

/* The structure lineage is drawn with page text, not rendered to a PNG: it is almost entirely words, and a
   raster of it shown at column width put every label under 9 CSS px. As markup it stays sharp, selectable and
   readable at any width. */
.lin{background:var(--surface);border:1px solid var(--line);border-radius:10px;overflow:hidden;
  box-shadow:var(--shadow)}
.lin-root{display:flex;flex-wrap:wrap;gap:6px 16px;align-items:baseline;padding:16px 22px;
  background:var(--surface-2);border-bottom:1px solid var(--line-2)}
.lin-root b{font-size:19px}
.lin-root span{font-size:14px;color:var(--ink-2)}
.lin-route{border-top:1px solid var(--line-2)}
.lin-route[data-r="a"]{border-left:5px solid var(--gold-soft)}
.lin-route[data-r="b"]{border-left:5px solid #b9d3de}
.lin-rhead{padding:13px 22px;font-size:15.5px;font-weight:600}
.lin-rhead span{font-weight:400;color:var(--ink-2)}
.lin-row{display:grid;grid-template-columns:minmax(290px,1.1fr) 26px minmax(190px,220px) minmax(320px,1.45fr);
  gap:0 18px;align-items:start;padding:17px 22px;border-top:1px dashed var(--line)}
.lin-op{font-size:15px;line-height:1.72}
.lin-arrow{font-size:20px;color:var(--muted);text-align:center;line-height:1.5}
.lin-fam{border:1px solid var(--line-2);border-radius:8px;padding:9px 13px;background:var(--surface-2)}
.lin-fam b{display:block;font-size:15.5px;line-height:1.4}
.lin-fam span{font-size:13px;color:var(--ink-2)}
.lin-mem{font-size:14.5px;line-height:1.75;color:var(--ink)}
.lin-mem em{display:block;margin-top:6px;font-size:13px;color:var(--muted);font-style:normal}
.lin-cfg{display:grid;grid-template-columns:minmax(230px,1fr) minmax(300px,1.6fr);gap:0 20px;
  padding:16px 22px;border-top:1px dashed var(--line);align-items:start}
.lin-cfg h4{margin:0;font-size:15.5px;font-weight:600}
.lin-cfg p{margin:4px 0 0;font-size:14.5px;line-height:1.72;color:var(--ink-2)}
.lin-cfg.sub{padding-left:60px}
.lin-cfg.sub h4{font-weight:500}
@media (max-width:900px){
  .lin-row{grid-template-columns:1fr;gap:10px}
  .lin-arrow{display:none}
  .lin-cfg{grid-template-columns:1fr;gap:6px}
  .lin-cfg.sub{padding-left:34px}
}
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
  .chartbox::after{content:"← 左右滑动查看完整图表";display:block;font-size:11.5px;color:var(--muted);
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
  <div class="eyebrow">恒电势 DFT · VASP + VASPsol++ · 1 M 隐式电解质</div>
  <h1>Au 表面形貌与充电响应</h1>
  <p class="lede" style="margin-top:14px">以 Au(111) 台面及相关邻晶面（Au(211)、(221)、(332)、(554)）为对象：表面形貌如何改变
  金属容纳电荷的能力，以及这种原子尺度的差别有多少真正传递到液相中的离子。两条主要结果：
  <b>金属上最容易积累正电荷的位置，不一定是阴离子富集最强的位置</b>；
  <b>原子数不变、只把一颗 Au 换个位置，就能明显平移充电曲线</b>——所以"缺陷种类"这个标签不足以概括充电响应。</p>
  <div class="kicker" id="kicker"></div>
  <div class="note" style="margin-top:20px;max-width:none"><b>数据范围（所有图表共用同一口径）。</b><span id="scope"></span></div>
</div></header>

<nav class="sticky"><div class="wrap">
  <a href="#origin">结构从哪来</a>
  <a href="#f-local">整体充电 ≠ 局部富集</a>
  <a href="#f-move">一颗原子换位置</a>
  <a href="#pzc">一 · 变的是零电荷点</a>
  <a href="#which">二 · 哪些形貌，往哪个方向</a>
  <a href="#where">三 · 阴离子究竟在哪里增多</a>
  <a href="#filter">四 · 电解质这个低通滤波器</a>
  <a href="#omega">五 · 电势与相对稳定性</a>
  <a href="#gallery">结构图谱</a>
  <a href="#maps">阴离子空间图</a>
  <a href="#method">方法与边界</a>
  <a href="#revlog">修订记录</a>
</div></nav>

<div class="wrap">

<section id="origin"><div class="finding">
  <div><div class="eyebrow">结构来源</div><h2>从 Au 晶体到缺陷构型：结构库的生成与采样</h2></div>
  <p class="lede">下面这张谱系回答的是"这些结构为什么存在、怎么构建、哪些之间才是受控对照"。
  每一行左边写的是<b>几何操作</b>，不是结构名称，操作取自实际写出 POSCAR 的建构脚本。
  单看后面的 33 张结构图，Island-7-compact 与 Island-7-elongated 像是两个无关的岛；
  谱系说明它们是<b>同一个胞里同样七个原子的两种排布</b>。</p>
  <div class="bleed">__LINEAGE__</div>
  <div class="note"><b>两件容易误会的事。</b>
  <br>· <b>四个高指数面不是"在 Au(111) 上加一条带"得到的。</b>Au(221)/(332)/(554) 由立方胞直接按 (hkl) 切割
  （面间距 a/(2√(h²+k²+l²)) = 0.693 / 0.443 / 0.256&nbsp;&Aring;），构成同一条路线上的台面宽度系列；
  Au(211) 由 <code>ase.build.fcc211</code> 生成，不是同一个生成器，也不在那个宽度系列里。
  <br>· <b><span data-n="ngeom2"></span> 个几何不是 <span data-n="ngeom2"></span> 个独立设计的缺陷。</b>
  其中 33 个是人工建构的；其余是它们在公共参考电势下的弛豫几何，以及由弛豫几何生成的随机位移、集体形变与路径像。
  <br><br><b>这套采样是什么，不是什么。</b>结构类别由晶体学与几何操作构建；弛豫提供局部参考态；扰动与路径用于覆盖
  非平衡构型。它们<b>不是</b>分子动力学中自然出现频率的统计，也<b>不是</b> DFT 自动找出的全部稳定形貌。
  同一个几何再对应若干电势点，但并非每个几何都有五个电势点。</div>
</div></section>

<section id="f-local"><div class="finding">
  <div><div class="eyebrow">主要结果 A</div><h2>金属上最容易积累正电荷的位置，不一定是阴离子富集最强的位置</h2></div>
  <p class="claim">A1-hcp（平板加一个吸附原子）在 U&nbsp;=&nbsp;<span data-n="fl_U"></span> 时，
  吸附原子所在低配位区的平均阴离子富集 <span data-n="fl_Kad"></span>，<b>低于</b>它周围台面的
  <span data-n="fl_Kte"></span>。与此同时整个电极更正：<span data-n="fl_sig"></span>。
  在条带台阶上这种错位更直接——<span data-n="fl_step"></span>。</p>
  <div class="card pad bleed" id="localfig"></div>
  <div class="note"><b>这里不能说"吸附原子排斥阴离子"。</b>K&nbsp;&gt;&nbsp;1 说明该区仍然富集，只是不如台面强。
  准确的说法是：<b>在本模型与本积分口径下，凸出的低配位位点可以增强整体正向充电，却不一定增强其正上方区域的平均阴离子富集。</b>
  这把三件事分开了——整个电极带多少电、金属电荷集中在哪里、离子在哪片可达液体中最多。
  <br><br><b>不是孤例。</b>整个数据集里，低配位区的富集低于同一结构台面的有 <span data-n="corr2"></span>。
  台阶上的错位在五个条带胞（三种台面宽度、两种沿边周期）里都出现，位移 4.3&ndash;4.6&nbsp;&Aring;，方向一致朝下台面。
  <br><br><b>机制上还不能定论。</b>自洽电势、介电屏蔽与离子可达空间在连续电解质模型里是耦合的，仅凭上面的排序
  不足以判定哪一项占主导。本数据也<b>不支持</b>"台阶两条边的金属响应显著不同"：
  Step-16x1 两条边的金属峰是 1.67&times; 与 1.69&times;，差约 1%。</div>
</div></section>

<section id="f-move"><div class="finding">
  <div><div class="eyebrow">主要结果 B</div><h2>原子数不变，只把一颗 Au 换个位置，充电曲线就整体平移</h2></div>
  <p class="claim">Step-8x2 与 Step-8x2_edge-vacancy_plus_foot-adatom 在<b>同一个胞</b>里、同样 72 个 Au，
  差别只是一颗台阶边缘原子移到了脚部空位。零电荷电势 <span data-n="mv_pzc"></span>，
  而割线电容 <span data-n="mv_C"></span>。两条 &sigma;(U) 斜率接近，横向位置不同。</p>
  <div class="card pad bleed" id="movefig"></div>
  <div class="note"><b>有意思的不是"多一个缺陷所以充电变了"。</b>原子数量没有改变，只是放置位置变了，
  充电曲线就明显平移。更有用的表述是：<b>在当前电势区间内，某些局部重排首先改变的是界面充电的零点，
  而不是显著改变整条曲线的斜率。</b>这里不是在比较"毫伏"和"百分比"哪个大，而是在描述两条曲线的形状与位置。
  <br><br><b>弛豫本身也能移动几十毫伏。</b>在 <span data-n="rx_n"></span> 个同时有理想与弛豫几何的结构里，
  弛豫使零电荷点移动的中位幅度 <span data-n="rx_med"></span>，最大 <span data-n="rx_max"></span>
  （<span data-n="rx_top"></span>）。相比之下跨形貌的零电荷点总跨度是 <span data-n="pzcspread2"></span>。
  所以"A1""Au211""R2"这些<b>名字并不足以唯一决定它们的电化学响应</b>；只用各类理想缺陷的一张 POSCAR 做排名，
  容易把理想几何的偶然特征当成整个结构族的规律。这正是采集弛豫、位移与路径构型的理由。
  <br><br><b>不能据此说什么。</b>不能说"正电势使这个原子更容易移动"，也不能说"这个终态更稳定"——
  那需要相对巨势基准与合适的路径能量。当前能确认的是：不同的原子排列具有不同的充电响应。</div>
</div></section>

<section id="pzc"><div class="finding">
  <div><div class="eyebrow">结论一</div><h2>同一电势下的电荷差异，更多体现为零电荷位置的变化</h2></div>
  <p class="claim">在采样点能夹住 &sigma;&nbsp;=&nbsp;0 的几何里，零电荷点跨度约 <span data-n="pzcspread"></span>，
  而割线电容跨度约 <span data-n="cspread"></span>。在本样本和本电势区间内，零电荷位置变化所对应的电荷尺度，
  比割线电容变化所对应的尺度大约一个量级。</p>
  <div class="stats" id="s1"></div>
  <div class="card pad bleed"><div class="chartbox"><div id="c_sigma"></div></div>
    <div class="caption">每个几何的表面电荷密度对内部电势 U&nbsp;=&nbsp;&mu;<sub>0</sub>&nbsp;&minus;&nbsp;&mu;<sub>e</sub>
    的曲线，全部取自实际收敛的 &mu;<sub>e</sub> 与 N<sub>e</sub>。实线是三个基准电势（统计口径）；金色虚线之外的浅色点线是
    已完成的 &plusmn;0.5&nbsp;V 扩展，仅供参看，不进入任何数值。曲线彼此接近平行：族与族之间的差别主要是横向平移。</div>
    <div class="legend" id="leg1"></div></div>
  <div class="note"><b>这个比值是什么，不是什么。</b>某个缺陷在给定电势下带更多正电，并不等于它更容易被极化。把
  &sigma;<sub>i</sub>&nbsp;=&nbsp;C<sub>i</sub>&thinsp;(U&nbsp;&minus;&nbsp;Z<sub>i</sub>) 在参考值附近展开，
  &delta;&sigma;<sub>i</sub>&nbsp;=&nbsp;&minus;C<sub>0</sub>&delta;Z<sub>i</sub>&nbsp;+&nbsp;(U&minus;Z<sub>0</sub>)&delta;C<sub>i</sub>&nbsp;&minus;&nbsp;&delta;C<sub>i</sub>&delta;Z<sub>i</sub>，
  除前两项外还有交叉项，而且跨几何时 &delta;C 与 &delta;Z 本身相关。下表给出的只是前两项按<b>极差</b>估计的<b>尺度</b>
  及其比值，<b>不是方差分解，也不能读成"某一项解释了百分之多少"</b>。
  <br><br>同样地，电容并非"几乎不变"：A3（三吸附原子）相对同胞平板就高约 15%。准确的说法是，在这个电势区间里，
  电容变化带来的电荷差异小于零电荷点平移带来的电荷差异。</div>
  <div class="tablewrap"><table id="t_decomp"></table></div>
  <div class="caption">跨胞比较带有数值偏移：审计记录过两个<em>平板</em>胞仅因胞形、k 点网格与 PREC 就相差 72&nbsp;meV 的中性
  &mu;<sub>e</sub>。同一个胞内部这项偏移对每一行是共同的，所以按胞分组更有利于隔离形貌效应，但这不等于所有数值误差都完全抵消。此外，全部
  <span data-n="ngeom"></span> 个几何里包含人为扰动、集体形变和不同胞／采样设置，整体跨度不应全部归因于"缺陷种类本身"。</div>
</div></section>

<section id="which"><div class="finding">
  <div><div class="eyebrow">结论二</div><h2>在本组匹配对照中，凸起把零电荷点推负，空位把它推正</h2></div>
  <p class="claim">以<em>同一个胞</em>里的平整平台为基准：本数据集中的吸附原子与台阶把零电荷点压低 29&ndash;97&nbsp;mV，
  三个空位构型则抬高 8&ndash;13&nbsp;mV，两种堆垛重构几乎不动。符号与经典的功函数论证一致——凸起抹平电子溢出、
  降低功函数，凹陷相反。这是对这批构型的描述，不是对所有凸／凹形貌的普遍规律。</p>
  <div class="card pad"><div class="chartbox"><div id="c_dpzc"></div></div>
    <div class="caption">只取理想构型，各自对照同一个胞里的平板成员。柱长是零电荷点的位移，右侧圆点标的是同一对照下的电容变化。</div></div>
  <div class="note"><b>电容同向变化，但幅度小得多。</b>加上去的原子比挖掉的影响大：三个吸附原子 +15%，单个吸附原子 +8%，
  一条台阶 +8%，而单个空位只有 +1%，两种重构在 1% 以内。凸向电解质的粗糙度既压低零电荷点又抬高电容。
  <br><br><b>注意符号。</b>由 &sigma;&nbsp;=&nbsp;C&thinsp;(U&nbsp;&minus;&nbsp;U<sub>pzc</sub>)，固定 U 时
  &Delta;&sigma;&nbsp;=&nbsp;&minus;C&thinsp;&Delta;U<sub>pzc</sub>：<b>零电荷点负移意味着同一电势下整体带更多正电</b>。
  实际数据正是如此——A1-fcc 的零电荷点比同胞 T-4x4 低 65&nbsp;mV，在 U&nbsp;&asymp;&nbsp;+0.2&nbsp;V 时
  &sigma; 是 +3.18 对 +2.22&nbsp;&micro;C/cm²。</div>
</div></section>

<section id="where"><div class="finding">
  <div><div class="eyebrow">结论三</div><h2>整体充电变强，与某处离子变多，不是同一件事</h2></div>
  <p class="claim">按最近表面原子的配位数给阴离子过量分区后，低配位位点（CN&nbsp;&le;&nbsp;8，按面积加权）的区域平均富集
  <em>低于</em>同一结构的平台区：<span data-n="uc_main"></span>。取每个几何最正的采样电势，只用理想与弛豫几何。例外见
  下面的例外表；两种积分口径给出相同的排序。</p>
  <div class="note" style="margin-top:0"><b>例外。</b><span data-n="uc_exc"></span>
  这个计数以前是写死在正文里的 37/41、四个例外；现在由 <code>regions.json</code> 在建站时重算并写入，
  定义（CN&nbsp;&le;&nbsp;8 按面积加权、取最正采样电势、只用理想与弛豫几何）也一并写出。</div>
  <p class="claim" style="border-left-color:var(--teal)">逐柱的空间相关进一步支持这一点：把金属的正电荷增量
  (&minus;&Delta;n<sub>e</sub>) 与阴离子增量 (&Delta;&Gamma;<sub>&minus;</sub>) 按柱子求相关，
  <span data-n="corr"></span>。这描述的是整张图上的<b>空间共变关系</b>：两者的起伏总体反号。它不直接给出峰位，也不保证每个结构的极大值都错开；具体位置请看下方各结构的空间图。</p>
  <div class="card pad"><div class="chartbox"><div id="c_regions"></div></div>
    <div class="caption">边界相对窗口下的富集比
    K<sub>&Omega;</sub>&nbsp;=&nbsp;&int;n<sub>&minus;</sub>&thinsp;/&thinsp;(n<sub>b</sub>&int;S<sub>ion</sub>)，
    取各结构最正的那个采样电势。柱按配位类分组，虚线是该结构的全胞平均。此图每种结构只画一个代表几何（理想或弛豫），
    共 33 条；上文 41 个的计数则覆盖同一结构的理想与弛豫两种几何。</div>
    <div class="legend" id="leg3"></div></div>
  <div class="note"><b>整体充电与局部离子分布并非简单对应。</b>即使整体零电荷点负移、同一电势下金属总正电荷更多，
  低配位原子对应的液相区域也可能没有更高的区域平均阴离子浓度。局部分布取决于自洽电势、离子可达性与几何分配三者共同作用。
  这比原来的错误解释更值得注意。
  <br><br>K 与 &Gamma; 回答的是不同问题，两者都给出：小区域可以浓度很高却只容纳很少额外离子，大区域可以只略微富集却贡献
  大部分总过量。图中同时给出各区面积占比，避免混淆。区域平均本身只能说明"该区平均浓度较低"；"热点在哪里"由上面的
  逐柱相关和下方的空间图支持。</div>
</div></section>

<section id="filter"><div class="finding">
  <div><div class="eyebrow">结论四</div><h2>离子响应里的短波空间起伏明显减弱</h2></div>
  <p class="claim">把金属电荷响应的每一个傅里叶模式与阴离子响应的同一模式相比，相对谱幅随波长<b>缩短</b>而单调降低：原子尺度
  （&lambda;&nbsp;&asymp;&nbsp;1&nbsp;&Aring;）只剩参考值的百分之零点几，缺陷尺度（&lambda;&nbsp;&asymp;&nbsp;9&nbsp;&Aring;）
  还有三成以上。这与空间屏蔽／平滑的图像相容。</p>
  <div class="stats" id="s4"></div>
  <div class="card pad"><div class="chartbox"><div id="c_tf"></div></div>
    <div class="caption">相对谱幅 |F<sub>ion</sub>(k)|&thinsp;/&thinsp;|F<sub>metal</sub>(k)| 对横向波长，汇总了所有
    两端电势齐全的几何。电子数面密度与离子数面密度量纲相同，但它们是不同的物理量，其比值没有绝对标度，所以每条曲线按<b>自身最长波长那一档</b>归一——这是分析选择，不是量纲所迫；因此长波端的
    1.0 是归一化的定义，<b>不代表物理上 100% 传递</b>，而且不同晶胞可用的最小波数不同，各条曲线的参考档并不完全一致。
    只有曲线形状有意义。阴影带是四分位距。</div></div>
  <div class="note"><b>两个数字是指标，不是材料常数。</b>下面两项都依赖上面说的归一化和分箱方式，应按定义读：
  <br>· <b><span data-n="half"></span></b>：相对谱幅降到各自参考值一半处的特征波长。这是"在当前傅里叶分箱、筛选与
  逐几何归一化定义下"的经验指标，<b>不是电解质的普适分辨率</b>。
  <br>· <b><span data-n="ratio"></span></b>：两个归一化对比度
  (max&minus;min)/&lang;|f|&rang; 的比值。它<b>不是</b>百分之多少的电荷、多少个傅里叶模式或多少信息被保留下来。
  <br><br><b>不把短波衰减唯一归因于某一项机制。</b>模型同时含非局域空腔、非线性介电与离子响应，三者耦合。观察到的谱幅随
  波长下降与空间屏蔽相容，但仅凭一套参数下的比值曲线，不能断言短波段由 4&nbsp;&Aring; 空腔而非德拜屏蔽决定。
  线性化 Poisson&ndash;Boltzmann 平板只作趋势参照，不是拟合。</div>
</div></section>

<section id="omega"><div class="finding">
  <div><div class="eyebrow">结论五</div><h2>电势对同组成构型相对巨势的贡献</h2></div>
  <p class="claim">对电子数曲线积分，可以得到电势<em>诱导</em>的相对巨势变化，全程不需要把两个总能相减。在
  &plusmn;0.2&nbsp;V 窗口内，跨形貌的最大值是 <span data-n="dom"></span>，与同一形貌自身各采样构型之间的散布
  （最大 <span data-n="domsame"></span>）同量级。<b>这不足以判断电势是否改变了稳定性排序</b>：那还需要参考电势下的
  相对巨势基准，本数据集尚未确定。</p>
  <div class="card pad"><div class="chartbox"><div id="c_omega"></div></div>
    <div class="caption">同成分、同胞、不同结构的配对，
    D&nbsp;=&nbsp;+&int;<sub>&mu;₀&minus;w</sub><sup>&mu;₀+w</sup>&thinsp;[N<sub>A</sub>&minus;N<sub>B</sub>]&thinsp;d&mu;，
    w&nbsp;=&nbsp;0.2&nbsp;V。负值表示电势往正方向移动时 A 相对被稳定。灰带是同一形貌自身采样构型之间配对的散布范围。
    符号：因 &part;&Omega;/&part;&mu;<sub>e</sub>&nbsp;=&nbsp;&minus;N<sub>e</sub> 且
    U&nbsp;=&nbsp;&mu;₀&minus;&mu;<sub>e</sub>，U&nbsp;=&nbsp;+w 对应<em>较低</em>的 &mu;。</div></div>
  <div class="note"><b>这是什么，不是什么。</b>它只是电势<em>诱导</em>的那一部分。它不回答在参考电势下谁更稳定，
  那需要数据集尚未确定的统一能量基准；没有这个基准，就无法判断排序是否发生了变化——若两个构型原本只差几 meV，
  几十 meV 足以改变排序；若原本差 1&nbsp;eV，则未必。它也从不跨越不同的 Au 原子数，那需要引入储库项。
  <br><br>&plusmn;0.5&nbsp;V 扩展完成后，应当用扩展后的实际曲线重新积分，而不是把现在的数值按窗口比例外推。</div>
</div></section>

<section id="gallery"><div class="finding">
  <div><div class="eyebrow">结构图谱</div><h2>全部三十三种结构，按配位数着色</h2></div>
  <p class="lede">每张图分三块：最左是<b>不画原子的简笔轮廓</b>，一眼看出缺陷长什么样；中间是按配位数着色的俯视图；
  右边是侧视图。三块全部由实际参与计算的几何生成，不是重新生成的理想结构，也不是手绘，所以简笔图不会和真实结构脱节。
  <b>点击任意图片可放大。</b></p>
  <div class="note" style="max-width:none"><b>图上各个标记的含义。</b>
  <br>· <b>颜色</b>＝配位数（3.4&nbsp;&Aring; 内的 Au 近邻数），括号里只是该配位数<em>常见</em>的环境，配位数不等于形貌身份。
  <br>· <b>亮度与描边</b>＝是否属于上表面，判据与区域分析一致（上方 2.35&nbsp;&Aring; 内有三个以上更高近邻即算被埋住）。
  <br>· <b>绿色圆环</b>＝加上去的 Au；<b>紫色虚圆</b>＝移走的位点。只在结构族本身由"加原子"或"去原子"定义时才标，
  避免把台阶上台面或邻晶面的台阶边误标成吸附原子。
  <br>· <b>蓝色 ×</b>＝类 hcp 配准（正上方两层处有原子）；<b>棕色空心方块</b>＝过渡带／畴界（配准介于类 fcc 与
  类 hcp 之间）。阈值按面内偏移量取 0.35 与 0.75 倍 fcc 偏移。这是一套<b>阈值相关的显示分类</b>，不是畴区的唯一定义，
  与建构时另一套阈值给出的计数不必相同，也不应为了复现旧计数去调阈值。参照的 fcc 偏移取自同层最近邻间距
  （R2 的压缩层量得 2.753&nbsp;&Aring;，而不是跨层投影给出的 4.40&nbsp;&Aring;）。
  这些标记<b>只画在 (111) 法向的结构上</b>：判据用的是全局 z 与 (111) 层间距，在 Au(211)/(221)/(332)/(554) 上没有意义，
  曾经在那里点出过几个无法解释的标记，现已撤下。
  <br>· <b>A–A′</b>（必要时加 <b>B–B′</b>）＝侧面剖线的实际位置。标题注明是真实剖线还是沿视线的投影包络；
  剖面以<b>平台高度</b>为基准，凸起在基准线上、凹陷在基准线下。
  <br>· <b>棕色粗虚线</b>＝母结构的原边缘。只画在由某个已存在结构改出来的结构上，取自建构时的原子映射。
  <br>· 有限特征（岛、坑、点缺陷、复合）的原子侧视图只画 A–A′ 附近的<b>一条窄带</b>。整胞投影会把坑前后两侧的顶层原子
  叠在坑上面，看起来像被填平了。窄带宽度是<b>真实垂直距离</b>（先前按"分数坐标差×另一根基矢长度"算，
  在 60° 胞里比标注值窄 1/sin60，标 ±3.23&nbsp;&Aring; 实为 ±2.80&nbsp;&Aring;）。
  <br><br><b>轮廓来自原子邻接，不是高度分区，尺度取同层间距。</b>高度分区不是结构的外轮廓：较低的原子会在两个上层
  原子之间抢到区域，把连通的团簇切成两半。现在轮廓是特征原子的圆盘并集，半径 0.62&nbsp;×&nbsp;同层最近邻间距，
  <b>同层</b>三个字是关键：先前那一步在投影里取最近邻，而吸附原子与其下方台面原子的投影距离只有
  a₀/√3&nbsp;=&nbsp;1.70&nbsp;&Aring;，在 A3、C1 和小坑上过半原子的"最近邻"都是这种跨层对，中位数因此给出
  1.70 而不是 2.94&nbsp;&Aring;，圆盘直径 2.10&nbsp;&Aring; 跨不过 2.94&nbsp;&Aring; 的键——A3 的三原子团被画成三个圆，
  Pit-7-compact 的坑底被画成十二个。逐结构复核过：圆盘连通分量数与原子成键连通分量数<b>逐一相等</b>。
  <br><br><b>周期副本是平移复制的。</b>距离场只在基本胞上算一次，再原样平移；先前给每个显示副本重新调用一次搜索范围
  只有 ±1 个像的距离函数，离中心两个胞远的副本就丢掉了特征——Step-8x2 三个副本的高区像素比例是
  53.18%／51.11%／0.40%，同一个周期结构的副本不应有这种差别。
  <br><br><b>剖面与 A–A′ 逐点对应。</b>先前只对剖面做循环移位而不动图上的 A、A′，两者因此脱节：Pit-19-8x8 沿所标
  A–A′ 的采样顺序是平台 16 → 坑底 75 → 平台 8，移位后成了坑底 38 → 平台 24 → 坑底 37，这正是"坑看起来像中央凸起"
  的原因。现在要把特征居中时，平移的是整套显示坐标（起点 A、剖线、剖面、原子窄带同时移动），不是单独一条曲线。
  原子侧视图也改用沿剖线／垂直剖线的同一坐标系，不再一律投影到全局 x 轴。
  <br><br><b>一个原子只取一个周期镜像。</b>沿线坐标与垂直距离曾各自独立折叠：垂直方向折回最近的周期带，沿线方向
  却仍用原始位置。斜胞里跨线的晶格矢量本身带有沿线分量，这样做不成立。取
  <b>a</b>&nbsp;=&nbsp;(10,&nbsp;0)、<b>b</b>&nbsp;=&nbsp;(5,&nbsp;8.6603)，同一个原子写成 <b>r</b> 与
  <b>r</b>+<b>b</b> 时，沿线坐标分别是 9.5 与 4.5&nbsp;&Aring;——仅仅换了写法就差半个周期。现在先按垂直距离选定镜像
  <b>r</b>−m<b>b</b>，两个坐标都从这一个镜像读出。
  <br><br><b>绘图有回归测试，并且测试本身被变异测试验证过。</b><code>scripts/check_gallery_drawing.py</code> 对 33 个结构
  断言五件事：圆盘轮廓的连通分量数等于特征原子成键图的连通分量数（键长阈值取胞内最短 Au–Au 距离的 1.15 倍，
  <b>不</b>取自被测的间距函数）；每个周期副本都是基本胞数组的平移且水平集覆盖率完全相同；剖面等于在
  A&nbsp;+&nbsp;s·t̂ 独立重采样的场；侧视窄带的宽度与坐标与一次独立的镜像搜索一致；整体平移一个面内晶格矢量后显示不变。
  失败时<b>非零退出</b>。<code>scripts/check_drawing_mutations.py</code> 把上述四个 bug 逐个放回去，要求回归全部报错——
  四个都被抓到，其中"圆盘半径取投影间距"这一项报出 8 条连通性不符，包括 A3（1 个成键团簇、3 个圆）与
  Pit-7-compact（1 个、12 个）。
  <br><br><b>图中没有画"离子可达边界"。</b>可达区域由 SION 决定，会随台面、岛、坑变化，用一条固定高度的水平线代表它
  与本研究的结论冲突，因此只标"电解液侧"。
  <br><br><b>几类结构的专门标注。</b>
  <br>· <b>四个邻晶面</b>（Au(211)/(221)/(332)/(554)）：俯视的渐变图只是高度分布，承担不了"台阶示意"，
  所以剖面改成沿落差最大的晶格方向的真实剖线，并<b>定性</b>标出局部 (111) 台面与周期边界处的台阶。
  唯一给出的数字是<b>宏观 (hkl) 与 (111) 的名义夹角</b>，由 arccos[(h+k+l)/(√3·√(h²+k²+l²))] 从晶面法向算出，
  <b>不是本图的测量值</b>。台面宽度与台阶高度<b>不标</b>：那需要先指定原子行并定义测量方向，本图没有做这件事。
  <br>　　（先前这里标过一个"倾角"，取 ptp(剖面)/(L−网格步长)，即高度分区场沿单一晶格方向、跨近一个周期的表观斜率。
  当剖线方向恰是最陡下降且台面很宽时它接近真值——Au554 给 5.6° 对名义 5.77°——但这不是同一个量，
  也不普遍接近：Au211 给 25.8° 对名义 19.47°。已撤下。）
  <br>· <b>拐角与边缘脱离终态</b>：叠加母结构的原直边（棕色虚线），并圈出改变的原子。Kink-edge1/2 的母结构是同一条
  无拐角的条带，按建构规则复原（条带是四整行各三个原子，外加独立一行里的一个原子，那个就是拐角原子）；
  Step-8x2_edge-vacancy_plus_foot-adatom 直接与 Step-8x2 逐原子比对，得到 1 个新增、1 个移走。
  <br>· <b>C1</b>：七个岛原子与台阶上台面<b>同高</b>，任何"高于多数原子"的几何判据都看不见它们，因此改用建构时的
  原子映射：与母结构 Step-8x4 逐原子比对，得 7 个新增、0 个移走。剖线也被拉到穿过这七个原子的那一行。
  <br>· <b>C2</b>：岛与坑不在同一条晶格线上，沿晶格方向的单条直线只能擦过其中一个的边缘（会把一个坑画成三个浅坑），
  所以画 A–A′ 过岛、B–B′ 过坑<b>两条</b>剖面，各配一张原子窄带侧视图。
  <br><br><b>简笔轮廓用于辅助辨认形貌，不要用于定量读取边界位置、宽度或峰位。</b>形貌类别由冻结计划的结构族加周期
  连通性判定，不由高低区面积比判定。</div>
  <div class="filters" id="filters"></div>
  <div class="gal bleed" id="gal"></div>
</div></section>

<section id="maps"><div class="finding">
  <div><div class="eyebrow">阴离子空间图</div><h2>逐个结构的阴离子分布</h2></div>
  <p class="lede">四张一组：最正采样电势下的逐柱阴离子过量、它在采样窗口两端之间的变化、同一电势步长下金属自身的正电荷变化
  (&minus;&Delta;n<sub>e</sub>)，以及切分区域所依据的配位数图。按周期平铺以便看清重复。<b>点击可放大。</b></p>
  <div class="note" style="max-width:none"><b>怎么读色标。</b>每幅图有<b>自己的</b>色标范围，<b>不同图之间不能直接比较幅度</b>，
  只能比较图内的空间分布。全为同号的量用顺序色（深=大），跨正负的量用以 0 为中心的发散色。单位：
  &Gamma;<sub>&minus;</sub> 与 &Delta;&Gamma;<sub>&minus;</sub> 是每投影面积的离子数（&Aring;<sup>&minus;2</sup>）；
  金属一侧画的是<b>正电荷</b>变化 &minus;&Delta;n<sub>e</sub>（e/&Aring;<sup>2</sup>，正值=失去电子）；配位数是离散标度。</div>
  <div class="gal bleed" id="mapgrid"></div>
</div></section>

<section id="method"><div class="finding">
  <div><div class="eyebrow">方法</div><h2>约定，以及这些数字不覆盖的范围</h2></div>
  <div class="grid2">
    <div class="card pad"><h3>电势与电荷</h3><p class="small" style="margin-top:8px">
      U&nbsp;=&nbsp;&mu;<sub>0</sub>&nbsp;&minus;&nbsp;&mu;<sub>e</sub>，&mu;<sub>0</sub>&nbsp;=&nbsp;&minus;4.9071&nbsp;eV
      是项目内部参考，既不是相对 RHE 的电势，也不是各结构自己的零电荷电势。
      &sigma;&nbsp;=&nbsp;&minus;e&thinsp;(N<sub>e</sub>&nbsp;&minus;&nbsp;N<sub>e</sub><sup>0</sup>)&thinsp;/&thinsp;A<sub>proj</sub>。
      所有数值都取自运行中实际收敛的 &mu;<sub>e</sub> 与 N<sub>e</sub>，从不用目标值代替，因为收敛判据允许运行在离目标
      10&nbsp;meV 处停下。</p></div>
    <div class="card pad"><h3>离子密度</h3><p class="small" style="margin-top:8px">
      由收敛后的电势和离子可达性掩码，通过模型自身的本构关系重建，含有限尺寸饱和项。没有加入任何显式离子，体系里也没有氯：
      这是连续介质电解质中非特异性的阴离子响应，不等同于氯的化学吸附偏好。</p></div>
    <div class="card pad"><h3>积分窗口</h3><p class="small" style="margin-top:8px">
      绝对口径取 z&nbsp;&lt;&nbsp;31&nbsp;&Aring; 以下的全部空间，各柱之和精确等于全胞；边界相对口径从每根柱自己的可达性
      边界往上取 10.4&nbsp;&Aring;，使抬高的岛与旁边的平台保留同样比例的衰减尾部。两种口径都给出，上述结论在两者下都成立。</p></div>
    <div class="card pad"><h3>区域怎么切（一次已修正的定义）</h3><p class="small" style="margin-top:8px">
      每根柱子指派给最近的<b>未被埋住</b>的表面原子，按该原子的配位数归类。最初的版本把"顶端 3&nbsp;&Aring; 内的所有原子"
      都当作候选，但 (111) 层间距只有 2.4&nbsp;&Aring;，第二层因此也成了候选，而它恰好位于空位点的正下方——在完全平整的
      T-4x4 上有 48% 的柱子被指派给 CN&nbsp;12 的第二层原子，把一个平整平台劈成了"平台"和"次表面"两类。改用未被埋住的
      原子后，平整面回到 100% 平台；上面的结论在修正前后都成立（修正前 36/41、中位 &minus;1.7%）。</p></div>
    <div class="card pad"><h3>数据范围的处理</h3><p class="small" style="margin-top:8px">
      每一个跨结构数值都只用三个基准电势（&minus;5.1071 / &minus;4.9071 / &minus;4.7071&nbsp;eV）。注意这三点并非每个几何都齐全：按计划，部分参考结构本来就只算了一或两个电势，需要三点的量（零电荷点、两段割线、巨势积分）只在点数够的几何上计算，点数不足的不进入同一区间的比较，各图下方给出实际参与的个数。
      &plusmn;0.5&nbsp;V 扩展仍在计算、覆盖不全，只在 &sigma;(U) 曲线上单独叠加显示，<b>不进入任何统计量</b>——
      否则会拿采样到 &plusmn;0.5&nbsp;V 的几何去和只采样到 &plusmn;0.2&nbsp;V 的几何比较，而且每完成一个任务
      页面上的数字就会变一次。页首的快照给出本页实际使用的态数、几何数与构建时间。</p></div>
    <div class="card pad"><h3>修订记录</h3><p class="small" style="margin-top:8px">
      本页经过两轮外部审查并据此修改。正文只写当前正确的解释；改动了什么、为什么改，收在页尾可展开的
      <a href="#revlog">修订记录</a>里。原始数据与脚本都在仓库中可追溯。</p></div>
    <div class="card pad"><h3>已知限制</h3><p class="small" style="margin-top:8px">
      生产参数下的力带有约 0.02&nbsp;eV/&Aring; 每原子的 egg-box 误差，因此这里没有任何结论建立在细小的力差异上。
      扰动构型与集体形变是人为设计的采样，不是热力学系综，只用于给出敏感性范围而不做平均。路径像很稀疏，不等于最小能量路径。</p></div>
  </div>
</div></section>

<section id="revlog"><div class="finding">
  <div><div class="eyebrow">修订记录</div><h2>改过什么，为什么改</h2></div>
  <p class="lede">正文只保留当前正确的解释。这里记录历次外部审查之后的实质改动，便于追溯；每一项在改之前都用数据或
  单元检验复核过。</p>
  <details open><summary style="cursor:pointer;font-weight:500;padding:10px 0">最近一轮（结构来源与主要结果）</summary>
  <ul class="small" style="line-height:1.9;max-width:62em">
    <li><b>新增结构生成谱系。</b>此前 33 张结构图说明了每个结构长什么样，却没有说明哪些之间是受控对照。
      谱系把每一步几何操作写出来，并把四个高指数面从"Au(111) 加条带"那一支里分出来——它们由不同的生成器
      直接按晶面切割。它用<b>页面文字</b>排版，不是一张图片：内容几乎全是字，栅格化之后被页面缩到栏宽就看不清了。</li>
    <li><b>所有图改为整行显示，并放宽到栏宽之外。</b>结构图与空间图原本两张一行、各约 560&nbsp;px，
      12 英寸画布上的 9&nbsp;pt 标签在屏幕上只有 7&nbsp;px，等于看不清。现在每行一张、最宽 1780&nbsp;px，
      同样的标签约 18&nbsp;px。空间图由 1×4 改为 2×2（画布 17.2 → 10 英寸），字号同时调大；
      SVG 图表也改为随容器放大。</li>
    <li><b>低配位区与平台的比较由写死改为重算。</b>正文原来写"41 个几何中 37 个、四个例外全是空位边缘，
      超出量只有 0.1%–0.3%"，而这个数字没有记录它的定义，也无法从 <code>regions.json</code> 复现。
      现在建站时按明确定义重算：<span data-n="uc_main"></span>；例外是
      <span data-n="uc_exc"></span>旧说法漏掉了其中两个，并低估了最大超出量。</li>
    <li><b>邻晶面剖面上的"台面倾角"撤下。</b>原先标的是高度分区剖面的 ptp/(L−网格步长)，即沿单一晶格方向、
      跨近一个周期的表观斜率。Au554 给 5.6° 对名义 5.77° 只是巧合；Au211 给 25.8° 对名义 19.47°。
      现在只标由晶面法向算出的名义夹角，并注明不是本图的测量值。</li>
  </ul></details>
  <details><summary style="cursor:pointer;font-weight:500;padding:10px 0">第二轮（图示表达）</summary>
  <ul class="small" style="line-height:1.9;max-width:62em">
    <li><b>删除图中人为画的"离子可达边界"。</b>原来在原子侧视图上于最高原子上方固定 4.2&nbsp;&Aring;、在简笔剪影上固定
      3.2&nbsp;&Aring; 各画一条水平虚线，两个偏移不同，都不是从 SION 读出来的。可达区域随台面、岛、坑变化，正是本研究的
      结论之一，用一条固定水平线代表它自相矛盾。现在只标"电解液侧"。</li>
    <li><b>形貌类别不再由面积比判定。</b>原规则把"高区占 30–70%"一律叫台阶，于是 Pit-19-8x8 这个有限浅坑被标成
      "两层平台（台阶）"。现在由冻结计划中的结构族加周期连通性判定：贯穿整胞的才是条带台阶，闭合的是岛或坑。</li>
    <li><b>剖线不再冒充。</b>过特征中心取不到时会退回投影包络，标题现在如实区分"过 A–A′ 的真实剖线"与"沿视线方向的
      投影包络"，并在俯视图上画出 A–A′ 的实际位置；剖面下方的实体注明只是基底示意。</li>
    <li><b>补上高度看不见的信息。</b>新增 hcp 配准标记（蓝色 ×），R1-hcp 终止面与平整 fcc 面、A1-hcp 与 A1-fcc
      此前在图上完全无法区分；新增"加上去的 Au"与"移走的位点"标注，按结构族门控以免把台阶上台面误标成吸附原子。</li>
    <li><b>配位数统计改为只针对上表面。</b>整块 slab 的统计把背面那一层也算进去，四层平整 slab 会写成"平台 32"，
      而面向电解液的只有 16 个。</li>
    <li>傅里叶段落的方向写反（应为"随波长缩短而降低"）；"不同量纲"改为"不同物理量"；逐柱相关不再被说成峰位定位；
      "同胞内数值偏移完全相同"改为"更有利于隔离形貌效应，但不保证完全抵消"；"三个基准电势每个几何都齐全"更正为
      部分参考结构本来就只有一两个电势。</li>
  </ul></details>
  <details><summary style="cursor:pointer;font-weight:500;padding:10px 0">第一轮（科学与代码）</summary>
  <ul class="small" style="line-height:1.9;max-width:62em">
    <li><b>巨势积分符号反了。</b>因 &part;&Omega;/&part;&mu;&nbsp;=&nbsp;&minus;N 且 U&nbsp;=&nbsp;&mu;₀&minus;&mu;，
      U&nbsp;=&nbsp;+w 对应更低的 &mu;，积分前不应再加负号。单元检验：令 N<sub>A</sub>&minus;N<sub>B</sub>&nbsp;=&nbsp;1、
      窗口 &plusmn;0.2&nbsp;V，应得 +0.4&nbsp;eV，原代码给 &minus;0.4&nbsp;eV，即把"谁被正电势稳定"说反了。</li>
    <li><b>结论三的机制解释反了，且混淆了全局与局部。</b>由 &sigma;&nbsp;=&nbsp;C(U&minus;U<sub>pzc</sub>)，零电荷点
      负移对应同一电势下<em>更正</em>，不是更不正；而且这里的零电荷点是整个周期体系的标量，不能赋给某个分区。已撤回，
      改为数据真正支持的说法。</li>
    <li><b>金属—离子相位比较的符号。</b>正向充电时电子减少、阴离子增加，直接比较会给空间上同位的图案凭空加上
      &pi; 相位。改为比较 &minus;&Delta;n<sub>e</sub> 与 &Delta;&Gamma;<sub>&minus;</sub>。</li>
    <li><b>统计口径被未完成的扩展污染。</b>充电分析在陆续读入 &plusmn;0.5&nbsp;V 的态，107 个几何里已有 68 个是五点，
      而页面仍写"三个电势"。现在跨结构数值只用三个基准电势，扩展单独叠加显示。</li>
    <li><b>"八成来自零电荷点"不成立。</b>极差之比是尺度比较，不是方差分解（展开式有交叉项，&delta;C 与 &delta;Z 相关）。
      同时删去"几乎不改变电容"（A3 相对同胞平板高 15%）。</li>
    <li>撤回"电势没有改变稳定性排序"（缺参考电势基准）与"&plusmn;0.5 将同比放大 2.5 倍"；半衰波长与对比度比降级为
      定义依赖的经验指标；短波衰减不再唯一归因于 4&nbsp;&Aring; 空腔。</li>
    <li><b>区域划分把平整平台劈成两半。</b>柱体原按"顶端 3&nbsp;&Aring; 内最近原子"指派，而层间距只有 2.4&nbsp;&Aring;，
      第二层恰在空位点正下方：完全平整的 T-4x4 有 48% 的柱子被标成 CN&nbsp;12。改为指派给未被埋住的原子后，
      平整面回到 100% 平台，结论由 36/41 变为 37/41。</li>
    <li>&sigma; 图纵轴写死 &plusmn;7&nbsp;&micro;C/cm² 而数据到 7.32，曲线被裁；改为按数据推导。</li>
  </ul></details>
</div></section>

</div>
<footer><div class="wrap">
  由 <span class="mono">dataset_v1/states.json</span> 与逐柱场归约生成。<span id="stamp"></span>
</div></footer>

<dialog class="lb" id="lb">
  <div class="bar"><span id="lbcap"></span><button type="button" id="lbclose">关闭 (Esc)</button></div>
  <img id="lbimg" alt="">
</dialog>

<script type="application/json" id="D">__DATA__</script>
<script>
const D = JSON.parse(document.getElementById("D").textContent);
const FAM = ["flat Au(111)","point defect","reconstruction-related","strip step","vicinal step face",
             "kink / edge rearrangement","single-layer island","single-layer pit","composite"];
const FZH = {"flat Au(111)":"平整 Au(111)","point defect":"点缺陷","reconstruction-related":"重构相关",
  "strip step":"条带台阶","vicinal step face":"邻晶面台阶","kink / edge rearrangement":"拐角 / 边缘重排",
  "single-layer island":"单层岛","single-layer pit":"单层坑","composite":"复合形貌"};
const CZH = {"kink/adatom":"拐角 / 吸附原子","edge/rim":"边缘","terrace":"平整平台","sub-surface/foot":"台阶脚 / 次表面"};
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
  const b=el("button",{type:"button",class:"zoom","aria-label":`放大 ${cap}`},
    el("img",{src:src,alt:alt,loading:"lazy"}));
  b.addEventListener("click",()=>{LBI.src=src;LBI.alt=alt;LBC.textContent=cap;
    if(LB.showModal) LB.showModal(); else LB.setAttribute("open","");});
  return b;}

/* ---- kicker + data snapshot, all generated from the same object the charts read ---- */
const SN=D.snapshot||{};
const NG=Object.keys(D.geometries).length, NS=Object.keys(D.structures).length;
const NPT=Object.values(D.geometries).reduce((a,g)=>a+g.pts.length,0);
const NEXT=Object.values(D.geometries).reduce((a,g)=>a+(g.ext?g.ext.length:0),0);
document.getElementById("kicker").append(
  ...[[NS,"种结构"],[NG,"个独立几何"],[NPT,"个基准态"],["3","个基准电势"],[NEXT,"个 ±0.5 V 扩展态（仅叠加显示）"]]
     .map(([v,l])=>el("span",{class:"chip"},`${v} ${l}`)));
document.getElementById("scope").append(
  `统计口径固定为三个基准电势 −5.1071 / −4.9071 / −4.7071 eV，共 ${NPT} 个态、`
  +`${NG} 个几何、${NS} 种结构。±0.5 V 扩展目前有 ${NEXT} 个态，覆盖不全，只在 σ(U) 图上以浅色叠加，不进入任何数值。`
  +`区域与空间分析覆盖 ${SN.n_region_states||"?"} 个态（${SN.n_region_structures||"?"} 种结构），`
  +`两端电势齐全、可做金属—离子对照的几何 ${Object.keys(D.transmission).length} 个。`,
  el("br"),
  el("span",{class:"mono",style:"font-size:12px"},
    `构建于 ${SN.built||"?"}${SN.commit?"　·　commit "+SN.commit:""}`));
document.getElementById("stamp").textContent = ` · 快照 ${SN.built||""}${SN.commit?" / "+SN.commit:""}`;

/* ---- Finding 1 stats ---- */
const zs=Object.values(D.geometries).filter(g=>g.z!==null).map(g=>g.z);
const cs=Object.values(D.geometries).filter(g=>g.C!==null).map(g=>g.C);
const dAll=D.decomposition["all geometries"];
const PZCSPREAD=1000*(Math.max(...zs)-Math.min(...zs)), CSPREAD=100*(Math.max(...cs)-Math.min(...cs))/median(cs);
document.getElementById("s1").append(...[
  [`${PZCSPREAD.toFixed(0)} mV`,`零电荷点跨度（${zs.length} 个能被采样点夹住的几何）`],
  [`${CSPREAD.toFixed(0)}%`,"同一批几何的割线电容跨度"],
  [`${median(cs).toFixed(1)} µF/cm²`,"割线电容中位数"],
  [`${dAll["U=+0.2V"].scale_ratio.toFixed(1)}倍`,"U = +0.2 V 处，两项的尺度之比（非贡献率）"],
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
 fill("corr",`${rs.length} 个几何里有 ${neg} 个为负，中位 r = ${median(rs).toFixed(2)}`);
 const hs=T.map(v=>v.tf&&v.tf.half).filter(x=>x);
 fill("half",`半衰波长 ${median(hs).toFixed(1)} Å`);
 fill("ratio",`对比度比 ${(100*median(T.map(v=>v.ratio).filter(x=>x))).toFixed(0)}%`);
}
{const U=D.undercoord;
 if(U){
  const txt=`${U.n} 个同时含低配位区与平台区的几何里有 ${U.n_lower} 个，`
    +`中位低 ${Math.abs(U.median_pct).toFixed(1)}%，最多低 ${Math.abs(U.min_pct).toFixed(1)}%（${U.min_structure}）`;
  fill("uc_main",txt); fill("corr2",txt);
  const v=U.exceptions.filter(e=>/^V[0-9]/.test(e.structure)).length;
  fill("uc_exc",`共 ${U.exceptions.length} 个：`
    +U.exceptions.map(e=>`${e.structure}（${e.config}，高 ${e.pct.toFixed(1)}%）`).join("、")
    +`。其中 ${v} 个是空位结构，另外的是 8×8 稀疏胞里的小岛与小坑。`);}}
fill("ngeom2",String(NG));
fill("pzcspread2",`${PZCSPREAD.toFixed(0)} mV`);

/* ---- the three lead figures, and every number quoted beside them ---- */
(function(){
  const FD=D.findings; if(!FD) return;
  const L=FD.local, M=FD.move, R=FD.relax;
  const put=(id,src,cap,alt)=>{const n=document.getElementById(id); if(n)
    n.append(zoomable(src,cap,alt), el("div",{class:"caption"},cap), el("div",{class:"hint"},"点击放大"));};
  put("localfig",`${D.imgbase}gallery/_finding_local.png`,
      "左：A1-hcp 的区域平均阴离子富集，吸附原子区低于周围台面。右：Step-16x1 上金属响应与阴离子响应的剖面，"
      +"各自对自身整胞平均归一化；金属的峰在台阶边，阴离子的峰在下台面里几个埃处。",
      "区域富集柱状图与台阶剖面");
  put("movefig",`${D.imgbase}gallery/_finding_move.png`,
      "左、中：Step-8x2 与把一颗边缘 Au 移到脚部之后的结构，紫圈是原位点、绿圈是新位置。"
      +"右：两条 σ(U)，斜率接近而横向位置不同。",
      "同组成对照的结构与两条充电曲线");
  fill("fl_U",`${L.U.toFixed(4)} V`);
  fill("fl_Kad",`K = ${L.K_adatom.toFixed(4)}（占面积 ${(100*L.area_adatom).toFixed(0)}%）`);
  fill("fl_Kte",`K = ${L.K_terrace.toFixed(4)}`);
  fill("fl_sig",`A1-hcp 的整胞面电荷 ${L.sigma_A1hcp>0?"+":""}${L.sigma_A1hcp.toFixed(2)} μC/cm²，`
               +`同胞平板 T-4x4 是 ${L.sigma_T4x4>0?"+":""}${L.sigma_T4x4.toFixed(2)} μC/cm²`
               +`（U = ${L.U_A1hcp.toFixed(4)} 与 ${L.U_T4x4.toFixed(4)} V，接近但不完全相等）`);
  fill("fl_step",`${L.step} 上金属响应的峰在台阶边、达整胞平均的 ${L.metal_peak.toFixed(2)} 倍，`
                +`而阴离子响应的峰移到下台面里 ${Math.min(...L.offsets_A.map(Math.abs)).toFixed(1)}–`
                +`${Math.max(...L.offsets_A.map(Math.abs)).toFixed(1)} Å 处，只有 ${L.ion_peak.toFixed(2)} 倍`);
  const A="Step-8x2", B="Step-8x2_edge-vacancy_plus_foot-adatom";
  fill("mv_pzc",`从 ${M[A].U_pzc_mV>0?"+":""}${M[A].U_pzc_mV.toFixed(1)} mV 移到 `
               +`${M[B].U_pzc_mV>0?"+":""}${M[B].U_pzc_mV.toFixed(1)} mV，共 ${M.dU_pzc_mV.toFixed(1)} mV`);
  fill("mv_C",`只从 ${M[A].C.toFixed(2)} 变到 ${M[B].C.toFixed(2)} μF/cm²（${M.dC_pct>0?"+":""}`
             +`${M.dC_pct.toFixed(1)}%）`);
  fill("rx_n",String(R.n));
  fill("rx_med",`${R.median_abs_dU_mV.toFixed(1)} mV`);
  fill("rx_max",`${R.max_abs_dU_mV.toFixed(1)} mV`);
  fill("rx_top",R.top.slice(0,3).map(t=>`${t.structure} ${t.dU_mV>0?"+":""}${t.dU_mV.toFixed(1)} mV`).join("、"));
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
  // the +-0.5 V points are drawn faint and dotted: visible, but plainly not part of the frozen statistics
  for(const g of Object.values(D.geometries)){
    if(!g.ext||!g.ext.length) continue;
    const all=g.pts.concat(g.ext).sort((a,b)=>a[0]-b[0]);
    F.s.append(el("polyline",{points:all.map(p=>`${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join(" "),
      fill:"none",stroke:FC[g.f]||"#888","stroke-width":.9,"stroke-opacity":.22,"stroke-dasharray":"2 3"}));}
  for(const g of Object.values(D.geometries))
    F.s.append(el("polyline",{points:g.pts.map(p=>`${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join(" "),
      fill:"none",stroke:FC[g.f]||"#888","stroke-width":1.2,"stroke-opacity":.7,"stroke-linecap":"round"}));
  for(const u of [-0.2,0.2]) F.s.append(el("line",{x1:X(u),x2:X(u),y1:m.t,y2:h-m.b,
    stroke:CSS("--gold"),"stroke-width":1,"stroke-dasharray":"4 3","stroke-opacity":.7}));
  F.s.append(el("text",{x:X(0.2)+5,y:m.t+13,"font-size":10,fill:CSS("--gold")},"统计口径边界 ±0.2 V"));
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
  t.append(el("thead",{},el("tr",{},...["范围","n","零电荷点跨度 (mV)","电容跨度 (%)","σ 散布：零电荷点项","电容项","比值"]
    .map(h=>el("th",{},h)))));
  const tb=el("tbody");
  const zh=k=>k==="all geometries"?"全部几何":k.replace(/^within cell A=(\d+) A\^2 \((.*)\)$/,"同一胞内 A=$1 Å² ($2)");
  for(const [k,v] of Object.entries(D.decomposition)){
    const w=v["U=+0.2V"];
    tb.append(el("tr",{},el("td",{},zh(k)),el("td",{class:"num"},String(v.n)),
      el("td",{class:"num"},fmt(v.U_pzc_spread_mV,0)),el("td",{class:"num"},fmt(v.C_spread_pct,1)),
      el("td",{class:"num"},fmt(w.from_pzc_shift_uC_per_cm2,2)),el("td",{class:"num"},fmt(w.from_capacitance_uC_per_cm2,2)),
      el("td",{class:"num"},fmt(w.ratio_pzc_over_capacitance,1)+"倍")));}
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
    "相对同胞平整平台的零电荷点位移 (mV)"));
  F.s.append(el("text",{x:w-m.r+12,y:m.t-13,"font-size":10,fill:CSS("--gold")},"电容变化"));
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
    F.s.append(el("text",{x:X(1)+5,y:m.t-2,"font-size":10.5,fill:CSS("--teal")},"体相浓度 K = 1"));}
  else{
    F.s.append(el("text",{x:m.l,y:m.t-2,"font-size":10.5,fill:CSS("--teal")},
      `注意：横轴从 ${lo.toFixed(2)} 起，未含体相基线 K = 1`));}
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
        r.K.toFixed(3)+"  (占面积 "+(100*r.a).toFixed(0)+"%)"));});
    F.s.append(el("line",{x1:X(v.regions.all.K),x2:X(v.regions.all.K),y1:y0+5,y2:y0+RH-7,stroke:CSS("--ink"),
      "stroke-width":1.2,"stroke-dasharray":"3 2"}));});
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-3,"text-anchor":"middle","font-size":11,fill:CSS("--ink-2")},
    "KΩ = 该区域内阴离子浓度相对体相的倍数"));
  document.getElementById("c_regions").append(F.s);
  document.getElementById("leg3").append(...CLS.map(([c,col])=>el("span",{},el("i",{class:"sw",style:`background:${col}`}),CZH[c])),
    el("span",{},el("i",{class:"sw",style:`background:${CSS("--ink")}`}),"全胞平均"));
})();

/* ---- transfer function ---- */
(function(){
  const curves=Object.values(D.transmission).filter(v=>v.tf&&v.tf.b.length);
  if(!curves.length) return;
  const halves=curves.map(v=>v.tf.half).filter(x=>x);
  document.getElementById("s4").append(...[
    [`${median(halves).toFixed(1)} Å`,"振幅衰减一半的波长（各几何的中位）"],
    [`${curves.length}`,"两端电势齐全、参与分析的几何数"],
    ["2.94 Å","Au–Au 最近邻间距，作为尺度参照"],
    [`${(100*median(Object.values(D.transmission).map(v=>v.ratio).filter(x=>x))).toFixed(0)}%`,
     "金属横向对比度存活到阴离子图的比例"],
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
  F.s.append(el("text",{x:X(hm)+6,y:m.t+14,"font-size":11,fill:CSS("--gold")},`半传输 ≈ ${hm.toFixed(1)} Å`));
  F.s.append(el("line",{x1:X(2.94),x2:X(2.94),y1:m.t,y2:h-m.b,stroke:CSS("--red"),"stroke-width":1.2,"stroke-dasharray":"2 3"}));
  F.s.append(el("text",{x:X(2.94)+5,y:h-m.b-9,"font-size":10,fill:CSS("--red")},"Au–Au 2.94 Å"));
  F.s.append(el("text",{x:(m.l+w-m.r)/2,y:h-5,"text-anchor":"middle","font-size":11,fill:CSS("--ink-2")},
    "横向波长 λ  (Å，对数轴)"));
  F.s.append(el("text",{x:14,y:(m.t+h-m.b)/2,"text-anchor":"middle","font-size":11,fill:CSS("--ink-2"),
    transform:`rotate(-90 14 ${(m.t+h-m.b)/2})`},"归一化传输率（对数）"));
  document.getElementById("c_tf").append(F.s);
})();

/* ---- dOmega ---- */
(function(){
  const seen=new Set(), rows=[];
  for(const p of D.pairs){
    const a=p.a.split("__"), b=p.b.split("__");
    if(!["ideal","relaxed"].includes(a[1])||!["ideal","relaxed"].includes(b[1])) continue;
    const key=[a[0],b[0]].sort().join("|"); if(seen.has(key)) continue; seen.add(key);
    rows.push({l:`${a[0]} ${a[1]}  对  ${b[0]} ${b[1]}`,v:p.d_relative_Omega_eV,z:p.dU_pzc_mV});}
  rows.sort((x,y)=>Math.abs(y.v)-Math.abs(x.v));
  const R=rows.slice(0,14); if(!R.length) return;
  const H=31,w=880,m={l:310,r:74,t:32,b:38},h=m.t+m.b+R.length*H;
  const lim=Math.max(0.035,...R.map(r=>Math.abs(r.v)))*1.12;
  const X=x=>m.l+(x+lim)/(2*lim)*(w-m.l-m.r);
  const F=frame(w,h,m), sc=D.same_structure_pair_scale;
  F.s.append(el("rect",{x:X(-sc.max),y:m.t-8,width:X(sc.max)-X(-sc.max),height:h-m.b-m.t+8,
    fill:CSS("--ink-2"),"fill-opacity":.09}));
  F.s.append(el("text",{x:X(0),y:m.t-14,"text-anchor":"middle","font-size":10.5,fill:CSS("--ink-2")},
    `同一形貌内部：|ΔΩ| 最大 ${(1000*sc.max).toFixed(0)} meV`));
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
    "U = ±0.2 V 上的 Δ(Ωₐ − Ωᵦ)   (meV)"));
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
  bar.append(mk("all",`全部 ${Object.keys(D.structures).length} 种`),
    ...fams.map(f=>mk(f,`${FZH[f]} (${Object.values(D.structures).filter(s=>s.f===f).length})`)));
  function render(){
    gal.textContent="";
    for(const [k,v] of Object.entries(D.structures).filter(([k,v])=>active==="all"||v.f===active)
        .sort((a,b)=>FAM.indexOf(a[1].f)-FAM.indexOf(b[1].f)||a[0].localeCompare(b[0]))){
      const tags=[el("span",{class:"tag"},`${v.n} Au`),el("span",{class:"tag"},`面积 ${v.A} Å²`),
        el("span",{class:"tag"},FZH[v.f]||v.f)];
      if(v.cn.kink) tags.push(el("span",{class:"tag k"},`${v.cn.kink} 个 CN≤6`));
      if(v.cn.edge) tags.push(el("span",{class:"tag e"},`${v.cn.edge} 个 CN 7–8`));
      if(v.C!=null) tags.push(el("span",{class:"tag"},`C ${v.C.toFixed(1)} µF/cm²`));
      if(v.z!=null) tags.push(el("span",{class:"tag"},`U₀ ${v.z>0?"+":""}${(1000*v.z).toFixed(0)} mV`));
      if(v.dz!=null) tags.push(el("span",{class:"tag"},`ΔU₀ ${v.dz>0?"+":""}${v.dz.toFixed(0)} mV`));
      gal.append(el("figure",{class:"gcard",style:"margin:0"},
        el("div",{class:"ghead"},el("h3",{},k),el("span",{class:"small"},v.cfg==="ideal"?"理想构型":"弛豫构型")),
        zoomable(`${D.imgbase}gallery/${k}.png`,`${k} · 简笔轮廓、俯视图与侧视图`,
          `${k}：简笔轮廓、俯视图与侧视图，原子按配位数着色`),
        v.sch?el("div",{class:"hint",style:"padding-top:10px;padding-bottom:0"},"轮廓："+v.sch):null,
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
      zoomable(`${D.imgbase}maps/${k}.png`,`${k} · 阴离子空间图`,
        `${k}：阴离子过量、其随电势的变化、金属电子数变化与配位数图`),
      el("div",{class:"hint"},"点击放大")));}
})();
</script>
"""


# The Artifact platform wraps the file in its own document skeleton; GitHub Pages does not, so a standalone build
# has to supply it. These are the parts of that skeleton the page actually relies on.
SKELETON_HEAD = """<!doctype html>
<html lang="zh-Hans">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<meta name="color-scheme" content="light dark">
<meta name="description" content="33 种 Au(111) 形貌在三个电子化学势下的恒电势 DFT 充电响应，以及有多少传递到了离子。">
<style>
:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}
body{margin:0}
img{max-width:100%}
[hidden]{display:none!important}
</style>
"""
SKELETON_TAIL = "\n</body>\n</html>\n"


# --------------------------------------------------------------------------------------------------------
# The structure-generation lineage. Content lives here, as data, and is rendered to markup rather than to an
# image: it is almost entirely text, and a raster of it is unreadable once the page scales it to the column.
# Every operation is taken from the build script that actually wrote the POSCAR.
LINEAGE_ROUTES = [
    ("a", "路线一 · 切出 Au(111) 四层平板", "再做局部增删 / 重排", [
        ("切出四层 (111) 平板（层高 5.0 / 7.4 / 9.8 / 12.2&nbsp;&Aring;）；只改面内周期胞，表面本身不动",
         "平整 Au(111)", 3,
         "T-4x4（4×4）　Flat-8x2（8×2）　Flat-16x1（16×1）",
         "共同基底与尺寸／镜像间距参考，不是三种缺陷"),
        ("顶层<b>移走</b> 1 / 2 / 3 个 Au；或在三重空位上<b>加</b> 1 或 3 个 Au",
         "点缺陷 / 极小团簇", 6,
         "V1 · V2 · V3（移走的位点）<br>A1-fcc · A1-hcp（同一个加原子，落在不同空位）　A3（三原子团簇）",
         "V 是移走，A 是加上；A1-fcc 与 A1-hcp 只差落在哪个三重空位"),
        ("按 fcc 延续注册，在顶层加一条<b>有限宽的单层条带</b>（条带恒占胞长的一半）",
         "条带直台阶", 5,
         "Step-8x1 · Step-16x1 · Step-24x1（上下台面各 4 / 8 / 12 行）<br>"
         "Step-8x2 · Step-16x2（沿台阶方向周期 ×2）",
         "一次同时产生上台面、下台面与两条不等价边缘"),
        ("从直台阶出发：在一条边<b>多放一个</b>原子；或把一个边缘原子<b>移到</b>脚部空位",
         "拐角 / 边缘重排", 3,
         "Kink-edge1 · Kink-edge2（8×3 条带，edge1 / edge2 各加一个）<br>"
         "Step-8x2_edge-vacancy_plus_foot-adatom（总 Au 数不变）",
         "后者与 Step-8x2 原子数、原子序都相同，只是一颗 Au 换了位置"),
        ("在顶层之上<b>加</b> 7 或 19 个 Au；或把同样 7 个 Au <b>改排</b>成 3+4 两排",
         "单层岛", 4,
         "Island-7-compact · Island-19-8x8（改尺寸）　Island-7-elongated（改形状，原子数不变）<br>"
         "Island-7-8x8（同一缺陷，扩大周围周期胞）",
         "尺寸、形状、镜像间距是三组分开的对照"),
        ("只从顶层<b>移走</b> 7 或 19 个 Au；或把 7 个缺失位改成 3+4 沟槽",
         "单层坑", 4,
         "Pit-7-compact · Pit-19-8x8（改宽度）　Pit-7-trench（改边缘形状）<br>"
         "Pit-7-8x8（同一缺陷，扩大周围周期胞）",
         "坑深恒为一层，比较的是宽度与边缘形状"),
        ("把<b>整个顶层</b>平移到 hcp 注册；或在 16 个位点的长度内放 <b>17</b> 个顶层 Au",
         "重构相关堆垛", 2,
         "R1-hcp-terminated（整层改注册）<br>R2-stripe-wall（多一个原子，形成注册过渡带）",
         "R1 不是一个 hcp 位吸附原子；R2 不是完整鱼骨重构"),
        ("在台阶<b>脚部</b>加七原子岛；或把坑里移走的 Au <b>就地</b>用来堆岛",
         "复合形貌", 2,
         "C1-island-near-step（Step-8x4 + 7 个 Au，实际与台阶相连）<br>"
         "C2-island+pit（8×8 胞内移走 7 个、加回 7 个，总 Au 数不变）",
         "C1 的岛与台阶连通，不是孤立的岛"),
    ]),
    ("b", "路线二 · 按高指数晶面直接切割", "得到规则台阶面，<b>不是</b>在 Au(111) 上加条带", [
        ("由立方胞直接按 (hkl) 切割（<code>build_vicinal_fixed.py</code>）；"
         "面间距 a/(2&radic;(h²+k²+l²)) = 0.693 / 0.443 / 0.256&nbsp;&Aring;",
         "台面宽度系列", 3,
         "Au221 · Au332 · Au554<br>与 (111) 的名义夹角 15.79° → 10.02° → 5.77°（台面渐宽）",
         "同一条构建路线上的宽度系列"),
        ("由 <code>ase.build.fcc211</code> 直接切割，<b>与上面三者不是同一个生成器</b>",
         "另一类台阶环境", 1,
         "Au211（与 (111) 的名义夹角 19.47°）",
         "不在上面的宽度系列里"),
    ]),
]
LINEAGE_CFG = [
    (0, "理想几何 · ideal", "33 个几何。直接作为静态参考态计算，不弛豫。"),
    (0, "在公共参考电势 &mu;<sub>0</sub> 下弛豫 · relax / relaxed",
     "16 个几何。给出每个结构的局部参考态；下面三类都由它生成。"),
    (1, "随机位移 · pert05 / pert10", "16 + 16 个几何。可动原子随机位移 0.05 / 0.10&nbsp;&Aring;。"),
    (1, "集体变形 · coll", "14 个几何。顶层间距 &minus;3%、面内应变 +1%、台阶边缘弯曲 0.15&nbsp;&Aring;。"),
    (1, "路径构型 · path",
     "12 个几何。吸附原子过桥位、边缘原子脱离到脚部、拐角原子沿边移动、岛／坑边原子进出。"),
]


def lineage_html():
    out = ['<div class="lin">',
           '<div class="lin-root"><b>fcc Au 晶体</b><span>建构用 a<sub>0</sub> = 4.158&nbsp;&Aring;</span>'
           '<span>33 种结构由下面两条路线生成</span></div>']
    for rid, head, sub, rows in LINEAGE_ROUTES:
        out.append(f'<div class="lin-route" data-r="{rid}">')
        out.append(f'<div class="lin-rhead">{head}　<span>{sub}</span></div>')
        for op, fam, n, mem, note in rows:
            out.append('<div class="lin-row">'
                       f'<div class="lin-op">{op}</div>'
                       '<div class="lin-arrow">&rarr;</div>'
                       f'<div class="lin-fam"><b>{fam}</b><span>{n} 个结构</span></div>'
                       f'<div class="lin-mem">{mem}<em>{note}</em></div>'
                       '</div>')
        out.append('</div>')
    out.append('<div class="lin-route" data-r="c">')
    out.append('<div class="lin-rhead">一个建构好的几何，如何变成多个被计算的构型　'
               '<span>33 个人工建构，其余都是它们的像；合计 107 个几何</span></div>')
    for depth, name, note in LINEAGE_CFG:
        cls = "lin-cfg sub" if depth else "lin-cfg"
        out.append(f'<div class="{cls}"><h4>{name}</h4><p>{note}</p></div>')
    out.append('<div class="lin-cfg"><h4>再乘上电势</h4><p>同一个几何在若干 TARGETMU 下各算一个电子态，'
               'U&nbsp;=&nbsp;&mu;<sub>0</sub>&nbsp;&minus;&nbsp;&mu;<sub>e</sub>。'
               '&plusmn;0.2&nbsp;V 窗口 291 个电子态已完成；&plusmn;0.5&nbsp;V 扩展进行中，本页结论不使用。'
               '<b>并非每个几何都有五个电势点。</b></p></div>')
    out.append('</div>')
    return "\n".join(out)


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
    page = HTML.replace("__LINEAGE__", lineage_html())
    page = page.replace("__DATA__", json.dumps(data, separators=(",", ":")).replace("</", "<\\/"))
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
    fpath = f"{ROOT}/analysis/gallery/_findings.json"
    if os.path.exists(fpath): data["findings"] = json.load(open(fpath))
    else: print("WARNING: no _findings.json; run scripts/render_findings.py first")
    data["undercoord"] = undercoordinated_vs_terrace()
    mapped = sorted(f[:-4] for f in os.listdir(f"{ROOT}/analysis/maps")) if os.path.exists(f"{ROOT}/analysis/maps") else []
    data["mapped"] = [m for m in mapped if m in data["structures"]]
    out = a.out or f"{WEB}/index.html"
    n = build(data, standalone=a.standalone, imgbase=a.imgbase, out_path=out)
    print(f"{n/1024:.0f} kB -> {out}  ({len(data['structures'])} structures, {len(data['mapped'])} maps"
          f"{', standalone, imgbase=' + repr(a.imgbase) if a.standalone else ''})")


if __name__ == "__main__":
    main()
