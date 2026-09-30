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

.gal{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:20px}
.gal.wide{grid-template-columns:1fr;gap:26px}
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
@media (max-width:560px){.gal{grid-template-columns:1fr}section{padding-block:38px}}
</style>

<header class="top"><div class="wrap">
  <div class="eyebrow">恒电势 DFT · VASP + VASPsol++ · 1 M 隐式电解质</div>
  <h1>Au 表面形貌与充电响应</h1>
  <p class="lede" style="margin-top:14px">以 Au(111) 台面及相关邻晶面（Au(211)、(221)、(332)、(554)）为对象：表面形貌如何改变
  金属容纳电荷的能力，以及这种原子尺度的差别有多少真正传递到液相中的离子。</p>
  <div class="kicker" id="kicker"></div>
  <div class="note" style="margin-top:20px;max-width:none"><b>数据范围（所有图表共用同一口径）。</b><span id="scope"></span></div>
</div></header>

<nav class="sticky"><div class="wrap">
  <a href="#pzc">一 · 变的是零电荷点</a>
  <a href="#which">二 · 哪些形貌，往哪个方向</a>
  <a href="#where">三 · 阴离子究竟在哪里增多</a>
  <a href="#filter">四 · 电解质这个低通滤波器</a>
  <a href="#omega">五 · 电势与相对稳定性</a>
  <a href="#gallery">结构图谱</a>
  <a href="#maps">阴离子空间图</a>
  <a href="#method">方法与边界</a>
</div></nav>

<div class="wrap">

<section id="pzc"><div class="finding">
  <div><div class="eyebrow">结论一</div><h2>零电荷点的移动幅度，明显大于电容的变化幅度</h2></div>
  <p class="claim">在采样点能夹住 &sigma;&nbsp;=&nbsp;0 的几何里，零电荷点跨度约 <span data-n="pzcspread"></span>，
  而割线电容跨度约 <span data-n="cspread"></span>。在本样本和本电势区间内，零电荷位置变化所对应的电荷尺度，
  比割线电容变化所对应的尺度大约一个量级。</p>
  <div class="stats" id="s1"></div>
  <div class="card pad"><div class="chartbox"><div id="c_sigma"></div></div>
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
  &mu;<sub>e</sub>。同一个胞内部这项偏移对每一行都相同，所以按胞分组的那几块才是形貌效应。此外，全部
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
  &sigma; 是 +3.18 对 +2.22&nbsp;&micro;C/cm²。本页此前把方向写反了，已更正。</div>
</div></section>

<section id="where"><div class="finding">
  <div><div class="eyebrow">结论三</div><h2>整体充电变强，与某处离子变多，不是同一件事</h2></div>
  <p class="claim">按最近表面原子的配位数给阴离子过量分区后，低配位位点的区域平均富集<em>低于</em>周围平台：在全部 41 个
  同时含低配位区与平台区的理想／弛豫几何中有 37 个如此，中位低 1.5%，最多低 6.1%（A3 的三吸附原子）。四个例外全是
  空位边缘，超出量只有 0.1%–0.3%。两种积分口径给出相同的排序。</p>
  <p class="claim" style="border-left-color:var(--teal)">逐柱的空间相关进一步支持这一点：把金属的正电荷增量
  (&minus;&Delta;n<sub>e</sub>) 与阴离子增量 (&Delta;&Gamma;<sub>&minus;</sub>) 按柱子求相关，
  <span data-n="corr"></span>。也就是说，<b>阴离子增多的地方，往往不是金属正电荷增加最多的地方</b>。</p>
  <div class="card pad"><div class="chartbox"><div id="c_regions"></div></div>
    <div class="caption">边界相对窗口下的富集比
    K<sub>&Omega;</sub>&nbsp;=&nbsp;&int;n<sub>&minus;</sub>&thinsp;/&thinsp;(n<sub>b</sub>&int;S<sub>ion</sub>)，
    取各结构最正的那个采样电势。柱按配位类分组，虚线是该结构的全胞平均。此图每种结构只画一个代表几何（理想或弛豫），
    共 33 条；上文 41 个的计数则覆盖同一结构的理想与弛豫两种几何。</div>
    <div class="legend" id="leg3"></div></div>
  <div class="note"><b>此前的机制解释是错的，已撤回。</b>本页曾写"吸附原子压低局部零电荷点，于是同一电势下那一小块更不正，
  吸引的阴离子更少"。这有两层问题：一是方向反了，零电荷点负移对应同一电势下<b>更正</b>（见结论二的更正说明）；
  二是这里算出的零电荷点是<b>整个周期体系</b>的一个标量，不是"吸附原子那一小块的局部零电荷点"，不能未经局部定义就
  当成分区属性。
  <br><br>正确的表述是：<b>整体充电与局部离子分布并非简单对应。</b>即使整体零电荷点负移、同一电势下金属总正电荷更多，
  低配位原子对应的液相区域也可能没有更高的区域平均阴离子浓度。局部分布取决于自洽电势、离子可达性与几何分配三者共同作用。
  这比原来的错误解释更值得注意。
  <br><br>K 与 &Gamma; 回答的是不同问题，两者都给出：小区域可以浓度很高却只容纳很少额外离子，大区域可以只略微富集却贡献
  大部分总过量。图中同时给出各区面积占比，避免混淆。区域平均本身只能说明"该区平均浓度较低"；"热点在哪里"由上面的
  逐柱相关和下方的空间图支持。</div>
</div></section>

<section id="filter"><div class="finding">
  <div><div class="eyebrow">结论四</div><h2>离子响应里的短波空间起伏明显减弱</h2></div>
  <p class="claim">把金属电荷响应的每一个傅里叶模式与阴离子响应的同一模式相比，相对谱幅随波长单调下降：原子尺度
  （&lambda;&nbsp;&asymp;&nbsp;1&nbsp;&Aring;）只剩参考值的百分之零点几，缺陷尺度（&lambda;&nbsp;&asymp;&nbsp;9&nbsp;&Aring;）
  还有三成以上。这与空间屏蔽／平滑的图像相容。</p>
  <div class="stats" id="s4"></div>
  <div class="card pad"><div class="chartbox"><div id="c_tf"></div></div>
    <div class="caption">相对谱幅 |F<sub>ion</sub>(k)|&thinsp;/&thinsp;|F<sub>metal</sub>(k)| 对横向波长，汇总了所有
    两端电势齐全的几何。该比值联系的是两个不同量纲的量，所以每条曲线按<b>自身最长波长那一档</b>归一；因此长波端的
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
    <b>符号更正</b>：因 &part;&Omega;/&part;&mu;<sub>e</sub>&nbsp;=&nbsp;&minus;N<sub>e</sub> 且
    U&nbsp;=&nbsp;&mu;₀&minus;&mu;<sub>e</sub>，U&nbsp;=&nbsp;+w 对应<em>较低</em>的 &mu;，积分前不应再加负号；
    此前的版本多加了一次负号，把"谁被相对稳定"说反了，现已改正。</div></div>
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
  <div class="note" style="max-width:none"><b>简笔轮廓是怎么来的。</b>先从真实原子坐标算出表面高度场（每根柱子指派给最近的
  未被埋住的原子），只在高度分布真的出现 &ge;1.2&nbsp;&Aring; 断层的地方切一刀画粗轮廓；没有断层的面（四个邻晶面）改用
  平滑渐变加下坡箭头。侧面剪影是过特征中心的一条真实剖线，不是投影最大值——否则一个紧凑的岛会被拉成和它 footprint
  一样宽的台阶。等高线之前做过一次小尺度高斯平滑。
  <br><br><b>简笔轮廓用于辅助辨认形貌，不要用于定量读取缺陷的边界位置、宽度或峰位</b>：平滑会移动等高线，也会改变窄特征的
  形状。定量的东西请看上面的分区数值和下面的空间图。<br>
  高度均匀的结构（两种堆垛重构）在这里看起来就是平的，这是实话：它们与平板的差别在层序和面内配准，不在高度。</div>
  <div class="note" style="max-width:none"><b>颜色的含义，以及为什么同一层里会出现不同颜色。</b>
  颜色编码的是<b>配位数</b>（3.4&nbsp;&Aring; 内的 Au 近邻数），不是原子所在的层。同一层里配位数本来就会不同，这正是要看的信息：
  台阶脚那一排原子虽然和平台同高，但上层平台压在它旁边，近邻数达到 10 以上，因此显灰色；而普通平台原子是 9，显黄色。
  <br><br>亮度和描边编码的是另一件事：<b>是否属于表面</b>。判据与区域分析<b>完全一致</b>——一个原子若上方
  2.35&nbsp;&Aring; 内有三个以上更高的近邻就算被埋住，否则算表面。实色深描边的正是区域分析给它分配柱体的那批原子；
  淡色无描边的是被整层压住的次表面原子。台阶的上下两层平台都算表面，所以都是实色。</div>
  <div class="filters" id="filters"></div>
  <div class="gal" id="gal"></div>
</div></section>

<section id="maps"><div class="finding">
  <div><div class="eyebrow">阴离子空间图</div><h2>逐个结构的阴离子分布</h2></div>
  <p class="lede">四张一组：最正采样电势下的逐柱阴离子过量、它在采样窗口两端之间的变化、同一电势步长下金属自身的正电荷变化
  (&minus;&Delta;n<sub>e</sub>)，以及切分区域所依据的配位数图。按周期平铺以便看清重复。<b>点击可放大。</b></p>
  <div class="note" style="max-width:none"><b>怎么读色标。</b>每幅图有<b>自己的</b>色标范围，<b>不同图之间不能直接比较幅度</b>，
  只能比较图内的空间分布。全为同号的量用顺序色（深=大），跨正负的量用以 0 为中心的发散色。单位：
  &Gamma;<sub>&minus;</sub> 与 &Delta;&Gamma;<sub>&minus;</sub> 是每投影面积的离子数（&Aring;<sup>&minus;2</sup>）；
  金属一侧画的是<b>正电荷</b>变化 &minus;&Delta;n<sub>e</sub>（e/&Aring;<sup>2</sup>，正值=失去电子）；配位数是离散标度。</div>
  <div class="gal wide" id="mapgrid"></div>
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
      每一个跨结构数值都只用三个基准电势（&minus;5.1071 / &minus;4.9071 / &minus;4.7071&nbsp;eV），因为每个几何都齐全。
      &plusmn;0.5&nbsp;V 扩展仍在计算、覆盖不全，只在 &sigma;(U) 曲线上单独叠加显示，<b>不进入任何统计量</b>——
      否则会拿采样到 &plusmn;0.5&nbsp;V 的几何去和只采样到 &plusmn;0.2&nbsp;V 的几何比较，而且每完成一个任务
      页面上的数字就会变一次。页首的快照给出本页实际使用的态数、几何数与构建时间。</p></div>
    <div class="card pad"><h3>本轮更正过的地方</h3><p class="small" style="margin-top:8px">
      巨势积分的符号（此前多加一次负号，结论反向）；结论三的机制解释（零电荷点负移对应同一电势下更正，且整体零电荷点
      不能当作局部属性）；金属—离子相位比较的符号（改用 &minus;&Delta;n<sub>e</sub>）；区域划分把平整平台的一半
      误判为次表面；&sigma; 图纵轴被写死而裁掉曲线。原始数据与脚本都在仓库中可追溯。</p></div>
    <div class="card pad"><h3>已知限制</h3><p class="small" style="margin-top:8px">
      生产参数下的力带有约 0.02&nbsp;eV/&Aring; 每原子的 egg-box 误差，因此这里没有任何结论建立在细小的力差异上。
      扰动构型与集体形变是人为设计的采样，不是热力学系综，只用于给出敏感性范围而不做平均。路径像很稀疏，不等于最小能量路径。</p></div>
  </div>
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
function frame(w,h,m){return {s:el("svg",{viewBox:`0 0 ${w} ${h}`,width:w,height:h,role:"img"}),w,h,m};}

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
  `统计口径固定为三个基准电势 −5.1071 / −4.9071 / −4.7071 eV，每个几何都齐全，共 ${NPT} 个态、`
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
 fill("ratio",`对比度比 ${(100*median(T.map(v=>v.ratio).filter(x=>x))).toFixed(0)}%`);}

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
    mapped = sorted(f[:-4] for f in os.listdir(f"{ROOT}/analysis/maps")) if os.path.exists(f"{ROOT}/analysis/maps") else []
    data["mapped"] = [m for m in mapped if m in data["structures"]]
    out = a.out or f"{WEB}/index.html"
    n = build(data, standalone=a.standalone, imgbase=a.imgbase, out_path=out)
    print(f"{n/1024:.0f} kB -> {out}  ({len(data['structures'])} structures, {len(data['mapped'])} maps"
          f"{', standalone, imgbase=' + repr(a.imgbase) if a.standalone else ''})")


if __name__ == "__main__":
    main()
