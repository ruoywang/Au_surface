#!/usr/bin/env python3
"""Build the standalone HTML report (Chinese) from analysis/web/data.json, inlined because the artifact CSP blocks fetch.

Usage (from Au_Cl/):  scripts/pyrun.sh scripts/build_web_report.py
Output: analysis/web/index.html  (publish with the gallery and map PNGs as supporting files)
"""
import json
import os

ROOT = "/anvil/scratch/x-rywang/Au_Cl"
WEB = f"{ROOT}/analysis/web"

HTML = r"""<title>Au(111) 充电图谱</title>
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
.chartbox{overflow-x:auto}
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
  <h1>Au(111) 充电图谱</h1>
  <p class="lede" style="margin-top:14px">表面形貌如何改变金属容纳电荷的能力，以及这种原子尺度的差别有多少真正传递到液相中的离子。
  33 种 Au(111) 结构、107 个独立几何、291 个在三个电子化学势下收敛的恒电势态。</p>
  <div class="kicker" id="kicker"></div>
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
  <div><div class="eyebrow">结论一</div><h2>形貌主要平移零电荷点，几乎不改变电容</h2></div>
  <p class="claim">在 104 个采样点能夹住 &sigma;&nbsp;=&nbsp;0 的几何里，零电荷点的跨度是 254&nbsp;mV，而割线电容的跨度只有
  15%。同一电势下两种形貌之间的电荷差异，约有八成来自零电荷点的平移，而不是电容的差别。</p>
  <div class="stats" id="s1"></div>
  <div class="card pad"><div class="chartbox"><div id="c_sigma"></div></div>
    <div class="caption">每个几何的表面电荷密度对内部电势 U&nbsp;=&nbsp;&mu;<sub>0</sub>&nbsp;&minus;&nbsp;&mu;<sub>e</sub>
    的曲线，全部取自实际收敛的 &mu;<sub>e</sub> 与 N<sub>e</sub>。曲线几乎彼此平行：族与族之间的差别是横向平移，不是斜率变化。</div>
    <div class="legend" id="leg1"></div></div>
  <div class="note"><b>为什么必须把两者分开。</b>某个缺陷在给定电势下带更多正电，并不等于它更容易被极化。把
  &sigma;(U)&nbsp;=&nbsp;C&thinsp;(U&nbsp;&minus;&nbsp;U<sub>pzc</sub>) 拆开，跨几何的 &sigma; 散布可分为
  &lang;C&rang;&thinsp;&times;&thinsp;零电荷点跨度 与 电容跨度&thinsp;&times;&thinsp;|U&nbsp;&minus;&nbsp;&lang;U<sub>pzc</sub>&rang;|
  两项，在整个采样窗口里第一项都占主导。</div>
  <div class="tablewrap"><table id="t_decomp"></table></div>
  <div class="caption">跨胞比较带有数值偏移：审计记录过两个<em>平板</em>胞仅因胞形、k 点网格与 PREC 就相差 72&nbsp;meV 的中性
  &mu;<sub>e</sub>。同一个胞内部这项偏移对每一行都相同，所以按胞分组的那几块是纯粹的形貌效应。</div>
</div></section>

<section id="which"><div class="finding">
  <div><div class="eyebrow">结论二</div><h2>凸起把零电荷点推负，空位把它推正</h2></div>
  <p class="claim">以<em>同一个胞</em>里的平整平台为基准：吸附原子与台阶把零电荷点压低 29&ndash;97&nbsp;mV，空位则抬高
  8&ndash;13&nbsp;mV，两种堆垛重构几乎不动。符号符合经典的功函数论证——凸起抹平电子溢出、降低功函数，凹陷相反。</p>
  <div class="card pad"><div class="chartbox"><div id="c_dpzc"></div></div>
    <div class="caption">只取理想构型，各自对照同一个胞里的平板成员。柱长是零电荷点的位移，右侧圆点标的是同一对照下的电容变化。</div></div>
  <div class="note"><b>电容同向变化，但幅度小得多。</b>加上去的原子比挖掉的影响大：三个吸附原子 +15%，单个吸附原子 +8%，
  一条台阶 +8%，而单个空位只有 +1%，两种重构在 1% 以内。凸向电解质的粗糙度既压低零电荷点又抬高电容，但在给定电势下真正
  改变电荷量的是前者。</div>
</div></section>

<section id="where"><div class="finding">
  <div><div class="eyebrow">结论三</div><h2>阴离子的热点并不在缺陷正上方</h2></div>
  <p class="claim">按最近表面原子的配位数给阴离子过量分区后，低配位位点的富集在 41 个几何里有 36 个<em>低于</em>周围平台，
  中位低 1.7%，最多低 6.1%。五个例外全是空位边缘，那是凹陷而不是凸起。两种积分口径给出相同的排序。</p>
  <div class="card pad"><div class="chartbox"><div id="c_regions"></div></div>
    <div class="caption">边界相对窗口下的富集比
    K<sub>&Omega;</sub>&nbsp;=&nbsp;&int;n<sub>&minus;</sub>&thinsp;/&thinsp;(n<sub>b</sub>&int;S<sub>ion</sub>)，
    取各结构最正的那个采样电势。柱按配位类分组，虚线是该结构的全胞平均。</div>
    <div class="legend" id="leg3"></div></div>
  <div class="note"><b>这就是充电结论从液相一侧看到的样子。</b>吸附原子压低了局部零电荷点，于是在固定的电极电势下，那一小块
  表面相对周围平台更不正，吸引到的阴离子反而更少。空位把零电荷点抬高，而空位边缘正是打破这一趋势的那几个几何。结论二和结论三
  是同一件物理，从界面的两侧各测了一次。
  <br><br>K 与 &Gamma; 回答的是不同问题，两者都给出：小区域可以浓度很高却只容纳很少额外离子，大区域可以只略微富集却贡献
  大部分总过量。表中同时给出各区面积占比，避免混淆。</div>
</div></section>

<section id="filter"><div class="finding">
  <div><div class="eyebrow">结论四</div><h2>电解质是一个低通滤波器，半传输波长约 8&nbsp;&Aring;</h2></div>
  <p class="claim">把金属电荷响应的每一个傅里叶模式与阴离子响应的同一模式相比，传输率从长波处的完全传递一路降到原子尺度的
  1% 以下。振幅在 &lambda;&nbsp;&asymp;&nbsp;7.7&nbsp;&Aring; 处衰减一半，约合两个半最近邻间距；金属横向对比度最终只有
  16% 存活到阴离子图上。</p>
  <div class="stats" id="s4"></div>
  <div class="card pad"><div class="chartbox"><div id="c_tf"></div></div>
    <div class="caption">归一化传输率 |F<sub>ion</sub>(k)|&thinsp;/&thinsp;|F<sub>metal</sub>(k)| 对横向波长，汇总了所有
    两端电势齐全的几何。该比值联系的是两个不同量纲的量，所以每条曲线按自身最长波长的一档归一，只有形状有意义；阴影带是四分位距。</div></div>
  <div class="note"><b>什么能传过去，什么传不过去。</b>Au&ndash;Au 间距是 2.94&nbsp;&Aring;：金属电荷的原子级起伏传到离子
  处只剩不到 2% 的相对振幅。缺陷尺度的特征——约 10&nbsp;&Aring; 的拐角重复周期、同量级宽度的岛——能传过去三分之一以上。
  线性化 Poisson&ndash;Boltzmann 平板配上 1&nbsp;M 的德拜长度给出同样的 k 依赖趋势，但实测曲线在短波端更陡，那一段由
  4&nbsp;&Aring; 的离子排斥空腔而非德拜屏蔽决定间距。这个理论式只作趋势参照，不是拟合。</div>
</div></section>

<section id="omega"><div class="finding">
  <div><div class="eyebrow">结论五</div><h2>在 &plusmn;0.2&nbsp;V 范围内，电势并没有把这些形貌重新排序</h2></div>
  <p class="claim">对电子数曲线积分，可以得到电势诱导的相对巨势变化，全程不需要把两个总能相减。在采样窗口内，跨形貌的最大值只有
  29&nbsp;meV，与同一形貌自身各采样构型之间的散布同量级。</p>
  <div class="card pad"><div class="chartbox"><div id="c_omega"></div></div>
    <div class="caption">同成分、同胞、不同结构的配对在 U&nbsp;=&nbsp;&plusmn;0.2&nbsp;V 上的
    &Delta;(&Omega;<sub>A</sub>&nbsp;&minus;&nbsp;&Omega;<sub>B</sub>)。负值表示电势往正方向移动时 A 相对被稳定。
    灰带是同一形貌自身采样构型之间配对的散布范围。</div></div>
  <div class="note"><b>这是什么，不是什么。</b>它只是电势<em>诱导</em>的那一部分。它不回答在参考电势下谁更稳定，那需要数据集
  尚未确定的统一能量基准；它也从不跨越不同的 Au 原子数，那需要引入储库项。&plusmn;0.5&nbsp;V 扩展会把窗口放大 2.5 倍，
  这些数值将同比放大。</div>
</div></section>

<section id="gallery"><div class="finding">
  <div><div class="eyebrow">结构图谱</div><h2>全部三十三种结构，按配位数着色</h2></div>
  <p class="lede">每张图用的都是实际参与计算的几何，不是重新生成的理想结构。着色分类与上文切分阴离子区域所用的分类完全一致，
  所以图与图可以直接对读。<b>点击任意图片可放大。</b></p>
  <div class="note" style="max-width:none"><b>颜色的含义，以及为什么同一层里会出现不同颜色。</b>
  颜色编码的是<b>配位数</b>（3.4&nbsp;&Aring; 内的 Au 近邻数），不是原子所在的层。同一层里配位数本来就会不同，这正是要看的信息：
  台阶脚那一排原子虽然和平台同高，但上层平台压在它旁边，近邻数达到 10 以上，因此显灰色；而普通平台原子是 9，显黄色。
  <br><br>亮度和描边编码的是另一件事：<b>是否属于表面</b>。实色深描边的原子在最外 3&nbsp;&Aring; 之内，也就是区域分析真正
  给它分配柱体的那批原子；淡色无描边的是更深的次表面原子。台阶的上下两层平台都算表面，所以都是实色。</div>
  <div class="filters" id="filters"></div>
  <div class="gal" id="gal"></div>
</div></section>

<section id="maps"><div class="finding">
  <div><div class="eyebrow">阴离子空间图</div><h2>逐个结构的阴离子分布</h2></div>
  <p class="lede">四张一组：最正采样电势下的逐柱阴离子过量、它在采样窗口两端之间的变化、同一电势步长下金属自身的电子数变化，
  以及切分区域所依据的配位数图。按周期平铺以便看清重复。<b>点击可放大。</b></p>
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

/* ---- kicker ---- */
const NG=Object.keys(D.geometries).length, NS=Object.keys(D.structures).length;
const NPT=Object.values(D.geometries).reduce((a,g)=>a+g.pts.length,0);
document.getElementById("kicker").append(
  ...[[NS,"种结构"],[NG,"个独立几何"],[NPT,"个收敛态"],["3","个电势：−5.1071 / −4.9071 / −4.7071 eV"]]
     .map(([v,l])=>el("span",{class:"chip"},`${v} ${l}`)));
document.getElementById("stamp").textContent = " · 其中 " + Object.keys(D.transmission).length +
  " 个几何具备两端电势的场分析。";

/* ---- Finding 1 stats ---- */
const zs=Object.values(D.geometries).filter(g=>g.z!==null).map(g=>g.z);
const cs=Object.values(D.geometries).filter(g=>g.C!==null).map(g=>g.C);
const dAll=D.decomposition["all geometries"];
document.getElementById("s1").append(...[
  [`${(1000*(Math.max(...zs)-Math.min(...zs))).toFixed(0)} mV`,"全部几何的零电荷点跨度"],
  [`${(100*(Math.max(...cs)-Math.min(...cs))/median(cs)).toFixed(0)}%`,"同一批几何的割线电容跨度"],
  [`${median(cs).toFixed(1)} µF/cm²`,"割线电容中位数"],
  [`${dAll["U=+0.2V"].ratio_pzc_over_capacitance.toFixed(1)}倍`,"U = +0.2 V 处，零电荷点项对电容项的比值"],
].map(([v,l])=>el("div",{class:"stat"},el("div",{class:"v"},v),el("div",{class:"l"},l))));

/* ---- sigma(U) ---- */
(function(){
  const F=frame(880,420,{l:56,r:14,t:14,b:42}), m=F.m, {w,h}=F;
  const X=x=>m.l+(x+0.56)/1.12*(w-m.l-m.r), Y=y=>h-m.b-(y+7)/14*(h-m.t-m.b);
  for(const t of [-0.5,-0.25,0,0.25,0.5]){
    F.s.append(el("line",{x1:X(t),x2:X(t),y1:m.t,y2:h-m.b,stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:X(t),y:h-m.b+16,"text-anchor":"middle","font-size":11,fill:CSS("--muted")},String(t)));}
  for(const t of [-6,-4,-2,0,2,4,6]){
    F.s.append(el("line",{x1:m.l,x2:w-m.r,y1:Y(t),y2:Y(t),stroke:CSS("--line"),"stroke-width":1}));
    F.s.append(el("text",{x:m.l-8,y:Y(t)+3.8,"text-anchor":"end","font-size":11,fill:CSS("--muted")},String(t)));}
  F.s.append(el("line",{x1:X(-0.56),x2:X(0.56),y1:Y(0),y2:Y(0),stroke:CSS("--line-2"),"stroke-width":1.4}));
  for(const g of Object.values(D.geometries))
    F.s.append(el("polyline",{points:g.pts.map(p=>`${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join(" "),
      fill:"none",stroke:FC[g.f]||"#888","stroke-width":1.1,"stroke-opacity":.62,"stroke-linecap":"round"}));
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
        zoomable(`gallery/${k}.png`,`${k} · 俯视图与侧视图`,`${k}：俯视图与侧视图，原子按配位数着色`),
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
      zoomable(`maps/${k}.png`,`${k} · 阴离子空间图`,
        `${k}：阴离子过量、其随电势的变化、金属电子数变化与配位数图`),
      el("div",{class:"hint"},"点击放大")));}
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
