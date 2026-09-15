import html, os, subprocess
OUTDIR = "/home/danielesalvi/OpenFOAM/Aerospace-Module-2/explainations/dipendenze"
os.makedirs(OUTDIR, exist_ok=True)
FONT = "DejaVu Sans"
STY = {
 "prog":('#DCEBFA','#2E6DB4','solid'), "core":('#EFEFEF','#6E6E6E','solid'), "mpp":('#FDE8D0','#C66A00','solid'),
 "case":('#E3F4E1','#3C8D3C','solid'), "bin":('#FFF4B8','#A68A00','solid'), "build":('#E8E4F6','#5E4FA2','solid'),
 "env":('#E0F2F1','#00796B','solid'), "miss":('#FDE0E0','#C62828','dashed'), "unused":('#DCEBFA','#C62828','dashed'),
 "gen":('#F7F7F7','#3C8D3C','dashed'),
}
GROUP = {"prog":"progetto","unused":"progetto","core":"core","mpp":"mpp","case":"caso","gen":"caso","bin":"bin",
         "build":"build","env":"env","miss":None}
GTITLE = {"progetto":("Il tuo codice (src/, applications/)","#2E6DB4"), "core":("OpenFOAM-13 (core, non modificato)","#6E6E6E"),
          "mpp":("Mutation++ (thirdParty/)","#C66A00"), "caso":("Casi e test","#3C8D3C"), "bin":("Librerie ed eseguibili prodotti","#A68A00"),
          "build":("Script di compilazione","#5E4FA2"), "env":("Ambiente (source etc/bashrc)","#00796B")}
N = {}
def n(i, kind, title, path="", note=""): N[i] = (kind, title, path, note)

# ------------------------------------------------------------------ nodi (tutti verificati)
n("env_bashrc","env","etc/bashrc","","da sorgentare dopo OpenFOAM-13")
n("env_settings","env","etc/config.sh/settings","","POLIMI_SRC, POLIMI_MODULES")
n("env_mpp","env","etc/config.sh/mutationpp","","MPP_DIRECTORY, MPP_DATA_DIRECTORY")
n("env_wl","env","etc/codeTemplates/dynamicCode/","fluidMulticomponentThermo, psiThermo","whitelist dynamicCode (contiene rrho)")
n("env_user","env","~/.OpenFOAM/13/codeTemplates/","","copia fatta da etc/bashrc\ncercata PRIMA di $FOAM_ETC")
n("b_all","build","Allwmake","radice")
n("b_src","build","src/Allwmake","","wmake all thermophysicalModels")
n("b_tm","build","src/thermophysicalModels/Allwmake")
n("b_specie","build","specie/Allwmake","","solo wmakeLnInclude: nessuna libreria")
n("b_mc","build","multicomponentThermo/Make/","files + options")
n("b_app","build","applications/Allwmake","","wmake all modules")
n("b_mod","build","applications/modules/Allwmake")
n("b_sf","build","shockFluid/Make/","files + options")
n("b_st","build","shockThermo/Make/","files + options")
n("b_tp","build","thirdParty/Allwmake","makeMutationpp","riga commentata in Allwmake:\nsi lancia a mano")
n("b_n2","build","nonEqTTv/Allwmake  +  N2/Make/","","non chiamato da Allwmake")
n("lib_he","bin","libhighEnthalpyThermophysicalModels.so","$FOAM_USER_LIBBIN")
n("lib_sf","bin","libshockFluid.so","$FOAM_USER_LIBBIN")
n("lib_st","bin","libshockThermo.so","$FOAM_USER_LIBBIN")
n("exe_n2","bin","Test-N2","$FOAM_USER_APPBIN")
n("lib_mpp","bin","libmutation++.so","thirdParty/Mutationpp/install/lib/")
n("s_heH","prog","highEnthalpyMulticomponentThermo.H","src/.../highEnthalpyMulticomponentThermo/",
  "classe astratta: Tve(), eve(), computeSourceVT()\nclasse composite: campi Tve_, eve_")
n("s_heC","prog","highEnthalpyMulticomponentThermo.C","","New(): chiama basicThermo::New")
n("s_heS","prog","highEnthalpyMulticomponentThermos.C","","istanziazioni statiche forCoeffGases/forGases\n(janaf, eConst, hConst: MAI rrho)")
n("s_HE","prog","HighEnthalpyMulticomponentThermo.H","IL BRIDGE",
  "costruttore, computeSourceVT(), correct(),\ninitMutation()")
n("s_rrho","unused","rrhoThermo.H / I.H / .C","src/.../specie/thermo/rrho/","copia rinominata di janafThermo\nnon compilata da nessun target")
n("s_sf","prog","shockFluid.H / shockFluid.C","applications/modules/shockFluid/","copia modificata del core:\nfluidThermo invece di psiThermo")
n("s_stH","prog","shockThermo.H","applications/modules/shockThermo/")
n("s_stC","prog","shockThermo.C","","costruttore + preSolve()\nregistra 'shockThermo' nella tabella solver")
n("s_tp","prog","thermophysicalPredictor.C","","specie, energia e, equazione di eve")
n("s_rdt","prog","setRDeltaT.C","","passo locale (LTS)")
n("c_foamRun","core","foamRun.C","applications/solvers/foamRun/","ciclo: prePredictor, momentumPredictor,\nthermophysicalPredictor, pressureCorrector")
n("c_solverNew","core","solverNew.C","src/finiteVolume/solver/",'libs.open("lib" + solverName + ".so")')
n("c_solver","core","solver.H","src/finiteVolume/solver/")
n("c_fluidSolver","core","fluidSolver.H","applications/modules/fluidSolver/")
n("c_sfcore","core","shockFluid/*.C del core","applications/modules/shockFluid/",
  "correctDensity, fluxPredictor, momentumPredictor,\nmoveMesh, pressureCorrector, setRDeltaT,\nthermophysicalPredictor + derivedFvPatchFields")
n("c_basic","core","basicThermo.H","src/thermophysicalModels/basic/")
n("c_fluid","core","fluidThermo.H","src/thermophysicalModels/basic/")
n("c_psi","core","psiThermo.H","src/thermophysicalModels/basic/")
n("c_mc","core","multicomponentThermo.H","src/thermophysicalModels/multicomponentThermo/")
n("c_fmc","core","fluidMulticomponentThermo.H","src/thermophysicalModels/multicomponentThermo/")
n("c_BT","core","BasicThermo.H","template <MixtureType, BasicThermoType>")
n("c_PT","core","PsiThermo.H","template <Thermo>")
n("c_MT","core","MulticomponentThermo.H","template <BaseThermo>")
n("c_FMT","core","FluidMulticomponentThermo.H","template <BaseThermo>")
n("c_btNew","core","basicThermoTemplates.C","src/thermophysicalModels/basic/basicThermo/",
  "basicThermo::New:\n1) cerca il tipo nella tabella statica\n2) altrimenti dynamicCode, con il template\n    che ha il NOME DELLA CLASSE")
n("c_make","core","makeFluidMulticomponentThermo.H","+ makeThermo.H, forGases.H","macro che registrano i tipi nelle tabelle")
n("c_janaf","core","janafThermo.H / I.H / .C","src/thermophysicalModels/specie/thermo/janaf/")
n("c_comb","core","combustionModel.H / combustionModelNew.C","src/combustionModels/")
n("c_tr","core","header di trasporto e schemi","",
  "compressibleMomentumTransportModel.H\nfluidThermoThermophysicalTransportModel.H\nfluidMulticomponentThermophysicalTransportModel.H\nmultivariateScheme.H")
n("miss_tmpl","miss","codeTemplates/dynamicCode/highEnthalpyMulticomponentThermo","NON ESISTE",
  "ne' nel progetto, ne' in ~/.OpenFOAM, ne' nel core")
n("m_h","mpp","mutation++.h","thirdParty/Mutationpp/install/include/mutation++/")
n("m_data","mpp","nonEqTTv/mutation-data/","","mixtures/air_5.xml, thermo/species.xml,\nthermo/elements.xml, transfer/VT.xml, ...")
n("k_allrun","case","solverHeatBath/Allrun","applications/test/nonEqTTv/")
n("k_ctrl","case","system/controlDict","","solver shockThermo; deltaT 1e-9; endTime 2e-5\nprobes T Tve p")
n("k_bmd","case","system/blockMeshDict","","una sola cella")
n("k_fv","case","system/fvSchemes, fvSolution","","div(phi,Yi_h); solver (Tve|eve).*, Yi.*")
n("k_pp","case","constant/physicalProperties","+ #include speciesThermo.janaf",
  "thermoType: highEnthalpyThermo ... thermo rrho\nhighEnthalpyMutation { mixture air_5 ... }")
n("k_mt","case","constant/momentumTransport","","laminar")
n("k_0","case","0/T  0/p  0/U  0/Tve  0/N2  0/Ydefault","","T=10000, Tve=1000, p=101325, U=0")
n("miss_comb","miss","constant/combustionProperties","ASSENTE","-> noCombustion")
n("k_cmp","case","compare-heatbath.py","","errori massimi + grafico")
n("k_ref","gen","reference.csv","generato da Allrun")
n("k_probes","gen","postProcessing/probes/0/","generato da foamRun")
n("t_n2","case","Test-N2.C","applications/test/nonEqTTv/N2/","heat bath 0D con solo Mutation++")
n("t_out","gen","output/results-N2.csv","generato da Test-N2")
n("tut","case","tutorials/shockThermo/shockTube/","","stesso heat bath, thermo janaf")
n("x_park","unused","Test-thermoMixturePark2T.C","applications/test/thermoMixturePark2T/",
  "include RRHOThermo.H che non esiste\nnon presente in nessun Allwmake")

E = {
 "inc": 'color="#333333", arrowhead=normal', "inh": 'color="#8E24AA", arrowhead=empty, penwidth=2.2',
 "bld": 'color="#1565C0", style=dashed, arrowhead=vee, penwidth=1.4', "run": 'color="#E65100", penwidth=2.4, arrowhead=normal',
 "io": 'color="#2E7D32", style=dotted, penwidth=2, arrowhead=open', "env": 'color="#00796B", style=dashed, penwidth=1.6, arrowhead=odiamond',
 "copy": 'color="#9E9E9E", style=dashed, arrowhead=onormal, penwidth=1.4',
}

def lab(title, path, note):
    t = '<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" CELLPADDING="1">'
    t += '<TR><TD ALIGN="LEFT"><B>%s</B></TD></TR>' % html.escape(title)
    if path: t += '<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="10" COLOR="#555555">%s</FONT></TD></TR>' % html.escape(path)
    for line in (note.split("\n") if note else []):
        t += '<TR><TD ALIGN="LEFT"><FONT POINT-SIZE="11">%s</FONT></TD></TR>' % html.escape(line)
    return t + '</TABLE>>'

def panel(fname, title, edges, rankdir="TB", extra_nodes=(), clusters=True):
    used = []
    for a, b, *_ in edges:
        for x in (a, b):
            if x not in used: used.append(x)
    for x in extra_nodes:
        if x not in used: used.append(x)
    L = ['digraph G {',
         '  graph [rankdir=%s, fontname="%s", labelloc=t, fontsize=30, pad=0.5, nodesep=0.45, ranksep=0.8, '
         'label=<<B>%s</B><BR/> >];' % (rankdir, FONT, html.escape(title)),
         '  node [fontname="%s", fontsize=13, margin="0.14,0.07"];' % FONT,
         '  edge [fontname="%s", fontsize=12];' % FONT]
    groups = {}
    for x in used:
        g = GROUP[N[x][0]]
        groups.setdefault(g, []).append(x)
    def nodeline(x):
        kind, t, p, no = N[x]; f, c, s = STY[kind]
        shape = "component" if kind == "bin" else "box"
        return '    %s [shape=%s, style="rounded,filled,%s", fillcolor="%s", color="%s", penwidth=1.8, label=%s];' % (
            x, shape, s, f, c, lab(t, p, no))
    for g, xs in groups.items():
        if g is None or not clusters:
            L += [nodeline(x) for x in xs]; continue
        gt, gc = GTITLE[g]
        L.append('  subgraph cluster_%s { label=<<B>%s</B>>; fontsize=17; style="rounded"; color="%s"; penwidth=2; margin=16;' % (g, html.escape(gt), gc))
        L += [nodeline(x) for x in xs]
        L.append('  }')
    for e in edges:
        a, b, kind = e[0], e[1], e[2]
        label = e[3] if len(e) > 3 else ""
        extra = e[4] if len(e) > 4 else ""
        lb = ', label=" %s "' % label.replace('"', '\\"') if label else ""
        L.append('  %s -> %s [%s%s%s];' % (a, b, E[kind], lb, (", " + extra) if extra else ""))
    L.append('}')
    dot = os.path.join(OUTDIR, fname + ".dot")
    open(dot, "w").write("\n".join(L) + "\n")
    for fmt in ("png", "svg"):
        subprocess.run(["dot", "-T" + fmt, "-Gdpi=100", dot, "-o", os.path.join(OUTDIR, fname + "." + fmt)], check=True)
    return os.path.join(OUTDIR, fname + ".png")

# ------------------------------------------------------------------ pannello 1: compilazione
P1 = [
 ("env_bashrc","env_settings","env","sorgenta"), ("env_bashrc","env_mpp","env","sorgenta"),
 ("env_bashrc","env_user","env","copia etc/codeTemplates"), ("env_wl","env_user","copy","copiato in"),
 ("env_settings","b_mc","env","POLIMI_SRC"), ("env_settings","b_st","env","POLIMI_SRC, POLIMI_MODULES"),
 ("env_mpp","b_mc","env","MPP_DIRECTORY"), ("env_mpp","b_st","env","MPP_DIRECTORY"), ("env_mpp","b_n2","env","MPP_DIRECTORY"),
 ("b_all","b_src","bld","1"), ("b_all","b_app","bld","2"), ("b_all","b_tp","copy","commentato"),
 ("b_src","b_tm","bld"), ("b_tm","b_specie","bld"), ("b_tm","b_mc","bld"),
 ("b_specie","s_rrho","bld","lnInclude"),
 ("b_mc","s_heC","bld","compila"), ("b_mc","s_heS","bld","compila"), ("b_mc","lib_he","bld","produce"),
 ("b_app","b_mod","bld"), ("b_mod","b_sf","bld"), ("b_mod","b_st","bld"),
 ("b_sf","s_sf","bld","compila"), ("b_sf","c_sfcore","bld","compila dal core"), ("b_sf","lib_sf","bld","produce"),
 ("b_st","s_stC","bld","compila"), ("b_st","s_tp","bld","compila"), ("b_st","s_rdt","bld","compila"), ("b_st","lib_st","bld","produce"),
 ("lib_he","lib_mpp","bld","linka"), ("lib_st","lib_sf","bld","linka"), ("lib_st","lib_he","bld","linka"),
 ("b_n2","t_n2","bld","compila"), ("b_n2","exe_n2","bld","produce"), ("exe_n2","lib_mpp","bld","linka"),
 ("b_tp","lib_mpp","bld","produce"),
]
# ------------------------------------------------------------------ pannello 2: include ed ereditarieta'
P2 = [
 ("s_stC","s_stH","inc"), ("s_tp","s_stH","inc"), ("s_rdt","s_stH","inc"),
 ("s_stH","s_sf","inh","shockThermo : shockFluid"), ("s_sf","c_fluidSolver","inh","shockFluid : fluidSolver"),
 ("c_fluidSolver","c_solver","inh","fluidSolver : solver"),
 ("s_stH","s_heH","inc","#include"), ("s_stH","c_comb","inc"), ("s_stH","c_tr","inc"), ("s_sf","c_tr","inc"),
 ("s_heC","s_heH","inc"), ("s_heS","s_heH","inc"), ("s_heS","c_make","inc","#include"),
 ("s_heH","s_HE","inc","#include"), ("s_HE","m_h","inc","#include"), ("t_n2","m_h","inc","#include"),
 ("s_heH","c_psi","inh","virtual"), ("s_heH","c_fmc","inh","virtual"),
 ("c_psi","c_fluid","inh","virtual"), ("c_fmc","c_fluid","inh","virtual"), ("c_fmc","c_mc","inh","virtual"),
 ("c_fluid","c_basic","inh","virtual"), ("c_mc","c_basic","inh","virtual"),
 ("s_HE","c_FMT","inh","HighEnthalpyMulticomponentThermo<T>\n: FluidMulticomponentThermo<T>"),
 ("c_FMT","c_MT","inh","T = MulticomponentThermo<...>"), ("c_MT","c_PT","inh","BaseThermo = PsiThermo<...>"),
 ("c_PT","c_BT","inh","Thermo = BasicThermo<Mixture, composite>"),
 ("c_BT","s_heH","inh","BasicThermoType = composite\n: highEnthalpyMulticomponentThermo", "constraint=false"),
 ("s_rrho","c_janaf","copy","copia rinominata"), ("x_park","s_rrho","inc","#include"),
]
# ------------------------------------------------------------------ pannello 3: esecuzione del test 3a
P3 = [
 ("k_allrun","exe_n2","run","1  Test-N2 10000 1000"), ("exe_n2","m_data","io","legge"), ("exe_n2","t_out","io","scrive"),
 ("t_out","k_ref","io","copiato"), ("k_allrun","m_data","env","export MPP_DATA_DIRECTORY"),
 ("k_allrun","k_bmd","run","2  blockMesh"), ("k_allrun","c_foamRun","run","3  foamRun"),
 ("c_foamRun","k_ctrl","io","legge 'solver'"), ("c_foamRun","c_solverNew","run","4"), ("c_solverNew","lib_st","run","carica"),
 ("c_solverNew","s_stC","run","5  crea shockThermo"), ("s_stC","s_heC","run","6  ::New"),
 ("s_heC","c_btNew","run","7  basicThermo::New"), ("c_btNew","k_pp","io","legge thermoType"),
 ("c_btNew","s_heS","run","8a  janaf: nella tabella statica"),
 ("c_btNew","miss_tmpl","run","8b  rrho: non in tabella,\ntemplate assente -> FatalError"),
 ("s_heS","s_HE","run","9  costruisce il bridge"), ("s_HE","lib_mpp","run","10  initMutation"),
 ("s_HE","k_pp","io","legge highEnthalpyMutation"), ("lib_mpp","m_data","io","legge"),
 ("s_heH","k_0","io","0/Tve"), ("c_basic","k_0","io","0/T, 0/p"), ("c_mc","k_0","io","0/N2, 0/Ydefault"), ("s_sf","k_0","io","0/U"),
 ("s_sf","k_mt","io","legge"), ("s_stC","c_comb","run","combustionModel::New"), ("c_comb","miss_comb","io","cerca"),
 ("s_stC","s_rdt","run","preSolve: solo se LTS\n(test 3a: Euler, non usato)","style=dashed"),
 ("c_foamRun","c_sfcore","run","a ogni passo:\nflussi, velocita', pressione"), ("c_foamRun","s_tp","run","11  a ogni passo"),
 ("s_tp","s_HE","run","12  computeSourceVT(), correct()"), ("s_tp","k_fv","io","schemi e solver"),
 ("c_foamRun","k_probes","io","scrive"), ("k_allrun","k_cmp","run","13  confronto"),
 ("k_cmp","k_ref","io","legge"), ("k_cmp","k_probes","io","legge"),
 ("tut","s_heS","run","thermo janaf: tabella statica","style=dashed"),
]
p1 = panel("pannello-1-compilazione", "1 · Come si compila: ambiente, script, librerie", P1)
p2 = panel("pannello-2-classi", "2 · #include ed ereditarieta' delle classi", P2)
p3 = panel("pannello-3-esecuzione-test-3a", "3 · Cosa succede quando lanci il test 3a (solverHeatBath/Allrun)", P3, clusters=False)

# ------------------------------------------------------------------ legenda + titolo
def cell(kind, text):
    f, c, st = STY[kind]
    return '<TD BGCOLOR="%s" BORDER="1" COLOR="%s" ALIGN="LEFT">%s</TD>' % (f, c, html.escape(text))
rows = [("prog","tuo codice"),("core","OpenFOAM-13 (core)"),("mpp","Mutation++"),("case","caso / test"),("bin","libreria o eseguibile"),
        ("build","script di compilazione"),("env","ambiente"),("gen","file generato"),("miss","file ASSENTE"),("unused","presente ma non compilato")]
arrows = [("#333333","━━▶","#include"),("#8E24AA","━━▷","eredita da (l'#include e' sottinteso)"),("#1565C0","╍╍▶","compila / linka / produce"),
          ("#E65100","━━▶","chiama o carica a runtime (numeri = ordine)"),("#2E7D32","┈┈▷","legge o scrive un file"),
          ("#00796B","╍╍◇","variabili d'ambiente"),("#9E9E9E","╍╍▷","copia / commentato")]
t = '<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="5" CELLPADDING="5">'
t += '<TR><TD COLSPAN="2" ALIGN="LEFT"><FONT POINT-SIZE="34"><B>Aerospace-Module-2 (branch revised): dipendenze fra i file</B></FONT></TD></TR>'
t += '<TR><TD COLSPAN="2" ALIGN="LEFT"><FONT POINT-SIZE="16">Verificato sui sorgenti del progetto, sui sorgenti di OpenFOAM-13 e sul Doxygen cpp.openfoam.org/v13</FONT></TD></TR>'
for i in range(max(len(rows), len(arrows))):
    left = cell(*rows[i]) if i < len(rows) else '<TD></TD>'
    right = ('<TD ALIGN="LEFT"><FONT COLOR="%s"><B>%s</B></FONT>  %s</TD>' % (arrows[i][0], arrows[i][1], html.escape(arrows[i][2]))) if i < len(arrows) else '<TD></TD>'
    t += '<TR>%s%s</TR>' % (left, right)
t += '</TABLE>>'
ld = os.path.join(OUTDIR, "legenda.dot")
open(ld, "w").write('digraph L { graph [pad=0.4]; node [fontname="%s", fontsize=15]; legend [shape=plaintext, label=%s]; }\n' % (FONT, t))
subprocess.run(["dot", "-Tpng", "-Gdpi=100", ld, "-o", os.path.join(OUTDIR, "legenda.png")], check=True)

# ------------------------------------------------------------------ montaggio in un'unica immagine
from PIL import Image
imgs = [Image.open(os.path.join(OUTDIR, "legenda.png"))] + [Image.open(x) for x in (p1, p2, p3)]
for x, im in zip(["legenda", "p1", "p2", "p3"], imgs): print(x, im.size)
W = max(im.width for im in imgs); gap = 40
H = sum(im.height for im in imgs) + gap * (len(imgs) + 1)
canvas = Image.new("RGB", (W + 2 * gap, H), "white")
y = gap
for im in imgs:
    canvas.paste(im.convert("RGB"), ((W + 2 * gap - im.width) // 2, y)); y += im.height + gap
out = "/home/danielesalvi/OpenFOAM/Aerospace-Module-2/explainations/dipendenze-progetto.png"
canvas.save(out); print("finale", canvas.size, out)
