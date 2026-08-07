# Milestone 8 — Diffusione del pool vibro-elettronico (κ_ve)

Tassello di fisica rimandato fin da **M3** e ricomparso come limite in ogni
caso multi-D (M5/M6 cono, M7 cilindro): la EveEqn del solver conduceva solo
convezione + rilassamento Landau-Teller + accoppiamento chimico, **senza il
termine di conduzione termica del modo vibro-elettronico**. Conseguenza
pratica: il flusso termico a parete era privo del contributo vibrazionale,
e la parete vibrazionale non poteva avere una BC con gradiente (un
`fixedValue` Tve senza diffusione non ha modo di propagarsi e innesca una
scacchiera — è esattamente ciò che avevamo visto in M7).

## 1. Il termine, e perché è semplicemente `laplacian(μ, eve)`

Nel modello two-temperature il modo vibro-elettronico conduce calore con
la sua conduttività κ_ve. Il flusso è q_ve = −κ_ve ∇Tve e nell'equazione
dell'energia vibrazionale compare −∇·q_ve = ∇·(κ_ve ∇Tve).

Con la **Eucken modificata per i modi interni**, κ_ve = μ·cv_ve. E poiché
eve = e_ve(Tve) è funzione della sola Tve (a composizione fissata),
puntualmente ∇eve = cv_ve ∇Tve. Quindi:

    κ_ve ∇Tve = μ · cv_ve ∇Tve = μ ∇eve      (esatto, non un'approssimazione)

e la conduzione vibrazionale è **∇·(μ ∇eve) = laplacian(μ, eve)**, implicita
nella variabile risolta eve (quindi stabile). È la stessa μ molecolare con
cui la EEqn principale conduce l'energia totale via
`thermophysicalTransport->divq(e)` — coerente.

Implementazione ([thermophysicalPredictor.C](../applications/modules/shockThermo/thermophysicalPredictor.C)), sotto lo stesso
guard `if (!inviscid)` della EEqn:

```cpp
if (!inviscid)
{
    EveEqn -= fvm::laplacian(thermo.mu(), eve);
}
```

Il segno rispecchia `EEqn += thermophysicalTransport->divq(e)` (divq porta
già dentro il −laplacian; qui `fvm::laplacian` è +laplacian, quindi `-=`).

## 2. Validazione: no-op dove deve, stabile e regressione-safe dove agisce

| Test | Effetto atteso | Risultato |
|---|---|---|
| **0D heat bath** (1 cella) | `laplacian` ≡ 0 → no-op esatto | T_tr max 0.02%, T_ve max 1.87% — **identici** ai valori M1/M3 |
| **1D relaxing** (gradienti reali, `inviscid=false`) | attivo ma piccolo (Peclet alto nella zona di rilassamento) | stabile; T max 1.03%/media 0.69%, Tve media 0.83% — **identici** a M4 |

Il caso 0D è a cella singola, quindi il termine è nullo per costruzione: è
la prova che non ho rotto la fisica di base. Il caso 1D ha viscosità reale
(`inviscid` è false quando max(μ)>0, `shockFluid.C:150`) quindi il termine
gira davvero, ma nella zona di rilassamento la convezione domina e la
diffusione vibrazionale è piccola — le metriche restano quelle di M4. Il
termine è quindi implementato, attivo e innocuo dove deve esserlo.

## 3. Prossimo passo: il payoff a parete (2D)

Il beneficio vero è a parete, dove finora mancava. Con `zeroGradient` su
Tve il flusso vibrazionale a parete è nullo per costruzione (dTve/dn=0);
per averlo serve una **BC di parete con gradiente** su Tve — `fixedValue`
Tw o `smoluchowskiJumpT` come per T. Prima di M8 questa BC era impraticabile
(scacchiera senza diffusione); **ora la diffusione la sostiene**. Il test:
- rimettere la parete vibrazionale a `smoluchowskiJumpT` (Tw=1000 sul
  cilindro, 297 sul cono) come nel paper;
- rilanciare cono/cilindro sul cluster e verificare che (a) la scacchiera
  NON ritorni grazie alla diffusione, (b) il flusso termico a parete ora
  includa il contributo vibrazionale (Stanton/C_H più vicini al paper).

Questo chiude il cerchio con il gap di flusso termico rimasto aperto in
M6/M7. È lavoro 2D su cluster, quindi separato da questa parte solver.

## 4. Rimandati (invariati)

Fig 9 (aria 5 specie), CVDV-QK, Tve multiple per specie nel solver, mesh
fine convergente del cilindro (re-grading near-wall, vedi M7 §6).
