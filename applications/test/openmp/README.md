# Parallelizzazione OpenMP del thermo

I cicli su celle del ponte con Mutation++
(`src/thermophysicalModels/multicomponentThermo/highEnthalpyMulticomponentThermo/`)
sono distribuiti sui thread con OpenMP. Sono cinque, tutti con la stessa
struttura — una cella, uno stato del gas, nessuna comunicazione con le
vicine — e valgono la maggior parte del tempo di calcolo del solver:

| funzione | cosa calcola |
|---|---|
| `updateSources()` | in un giro: ω_s di tutte le specie (eq. 27), calore di reazione (eq. 1), sorgente di `eve` V-T + chimica (eq. 8 e 30) |
| `correct()` | dalle energie a `T` e `Tve`, Newton di Mutation++ |
| `updatePsi()` | `psi = rho/p`, celle e facce di bordo |

Le `computeSourceY/E/Ve` che il solver chiama non girano più sulle celle:
leggono quello che `updateSources()` ha messo da parte (vedi
[Il consolidamento delle sorgenti](#il-consolidamento-delle-sorgenti)).

## Come, e perché così

`Mutation::Mixture` **non è thread-safe**: `setState()` muta l'oggetto e tutto
quello che si chiede dopo (`T()`, `Tv()`, `netProductionRates()`, …) legge
quello stato. Ogni thread ha quindi la sua mixture, i suoi vibratori e i suoi
buffer, raccolti in un `MutationWorkspace` e costruiti una volta sola
all'avvio in `initMutation()`.

Non è bastato. La stessa `Mutation++` teneva i buffer del Newton che inverte
le energie in **variabili `static` di funzione**, una copia per processo
condivisa fra tutte le miscele e tutti i thread. Vedi
[`docker/mutationpp-thread-safety.patch`](../../../docker/mutationpp-thread-safety.patch):
senza quella patch questo test dà temperature che differiscono di ordini di
grandezza fra celle identiche, e cambiano da una corsa all'altra.

Altri dettagli che contano:

- **Numero di thread limitato dalla taglia del ciclo** (`minCellsPerThread`).
  Senza, i casi heat bath a una cella e le patch di bordo da due facce che
  `updatePsi()` percorre a ogni passo lancerebbero otto thread per niente, e i
  thread senza lavoro poi girano a vuoto: sul test `fig3a` a una cella il
  tempo di CPU passava da 16 s a 2 minuti.
- **`schedule(static)`**, non `runtime`: il default di libgomp per `runtime` è
  `dynamic` con fette da **una** cella, il peggiore dei mondi.
- **Buffer `thread_local`** in `mutationSources.H`: `speciesEv` e `speciesEve`
  allocavano due vettori sullo heap a ogni chiamata, cioè per ogni specie di
  ogni cella, mettendo i thread in coda sull'allocatore.
- **Le eccezioni non escono dalla regione parallela**: un'eccezione che scappa
  da un `omp for` termina il programma senza dire niente, e il Newton di
  Mutation++ può non convergere. `parallelFor()` le cattura e le rilancia come
  `FatalError` dopo la regione.

## Il test

```sh
./Allrun                          # 20000 celle, 200 passi, 1 2 4 8 thread
./Allrun 50000 1e-7 "1 4"         # celle, tempo finale, lista thread
```

Il caso è una fila di celle **tutte identiche**: ognuna è l'heat bath della
fig. 7 del paper (N2 + N reagente, il più caro perché accende anche la
chimica), il gas è fermo e il campo uniforme, quindi nessuna cella parla con
le vicine. Le patch `empty` su *y* e *z* impongono una cella sola in quelle
direzioni, per questo si raffina solo in *x*.

Due controlli distinti:

- **accordo fra numeri di thread diversi.** Non è atteso bit a bit:
  `solveEnergies()` è un Newton che parte dalle temperature correnti della
  miscela, cioè da quelle dell'ultima cella trattata da *quel* thread. Con un
  thread la catena è 0, 1, 2, …; con otto ogni thread ha la sua. Stessa
  soluzione a meno della tolleranza del Newton, 1e-12.
- **riproducibilità a parità di thread.** Questa sì, bit a bit: con
  `schedule(static)` la ripartizione è fissa. È il controllo che smaschera le
  corse critiche, ed è quello che ha trovato il bug di Mutation++.

## Risultati

MacBook Air M2, 8 core (4 performance + 4 efficiency), container arm64
nativo, 20000 celle, 100 passi. Migliore di tre ripetizioni alternate fra i
due casi, a macchina fredda:

| | 1 thread | 8 thread | speedup OpenMP |
|---|---:|---:|---:|
| prima del consolidamento | 10.80 s | 5.12 s | 2.11× |
| dopo il consolidamento | 4.88 s | 3.27 s | 1.49× |

Rispetto al punto di partenza (10.80 s, un thread, senza consolidamento):

| | tempo | speedup |
|---|---:|---:|
| solo OpenMP, 8 thread | 5.12 s | 2.11× |
| solo consolidamento, 1 thread | 4.88 s | 2.21× |
| consolidamento + 8 thread | 3.27 s | **3.30×** |

Accordo fra 1 e 8 thread: 1.7e-14 relativo. Riproducibile bit a bit.

> **Attenzione ai tempi su un portatile senza ventola.** Su questo MacBook Air
> lo stesso caso misura 4.88 s a freddo e 18.4 s dopo mezz'ora di carico
> sostenuto, e il caso a 8 thread degrada più di quello a 1 thread, quindi lo
> *speedup* apparente crolla (1.27× invece di 3.30×). Le misure qui sopra sono
> il migliore di tre ripetizioni con i due casi alternati, a macchina fredda;
> chi rifà il test deve fare lo stesso, altrimenti misura il dissipatore.

Dopo il consolidamento il pavimento seriale è circa 2.5 s su 4.88, cioè **più
della metà del tempo**: sono i solutori lineari di `Yi`, `e` ed `eve`, la
ricostruzione dei flussi e la scrittura dei campi, che OpenMP qui non tocca e
che crescono con le celle esattamente come la parte parallela. Aumentare la
mesh non sposta il rapporto.

Oltre i 4 thread si guadagna poco perché gli ultimi quattro finiscono sui core
*efficiency*, e con le fette fisse di `schedule(static)` tutti aspettano loro.
Passando a `schedule(dynamic, 64)` a 8 thread si recupera un 10-15%, ma si
perde la riproducibilità bit a bit: le fette vanno a chi si libera, quindi in
ordine diverso a ogni corsa. Per provarlo va cambiato il pragma in
`HighEnthalpyMulticomponentThermo.H` e ricompilato.

## Il consolidamento delle sorgenti

Prima esisteva un secondo problema, più grosso del parallelismo: ogni
`computeSource*` faceva il suo giro sulle celle e ricalcolava
`productionRates()`. Quattro volte per le specie risolte, una per il calore di
reazione, una dentro la sorgente di `eve`: **sei calcoli identici della stessa
cosa, il 65% del tempo dell'intera corsa.**

Misurato con timer sulle cinque funzioni (20000 celle, 100 passi, 1 thread,
`foamRun` 11.46 s in totale):

| funzione | chiamate | tempo | per chiamata |
|---|---:|---:|---:|
| `computeSourceY` | 800 (4 specie × 200) | 4.95 s | 6.18 ms |
| `computeSourceE` | 200 | 1.23 s | 6.16 ms |
| `computeSourceVe` | 200 | 2.01 s | 10.05 ms |
| `correct` (Newton) | 200 | 0.46 s | 2.32 ms |
| `updatePsi` | 201 | 0.12 s | 0.59 ms |
| **thermo** | | **8.77 s** | **77% del totale** |

`computeSourceE` costava 6.16 ms per chiamata e `computeSourceY` 6.18 ms:
identici, perché facevano lo stesso lavoro. `computeSourceVe` costava 10.05 ms,
cioè gli stessi ~6.2 ms più 3.9 ms di `sourceVT`/`sourceCV`. Lo stesso dice
`perf`: 43% in `libmutation++` più 26% in `libm`, e il simbolo più caro è
`RrhoDB::gibbs` (11.8%), le costanti di equilibrio, che si raggiungono solo
attraverso `productionRates`.

Ora `updateSources()` fa un giro solo e ne ricava insieme le ω_s di tutte le
specie, il calore di reazione e l'accoppiamento chimica-vibrazione. Il solver
la chiama all'inizio del predictor.

Non è solo più veloce, è anche più corretto: prima la sorgente della specie *i*
vedeva le specie *0…i-1* già avanzate dalle loro `YiEqn`, perché `YiEqn.solve()`
sta **dentro** il ciclo sulle specie, e il calore di reazione le vedeva tutte
avanzate e rinormalizzate. Le eq. 27 e 30 del paper valutano tutti gli ω_s allo
stesso stato.

Effetto sui risultati: nullo in pratica. I cinque casi senza chimica sono
identici bit a bit (senza reazioni non c'era nessuna deriva da correggere); la
fig. 7, l'unico caso reagente, si sposta di **0.011 K su Tv** (0.0002%) e di
1e-7 sulle densità numeriche — il passo di 1 ns è così corto che la
composizione cambia pochissimo dentro un singolo predictor.

## Cosa resta da guadagnare

Dopo il consolidamento OpenFOAM è più della metà del tempo a un thread. Da lì
in poi aggiungere thread rende poco: la leva successiva è MPI sui casi con una
mesh vera, non OpenMP.
