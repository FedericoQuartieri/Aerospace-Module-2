# A cosa serve `janaf`/`rrho` (e a cosa NON serve)

Stato del codice: branch `revised` (20 settembre 2026). Ogni riferimento sotto
è stato verificato leggendo il codice, non a memoria.

---

## 0. In una frase

`rrho` è **una copia rinominata della classe `janaf` di OpenFOAM** (stesso
polinomio NASA-7). Non calcola la fisica ad alta temperatura: quella la fa
**Mutation++**. `rrho`/`janaf` serve solo per tre cose di contorno — equazione
di stato, controllo del passo temporale, avvio — più **una cosa vera**:
**viscosità e conducibilità del gas**, che Mutation++ non calcola mai.

| # | cosa | chi la usa | dove |
|---|---|---|---|
| 1 | peso molare per specie | equazione di stato `perfectGas` | `specie{molWeight}` |
| 2 | polinomio Cp/h (7 coeff.) | `Cp()` nel controllo del `deltaT`; decode T←e all'avvio | `thermodynamics{...CpCoeffs}` |
| 3 | **coefficienti di trasporto** | **μ e κ veri del solver** (sforzo viscoso, conduzione termica in `e`) | `transport{As,Ts,Amu,...}` |
| 4 | dati spettroscopici RRHO (θrot, θvib, θdiss, g, hf) | **nessuno**: il codice non li legge | `thermodynamics{thetaVib,...}` |

Tutto il resto — le due temperature, la densità, le velocità di reazione, lo
scambio V-T — viene da Mutation++, che ha un **proprio** database interno
anch'esso chiamato "RRHO" ma che non c'entra nulla con la classe OpenFOAM
`rrho`. Il nome è lo stesso (rigid-rotor harmonic-oscillator) per pura
coincidenza terminologica, non perché sia lo stesso codice o lo stesso dato.

---

## 1. Perché si chiama `rrho` se è `janaf`?

`janaf` è già una classe del core di OpenFOAM (thermo con polinomio NASA-7 a
7 coefficienti, due intervalli di temperatura). Serve tuttavia poterla
istanziare in combinazione con l'`equationOfState` e il `transport` che il
progetto usa, ed `janaf` non è fra le combinazioni pre-generate dalle macro
`forGases`/`forCoeffGases` del core. Anziché toccare OpenFOAM, si è copiata la
classe sotto un altro nome e la si è istanziata a mano:

```cpp
// src/thermophysicalModels/multicomponentThermo/highEnthalpyMulticomponentThermo/
// highEnthalpyMulticomponentThermos.C, righe 39-53

// rrho e' una copia di janaf, ma non e' fra i thermo delle macro forGases e
// forCoeffGases del core: lo si istanzia qui con le stesse combinazioni
// (gas perfetto, trasporto const o sutherland, energia h o e)
#define forRrhoGasEqns(Mu, He, Macro, Args...)                                 \
    forThermo(Mu, He, rrhoThermo, perfectGas, specie, Macro, Args)
```

Confronta [`rrhoThermo.H`](../src/thermophysicalModels/specie/thermo/rrho/rrhoThermo.H)
con `janafThermo.H` di OpenFOAM-13: stessa formula, stessi campi
(`Tlow`, `Thigh`, `Tcommon`, `highCpCoeffs`, `lowCpCoeffs`), nomi di file
diversi. **`rrho` non è un modello fisico diverso da `janaf`: è lo stesso
identico modello, con un altro nome**, perché il nome `janaf` era già preso e
serviva una seconda istanza indipendente da modificare/estendere senza
toccare il core.

---

## 2. Chi lo seleziona

`thermoType.thermo` nei `physicalProperties` di ogni caso:

```
// tutorials/shockThermo/shockTube/constant/physicalProperties, righe 17-26
thermoType
{
    type            highEnthalpyThermo;
    mixture         coefficientWilkeMulticomponentMixture;
    transport       sutherland;
    thermo          rrho;               // <- qui
    energy          sensibleInternalEnergy;
    equationOfState perfectGas;
    specie          specie;
}
```

`highEnthalpyThermo` è l'alias di `HighEnthalpyMulticomponentThermo<rrhoThermo<...>>`:
la classe che fa da ponte con Mutation++ **eredita** da `rrho`, ma ne
**riscrive** (override) i metodi che contano davvero — `correct()`,
`updateSources()`, `updatePsi()` — vedi
[`HighEnthalpyMulticomponentThermo.H`](../src/thermophysicalModels/multicomponentThermo/highEnthalpyMulticomponentThermo/HighEnthalpyMulticomponentThermo.H).
`rrho` resta come classe base solo per le funzioni che **non** vengono
riscritte.

---

## 3. Il file dati: `speciesThermo.janaf`

Il file (formato JANAF/NASA-9, dati RRHO da Vincenti-Kruger e Gurvich) ha,
per ogni specie, quattro blocchi. Esempio (N2+):

```
// tutorials/shockThermo/shockTube/constant/speciesThermo.janaf

N2+
{
    specie
    {
        molWeight       28.012855;
    }

    thermodynamics
    {
        Tlow            298.15;
        Thigh           20000;
        Tcommon         10149.075;
        highCpCoeffs    ( -13.61183854 0.005568066047 ... );
        lowCpCoeffs     ( 3.188966155 0.0008739140687 ... );

        //RRHO model
        thetaRot        2.9;
        thetaVib        3390;
        thetaDiss       113000;
        g               ( 2 4 2 4 8 8 4 4 4 4 8 8 4 4 2 2 4 );
        thetaElec       ( 0 1.3189972e4 ... );
        hf              5.3886e7;
    }

    transport
    {
        As              1.907238654e-06;
        Ts              403.2298446;
        Amu             0;   Bmu    2.5;   Cmu   -32.0827;
        Ak              0;   Bk    -0.03723; Ck  0.84192; Dk -3.59040; Ek -18.65620;
    }

    elements { e -1; N 2; }
}
```

Di questi quattro blocchi, **`rrhoThermo.C` (righe 67-71) legge solo**:

```cpp
Tlow_(dict.subDict("thermodynamics").lookup<scalar>("Tlow")),
Thigh_(dict.subDict("thermodynamics").lookup<scalar>("Thigh")),
Tcommon_(dict.subDict("thermodynamics").lookup<scalar>("Tcommon")),
highCpCoeffs_(dict.subDict("thermodynamics").lookup("highCpCoeffs")),
lowCpCoeffs_(dict.subDict("thermodynamics").lookup("lowCpCoeffs"))
```

`thetaRot`, `thetaVib`, `thetaDiss`, `g`, `thetaElec`, `hf` **non sono letti
da nessuna parte nel codice C++** del progetto (verificato con grep su
`applications/` e `src/`, zero occorrenze fuori dai file dati). Sono lì come
documentazione del modello RRHO originale da cui derivano i `CpCoeffs`, non
perché servano a runtime.

`molWeight` lo legge invece l'`equationOfState` (`perfectGas`/`specie`), e
`transport{...}` lo legge `sutherlandTransport`. Nessuno di questi due è
`rrho`, ma vivono nello stesso file perché OpenFOAM organizza i thermo file
per specie, non per modello.

---

## 4. I quattro usi reali a runtime

### 4.1 Equazione di stato — `molWeight`

Dà la costante specifica R = Ru/M di ogni specie a `perfectGas`. È
scaffolding richiesto dal template `FluidMulticomponentThermo<Thermo>`: la
densità/`psi` **nelle celle**, a runtime, viene comunque sovrascritta da
Mutation++ (`updatePsi()`, §5), ma il framework OpenFOAM ne ha comunque
bisogno per compilare ed esistere.

### 4.2 `Cp()` — controllo del passo temporale

```cpp
// applications/modules/shockThermo/setRDeltaT.C, riga 79
mag(reaction->Qdot())/(alphaTemp*rho*thermo.Cp()*thermo.T())
```

Unico posto nel solver dove il `Cp()` polinomiale di `rrho` (non quello di
Mutation++) viene chiamato direttamente.

### 4.3 Decode iniziale T←e, poi buttato via

Il costruttore di `FluidMulticomponentThermo` (classe base, prima che
`HighEnthalpyMulticomponentThermo` prenda controllo) fa il suo normale
decode "energia → temperatura" usando il polinomio `rrho`, limitato al
`Thigh` del file dati (20000 K). Siccome per il caso reale serve poter
partire anche sopra quella soglia, il costruttore del bridge Mutation++
**rilegge T dal caso** subito dopo, buttando via quel valore:

```cpp
// HighEnthalpyMulticomponentThermo.H, righe 330-333
// ---- energie di partenza di ogni cella da p, T e Tve letti dal caso
// il thermo base (janaf) ha gia' ricalcolato T dall'energia limitandola
// al suo intervallo di validita' (20000 K): rileggo T dal file del caso
```

### 4.4 Viscosità e conducibilità — quello che conta davvero

Qui `rrho`/`janaf` non c'entra come classe, ma il **file dati sì**: i
coefficienti `As, Ts` (Sutherland) e `Amu,Bmu,Cmu / Ak,Bk,Ck,Dk,Ek`
(Blottner-style) in `transport{...}` alimentano `sutherlandTransport`
(selezionato in `thermoType.transport`). Ho verificato che **Mutation++ non
viene mai interrogato per μ o κ** (zero chiamate a `mix->viscosity()` o
equivalenti in tutto `applications/` e `src/`). Il μ realmente usato nelle
equazioni del moto e dell'energia è quello di questo file:

```cpp
// applications/modules/shockFluid/shockFluid.C, riga 160
max(thermo_.mu().primitiveField()) > 0
```

```cpp
// applications/modules/shockThermo/thermophysicalPredictor.C, righe 148-154
const surfaceScalarField devTauDotU(
    "devTauDotU",
    devTau() & (a_pos()*U_pos() + a_neg()*U_neg())   // usa mu()
);
EEqn += thermophysicalTransport->divq(e) + fvc::div(devTauDotU);  // usa kappa/alphaEff
```

**Nota a margine**, emersa controllando questo file: l'equazione di `eve`
oggi (righe 173-178 dello stesso file) è solo

```cpp
fvScalarMatrix EveEqn(fvm::ddt(rho, eve) == Q_ve());
```

senza nessun termine di diffusione (`divq`/`laplacian`). Se la memoria del
progetto parla di una diffusione di `eve` aggiunta in una milestone
precedente (κ_ve via Eucken), sul branch `revised` attuale non è presente:
i coefficienti di trasporto di questo file scaldano solo l'equazione di `e`,
non quella di `eve`. Da verificare se sia un regresso o una scelta
consapevole.

---

## 5. Cosa fa Mutation++ al posto di `rrho`

Tutto quello che conta per la fisica a due temperature **non** passa da
`rrho`/`janaf`, ma da chiamate dirette a Mutation++ dentro
`HighEnthalpyMulticomponentThermo`:

| cosa | funzione che lo fa | riga |
|---|---|---|
| T_tr, T_ve da (ρ_s, e, e_ve) | `correct()` → `mix->setState(...,0)`, `mix->T()`, `mix->Tv()` | 542-589 |
| densità/`psi` da (p, T, Tve, Y) | `updatePsi()` → `mix->density()` | 690-736 |
| produzione chimica ω_s | `updateSources()` → `productionRates(...)` | 410-469 |
| scambio V-T (Landau-Teller) | `updateSources()` → `sourceVT(...)` | 441 |
| entalpie di formazione | `initMutation()` → `mix0.speciesHOverRT(...)` | 671-683 |

Il database termodinamico che Mutation++ usa per tutto questo è impostato
qui:

```
// tutorials/shockThermo/shockTube/constant/physicalProperties, righe 28-38
highEnthalpyMutation
{
    mixture                 air_11;
    stateModel              ChemNonEqTTv;
    thermodynamicDatabase   RRHO;     // <- database INTERNO di Mutation++
    mechanism               none;
    ...
}
```

Questo `RRHO` è un database Mutation++ (file XML dentro
`thirdParty/Mutationpp/`, gitignorato), **completamente separato** dalla
classe OpenFOAM `rrho` del punto 1. Stesso nome per la stessa idea fisica
(rigid-rotor harmonic-oscillator), due implementazioni indipendenti che non
si parlano.

---

## 6. Riepilogo: chi legge cosa

| entry nel file dati | letta da | usata per | dove nel codice |
|---|---|---|---|
| `specie.molWeight` | `perfectGas`/`specie` | R = Ru/M | scaffolding EOS |
| `thermodynamics.Tlow/Thigh/Tcommon/*CpCoeffs` | `rrhoThermo.C:67-71` | `Cp()` in `setRDeltaT`; decode T iniziale (poi sovrascritto) | [setRDeltaT.C:79](../applications/modules/shockThermo/setRDeltaT.C#L79) |
| `thermodynamics.thetaRot/thetaVib/thetaDiss/g/thetaElec/hf` | **nessuno** | niente (solo documentazione) | — |
| `transport.As/Ts/Amu/.../Ek` | `sutherlandTransport` | μ e κ **realmente usati** nel solver | [shockFluid.C:160](../applications/modules/shockFluid/shockFluid.C#L160), [thermophysicalPredictor.C:148-154](../applications/modules/shockThermo/thermophysicalPredictor.C#L148-L154) |
| `elements{}` | nessuno (nel nostro codice) | documentazione stechiometria | — |
| — (dati equivalenti dentro Mutation++) | Mutation++ | T_tr, T_ve, densità, ω_s, scambio V-T | [HighEnthalpyMulticomponentThermo.H](../src/thermophysicalModels/multicomponentThermo/highEnthalpyMulticomponentThermo/HighEnthalpyMulticomponentThermo.H) |
