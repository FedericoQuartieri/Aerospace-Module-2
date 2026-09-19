# Ambiente di sviluppo

OpenFOAM-13 + Mutation++ + OpenMP in un container Linux, con il progetto
montato da macOS.

Su Apple Silicon l'immagine è **arm64 nativa**: `dl.openfoam.org` pubblica
`openfoam13` anche per arm64, quindi niente emulazione QEMU. Conta perché i
tempi misurati per la parallelizzazione OpenMP sarebbero altrimenti inutili.

## Uso

```sh
./docker/build.sh                 # una volta sola (~15 min: compila Mutation++)
./docker/run.sh                   # shell nel container, progetto su /project
./docker/run.sh ./Allwmake        # compila e esce
OMP_NUM_THREADS=4 ./docker/run.sh # con 4 thread invece di tutti
```

Il progetto è un bind mount, non una copia: quello che si compila dentro
finisce in `platforms/` sul Mac, e le modifiche fatte sul Mac si vedono subito
dentro. Il Dockerfile va ricostruito solo se cambia lui.

`/root` sta invece in un volume Docker (`aero-m2-home`). È lì che OpenFOAM
installa le librerie dell'utente (`$FOAM_USER_LIBBIN`, cioè `libshockThermo.so`
e `libhighEnthalpyThermophysicalModels.so`): senza il volume sparirebbero a
ogni `--rm` e si ricompilerebbe da capo ogni volta. È anche un filesystem
nativo della VM, molto più veloce di quello condiviso con macOS.

## Com'è fatta

- `ubuntu:24.04` + `openfoam13` dal repository della OpenFOAM Foundation
- Mutation++ compilato in `/opt/Mutationpp`, **dentro l'immagine**: il progetto
  lo cercherebbe in `thirdParty/Mutationpp` (vedi `etc/config.sh/mutationpp`),
  ma l'entrypoint reindirizza `MPP_DIRECTORY` su `/opt`. Così la build non
  passa dal filesystem condiviso con macOS, che è lento, e l'immagine resta
  autosufficiente.
- `entrypoint.sh` sorgente prima `/opt/openfoam13/etc/bashrc` (definisce le
  `WM_*`, `FOAM_USER_LIBBIN`, `WM_THIRD_PARTY_DIR`) e poi `/project/etc/bashrc`
  (definisce le `POLIMI_*`), che si appoggia alle prime: l'ordine è obbligato.

## Risorse

La VM di Docker Desktop ha 8 CPU e 4 GB di RAM. I 4 GB sono il motivo per cui
Mutation++ si compila con `-j2`: le unità che includono Eigen prendono oltre
1 GB ciascuna e con `-j8` la build muore con `cannot allocate memory`.

Le 8 CPU sono i core del M2 (4 performance + 4 efficiency): oltre i 4 thread
lo scaling OpenMP resta sotto il lineare perché gli ultimi thread finiscono
sui core efficiency, che sono più lenti. Non è un difetto della
parallelizzazione.

Per cambiare le risorse: Docker Desktop → Settings → Resources.
