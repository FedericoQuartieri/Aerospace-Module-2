"""Confronta i campi finali delle corse a diverso numero di thread.

Due controlli distinti, perche' misurano cose diverse.

1. RIPRODUCIBILITA' a parita' di thread. Con schedule(static) la ripartizione
   delle celle fra i thread e' fissa, quindi due corse con lo stesso numero di
   thread devono dare gli stessi bit. Se non li danno c'e' una corsa critica:
   e' cosi' che si e' scoperto che Mutation++ teneva i buffer del Newton in
   variabili static di funzione, condivise fra tutti i thread (vedi
   docker/mutationpp-thread-safety.patch).

2. ACCORDO fra 1 thread e N thread. Qui l'uguaglianza bit a bit NON e'
   attesa, ed e' importante capire perche' non lo e': solveEnergies() di
   Mutation++ e' un Newton che parte dalle temperature correnti della
   miscela, cioe' da quelle dell'ultima cella che quel thread ha trattato.
   Con un thread la catena e' cella 0, 1, 2, ...; con otto ogni thread ha la
   sua. Punti di partenza diversi, stessa soluzione a meno della tolleranza
   del Newton (rtol 1e-12). Lo scarto atteso e' quello, non zero.

uso: python3 check.py scalingBath "1 2 4 8"
"""

import re
import sys
from pathlib import Path

# soglia di accordo fra numeri di thread diversi: molto sopra la tolleranza
# del Newton di Mutation++ (1e-12) e molto sotto qualunque effetto fisico
TOLLERANZA = 1.0e-9


def read_field(path, n_expected=None):
    """Valori interni di un volScalarField, anche se scritto come uniform."""
    text = path.read_text()
    body = text[text.index("internalField"):]

    uniform = re.match(r"internalField\s+uniform\s+([^;]+);", body)
    if uniform:
        value = float(uniform.group(1))
        return [value] * (n_expected or 1)

    n = int(re.search(r"nonuniform\s+List<scalar>\s*\n?\s*(\d+)", body).group(1))
    start = body.index("(", body.index("List<scalar>"))
    end = body.index(")", start)
    values = [float(v) for v in body[start + 1:end].split()]
    assert len(values) == n, f"{path}: attesi {n} valori, trovati {len(values)}"
    return values


def relative_difference(values, reference):
    """Massimo scarto relativo fra due campi."""
    scale = max(abs(v) for v in reference) or 1.0
    return max(abs(a - b) for a, b in zip(values, reference)) / scale


def main():
    case = Path(sys.argv[1])
    threads = [int(t) for t in sys.argv[2].split()]
    fields = ["T", "Tve", "N2", "N"]
    base = threads[0]

    # quante celle: la mesh puo' aver scritto un campo come uniform
    sizes = [len(read_field(case / f"result-{n}" / "T")) for n in threads]
    nCells = max(sizes)

    reference = {}
    worst = 0.0

    print(f"--- accordo con la corsa a {base} thread "
          f"(atteso: la tolleranza del Newton, ~1e-12) ---")
    for name in fields:
        reference[name] = read_field(case / f"result-{base}" / name, nCells)
        line = [f"{name:>4}: {reference[name][0]:.10g}"]
        for n in threads[1:]:
            values = read_field(case / f"result-{n}" / name, nCells)
            diff = relative_difference(values, reference[name])
            worst = max(worst, diff)
            line.append(f"{n} thread {diff:.1e}")
        print("   " + "   ".join(line))

    # riproducibilita': stessa corsa due volte con lo stesso numero di thread
    print()
    repeat = case / f"result-{threads[-1]}-bis"
    reproducible = None
    if repeat.exists():
        print(f"--- riproducibilita' a {threads[-1]} thread "
              f"(atteso: identico bit a bit) ---")
        reproducible = True
        for name in fields:
            a = read_field(case / f"result-{threads[-1]}" / name, nCells)
            b = read_field(repeat / name, nCells)
            same = a == b
            reproducible = reproducible and same
            print(f"   {name:>4}: {'identico' if same else 'DIVERSO'}")

    print()
    ok = worst < TOLLERANZA and reproducible is not False
    if ok:
        print(f"OK: scarto massimo fra numeri di thread diversi {worst:.1e}, "
              f"sotto la soglia {TOLLERANZA:.0e}")
        if reproducible:
            print("OK: a parita' di thread il risultato e' riproducibile "
                  "bit a bit")
        return 0

    if worst >= TOLLERANZA:
        print(f"ERRORE: scarto {worst:.1e} sopra la soglia {TOLLERANZA:.0e}: "
              "il risultato dipende dal numero di thread")
    if reproducible is False:
        print("ERRORE: due corse con lo stesso numero di thread danno "
              "risultati diversi -> corsa critica")
    return 1


if __name__ == "__main__":
    sys.exit(main())
