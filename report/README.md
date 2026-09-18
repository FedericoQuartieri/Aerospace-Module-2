# Relazione di progetto

Relazione sul modello a due temperature implementato in questo repository
(branch `revised`) e sulla riproduzione dei casi di verifica 0D di
Casseau et al., *Aerospace* **3** (2016) 34.

## Compilazione

```sh
make            # produce main.pdf
make pulito     # rimuove i file ausiliari
```

Serve una installazione TeX Live con `latexmk` e `bibtex`.

## Struttura

| file | contenuto |
|------|-----------|
| `main.tex` | preambolo, frontespizio, inclusione delle sezioni |
| `sezioni/00-sommario.tex` | sommario |
| `sezioni/01-introduzione.tex` | introduzione |
| `sezioni/02-lavoro.tex` | il lavoro svolto |
| `sezioni/03-risultati.tex` | verifica e risultati |
| `sezioni/04-discussione.tex` | discussione |
| `sezioni/05-limiti.tex` | limiti e sviluppi futuri |
| `bibliografia.bib` | riferimenti |

Le sezioni si aggiungono una alla volta: ogni nuova sezione e' un file in
`sezioni/` incluso con `\input` in `main.tex`.

## Figure

Le figure dei risultati **non** sono duplicate qui: `\graphicspath` in
`main.tex` punta direttamente a `../applications/test/nonEqTTv/output/` e a
`../explainations/dipendenze/`, cosi' il pdf usa sempre l'ultima versione
prodotta da `Allrun`/`compare.py`. La cartella `immagini/` e' per le figure
scritte apposta per la relazione.
