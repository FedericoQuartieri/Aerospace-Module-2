# Relazione di progetto

Relazione sul modello a due temperature implementato in questo repository
(branch `revised`) e sulla riproduzione dei casi di verifica 0D di
Casseau et al., *Aerospace* **3** (2016) 34.

## Compilazione

```sh
make            # produce main.pdf, la versione inglese
make ita        # produce main-ita.pdf, la versione italiana
make pulito     # rimuove i file ausiliari
```

Serve una installazione TeX Live con `latexmk` e `bibtex`.

## Struttura

La versione principale e' in inglese (`main.tex`, sezioni in `sections-en/`,
figure in `images-en/`); quella italiana e' `main-ita.tex`, con le sezioni in
`sezioni/` e le figure in `immagini/`. Le due versioni si aggiornano a mano:
una modifica a una sezione va riportata nel file corrispondente dell'altra
lingua.

| file | contenuto |
|------|-----------|
| `main.tex` | versione inglese: preambolo, frontespizio, inclusione delle sezioni di `sections-en/` |
| `main-ita.tex` | versione italiana: preambolo, frontespizio, inclusione delle sezioni di `sezioni/` |
| `sezioni/00-sommario.tex` | sommario |
| `sezioni/01-introduzione.tex` | introduzione |
| `sezioni/02-lavoro.tex` | il lavoro svolto |
| `sezioni/03-risultati.tex` | verifica e risultati |
| `sezioni/04-discussione.tex` | discussione |
| `sezioni/05-conclusioni.tex` | conclusioni e sviluppi futuri |
| `bibliografia.bib` | riferimenti |

Le sezioni si aggiungono una alla volta: ogni nuova sezione e' un file in
`sections-en/` incluso con `\input` in `main.tex`, e il suo corrispondente
italiano in `sezioni/`, incluso in `main-ita.tex`.

## Figure

I grafici della sezione 3 si rigenerano con

```sh
python3 immagini/genera-grafici.py --en   # legenda in inglese, in images-en/ (main.tex)
python3 immagini/genera-grafici.py        # legenda in italiano, in immagini/ (main-ita.tex)
```

Lo script legge i risultati del progetto (`../applications/test/nonEqTTv/output/`)
e le curve estratte dal paper (`../applications/test/nonEqTTv/paper-data/`):
sono le stesse curve dei grafici di `compare.py`, ma senza titolo, con la
legenda in italiano e con nomi che non richiamano la numerazione del paper.
