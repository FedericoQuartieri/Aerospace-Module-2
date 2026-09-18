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
| `sezioni/05-conclusioni.tex` | conclusioni e sviluppi futuri |
| `bibliografia.bib` | riferimenti |

Le sezioni si aggiungono una alla volta: ogni nuova sezione e' un file in
`sezioni/` incluso con `\input` in `main.tex`.

## Figure

I grafici della sezione 3 sono in `immagini/` e si rigenerano con

```sh
python3 immagini/genera-grafici.py
```

Lo script legge i risultati del progetto (`../applications/test/nonEqTTv/output/`)
e le curve estratte dal paper (`../applications/test/nonEqTTv/paper-data/`):
sono le stesse curve dei grafici di `compare.py`, ma senza titolo, con la
legenda in italiano e con nomi che non richiamano la numerazione del paper.
