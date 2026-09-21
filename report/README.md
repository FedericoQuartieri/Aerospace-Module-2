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

La relazione e' in inglese: `main.tex`, con le sezioni in `sections-en/` e le
figure in `images-en/`.

| file | contenuto |
|------|-----------|
| `main.tex` | preambolo, frontespizio, inclusione delle sezioni di `sections-en/` |
| `sections-en/00-summary.tex` | abstract |
| `sections-en/01-introduction.tex` | introduzione |
| `sections-en/02-work.tex` | modello e implementazione |
| `sections-en/03-results.tex` | verifica e risultati |
| `sections-en/04-discussion.tex` | discussione |
| `sections-en/05-conclusions.tex` | conclusioni e sviluppi futuri |
| `bibliografia.bib` | riferimenti |

Le sezioni si aggiungono una alla volta: ogni nuova sezione e' un file in
`sections-en/` incluso con `\input` in `main.tex`.

## Figure

I grafici della sezione 3 si rigenerano con

```sh
python3 images-en/genera-grafici.py
```

Lo script legge i risultati del progetto (`../applications/test/nonEqTTv/output/`)
e le curve estratte dal paper (`../applications/test/nonEqTTv/paper-data/`):
sono le stesse curve dei grafici di `compare.py`, ma senza titolo e con nomi
che non richiamano la numerazione del paper.
