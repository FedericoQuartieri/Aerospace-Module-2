# Relazione di progetto - compilazione
#
#   make        genera main.pdf, la versione inglese (latexmk gestisce da solo
#               le passate e bibtex; figure in images-en/)
#   make ita    genera main-ita.pdf, la versione italiana (figure in immagini/)
#   make pulito rimuove i file ausiliari, tiene i pdf
#   make pulisci rimuove anche i pdf

MAIN = main

.PHONY: all ita pulito pulisci

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex sections-en/*.tex bibliografia.bib assets/logo-polimi.pdf
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex

ita: $(MAIN)-ita.pdf

$(MAIN)-ita.pdf: $(MAIN)-ita.tex sezioni/*.tex bibliografia.bib assets/logo-polimi.pdf
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(MAIN)-ita.tex

pulito:
	latexmk -c $(MAIN).tex
	latexmk -c $(MAIN)-ita.tex

pulisci:
	latexmk -C $(MAIN).tex
	latexmk -C $(MAIN)-ita.tex
