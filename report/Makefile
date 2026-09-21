# Relazione di progetto - compilazione
#
#   make        genera main.pdf (latexmk gestisce da solo le passate e bibtex;
#               figure in images-en/)
#   make pulito rimuove i file ausiliari, tiene i pdf
#   make pulisci rimuove anche i pdf

MAIN = main

.PHONY: all pulito pulisci

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex sections-en/*.tex bibliografia.bib assets/logo-polimi.pdf
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex

pulito:
	latexmk -c $(MAIN).tex

pulisci:
	latexmk -C $(MAIN).tex
