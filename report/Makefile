# Relazione di progetto - compilazione
#
#   make        genera main.pdf (latexmk gestisce da solo le passate e bibtex)
#   make pulito rimuove i file ausiliari, tiene il pdf
#   make pulisci rimuove anche il pdf

MAIN = main

.PHONY: all pulito pulisci

all: $(MAIN).pdf

$(MAIN).pdf: $(MAIN).tex sezioni/*.tex bibliografia.bib
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex

pulito:
	latexmk -c

pulisci:
	latexmk -C
