############################ -*- Mode: Makefile -*- ###########################
SHELL 		= /bin/bash

figuras = 			\
	DifIndicadoresCLI	\
	IndicadoresCLI		\
	DifIndicadoresCCI	\
	IndicadoresCCI		\
	Graph-CLI-CCI		\
	Graph-CLI-CCI-tip	\

#FIGURAS :=$(addsuffix .pdf,   $(figuras))

xfig2pdflatex: $(figuras)

$(figuras):
	fig2dev -L pdftex   -m .8 -n $(addprefix figuras/, $@.pdf) $@.fig $@.pdf
	fig2dev -L pdftex_t -m .8 -p $(addprefix figuras/, $@.pdf) $@.fig $@.tex

