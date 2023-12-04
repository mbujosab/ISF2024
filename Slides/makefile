############################ -*- Mode: Makefile -*- ###########################
PAPER	   = ISF2024-AGF-MB

NOMBREDIR  = Version0

############
LATEXFILES    = $(VPATH)/DatosDocumento.tex $(VPATH)/texto.tex $(VPATH)/Documento.tex $(VPATH)/referencias/ 
OCGAUXFILES   = $(OCGPATH)/ocg-p.sty $(OCGPATH)/ocgtools.sty # $(OCGPATH)/ocgtools-preview.sh

PREDIR   = $(VPATH)/pre
PREFILES = \
	01preBeamer	\
	01slides	\
	01packages	\
	01AMS		\
	01colores	\
	notacion	\
	notacion2	\
	01notasOCG	\

PREFILESLATEX = $(addsuffix .tex, $(PREFILES))
PREAMBLELATEX = $(addprefix $(PREDIR)/, $(PREFILESLATEX))

VPATH	   = $(PWD)/$(NOMBREDIR)
OCGPATH    = $(PWD)/accesorios
DIRLOCAL   = $(subst $(HOME)/,$(empty),$(PWD))
DIRTEMP	   = $(HOME)/temporal/$(DIRLOCAL)
SHELL	   = /bin/bash


sinocg: ocg = $(DIRTEMP)
sinocg:
	make screen
	make print

all: ocg = $(DIRTEMP)/preview.labels
all: 
	make screen
	make print

screen: $(DIRTEMP) $(LATEXFILES) $(PREAMBLELATEX)
	make $(PAPER)-screen.pdf tipo=screen
	touch $@

print: $(DIRTEMP)  $(LATEXFILES) $(PREAMBLELATEX)
	make $(PAPER)-print.pdf  tipo=print
	touch $@

######################################################

batch          :=                -interaction batchmode
esc_bat        := --shell-escape $(batch)
esc_err        := --shell-escape -interaction errorstopmode

$(PAPER)-screen.pdf: $(DIRTEMP) $(LATEXFILES) $(PAPER)-screen.tex $(ocg) # $(DIRTEMP)/preview.labels # $(ocg) #
$(PAPER)-print.pdf:  $(DIRTEMP) $(LATEXFILES) $(PAPER)-print.tex
$(PAPER)-$(tipo).pdf: %.pdf: %.tex   
	cd $(DIRTEMP) && latexmk -pdf -pdflatex='pdflatex $(esc_err) $(<F)' -use-make $(<F)
#	cd $(DIRTEMP) && latexmk -pdf -pdflatex='pdflatex $(esc_bat) $(<F)' -use-make $(<F)
	cp $(DIRTEMP)/$@ .

$(PAPER)-$(tipo).tex: $(PREDIR)/$(tipo).tex $(PREAMBLELATEX) $(VPATH)/DatosDocumento.tex $(VPATH)/Documento.tex
	cat $? >> $(DIRTEMP)/$@

$(DIRTEMP)/preview.labels: $(LATEXFILES) $(PAPER)-screen.tex
	cd $(DIRTEMP) && pdflatex $(esc_err) "\AtBeginDocument{\usepackage[pdftex,displaymath,floats,active,tightpage]{preview}}\PassOptionsToPackage{noocg}{ocgtools} \def\ocgpreview#1#2{#2}\input $(PAPER)-screen.tex "
#	cd $(DIRTEMP) && pdflatex $(esc_bat) "\AtBeginDocument{\usepackage[pdftex,displaymath,floats,active,tightpage]{preview}}\PassOptionsToPackage{noocg}{ocgtools} \def\ocgpreview#1#2{#2}\input $(PAPER)-screen.tex "
	cd $(DIRTEMP) && cp $(PAPER)-screen.pdf  preview-temp.pdf
	cd $(DIRTEMP) && grep '\\newlabel{' ocgtools-preview.aux | sed 's/newlabel/maplabelstoall/g' > preview.labels

######################################################

$(DIRTEMP): $(LATEXFILES) $(OCGAUXFILES);
	mkdir -p $@
	ln -s -f $? $@
	ln -s -f $(VPATH)/figuras/ $@


clean:          ; rm -r print screen

cleanall: 
	rm -r $(DIRTEMP); rm -r $(PAPER)-*
	make clean
