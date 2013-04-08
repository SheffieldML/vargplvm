FILE = 2012_CVPR_del

tex: *.tex
	pdflatex $(FILE)
bib: *.bib *.tex
	bibtex $(FILE)
all: *.tex *.bib
	$(bib) $(tex)

