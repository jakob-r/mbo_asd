view:
	open main.pdf

all: plots pdf
.PHONY: all

plots:
	cd generator; Rscript main.R

pdf:
	latexmk -pdf main.tex

clean:  ## Clean output files
	latexmk -c main.tex
	rm -rf generator/memoise
