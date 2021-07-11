view:
	open main.pdf

all: plots pdf
.PHONY: all arxiv

plots:
	cd generator; Rscript main.R

pdf:
	latexmk -pdf main.tex

arxiv:
	## cp main.tex arxiv
	cp -r generated arxiv
	rm arxiv/generated/README.txt
	cp literature_zotero.bib arxiv
	cd arxiv & latexmk -pdf main.tex	

diff:
	latexdiff-vc --disable-citation-markup --flatten --git --force -r 0e43d77a834beb75ea01cecf8c5726061867b520 main.tex ## needs manual labour: remove bib and put normal bib command + figure captions

diffrev:
	latexdiff-vc --disable-citation-markup --flatten --git --force -r 57820aacbbeddddb8c5523606b80d71f8df23f6f main.tex ## needs manual labour: remove bib and put normal bib command + figure captions

review:
	pandoc response_to_the_reviewers.md --pdf-engine=pdflatex -o response.pdf

clean:  ## Clean output files
	latexmk -c main.tex
	rm -rf generator/memoise


