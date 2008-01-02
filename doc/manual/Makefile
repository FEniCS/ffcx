FILENAME=ufc-user-manual

all:
	latex $(FILENAME).tex

final:
	latex $(FILENAME).tex
	bibtex $(FILENAME)
	latex $(FILENAME).tex
	makeindex $(FILENAME)
	latex $(FILENAME)
	dvips -P pdf -o $(FILENAME).ps $(FILENAME).dvi
	ps2pdf $(FILENAME).ps $(FILENAME).pdf

clean:
	rm -f $(FILENAME).aux
	rm -f $(FILENAME).idx
	rm -f $(FILENAME).log
	rm -f $(FILENAME).out
	rm -f $(FILENAME).toc
	rm -f $(FILENAME).bbl
	rm -f $(FILENAME).blg
	rm -f $(FILENAME).ilg
	rm -f $(FILENAME).ind

purge: clean
	rm -f $(FILENAME).dvi
	rm -f $(FILENAME).ps
	rm -f $(FILENAME).pdf

