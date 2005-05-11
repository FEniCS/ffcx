FILENAME='ffc-user-manual'

all:
	latex $(FILENAME).tex
	dvips -P pdf -o $(FILENAME).ps $(FILENAME).dvi
	ps2pdf $(FILENAME).ps $(FILENAME).pdf
