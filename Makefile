ODIR = ./obj/
SDIR = ./src/
BDIR = ./bin/

all: cpp

experiments : cpp
	cd data && ../bin/main experiment

cpp: 
	$(MAKE) -C ./src all 

img: cpp
	$(MAKE) -C ./img all

test: cpp 
	$(BDIR)test

#latex: main.pdf presentation.pdf
#main.pdf: main.tex
	#pdflatex main.tex > /dev/null

clear: clean
clean:
	$(MAKE) -C ./src clean
	$(MAKE) -C ./img clean
	\rm *.pdf || true
clean-data:
	\rm ./data/*.dat
