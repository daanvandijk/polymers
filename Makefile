ODIR = ./obj/
SDIR = ./src/
BDIR = ./bin/

all: cpp

publish: cpp
	rm publish -rf
	mkdir publish/ || true
	mkdir publish/bin || true
	mkdir publish/img || true
	mkdir publish/data || true
	cp bin/main publish/bin/main
	cp img/ publish/ -r
	rm "publish/img/__pycache__" -rf
	rm publish/img/*.png
	rm publish/img/coparisons -rf
	zip -r publish.zip publish/

experiments : cpp
	cd data && ../bin/main "../experiments.json" experiment

list:
	bin/main "experiments.json" list

cpp: 
	$(MAKE) -C ./src all 

img: cpp
	$(MAKE) -C ./img all

test: cpp 
	$(BDIR)test

clear: clean
clean:
	$(MAKE) -C ./src clean
	$(MAKE) -C ./img clean
	\rm *.pdf || true
	rm publish/ || true -rf
	rm publish.zip || true
clean-data:
	\rm ./data/*.dat

#latex: main.pdf presentation.pdf
#main.pdf: main.tex
	#pdflatex main.tex > /dev/null
