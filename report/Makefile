all: report.pdf


report.pdf: ./comparisons.tex
	pdflatex report.tex

comparisons.tex: ../comparison.json
	./comparisons.py


clear: clean

clean:
	\rm *.pdf || true
	\rm *.aux || true
	\rm *.log || true
