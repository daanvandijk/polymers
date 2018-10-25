all: report.pdf


report.pdf: ./comparisons.tex
	pdflatex report.tex >> /dev/null

comparisons.tex: ../comparison.json
	./comparisons.py


clear: clean

clean:
	\rm *.pdf || true
	\rm *.aux || true
	\rm *.log || true
	rm comparisons.tex
