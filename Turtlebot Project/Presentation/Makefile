all: presentation clean

presentation: slides
	latexmk -pdf ./latex/presentation.tex

slides:
	pandoc -t beamer ./markdown/presentation.md -o ./latex/slides.tex

clean:
	rm *.aux *.fdb_latexmk *.fls *.log *.nav *.out *.snm *.toc 