#!/bin/bash
# Build the manuscript (PDF and DOCX) with pandoc. Run from the repository root or from paper/manuscript.
set -e
cd "$(dirname "$0")"
cat 00_frontmatter.md 01_intro.md 02_results.md 03_discussion.md 04_methods.md 05_legends.md 06_backmatter.md > manuscript.md
pandoc manuscript.md --citeproc --bibliography=references.bib --csl=plos.csl -o manuscript.docx --resource-path=.:../figures 2>/dev/null || \
pandoc manuscript.md --citeproc --bibliography=references.bib -o manuscript.docx --resource-path=.:../figures
pandoc manuscript.md --citeproc --bibliography=references.bib -o manuscript.tex --standalone --resource-path=.:../figures -H header.tex
pandoc manuscript.md --citeproc --bibliography=references.bib -o manuscript.pdf --pdf-engine=pdflatex --resource-path=.:../figures -H header.tex
pandoc supplement.md --citeproc --bibliography=references.bib -o supplement.pdf --pdf-engine=pdflatex --resource-path=.:../figures -H header.tex
echo "built manuscript.pdf, manuscript.docx, supplement.pdf"
