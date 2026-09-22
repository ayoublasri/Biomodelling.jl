#!/bin/bash
# Build every document with pandoc. Run from the repository root or from paper/manuscript.
#
# Three documents share one bibliography, one CSL and one number source:
#   .        the integrated preprint (bioRxiv), which carries all eight figures
#   paper1/  the framework paper: Figs 1-7, for PLOS Computational Biology
#   paper2/  the schedule-optimisation paper: the clinical case studies
# fill_numbers.py fills the templates of all three from a single pass over paper/output,
# so a number cannot disagree between them.
set -e
cd "$(dirname "$0")"
ROOT=$PWD
python fill_numbers.py

# $1 = directory, $2 = output stem
build_one() {
  local dir=$1 stem=$2
  ( cd "$dir"
    for f in 00_frontmatter.md 01_intro.md 02_results.md 03_discussion.md 04_methods.md 05_legends.md 06_backmatter.md; do
      [ -f "$f" ] && { cat "$f"; printf '\n\n'; }
    done > "$stem.md"
    local common=(--citeproc --bibliography="$ROOT/references.bib" --csl="$ROOT/plos.csl" --resource-path=".:$ROOT:$ROOT/../figures")
    pandoc "$stem.md" "${common[@]}" -o "$stem.docx" 2>/dev/null || pandoc "$stem.md" "${common[@]}" -o "$stem.docx"
    pandoc "$stem.md" "${common[@]}" -o "$stem.tex" --standalone -H "$ROOT/header.tex"
    pandoc "$stem.md" "${common[@]}" -o "$stem.pdf" --pdf-engine=pdflatex -H "$ROOT/header.tex"
    [ -f supplement.md ] && pandoc supplement.md "${common[@]}" -o "${stem}_supplement.pdf" --pdf-engine=pdflatex -H "$ROOT/header.tex"
    echo "built $dir/$stem.pdf, $dir/$stem.docx${SUPP:+, }$([ -f supplement.md ] && echo "$dir/${stem}_supplement.pdf")"
  )
}

build_one . manuscript
# the preprint keeps its historical supplement filename
[ -f manuscript_supplement.pdf ] && mv manuscript_supplement.pdf supplement.pdf
for d in paper1 paper2; do
  [ -d "$d" ] && [ -f "$d/02_results.md" ] && build_one "$d" "$d"
done
