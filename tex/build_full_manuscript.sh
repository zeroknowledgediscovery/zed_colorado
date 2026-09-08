#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

latexmk \
  -pdf \
  -interaction=nonstopmode \
  -halt-on-error \
  zebra_genomics_boundary_theorem_empirical.tex

latexmk \
  -pdf \
  -interaction=nonstopmode \
  -halt-on-error \
  zebra_genomics_boundary_theorem_empirical_supplement.tex

echo "Built: $(pwd)/zebra_genomics_boundary_theorem_empirical.pdf"
echo "Built: $(pwd)/zebra_genomics_boundary_theorem_empirical_supplement.pdf"
