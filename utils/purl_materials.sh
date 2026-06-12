#!/usr/env bash

# check if we are in the project directory (contains `_quarto.yml` file)
if [ ! -f "_quarto.yml" ]; then
    echo "Error: This script must be run from the project root directory (where _quarto.yml is located)."
    exit 1
fi

# create output directory
mkdir -p course_files/scripts

# extract R code from all .qmd files in the materials directory and save to course_files/scripts
for f in materials/*.qmd; do
    # if file does not contain `{r}` chunks, skip it
    if ! grep -q "{r" "$f"; then
        echo "Skipping $f (no R chunks found)"
        continue
    fi

    # otherwise purl
    out="course_files/scripts/$(basename "${f%.qmd}.R")"

    Rscript -e "
        knitr::purl(
            '$f',
            output = '$out',
            documentation = 0
        )
    "
done