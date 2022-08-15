#!/usr/bin/env bash
echo "self, scpt=${0}"; scpt=${0}
rmd=${1:-uni_grp.rmd}; echo "arg1, rmd=$rmd"

# R -e "rmarkdown::render(\"${rmd}\",output_format=\"html_document\")"
R -e "rmarkdown::render(\"${rmd}\",output_format=\"word_document\")"