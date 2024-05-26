#!/usr/bin/env bash
echo "self, scpt=${0}"; scpt=${0}

apt update
# apt install -y pandoc
apt install -y xfonts-utils \
  texlive-latex-base \
  texlive-latex-extra \
  texlive-fonts-recommended \
  lmodern \
  procps && apt clean

R -e 'BiocManager::install(c("ggpubr","jvm"),update=FALSE)'
R -e 'BiocManager::install(c("rmarkdown"),update=FALSE)'
R -e 'BiocManager::install(c("DT","plotly"),update=FALSE)'
