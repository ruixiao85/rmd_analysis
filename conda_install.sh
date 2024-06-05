
conda create -y -n rmd -c conda-forge r-biocmanager r-tidyverse r-rmarkdown r-ggpubr r-boot r-plotly r-dt

conda activate rmd

R -e 'BiocManager::install(c("jvm"),update=FALSE)'


