#!/usr/bin/env Rscript
gr_arg<-function(a,k) sub(k,"",a[grep(k,a)])
cat("Running ",scpt<-gr_arg(commandArgs(),"--file="),"\n")
args=commandArgs(trailingOnly=TRUE); nargs=length(args) # print(args)
cat("arg1: fin [",fin<-if (nargs>0) args[1] else "cl0905.tsv","]\n")
cat("arg2: grp [",grp<-if (nargs>1) args[2] else "GCW,GWW,GWE,GWI,LAEF,LAVImax,LAVImin","]\n")

# docker run -it --rm -v $PWD:$PWD -w $PWD rocker/verse bash
library(tidyverse)
grps<-unlist(strsplit(grp,","))

# fbs<-tools::file_path_sans_ext(fin)
# ext1<-tolower(substr(tools::file_ext(fin),1,1))
# if (ext1=='c') { df<-read_csv(fin)
# } else if (ext1=='x') { df<-readxl::read_excel(fin)
# } else { df<-read_tsv(fin) }

# dfs<-read_tsv(fin,na=c("","NA","0"))
dfs<-read_tsv(fin,na=c("","NA"))



# BiocManager::install("PerformanceAnalytics",update=F); library(PerformanceAnalytics)
# df<-dfs[,grps]
# pdf("corplot_Pearson_Spearman.pdf")
# chart.Correlation(df, histogram=TRUE, method="pearson", pch=19)
# chart.Correlation(df, histogram=TRUE, method="spearman", pch=19)
# dev.off()


BiocManager::install("ggpubr",update=F); library(ggpubr)

for (sm in c("spearman")) {
# for (sm in c("pearson","kendall")) {
  # sink(paste0("ggpubr_",sm,".tsv"))
  # cat("GroupSelected","\t","Pair","\t","CorCoef_",sm,"\t","PValue","\n",sep="")
  for (ss in c("Male","Female")) {
    df<-dfs[dfs$gender==ss,]
    # for (sg in c("control","PAF","all")) {
    for (sg in c("PAF")) {
      if (sg!="all") { df<-df[df$group==sg,] }
      df<-df[,grps]
      pdf(paste0("ggpubr_",sm,"_",sg,"_",ss,".pdf"), width=5, height=4, family="Times")
      for (x in grps) {
        for (y in grps) {
          # if (x!=y) {
          if (substr(x,1,2)=="LA" && substr(y,1,1)=="G") {
            # pdf(paste0("ggpubr_",sm,"_",sg,"_",x,"-",y,".pdf"), width=5, height=4, family="Times")
            # png(paste0("ggpubr_",sm,"_",sg,"_",x,"-",y,".png"), width=5, height=4, units="in", res=200)
            # prm=list(data=df, x=x, y=y)
            # do.call(ggscatter, prm)
            print(ggscatter(df, x=x, y=y, add="reg.line", conf.int=TRUE,
              color="darkblue", shape=1, add.params=list(color="darkblue", fill="lightblue"),
              cor.coef=TRUE, cor.method=sm)+border(color="darkblue")) # pearson
            res<-cor.test(df[[x]], df[[y]], method=sm, conf.level=0.95)
            cat(sg,"\t",x,"-",y,"\t",res$estimate,"\t",res$p.value,"\n",sep="")
            # break
          }
        }
      }
      dev.off()
    }
  }
  # sink()
}

