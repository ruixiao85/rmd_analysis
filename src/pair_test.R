#!/usr/bin/env Rscript
gr_arg<-function(a,k) sub(k,"",a[grep(k,a)])
cat("Running ",scpt<-gr_arg(commandArgs(),"--file="),"\n")
args=commandArgs(trailingOnly=TRUE); nargs=length(args) # print(args)
cat("arg1: fin [",fin<-if (nargs>0) args[1] else NULL,"]\n")
cat("arg2: grp [",grp<-if (nargs>1) args[2] else 1,"]\n")

# docker run -it --rm -v $PWD:$PWD -w $PWD rocker/tidyverse bash
# library(tidyverse)
# BiocManager::install("ggpubr"); library(ggpubr)
# BiocManager::install("jmv"); library(jmv)


fbs<-tools::file_path_sans_ext(fin)
ext1<-tolower(substr(tools::file_ext(fin),1,1))
if (ext1=='c') { df<-read_csv(fin)
} else if (ext1=='x') { df<-readxl::read_excel(fin)
} else { df<-read_tsv(fin) }

# dfs<-read_tsv(fin,na=c("","NA","0"))

df<-dfs[,c(1,8:29)]


fac<-unique(df[[grp]])
# df[[grp]]<-factor(df[[grp]],levels=c("Control","PAF"))
nf<-length(fac)
hf<-colnames(df)[grp]

# fac<-c("A","B","C")
pair<-list()
for (l in seq(1,nf,1)) {
  if (l<nf) {
    for (r in seq(l+1,nf,1)) {
      cat(l,fac[l],"\t",r,fac[r],"\n")
      # pair<-append(pair,list(c(fac[l],fac[r])))
      pr<-list(c(fac[l],fac[r])); names(pr)<-paste(l,r,sep="-")
      pair<-append(pair,pr)
    }
  }
}

sink(paste0(fbs,"_res.csv")); sep=","
cat(paste(hf,nf,sep=sep,end="\n"))
for (pn in names(pair)) {
  cat(pn,sep="",end=",")
  cat(pair[[pn]],sep=sep,end="\n")
}
pdf(paste0(fbs,".pdf"))
for (ci in seq(1,ncol(df))){
  if (ci!=grp) { # escape grp_idx
    hy<-colnames(df)[ci]
    p<-ggboxplot(df, x=hf, y=hy,
      color=hf, add="jitter", shape=hf)
    p<-p+stat_compare_means(comparisons=pair, method="t.test")
    if (nf>2) p<-p+stat_compare_means(label.y=50, method="anova")
    print(p)
    cat(hy,sep)
    for (f in fac) {
      vec<-df[[ci]][df[[grp]]==f]
      vec<-vec[!is.na(vec)] # remove NA
      cat(f,length(vec),mean(vec),sd(vec),median(vec),"",sep=sep)
    }
    for (pn in names(pair)) {
      ps<-pair[[pn]]
      t_welch<-t.test(df[[ci]][df[[hf]]==ps[1]],df[[ci]][df[[hf]]==ps[2]])
      cat(pn,t_welch$p.value,"",sep=sep)
    }
    cat("\n")
  }
}
sink()
dev.off()


if (FALSE) {
dgt=2
param=list(Group=c("control","disease"),gender=c("F","M"))
for (pn in names(param)){
  df[[pn]]<-factor(df[[pn]],levels=param[[pn]])
  p1=param[[pn]][1]; p2=param[[pn]][2]; cat(pn,":",p1,p2,end="\n")
  pdf(paste0(pn,"_paired_ttest2.pdf"),height=5,width=5)
  for (ci in 3:ncol(df)) {
  # for (ci in 3:4) {
    p<-ggboxplot(df, x=pn, y=colnames(df)[ci],
        color=pn, palette =c("#00AFBB", "#FC4E07"),
        add="jitter", shape=pn,
        xlab=paste0(round(mean(df[df[[pn]]==p1,ci,drop=T],na.rm=T),dgt),
          "±",round(sd(df[df[[pn]]==p1,ci,drop=T],na.rm=T),dgt)," \t",
      " ",round(mean(df[df[[pn]]==p2,ci,drop=T],na.rm=T),dgt),"±",
        round(sd(df[df[[pn]]==p2,ci,drop=T],na.rm=T),dgt),"\n")
      )
    # print(p+stat_compare_means(method="t.test"))
    print(p+stat_compare_means(method="t.test",paired=T))
    # print(p+stat_compare_means(comparisons=list(p=param[[p]]))
  }
  dev.off()
}
library(ggpubr)
for (ci in 3:ncol(p)) {
  d<-data.frame(pre=p[p$GRP=="PAF1",ci,drop=T],post=p[p$GRP=="PAF2",ci,drop=T])
  tp<-t.test(d$pre, d$post, var.equal = TRUE, paired = TRUE)
  print(ggpaired(d,cond1="pre",cond2="post",fill="condition",line.color="gray",pallette="jco",
  xlab=paste0("paired t-test P=",round(tp$p.value,3)),ylab=colnames(p)[ci]))
}


install.packages("pwr")
library(pwr)
pwr_rc=function(r1,r2,c){
  delta<-abs(pa[r1,c]-pa[r2,c])
  sigma<-mean(ps[r1,c],ps[r2,c])
  ratio<-delta/sigma
  # cat(pwr.t.test(d=ratio, sig.level=.05, power = .90, type = 'two.sample')$n,endl="\n")
  cat(pwr.t.test(d=ratio, sig.level=.05, power = .80, type = 'two.sample')$n,endl="\n")
}
for (c in 2:5) {
  cat(colnames(p)[c],endl="\n")
  cat(" cAF vs Ctrl N = "); pwr_rc(1,2,c) # cAF vs Ctrl GWI
  cat(" pAF vs Ctrl N = "); pwr_rc(2,3,c) # pAF vs Ctrl GWI N=25
  cat(" cAF vs pAF N = "); pwr_rc(1,3,c) # pAF vs cAF GWI
}

BiocManager::install("ggpubr")
library("ggpubr")
library("tidyverse")
df<-read_tsv("ttr2.tsv")
dgt=2
param=list(Group=c("control","disease"),gender=c("F","M"))
for (pn in names(param)){
  df[[pn]]<-factor(df[[pn]],levels=param[[pn]])
  p1=param[[pn]][1]; p2=param[[pn]][2]; cat(pn,":",p1,p2,end="\n")
  pdf(paste0(pn,"_paired_ttest2.pdf"),height=5,width=5)
  for (ci in 3:ncol(df)) {
  # for (ci in 3:4) {
    p<-ggboxplot(df, x=pn, y=colnames(df)[ci],
        color=pn, palette =c("#00AFBB", "#FC4E07"),
        add="jitter", shape=pn,
        xlab=paste0(round(mean(df[df[[pn]]==p1,ci,drop=T],na.rm=T),dgt),
          "±",round(sd(df[df[[pn]]==p1,ci,drop=T],na.rm=T),dgt)," \t",
      " ",round(mean(df[df[[pn]]==p2,ci,drop=T],na.rm=T),dgt),"±",
        round(sd(df[df[[pn]]==p2,ci,drop=T],na.rm=T),dgt),"\n")
      )
    # print(p+stat_compare_means(method="t.test"))
    print(p+stat_compare_means(method="t.test",paired=T))
    # print(p+stat_compare_means(comparisons=list(p=param[[p]]))
  }
  dev.off()
}
}