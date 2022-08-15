read_tb<-function(file,...){
  ext<-tolower(tools::file_ext(file)); ext1<-substr(ext,1,1)
  # readRDS(file)
  if (ext1=='x') return(read_excel(file,...))
  if (ext1=='c') return(read_csv(file,...))
  return(read_tsv(file,...))
}
is_valid<-function(x) !is.na(x) & trimws(x)!=""
p_sym<-function(p) {
  if (is.na(p)||p>0.1) return("ns")
  if (p>0.05) return(".")
  if (p>0.01) return("*")
  if (p>0.001) return("**")
  return("***")
}
ps_sym<-function(ps) sapply(ps, FUN=p_sym)
ps_fmt<-function(ps) paste(format.pval(ps, eps=0.001, digits=2), ps_sym(ps))
str_wrap<-function(ss, lim=12) {
  sapply(ss, function(x) paste(strwrap(gsub(" |_|/"," ",x),width=lim),collapse="\n"))
}
fac_rev<-function(x) factor(x, levels=sort(unique(x),decreasing=T))
fac_wrap_rev<-function(x) fac_rev(str_wrap(x))
is_cont<-function(vec) {
  type<-purrr::map_chr(vec, class)[1]
  return (type=="numeric" && length(unique(vec))>4)
}
pct_true<-function(x,i) { sum(x[i])/length(x) }
freq_pct<-function(list) {
  tab<-as.data.frame(table(list,dnn=""),row.names=1)
  if (nrow(tab)<1) return(NULL)
  if (nrow(tab)==1) { row.names(tab)<-tab[,1]; tab<-tab[-1] }
  tab$Pct<-round(100.0*tab$Freq/sum(tab$Freq),1)
  return(tab)
}
chisq_multi<-function(mat, p.adj="fdr", verbose=F) {
  mat<-mat %>% select(Group, Category, n)
  catu<-unique(mat$Category)
  grpc<-combn(unique(mat$Group), 2)
  resm<-tibble(Comparison=character(), Category=character(), p=numeric())
  mat<-mat %>% group_by(Group) %>%
    summarize(across(), other=sum(n), .groups='keep') %>% mutate(other=other-n)
  for (ci in 1:length(catu)) {
    matc<-mat %>% filter(Category==catu[ci]) %>% select(-Category)
    for (gi in 1:ncol(grpc)) {
      matcg<-matc %>% filter(Group %in% grpc[,gi]) %>%
        column_to_rownames(var="Group")
      cst<-chisq.test(matcg) # print(matcg)
      if (verbose) {
        cat(paste(paste(grpc[1,gi],grpc[2,gi],sep="-"),catu[ci],sep=" @ "),end="\n")
        print(cst)
      }
      resm<-resm %>% add_row(Comparison=paste(grpc[1,gi],grpc[2,gi],sep="-"),
        Category=catu[ci],p=cst$`p.value`)
    }
  }
  if (is.null(p.adj)) {
    resm$signif<-ps_sym(resm$p)
  } else {
    resm[[p.adj]]<-p.adjust(resm$p, method=p.adj)
    resm$signif<-ps_sym(resm[[p.adj]])
  }
  return(resm)
}


